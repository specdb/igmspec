""" Module to ingest HD-LLS Survey data

Prochaska et al. 2015
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb
import os

from astropy.table import Table, Column
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
from astropy.io import fits


def meta_for_build():
    """ Generates the meta data needed for the IGMSpec build
    Returns
    -------
    meta : Table
    """
    hdlls_meta = Table.read(os.getenv('RAW_IGMSPEC')+'/HD-LLS_DR1.fits.gz')
    # Cut down to unique QSOs
    names = np.array([name[0:26] for name in hdlls_meta['Name']])
    uni, uni_idx = np.unique(names, return_index=True)
    hdlls_meta = hdlls_meta[uni_idx]
    nqso = len(hdlls_meta)
    #
    meta = Table()
    meta['RA'] = hdlls_meta['RA']
    meta['DEC'] = hdlls_meta['DEC']
    meta['EPOCH'] = [2000.]*nqso
    meta['zem'] = hdlls_meta['Z_QSO']
    meta['sig_zem'] = [0.]*nqso
    meta['flag_zem'] = [str('UNKN')]*nqso
    # Return
    return meta


def hdf5_adddata(hdf, IDs, sname, debug=False):
    """ Append HD-LLS data to the h5 file

    Parameters
    ----------
    hdf : hdf5 pointer
    IDs : ndarray
      int array of IGMsp_ID values in mainDB
    sname : str
      Survey name

    Returns
    -------

    """
    #import tarfile
    # Add Survey
    hdlls_grp = hdf.create_group(sname)
    # Checks
    if sname != 'HD-LLS_DR1':
        raise IOError("Not expecting this survey..")
    if np.sum(IDs < 0) > 0:
        raise ValueError("Bad ID values")
    # Open Meta tables
    cut_meta = meta_for_build()
    if len(cut_meta) != len(IDs):
        raise ValueError("Wrong sized table..")
    hdlls_meta = Table.read(os.getenv('RAW_IGMSPEC')+'/HD-LLS_DR1.fits.gz')
    nlls = len(hdlls_meta)

    # Generate ID array from RA/DEC
    c_cut = SkyCoord(ra=cut_meta['RA'], dec=cut_meta['DEC'], unit='deg')
    c_all = SkyCoord(ra=hdlls_meta['RA'], dec=hdlls_meta['DEC'], unit='deg')
    # Find new sources
    idx, d2d, d3d = match_coordinates_sky(c_all, c_cut, nthneighbor=1)
    if np.sum(d2d > 0.1*u.arcsec):
        raise ValueError("Bad matches in HD-LLS")
    meta_IDs = IDs[idx]

    # Kludgy Table judo
    hdlls_full = hdlls_meta[0:1]
    spec_files = []

    full_IDs = []
    # Loop me to bid the full survey catalog
    for kk,row in enumerate(hdlls_meta):
        for spec_file in row['SPEC_FILES']:
            if spec_file == 'NULL':
                continue
            # Add to full table
            hdlls_full.add_row(row)
            spec_files.append(spec_file)
            full_IDs.append(meta_IDs[kk])
    # Trim down
    hdlls_full = hdlls_full[1:]
    hdlls_full.remove_column('SPEC_FILES')
    hdlls_full.add_column(Column(spec_files,name='SPEC_FILE'))
    hdlls_full.add_column(Column(full_IDs, name='IGMsp_ID'))
    uni, uni_idx = np.unique(np.array(spec_files), return_index=True)
    hdlls_full = hdlls_full[uni_idx]


    # Build spectra
    if debug:
        nspec = 10
    else:
        nspec = len(hdlls_full)
    max_npix = 210000
    data = np.ma.empty((1,), #self.npix),  # THIS WILL NOT SCALE FOR BOSS!!
                       dtype=[(str('wave'), 'float64', (max_npix)),
                              (str('flux'), 'float32', (max_npix)),
                              (str('sig'),  'float32', (max_npix)),
                              #(str('co'),   'float32', (max_npix)),
                             ])
    # Init
    full_idx = np.zeros(len(hdlls_full), dtype=int)
    spec_set = hdf[sname].create_dataset('spec', data=data, chunks=True,
                                         maxshape=(None,), compression='gzip')
    spec_set.resize((nspec,))

    import glob
    members = glob.glob(os.getenv('RAW_IGMSPEC')+'/{:s}/*fits'.format(sname))
    for jj,member in enumerate(members):
    #tar = tarfile.open(os.getenv('RAW_IGMSPEC')+'/HD-LLS_DR1_spectra.tar.gz')
    #for jj,member in enumerate(tar.getmembers()):
        kk = jj
        """
        if '.' not in member.name:
            print('Skipping a likely folder: {:s}'.format(member.name))
            continue
        if jj == 0:
            raise ValueError("HD-LLS: Did not expect to get here")
        """
        # Extract
        #f = tar.extractfile(member)
        f = member
        hdu = fits.open(f)
        # Parse name
        #fname = f.name.split('/')[-1]
        fname = f.split('/')[-1]
        mt = np.where(hdlls_full['SPEC_FILE'] == fname)[0]
        if len(mt) != 1:
            pdb.set_trace()
            raise ValueError("HD-LLS: No match to spectral file?!")
        else:
            print('loading {:s}'.format(fname))
            full_idx[kk] = mt[0]
        # npix
        npix = hdu[0].header['NAXIS1']
        if npix > max_npix:
            raise ValueError("Too many pixels... ({:d})".format(npix))
        # Some fiddling about
        for key in ['wave','flux','sig']:
            data[key] = 0.  # Important to init (for compression too)
        data['flux'][0][:npix] = hdu[0].data
        data['sig'][0][:npix] = hdu[1].data
        data['wave'][0][:npix] = hdu[2].data
        # Only way to set the dataset correctly
        spec_set[kk] = data

    #tar.close()

    # Add HDLLS meta to hdf5
    hdf[sname]['meta'] = hdlls_full[full_idx]  # No special dataset necessary?
    #
    return


