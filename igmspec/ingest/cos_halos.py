""" Module to ingest COS-Halos

Tumlinson et al. 2013
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import pdb
import warnings
import os, json, glob, imp

from astropy.table import Table, Column
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
from astropy.io import fits
from astropy.time import Time

from linetools.spectra import io as lsio
from linetools import utils as ltu

from igmspec.ingest import utils as iiu

#igms_path = imp.find_module('igmspec')[1]


def grab_meta():
    """ Grab COS-Halos meta table
    Returns
    -------

    """
    from igmspec.cat_utils import zem_from_radec
    from igmspec.igmspec import IgmSpec
    igmsp = IgmSpec(db_file=os.getenv('IGMSPEC_DB')+'/IGMspec_DB_v01.hdf5', skip_test=True)

    summ_file = os.getenv('RAW_IGMSPEC')+'/COS-Halos/cos_halos_obs.ascii'
    chalos_meta = Table.read(summ_file, format='ascii')
    # RA/DEC, DATE
    ra = []
    dec = []
    for row in chalos_meta:
        coord = ltu.radec_to_coord(row['QSO'])
        ra.append(coord.ra.value)
        dec.append(coord.dec.value)
    chalos_meta.add_column(Column(ra, name='RA'))
    chalos_meta.add_column(Column(dec, name='DEC'))
    chalos_meta.add_column(Column(['2010-01-01']*len(chalos_meta), name='DATE-OBS'))
    # Others
    chalos_meta.add_column(Column(['HST']*len(chalos_meta), name='TELESCOPE'))
    chalos_meta.add_column(Column(['COS']*len(chalos_meta), name='INSTR'))
    chalos_meta.add_column(Column(['G130M/G160M']*len(chalos_meta), name='GRATING'))
    chalos_meta.add_column(Column([20000.]*len(chalos_meta), name='R'))
    chalos_meta.add_column(Column([2000.]*len(chalos_meta), name='EPOCH'))
    # Myers for zem
    zem, zsource = zem_from_radec(chalos_meta['RA'], chalos_meta['DEC'], igmsp.idb.hdf)
    badz = zem <= 0.
    if np.sum(badz) > 0:
        raise ValueError("Bad zem in COS-Halos")
    chalos_meta['zem'] = zem
    chalos_meta['sig_zem'] = 0.  # Need to add
    chalos_meta['flag_zem'] = zsource
    # Done
    return chalos_meta


def meta_for_build():
    """ Generates the meta data needed for the IGMSpec build
    Returns
    -------
    meta : Table
    """
    # Cut down to unique QSOs
    chalos_meta = grab_meta()
    #
    meta = Table()
    for key in ['RA', 'DEC', 'zem', 'sig_zem', 'flag_zem']:
        meta[key] = chalos_meta[key]
    meta['STYPE'] = str('QSO')
    # Return
    return meta


def hdf5_adddata(hdf, IDs, sname, debug=False, chk_meta_only=False,
                 mk_test_file=False):
    """ Append COS-Halos data to the h5 file

    Parameters
    ----------
    hdf : hdf5 pointer
    IDs : ndarray
      int array of IGM_ID values in mainDB
    sname : str
      Survey name
    chk_meta_only : bool, optional
      Only check meta file;  will not write
    mk_test_file : bool, optional
      Generate the debug test file for Travis??

    Returns
    -------

    """
    # Add Survey
    print("Adding {:s} survey to DB".format(sname))
    chalos_grp = hdf.create_group(sname)
    # Load up
    meta = grab_meta()
    bmeta = meta_for_build()
    # Checks
    if sname != 'COS-Halos':
        raise IOError("Not expecting this survey..")
    if np.sum(IDs < 0) > 0:
        raise ValueError("Bad ID values")
    # Open Meta tables
    if len(bmeta) != len(IDs):
        raise ValueError("Wrong sized table..")

    # Generate ID array from RA/DEC
    c_cut = SkyCoord(ra=bmeta['RA'], dec=bmeta['DEC'], unit='deg')
    c_all = SkyCoord(ra=meta['RA'], dec=meta['DEC'], unit='deg')
    # Find new sources
    idx, d2d, d3d = match_coordinates_sky(c_all, c_cut, nthneighbor=1)
    if np.sum(d2d > 0.1*u.arcsec):
        raise ValueError("Bad matches in COS-Halos")
    meta_IDs = IDs[idx]

    # Loop me to bid the full survey catalog
    meta.add_column(Column(meta_IDs, name='IGM_ID'))

    # Build spectra (and parse for meta)
    nspec = len(meta)
    max_npix = 20000  # Just needs to be large enough
    data = np.ma.empty((1,),
                       dtype=[(str('wave'), 'float64', (max_npix)),
                              (str('flux'), 'float32', (max_npix)),
                              (str('sig'),  'float32', (max_npix)),
                              #(str('co'),   'float32', (max_npix)),
                             ])
    # Init
    spec_set = hdf[sname].create_dataset('spec', data=data, chunks=True,
                                         maxshape=(None,), compression='gzip')
    spec_set.resize((nspec,))
    wvminlist = []
    wvmaxlist = []
    npixlist = []
    speclist = []
    # Loop
    path = os.getenv('RAW_IGMSPEC')+'/COS-Halos/'
    maxpix = 0
    for jj,row in enumerate(meta):
        # Generate full file
        coord = ltu.radec_to_coord((row['RA'],row['DEC']))
        full_file = path+'/J{:s}{:s}_nbin3_coadd.fits.gz'.format(coord.ra.to_string(unit=u.hour,sep='',pad=True)[0:4],
                                               coord.dec.to_string(sep='',pad=True,alwayssign=True)[0:5])
        # Extract
        print("COS-Halos: Reading {:s}".format(full_file))
        spec = lsio.readspec(full_file)
        # Parse name
        fname = full_file.split('/')[-1]
        # npix
        npix = spec.npix
        if npix > max_npix:
            raise ValueError("Not enough pixels in the data... ({:d})".format(npix))
        else:
            maxpix = max(npix,maxpix)
        # Some fiddling about
        for key in ['wave','flux','sig']:
            data[key] = 0.  # Important to init (for compression too)
        data['flux'][0][:npix] = spec.flux.value
        data['sig'][0][:npix] = spec.sig.value
        data['wave'][0][:npix] = spec.wavelength.value
        # Meta
        speclist.append(str(fname))
        wvminlist.append(np.min(data['wave'][0][:npix]))
        wvmaxlist.append(np.max(data['wave'][0][:npix]))
        npixlist.append(npix)
        if chk_meta_only:
            continue
        # Only way to set the dataset correctly
        spec_set[jj] = data

    #
    print("Max pix = {:d}".format(maxpix))
    # Add columns
    meta.add_column(Column(speclist, name='SPEC_FILE'))
    meta.add_column(Column(npixlist, name='NPIX'))
    meta.add_column(Column(wvminlist, name='WV_MIN'))
    meta.add_column(Column(wvmaxlist, name='WV_MAX'))
    meta.add_column(Column(np.arange(nspec,dtype=int), name='SURVEY_ID'))

    # Add HDLLS meta to hdf5
    if iiu.chk_meta(meta):
        if chk_meta_only:
            pdb.set_trace()
        hdf[sname]['meta'] = meta
    else:
        raise ValueError("meta file failed")
    # References
    refs = [dict(url='http://adsabs.harvard.edu/abs/2013ApJ...777...59T',
                 bib='tumlinson+13')
            ]
    jrefs = ltu.jsonify(refs)
    hdf[sname]['meta'].attrs['Refs'] = json.dumps(jrefs)
    #
    return

