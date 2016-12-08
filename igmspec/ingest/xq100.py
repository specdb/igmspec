""" Module to ingest XQ-100 Survey data

Lopez et al. 2016
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import pdb
import os, glob
import imp
import json

from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.table import Table, Column, vstack
from astropy.time import Time
from astropy.io import fits
from astropy import units as u

from linetools.spectra import io as lsio
from linetools import utils as ltu

from specdb.build.utils import chk_meta

igms_path = imp.find_module('igmspec')[1]


def grab_meta():
    """ Grab XQ-100 meta Table

    Returns
    -------

    """
    #
    xq100_table = Table.read(os.getenv('RAW_IGMSPEC')+'/XQ-100/XQ100_v1_2.fits.gz')
    nqso = len(xq100_table)
    # ESO meta
    eso_tbl = Table.read(os.getenv('RAW_IGMSPEC')+'/XQ-100/metadata_eso_XQ100.csv', format='ascii.csv')
    ar_files = eso_tbl['ARCFILE'].data
    # Spectral files
    spec_files = glob.glob(os.getenv('RAW_IGMSPEC')+'/XQ-100/ADP.*')
    # Dummy column
    xq100_coords = SkyCoord(ra=xq100_table['RA'], dec=xq100_table['DEC'], unit='deg')
    matches = []
    sv_spec_files = []
    sv_orig_files = []
    sv_rescale_files = []
    for spec_file in spec_files:
        if 'ADP.2016-07-15T08:22:40.682.fits' in spec_file:
            print("XQ-100: Skipping summary file")
            continue
        # ESO file
        ssfile = spec_file[spec_file.rfind('/')+1:-5]
        eso_mt = np.where(ar_files == ssfile)[0]
        try:
            ofile = eso_tbl['ORIGFILE'][eso_mt][0]
        except IndexError:
            print("XQ-100: File {:s} not really in XQ100!".format(spec_file))
            continue
        if ('_1' in ofile) or ('_2' in ofile) or ('_3' in ofile) or ('_4' in ofile):
            print("XQ-100: Skipping additional file: {:s}".format(ofile))
            continue
        # Match
        hdu = fits.open(spec_file)
        head0 = hdu[0].header
        if head0['DISPELEM'] == 'UVB,VIS,NIR':
            print("XQ-100: Skipping merged spectrum file")
            if 'rescale' not in ofile:
                print('no rescale')
                pdb.set_trace()
            continue
        try:
            coord = SkyCoord(ra=head0['RA'], dec=head0['DEC'], unit='deg')
        except KeyError:
            pdb.set_trace()
        sep = coord.separation(xq100_coords)
        imt = np.argmin(sep)
        if sep[imt] > 0.1*u.arcsec:
            pdb.set_trace()
            raise ValueError("Bad offset")
        # Save
        matches.append(imt)
        sv_spec_files.append(spec_file)
        sv_orig_files.append(ofile)
    # Finish up
    xq100_meta = xq100_table[np.array(matches)]
    nspec = len(xq100_meta)
    # Add spec_files
    xq100_meta['SPEC_FILE'] = sv_spec_files
    xq100_meta['ORIG_FILE'] = sv_orig_files
    # Add zem
    xq100_meta['zem_GROUP'] = xq100_meta['Z_QSO']
    xq100_meta['sig_zem'] = xq100_meta['ERR_ZQSO']
    xq100_meta['flag_zem'] = [str('XQ-100')]*nspec
    # Rename
    xq100_meta.rename_column('RA','RA_GROUP')
    xq100_meta.rename_column('DEC','DEC_GROUP')
    # DATE-OBS
    meanmjd = []
    for row in xq100_meta:
        gdm = row['MJD_OBS'] > 0.
        meanmjd.append(np.mean(row['MJD_OBS'][gdm]))
    t = Time(meanmjd, format='mjd', out_subfmt='date')  # Fixes to YYYY-MM-DD
    xq100_meta.add_column(Column(t.iso, name='DATE-OBS'))
    #
    xq100_meta.add_column(Column([2000.]*nspec, name='EPOCH'))
    xq100_meta['STYPE'] = str('QSO')
    # Sort
    xq100_meta.sort('RA_GROUP')
    # Check
    assert chk_meta(xq100_meta, chk_cat_only=True)
    #
    return xq100_meta


'''
def meta_for_build(xq100_meta=None):
    """ Generates the meta data needed for the IGMSpec build
    Returns
    -------
    meta : Table
    """
    if xq100_meta is None:
        xq100_meta = grab_meta()
    # Cut down to unique QSOs
    names = np.array([name[0:20] for name in xq100_meta['OBJ_NAME']])
    uni, uni_idx = np.unique(names, return_index=True)
    xq100_meta = xq100_meta[uni_idx]
    nqso = len(xq100_meta)
    #
    meta = Table()
    for key in ['RA', 'DEC', 'zem', 'sig_zem', 'flag_zem']:
        meta[key] = xq100_meta[key]
    meta['STYPE'] = [str('QSO')]*nqso
    # Return
    return meta
'''


def hdf5_adddata(hdf, IDs, sname, debug=False, chk_meta_only=False):
    """ Append XQ-100 data to the h5 file

    Parameters
    ----------
    hdf : hdf5 pointer
    IDs : ndarray
      int array of IGM_ID values in mainDB
    sname : str
      Survey name
    chk_meta_only : bool, optional
      Only check meta file;  will not write

    Returns
    -------

    """
    # Add Survey
    print("Adding {:s} survey to DB".format(sname))
    xq100_grp = hdf.create_group(sname)
    # Load up
    meta = grab_meta()
    bmeta = meta_for_build(xq100_meta=meta)
    eso_meta = Table.read(os.getenv('RAW_IGMSPEC')+'/XQ-100/metadata_eso_XQ100.csv', format='ascii.csv')
    if len(meta) != 300:
        pdb.set_trace()
    # Checks
    if sname != 'XQ-100':
        raise IOError("Expecting XQ-100!!")
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
        raise ValueError("Bad matches in XQ-100")
    meta_IDs = IDs[idx]
    meta.add_column(Column(meta_IDs, name='IGM_ID'))

    # Build spectra (and parse for meta)
    nspec = len(meta)
    max_npix = 20000  # Just needs to be large enough
    data = np.ma.empty((1,),
                       dtype=[(str('wave'), 'float64', (max_npix)),
                              (str('flux'), 'float32', (max_npix)),
                              (str('sig'),  'float32', (max_npix)),
                              (str('co'),   'float32', (max_npix)),
                             ])
    # Init
    spec_set = hdf[sname].create_dataset('spec', data=data, chunks=True,
                                         maxshape=(None,), compression='gzip')
    spec_set.resize((nspec,))
    Rlist = []
    wvminlist = []
    wvmaxlist = []
    npixlist = []
    speclist = []
    gratinglist = []
    telelist = []
    instrlist = []
    # Loop
    maxpix = 0
    for jj,row in enumerate(meta):
        #
        print("XQ-100: Reading {:s}".format(row['SPEC_FILE']))
        spec = lsio.readspec(row['SPEC_FILE'])
        # Parse name
        fname = row['SPEC_FILE'].split('/')[-1]
        # npix
        npix = spec.npix
        if npix > max_npix:
            raise ValueError("Not enough pixels in the data... ({:d})".format(npix))
        else:
            maxpix = max(npix,maxpix)
        # Continuum
        # Some fiddling about
        for key in ['wave','flux','sig']:
            data[key] = 0.  # Important to init (for compression too)
        data['flux'][0][:npix] = spec.flux.value
        data['sig'][0][:npix] = spec.sig.value
        data['wave'][0][:npix] = spec.wavelength.to('AA').value
        data['co'][0][:npix] = spec.co.value
        # Meta
        head = spec.header
        speclist.append(str(fname))
        wvminlist.append(np.min(data['wave'][0][:npix]))
        wvmaxlist.append(np.max(data['wave'][0][:npix]))
        telelist.append(head['TELESCOP'])
        instrlist.append(head['INSTRUME'])
        gratinglist.append(head['DISPELEM'])
        npixlist.append(npix)
        if gratinglist[-1] == 'NIR':  # From Lopez+16
            Rlist.append(4350.)
        elif gratinglist[-1] == 'VIS':
            Rlist.append(7450.)
        elif gratinglist[-1] == 'UVB':
            Rlist.append(5300.)
        else:
            pdb.set_trace()
            raise ValueError("UH OH")
        # Only way to set the dataset correctly
        if chk_meta_only:
            continue
        spec_set[jj] = data

    #
    print("Max pix = {:d}".format(maxpix))
    # Add columns
    meta.add_column(Column(gratinglist, name='GRATING'))
    meta.add_column(Column(telelist, name='TELESCOPE'))
    meta.add_column(Column(instrlist, name='INSTR'))
    meta.add_column(Column(npixlist, name='NPIX'))
    meta.add_column(Column(wvminlist, name='WV_MIN'))
    meta.add_column(Column(wvmaxlist, name='WV_MAX'))
    meta.add_column(Column(Rlist, name='R'))
    meta.add_column(Column(np.arange(nspec,dtype=int),name='SURVEY_ID'))

    # Add HDLLS meta to hdf5
    if chk_meta(meta):
        if chk_meta_only:
            pdb.set_trace()
        hdf[sname]['meta'] = meta
    else:
        raise ValueError("meta file failed")

    # References
    refs = [dict(url='http://adsabs.harvard.edu/abs/2016arXiv160708776L',
                 bib='lopez+16')]
    jrefs = ltu.jsonify(refs)
    hdf[sname]['meta'].attrs['Refs'] = json.dumps(jrefs)
    #
    return


