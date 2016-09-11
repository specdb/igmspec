""" Module to ingest HST+FUSE AGN spectra

Cooksey et al. 2010
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import pdb
import warnings
import os, json

from astropy.table import Table, Column
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
from astropy.io.fits import Header

from linetools.spectra import io as lsio
from linetools import utils as ltu

from igmspec.ingest import utils as iiu
from igmspec import defs

#igms_path = imp.find_module('igmspec')[1]


def grab_meta():
    """ Grab KODIAQ meta Table
    Returns
    -------

    """
    hstc_file = os.getenv('RAW_IGMSPEC')+'/HST_Cooksey/HSTQSO_pre-SM4.lst'
    hstc_meta = Table.read(hstc_file, format='ascii')
    #
    hstc_meta.add_column(Column([2000.]*len(hstc_meta), name='EPOCH'))
    return hstc_meta

def meta_for_build():
    """ Generates the meta data needed for the IGMSpec build
    Returns
    -------
    meta : Table
    """
    # Cut down to unique QSOs
    hstc_meta = grab_meta()
    names = hstc_meta['QSO'].data
    uni, uni_idx = np.unique(names, return_index=True)
    hstc_meta = hstc_meta[uni_idx]
    nqso = len(hstc_meta)
    #
    meta = Table()
    meta['RA'] = hstc_meta['RA']
    meta['DEC'] = hstc_meta['DEC']
    meta['zem'] = hstc_meta['zem']
    meta['sig_zem'] = [0.]*nqso
    meta['flag_zem'] = [str('UNKWN')]*nqso
    meta['STYPE'] = [str('QSO')]*nqso
    # Return
    return meta


def hdf5_adddata(hdf, IDs, sname, debug=False, chk_meta_only=False,
                 mk_test_file=False):
    """ Append HST/FUSE data to the h5 file

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
    Rdicts = defs.get_res_dicts()
    # Add Survey
    print("Adding {:s} survey to DB".format(sname))
    hstc_grp = hdf.create_group(sname)
    # Load up
    meta = grab_meta()
    bmeta = meta_for_build()
    # Checks
    if sname != 'HST_Cooksey':
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
        raise ValueError("Bad matches in HST_Cooksey")
    meta_IDs = IDs[idx]

    # Loop me to bid the full survey catalog
    meta.add_column(Column(meta_IDs, name='IGM_ID'))

    # Build spectra (and parse for meta)
    nspec = len(meta)
    max_npix = 30000  # Just needs to be large enough
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
    Rlist = []
    wvminlist = []
    wvmaxlist = []
    npixlist = []
    gratinglist = []
    # Loop
    #path = os.getenv('RAW_IGMSPEC')+'/KODIAQ_data_20150421/'
    path = os.getenv('RAW_IGMSPEC')+'/HST_Cooksey/'
    maxpix = 0
    for jj,row in enumerate(meta):
        # Generate full file
        full_file = path+'{:s}/{:s}/{:s}'.format(
                         row['QSO'],row['INSTR'],row['SPEC_FILE'])
        # Extract
        if row['INSTR'] == 'FUSE':
            hext = 1
        else:
            hext = 0
        print("HST_Cooksey: Reading {:s}".format(full_file))
        spec = lsio.readspec(full_file, head_ext=hext)
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
        if row['INSTR'] == 'FUSE':
            for ss in range(len(spec.header['HISTORY'])):
                card = Header(spec.header['HISTORY'][ss])
                pdb.set_trace()
                if card.keys()[0] == 'APERTURE':
                    pdb.set_trace()
                    aper = Header(spec.header['HISTORY'][ss])
        gratinglist.append(spec.header['GRATING'])
        Rlist.append(Rdicts[row['INSTR']][gratinglist[-1]])
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
    meta.add_column(Column(npixlist, name='NPIX'))
    meta.add_column(Column(wvminlist, name='WV_MIN'))
    meta.add_column(Column(wvmaxlist, name='WV_MAX'))
    meta.add_column(Column(np.arange(nspec,dtype=int),name='SURVEY_ID'))

    # Add HDLLS meta to hdf5
    if iiu.chk_meta(meta):
        if chk_meta_only:
            pdb.set_trace()
        hdf[sname]['meta'] = meta
    else:
        raise ValueError("meta file failed")
    # References
    refs = [dict(url='http://adsabs.harvard.edu/abs/2010ApJ...708..868C',
                 bib='cooksey10')
            ]
    jrefs = ltu.jsonify(refs)
    hdf[sname]['meta'].attrs['Refs'] = json.dumps(jrefs)
    #
    return

