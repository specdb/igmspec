""" Module to ingest HD-LLS Survey data

Prochaska et al. 2015
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import pdb
import warnings
import os, json

from astropy.table import Table, Column
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
from astropy.io import fits

from linetools.spectra import io as lsio
from linetools import utils as ltu

from specdb.build.utils import chk_meta

#igms_path = imp.find_module('igmspec')[1]


def grab_meta():
    """ Grab KODIAQ meta Table
    Returns
    -------

    """
    hstz2_meta = Table.read(os.getenv('RAW_IGMSPEC')+'/HST_z2/hst_z2.ascii', format='ascii')
    # RA/DEC, DATE
    ra = []
    dec = []
    for row in hstz2_meta:
        # Fix DEC
        # Get RA/DEC
        coord = ltu.radec_to_coord((row['ra'],row['dec']))
        ra.append(coord.ra.value)
        dec.append(coord.dec.value)
    hstz2_meta.add_column(Column(ra, name='RA'))
    hstz2_meta.add_column(Column(dec, name='DEC'))
    # RENAME
    hstz2_meta.rename_column('obsdate','DATE-OBS')
    hstz2_meta.rename_column('tel','TELESCOPE')
    hstz2_meta.rename_column('inst','INSTR')
    hstz2_meta.rename_column('grating','GRATING')
    hstz2_meta.rename_column('resolution','R')
    return hstz2_meta

def meta_for_build():
    """ Generates the meta data needed for the IGMSpec build
    Returns
    -------
    meta : Table
    """
    # Cut down to unique QSOs
    hstz2_meta = grab_meta()
    names = np.array([name[0:26] for name in hstz2_meta['qso']])
    uni, uni_idx = np.unique(names, return_index=True)
    hstz2_meta = hstz2_meta[uni_idx]
    nqso = len(hstz2_meta)
    #
    meta = Table()
    meta['RA'] = hstz2_meta['RA']
    meta['DEC'] = hstz2_meta['DEC']
    meta['zem'] = hstz2_meta['zem']
    meta['sig_zem'] = [0.]*nqso
    meta['flag_zem'] = [str('SDSS_PIPE')]*nqso
    meta['STYPE'] = [str('QSO')]*nqso
    # Return
    return meta


def hdf5_adddata(hdf, IDs, sname, debug=False, chk_meta_only=False,
                 mk_test_file=False):
    """ Append HST_z2 data to the h5 file

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
    hstz2_grp = hdf.create_group(sname)
    # Load up
    meta = grab_meta()
    bmeta = meta_for_build()
    # Checks
    if sname != 'HST_z2':
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
        raise ValueError("Bad matches in HST_z2")
    meta_IDs = IDs[idx]

    # Loop me to bid the full survey catalog
    meta.add_column(Column(meta_IDs, name='IGM_ID'))

    # Build spectra (and parse for meta)
    nspec = len(meta)
    max_npix = 300  # Just needs to be large enough
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
    gratinglist = []
    npixlist = []
    speclist = []
    # Loop
    #path = os.getenv('RAW_IGMSPEC')+'/KODIAQ_data_20150421/'
    path = os.getenv('RAW_IGMSPEC')+'/HST_z2/'
    maxpix = 0
    for jj,row in enumerate(meta):
        # Generate full file
        if row['INSTR'] == 'ACS':
            full_file = path+row['qso']+'.fits.gz'
        elif row['INSTR'] == 'WFC3':
            coord = ltu.radec_to_coord((row['RA'],row['DEC']))
            full_file = path+'/J{:s}{:s}_wfc3.fits.gz'.format(coord.ra.to_string(unit=u.hour,sep='',precision=2,pad=True),
                                               coord.dec.to_string(sep='',pad=True,alwayssign=True,precision=1))
        # Extract
        print("HST_z2: Reading {:s}".format(full_file))
        hduf = fits.open(full_file)
        head = hduf[0].header
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
    meta.add_column(Column([2000.]*nspec, name='EPOCH'))
    meta.add_column(Column(speclist, name='SPEC_FILE'))
    meta.add_column(Column(npixlist, name='NPIX'))
    meta.add_column(Column(wvminlist, name='WV_MIN'))
    meta.add_column(Column(wvmaxlist, name='WV_MAX'))
    meta.add_column(Column(np.arange(nspec,dtype=int),name='SURVEY_ID'))

    # Add HDLLS meta to hdf5
    if chk_meta(meta):
        if chk_meta_only:
            pdb.set_trace()
        hdf[sname]['meta'] = meta
    else:
        raise ValueError("meta file failed")
    # References
    refs = [dict(url='http://adsabs.harvard.edu/abs/2011ApJS..195...16O',
                 bib='omeara11')
            ]
    jrefs = ltu.jsonify(refs)
    hdf[sname]['meta'].attrs['Refs'] = json.dumps(jrefs)
    #
    return

