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
from specdb.build.utils import init_data

#igms_path = imp.find_module('igmspec')[1]


def grab_meta():
    """ Grab KODIAQ meta Table
    Returns
    -------

    """
    hstz2_meta = Table.read(os.getenv('RAW_IGMSPEC')+'/HST_z2/hst_z2.ascii', format='ascii')
    nspec = len(hstz2_meta)
    # RA/DEC, DATE
    ra = []
    dec = []
    for row in hstz2_meta:
        # Fix DEC
        # Get RA/DEC
        coord = ltu.radec_to_coord((row['ra'],row['dec']))
        ra.append(coord.ra.value)
        dec.append(coord.dec.value)
    hstz2_meta.add_column(Column(ra, name='RA_GROUP'))
    hstz2_meta.add_column(Column(dec, name='DEC_GROUP'))
    # z
    hstz2_meta.rename_column('zem', 'zem_GROUP')
    hstz2_meta['sig_zem'] = [0.]*nspec
    hstz2_meta['flag_zem'] = [str('SDSS_PIPE')]*nspec
    hstz2_meta['STYPE'] = [str('QSO')]*nspec
    #
    hstz2_meta.rename_column('obsdate','DATE-OBS')
    hstz2_meta.rename_column('tel','TELESCOPE')
    hstz2_meta.rename_column('inst','INSTR')
    hstz2_meta.rename_column('grating','GRATING')
    hstz2_meta.rename_column('resolution','R')
    # Check
    assert chk_meta(hstz2_meta, chk_cat_only=True)
    # Return
    return hstz2_meta

'''
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
'''


def hdf5_adddata(hdf, sname, meta, debug=False, chk_meta_only=False,
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
    # Checks
    if sname != 'HST_z2':
        raise IOError("Not expecting this survey..")

    # Build spectra (and parse for meta)
    nspec = len(meta)
    max_npix = 300  # Just needs to be large enough
    data =init_data(max_npix)
    # Init
    spec_set = hdf[sname].create_dataset('spec', data=data, chunks=True,
                                         maxshape=(None,), compression='gzip')
    spec_set.resize((nspec,))
    wvminlist = []
    wvmaxlist = []
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
        #head = hduf[0].header
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
    meta.add_column(Column(np.arange(nspec,dtype=int),name='GROUP_ID'))

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

