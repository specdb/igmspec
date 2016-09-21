""" Module to ingest HD-LLS Survey data

Prochaska et al. 2015
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import pdb
import warnings
import os, json, glob, imp
import datetime

from astropy.table import Table, Column, vstack
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
from astropy.io import fits

from linetools.spectra import io as lsio
from linetools import utils as ltu

from specdb.build.utils import chk_meta

#igms_path = imp.find_module('igmspec')[1]


def grab_meta():
    """ Grab HSTQSO meta Table
    Returns
    -------

    """
    summ_file = os.getenv('RAW_IGMSPEC')+'/HSTQSO/hstqso.lst'
    hstqso_meta = Table.read(summ_file, format='ascii')
    nspec = len(hstqso_meta)
    # RA/DEC
    radec_file = os.getenv('RAW_IGMSPEC')+'/HSTQSO/all_qso_table.txt'
    radec = Table.read(radec_file, format='ascii')
    # DATE-OBS
    date_files = glob.glob(os.getenv('RAW_IGMSPEC')+'/HSTQSO/date_obs*')
    for ss,date_file in enumerate(date_files):
        if ss == 0:
            tab_date = Table.read(date_file, format='ascii')
        else:
            tab_date = vstack([tab_date, Table.read(date_file, format='ascii')])
    # RA/DEC, DATE
    hstqso_meta.add_column(Column(['2000-01-01']*nspec, name='DATE-OBS'))
    for jj,row in enumerate(hstqso_meta):
        if row['INST'] == 'COS':
            continue
        # DATE
        spec = row['SPEC_FILE'].split('.')[0]
        mt1 = np.where(tab_date['SPEC'] == spec)[0]
        if len(mt1) == 0:
            print("NO DATE MATCH for {:s}!".format(spec))
            pdb.set_trace()
        else:
            mt1 = mt1[0] # TAKING THE FIRST ONE
        joe_date = tab_date['DATE-OBS'][mt1].split('-')
        hstqso_meta[jj]['DATE-OBS'] = '{:s}-{:02d}-{:02d}'.format(joe_date[0], int(joe_date[1]), int(joe_date[2]))
        if int(joe_date[1]) > 12:
            pdb.set_trace()
        # RA/DEC
        if row['INST'] != 'FOS':
            continue
        mt = np.where(radec['File_ID'] == row['QSO_ALT_NAME'])[0]
        if len(mt) == 0:
            mt = np.where(radec['File_ID'] == row['QSO_NAME'])[0]
            if len(mt) == 0:
                print("NO RA/DEC MATCH!")
                pdb.set_trace()
            else:
                mt = mt[0]
        else:
            mt = mt[0]
        hstqso_meta[jj]['RA'] = radec['RA'][mt]
        hstqso_meta[jj]['DEC'] = radec['DEC'][mt]
    # TEST
    from astropy.time import Time
    try:
        tval = Time(list(hstqso_meta['DATE-OBS'].data), format='iso')
    except:
        pdb.set_trace()
    # RENAME
    hstqso_meta.rename_column('GRATE', 'GRATING')
    hstqso_meta.rename_column('QSO_ZEM', 'zem')
    hstqso_meta.rename_column('INST', 'INSTR')
    # ADD
    hstqso_meta.add_column(Column(['HST']*nspec, name='TELESCOPE'))
    return hstqso_meta


def meta_for_build():
    """ Generates the meta data needed for the IGMSpec build
    Returns
    -------
    meta : Table
    """
    # Meta
    hstqso_meta = grab_meta()
    # Cut down to unique sources
    coord = SkyCoord(ra=hstqso_meta['RA'], dec=hstqso_meta['DEC'], unit='deg')
    idx, d2d, d3d = match_coordinates_sky(coord, coord, nthneighbor=2)
    dups = np.where(d2d < 1.5*u.arcsec)[0]  # Closest lens is ~2"
    keep = np.array([True]*len(hstqso_meta))
    for idup in dups:
        dcoord = SkyCoord(ra=hstqso_meta['RA'][idup], dec=hstqso_meta['DEC'][idup], unit='deg')
        sep = dcoord.separation(coord)
        isep = np.where(sep < 1.5*u.arcsec)[0]
        keep[isep] = False
        keep[np.min(isep)] = True  # Only keep 1
    hstqso_meta = hstqso_meta[keep]
    nqso = len(hstqso_meta)
    #
    meta = Table()
    meta['RA'] = hstqso_meta['RA']
    meta['DEC'] = hstqso_meta['DEC']
    meta['zem'] = hstqso_meta['zem']
    meta['sig_zem'] = [0.]*nqso
    meta['flag_zem'] = [str('UNKWN')]*nqso
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
    if sname != 'HSTQSO':
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
    if np.sum(d2d > 1.5*u.arcsec):
        pdb.set_trace()
        raise ValueError("Bad matches in HSTQSO")
    meta_IDs = IDs[idx]

    # Loop me to bid the full survey catalog
    meta.add_column(Column(meta_IDs, name='IGM_ID'))

    # Build spectra (and parse for meta)
    nspec = len(meta)
    max_npix = 61000  # Just needs to be large enough
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
    Rlist = []
    # Loop
    #path = os.getenv('RAW_IGMSPEC')+'/KODIAQ_data_20150421/'
    path = os.getenv('RAW_IGMSPEC')+'/HSTQSO/'
    maxpix = 0
    for jj,row in enumerate(meta):
        # Generate full file
        full_file = path+row['SPEC_FILE']+'.gz'
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
        if 'FOS-L' in fname:
            Rlist.append(300.)
        elif 'FOS-H' in fname:
            Rlist.append(14000.)
        elif 'STIS' in fname:
            if row['GRATING'] == 'G230L':
                Rlist.append(700.)
            elif row['GRATING'] == 'G140L':
                Rlist.append(1200.)
            else:
                raise ValueError("Bad STIS grating")
        elif 'COS' in fname:
            Rlist.append(18000.)  # Approximate, depends on G140L vs G230L
        else:
            raise ValueError("Missing instrument!")
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
    refs = [dict(url='http://adsabs.harvard.edu/abs/2011ApJS..195...16O',
                 bib='omeara11')
            ]
    jrefs = ltu.jsonify(refs)
    hdf[sname]['meta'].attrs['Refs'] = json.dumps(jrefs)
    #
    return

