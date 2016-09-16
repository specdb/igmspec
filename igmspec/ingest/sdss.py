""" Module to ingest SDSS II (aka SDSS) data products
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import os, json
import pdb

import datetime

from astropy.table import Table, Column
from astropy.time import Time
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u

from linetools.spectra import io as lsio
from linetools import utils as ltu

from igmspec.ingest import utils as iiu


def get_specfil(row, dr7=False):
    """Parse the SDSS spectrum file
    Requires a link to the database Class
    """
    if dr7:
        path = os.getenv('RAW_IGMSPEC')+'/SDSS/Schneider/'
    else:
        path = os.getenv('RAW_IGMSPEC')+'/SDSS/spectro_DR7/1d_26/'
    # Generate file name (DR4 is different)
    pnm = '{0:04d}'.format(row['PLATE'])
    fnm = '{0:03d}'.format(row['FIBERID'])
    mjd = str(row['MJD'])
    if dr7:
        sfil = path+'spSpec-'
    else:
        sfil = path+pnm+'/1d/'+'spSpec-'
    # Finish
    specfil = sfil+mjd+'-'+pnm+'-'+fnm+'.fit.gz'  # Is usually gzipped
    return specfil


def grab_meta():
    """ Grab SDSS meta Table

    Returns
    -------
    meta
    """
    old = False
    #sdss_meta = Table.read(os.getenv('RAW_IGMSPEC')+'/SDSS/SDSS_DR7_qso.fits.gz')
    sdss_meta = Table.read(os.getenv('RAW_IGMSPEC')+'/SDSS/dr7qso.fit.gz')
    nspec = len(sdss_meta)
    # DATE
    #t = Time(list(sdss_meta['MJD'].data), format='mjd', out_subfmt='date')  # Fixes to YYYY-MM-DD
    t = Time(list(sdss_meta['SMJD'].data), format='mjd', out_subfmt='date')  # Fixes to YYYY-MM-DD
    sdss_meta.add_column(Column(t.iso, name='DATE-OBS'))
    # Add a few columns
    sdss_meta.add_column(Column([2000.]*nspec, name='EPOCH'))
    sdss_meta.add_column(Column([2000.]*nspec, name='R'))
    sdss_meta.add_column(Column(['SDSS']*nspec, name='INSTR'))
    sdss_meta.add_column(Column(['BOTH']*nspec, name='GRATING'))
    sdss_meta.add_column(Column(['SDSS 2.5-M']*nspec, name='TELESCOPE'))
    # Rename
    if old:
        sdss_meta.rename_column('RAOBJ', 'RA')
        sdss_meta.rename_column('DECOBJ', 'DEC')
        sdss_meta.rename_column('Z_ERR', 'sig_zem')
    else:
        sdss_meta.rename_column('z', 'zem')          # Some of these were corrected by QPQ
        sdss_meta['sig_zem'] = 0.
        sdss_meta['flag_zem'] = 'SDSS-DR7'
    # Sort
    sdss_meta.sort('RA')
    # Return
    return sdss_meta


def meta_for_build(old=False):
    """ Load the meta info
    
    JXP made DR7 -- Should add some aspect of the official list..
      Am worried about the coordinates some..

    Returns
    -------

    """
    sdss_meta = grab_meta()
    # Cut down to unique sources
    coord = SkyCoord(ra=sdss_meta['RA'], dec=sdss_meta['DEC'], unit='deg')
    idx, d2d, d3d = match_coordinates_sky(coord, coord, nthneighbor=2)
    dups = np.where(d2d < 0.5*u.arcsec)[0]
    keep = np.array([True]*len(sdss_meta))
    for idup in dups:
        dcoord = SkyCoord(ra=sdss_meta['RA'][idup], dec=sdss_meta['DEC'][idup], unit='deg')
        sep = dcoord.separation(coord)
        isep = np.where(sep < 0.5*u.arcsec)[0]
        keep[isep] = False
        keep[np.min(isep)] = True  # Only keep 1
    sdss_meta = sdss_meta[keep]
    # Cut one more (pair of QSOs)
    if old:
        bad_dup_c = SkyCoord(ra=193.96678*u.deg, dec=37.099741*u.deg)
        coord = SkyCoord(ra=sdss_meta['RA'], dec=sdss_meta['DEC'], unit='deg')
        sep = bad_dup_c.separation(coord)
        assert np.sum(sep < 2*u.arcsec) == 2
        badi = np.argmin(bad_dup_c.separation(coord))
        keep = np.array([True]*len(sdss_meta))
        keep[badi] = False
        sdss_meta = sdss_meta[keep]
    #
    nqso = len(sdss_meta)
    meta = Table()
    for key in ['RA', 'DEC', 'zem', 'sig_zem']:
        meta[key] = sdss_meta[key]
    meta['flag_zem'] = [str('SDSS-DR7')]*nqso
    meta['STYPE'] = [str('QSO')]*nqso
    # Return
    return meta


def hdf5_adddata(hdf, IDs, sname, debug=False, chk_meta_only=False, sdss_hdf=None, **kwargs):
    """ Add SDSS data to the DB

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
    from igmspec.cat_utils import zem_from_radec
    # Add Survey
    print("Adding {:s} survey to DB".format(sname))
    if sdss_hdf is not None:
        print("Using previously generated {:s} dataset...".format(sname))
        sdss_hdf.copy(sname, hdf)
        return
    sdss_grp = hdf.create_group(sname)
    # Load up
    meta = grab_meta()
    bmeta = meta_for_build()
    # Checks
    if sname != 'SDSS_DR7':
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
    if np.sum(d2d > 1.2*u.arcsec):  # There is one system offset by 1.1"
        raise ValueError("Bad matches in SDSS")
    meta_IDs = IDs[idx]
    meta.add_column(Column(meta_IDs, name='IGM_ID'))

    # Fix zem
    zem, zsource = zem_from_radec(meta['RA'], meta['DEC'], hdf)
    gdz = zem > 0.
    meta['zem'][gdz] = zem[gdz]
    meta['flag_zem'][gdz] = zem[gdz]


    # Build spectra (and parse for meta)
    nspec = len(meta)
    max_npix = 4000  # Just needs to be large enough
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
    maxpix = 0
    for jj,row in enumerate(meta):
        full_file = get_specfil(row)
        if not os.path.isfile(full_file):
            full_file = get_specfil(row, dr7=True)
        # Extract
        print("SDSS: Reading {:s}".format(full_file))
        # Parse name
        fname = full_file.split('/')[-1]
        # Generate full file
        spec = lsio.readspec(full_file)
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
        # Only way to set the dataset correctly
        if chk_meta_only:
            continue
        spec_set[jj] = data

    #
    print("Max pix = {:d}".format(maxpix))
    # Add columns
    meta.add_column(Column(speclist, name='SPEC_FILE'))
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
    refs = [dict(url='http://adsabs.harvard.edu/abs/2010AJ....139.2360S',
                 bib='sdss_qso_dr7'),
            ]
    jrefs = ltu.jsonify(refs)
    hdf[sname]['meta'].attrs['Refs'] = json.dumps(jrefs)
    #
    return
