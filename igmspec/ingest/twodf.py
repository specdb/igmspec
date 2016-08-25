""" Module to ingest 2dF/6dF quasars
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import os, json
import pdb

import datetime

from astropy.table import Table, Column
from astropy.io import fits
from astropy.time import Time

from linetools.spectra import io as lsio
from linetools.spectra.xspectrum1d import XSpectrum1D
from linetools import utils as ltu

from igmspec.ingest import utils as iiu

def get_specfil(row):
    """Parse the 2QZ spectrum file
    Requires a link to the database Class
    """
    path = os.getenv('RAW_IGMSPEC')+'/2QZ/2df/fits/'
    # RA/DEC folder
    path += 'ra{:02d}_{:02d}/'.format(row['RAh00'], row['RAh00']+1)
    # File name
    sfil = path+row['Name']
    if row['ispec'] == 1:
        sfil += 'a'
    elif row['ispec'] == 2:
        sfil += 'b'
    else:
        raise ValueError("Bad ispec value")
    # Finish
    specfil = sfil+'.fits.gz'
    return specfil


def grab_meta():
    """ Grab GGG meta Table
    Catalog -- http://www.2dfquasar.org/Spec_Cat/catalogue.html

    Returns
    -------
    meta
    """
    catfil = os.getenv('RAW_IGMSPEC')+'/2QZ/2QZ393524355693.out'
    tdf_meta = Table.read(catfil, format='ascii')
    # Rename columns
    clms = ['Name', 'RAh00', 'RAm00', 'RAs00', 'DECd00', 'DECm00', 'DECs00',
       'ID','cat_name', 'Sector', 'RAh50', 'RAm50', 'RAs50', 'DECd50', 'DECm50', 'DECs50',
       'UKST', 'XAPM','YAPM','RA50','DEC50','bj','u-b','b-r','Nobs',
       'z1','q1','ID1','date1','fld1','fiber1','SN1',
       'z2','q2','ID2','date2','fld2','fiber2','SN2',
        'zprev','rflux','Xray','EBV','comm1','comm2']
    for ii in range(1,46):
        tdf_meta.rename_column('col{:d}'.format(ii), clms[ii-1])
    # Cut down to QSOs and take only 1 spectrum
    ispec = []
    zspec = []
    datespec = []
    for row in tdf_meta:
        if 'QSO' in row['ID1']:
            ispec.append(1)
            zspec.append(row['z1'])
            sdate = str(row['date1'])
            datespec.append('{:s}-{:s}-{:s}'.format(sdate[0:4],sdate[4:6], sdate[6:8]))
        elif 'QSO' in row['ID2']:
            ispec.append(2)
            zspec.append(row['z2'])
            sdate = str(row['date2'])
            datespec.append('{:s}-{:s}-{:s}'.format(sdate[0:4],sdate[4:6], sdate[6:8]))
        else:
            ispec.append(0)
            zspec.append(-1.)
            datespec.append('')
    tdf_meta['ispec'] = ispec
    tdf_meta['zem'] = zspec
    tdf_meta['DATE'] = datespec
    cut = tdf_meta['ispec'] > 0
    tdf_meta = tdf_meta[cut]
    nspec = len(tdf_meta)
    # DATE
    t = Time(list(tdf_meta['DATE'].data), format='iso', out_subfmt='date')  # Fixes to YYYY-MM-DD
    tdf_meta.add_column(Column(t.iso, name='DATE-OBS'))
    # Add a few columns
    tdf_meta.add_column(Column([2000.]*nspec, name='EPOCH'))
    # Resolution
    #  2df 8.6A FWHM
    #  R=580 at 5000A
    #  6df??
    tdf_meta.add_column(Column([580.]*nspec, name='R'))
    #
    tdf_meta.add_column(Column(['2dF']*nspec, name='INSTR'))
    tdf_meta.add_column(Column(['300B']*nspec, name='GRATING'))
    tdf_meta.add_column(Column(['UKST']*nspec, name='TELESCOPE'))
    # Rename
    rad = (tdf_meta['RAh00']*3600 + tdf_meta['RAm00']*60 + tdf_meta['RAs00'])*360./86400.
    decd = tdf_meta['DECd00'] + tdf_meta['DECm00']/60 + tdf_meta['DECs00']/3600.
    tdf_meta['RA'] = rad
    tdf_meta['DEC'] = decd
    tdf_meta['sig_zem'] = [0.]*nspec
    # Sort
    tdf_meta.sort('RA')
    # Return
    return tdf_meta


def meta_for_build():
    """ Load the meta info
    JXP made DR7 -- Should add some aspect of the official list..
      Am worried about the coordinates some..

    Returns
    -------

    """
    tdf_meta = grab_meta()
    nqso = len(tdf_meta)
    #
    #
    meta = Table()
    for key in ['RA', 'DEC', 'zem', 'sig_zem']:
        meta[key] = tdf_meta[key]
    meta['flag_zem'] = [str('2QZ')]*nqso  # QPQ too
    meta['STYPE'] = [str('QSO')]*nqso  # QPQ too
    # Return
    return meta


def hdf5_adddata(hdf, IDs, sname, debug=False, chk_meta_only=False):
    """ Add 2QZ data to the DB

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
    tqz_grp = hdf.create_group(sname)
    # Load up
    meta = grab_meta()
    bmeta = meta_for_build()
    # Checks
    if sname != '2QZ':
        raise IOError("Not expecting this survey..")
    if np.sum(IDs < 0) > 0:
        raise ValueError("Bad ID values")
    # Open Meta tables
    if len(bmeta) != len(IDs):
        raise ValueError("Wrong sized table..")

    # Generate ID array from RA/DEC
    meta_IDs = IDs
    meta.add_column(Column(meta_IDs, name='IGM_ID'))

    # Add zem

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
        # Extract
        print("2QZ: Reading {:s}".format(full_file))
        # Parse name
        fname = full_file.split('/')[-1]
        if debug:
            if jj > 500:
                speclist.append(str(fname))
                if not os.path.isfile(full_file):
                    raise IOError("2QZ file {:s} does not exist".format(full_file))
                wvminlist.append(np.min(data['wave'][0][:npix]))
                wvmaxlist.append(np.max(data['wave'][0][:npix]))
                npixlist.append(npix)
                continue
        # Read
        hdu = fits.open(full_file)
        head0 = hdu[0].header
        wave = lsio.setwave(head0)
        flux = hdu[0].data
        var = hdu[2].data
        sig = np.zeros_like(flux)
        gd = var > 0.
        sig[gd] = np.sqrt(var)
        # npix
        spec = XSpectrum1D.from_tuple((wave,flux,sig))
        npix = spec.npix
        spec.header = head0
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
    #
    # References
    refs = [dict(url='http://adsabs.harvard.edu/abs/2004MNRAS.349.1397C',
                 bib='2QZ')
            ]
    jrefs = ltu.jsonify(refs)
    hdf[sname]['meta'].attrs['Refs'] = json.dumps(jrefs)
    #
    return
