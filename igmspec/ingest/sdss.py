""" Module to ingest SDSS II (aka SDSS) data products
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import os
import pdb

import datetime

from astropy.table import Table, Column
from astropy.io import fits
from astropy.time import Time

from linetools.spectra import io as lsio

from igmspec.ingest import utils as iiu

def get_specfil(row):
    """Parse the SDSS spectrum file
    Requires a link to the database Class
    """
    path = os.getenv('SDSSPATH')+'/DR7_QSO/spectro/1d_26/'
    # Generate file name (DR4 is different)
    pnm = '{0:04d}'.format(row['PLATE'])
    fnm = '{0:03d}'.format(row['FIBERID'])
    mjd = str(row['MJD'])
    sfil = path+pnm+'/1d/'+'spSpec-'
    # Finish
    specfil = sfil+mjd+'-'+pnm+'-'+fnm+'.fit'  # Is usually gzipped
    return specfil


def grab_meta():
    """ Grab GGG meta Table
    Returns
    -------
    meta
    """
    sdss_meta = Table.read(os.getenv('RAW_IGMSPEC')+'/SDSS/SDSS_DR7_qso.fits.gz')
    nqso = len(sdss_meta)
    # DATE
    dateobs = []
    for row in sdss_meta:
        t = Time(row['MJD'], format='mjd')
        tymd = str(t.iso.split(' ')[0])
        tval = datetime.datetime.strptime(tymd, '%Y-%m-%d')
        dateobs.append(datetime.datetime.strftime(tval,'%Y-%b-%d'))
    sdss_meta.add_column(Column(dateobs, name='DATE-OBS'))
    # Return
    return sdss_meta


def meta_for_build():
    """ Load the meta info
    JXP made DR7 -- Should add some aspect of the official list..
      Am worried about the coordinates some..

    Returns
    -------

    """
    sdss_meta = grab_meta()
    nqso = len(sdss_meta)
    #
    #
    meta = Table()
    meta['RA'] = sdss_meta['RAOBJ']
    meta['DEC'] = sdss_meta['DECOBJ']
    meta['zem'] = sdss_meta['Z'] # Some of these were corrected by QPQ
    meta['sig_zem'] = sdss_meta['Z_ERR']
    meta['flag_zem'] = [str('SDSS')]*nqso  # QPQ too
    # Return
    return meta


def hdf5_adddata(hdf, IDs, sname, debug=False, chk_meta_only=False):
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
       # Add Survey
    print("Adding {:s} survey to DB".format(sname))
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
    meta_IDs = IDs
    meta.add_column(Column(meta_IDs, name='IGM_ID'))

    # Add zem
    meta.rename_column('Z', 'zem')
    meta.rename_column('Z_ERR', 'sig_zem')

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
    Rlist = []
    wvminlist = []
    wvmaxlist = []
    npixlist = []
    speclist = []
    # Loop
    maxpix = 0
    for jj,row in enumerate(meta):
        # Generate full file
        full_file = get_specfil(row)
        # Extract
        print("SDSS: Reading {:s}".format(full_file))
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
        Rlist.append(2000.)
        # Only way to set the dataset correctly
        if chk_meta_only:
            continue
        spec_set[jj] = data

    #
    print("Max pix = {:d}".format(maxpix))
    # Add columns
    meta.add_column(Column([2000.]*nspec, name='EPOCH'))
    meta.add_column(Column(speclist, name='SPEC_FILE'))
    meta.add_column(Column(npixlist, name='NPIX'))
    meta.add_column(Column(wvminlist, name='WV_MIN'))
    meta.add_column(Column(wvmaxlist, name='WV_MAX'))
    meta.add_column(Column(Rlist, name='R'))
    meta.sort('RA')
    meta.add_column(Column(np.arange(nspec,dtype=int),name='SURVEY_ID'))

    # Add HDLLS meta to hdf5
    if iiu.chk_meta(meta):
        if chk_meta_only:
            pdb.set_trace()
        hdf[sname]['meta'] = meta
    else:
        raise ValueError("meta file failed")
    #
    return
