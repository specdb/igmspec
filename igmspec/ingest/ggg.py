""" Module to ingest GGG Survey data

Worseck et al. 2014
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb
import os
import glob
import imp
import datetime

from astropy.table import Table, Column, vstack
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
from astropy.io import fits
from astropy.time import Time

from linetools import utils as ltu
from linetools.spectra import io as lsio

from igmspec.ingest import utils as iiu

igms_path = imp.find_module('igmspec')[1]


def grab_meta():
    """ Grab GGG meta Table
    Returns
    -------

    """
    # This table has units in it!
    ggg_meta = Table.read(os.getenv('RAW_IGMSPEC')+'/GGG/GGG_catalog.fits')
    # RA/DEC, DATE
    dateobs = []
    for row in ggg_meta:
        # DATE
        t = Time(row['RMJD'], format='mjd')
        tymd = str(t.iso.split(' ')[0])
        tval = datetime.datetime.strptime(tymd, '%Y-%m-%d')
        dateobs.append(datetime.datetime.strftime(tval,'%Y-%b-%d'))
    ggg_meta.add_column(Column(dateobs, name='DATE-OBS'))
    # Turn off units
    for key in ['RA', 'DEC']:
        ggg_meta[key].unit = None
    #
    return ggg_meta


def meta_for_build():
    """ Generates the meta data needed for the IGMSpec build
    Returns
    -------
    meta : Table
    """
    ggg_meta = grab_meta()
    # Cut down to unique QSOs
    names = np.array([name[0:26] for name in ggg_meta['SDSSJ']])
    uni, uni_idx = np.unique(names, return_index=True)
    ggg_meta = ggg_meta[uni_idx]
    nqso = len(ggg_meta)
    #
    meta = Table()
    meta['RA'] = ggg_meta['RA']
    meta['DEC'] = ggg_meta['DEC']
    meta['zem'] = ggg_meta['z_gmos']
    meta['sig_zem'] = ggg_meta['zerror_gmos']
    meta['flag_zem'] = [str('GGG')]*nqso
    # Return
    return meta


def hdf5_adddata(hdf, IDs, sname, debug=False, chk_meta_only=False):
    """ Append GGG data to the h5 file

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
    ggg_grp = hdf.create_group(sname)
    # Load up
    meta = grab_meta()
    bmeta = meta_for_build()
    # Checks
    if sname != 'GGG':
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
    meta['zem'] = meta['z_gmos']

    # Double up for the two gratings
    meta = vstack([meta,meta])

    # Build spectra (and parse for meta)
    nspec = len(meta)
    max_npix = 1600  # Just needs to be large enough
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
    path = os.getenv('RAW_IGMSPEC')+'/GGG/'
    maxpix = 0
    for jj,row in enumerate(meta):
        # Generate full file
        if jj >= nspec//2:
            full_file = path+row['name']+'_R400.fits.gz'
        else:
            full_file = path+row['name']+'_B600.fits.gz'
        # Extract
        print("GGG: Reading {:s}".format(full_file))
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
        if 'R400' in fname:
            Rlist.append(833.)
        else:
            Rlist.append(940.)
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


