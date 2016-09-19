""" Module for catalog utilities
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import h5py
import pdb

from astropy import units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord, match_coordinates_sky

from igmspec import defs as idefs


def zem_from_radec(ra, dec, qsos, qtoler=2*u.arcsec):
    """ Parse quasar catalog (Myers) for zem

    Parameters
    ----------
    ra : list or array
      RA in deg
    dec : list or array
      DEC in deg
    qsos : Table
      Must contain RA,DEC,ZEM_SOURCE

    Returns
    -------
    zem : array
      Redshifts
    zsource : array
      str array of sources
    """
    # Generate coordinates
    icoord = SkyCoord(ra=ra, dec=dec, unit='deg')
    # Quasar catalog
    qcoord = SkyCoord(ra=qsos['RA'], dec=qsos['DEC'], unit='deg')
    # Match
    idx, d2d, d3d = match_coordinates_sky(icoord, qcoord, nthneighbor=1)
    good = d2d < qtoler
    pdb.set_trace()
    # Finish
    zem = np.zeros_like(ra.data)
    try:
        zem[good] = qsos['ZEM'][idx[good]]
    except IndexError:
        pdb.set_trace()
    zsource = np.array(['NONENONE']*len(ra))
    zsource[good] = qsos['ZEM_SOURCE'][idx[good]]

    # Return
    return zem, zsource


def flag_to_surveys(flag, survey_dict=None):
    """ Convert flag_survey to list of surveys

    Parameters
    ----------
    flag : int

    Returns
    -------
    surveys : list

    """
    if survey_dict is None:
        survey_dict = idefs.get_survey_dict()
    #
    surveys = []
    for key,sflag in survey_dict.items():
        if flag % (2**sflag) >= sflag:
            surveys.append(key)
    # Return
    return surveys


def write_cat_to_fits(DB_file, cat_fits_file):
    """ Simple script to write the catalog file to a FITS file (mainly for others)
    Parameters
    ----------
    DB_file : str
      Full path to the DB file which contains the catalog
    cat_fits_file : str
      Filename for the FITS file

    Returns
    -------

    """
    if '.fits' not in cat_fits_file:
        raise IOError("Output file {:s} must have .fits extension".format(cat_fits_file))
    # Read
    hdf = h5py.File(DB_file, 'r')
    cat = Table(hdf['catalog'])
    # Write
    cat.write(cat_fits_file)
    # Finish
    hdf.close()
    return
