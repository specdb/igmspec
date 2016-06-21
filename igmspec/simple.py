""" Module for simple igmspec routines
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os, glob
import warnings
import pdb


from igmspec import query_catalog as iqcat
from igmspec import cat_utils as icu
from igmspec import interface_db as igidb

from astropy import units as u

def spec_from_coord(coord, toler=5.*u.arcsec, isurvey=None):
    """ Radial search for spectra around given coordinate
    Best for single searches (i.e. slower than other approaches)

    Parameters
    ----------
    coord : str, tuple, SkyCoord
      See linetools.utils.radec_to_coord
    toler : Angle or Quantity
      Search radius
    isurvey : str or list, optional
      One or more surveys to include
    show_meta : bool, optional
      Show meta data?

    Returns
    -------
    spec : XSpectrum1D
      One or more spectra satisfying the radial search
    meta : Table
      Meta data related to spec

    """
    # Catalog
    qcat = iqcat.QueryCatalog()
    ids = qcat.radial_search(coord, toler)
    if len(ids) == 0:
        warnings.warn("No sources found at your coordinate.  Returning none")
        return None
    elif len(ids) > 1:
        warnings.warn("Found multiple sources.  Hope you expected that.")

    # Surveys
    if isurvey is None:
        surveys = qcat.surveys
    else:
        surveys = qcat.in_surveys(isurvey)

    # Grab all the spectra in the survey(s)
    idb = igidb.InterfaceDB(verbose=False)

    # Load spectra
    spec = idb.grab_spec(surveys, ids)
    return spec, idb.meta
