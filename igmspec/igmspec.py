""" Module for IgmSpec Class
"""
from __future__ import print_function, absolute_import, division, unicode_literals

#import numpy as np
import pdb
import warnings

from astropy import units as u

from igmspec.query_catalog import QueryCatalog
from igmspec.interface_db import InterfaceDB


class IgmSpec(object):
    """ An uber-class for the package
    Ideally one-stop-shopping for most consumers

    Parameters
    ----------

    Attributes
    ----------
    qcat : QueryCatalog
    idb : InterfaceDB
    """

    def __init__(self, **kwargs):
        """
        """
        # Init
        self.qcat = QueryCatalog(**kwargs)
        self.idb = InterfaceDB(**kwargs)
        # Checks
        assert self.idb.db_file == self.qcat.db_file
        for survey in self.idb.surveys:
            assert survey in self.qcat.surveys

    def spec_from_coord(self, coord, tol=5.*u.arcsec, isurvey=None, **kwargs):
        """ Radial search for spectra around given coordinate
        Best for single searches (i.e. slower than other approaches)

        Parameters
        ----------
        coord : str, tuple, SkyCoord
          See linetools.utils.radec_to_coord
        tol : Angle or Quantity, optional
          Search radius
        isurvey : str or list, optional
          One or more surveys to include
        kwargs :
          fed to grab_spec


        Returns
        -------
        spec : XSpectrum1D or list of XSpectrum1D
          One or more spectra satisfying the radial search
        meta : Table or list of Tables
          Meta data related to spec

        """
        # Catalog
        ids = self.qcat.radial_search(coord, tol)
        if len(ids) == 0:
            warnings.warn("No sources found at your coordinate.  Returning none")
            return None
        elif len(ids) > 1:
            warnings.warn("Found multiple sources.  Hope you expected that.")

        # Surveys
        if isurvey is None:
            surveys = self.qcat.surveys
        else:
            surveys = self.qcat.in_surveys(isurvey)

        # Load spectra
        spec, meta = self.idb.grab_spec(surveys, ids, **kwargs)
        return spec, meta

    def __getattr__(self, k):
        """ Overload attributes using the underlying classes

        Parameters
        ----------
        k

        Returns
        -------

        """
        # Try DB first
        try:
            return getattr(self.idb, k)
        except AttributeError:
            # Try qcat last
            return getattr(self.qcat, k)

    def __repr__(self):
        txt = '<{:s}:  IGM_file={:s} with {:d} sources\n'.format(self.__class__.__name__,
                                            self.db_file, len(self.cat))
        # Surveys
        txt += '   Loaded surveys are {} \n'.format(self.surveys)
        txt += '>'
        return (txt)
