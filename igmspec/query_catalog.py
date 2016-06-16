""" Module to interface with hdf5 database for IGMspec
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import os, glob, imp
import psutil
import warnings
import h5py
import numpy as np
import pdb


from astropy.table import Table
from astropy import units as u
from astropy.units import Quantity
from astropy.coordinates import SkyCoord, match_coordinates_sky, Angle

from igmspec.defs import z_priority, survey_flag
from igmspec import db_utils as idbu

from linetools import utils as ltu

class QueryCatalog(object):
    """ A Class for querying the IGMspec catalog

    Parameters
    ----------

    Attributes
    ----------
    cat : Table
      Astropy Table holding the IGMspec catalog
    """

    def __init__(self, db_file=None, maximum_ram=10.):
        """
        Returns
        -------

        """
        # Init
        # Load catalog
        self.load_cat(db_file)
        # Setup
        self.setup()

    def load_cat(self, db_file=None):
        """ Open the DB file
        Parameters
        ----------
        db_file

        Returns
        -------

        """
        #
        if db_file is None:
            db_file = idbu.grab_dbfile()
        print("Using {:s} for the catalog file".format(db_file))
        hdf = h5py.File(db_file,'r')
        self.cat = Table(hdf['catalog'].value)
        self.db_file = db_file

    def match_coord(self, cat_coords, toler=0.5*u.arcsec, verbose=True):
        """ Match an input set of SkyCoords to the catalog within a given radius

        Parameters
        ----------
        coords : SkyCoord
          Single or array
        toler : Angle or Quantity, optional
          Tolerance for a match
        verbose : bool, optional

        Returns
        -------
        indices : bool array
          True = match

        """
        # Checks
        if not isinstance(toler, (Angle, Quantity)):
            raise IOError("Input radius must be an Angle type, e.g. 10.*u.arcsec")
        # Match
        idx, d2d, d3d = match_coordinates_sky(self.coords, cat_coords, nthneighbor=1)
        good = d2d < toler
        # Return
        if verbose:
            print("Your search yielded {:d} matches".format(np.sum(good)))
        return good

    def radial_search(self, inp, radius, verbose=True):
        """ Search for sources in a radius around the input coord

        Parameters
        ----------
        inp : tuple or SkyCoord
          See linetools.utils.radec_to_coord
        toler
        verbose

        Returns
        -------

        """
        # Convert to SkyCoord
        coord = ltu.radec_to_coord(inp)
        # Separation
        sep = coord.separation(self.coords)
        # Match
        good = sep < radius
        # Return
        if verbose:
            print("Your search yielded {:d} match[es]".format(np.sum(good)))
        return good

    def setup(self):
        """ Set up a few things, e.g. SkyCoord for the catalog
        Returns
        -------

        """
        self.coords = SkyCoord(ra=self.cat['RA'], dec=self.cat['DEC'], unit='deg')


    def __repr__(self):
        txt = '<{:s}:  DB_file={:s} with {:d} sources>'.format(self.__class__.__name__,
                                            self.db_file, len(self.cat))
        return (txt)
