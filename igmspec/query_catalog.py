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

from igmspec import defs as idefs
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
        return self.cat['IGMsp_ID'][good]

    def radial_search(self, inp, radius, verbose=True):
        """ Search for sources in a radius around the input coord

        Parameters
        ----------
        inp : str or tuple or SkyCoord
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
        return self.cat['IGMsp_ID'][good]

    def show_meta(self, good):
        """
        Parameters
        ----------
        good : bool array
          True means show

        Returns
        -------

        """
        # IGMspec catalog
        # Catalog keys
        cat_keys = ['IGMsp_ID', 'RA', 'DEC', 'zem', 'flag_survey']
        for key in self.cat.keys():
            if key not in cat_keys:
                cat_keys += [key]
        self.cat[cat_keys][good].pprint(max_width=120)
        # Print survey dict
        print("----------")
        print("Survey key:")
        for survey in self.surveys:
            print("    {:s}: {:d}".format(survey, idefs.get_survey_dict()[survey]))


    def setup(self):
        """ Set up a few things, e.g. SkyCoord for the catalog
        Returns
        -------

        """
        from igmspec import cat_utils as icu
        # SkyCoord
        self.coords = SkyCoord(ra=self.cat['RA'], dec=self.cat['DEC'], unit='deg')
        # Formatting the Table
        self.cat['RA'].format = '7.3f'
        self.cat['DEC'].format = '7.3f'
        self.cat['zem'].format = '6.3f'
        self.cat['sig_zem'].format = '5.3f'
        # Surveys
        surveys = idefs.get_survey_dict()
        unif = np.unique(self.cat['flag_survey'])
        all_surveys = []
        for ifs in unif:
            all_surveys += icu.flag_to_surveys(ifs)
        self.surveys = list(np.unique(all_surveys))

    def __repr__(self):
        txt = '<{:s}:  DB_file={:s} with {:d} sources\n'.format(self.__class__.__name__,
                                            self.db_file, len(self.cat))
        # Surveys
        txt += '   Loaded surveys are {} \n'.format(self.surveys)
        txt += '>'
        return (txt)
