""" Module to interface with hdf5 database for IGMspec
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os, glob, imp
import psutil
import warnings
import h5py
import pdb

import numpy as np

from astropy.table import Table, vstack, Column
from astropy import units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky

from linetools.spectra.xspectrum1d import XSpectrum1D

from igmspec.defs import z_priority, survey_flag
from igmspec.db_utils import grab_dbfile

class InterfaceDB(object):
    """ A Class for interfacing with the DB

    Parameters
    ----------

    Attributes
    ----------
    memory_used : float
      Used memory in Gb
    memory_warning : float
      Value at which a Warning is raised
    hdf : pointer to DB
    maximum_ram : float, optonal
      Maximum memory allowed for the Python session, in Gb
    """

    def __init__(self, db_file=None, maximum_ram=10.):
        """
        Returns
        -------

        """
        # Init
        self.memory_used = 0.
        self.memory_warning = 5.  # Gb
        self.memory_max = 10.  # Gb
        self.hdf = None
        # Open DB
        self.open_db(db_file)
        # Check memory
        self.update()

    def open_db(self, db_file=None):
        """ Open the DB file
        Parameters
        ----------
        db_file

        Returns
        -------

        """
        if db_file is None:
            db_file = grab_dbfile()
        #
        print("Using {:s} for the DB file".format(db_file))
        self.hdf = h5py.File(db_file,'r')
        self.db_file = db_file
        self.survey_IDs = None
        #
        surveys = self.hdf.keys()
        surveys.pop(surveys.index('catalog'))
        self.surveys = surveys
        print("Available surveys: {}".format(self.surveys))

    def grab_meta(self, survey):
        """ Grab meta data for survey
        Returns
        -------

        """

    def grab_spec(self, survey, IGMsp_IDs):
        """ Grab spectra using staged IDs

        Parameters
        ----------
        survey

        Returns
        -------

        """
        if self.stage_data(survey, IGMsp_IDs):
            data = self.hdf[survey]['spec'][self.survey_bool]
        else:
            print("Staging failed..  Not returning spectra")
            return
        # Deal with padding?
        tmp = np.sum(data['sig'],axis=0)
        notzero = np.where(tmp != 0.)[0]
        npix = np.max(notzero)
        # Generate XSpectrum1D
        spec = XSpectrum1D(data['wave'][:,:npix], data['flux'][:,:npix],
                           sig=data['sig'][:,:npix])
        # Return
        return spec

    def stage_data(self, survey, IGMsp_IDs):
        """ Stage the spectra for serving
        Mainly checks the memory

        Parameters
        ----------
        survey : str
          Name of the Survey
        IGMS_IDs : int or ndarray
          IGMS_ID values

        Returns
        -------
        check : bool
        If True, self.survey_IDs is filled
          Indices in the survey dataset

        """
        # Checks
        if survey not in self.hdf.keys():
            raise IOError("Survey {:s} not in your DB file {:s}".format(self.db_file))
        # Find IGMS_IDs indices in survey
        meta = self.hdf[survey]['meta'].value
        if 'IGMsp_ID' not in meta.dtype.names:
            raise ValueError("Meta table in {:s} survey is missing IGSM_ID column!".format(survey))
        in_survey = np.in1d(meta['IGMsp_ID'], IGMsp_IDs)
        nhits = np.sum(in_survey)
        # Memory check (approximate; ignores meta data)
        spec_bytes = self.hdf[survey]['spec'][0].nbytes/1e9  # Gb
        new_memory = spec_bytes*nhits
        if new_memory + self.memory_used > self.memory_max:
            warnings.warn("This request would exceed your maximum memory limit of {:g} Gb".format(self.memory_max))
            return False
        else:
            self.survey_bool = in_survey
            return True


    def update(self):
        """ Update key attributes

        Returns
        -------

        """
        # Memory
        process = psutil.Process(os.getpid())
        self.memory_used = process.memory_info().rss/1e9  # Gb
        if self.memory_used > self.memory_warning:
            warnings.warn("Your memory usage -- {:g} Gb -- is high".format(self.memory_used))

    def __repr__(self):
        txt = '<{:s}:  DB_file={:s} \n'.format(self.__class__.__name__,
                                            self.db_file)
        # Surveys
        txt += '   Loaded surveys are {} \n'.format(self.surveys)
        txt += '>'
        return (txt)
