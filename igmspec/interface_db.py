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
    survey_bool : bool array
      Used to grab data from a given survey
    """

    def __init__(self, db_file=None, maximum_ram=10.,verbose=True):
        """
        Returns
        -------

        """
        # Init
        self.memory_used = 0.
        self.memory_warning = 5.  # Gb
        self.memory_max = 10.  # Gb
        self.hdf = None
        self.verbose = verbose
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
        if self.verbose:
            print("Using {:s} for the DB file".format(db_file))
        self.hdf = h5py.File(db_file,'r')
        self.db_file = db_file
        self.survey_IDs = None
        #
        surveys = self.hdf.keys()
        surveys.pop(surveys.index('catalog'))
        self.surveys = surveys
        if self.verbose:
            print("Available surveys: {}".format(self.surveys))

    def grab_ids(self, survey, IGM_IDs, meta=None, match_meta=None):
        """ Grab the rows in a survey matching IGM_IDs

        Parameters
        ----------
        survey : str
          Name of the Survey
        IGM_IDs : int or ndarray
          IGM_ID values
        meta : Table, optional
          Meta data for the survey (usually read from hdf)
        match_meta : dict, optional
          key/value to filter spec with, e.g. {'INSTR': HIRES}

        Returns
        ------
        match_survey : ndarray
          bool array of meta rows from the survey matching IGM_IDs

        """
        # Check
        if survey not in self.hdf.keys():
            raise IOError("Survey {:s} not in your DB file {:s}".format(survey, self.db_file))
        # Find IGMS_IDs indices in survey
        if meta is None:
            meta = self.hdf[survey]['meta'].value
        match_survey = np.in1d(meta['IGM_ID'], IGM_IDs)
        # Match meta?
        if match_meta is not None:
            for key,value in match_meta.items():
                match_survey = match_survey & (meta[key] == value)
        # Store and return
        self.survey_bool = match_survey
        self.meta = meta[match_survey]
        return match_survey

    def grab_meta(self, survey, IGM_IDs, show=True):
        """ Grab meta data for survey
        Parameters
        ----------
        survey : str or list
        IGM_IDs : int or array
        show : bool, optional
          Show the Meta table (print) in addition to returning

        Returns
        -------
        meta : Table(s)

        """
        #
        if isinstance(survey, list):
            all_meta = []
            for isurvey in survey:
                all_meta.append(self.grab_meta(isurvey, IGM_IDs))
            return all_meta
        # Grab IDs then cut
        _ = self.grab_ids(survey, IGM_IDs)
        meta = Table(self.hdf[survey]['meta'].value)[self.survey_bool]
        self.meta = meta
        self.meta['RA'].format = '7.3f'
        self.meta['DEC'].format = '7.3f'
        self.meta['zem'].format = '6.3f'
        self.meta['WV_MIN'].format = '6.1f'
        self.meta['WV_MAX'].format = '6.1f'
        # Show?
        if show:
            meta_keys = ['IGM_ID', 'RA', 'DEC', 'zem', 'SPEC_FILE']
            for key in meta.keys():
                if key not in meta_keys:
                    meta_keys += [key]
            meta[meta_keys].pprint(max_width=120)
        # Load and return
        return meta

    def grab_spec(self, survey, IGM_IDs, verbose=None, **kwargs):
        """ Grab spectra using staged IDs

        Parameters
        ----------
        survey : str or list
        IGM_IDs : int or intarr

        Returns
        -------

        """
        if verbose is None:
            verbose = self.verbose
        if isinstance(survey, list):
            all_spec = []
            for isurvey in survey:
                all_spec.append(self.grab_spec(isurvey, IGM_IDs, **kwargs))
            return all_spec
        # Grab IDs
        if self.stage_data(survey, IGM_IDs, **kwargs):
            if verbose:
                print("Loaded spectra")
            data = self.hdf[survey]['spec'][self.survey_bool]
        else:
            print("Staging failed..  Not returning spectra")
            return
        # Generate XSpectrum1D
        spec = XSpectrum1D(data['wave'], data['flux'], sig=data['sig'], masking='edges')
        # Return
        return spec

    def stage_data(self, survey, IGM_IDs, verbose=None, **kwargs):
        """ Stage the spectra for serving
        Mainly checks the memory

        Parameters
        ----------
        survey : str
          Name of the Survey
        IGM_IDs : int or ndarray
          IGM_ID values

        Returns
        -------
        check : bool
        If True, self.survey_IDs is filled
          Indices in the survey dataset

        """
        if verbose is None:
            verbose = self.verbose
        # Checks
        if survey not in self.hdf.keys():
            raise IOError("Survey {:s} not in your DB file {:s}".format(survey, self.db_file))
        match_survey = self.grab_ids(survey, IGM_IDs, **kwargs)
        # Meta?
        nhits = np.sum(match_survey)
        # Memory check (approximate; ignores meta data)
        spec_Gb = self.hdf[survey]['spec'][0].nbytes/1e9  # Gb
        new_memory = spec_Gb*nhits
        if new_memory + self.memory_used > self.memory_max:
            warnings.warn("This request would exceed your maximum memory limit of {:g} Gb".format(self.memory_max))
            return False
        else:
            if verbose:
                print("Staged {:d} spectra totalling {:g} Gb".format(nhits, new_memory))
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
