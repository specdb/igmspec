""" Module for IgmSpec Class
"""
from __future__ import print_function, absolute_import, division, unicode_literals

#import numpy as np
import pdb
import warnings

from astropy import units as u

from specdb.specdb import SpecDB


class IgmSpec(SpecDB):
    """ Main class for using IGMSpec DB
    See SpecDB for full docs

    Parameters
    ----------
    skip_test : bool, optional
      Skip tests?  Highly *not* recommended

    Attributes
    ----------
    """

    def __init__(self, db_file=None, skip_test=False, **kwargs):
        """
        """
        # db_file
        SpecDB.__init__(self, db_file=db_file, skip_test=skip_test, **kwargs)

    def grab_dbfile(self):
        """ Grabs the DB file
        Returns
        -------
        db_file : str
          full path to the DB file

        """
        import os, glob
        if os.getenv('IGMSPEC_DB') is None:
            warnings.warn('Environmental variable IGMSPEC_DB not set.  Assuming this is a test')
            import igmspec
            db_dir = igmspec.__path__[0]+'/tests/files/'
        else:
            db_dir = os.getenv('IGMSPEC_DB')
        #
        fils = glob.glob(db_dir+'/IGMspec_DB_*hdf5')
        fils.sort()
        db_file = fils[-1]  # Should grab the latest
        # Return
        return db_file

    def __repr__(self):
        txt = '<{:s}:  IGM_file={:s} with {:d} sources\n'.format(self.__class__.__name__,
                                            self.db_file, len(self.cat))
        # Surveys
        txt += '   Loaded surveys are {} \n'.format(self.surveys)
        txt += '>'
        return (txt)
