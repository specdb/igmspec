""" Module for DB utilities
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os, glob
import pdb

def grab_dbfile():
    """ Grabs the DB file
    Returns
    -------
    db_file : str
      full path to the DB file

    """
    if os.getenv('IGMSPEC_DB') is None:
        raise IOError("You need to set the environmental variable IGMSPEC_DB")
    fils = glob.glob(os.getenv('IGMSPEC_DB')+'/IGMspec_DB_*hdf5')
    fils.sort()
    db_file = fils[-1] # Should grab the lateset
    # Return
    return db_file
