""" Module to ingest Myers' QSOs
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os
import pdb

from astropy.table import Table
from astropy.io import fits


def add_to_hdf(hdf):
    """ Add Meyers catalog to hdf file

    Parameters
    ----------
    hdf : HDF5 file
    """
    # Load
    ADM_qso, date = load()
    # Add
    hdf['quasars'] = ADM_qso
    hdf['quasars'].attrs['DATE'] = date
    #
    return


def load():
    """ Load catalog

    Parameters
    ----------

    Returns
    -------
    cat : Table
    date : str
      DATE of creation

    """
    ADM_file = os.getenv('RAW_IGMSPEC')+'/Myers/GTR-ADM-QSO-master-wvcv.fits.gz'
    ADM_qso = Table.read(ADM_file)
    # Grab header for DATE
    head1 = fits.open(ADM_file)[1].header
    # Return
    return ADM_qso, head1['DATE']


