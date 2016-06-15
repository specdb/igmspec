""" Module to ingest HD-LLS Survey data

Prochaska et al. 2015
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os
from astropy.table import Table


def meta_for_build():
    """ Generates the meta data needed for the IGMSpec build
    Returns
    -------
    meta : Table
    """
    hdlls_meta = Table.read(os.getenv('RAW_IGMSPEC')+'/HD-LLS_DR1.fits.gz')
    nqso = len(hdlls_meta)
    #
    meta = Table()
    meta['RA'] = hdlls_meta['RA']
    meta['DEC'] = hdlls_meta['DEC']
    meta['EPOCH'] = [2000.]*nqso
    meta['zem'] = hdlls_meta['Z_QSO']
    meta['flag_z'] = ['UNKN']*nqso
    # Return
    return meta

def append(h5):
    """ Append HD-LLS to the h5 file

    Parameters
    ----------
    h5 : hdf5 pointer

    Returns
    -------

    """
    # Open Meta table



