""" Module to ingest SDSS III (aka BOSS) data products
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os
from astropy.table import Table
from astropy.io import fits


def meta_for_build():
    """ Load the meta info
    DR12 quasars : https://data.sdss.org/datamodel/files/BOSS_QSO/DR12Q/DR12Q.html

    Returns
    -------

    """
    boss_meta = Table.read(os.getenv('RAW_IGMSPEC')+'/DR12Q.fits')
    nqso = len(boss_meta)
    #
    meta = Table()
    meta['RA'] = boss_meta['RA']
    meta['DEC'] = boss_meta['DEC']
    meta['EPOCH'] = [2000.]*nqso
    meta['zem'] = boss_meta['Z_PCA']
    meta['sig_zem'] = boss_meta['ERR_ZPCA']
    meta['flag_zem'] = ['BOSS_PCA']*nqso
    # Return
    return meta


def qsos_dr12():
    """
    for the data model
    Returns
    -------

    """
