""" Module to ingest SDSS III (aka BOSS) data products
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import os
import pdb

from astropy.table import Table
from astropy.io import fits


def meta_for_build():
    """ Load the meta info
    DR12 quasars : https://data.sdss.org/datamodel/files/BOSS_QSO/DR12Q/DR12Q.html

    Returns
    -------

    """
    boss_meta = Table.read(os.getenv('RAW_IGMSPEC')+'/BOSS/DR12Q.fits')
    nqso = len(boss_meta)
    #
    #
    meta = Table()
    meta['RA'] = boss_meta['RA']
    meta['DEC'] = boss_meta['DEC']
    meta['zem'] = boss_meta['Z_PCA']
    meta['sig_zem'] = boss_meta['ERR_ZPCA']
    meta['flag_zem'] = [str('BOSS_PCA ')]*nqso
    # Fix bad redshifts
    bad_pca = boss_meta['Z_PCA'] < 0.
    meta['zem'][bad_pca] = boss_meta['Z_PIPE'][bad_pca]
    meta['sig_zem'][bad_pca] = boss_meta['ERR_ZPIPE'][bad_pca]
    meta['flag_zem'][bad_pca] = str('BOSS_PIPE')
    meta['flag_zem'] = str('BOSS_PCA')
    # Return
    return meta


def qsos_dr12():
    """
    for the data model
    Returns
    -------

    """
