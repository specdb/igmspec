""" Module to ingest SDSS III (aka BOSS) data products
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import os
import pdb
import datetime

from astropy.table import Table, Column
from astropy.time import Time

def grab_meta():
    """ Grab BOSS meta Table
    Returns
    -------

    """
    boss_meta = Table.read(os.getenv('RAW_IGMSPEC')+'/BOSS/DR12Q.fits')
    nboss = len(boss_meta)
    # DATE-OBS
    t = Time(list(boss_meta['MJD'].data), format='mjd', out_subfmt='date')  # Fixes to YYYY-MM-DD
    boss_meta.add_column(Column(t.iso, name='DATE-OBS'))
    # Add columns
    boss_meta.add_column(Column(['BOSS']*nboss, name='INSTR'))
    boss_meta.add_column(Column(['BOTH']*nboss, name='GRATING'))
    boss_meta.add_column(Column([2000.]*nboss, name='R'))  # RESOLUTION
    boss_meta.add_column(Column(['SDSS 2.5-M']*nboss, name='TELESCOPE'))
    # Redshift logic
    boss_meta['zem'] = boss_meta['Z_PCA']
    boss_meta['sig_zem'] = boss_meta['ERR_ZPCA']
    boss_meta['flag_zem'] = [str('BOSS_PCA ')]*nboss
    # Fix bad redshifts
    bad_pca = boss_meta['Z_PCA'] < 0.
    boss_meta['zem'][bad_pca] = boss_meta['Z_PIPE'][bad_pca]
    boss_meta['sig_zem'][bad_pca] = boss_meta['ERR_ZPIPE'][bad_pca]
    boss_meta['flag_zem'][bad_pca] = str('BOSS_PIPE')
    #
    return boss_meta

def meta_for_build():
    """ Load the meta info
    DR12 quasars : https://data.sdss.org/datamodel/files/BOSS_QSO/DR12Q/DR12Q.html

    Returns
    -------

    """
    boss_meta = grab_meta()
    #
    meta = Table()
    for key in ['RA', 'DEC', 'zem', 'sig_zem', 'flag_zem']:
        meta[key] = boss_meta[key]
    # Return
    return meta


def qsos_dr12():
    """
    for the data model
    Returns
    -------

    """
