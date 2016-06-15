""" Module to build the hdf5 database file for IGMspec
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import igmspec
import h5py
import pdb

from igmspec.defs import z_priority, survey_flag
from igmspec.ingest import boss

from astropy.table import Table, vstack, Column
from astropy import units as u

def chk_maindb_join(maindb, newdb):
    """ Check that new data is consistent with existing table

    Parameters
    ----------
    maindb : Table
    newdb : Table

    Returns
    -------
    answer : bool

    """
    for key in newdb.keys():
        try:
            assert key in maindb.keys()
        except AssertionError:
            pdb.set_trace()
            return False
    return True

def ver01():
    """ Build version 0.1
    Returns
    -------

    """
    outfil = igmspec.__path__[0]+'/../IGMspec_DB.hdf5'

    # Defs
    zpri = z_priority()
    lenz = [len(zpi) for zpi in zpri]
    dummyf = '#'*np.max(np.array(lenz))  # For the Table

    # Main DB Table  (WARNING: THIS MAY TURN INTO SQL)
    idict = dict(RA=0., DEC=0., IGMS_ID=0, EPOCH=2000.,
                 zem=0., sig_zem=0., flag_zem=dummyf, flag_survey=0)
    tkeys = idict.keys()
    lst = [[idict[tkey]] for tkey in tkeys]
    maindb = Table(lst, names=tkeys)

    # BOSS-DR12
    # Read
    boss_meta = boss.meta_for_build()
    nboss = len(boss_meta)
    # IDs
    boss_meta.add_column(Column(np.arange(nboss,dtype=int), name='IGMS_ID'))
    flag_s = survey_flag('BOSS-DR12')
    boss_meta.add_column(Column([flag_s]*nboss, name='flag_survey'))
    # Check
    assert chk_maindb_join(maindb, boss_meta)
    # Append
    maindb = vstack([maindb,boss_meta])
    maindb = maindb[1:]  # Eliminate dummy line
    pdb.set_trace()

if __name__ == '__main__':
    ver01()

