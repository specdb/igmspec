""" Module to build the hdf5 database file for IGMspec
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import igmspec
import h5py
import numbers
import pdb

from igmspec.defs import z_priority, survey_flag
from igmspec.ingest import boss, hdlls

from astropy.table import Table, vstack, Column
from astropy import units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky


def add_to_flag(cur_flag, add_flag):
    """ Add a bitwise flag to an existing flat
    Parameters
    ----------
    cur_flag : int or ndarray
    add_flag : int

    Returns
    -------
    new_flag : int or ndarray

    """
    if isinstance(cur_flag, numbers.Integral):
        if (cur_flag % add_flag) < (add_flag//2):
            cur_flag += add_flag
        return cur_flag
    else:  # Array
        mods = cur_flag % add_flag
        upd = mods < (add_flag//2)
        cur_flag[upd] += add_flag
        return cur_flag


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
    # One way
    for key in newdb.keys():
        try:
            assert key in maindb.keys()
        except AssertionError:
            pdb.set_trace()
            return False
    return True


def get_new_ids(maindb, newdb, toler=2*u.arcsec, chk=True):
    """ Generate new IGMsp_IDs for an input DB
    Parameters
    ----------
    maindb : Table
    newdb : Table
    chk : bool, optional
      Perform some checks

    Returns
    -------
    ids : ndarray (int)
      Old IDs are filled with negative their value
      New IDs are generated as needed

    """
    IDs = np.zeros(len(newdb), dtype=int)
    # Setup
    c_main = SkyCoord(ra=maindb['RA'], dec=maindb['DEC'], unit='deg')
    c_new = SkyCoord(ra=newdb['RA'], dec=newdb['DEC'], unit='deg')
    # Find new sources
    idx, d2d, d3d = match_coordinates_sky(c_new, c_main, nthneighbor=1)
    new = d2d > toler
    #
    # Old IDs
    IDs[~new] = -1 * maindb[idx[~new]]['IGMsp_ID']
    nnew = np.sum(new)
    # New IDs
    if nnew > 0:
        # Generate
        newID = np.max(maindb['IGMsp_ID']) + 1
        newIDs = newID + np.arange(nnew, dtype=int)
        # Insert
        IDs[new] = newIDs
    if chk:
        print("The following sources were previously in the DB")
        print(newdb[~new])
    # Return
    return IDs



def ver01():
    """ Build version 0.1
    Returns
    -------

    """
    # HDF5 file
    outfil = igmspec.__path__[0]+'/../DB/IGMspec_DB_ver01.hdf5'
    hdf = h5py.File(outfil,'w')

    # Defs
    zpri = z_priority()
    lenz = [len(zpi) for zpi in zpri]
    dummyf = str('#')*np.max(np.array(lenz))  # For the Table

    # Main DB Table  (WARNING: THIS MAY TURN INTO SQL)
    idict = dict(RA=0., DEC=0., IGMsp_ID=0, EPOCH=2000.,
                 zem=0., sig_zem=0., flag_zem=dummyf, flag_survey=0)
    tkeys = idict.keys()
    lst = [[idict[tkey]] for tkey in tkeys]
    maindb = Table(lst, names=tkeys)

    ''' BOSS_DR12 '''
    # Read
    boss_meta = boss.meta_for_build()
    nboss = len(boss_meta)
    # IDs
    boss_ids = np.arange(nboss,dtype=int)
    boss_meta.add_column(Column(boss_ids, name='IGMsp_ID'))
    # Survey flag
    flag_s = survey_flag('BOSS_DR12')
    boss_meta.add_column(Column([flag_s]*nboss, name='flag_survey'))
    # Check
    assert chk_maindb_join(maindb, boss_meta)
    # Append
    maindb = vstack([maindb,boss_meta], join_type='exact')
    maindb = maindb[1:]  # Eliminate dummy line
    # Update hf5 file (TBD)

    ''' HD-LLS '''
    sname = 'HD-LLS_DR1'
    # Read
    hdlls_meta = hdlls.meta_for_build()
    nhdlls = len(hdlls_meta)
    # IDs
    hdlls_ids = get_new_ids(maindb, hdlls_meta)
    # Truncate
    new = hdlls_ids > 0
    hdlls_meta = hdlls_meta[new]
    nnew = len(hdlls_meta)
    hdlls_meta.add_column(Column(hdlls_ids[new], name='IGMsp_ID'))
    # Reset IDs to all positive
    hdlls_ids = np.abs(hdlls_ids)
    # Survey flag
    flag_s = survey_flag(sname)
    hdlls_meta.add_column(Column([flag_s]*nnew, name='flag_survey'))
    midx = np.array(maindb['IGMsp_ID'][hdlls_ids[~new]])
    maindb['flag_survey'][midx] += flag_s   # ASSUMES NOT SET ALREADY
    # Append
    assert chk_maindb_join(maindb, hdlls_meta)
    maindb = vstack([maindb,hdlls_meta], join_type='exact')
    pdb.set_trace()
    # Update hf5 file
    hdlls.hdf5_adddata(hdf, hdlls_ids, sname)

    # Finish
    hdf['catalog'] = maindb
    hdf.close()
    print("Wrote {:s} DB file".format(outfil))

if __name__ == '__main__':
    ver01()

