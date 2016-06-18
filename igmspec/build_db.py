""" Module to build the hdf5 database file for IGMspec
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import igmspec
import h5py
import numbers
import pdb

from igmspec import defs
from igmspec.ingest import boss, hdlls, kodiaq, ggg, sdss

from astropy.table import Table, vstack, Column
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


def get_new_ids(maindb, newdb, chk=True):
    """ Generate new IGM_IDs for an input DB
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
    cdict = defs.get_cat_dict()
    IDs = np.zeros(len(newdb), dtype=int)
    # Setup
    c_main = SkyCoord(ra=maindb['RA'], dec=maindb['DEC'], unit='deg')
    c_new = SkyCoord(ra=newdb['RA'], dec=newdb['DEC'], unit='deg')
    # Find new sources
    idx, d2d, d3d = match_coordinates_sky(c_new, c_main, nthneighbor=1)
    new = d2d > cdict['match_toler']
    #
    # Old IDs
    IDs[~new] = -1 * maindb[idx[~new]]['IGM_ID']
    nnew = np.sum(new)
    # New IDs
    if nnew > 0:
        # Generate
        newID = np.max(maindb['IGM_ID']) + 1
        newIDs = newID + np.arange(nnew, dtype=int)
        # Insert
        IDs[new] = newIDs
    if chk:
        print("The following sources were previously in the DB")
        print(newdb[~new])
    # Return
    return IDs


def set_new_ids(maindb, newdb, chk=True):
    """ Set the new IDs
    Parameters
    ----------
    maindb
    newdb
    toler
    chk

    Returns
    -------
    cut_db : Table
      Cut to the new QSOs
    new : bool array
    ids : ID values

    """
    # IDs
    ids = get_new_ids(maindb, newdb)  # Includes new and old
    # Truncate
    new = ids > 0
    cut_db = newdb[new]
    cut_db.add_column(Column(ids[new], name='IGM_ID'))
    # Reset IDs to all positive
    ids = np.abs(ids)
    #
    return cut_db, new, ids


def ver01():
    """ Build version 0.1
    Returns
    -------

    """
    # HDF5 file
    outfil = igmspec.__path__[0]+'/../DB/IGMspec_DB_ver01.hdf5'
    hdf = h5py.File(outfil,'w')

    # Defs
    zpri = defs.z_priority()
    lenz = [len(zpi) for zpi in zpri]
    dummyf = str('#')*np.max(np.array(lenz))  # For the Table
    cdict = defs.get_cat_dict()

    # Main DB Table  (WARNING: THIS MAY TURN INTO SQL)
    idict = dict(RA=0., DEC=0., IGM_ID=0, zem=0., sig_zem=0.,
                 flag_zem=dummyf, flag_survey=0)
    tkeys = idict.keys()
    lst = [[idict[tkey]] for tkey in tkeys]
    maindb = Table(lst, names=tkeys)

    ''' BOSS_DR12 '''
    # Read
    boss_meta = boss.meta_for_build()
    nboss = len(boss_meta)
    # IDs
    boss_ids = np.arange(nboss,dtype=int)
    boss_meta.add_column(Column(boss_ids, name='IGM_ID'))
    # Survey flag
    flag_s = defs.survey_flag('BOSS_DR12')
    boss_meta.add_column(Column([flag_s]*nboss, name='flag_survey'))
    # Check
    assert chk_maindb_join(maindb, boss_meta)
    # Append
    maindb = vstack([maindb,boss_meta], join_type='exact')
    maindb = maindb[1:]  # Eliminate dummy line
    #maindb = maindb[1:3]  # For testing

    ''' GGG '''
    sname = 'GGG'
    ggg_meta = ggg.meta_for_build()
    # IDs
    ggg_cut, new, ggg_ids = set_new_ids(maindb, ggg_meta)
    nnew = np.sum(new)
    # Survey flag
    flag_s = defs.survey_flag(sname)
    ggg_cut.add_column(Column([flag_s]*nnew, name='flag_survey'))
    midx = np.array(maindb['IGM_ID'][ggg_ids[~new]])
    maindb['flag_survey'][midx] += flag_s   # ASSUMES NOT SET ALREADY
    # Append
    assert chk_maindb_join(maindb, ggg_cut)
    maindb = vstack([maindb,ggg_cut], join_type='exact')
    # Update hf5 file
    ggg.hdf5_adddata(hdf, ggg_ids, sname)

    ''' SDSS DR7'''
    sname = 'SDSS_DR7'
    sdss_meta = sdss.meta_for_build()
    # IDs
    sdss_cut, new, sdss_ids = set_new_ids(maindb, sdss_meta)
    nnew = np.sum(new)
    # Survey flag
    flag_s = defs.survey_flag(sname)
    sdss_cut.add_column(Column([flag_s]*nnew, name='flag_survey'))
    midx = np.array(maindb['IGM_ID'][sdss_ids[~new]])
    maindb['flag_survey'][midx] += flag_s   # ASSUMES NOT SET ALREADY
    # Append
    assert chk_maindb_join(maindb, sdss_cut)
    maindb = vstack([maindb, sdss_cut], join_type='exact')
    # Update hf5 file
    #sdss.hdf5_adddata(hdf, sdss_ids, sname)

    ''' KODIAQ DR1 '''
    sname = 'KODIAQ_DR1'
    kodiaq_meta = kodiaq.meta_for_build()
    # IDs
    kodiaq_cut, new, kodiaq_ids = set_new_ids(maindb, kodiaq_meta)
    nnew = np.sum(new)
    # Survey flag
    flag_s = defs.survey_flag(sname)
    kodiaq_cut.add_column(Column([flag_s]*nnew, name='flag_survey'))
    midx = np.array(maindb['IGM_ID'][kodiaq_ids[~new]])
    maindb['flag_survey'][midx] += flag_s   # ASSUMES NOT SET ALREADY
    # Append
    assert chk_maindb_join(maindb, kodiaq_cut)
    maindb = vstack([maindb,kodiaq_cut], join_type='exact')
    # Update hf5 file
    kodiaq.hdf5_adddata(hdf, kodiaq_ids, sname)

    ''' HD-LLS '''
    sname = 'HD-LLS_DR1'
    # Read
    hdlls_meta = hdlls.meta_for_build()
    # IDs
    hdlls_cut, new, hdlls_ids = set_new_ids(maindb, hdlls_meta)
    nnew = np.sum(new)
    # Survey flag
    flag_s = defs.survey_flag(sname)
    hdlls_cut.add_column(Column([flag_s]*nnew, name='flag_survey'))
    midx = np.array(maindb['IGM_ID'][hdlls_ids[~new]])
    maindb['flag_survey'][midx] += flag_s   # ASSUMES NOT SET ALREADY
    # Append
    assert chk_maindb_join(maindb, hdlls_cut)
    maindb = vstack([maindb,hdlls_cut], join_type='exact')
    # Update hf5 file
    hdlls.hdf5_adddata(hdf, hdlls_ids, sname)

    """
    ''' GGG '''
    sname = 'GGG'
    ggg_meta = ggg.meta_for_build()
    # IDs
    ggg_cut, new, ggg_ids = set_new_ids(maindb, ggg_meta)
    nnew = np.sum(new)
    # Survey flag
    flag_s = defs.survey_flag(sname)
    ggg_cut.add_column(Column([flag_s]*nnew, name='flag_survey'))
    midx = np.array(maindb['IGM_ID'][ggg_ids[~new]])
    maindb['flag_survey'][midx] += flag_s   # ASSUMES NOT SET ALREADY
    # Append
    assert chk_maindb_join(maindb, ggg_cut)
    maindb = vstack([maindb,ggg_cut], join_type='exact')
    # Update hf5 file
    ggg.hdf5_adddata(hdf, ggg_ids, sname)
    """

    # Finish
    hdf['catalog'] = maindb
    hdf['catalog'].attrs['EPOCH'] = 2000.
    hdf['catalog'].attrs['Z_PRIORITY'] = zpri
    #hdf['catalog'].attrs['CAT_DICT'] = cdict
    #hdf['catalog'].attrs['SURVEY_DICT'] = defs.get_survey_dict()
    hdf.close()
    print("Wrote {:s} DB file".format(outfil))
