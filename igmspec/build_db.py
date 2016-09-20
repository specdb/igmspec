""" Module to build the hdf5 database file for IGMspec
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import igmspec
import os

import h5py
import numbers, json
import pdb

from igmspec import defs
from igmspec.ingest import boss, hdlls, kodiaq, ggg, sdss, hst_z2, myers, twodf, xq100
from igmspec.ingest import hdla100
from igmspec.ingest import esidla
from igmspec.ingest import cos_halos

from astropy.table import Table, vstack, Column
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u

from linetools import utils as ltu


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
    """Check that new data is consistent with existing table

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


def chk_for_duplicates(maindb):
    """ Generate new IGM_IDs for an input DB

    Parameters
    ----------
    maindb : Table

    Return
    ------
    result : bool
      * True = pass
      * False = fail
    """
    c_main = SkyCoord(ra=maindb['RA'], dec=maindb['DEC'], unit='deg')
    # Find candidate dups
    idx, d2d, d3d = match_coordinates_sky(c_main, c_main, nthneighbor=2)
    cand_dups = d2d < 2*u.arcsec
    # Finish
    if np.sum(cand_dups) > 0:
        return False
    else:
        return True


def get_new_ids(maindb, newdb, chk=True, idkey='IGM_ID'):
    """ Generate new IGM_IDs for an input DB

    Parameters
    ----------
    maindb : Table
    newdb : Table
    chk : bool, optional
      Perform some checks
    idkey : str, optional
      Key for ID

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
    # Old IDs
    IDs[~new] = -1 * maindb[idkey][idx[~new]]
    nnew = np.sum(new)
    # New IDs
    if nnew > 0:
        # Generate
        newID = np.max(maindb[idkey]) + 1
        newIDs = newID + np.arange(nnew, dtype=int)
        # Insert
        IDs[new] = newIDs
    if chk:
        print("The following sources were previously in the DB")
        print(newdb[~new])
    # Return
    return IDs


def set_new_ids(maindb, newdb, chk=True, idkey='IGM_ID'):
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
    ids = get_new_ids(maindb, newdb, idkey=idkey)  # Includes new and old
    # Truncate
    new = ids > 0
    cut_db = newdb[new]
    cut_db.add_column(Column(ids[new], name=idkey))
    # Reset IDs to all positive
    ids = np.abs(ids)
    #
    return cut_db, new, ids


def start_maindb(private=False):
    """ Start the main DB catalog

    Returns
    -------
    maindb : Table
    tkeys : list
      List of columns in the table
    private : bool, optional
      Private DB?

    """
    idict = defs.get_db_table_format()
    #if private:
    #    idict['PRIV_ID'] = 0
        #idict.pop('IGM_ID')
    tkeys = idict.keys()
    lst = [[idict[tkey]] for tkey in tkeys]
    maindb = Table(lst, names=tkeys)
    # Return
    return maindb, tkeys


def ver01(test=False, mk_test_file=False, **kwargs):
    """ Build version 1.0

    Parameters
    ----------
    test : bool, optional
      Run test only
    mk_test_file : bool, optional
      Generate the test file for Travis tests?
      Writes catalog and HD-LLS dataset only

    Returns
    -------

    """
    version = 'v01'
    # HDF5 file
    if mk_test_file:
        outfil = igmspec.__path__[0]+'/tests/files/IGMspec_DB_{:s}_debug.hdf5'.format(version)
        print("Building debug file: {:s}".format(outfil))
        test = True
    else:
        outfil = igmspec.__path__[0]+'/../DB/IGMspec_DB_{:s}.hdf5'.format(version)
    hdf = h5py.File(outfil,'w')

    ''' Myers QSOs '''
    myers.add_to_hdf(hdf)

    # Main DB Table
    maindb, tkeys = start_maindb()

    ''' BOSS_DR12 '''
    # Read
    sname = 'BOSS_DR12'
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
    if mk_test_file:
        maindb = maindb[1:100]  # Eliminate dummy line
    else:
        maindb = maindb[1:]  # Eliminate dummy line
    tmp=chk_for_duplicates(maindb)
    if not test:
        boss.hdf5_adddata(hdf, boss_ids, sname, **kwargs)

    ''' SDSS DR7'''
    sname = 'SDSS_DR7'
    print('===============\n Doing {:s} \n===============\n'.format(sname))
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
    if mk_test_file:
        sdss_cut = sdss_cut[0:100]
    maindb = vstack([maindb, sdss_cut], join_type='exact')
    # Update hf5 file
    #if not test:
    sdss.hdf5_adddata(hdf, sdss_ids, sname, **kwargs)

    ''' KODIAQ DR1 '''
    sname = 'KODIAQ_DR1'
    print('==================\n Doing {:s} \n==================\n'.format(sname))
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
    if not test:
        kodiaq.hdf5_adddata(hdf, kodiaq_ids, sname)

    ''' HD-LLS '''
    sname = 'HD-LLS_DR1'
    print('===============\n Doing {:s} \n==============\n'.format(sname))
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
    if (not test) or mk_test_file:
        hdlls.hdf5_adddata(hdf, hdlls_ids, sname, mk_test_file=mk_test_file)

    ''' GGG '''
    sname = 'GGG'
    print('===============\n Doing {:s} \n==============\n'.format(sname))
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
    if not mk_test_file:
        ggg.hdf5_adddata(hdf, ggg_ids, sname)

    # Check for duplicates
    if not chk_for_duplicates(maindb):
        raise ValueError("Failed duplicates")

    # Check for junk
    zpri = defs.z_priority()

    # Finish
    hdf['catalog'] = maindb
    hdf['catalog'].attrs['EPOCH'] = 2000.
    hdf['catalog'].attrs['Z_PRIORITY'] = zpri
    hdf['catalog'].attrs['VERSION'] = version
    #hdf['catalog'].attrs['CAT_DICT'] = cdict
    hdf['catalog'].attrs['SURVEY_DICT'] = json.dumps(ltu.jsonify(defs.get_survey_dict()))
    hdf.close()
    print("Wrote {:s} DB file".format(outfil))


def ver02(test=False, mk_test_file=False, skip_copy=False):
    """ Build version 2.X

    Reads previous datasets from v1.X

    Parameters
    ----------
    test : bool, optional
      Run test only
    mk_test_file : bool, optional
      Generate the test file for Travis tests?
      Writes catalog and HD-LLS dataset only
    skip_copy : bool, optional
      Skip copying the data from v01

    Returns
    -------
    """
    import os
    # Read v1.X
    #v01file = igmspec.__path__[0]+'/../DB/IGMspec_DB_v01.hdf5'
    v01file = os.getenv('IGMSPEC_DB')+'/IGMspec_DB_v01.hdf5'
    v01file_debug = igmspec.__path__[0]+'/tests/files/IGMspec_DB_v01_debug.hdf5'
    print("Loading v01")
    v01hdf = h5py.File(v01file,'r')
    maindb = Table(v01hdf['catalog'].value)

    # Start new file
    version = 'v02'
    if mk_test_file:
        outfil = igmspec.__path__[0]+'/tests/files/IGMspec_DB_{:s}_debug.hdf5'.format(version)
        print("Building debug file: {:s}".format(outfil))
        test = True
    else:
        outfil = igmspec.__path__[0]+'/../DB/IGMspec_DB_{:s}.hdf5'.format(version)
    hdf = h5py.File(outfil,'w')

    # Copy over the old stuff
    if (not test) and (not skip_copy):
        for key in v01hdf.keys():
            if key == 'catalog':
                continue
            else:
                v01hdf.copy(key, hdf)
    if mk_test_file:
        v01hdf_debug = h5py.File(v01file_debug,'r')
        # Copy original
        for key in v01hdf_debug.keys():
            if key == 'catalog':
                dmaindb = v01hdf_debug[key].value
            else:
                v01hdf_debug.copy(key, hdf)
        # Add subset of quasars
        idx = np.array([False]*v01hdf['quasars'].size)
        idx[0:100] = True
        idx[161121] = True
        idx[161130] = True
        hdf['quasars'] = v01hdf['quasars'].value[idx]
        # Add some SDSS for script test
        bsdssi = np.where(maindb['flag_survey'] == 3)[0][0:10]
        sdss_meta = v01hdf['SDSS_DR7']['meta']
        sdssi = np.in1d(maindb['IGM_ID'][bsdssi], sdss_meta['IGM_ID'])
        hdf.create_group('SDSS_DR7')
        ibool = np.array([False]*len(sdss_meta))
        ibool[sdssi] = True
        # Generate
        hdf['SDSS_DR7']['meta'] = sdss_meta[ibool]
        hdf['SDSS_DR7']['spec'] = v01hdf['SDSS_DR7']['spec'][ibool]
        # Finish
        test = True
        maindb = dmaindb

    ''' 2QZ '''
    if not mk_test_file:
        sname = '2QZ'
        print('===============\n Doing {:s} \n==============\n'.format(sname))
        # Read
        tdf_meta = twodf.meta_for_build()
        # IDs
        tdf_cut, new, tdf_ids = set_new_ids(maindb, tdf_meta)
        nnew = np.sum(new)
        # Survey flag
        flag_s = defs.survey_flag(sname)
        tdf_cut.add_column(Column([flag_s]*nnew, name='flag_survey'))
        midx = np.array(maindb['IGM_ID'][tdf_ids[~new]])
        maindb['flag_survey'][midx] += flag_s
        # Append
        assert chk_maindb_join(maindb, tdf_cut)
        maindb = vstack([maindb,tdf_cut], join_type='exact')
        # Update hf5 file
        twodf.hdf5_adddata(hdf, tdf_ids, sname)

    ''' COS-Halos '''
    if not mk_test_file:
        sname = 'COS-Halos'
        print('===============\n Doing {:s} \n==============\n'.format(sname))
        # Read
        chalos_meta = cos_halos.meta_for_build()
        # IDs
        chalos_cut, new, chalos_ids = set_new_ids(maindb, chalos_meta)
        nnew = np.sum(new)
        # Survey flag
        flag_s = defs.survey_flag(sname)
        chalos_cut.add_column(Column([flag_s]*nnew, name='flag_survey'))
        midx = np.array(maindb['IGM_ID'][chalos_ids[~new]])
        maindb['flag_survey'][midx] += flag_s
        # Append
        assert chk_maindb_join(maindb, chalos_cut)
        maindb = vstack([maindb, chalos_cut], join_type='exact')
        # Update hf5 file
        cos_halos.hdf5_adddata(hdf, chalos_ids, sname)#, mk_test_file=mk_test_file)

    ''' HDLA100 '''
    if not mk_test_file:
        sname = 'HDLA100'
        print('===============\n Doing {:s} \n==============\n'.format(sname))
        # Read
        hdla100_meta, _ = hdla100.meta_for_build()
        # IDs
        hdla100_cut, new, hdla100_ids = set_new_ids(maindb, hdla100_meta)
        nnew = np.sum(new)
        # Survey flag
        flag_s = defs.survey_flag(sname)
        hdla100_cut.add_column(Column([flag_s]*nnew, name='flag_survey'))
        midx = np.array(maindb['IGM_ID'][hdla100_ids[~new]])
        maindb['flag_survey'][midx] += flag_s
        # Append
        assert chk_maindb_join(maindb, hdla100_cut)
        maindb = vstack([maindb, hdla100_cut], join_type='exact')
        # Update hf5 file
        hdla100.hdf5_adddata(hdf, hdla100_ids, sname)#, mk_test_file=mk_test_file)

    ''' ESI-DLA '''
    if not mk_test_file:
        sname = 'ESI_DLA'
        print('===============\n Doing {:s} \n==============\n'.format(sname))
        # Read
        esidla_meta = esidla.meta_for_build()
        # IDs
        esidla_cut, new, esidla_ids = set_new_ids(maindb, esidla_meta)
        nnew = np.sum(new)
        # Survey flag
        flag_s = defs.survey_flag(sname)
        esidla_cut.add_column(Column([flag_s]*nnew, name='flag_survey'))
        midx = np.array(maindb['IGM_ID'][esidla_ids[~new]])
        maindb['flag_survey'][midx] += flag_s
        # Append
        assert chk_maindb_join(maindb, esidla_cut)
        maindb = vstack([maindb, esidla_cut], join_type='exact')
        # Update hf5 file
        esidla.hdf5_adddata(hdf, esidla_ids, sname)#, mk_test_file=mk_test_file)

    ''' XQ-100 '''
    #if not mk_test_file:
    if False:
        sname = 'XQ-100'
        print('===============\n Doing {:s} \n==============\n'.format(sname))
        # Read
        xq100_meta = xq100.meta_for_build()
        # IDs
        xq100_cut, new, xq100_ids = set_new_ids(maindb, xq100_meta)
        nnew = np.sum(new)
        # Survey flag
        flag_s = defs.survey_flag(sname)
        xq100_cut.add_column(Column([flag_s]*nnew, name='flag_survey'))
        midx = np.array(maindb['IGM_ID'][xq100_ids[~new]])
        maindb['flag_survey'][midx] += flag_s
        # Append
        assert chk_maindb_join(maindb, xq100_cut)
        maindb = vstack([maindb, xq100_cut], join_type='exact')
        # Update hf5 file
        xq100.hdf5_adddata(hdf, xq100_ids, sname)#, mk_test_file=mk_test_file)


    ''' HST_z2 '''
    if not mk_test_file:
        sname = 'HST_z2'
        print('===============\n Doing {:s} \n==============\n'.format(sname))
        # Read
        hstz2_meta = hst_z2.meta_for_build()
        # IDs
        hstz2_cut, new, hstz2_ids = set_new_ids(maindb, hstz2_meta)
        nnew = np.sum(new)
        if nnew > 0:
            raise ValueError("All of these should be in SDSS")
        # Survey flag
        flag_s = defs.survey_flag(sname)
        midx = np.array(maindb['IGM_ID'][hstz2_ids[~new]])
        maindb['flag_survey'][midx] += flag_s
        # Update hf5 file
        if (not test):# or mk_test_file:
            hst_z2.hdf5_adddata(hdf, hstz2_ids, sname, mk_test_file=mk_test_file)

    # Finish
    hdf['catalog'] = maindb
    hdf['catalog'].attrs['EPOCH'] = 2000.
    zpri = v01hdf['catalog'].attrs['Z_PRIORITY']
    hdf['catalog'].attrs['Z_PRIORITY'] = zpri
    hdf['catalog'].attrs['VERSION'] = version
    hdf['catalog'].attrs['SURVEY_DICT'] = json.dumps(ltu.jsonify(defs.get_survey_dict()))
    #hdf['catalog'].attrs['CAT_DICT'] = cdict
    hdf.close()
    print("Wrote {:s} DB file".format(outfil))
