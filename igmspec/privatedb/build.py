""" Module to build a private DB
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, glob
import h5py
import pdb

from igmspec import defs

from astropy.table import Table, vstack, Column
from astropy.coordinates import SkyCoord
#from astropy import units as u

from linetools import utils as ltu


def grab_files(tree_root, skip_files=('c.fits', 'C.fits', 'e.fits', 'E.fits')):
    """ Generate a list of FITS files within the file tree

    Parameters
    ----------
    tree_root : str
      Top level path of the tree of FITS files

    Returns
    -------
    files : list
      List of FITS files
    skip_files : tuple
      List of file roots to skip as primary files when ingesting

    """
    walk = os.walk(tree_root)
    folders = ['.']
    pfiles = []
    while len(folders) > 0:
        # Search for fits files
        ofiles = []
        for folder in folders:
            ofiles += glob.glob(tree_root+folder+'/*.fits*')
            # Eliminate error and continua files
            for ofile in ofiles:
                flg = True
                # Ugly loop
                for skip_file in skip_files:
                    if skip_file in ofile:
                        flg = False
                if flg:
                    pfiles.append(ofile)
            # walk
        folders = next(walk)[1]
    # Return
    return pfiles


def mk_meta(files, fname=False, stype='QSO'):
    """ Generate a meta Table from an input list of files

    Parameters
    ----------
    files : list
      List of FITS files
    fname : bool, optional
      Attempt to parse RA/DEC from the file name
      Format must be
      SDSSJ######(.##)+/-######(.#)[x]
        where x cannot be a #. or +/-

    Returns
    -------
    meta : Table
      Meta table
    """
    from igmspec.igmspec import IgmSpec
    igmsp = IgmSpec(skip_test=True)
    #
    coordlist = []
    for ifile in files:
        if fname:
            # Starting index
            if 'SDSSJ' in ifile:
                i0 = ifile.find('SDSSJ')+4
            else:
                i0 = ifile.rfind('J')+1
            # Find end (ugly)
            for ii in range(i0+1,99999):
                if ifile[ii] in ('0','1','2','3','4','5','6','7','8','9',
                                 '.','+','-'):
                    continue
                else:
                    i1 = ii
                    break
        # Get coord
        coord = ltu.radec_to_coord(ifile[i0:i1])
        coordlist.append(coord)
    coords = SkyCoord(ra=[coord.ra.degree for coord in coordlist], dec=[coord.dec.degree for coord in coordlist], unit='deg')

    # Generate Meta Table
    idict = defs.get_db_table_format()
    idict['PRIV_ID'] = 0
    tkeys = idict.keys()
    lst = [[idict[tkey]] for tkey in tkeys]
    maindb = Table(lst, names=tkeys)

    # Fill
    meta = Table()
    meta['RA'] = coords.ra.deg
    meta['DEC'] = coords.dec.deg
    meta['STYPE'] = [stype]*len(coords)

    # Redshift
    pdb.set_trace()


    # Test
    maindb = vstack([maindb,meta], join_type='exact')
    maindb = maindb[1:]

    # Add other meta (as desired)
    return maindb




def ver01(test=False, mk_test_file=False):
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

    # Defs
    zpri = defs.z_priority()
    lenz = [len(zpi) for zpi in zpri]
    dummyf = str('#')*np.max(np.array(lenz))  # For the Table
    stypes = defs.list_of_stypes()
    lens = [len(stype) for stype in stypes]
    dummys = str('#')*np.max(np.array(lens))  # For the Table
    #cdict = defs.get_cat_dict()

    # Main DB Table  (WARNING: THIS MAY TURN INTO SQL)
    idict = dict(RA=0., DEC=0., IGM_ID=0, zem=0., sig_zem=0.,
                 flag_zem=dummyf, flag_survey=0, STYPE=dummys)
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
    if mk_test_file:
        maindb = maindb[1:100]  # Eliminate dummy line
    else:
        maindb = maindb[1:]  # Eliminate dummy line
    #if not test:
    #    boss.hdf5_adddata(hdf, sdss_ids, sname)


    # Finish
    hdf['catalog'] = maindb
    hdf['catalog'].attrs['EPOCH'] = 2000.
    hdf['catalog'].attrs['Z_PRIORITY'] = zpri
    hdf['catalog'].attrs['VERSION'] = version
    #hdf['catalog'].attrs['CAT_DICT'] = cdict
    #hdf['catalog'].attrs['SURVEY_DICT'] = defs.get_survey_dict()
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
    # Read v1.X
    v01file = igmspec.__path__[0]+'/../DB/IGMspec_DB_v01.hdf5'
    v01file_debug = igmspec.__path__[0]+'/tests/files/IGMspec_DB_v01_debug.hdf5'
    print("Loading v01")
    v01hdf = h5py.File(v01file,'r')
    maindb = v01hdf['catalog'].value

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
        # Copy orginal
        for key in v01hdf_debug.keys():
            if key == 'catalog':
                dmaindb = v01hdf_debug[key].value
            else:
                v01hdf_debug.copy(key, hdf)
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
    #hdf['catalog'].attrs['CAT_DICT'] = cdict
    #hdf['catalog'].attrs['SURVEY_DICT'] = defs.get_survey_dict()
    hdf.close()
    print("Wrote {:s} DB file".format(outfil))
