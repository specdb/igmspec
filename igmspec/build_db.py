""" Module to build the hdf5 database file for IGMspec
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os, warnings

import h5py
import json
import datetime
import pdb
from collections import OrderedDict

from specdb import defs
from specdb.build import utils as sdbbu

import igmspec
from igmspec.ingest import boss, hdlls, kodiaq, ggg, sdss, hst_z2, myers, twodf, xq100
from igmspec.ingest import hdla100
from igmspec.ingest import esidla
from igmspec.ingest import cos_halos
from igmspec.ingest import hst_qso
from igmspec.ingest import hst_cooksey as hst_c
from igmspec.ingest import cos_dwarfs
from igmspec.ingest import musodla
from igmspec.ingest import uves_dall

from astropy.table import Table, vstack, Column
from astropy import units as u

from linetools import utils as ltu

from specdb.specdb import IgmSpec

#from igmspec.defs import get_survey_dict
#survey_dict = get_survey_dict()


def ver01(test=False, clobber=False, publisher='J.X. Prochaska', **kwargs):
    """ Build version 1.0

    Parameters
    ----------
    test : bool, optional
      Run test only

    Returns
    -------

    """
    version = 'v01'
    # HDF5 file
    outfil = igmspec.__path__[0]+'/../DB/IGMspec_DB_{:s}.hdf5'.format(version)
    # Chk clobber
    if os.path.isfile(outfil):
        if clobber:
            warnings.warn("Overwriting previous DB file {:s}".format(outfil))
        else:
            warnings.warn("Not overwiting previous DB file.  Use clobber=True to do so")
            return
    # Begin
    hdf = h5py.File(outfil,'w')

    ''' Myers QSOs '''
    myers.add_to_hdf(hdf)

    # Main DB Table
    idkey = 'IGM_ID'
    maindb, tkeys = sdbbu.start_maindb(idkey)

    # Group dict
    group_dict = {}

    # Organize for main loop
    groups = get_build_groups(version)

    pair_groups = ['SDSS_DR7']

    meta_only = False
    # Loop over the groups
    for gname in groups:
        # Meta
        if gname == 'SDSS_DR7':
            meta = groups[gname].grab_meta(hdf)
        else:
            meta = groups[gname].grab_meta()
        # Survey flag
        flag_g = sdbbu.add_to_group_dict(gname, group_dict)
        # IDs
        maindb = sdbbu.add_ids(maindb, meta, flag_g, tkeys, idkey,
                               first=(flag_g==1), close_pairs=(gname in pair_groups))
        # Spectra
        if not meta_only:
            groups[gname].hdf5_adddata(hdf, gname, meta, idkey)

    # Check for duplicates -- There is 1 pair in SDSS (i.e. 2 duplicates)
    if not sdbbu.chk_for_duplicates(maindb, dup_lim=2):
        raise ValueError("Failed duplicates")

    # Check for junk
    zpri = defs.z_priority()

    # Finish
    sdbbu.write_hdf(hdf, str('igmspec'), maindb, zpri,
                    group_dict, version, publisher=str(publisher))
    print("Wrote {:s} DB file".format(outfil))
    print("Update DB info in specdb.defs.dbase_info !!")


def ver02(test=False, skip_copy=False, publisher='J.X. Prochaska', clobber=False):
    """ Build version 2.X

    Reads previous datasets from v1.X

    Parameters
    ----------
    test : bool, optional
      Run test only
    skip_copy : bool, optional
      Skip copying the data from v01

    Returns
    -------
    """
    import os
    from specdb.specdb import IgmSpec
    # Read v01
    v01file = os.getenv('IGMSPEC_DB')+'/IGMspec_DB_v01.hdf5'
    #v01file_debug = igmspec.__path__[0]+'/tests/files/IGMspec_DB_v01_debug.hdf5'
    print("Loading v01")
    igmsp_v01 = IgmSpec(db_file=v01file)
    v01hdf = igmsp_v01.hdf
    maindb = igmsp_v01.cat.copy()

    # Start new file
    version = 'v02'
    outfil = igmspec.__path__[0]+'/../DB/IGMspec_DB_{:s}.hdf5'.format(version)
    # Clobber?
    if not chk_clobber(outfil, clobber=clobber):
        return
    # Begin
    hdf = h5py.File(outfil,'w')

    # Copy over the old stuff
    #skip_copy = True
    if (not test) and (not skip_copy):
        for key in v01hdf.keys():
            if key == 'catalog':
                continue
            else:
                v01hdf.copy(key, hdf)
    # Setup
    new_groups = get_build_groups(version)

    pair_groups = []
    group_dict = igmsp_v01.qcat.group_dict
    # Set/Check keys (and set idkey internally for other checks)
    idkey = 'IGM_ID'
    _, tkeys = sdbbu.start_maindb(idkey)
    mkeys = list(maindb.keys())
    for key in tkeys:
        assert key in mkeys

    meta_only = False
    # Loop over the groups
    for gname in new_groups:
        print("Working on group: {:s}".format(gname))
        # Meta
        meta = new_groups[gname].grab_meta()
        # Survey flag
        flag_g = sdbbu.add_to_group_dict(gname, group_dict)
        # IDs
        maindb = sdbbu.add_ids(maindb, meta, flag_g, tkeys, idkey,
                               first=(flag_g==1), close_pairs=(gname in pair_groups))
        # Spectra
        if not meta_only:
            new_groups[gname].hdf5_adddata(hdf, gname, meta)

    # Check for duplicates -- There is 1 pair in SDSS (i.e. 2 duplicates)
    if not sdbbu.chk_for_duplicates(maindb, dup_lim=2):
        raise ValueError("Failed duplicates")

    # Finish
    zpri = v01hdf['catalog'].attrs['Z_PRIORITY']
    sdbbu.write_hdf(hdf, str('igmspec'), maindb, zpri,
                    group_dict, version, publisher=str(publisher))

    print("Wrote {:s} DB file".format(outfil))
    print("Update DB info in specdb.defs.dbase_info !!")


def chk_clobber(outfil, clobber=False):
    """ Simple clobber check
    outfil : str
    clobber : bool, optional
    """
    # Chk clobber
    if os.path.isfile(outfil):
        if clobber:
            warnings.warn("Overwriting previous DB file {:s}".format(outfil))
            return True
        else:
            warnings.warn("Not overwiting previous DB file.  Set clobber=True to do so")
            return False
    else:
        return True

def get_build_groups(version):
    """
    Parameters
    ----------
    version : str

    Returns
    -------
    build_groups : dict

    """

    groups = OrderedDict()
    if version == 'v01':
        groups['BOSS_DR12'] = boss
        groups['SDSS_DR7'] = sdss
        groups['KODIAQ_DR1'] = kodiaq
        groups['HD-LLS_DR1'] = hdlls
        groups['GGG'] = ggg
    elif version == 'v02':
        groups = OrderedDict()
        groups['HST_z2'] = hst_z2       # O'Meara et al. 2011
        groups['XQ-100'] = xq100        # Lopez et al. 2016
        groups['HDLA100'] = hdla100     # Neeleman et al. 2013
        groups['2QZ'] = twodf           # Croom et al.
        groups['ESI_DLA'] = esidla      # Rafelski et al. 2012, 2014
        groups['COS-Halos'] = cos_halos # Tumlinson et al. 2013
        groups['COS-Dwarfs'] = cos_dwarfs # Bordoloi et al. 2014
        groups['HSTQSO'] = hst_qso      # Ribaudo et al. 2011; Neeleman et al. 2016
        groups['MUSoDLA'] = musodla     # Jorgensen et al. 2013
        groups['UVES_Dall'] = uves_dall # Dall'Aglio et al. 2008
        groups['UVpSM4'] = hst_c        # Cooksey et al. 2010, 2011
    else:
        raise IOError("Not ready for this version")
    # Return
    return groups
