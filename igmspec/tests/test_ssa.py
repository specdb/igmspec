# Module to run tests on scripts

# TEST_UNICODE_LITERALS

import pytest
import os

from specdb.specdb import IgmSpec
from specdb import ssa as spdb_ssa
import igmspec

#version = 'v01'
version = 'v02'

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_sdss_ssa_querydata():
    db_file = igmspec.__path__[0]+'/../DB/IGMspec_DB_v02.hdf5'
    igmsp = IgmSpec(db_file=db_file)
    #
    ssa = spdb_ssa.SSAInterface(igmsp)
    votable = ssa.querydata('0.027228,0.515341', SIZE=1e-3)
    # Write
    votable.to_xml('sdss_querydata.xml')




