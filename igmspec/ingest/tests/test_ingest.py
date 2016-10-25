# Module to run tests on ingest scripts

import pytest
import numpy as np
import h5py

from astropy.table import Column

from specdb.build import utils as sdbbu

from igmspec.ingest import hst_z2
from igmspec.defs import get_survey_dict
survey_dict = get_survey_dict()

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)


def test_hstz2():
    # Setup
    outfil = 'unit_test.hdf5'
    hdf = h5py.File(outfil,'w')
    maindb, tkeys = sdbbu.start_maindb(extras=dict(IGM_ID=0))
    #
    sname = 'HST_z2'
    hstz2_meta = hst_z2.meta_for_build(unit_test=True)
    # IDs
    hstz2_cut, new, hstz2_ids = sdbbu.set_new_ids(maindb, hstz2_meta)
    nnew = np.sum(new)
    assert nnew == 69
    #if nnew > 0:
    #    raise ValueError("All of these should be in SDSS")
    # Survey flag
    flag_s = survey_dict[sname]
    hstz2_cut.add_column(Column([flag_s]*nnew, name='flag_survey'))
    midx = np.array(maindb['IGM_ID'][hstz2_ids[~new]])
    maindb['flag_survey'][midx] += flag_s
    # Update hf5 file
    hst_z2.hdf5_adddata(hdf, hstz2_ids, sname)
    # Test?
