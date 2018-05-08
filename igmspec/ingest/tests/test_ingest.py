# Module to run tests on ingest scripts

import os
import pytest

from igmspec.ingest.hdlls import grab_meta

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)


def test_hdlls():
    if os.getenv('RAW_IGMSPEC') is None:
        assert True
        return
    #
    meta = grab_meta()
    assert len(meta) == 145

