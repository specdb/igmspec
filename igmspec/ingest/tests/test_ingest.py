# Module to run tests on ingest scripts

import os
import pytest

from ..hdlls import meta_for_build

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)


def test_hdlls():
    if os.getenv('RAW_IGMSPEC') is None:
        assert True
        return
    #
    meta = meta_for_build()
    assert len(meta) == 129

