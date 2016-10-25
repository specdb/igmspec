# Module to run tests on scripts

import matplotlib
matplotlib.use('agg')  # For Travis

# TEST_UNICODE_LITERALS

import pytest
import os

from ..scripts import plot_igmspec, sdss_igmspec

#version = 'v01'
version = 'v02'

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


"""
def test_plot_spec():
    pargs = plot_igmspec.parser(['J000345.00-232346.5', '-s=HD-LLS_DR1'])
    db_file = data_path('IGMspec_DB_{:s}_debug.hdf5'.format(version))
    plot_igmspec.main(pargs, db_file=db_file, unit_test=True)


def test_sdss():
    pargs = sdss_igmspec.parser(['751', '354'])
    db_file = data_path('IGMspec_DB_{:s}_debug.hdf5'.format(version))
    sdss_igmspec.main(pargs, db_file=db_file, unit_test=True)
"""
