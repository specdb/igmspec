#!/usr/bin/env python
"""
Run a build of the DB
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import pdb

try:  # Python 3
    ustr = unicode
except NameError:
    ustr = str

def parser(options=None):
    import argparse
    # Parse
    parser = argparse.ArgumentParser(
        description='Build the igmspec DB')
    parser.add_argument("-v", "--version", help="DB version to generate")
    parser.add_argument("-t", "--test", default=False, action='store_true', help="Test?")
    parser.add_argument("--boss_hdf", help="HDF file with BOSS dataset [avoids repeating spectra ingestion]")
    parser.add_argument("--sdss_hdf", help="HDF file with SDSS dataset [avoids repeating spectra ingestion]")
    parser.add_argument("--clobber", default=False, action='store_true', help="Clobber existing file?")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args=None):
    """ Run
    Parameters
    ----------
    args

    Returns
    -------

    """
    from igmspec import build_db
    import h5py

    # Grab arguments
    pargs = parser(options=args)

    # BOSS
    if pargs.boss_hdf is not None:
        boss_hdf = h5py.File(pargs.boss_hdf,'r')
    else:
        boss_hdf = None
    if pargs.sdss_hdf is not None:
        if boss_hdf is not None:
            if pargs.boss_hdf==pargs.sdss_hdf:
                sdss_hdf = boss_hdf
            else:
                sdss_hdf = h5py.File(pargs.sdss_hdf,'r')
        else:
            sdss_hdf = h5py.File(pargs.sdss_hdf,'r')
    else:
        sdss_hdf = None

    # Run
    if pargs.version is None:
        print("Building v02 of the igmspec DB")
        build_db.ver02(test=pargs.test, clobber=pargs.clobber)
    elif pargs.version == 'v01':
        print("Building v01 of the igmspec DB")
        build_db.ver01(test=pargs.test,
                       boss_hdf=boss_hdf, sdss_hdf=sdss_hdf, clobber=pargs.clobber)
    elif pargs.version == 'v02':
        print("Building v02 of the igmspec DB")
        build_db.ver02(test=pargs.test, clobber=pargs.clobber)
    elif pargs.version == 'v02.1':
        print("Building v02.1 of the igmspec DB")
        build_db.ver02(test=pargs.test, clobber=pargs.clobber)
    else:
        raise IOError("Bad version number")

if __name__ == '__main__':
    main()
