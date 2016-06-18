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
    #parser.add_argument("-llist", default='ISM', action='store_true', help="Name of LineList:  ISM, HI, H2, CO, etc.")

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

    # Grab arguments
    pargs = parser(options=args)

    # Run
    if pargs.version is None:
        build_db.ver01()
    elif pargs.version == 'v01':
        build_db.ver01()
    else:
        raise IOError("Bad version number")

if __name__ == '__main__':
    main()
