#!/usr/bin/env python

"""
Runs igmspec related activities
  Mainly loads and plots a requested spectrum
"""

import pdb

def parser(options=None):

    import argparse

    parser = argparse.ArgumentParser(description='igmspec script')
    parser.add_argument("obj", type=str, help="Name of the Object,e.g. J081240+320808")
    parser.add_argument("--toler", default=5., type=float, help="Maximum offset in arcsec [default=5.]")
    parser.add_argument("--meta", default=True, help="Show meta data? [default: True]", action="store_true")
    parser.add_argument("--survey", help="Name of Survey to use")
    parser.add_argument("--select", default=0, type=int, help="Name of Survey to use [default: 0]")
    parser.add_argument("--mplot", default=False, help="Use simple matplotlib plot [default: False]")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args, unit_test=False):
    """ Run
    """
    import sys

    from astropy import units as u
    from igmspec import query_catalog as iqcat
    from igmspec import cat_utils as icu
    from igmspec import interface_db as igidb

    # Query the catalog
    qcat = iqcat.QueryCatalog()
    ids = qcat.radial_search(args.obj, args.toler*u.arcsec)
    if len(ids) == 0:
        print("Try another source...")
    elif len(ids) == 1:
        # Grab catalog entry
        row = qcat.get_cat(ids)
    else:
        print("Your query gave multiple hits.  Taking the first")
        row = qcat.get_cat(ids[0])

    # Choose survey
    if args.survey is not None:
        survey = args.survey
    else:
        fs = row['flag_survey']
        surveys = icu.flag_to_surveys(fs)
        survey = surveys[0]
        if survey == 'BOSS_DR12':
            print("BOSS NOT READY YET;  SPECIFY A DIFFERENT SURVEY")
            return
        print("Using survey={:s}.  Specify a different one with --survey".format(survey))
        print("Here is the full list {}".format(surveys))

    # Load DB
    idb = igidb.InterfaceDB(verbose=False)
    # Load Meta?
    if args.meta:
        meta = idb.grab_meta(survey, ids, show=True)
    # Load spectra
    spec = idb.grab_spec(survey, ids)
    spec.select = args.select
    # Show  [may transition to xspec]
    if args.mplot:
        spec.plot()
    else:
        spec.plot(xspec=True)
