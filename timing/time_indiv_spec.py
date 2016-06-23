""" Test time to load spectra one by one
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import cProfile, pstats

from astropy.coordinates import SkyCoord
from igmspec.igmspec import IgmSpec


def time_coord_to_spec(survey='HD-LLS_DR1', ntrials=1000, seed=123):
    """ Time the process of grabbing spectra given an input coordinate

    Parameters
    ----------
    survey : str, optional
      Survey to test on
    ntrials = int, optional

    Returns
    -------

    """
    # Init
    igmsp = IgmSpec()
    rstate = np.random.RandomState(seed)

    # Grab survey
    meta = igmsp.idb.grab_meta(survey)
    coords = SkyCoord(ra=meta['RA'], dec=meta['DEC'], unit='deg')

    rani = rstate.randint(0,len(meta),ntrials)
    # Loop
    for ii in rani:
        coord = coords[ii]
        # Grab
        speclist, meta = igmsp.spec_from_coord(coord, isurvey=[survey])


# Command line execution
if __name__ == '__main__':
    #cProfile.run('time_coord_to_spec(ntrials=100)')
    cProfile.run('time_coord_to_spec(ntrials=100)', 'coord_to_spec.stats')
    stats = pstats.Stats('coord_to_spec.stats')
    stats.strip_dirs()
    stats.sort_stats('cumulative')
    stats.print_stats()


