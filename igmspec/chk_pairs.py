""" Module to check for pairs in igmspec
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb

from specdb.specdb import IgmSpec

from astropy import units as u
from astropy.coordinates import match_coordinates_sky, SkyCoord
from astropy.table import Table

def skip_gd_pair():
    """ Currently up-to-date with v02
    Returns
    -------
    stbl : Table
      Table of good pairs

    """
    skip = [[10.9697,   4.4073],
        [15.3189,   2.0326],
        [23.7452,  24.5140],
        [28.4283,  20.9148],
        [35.1739,   1.1985],
        [40.0218,  -0.6527],
        [41.2998,  -1.2216],
        [41.9780,   0.6380],
        [123.3026,  54.2806],
        [126.6732,  45.7450],
        [131.5022,   7.0747],
        [150.3681,  50.4663],
        [150.3362,  55.8989],  # FOS lens
        [158.2551,  47.2532],
        [164.0395,  55.2669],
        [170.1281,  54.7426],
        [176.7206,  16.7400],
        [188.5052,   6.5367],
        [190.7380,  25.7174],
        [193.7286,   8.7812],
        [196.9841,   4.3710],
        [198.7737,  47.9047],
        [201.3239,  37.6164],
        [211.2581,  44.8000],
        [222.7320,  47.0272],
        [238.3773,  22.5040],
        [243.2571,   8.1350],
        [253.7555,  26.0882],
        [357.0800,   0.9549],
        [116.9959,  43.3015],
        [184.6687,  50.2621],
        [166.63912,   -18.35661],  # FOS lens
        [166.6396, -18.3567], # FOS lens
        [216.9947,  -1.3601],
        [9.9763, -27.4229],  # 2QZ pair
        [341.6578, -29.4963], # 2QZ pair
        ]
    # Table
    sa = np.array(skip)
    stbl = Table()
    stbl['RA'] = sa[:,0]
    stbl['DEC'] = sa[:,1]
    # Return
    return stbl


def chk_for_pairs(maindb, pair_sep=10*u.arcsec):
    """ Generate new IGM_IDs for an input DB

    Parameters
    ----------
    maindb : Table

    Return
    ------
    result : bool
      * True = pass
      * False = fail
    """
    c_main = SkyCoord(ra=maindb['RA'], dec=maindb['DEC'], unit='deg')
    # Find candidate dups
    idx, d2d, d3d = match_coordinates_sky(c_main, c_main, nthneighbor=2)
    cand_pairs = np.where(d2d < pair_sep)[0]
    # Finish
    print("There are {:d} potential pairs with separation theta<{:g}".format(len(cand_pairs)//2,pair_sep))
    return cand_pairs


def chk_v02(pair_sep=10*u.arcsec):
    """ Check v02 for pairs
    Returns
    -------

    """
    print("checking..")
    igmsp = IgmSpec()
    # Grab candidate pairs
    cpairs = chk_for_pairs(igmsp.qcat.cat, pair_sep=pair_sep)
    # Coords
    c_main = SkyCoord(ra=igmsp.qcat.cat['RA'], dec=igmsp.qcat.cat['DEC'], unit='deg')
    # Skip
    stbl = skip_gd_pair()

    # Loop
    flg_cp = np.array([False]*len(igmsp.qcat.cat))
    for qq, cpair in enumerate(cpairs):
        # Skip those already done
        if flg_cp[cpair]:
            continue
        # Find the matches
        sep = c_main[cpair].separation(c_main)
        pairs = sep < pair_sep
        flg_cp[pairs] = True
        # Skip pairs with very different zem
        if np.sum(pairs) == 2:
            zem = igmsp.qcat.cat['zem'][pairs]
            if np.abs(zem[0]-zem[1]) > 0.1:
                continue
            # Both BOSS?
            if (igmsp.qcat.cat['flag_survey'][pairs][0] == 1.) & (
                igmsp.qcat.cat['flag_survey'][pairs][1] == 1.):
                continue
        # Skip table?
        if np.min(np.abs(igmsp.qcat.cat['RA'][pairs][0]-stbl['RA'])) < 1e-4:
            continue
        # XQ-100?
        if igmsp.qcat.cat['flag_survey'][pairs][1] == 64.:
            print("Skipping XQ-100")
            pdb.set_trace()
            continue
        # Print
        print('qq = {:d}'.format(qq))
        print(igmsp.qcat.cat[['RA','DEC','IGM_ID','zem','flag_survey']][pairs])
        print(sep.to('arcsec')[pairs])
        pdb.set_trace()
    # All clear?
    print("All clear..")


# Command line execution
if __name__ == '__main__':
    import sys

    if len(sys.argv) == 1: #
        flg = 0
        flg += 2**0  # v02
    else:
        flg = sys.argv[1]

    # SLLS Ions
    if (flg % 2**1) >= 2**0:
        #mk_lls_dr1(wrspec=False)
        chk_v02()