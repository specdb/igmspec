""" Module to for ingest utilities
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import datetime
import pdb

def chk_meta(meta):
    """ Vettes a meta Table prior to its being ingested into the hdf

    Parameters
    ----------
    meta

    Returns
    -------
    chk : bool

    """
    chk = True
    # Required columns
    req_clms = ['IGM_ID', 'RA', 'DEC', 'EPOCH', 'zem', 'R', 'WV_MIN',
                'WV_MAX', 'DATE-OBS', 'SURVEY_ID', 'NPIX', 'SPEC_FILE',
                'INSTR', 'GRATING', 'TELESCOP']
    meta_keys = meta.keys()
    for clm in req_clms:
        if clm not in meta_keys:
            chk = False
            print("Missing column {:s} in meta".format(clm))
    # Check date formatting
    for row in meta:
        try:
            tval = datetime.datetime.strptime(row['DATE-OBS'], '%Y-%b-%d')
        except:
            print("Bad or missing DATE-OBS value")
            chk = False
    # Check for unicode
    # Return
    return chk


def set_resolution(head, instr=None):
    """ Sets resolution based on the instrument and header

    Parameters
    ----------
    head : FITS header
    instr : str, optional
      If not provided, attempt to grab from header

    Returns
    -------

    """
    from igmspec import defs
    # Dicts
    Rdicts = defs.get_res_dicts()
    # Grab instrument
    if instr is None:
        if 'CURRINST' in head.keys():  # ESI
            instr = head['CURRINST'].strip()
        elif 'INSTRUME' in head.keys():
            if 'HIRES' in head['INSTRUME']:
                instr = 'HIRES'
        else:
            pass
        if instr is None:
            raise ValueError("NEED MORE INFO FOR INSTR")

    # Grab resolution
    if instr == 'ESI':
        try:
            return Rdicts['ESI'][head['SLMSKNAM']]
        except KeyError:
            pdb.set_trace()
    elif instr == 'HIRES':
        try:
            return Rdicts['HIRES'][head['DECKNAME'].strip()]
        except KeyError:
            print("Need to add {:s}".format(head['DECKNAME']))
            pdb.set_trace()
    else:
        raise IOError("Not read for this instrument")
