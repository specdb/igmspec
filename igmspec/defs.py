""" Module for key definitions in the IGMspec database
"""
from __future__ import print_function, absolute_import, division, unicode_literals


def z_priority():
    """ List of redshift priorities for setting the DB redshift
    Returns
    -------
    zpri : list

    """
    zpri = [
        'BOSS_PCA',   # PCA analysis by Paris et al. 2015 on BOSS spectra
        'BOSS_PIPE',  # BOSS Pipeline redshifts
        'UNKN',       # Unknown
    ]
    return zpri

def survey_flag(survey, iflag=None):
    """ Defines bitwise survey flag for IGMspec
    Parameters
    ----------
    survey : str
      Name of the survey
    iflag : int, optional
    Returns
    -------
    flag_val : int

    """
    survey_dict = {'BOSS-DR12': 1}
    #
    return survey_dict[survey]

