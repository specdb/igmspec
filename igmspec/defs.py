""" Module for key definitions in the IGMspec database
"""
from __future__ import print_function, absolute_import, division, unicode_literals

from collections import OrderedDict
from astropy import units as u

def instruments():
    """ Dict of allowed instruments (and definitions)
    Includes allowed gratings

    Returns
    -------
    inst_dict

    """
    inst_dict = {
        # Spectrograph for SDSS-III (BOSS) survey ; https://www.sdss3.org/instruments/boss_spectrograph.php
        'BOSS': dict(gratings=['BLUE', 'RED', 'BOTH']),
        # Spectrograph for SDSS-I/II survey; http://classic.sdss.org/dr7/instruments/spectrographs/index.html
        'SDSS': dict(gratings=['BLUE', 'RED', 'BOTH']),
        # Keck/HIRES spectrometer -- BLUE/RED refer to the cross-disperser
        'HIRES': dict(gratings=['BLUE', 'RED', 'BOTH']),
    }

def z_priority():
    """ List of redshift priorities for setting the DB redshift
    Returns
    -------
    zpri : list

    """
    zpri = [
        str('GGG'),        # GGG redshifts
        str('BOSS_PCA'),   # PCA analysis by Paris et al. 2015 on BOSS spectra
        str('BOSS_PIPE'),  # BOSS Pipeline redshifts
        str('UNKN'),       # Unknown
    ]
    return zpri

def get_cat_dict():
    """ Definitions for the catalog
    Returns
    -------

    """
    cdict = dict(match_toler=2*u.arcsec)
    return cdict


def get_survey_dict():
    """ Return the survey dict
    Returns
    -------

    """
    survey_dict = OrderedDict()
    survey_dict['BOSS_DR12'] = 1
    survey_dict['SDSS_DR7'] = 2
    survey_dict['KODIAQ_DR1'] = 4   # O'Meara et al. 2016
    survey_dict['HD-LLS_DR1'] = 8   # Prochaska et al. 2015
    survey_dict['GGG'] = 16         # Worseck et al. 201X
    #
    return survey_dict


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
    survey_dict = get_survey_dict()
    #
    return survey_dict[survey]

def get_res_dicts():
    """ Resolution dicts

    Returns
    -------
    Rdicts : dict
      dict of R dicts

    """
    ESI_Rdict = {'0.50_arcsec': 4545./0.5, '0.75_arcsec': 4545./0.75}
    HIRES1 = 36000.*1.148  # https://koa.ipac.caltech.edu/UserGuide/deckname_detail.html
    HIRES_Rdict = {'C1': HIRES1/0.861,
                   'C5': HIRES1/1.148,
                   'B2': HIRES1/0.574,
                   'B5': HIRES1/0.861,
                   'E3': HIRES1/0.4,
                   }
    #
    Rdicts = dict(ESI=ESI_Rdict, HIRES=HIRES_Rdict)
    #
    return Rdicts
