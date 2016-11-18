""" Module for key definitions in the IGMspec database
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import pdb

from collections import OrderedDict
from astropy import units as u


def z_priority():
    """ List of redshift priorities for setting the DB redshift
    See also myers.zbest_myers

    Returns
    -------
    zpri : list

    """
    zpri = [
        str('GGG'),        # GGG redshifts
        str('SDSS-HW'),    # SDSS redshifts with Hewitt&Wild
        str('BOSS_PCA'),   # PCA analysis by Paris et al. 2015 on BOSS spectra
        str('XQ-100'),     # XQ-100 redshifts
        str('BOSS_PIPE'),  # BOSS Pipeline redshifts
        str('2QZ'),        #
        str('2SLAQ'),      #
        str('AUS'),
        str('AGES'),
        str('COSMOS'),
        str('FAN'),
        str('MMT'),
        str('PAPOVICH'),
        str('GLIKMAN'),
        str('MADDOX'),
        str('LAMOST'),
        str('MCGREER'),
        str('VCV'),
        str('ALLBOSS'),
        str('Dall08'),     # Dall'Aglio et al. 2008
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
    survey_dict['HST_z2'] = 2**5    # O'Meara et al. 2011
    survey_dict['XQ-100'] = 2**6    # Lopez et al. 2016
    survey_dict['HDLA100'] = 2**7   # Neeleman et al. 2013
    survey_dict['2QZ'] = 2**8       # Croom et al.
    survey_dict['ESI_DLA'] = 2**9   # Rafelski et al. 2012, 2014
    survey_dict['COS-Halos'] = 2**10 # Tumlinson et al. 2013
    survey_dict['COS-Dwarfs'] = 2**11 # Bordoloi et al. 2014
    survey_dict['HSTQSO'] = 2**12 # Ribaudo et al. 2011; Neeleman et al. 2016
    survey_dict['MUSoDLA'] = 2**13 # Jorgensen et al. 2013
    survey_dict['UVES_Dall'] = 2**14 # Dall'Aglio et al. 2008
    #
    return survey_dict

