""" Module for key definitions in the IGMspec database
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np

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
        # Keck/ESI spectrometer -- ECH
        'ESI': dict(gratings=['ECH']),
        # Keck/LRIS spectrometer
        'LRISb': dict(gratings=['ECH']),
        'LRISr': dict(gratings=['ECH']),
        # Magellan MIKE spectrometer
        'MIKE': dict(gratings=['BOTH']),   # HD-LLS spliced blue and red
        'MIKEb': dict(gratings=['BLUE']),
        'MIKEr': dict(gratings=['RED']),
        # Magellan MagE spectrometer
        'MagE': dict(gratings=['N/A']),
        # MMT BCS
        'MMT': dict(gratings=['??']),
        # Gemini GMOS spectrometer
        'GMOS-S': dict(gratings=['R400', 'B600']),
        'GMOS-N': dict(gratings=['R400', 'B600']),
        'GNIRS': dict(gratings=['ECH']),
        # UKST
        '2dF': dict(gratings=['300B']),
        # HST
        'ACS': dict(gratings=['PR200L']),
        'WFC3': dict(gratings=['G280']),
        # VLT
        'XSHOOTER': dict(gratings=['UVB,VIS,NIR']),
        'ISAAC': dict(gratings=['SW_MRes']),
    }
    return inst_dict

def list_of_stypes():
    """ List of source types
    Returns
    -------
    stypes : list

    """
    stypes = [
        str('QSO'),        # Quasars
        str('GRB'),        # Gamma ray burst
    ]
    return stypes

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


def get_db_table_format():
    """ Returns DB Table format

    Returns
    -------
    idict : dict
      Describes the table columns
    """
        # Defs
    zpri = z_priority()
    lenz = [len(zpi) for zpi in zpri]
    dummyf = str('#')*np.max(np.array(lenz))  # For the Table
    stypes = list_of_stypes()
    lens = [len(stype) for stype in stypes]
    dummys = str('#')*np.max(np.array(lens))  # For the Table

    # Dict for Table
    idict = dict(RA=0., DEC=0., IGM_ID=0, zem=0., sig_zem=0.,
                 flag_zem=dummyf, flag_survey=0, STYPE=dummys)
    # Return
    return idict

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
    survey_dict['2QZ'] = 2**7       # Croom et al.
    survey_dict['ESI_DLA'] = 2**8   # Rafelski et al. 2012, 2014
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
    ESI_Rdict = {'0.50_arcsec': 4545./0.5, '0.75_arcsec': 4545./0.75,
                 '1.00_arcsec': 4545./1., 'Unknown': 4545./1.}
    HIRES1 = 36000.*1.148  # https://koa.ipac.caltech.edu/UserGuide/deckname_detail.html
    HIRES_Rdict = {'C1': HIRES1/0.861,
                   'C5': HIRES1/1.148,
                   'B2': HIRES1/0.574,
                   'B5': HIRES1/0.861,
                   'E3': HIRES1/0.4,
                   }
    LRISb_Rdict = {'400/3400': 500.}      # Assumes 1" slit
    LRISr_Rdict = {'600/7500': 1595.}     # Assumes 1" slit
    MagE_Rdict = {'0.70': 4100./0.7}
    GMOS_Rdict = {'B600+_G5307': 844.,    # Assumes 1" slit
                  'B600+_G5323': 844.,
                  'B1200+_G5301': 1872.,
                  }
    GNIRS_Rdict = {'32/mm_G5506': 1700.,    # Assumes 0.3" slit
                  }
    #
    Rdicts = dict(ESI=ESI_Rdict, HIRES=HIRES_Rdict, MagE=MagE_Rdict,
                  GMOS=GMOS_Rdict, GNIRS=GNIRS_Rdict, LRISb=LRISb_Rdict,
                  LRISr=LRISr_Rdict,
                  )
    #
    return Rdicts


def slit_dict():
    """ Slit dict

    Returns
    -------
    sdict : dict
      Translates slit mask name to slit with in arcsec

    """
    sdict = {'long_1.0':1.0}
    return sdict


def get_req_clms():
    """
    Returns
    -------
    req_clms : list
      List of required columns for meta data
    """
    req_clms = ['RA', 'DEC', 'EPOCH', 'zem', 'R', 'WV_MIN',
            'WV_MAX', 'DATE-OBS', 'SURVEY_ID', 'NPIX', 'SPEC_FILE',
            'INSTR', 'GRATING', 'TELESCOPE', 'IGM_ID']
    return req_clms
