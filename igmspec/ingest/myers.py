""" Module to ingest Myers' QSOs
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os
import pdb

from astropy.table import Table
from astropy.io import fits


def add_to_hdf(hdf):
    """ Add Myers catalog to hdf file

    Parameters
    ----------
    hdf : HDF5 file
    """
    print("Adding Myers catalog")
    # Load
    ADM_qso, date = load()
    # Redshifts
    zbest_myers(ADM_qso)
    # Cut down
    ztrim = (ADM_qso['ZEM'] >= 0.1) & (ADM_qso['ZEM'] <= 7.0)
    coordtrim = (ADM_qso['RA'] >= 0.0) & (ADM_qso['RA'] <= 360.0) & (np.abs(
            ADM_qso['DEC']) <= 90.0)
    keep = ztrim & coordtrim
    ADM_qso = ADM_qso[keep]
    # Add
    hdf['quasars'] = ADM_qso
    hdf['quasars'].attrs['DATE'] = date
    #
    return


def load():
    """ Load catalog

    Parameters
    ----------

    Returns
    -------
    cat : Table
    date : str
      DATE of creation

    """
    ADM_file = os.getenv('RAW_IGMSPEC')+'/Myers/GTR-ADM-QSO-master-wvcv.fits.gz'
    ADM_qso = Table.read(ADM_file)
    # Grab header for DATE
    head1 = fits.open(ADM_file)[1].header
    # Return
    return ADM_qso, head1['DATE']


def zbest_myers(ADM_qso):
    """ Assign best redshift within the Myers catalog

    Parameters
    ----------
    ADM_qso : Table
      Myers catalog without ZEM, ZEM_SOURCE columns

    Returns
    -------
    Nothing; fills Myers catalog with ZEM, ZEM_SOURCE columns

     0 SDSS (Schneider et al. with Hewett and Wild redshifts)
     1 2QZ
     2 2SLAQ
     3 AUS
     4 AGES
     5 COSMOS
     6 FAN
     7 BOSS (Paris et al. through DR12+SEQUELS)
     8 MMT
     9 KDE (Photometric; Richards et al.)
    10 XDQSOZ (Photometric; Bovy et al.)
    11 PAPOVICH
    12 GLIKMAN
    13 MADDOX
    14 LAMOST
    15 VHS (Photometric; calculated using the Vista Hemisphere Survey IR-data)
    16 MCGREER
    17 VCV
    18 ALLBOSS
    """
    #nmyers = len(ADM_qso)
    #zstr = replicate(create_struct('ZEM', 0.0, 'ZEM_SOURCE', ''), nmyers)
    #myers = struct_addtags(a, zstr)
    #; Bits for Myers survey SOURCEBIT in order of redshift precedenece
    #;      HW  , BOSS , all the rest
    myers_pref = [0, 7, 1, 2, 3, 4, 5, 6, 8, 11, 12, 13, 14, 16, 17, 18]
    myers_binary = [2**ipref for ipref in myers_pref]
    #myers_binary = [2**0, 2**7, 2**1, 2**2, 2**3, 2**4, 2**5, 2**6, 2**8, 2**11,
    #                  2**12, 2**13, 2**14, 2**16, 2**17, 2**18]
    myers_source = ['SDSS-HW', 'BOSS_PCA', '2QZ', '2SLAQ', 'AUS', 'AGES', 'COSMOS', 'FAN', 'MMT', 'PAPOVICH',
                      'GLIKMAN', 'MADDOX', 'LAMOST', 'MCGREER', 'VCV', 'ALLBOSS']
    myers_source = [str(msrc) for msrc in myers_source]  # For hdf5
    #; Above gives top priority to HW, and second priority to BOSS

    # Assign the best redshift to Myers targets
    zem = []
    zem_source = []
    for row in ADM_qso:
        try:
            indx = min(np.where(row['SOURCEBIT'] & myers_binary)[0])
        except ValueError:
            indx = 0
        # Fill
        zem.append(row['ZBEST'][myers_pref[indx]])
        zem_source.append(myers_source[indx])
    # Add to Table
    ADM_qso['ZEM'] = zem
    ADM_qso['ZEM_SOURCE'] = zem_source

