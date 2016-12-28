""" Module to ingest Myers' QSOs
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import numpy as np
import os
#import pdb
from xastropy.xutils import xdebug as xdb
from astropy.table import Table
from astropy.io import fits
from astropy import units as u

# What should this matching tolerance be?? Set to 2.0" for now
def add_to_hdf(hdf, Z_MIN = 0.1, Z_MAX = 7.1, MATCH_TOL = 2.0*u.arcsec):

    from astropy import units as u
    from astropy.table import QTable, Table, Column, hstack, vstack
    from astropy.coordinates import SkyCoord, match_coordinates_sky, search_around_sky

    ## SDSS/BOSS data stuff
    specfile = os.getenv('RAW_IGMSPEC') + '/sdss/specObj-dr12_trim.fits'
    spec = Table.read(specfile)
    # Read in select columns from DR12 photometry. This and the file above are aligned
    posfile = os.getenv('RAW_IGMSPEC') + '/sdss/photoPosPlate-dr12_trim.fits'
    phot = Table.read(posfile)

    # Trim to QSO, Specprimary, spec2d called it a QSO, redshift flag cuts, sanity check on coords
    itrim = (spec['SPECPRIMARY'] == 1) & \
            [('QSO' in q) for q in spec['CLASS']] & \
            (spec['ZWARNING'] < 5) & \
            (spec['PLUG_RA'] >= 0.0) & (spec['PLUG_RA'] <= 360.0) & \
            (np.abs(spec['PLUG_DEC']) <= 90.0)
    spec = spec[itrim]
    phot = phot[itrim]
    sdss_boss1 = hstack([spec, phot], join_type='exact')
    # Add SDSS prefix to all SDSS tags
    for key in sdss_boss1.keys():
        sdss_boss1.rename_column(key, 'SDSS_' + key)

    # Read in the Myers file, match it to Myers sweeps photometry
    # Myers master QSO catalog
    ADM_file = os.getenv('RAW_IGMSPEC') + '/Myers/GTR-ADM-QSO-master-wvcv.fits'
    ADM_qso = Table.read(ADM_file)
    head1 = fits.open(ADM_file)[1].header
    DATE = head1['DATE']

    # Photometry for Myers QSO catalog. This file is not aligned with the catalog file, i.e. it is a
    # superset that includes the catalog file. For that reason we need to match and tack on photometry
    ADM_sweep_file = os.getenv('RAW_IGMSPEC') + '/Myers/GTR-ADM-QSO-master-sweeps-Feb5-2016.fits'
    ADM_sweep = Table.read(ADM_sweep_file)

    c_qso = SkyCoord(ra=ADM_qso['RA'] * u.deg, dec=ADM_qso['DEC'] * u.deg)
    c_swp = SkyCoord(ra=ADM_sweep['RA'] * u.deg, dec=ADM_sweep['DEC'] * u.deg)

    ## Create an aligned Table for matching photometry from sweeps
    nqso = len(ADM_qso)
    qso_phot = Table(np.repeat(np.zeros_like(ADM_sweep[0]), nqso))
    # Rename the RA and DEC
    qso_phot.rename_column('RA', 'RA_sweep')
    qso_phot.rename_column('DEC', 'DEC_sweep')
    # Cull out the keys which already exist in the ADM_qso Table (except
    # for the RA and DEC, which we renamed)
    dupe_keys = list(set(ADM_qso.keys()) & set(qso_phot.keys()))
    qso_phot.remove_columns(dupe_keys)

    idx, d2d, d3d = c_qso.match_to_catalog_sky(c_swp)
    # Currently using 1.0" for matching, as for the SDSS objects, these will mostly be the exact
    # same coordinates.
    itrim = (d2d <= 1.0 * u.arcsec)
    qso_phot[:][itrim] = ADM_sweep[:][idx[itrim]]
    ADM_qso = hstack([ADM_qso, qso_phot], join_type='exact')
    # Trim to only spectroscopic objects
    ispec = spectro_myers(ADM_qso)
    ADM_qso = ADM_qso[ispec]
    # assign best redshifts to ZEM tag
    zbest_myers(ADM_qso)
    # Add MYERS prefix to all MYERS tags
    for key in ADM_qso.keys():
        ADM_qso.rename_column(key, 'MYERS_' + key)

    # Now we meatch the SDSS/BOSS and Myers catalogs to create one master QSO catalog
    #
    # There are three groups of objects, 1) SDSS-MYERS match, 2) SDSS only, 3) Myers only.
    # Deal with each in turn.

    # 1) SDSS-MYERS match. Add Myers tags to the SDSS structure
    c_sdss = SkyCoord(ra=sdss_boss1['SDSS_PLUG_RA'] * u.deg, dec=sdss_boss1['SDSS_PLUG_DEC'] * u.deg)
    c_myers = SkyCoord(ra=ADM_qso['MYERS_RA'] * u.deg, dec=ADM_qso['MYERS_DEC'] * u.deg)
    isdss, imyers, d2d, _ = search_around_sky(c_sdss, c_myers, MATCH_TOL)
    sdss_myers = hstack([sdss_boss1[isdss], ADM_qso[imyers]], join_type='exact')
    sdss_myers['SDSS_MYERS_FLAG'] = 'SDSS_MYERS'
    sdss_myers['RA'] = sdss_myers['SDSS_PLUG_RA']
    sdss_myers['DEC'] = sdss_myers['SDSS_PLUG_DEC']
    sdss_myers['SOURCEBIT'] = sdss_myers['MYERS_SOURCEBIT']
    sdss_myers['ZEM'] = sdss_myers['MYERS_ZEM']
    sdss_myers['ZEM_SOURCE'] = sdss_myers['MYERS_ZEM_SOURCE']

    # 2) SDSS only
    # Find the SDSS objects that have no match in the Myers catalog
    inomatch = np.ones(len(c_sdss), dtype=bool)
    inomatch[isdss] = False
    sdss_only = sdss_boss1[inomatch]
    sdss_only['SDSS_MYERS_FLAG'] = 'SDSS_ONLY'
    sdss_only['RA'] = sdss_only['SDSS_PLUG_RA']
    sdss_only['DEC'] = sdss_only['SDSS_PLUG_DEC']
    sdss_only['SOURCEBIT'] = 2 ** 19  # New source bit for SDSS only objects
    sdss_only['ZEM'] = sdss_only['SDSS_Z']
    sdss_only['ZEM_SOURCE'] = 'SDSS_ONLY'

    # 3) Myers only
    # Find the Myers objects that have no match in SDSS/BOSS
    inomatch = np.ones(len(c_myers), dtype=bool)
    inomatch[imyers] = False
    myers_only = ADM_qso[inomatch]
    myers_only['SDSS_MYERS_FLAG'] = 'MYERS_ONLY'
    myers_only['RA'] = myers_only['MYERS_RA']
    myers_only['DEC'] = myers_only['MYERS_DEC']
    myers_only['SOURCEBIT'] = myers_only['MYERS_SOURCEBIT']
    myers_only['ZEM'] = myers_only['MYERS_ZEM']
    myers_only['ZEM_SOURCE'] = myers_only['MYERS_ZEM_SOURCE']

    sdss_myers_out = vstack([sdss_myers, sdss_only, myers_only])

    # Cut down
    ztrim = (sdss_myers_out['ZEM'] >= Z_MIN) & (sdss_myers_out['ZEM'] <= Z_MAX)
    coordtrim = (sdss_myers_out['RA'] >= 0.0) & (sdss_myers_out['RA'] <= 360.0) & (np.abs(
        sdss_myers_out['DEC']) <= 90.0)
    keep = ztrim & coordtrim
    sdss_myers_out = sdss_myers_out[keep]
    hdf['quasars'] = sdss_myers_out
    hdf['quasars'].attrs['DATE'] = DATE

    return


def add_to_hdf_old(hdf):
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

def spectro_myers(ADM_qso):
    """ Returns indices of objects in the Myers catalog which are real spectroscopic QSOs.

     Parameters
     ----------
     ADM_qso : Table
       Myers catalog without ZEM, ZEM_SOURCE columns

     Returns
     -------
     Aligned array of booleans with True indicating a spectroscopic QSO.

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

     Notes from discussion with Myers:

    For this reason we now exclude anyting with only bit 18 set
    --------
    I've also added, as bit 2L^18, a list of everything that
    was visually inspected for BOSS. The redshifts for these objects
    aren't necessarily correct, but if this bit is set for an
    object and the object does not have a redshift corresponding
    to bit 2L^0 or bit 2L^7 then this is not a quasar, as it was
    visually inspected and not ultimately
    included in a quasar catalog.


    """

    ispec = ((ADM_qso['SOURCEBIT'] & 2**0) != False) | \
            ((ADM_qso['SOURCEBIT'] & 2**1) != False) | \
            ((ADM_qso['SOURCEBIT'] & 2**2) != False) | \
            ((ADM_qso['SOURCEBIT'] & 2**3) != False) | \
            ((ADM_qso['SOURCEBIT'] & 2**4) != False) | \
            ((ADM_qso['SOURCEBIT'] & 2**5) != False) | \
            ((ADM_qso['SOURCEBIT'] & 2**6) != False) | \
            ((ADM_qso['SOURCEBIT'] & 2**7) != False) | \
            ((ADM_qso['SOURCEBIT'] & 2**8) != False) | \
            ((ADM_qso['SOURCEBIT'] & 2**11) != False) | \
            ((ADM_qso['SOURCEBIT'] & 2**12) != False) | \
            ((ADM_qso['SOURCEBIT'] & 2**13) != False) | \
            ((ADM_qso['SOURCEBIT'] & 2**14) != False) | \
            ((ADM_qso['SOURCEBIT'] & 2**16) != False) | \
            ((ADM_qso['SOURCEBIT'] & 2**17) != False) & \
            (ADM_qso['SOURCEBIT'] != 2**18)

    return ispec


