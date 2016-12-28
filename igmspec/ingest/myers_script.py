

from astropy import units as u
from astropy.table import QTable, Table, Column, hstack, vstack
from astropy.coordinates import SkyCoord, match_coordinates_sky, search_around_sky
from xastropy.xutils import xdebug as xdb
from igmspec.ingest.myers import spectro_myers, zbest_myers

MATCH_TOL=3.0 # What should this matching tolerance be?? Set to 3.0" for
## SDSS/BOSS data stuff
specfile = os.getenv('RAW_IGMSPEC') + '/sdss/specObj-dr12_trim.fits'
spec = Table.read(specfile)
# Read in select columns from DR12 photometry. This and the file above are aligned
posfile =  os.getenv('RAW_IGMSPEC') + '/sdss/photoPosPlate-dr12_trim.fits'
phot = Table.read(posfile)

# Trim to QSO, Specprimary, spec2d called it a QSO, redshift flag cuts, sanity check on coords
itrim = (spec['SPECPRIMARY'] == 1) &\
        [('QSO' in q) for q in spec['CLASS']] &\
        (spec['ZWARNING'] < 5) &\
        (spec['PLUG_RA'] >= 0.0) & (spec['PLUG_RA'] <= 360.0) &\
        (np.abs(spec['PLUG_DEC']) <= 90.0)
spec = spec[itrim]
phot = phot[itrim]
sdss_boss1 = hstack([spec,phot],join_type='exact')
# Add SDSS prefix to all SDSS tags
for key in sdss_boss1.keys():
    sdss_boss1.rename_column(key,'SDSS_' + key)


# Read in the Myers file, match it to Myers sweeps
ADM_file = os.getenv('RAW_IGMSPEC') + '/Myers/GTR-ADM-QSO-master-wvcv.fits.gz'
ADM_qso = Table.read(ADM_file)

ADM_sweep_file = os.getenv('RAW_IGMSPEC') + '/Myers/GTR-ADM-QSO-master-sweeps-Feb5-2016.fits'
ADM_sweep = Table.read(ADM_sweep_file)

c_qso = SkyCoord(ra=ADM_qso['RA']*u.deg,dec=ADM_qso['DEC']*u.deg)
c_swp = SkyCoord(ra=ADM_sweep['RA']*u.deg,dec=ADM_sweep['DEC']*u.deg)

## Create an aligned Table for matching photometry from sweeps
nqso=len(ADM_qso)
qso_phot = Table(np.repeat(np.zeros_like(ADM_sweep[0]), nqso))
# Rename the RA and DEC
qso_phot.rename_column('RA','RA_sweep')
qso_phot.rename_column('DEC','DEC_sweep')
# Cull out the keys which already exist in the ADM_qso Table (except
# for the RA and DEC, which we renamed)
dupe_keys = list(set(ADM_qso.keys()) & set(qso_phot.keys()))
qso_phot.remove_columns(dupe_keys)

idx, d2d, d3d = c_qso.match_to_catalog_sky(c_swp)
itrim = (d2d <= 1.0*u.arcsec)
qso_phot[:][itrim]=ADM_sweep[:][idx[itrim]]
ADM_qso = hstack([ADM_qso,qso_phot],join_type='exact')
# Trim to only spectro objects
ispec = spectro_myers(ADM_qso)
ADM_qso = ADM_qso[ispec]
# assign best redshifts to ZEM tag
zbest_myers(ADM_qso)
# Add MYERS prefix to all MYERS tags
for key in ADM_qso.keys():
    ADM_qso.rename_column(key,'MYERS_' + key)

# There are three groups of objects, 1) SDSS-MYERS match, 2) SDSS only, 3) Myers only.
# Deal with each in turn.

# 1) SDSS-MYERS match. Add Myers tags to the SDSS structure
c_sdss = SkyCoord(ra=sdss_boss1['SDSS_PLUG_RA']*u.deg, dec = sdss_boss1['SDSS_PLUG_DEC']*u.deg)
c_myers= SkyCoord(ra=ADM_qso['MYERS_RA']*u.deg,dec=ADM_qso['MYERS_DEC']*u.deg)
isdss, imyers, d2d, _ = search_around_sky(c_sdss,c_myers,MATCH_TOL*u.arcsec)
sdss_myers = hstack([sdss_boss1[isdss], ADM_qso[imyers]],join_type='exact')
sdss_myers['SDSS_MYERS_FLAG'] = 'SDSS_MYERS'
sdss_myers['RA'] = sdss_myers['SDSS_PLUG_RA']
sdss_myers['DEC'] = sdss_myers['SDSS_PLUG_DEC']
sdss_myers['SOURCEBIT'] = sdss_myers['MYERS_SOURCEBIT']
sdss_myers['ZEM'] = sdss_myers['MYERS_ZEM']
sdss_myers['ZEM_SOURCE'] = sdss_myers['MYERS_ZEM_SOURCE']

# 2) SDSS only
# Find the SDSS objects that have no match in the Myers catalog
inomatch = np.ones(len(c_sdss),dtype=bool)
inomatch[isdss]=False
sdss_only = sdss_boss1[inomatch]
sdss_only['SDSS_MYERS_FLAG']='SDSS_ONLY'
sdss_only['RA'] = sdss_only['SDSS_PLUG_RA']
sdss_only['DEC'] = sdss_only['SDSS_PLUG_DEC']
sdss_only['SOURCEBIT'] = 2**19 # New source bit for SDSS only objects
sdss_only['ZEM'] = sdss_only['SDSS_Z']
sdss_only['ZEM_SOURCE'] = 'SDSS_ONLY'

# 3) Myers only
# Find the Myers objects that have no match in SDSS/BOSS
inomatch = np.ones(len(c_myers),dtype=bool)
inomatch[imyers]=False
myers_only = ADM_qso[inomatch]
myers_only['SDSS_MYERS_FLAG']='MYERS_ONLY'
myers_only['RA'] = myers_only['MYERS_RA']
myers_only['DEC'] = myers_only['MYERS_DEC']
myers_only['SOURCEBIT'] = myers_only['MYERS_SOURCEBIT']
myers_only['ZEM'] = myers_only['MYERS_ZEM']
myers_only['ZEM_SOURCE'] = myers_only['MYERS_ZEM_SOURCE']

sdss_myers_out = vstack([sdss_myers, sdss_only, myers_only])




