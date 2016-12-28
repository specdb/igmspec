

from astropy import units as u
from astropy.table import QTable, Table, Column, hstack, vstack
from astropy.coordinates import SkyCoord, match_coordinates_sky
from xastropy.xutils import xdebug as xdb


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

## SDSS/BOSS data stuff
specfile = os.getenv('RAW_IGMSPEC') + '/sdss/specObj-dr12_trim.fits'
spec = Table.read(specfile)
# Read in select columns from DR12 photometry. This and the file above are aligned
posfile =  os.getenv('RAW_IGMSPEC') + '/photoPosPlate-dr12.fits'
phot = Table.read(posfile)

# Trim to QSO, Specprimary, spec2d called it a QSO, redshift flag cuts, sanity check on coords
itrim = (spec['SPECPRIMARY'] == 1) &\
        [('QSO' in q) for q in spec['CLASS']] &\
        (spec['ZWARNING'] < 5) &\
        (spec['PLUG_RA'] >= 0.0) & (spec['PLUG_RA'] <= 360.0) &\
        (spec['PLUG_DEC'] <= 90.0)
spec = spec[itrim]
phot = phot[itrim]
sdss_boss1 = hstack(spec,phot,join='exact')
