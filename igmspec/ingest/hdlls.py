""" Module to ingest HD-LLS Survey data

Prochaska et al. 2015
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import pdb
import os, json, glob, imp
import datetime

from astropy.table import Table, Column
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
from astropy.io import fits

from linetools import utils as ltu

from specdb.build.utils import chk_meta, set_resolution

igms_path = imp.find_module('igmspec')[1]


def grab_meta_mike():
    """ Grab MIKE meta Table
    Returns
    -------

    """
    mike_file = igms_path+'/data/meta/HD-LLS_DR1_MIKE.ascii'
    mike_meta = Table.read(mike_file, format='ascii', delimiter='&',
                           guess=False, comment='#')
    # RA/DEC, DATE
    ra = []
    dec = []
    dateobs = []
    for row in mike_meta:
        # Fix DEC
        if '--' in row['sDEC']:
            row['sDEC'] = row['sDEC'].replace('--','-')
        # Get RA/DEC
        coord = ltu.radec_to_coord((row['sRA'],row['sDEC']))
        ra.append(coord.ra.value)
        dec.append(coord.dec.value)
        # DATE
        dvals = row['DATE'].split(' ')
        dateobs.append(str('{:s}-{:s}-{:s}'.format(dvals[2],dvals[1],dvals[0])))
    mike_meta.add_column(Column(ra, name='RA'))
    mike_meta.add_column(Column(dec, name='DEC'))
    mike_meta.add_column(Column(dateobs, name='DATE-OBS'))
    #
    return mike_meta


def meta_for_build():
    """ Generates the meta data needed for the IGMSpec build
    Returns
    -------
    meta : Table
    """
    # Cut down to unique QSOs
    hdlls_meta = Table.read(os.getenv('RAW_IGMSPEC')+'/HD-LLS_DR1/HD-LLS_DR1.fits')
    names = np.array([name[0:26] for name in hdlls_meta['Name']])
    uni, uni_idx = np.unique(names, return_index=True)
    hdlls_meta = hdlls_meta[uni_idx]
    nqso = len(hdlls_meta)
    #
    meta = Table()
    meta['RA'] = hdlls_meta['RA']
    meta['DEC'] = hdlls_meta['DEC']
    meta['zem'] = hdlls_meta['Z_QSO']
    meta['sig_zem'] = [0.]*nqso
    meta['flag_zem'] = [str('UNKN')]*nqso
    meta['STYPE'] = [str('QSO')]*nqso
    # Return
    return meta


def hdf5_adddata(hdf, IDs, sname, debug=False, chk_meta_only=False,
                 mk_test_file=False):
    """ Append HD-LLS data to the h5 file

    Parameters
    ----------
    hdf : hdf5 pointer
    IDs : ndarray
      int array of IGM_ID values in mainDB
    sname : str
      Survey name
    chk_meta_only : bool, optional
      Only check meta file;  will not write
    mk_test_file : bool, optional
      Generate the debug test file for Travis??

    Returns
    -------

    """
    from specdb import defs
    # Add Survey
    print("Adding {:s} survey to DB".format(sname))
    hdlls_grp = hdf.create_group(sname)
    # Load up
    Rdicts = defs.get_res_dicts()
    mike_meta = grab_meta_mike()
    mike_coord = SkyCoord(ra=mike_meta['RA'], dec=mike_meta['DEC'], unit='deg')
    # Checks
    if sname != 'HD-LLS_DR1':
        raise IOError("Not expecting this survey..")
    if np.sum(IDs < 0) > 0:
        raise ValueError("Bad ID values")
    # Open Meta tables
    cut_meta = meta_for_build()
    if len(cut_meta) != len(IDs):
        raise ValueError("Wrong sized table..")
    # DR1 Table by LLS, not spectrum
    hdlls_meta = Table.read(os.getenv('RAW_IGMSPEC')+'/HD-LLS_DR1/HD-LLS_DR1.fits')
    nlls = len(hdlls_meta)

    # Generate ID array from RA/DEC
    c_cut = SkyCoord(ra=cut_meta['RA'], dec=cut_meta['DEC'], unit='deg')
    c_all = SkyCoord(ra=hdlls_meta['RA'], dec=hdlls_meta['DEC'], unit='deg')
    # Find new sources
    idx, d2d, d3d = match_coordinates_sky(c_all, c_cut, nthneighbor=1)
    if np.sum(d2d > 0.1*u.arcsec):
        raise ValueError("Bad matches in HD-LLS")
    meta_IDs = IDs[idx]

    # Kludgy Table judo
    hdlls_full = hdlls_meta[0:1]
    spec_files = []

    full_IDs = []
    # Loop me to bid the full survey catalog
    for kk,row in enumerate(hdlls_meta):
        for spec_file in row['SPEC_FILES']:
            if spec_file == 'NULL':
                continue
            # Add to full table
            hdlls_full.add_row(row)
            spec_files.append(spec_file)
            full_IDs.append(meta_IDs[kk])
    # Build
    hdlls_full = hdlls_full[1:]
    hdlls_full.remove_column('SPEC_FILES')
    hdlls_full.add_column(Column(spec_files,name='SPEC_FILE'))
    hdlls_full.add_column(Column(full_IDs, name='IGM_ID'))
    # Cut on unique SPEC_FILEs
    uni, uni_idx = np.unique(np.array(spec_files), return_index=True)
    if sname == 'HD-LLS_DR1':
        # REMOVE ONE FILE (A DUPLICATE) BY HAND
        mt = uni != 'HD-LLS_J130756.70+042215.0_MIKE.fits'
        uni_idx = uni_idx[mt]
    hdlls_full = hdlls_full[uni_idx]
    full_coord = SkyCoord(ra=hdlls_full['RA'], dec=hdlls_full['DEC'], unit='deg')

    # Build spectra (and parse for meta)
    if mk_test_file:
        hdlls_full = hdlls_full[0:3]
    nspec = len(hdlls_full)
    max_npix = 210000  # Just needs to be large enough
    data = np.ma.empty((1,),
                       dtype=[(str('wave'), 'float64', (max_npix)),
                              (str('flux'), 'float32', (max_npix)),
                              (str('sig'),  'float32', (max_npix)),
                              #(str('co'),   'float32', (max_npix)),
                             ])
    # Init
    full_idx = np.zeros(len(hdlls_full), dtype=int)
    spec_set = hdf[sname].create_dataset('spec', data=data, chunks=True,
                                         maxshape=(None,), compression='gzip')
    spec_set.resize((nspec,))
    Rlist = []
    wvminlist = []
    wvmaxlist = []
    dateobslist = []
    npixlist = []
    instrlist = []
    gratinglist = []
    telelist = []
    # Loop
    members = glob.glob(os.getenv('RAW_IGMSPEC')+'/{:s}/*fits'.format(sname))
    kk = -1
    for jj,member in enumerate(members):
        if 'HD-LLS_DR1.fits' in member:
            continue
        kk += 1
        # Extract
        f = member
        hdu = fits.open(f)
        # Parse name
        fname = f.split('/')[-1]
        mt = np.where(hdlls_full['SPEC_FILE'] == fname)[0]
        if mk_test_file and (jj>=3):
            continue
        if len(mt) != 1:
            pdb.set_trace()
            raise ValueError("HD-LLS: No match to spectral file?!")
        else:
            print('loading {:s}'.format(fname))
            full_idx[kk] = mt[0]
        # npix
        head = hdu[0].header
        npix = head['NAXIS1']
        if npix > max_npix:
            raise ValueError("Not enough pixels in the data... ({:d})".format(npix))
        # Some fiddling about
        for key in ['wave','flux','sig']:
            data[key] = 0.  # Important to init (for compression too)
        # Double check
        if kk == 0:
            assert hdu[1].name == 'ERROR'
            assert hdu[2].name == 'WAVELENGTH'
        # Write
        data['flux'][0][:npix] = hdu[0].data
        data['sig'][0][:npix] = hdu[1].data
        data['wave'][0][:npix] = hdu[2].data
        # Meta
        wvminlist.append(np.min(data['wave'][0][:npix]))
        wvmaxlist.append(np.max(data['wave'][0][:npix]))
        npixlist.append(npix)
        if 'HIRES' in fname:
            instrlist.append('HIRES')
            telelist.append('Keck-I')
            gratinglist.append('BOTH')
            try:
                Rlist.append(set_resolution(head))
            except ValueError:
                # A few by hand (pulled from Table 1)
                if 'J073149' in fname:
                    Rlist.append(Rdicts['HIRES']['C5'])
                    tval = datetime.datetime.strptime('2006-01-04', '%Y-%m-%d')
                elif 'J081435' in fname:
                    Rlist.append(Rdicts['HIRES']['C1'])
                    tval = datetime.datetime.strptime('2006-12-26', '%Y-%m-%d') # 2008 too
                elif 'J095309' in fname:
                    Rlist.append(Rdicts['HIRES']['C1'])
                    tval = datetime.datetime.strptime('2005-03-18', '%Y-%m-%d')
                elif 'J113418' in fname:
                    Rlist.append(Rdicts['HIRES']['C5'])
                    tval = datetime.datetime.strptime('2006-01-05', '%Y-%m-%d')
                elif 'J135706' in fname:
                    Rlist.append(Rdicts['HIRES']['C5'])
                    tval = datetime.datetime.strptime('2007-04-28', '%Y-%m-%d')
                elif 'J155556.9' in fname:
                    Rlist.append(Rdicts['HIRES']['C5'])
                    tval = datetime.datetime.strptime('2005-04-15', '%Y-%m-%d')
                elif 'J212329' in fname:
                    Rlist.append(Rdicts['HIRES']['E3'])
                    tval = datetime.datetime.strptime('2006-08-20', '%Y-%m-%d')
                else:
                    pdb.set_trace()
            else:
                tval = datetime.datetime.strptime(head['DATE-OBS'], '%Y-%m-%d')
            dateobslist.append(datetime.datetime.strftime(tval,'%Y-%m-%d'))
        elif 'ESI' in fname:
            instrlist.append('ESI')
            telelist.append('Keck-II')
            gratinglist.append('ECH')
            try:
                Rlist.append(set_resolution(head))
            except ValueError:
                print("Using R=6,000 for ESI")
                Rlist.append(6000.)
            try:
                tval = datetime.datetime.strptime(head['DATE'], '%Y-%m-%d')
            except KeyError:
                if ('J223438.5' in fname) or ('J231543' in fname):
                    tval = datetime.datetime.strptime('2004-09-11', '%Y-%m-%d')
                else:
                    pdb.set_trace()
            dateobslist.append(datetime.datetime.strftime(tval,'%Y-%m-%d'))
        elif 'MIKE' in fname:  # APPROXIMATE
            if 'MIKEr' in fname:
                instrlist.append('MIKEr')
                gratinglist.append('RED')
            elif 'MIKEb' in fname:
                instrlist.append('MIKEb')
                gratinglist.append('BLUE')
            else:
                instrlist.append('MIKE')
                gratinglist.append('BOTH')
            telelist.append('Magellan')
            sep = full_coord[mt[0]].separation(mike_coord)
            imin = np.argmin(sep)
            if sep[imin] > 1.*u.arcsec:
                pdb.set_trace()
                raise ValueError("Bad separation in MIKE")
            # R and Date
            Rlist.append(25000. / mike_meta['Slit'][imin])
            tval = datetime.datetime.strptime(mike_meta['DATE-OBS'][imin], '%Y-%b-%d')
            dateobslist.append(datetime.datetime.strftime(tval,'%Y-%m-%d'))
        elif 'MAGE' in fname:  # APPROXIMATE
            instrlist.append('MagE')
            if 'Clay' in head['TELESCOP']:
                telelist.append('Magellan/Clay')
            else:
                telelist.append('Magellan/Baade')
            gratinglist.append('N/A')
            Rlist.append(set_resolution(head))
            dateobslist.append(head['DATE-OBS'])
        else:  # MagE
            raise ValueError("UH OH")
        # Only way to set the dataset correctly
        if chk_meta_only:
            continue
        spec_set[kk] = data

    # Add columns
    meta = hdlls_full[full_idx]
    nmeta = len(meta)
    meta.add_column(Column([2000.]*nmeta, name='EPOCH'))
    meta.add_column(Column(npixlist, name='NPIX'))
    meta.add_column(Column([str(date) for date in dateobslist], name='DATE-OBS'))
    meta.add_column(Column(wvminlist, name='WV_MIN'))
    meta.add_column(Column(wvmaxlist, name='WV_MAX'))
    meta.add_column(Column(Rlist, name='R'))
    meta.add_column(Column(np.arange(nmeta,dtype=int),name='SURVEY_ID'))
    meta.add_column(Column(gratinglist, name='GRATING'))
    meta.add_column(Column(instrlist, name='INSTR'))
    meta.add_column(Column(telelist, name='TELESCOPE'))
    meta.rename_column('Z_QSO', 'zem')

    # Add HDLLS meta to hdf5
    if chk_meta(meta):
        if chk_meta_only:
            pdb.set_trace()
        hdf[sname]['meta'] = meta
    else:
        raise ValueError("meta file failed")
    # References
    refs = [dict(url='http://adsabs.harvard.edu/abs/2015ApJS..221....2P',
                 bib='prochaska+15'),
            ]
    jrefs = ltu.jsonify(refs)
    hdf[sname]['meta'].attrs['Refs'] = json.dumps(jrefs)
    #
    return


