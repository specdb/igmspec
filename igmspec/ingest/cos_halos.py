""" Module to ingest COS-Halos

Tumlinson et al. 2013
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import pdb
import warnings
import os, json, glob, imp

from astropy.table import Table, Column, vstack
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
from astropy.time import Time

from linetools.spectra import io as lsio
from linetools import utils as ltu

from specdb.build.utils import chk_meta
from specdb.build.utils import init_data
#igms_path = imp.find_module('igmspec')[1]


def grab_meta():
    """ Grab COS-Halos meta table
    Returns
    -------

    """
    from time import strptime
    from specdb.zem.utils import zem_from_radec
    from specdb.specdb import IgmSpec
    from specdb.defs import get_res_dicts
    Rdicts = get_res_dicts()
    igmsp = IgmSpec(db_file=os.getenv('SPECDB')+'/IGMspec_DB_v01.hdf5', skip_test=True)

    summ_file = os.getenv('RAW_IGMSPEC')+'/COS-Halos/cos_halos_obs.ascii'
    chalos_meta = Table.read(summ_file, format='ascii')
    # RA/DEC, DATE
    # Visits from this page: http://www.stsci.edu/cgi-bin/get-visit-status?id=11598&markupFormat=html
    visit_file = os.getenv('RAW_IGMSPEC')+'/COS-Halos/cos_halos_visits.ascii'
    ch_visits = Table.read(visit_file,format='ascii')
    ra = []
    dec = []
    datet = []
    for row in chalos_meta:
        coord = ltu.radec_to_coord(row['QSO'])
        ra.append(coord.ra.value)
        dec.append(coord.dec.value)
        #
        visit = row['Visit']
        mtv = np.where(ch_visits['Visit'] == visit)[0]
        if len(mtv) != 1:
            pdb.set_trace()
        else:
            chv = ch_visits['Start_UT'][mtv].data[0]
        icmma = chv.find(',')
        datet.append('{:s}-{:02d}-{:02d}'.format(
                chv[icmma+1:icmma+5], strptime(chv[:3],'%b').tm_mon,
                int(chv[3:icmma])))
    chalos_meta.add_column(Column(ra, name='RA'))
    chalos_meta.add_column(Column(dec, name='DEC'))
    chalos_meta.add_column(Column(datet, name='DATE-OBS'))
    # Others
    chalos_meta.add_column(Column(['      ']*len(chalos_meta), name='TELESCOPE')) # Padding
    chalos_meta.add_column(Column(['     ']*len(chalos_meta), name='INSTR')) # Padding for HIRES
    chalos_meta.add_column(Column(['G130M/G160M']*len(chalos_meta), name='DISPERSER'))
    chalos_meta.add_column(Column([20000.]*len(chalos_meta), name='R'))
    chalos_meta.add_column(Column([2000.]*len(chalos_meta), name='EPOCH'))
    chalos_meta['INSTR'] = 'COS' # Deals with padding
    chalos_meta['TELESCOPE'] = 'HST'
    # Myers for zem
    zem, zsource = zem_from_radec(chalos_meta['RA'], chalos_meta['DEC'], Table(igmsp.hdf['quasars'].value))
    badz = zem <= 0.
    if np.sum(badz) > 0:
        raise ValueError("Bad zem in COS-Halos")
    chalos_meta['zem'] = zem
    chalos_meta['sig_zem'] = 0.  # Need to add
    chalos_meta['flag_zem'] = zsource
    # HIRES
    hires_files = glob.glob(os.getenv('RAW_IGMSPEC')+'/COS-Halos/HIRES/J*f.fits.gz')
    hires_tab = chalos_meta[0:0]
    subnm = np.array([row['QSO'][4:9] for row in chalos_meta])
    signs = np.array([row['QSO'][14] for row in chalos_meta])
    for ifile in hires_files:
        print(ifile)
        fname = ifile.split('/')[-1]
        mt = np.where((subnm == fname[0:5]) & (signs == fname[5]))[0]
        if len(mt) != 1:
            pdb.set_trace()
        # Add row
        hires_tab.add_row(chalos_meta[mt[0]])
        hires_tab[-1]['INSTR'] = 'HIRES'
        hires_tab[-1]['TELESCOPE'] = 'Keck I'
        hires_tab[-1]['DISPERSER'] = 'Red'
        hires_tab[-1]['R'] = Rdicts['HIRES']['C1']
    # Combine
    chalos_meta = vstack([chalos_meta, hires_tab])
    chalos_meta['STYPE'] = str('QSO')
    # Rename
    chalos_meta.rename_column('RA', 'RA_GROUP')
    chalos_meta.rename_column('DEC', 'DEC_GROUP')
    chalos_meta.rename_column('zem', 'zem_GROUP')
    # Check
    assert chk_meta(chalos_meta, chk_cat_only=True)
    # Done
    return chalos_meta


def hdf5_adddata(hdf, sname, meta, debug=False, chk_meta_only=False,
                 mk_test_file=False):
    """ Append COS-Halos data to the h5 file

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
    # Add Survey
    print("Adding {:s} survey to DB".format(sname))
    chalos_grp = hdf.create_group(sname)
    # Load up
    # Checks
    if sname != 'COS-Halos':
        raise IOError("Not expecting this survey..")

    # Build spectra (and parse for meta)
    nspec = len(meta)
    max_npix = 160000  # Just needs to be large enough
    data = init_data(max_npix, include_co=False)
    # Init
    spec_set = hdf[sname].create_dataset('spec', data=data, chunks=True,
                                         maxshape=(None,), compression='gzip')
    spec_set.resize((nspec,))
    wvminlist = []
    wvmaxlist = []
    npixlist = []
    speclist = []
    # Loop
    path = os.getenv('RAW_IGMSPEC')+'/COS-Halos/'
    maxpix = 0
    for jj,row in enumerate(meta):
        # Generate full file
        coord = ltu.radec_to_coord((row['RA_GROUP'],row['DEC_GROUP']))
        if row['INSTR'].strip() == 'COS':
            full_file = path+'/J{:s}{:s}_nbin3_coadd.fits.gz'.format(coord.ra.to_string(unit=u.hour,sep='',pad=True)[0:4],
                                               coord.dec.to_string(sep='',pad=True,alwayssign=True)[0:5])
        else: # HIRES
            full_file = path+'/HIRES/J{:s}{:s}_f.fits.gz'.format(coord.ra.to_string(unit=u.hour,sep='',pad=True)[0:4], coord.dec.to_string(sep='',pad=True,alwayssign=True)[0:5])
        # Extract
        print("COS-Halos: Reading {:s}".format(full_file))
        spec = lsio.readspec(full_file)
        # Parse name
        fname = full_file.split('/')[-1]
        # npix
        npix = spec.npix
        if npix > max_npix:
            raise ValueError("Not enough pixels in the data... ({:d})".format(npix))
        else:
            maxpix = max(npix,maxpix)
        # Some fiddling about
        for key in ['wave','flux','sig']:
            data[key] = 0.  # Important to init (for compression too)
        data['flux'][0][:npix] = spec.flux.value
        data['sig'][0][:npix] = spec.sig.value
        data['wave'][0][:npix] = spec.wavelength.value
        # Meta
        speclist.append(str(fname))
        wvminlist.append(np.min(data['wave'][0][:npix]))
        wvmaxlist.append(np.max(data['wave'][0][:npix]))
        npixlist.append(npix)
        if chk_meta_only:
            continue
        # Only way to set the dataset correctly
        spec_set[jj] = data

    #
    print("Max pix = {:d}".format(maxpix))
    # Add columns
    meta.add_column(Column(speclist, name='SPEC_FILE'))
    meta.add_column(Column(npixlist, name='NPIX'))
    meta.add_column(Column(wvminlist, name='WV_MIN'))
    meta.add_column(Column(wvmaxlist, name='WV_MAX'))
    meta.add_column(Column(np.arange(nspec,dtype=int), name='GROUP_ID'))

    # Add HDLLS meta to hdf5
    if chk_meta(meta):
        if chk_meta_only:
            pdb.set_trace()
        hdf[sname]['meta'] = meta
    else:
        raise ValueError("meta file failed")
    # References
    refs = [dict(url='http://adsabs.harvard.edu/abs/2013ApJ...777...59T',
                 bib='tumlinson+13'),
            dict(url='http://adsabs.harvard.edu/abs/2013ApJS..204...17W',
                         bib='werk+13')
            ]
    jrefs = ltu.jsonify(refs)
    hdf[sname]['meta'].attrs['Refs'] = json.dumps(jrefs)
    #
    return


def add_ssa(hdf, dset):
    """  Add SSA info to meta dataset
    Parameters
    ----------
    hdf
    dset : str
    """
    from specdb.ssa import default_fields
    Title = '{:s}: Quasar Spectra from the COS-Halos Survey'.format(dset)
    ssa_dict = default_fields(Title, flux='flambda', fxcalib='ABSOLUTE')
    hdf[dset]['meta'].attrs['SSA_COS'] = json.dumps(ltu.jsonify(ssa_dict))
    # HIRES
    ssa_dict = default_fields(Title, flux='normalized')
    hdf[dset]['meta'].attrs['SSA_HIRES'] = json.dumps(ltu.jsonify(ssa_dict))
