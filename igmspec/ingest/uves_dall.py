""" Module to ingest UVES data from Dall'Aglio

Dall'Aglio et al. 2008, A&A, 491, 465
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import pdb
import os, glob
import imp
import json

from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.table import Table, Column, vstack
from astropy.time import Time
from astropy import units as u

from linetools.spectra import io as lsio
from linetools import utils as ltu

from specdb.specdb import IgmSpec
from specdb.build.utils import chk_meta
from specdb.build.utils import init_data
from specdb.zem.utils import zem_from_radec

igms_path = imp.find_module('igmspec')[1]


def grab_meta():
    """ Grab UVES Dall'Aglio meta table

    Returns
    -------

    """
    #
    uvesdall_meta = Table.read(os.getenv('RAW_IGMSPEC')+'/UVES_Dall/uves_dall_summ.dat', format='ascii')
    nspec = len(uvesdall_meta)
    # DATE
    #datearr = [day.split('/') for day in list(uvesdall_meta['ObsDate'])]
    #ndate = ['20'+str(day[2])+'-'+str(day[0])+'-'+str(day[1]) for day in datearr]
    t = Time(uvesdall_meta['OBS-DATE'], out_subfmt='date')  # Fixes to YYYY-MM-DD
    uvesdall_meta.add_column(Column(t.iso, name='DATE-OBS'))
    # RA/DEC
    coord = SkyCoord(ra=uvesdall_meta['RA'], dec=uvesdall_meta['DEC'], unit=(u.hour,u.deg))
    rad = [icoord.ra.value for icoord in coord]
    decd = [icoord.dec.value for icoord in coord]
    uvesdall_meta.rename_column('RA', 'RA_STR')
    uvesdall_meta.rename_column('DEC', 'DEC_STR')
    uvesdall_meta['RA_GROUP'] = rad
    uvesdall_meta['DEC_GROUP'] = decd
    # Add zem
    igmsp = IgmSpec()
    ztbl = Table(igmsp.hdf['quasars'].value)
    zem, zsource = zem_from_radec(rad, decd, ztbl)
    badz = np.where(zem < 0.1)[0]
    for ibadz in badz:
        if uvesdall_meta['NAME'][ibadz] == 'HE2243-6031':
            zem[ibadz] = 3.005
            zsource[ibadz] = 'FOP13'  # Fumagalli+13
        elif uvesdall_meta['NAME'][ibadz] == 'HE1341-1020':
            zem[ibadz] = 2.137
            zsource[ibadz] = 'Dall08'  # Dall'Aglio+08
        elif uvesdall_meta['NAME'][ibadz] == 'Q0002-422':
            zem[ibadz] = 2.769
            zsource[ibadz] = 'Dall08'  # Dall'Aglio+08
        elif uvesdall_meta['NAME'][ibadz] == 'PKS2000-330':
            zem[ibadz] = 3.786
            zsource[ibadz] = 'Dall08'  # Dall'Aglio+08
        else:
            raise ValueError("Should not be here")

    uvesdall_meta['zem_GROUP'] = zem
    uvesdall_meta['sig_zem'] = [0.]*nspec
    uvesdall_meta['flag_zem'] = zsource
    #
    uvesdall_meta.add_column(Column([2000.]*nspec, name='EPOCH'))
    uvesdall_meta.add_column(Column(['VLT']*nspec, name='TELESCOPE'))
    uvesdall_meta.add_column(Column(['UVES']*nspec, name='INSTR'))
    uvesdall_meta.add_column(Column(['BOTH']*nspec, name='GRATING'))
    uvesdall_meta.add_column(Column([45000.]*nspec, name='R'))
    uvesdall_meta['STYPE'] = str('QSO')
    # Sort
    uvesdall_meta.sort('RA_GROUP')
    # Check
    assert chk_meta(uvesdall_meta, chk_cat_only=True)
    return uvesdall_meta

'''
def meta_for_build(uvesdall_meta=None):
    """ Generates the meta data needed for the IGMSpec build
    Returns
    -------
    meta : Table
    """
    if uvesdall_meta is None:
        uvesdall_meta = grab_meta()
    nqso = len(uvesdall_meta)
    #
    meta = Table()
    for key in ['RA', 'DEC', 'zem', 'sig_zem', 'flag_zem']:
        meta[key] = uvesdall_meta[key]
    meta['STYPE'] = [str('QSO')]*nqso
    # Return
    return meta
'''


def hdf5_adddata(hdf, sname, meta, debug=False, chk_meta_only=False):
    """ Append UVES_Dall data to the h5 file

    Parameters
    ----------
    hdf : hdf5 pointer
    IDs : ndarray
      int array of IGM_ID values in mainDB
    sname : str
      Survey name
    chk_meta_only : bool, optional
      Only check meta file;  will not write

    Returns
    -------

    """
    from specdb import defs
    # Add Survey
    print("Adding {:s} survey to DB".format(sname))
    uvesdall_grp = hdf.create_group(sname)
    # Load up
    Rdicts = defs.get_res_dicts()
    # Checks
    if sname != 'UVES_Dall':
        raise IOError("Expecting UVES_Dall!!")

    # Build spectra (and parse for meta)
    nspec = len(meta)
    max_npix = 150000  # Just needs to be large enough
    data = init_data(max_npix, include_co=True)
    # Init
    spec_set = hdf[sname].create_dataset('spec', data=data, chunks=True,
                                         maxshape=(None,), compression='gzip')
    spec_set.resize((nspec,))
    wvminlist, wvmaxlist, npixlist, speclist = [], [], [], []
    # Loop
    maxpix = 0
    for jj,row in enumerate(meta):
        # Read
        specfile = os.getenv('RAW_IGMSPEC')+'/UVES_Dall/{:s}_flux.dat'.format(row['NAME'])
        print("UVES_Dall: Reading {:s}".format(specfile))
        spec = Table.read(specfile,format='ascii.fast_no_header',guess=False)#, data_start=1)
        # Parse name
        fname = specfile.split('/')[-1]
        # npix
        npix = len(spec['col1'])
        if npix > max_npix:
            raise ValueError("Not enough pixels in the data... ({:d})".format(npix))
        else:
            maxpix = max(npix,maxpix)
        # Continuum
        # Some fiddling about
        for key in ['wave','flux','sig','co']:
            data[key] = 0.  # Important to init (for compression too)
        data['flux'][0][:npix] = spec['col2']
        data['sig'][0][:npix] = spec['col3']
        data['wave'][0][:npix] = spec['col1']
        data['co'][0][:npix] = spec['col4']
        # Meta
        speclist.append(str(fname))
        wvminlist.append(np.min(data['wave'][0][:npix]))
        wvmaxlist.append(np.max(data['wave'][0][:npix]))
        npixlist.append(npix)
        # Only way to set the dataset correctly
        if chk_meta_only:
            continue
        spec_set[jj] = data

    #
    print("Max pix = {:d}".format(maxpix))
    # Add columns
    meta.add_column(Column(speclist, name='SPEC_FILE'))
    meta.add_column(Column(npixlist, name='NPIX'))
    meta.add_column(Column(wvminlist, name='WV_MIN'))
    meta.add_column(Column(wvmaxlist, name='WV_MAX'))
    meta.add_column(Column(np.arange(nspec,dtype=int),name='GROUP_ID'))

    # Add HDLLS meta to hdf5
    if chk_meta(meta):
        if chk_meta_only:
            pdb.set_trace()
        hdf[sname]['meta'] = meta
    else:
        raise ValueError("meta file failed")

    # References
    refs = [dict(url='http://adsabs.harvard.edu/abs/2008A%26A...491..465D',
                 bib='dallaglio+08'),
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
    ssa_dict = default_fields(flux='flambda')
    ssa_dict['FluxCalib']='RELATIVE'
    ssa_dict['Title'] = 'Dall''Aglio et al. (2008) compilation of VLT/UVES spectra'
    hdf[dset]['meta'].attrs['SSA'] = json.dumps(ltu.jsonify(ssa_dict))
