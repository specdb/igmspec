""" Module to ingest High z ESI DLA

Rafelski et al. 2012, 2014
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
from astropy.io import fits
from astropy import units as u

from linetools.spectra import io as lsio
from linetools import utils as ltu

from specdb.build.utils import chk_meta

igms_path = imp.find_module('igmspec')[1]


def grab_meta():
    """ Grab High-z ESI meta Table

    Returns
    -------

    """
    #
    esidla_meta = Table.read(os.getenv('RAW_IGMSPEC')+'/HighzESIDLA/ascii_highz_rafelski.list', format='ascii')
    nspec = len(esidla_meta)
    # DATE
    datearr = [day.split('/') for day in list(esidla_meta['ObsDate'])]
    ndate = ['20'+str(day[2])+'-'+str(day[0])+'-'+str(day[1]) for day in datearr]
    t = Time(ndate, out_subfmt='date')  # Fixes to YYYY-MM-DD
    esidla_meta.add_column(Column(t.iso, name='DATE-OBS'))
    # Add zem
    esidla_meta['sig_zem'] = [0.]*nspec
    esidla_meta['flag_zem'] = [str('SDSS')]*nspec
    #
    esidla_meta.add_column(Column([2000.]*nspec, name='EPOCH'))
    esidla_meta.add_column(Column(['KeckII']*nspec, name='TELESCOPE'))
    esidla_meta.add_column(Column(['ESI']*nspec, name='INSTR'))
    esidla_meta.add_column(Column(['ECH']*nspec, name='GRATING'))
    # Rename
    esidla_meta.rename_column('RA', 'RA_GROUP')
    esidla_meta.rename_column('DEC', 'DEC_GROUP')
    esidla_meta.rename_column('zem', 'zem_GROUP')
    esidla_meta['STYPE'] = str('QSO')
    # Sort
    esidla_meta.sort('RA_GROUP')
    # Check
    assert chk_meta(esidla_meta, chk_cat_only=True)
    #
    return esidla_meta

'''
def meta_for_build(esidla_meta=None):
    """ Generates the meta data needed for the IGMSpec build
    Returns
    -------
    meta : Table
    """
    if esidla_meta is None:
        esidla_meta = grab_meta()
    nqso = len(esidla_meta)
    #
    meta = Table()
    for key in ['RA', 'DEC', 'zem', 'sig_zem', 'flag_zem']:
        meta[key] = esidla_meta[key]
    meta['STYPE'] = [str('QSO')]*nqso
    # Return
    return meta
'''


def hdf5_adddata(hdf, IDs, sname, debug=False, chk_meta_only=False):
    """ Append ESI data to the h5 file

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
    esidla_grp = hdf.create_group(sname)
    # Load up
    Rdicts = defs.get_res_dicts()
    meta = grab_meta()
    bmeta = meta_for_build()
    if len(meta) != len(bmeta):
        raise ValueError("Should be the same size")
    # Checks
    if sname != 'ESI_DLA':
        raise IOError("Expecting ESI_DLA!!")
    if np.sum(IDs < 0) > 0:
        raise ValueError("Bad ID values")
    # Open Meta tables
    if len(bmeta) != len(IDs):
        raise ValueError("Wrong sized table..")

    # Generate ID array from RA/DEC
    c_cut = SkyCoord(ra=bmeta['RA'], dec=bmeta['DEC'], unit='deg')
    c_all = SkyCoord(ra=meta['RA'], dec=meta['DEC'], unit='deg')
    # Find new sources
    idx, d2d, d3d = match_coordinates_sky(c_all, c_cut, nthneighbor=1)
    if np.sum(d2d > 0.1*u.arcsec):
        raise ValueError("Bad matches in ESI_DLA")
    meta_IDs = IDs[idx]
    meta.add_column(Column(meta_IDs, name='IGM_ID'))

    # Build spectra (and parse for meta)
    nspec = len(meta)
    max_npix = 50000  # Just needs to be large enough
    data = np.ma.empty((1,),
                       dtype=[(str('wave'), 'float64', (max_npix)),
                              (str('flux'), 'float32', (max_npix)),
                              (str('sig'),  'float32', (max_npix)),
                              #(str('co'),   'float32', (max_npix)),
                             ])
    # Init
    spec_set = hdf[sname].create_dataset('spec', data=data, chunks=True,
                                         maxshape=(None,), compression='gzip')
    spec_set.resize((nspec,))
    Rlist = []
    wvminlist = []
    wvmaxlist = []
    npixlist = []
    speclist = []
    # Loop
    maxpix = 0
    for jj,row in enumerate(meta):
        #
        specfile = os.getenv('RAW_IGMSPEC')+'/HighzESIDLA/{:s}a_xF.fits'.format(
            row['Name'])
        print("ESI_DLA: Reading {:s}".format(specfile))
        spec = lsio.readspec(specfile)
        # Parse name
        fname = specfile.split('/')[-1]
        # npix
        npix = spec.npix
        if npix > max_npix:
            raise ValueError("Not enough pixels in the data... ({:d})".format(npix))
        else:
            maxpix = max(npix,maxpix)
        # Continuum
        # Some fiddling about
        for key in ['wave','flux','sig']:
            data[key] = 0.  # Important to init (for compression too)
        data['flux'][0][:npix] = spec.flux.value
        data['sig'][0][:npix] = spec.sig.value
        data['wave'][0][:npix] = spec.wavelength.to('AA').value
        #data['co'][0][:npix] = spec.co.value
        # Meta
        head = spec.header
        speclist.append(str(fname))
        wvminlist.append(np.min(data['wave'][0][:npix]))
        wvmaxlist.append(np.max(data['wave'][0][:npix]))
        npixlist.append(npix)
        try:
            Rlist.append(Rdicts['ESI'][head['SLMSKNAM']])
        except KeyError:
            if row['Slit'] == 0.75:
                Rlist.append(Rdicts['ESI']['0.75_arcsec'])
            elif row['Slit'] == 0.5:
                Rlist.append(Rdicts['ESI']['0.50_arcsec'])
            else:
                pdb.set_trace()
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
    meta.add_column(Column(Rlist, name='R'))
    meta.add_column(Column(np.arange(nspec,dtype=int),name='SURVEY_ID'))

    # Add HDLLS meta to hdf5
    if chk_meta(meta):
        if chk_meta_only:
            pdb.set_trace()
        hdf[sname]['meta'] = meta
    else:
        raise ValueError("meta file failed")

    # References
    refs = [dict(url='http://adsabs.harvard.edu/abs/2012ApJ...755...89R',
                 bib='rafelski+12'),
            dict(url='http://adsabs.harvard.edu/abs/2014ApJ...782L..29R',
                         bib='rafelski+14'),
            ]
    jrefs = ltu.jsonify(refs)
    hdf[sname]['meta'].attrs['Refs'] = json.dumps(jrefs)
    #
    return


