""" Module to ingest GGG Survey data

Worseck et al. 2014
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import pdb
import os
import imp
import json

from astropy.table import Table, Column, vstack
from astropy.time import Time

from linetools.spectra import io as lsio
from linetools import utils as ltu

from specdb.build.utils import chk_meta
from specdb.build.utils import init_data

igms_path = imp.find_module('igmspec')[1]


def grab_meta():
    """ Grab GGG meta Table
    Returns
    -------

    """
    # This table has units in it!
    meta = Table.read(os.getenv('RAW_IGMSPEC')+'/GGG/GGG_catalog.fits.gz')
    nqso = len(meta)
    # Turn off RA/DEC units
    for key in ['RA', 'DEC']:
        meta[key].unit = None
    meta.rename_column('RA', 'RA_GROUP')
    meta.rename_column('DEC', 'DEC_GROUP')
    #
    # Add zem
    meta['zem_GROUP'] = meta['z_gmos']
    meta['sig_zem'] = meta['zerror_gmos']
    meta['flag_zem'] = [str('GGG')]*nqso
    meta.add_column(Column([2000.]*nqso, name='EPOCH'))
    #
    meta['STYPE'] = [str('QSO')]*nqso
    # Double up for the two gratings
    ggg_meta = vstack([meta,meta])
    # Check
    assert chk_meta(ggg_meta, chk_cat_only=True)
    # Return
    return ggg_meta

'''
def meta_for_build():
    """ Generates the meta data needed for the IGMSpec build
    Returns
    -------
    meta : Table
    """
    ggg_meta = grab_meta()
    # Cut down to unique QSOs
    names = np.array([name[0:26] for name in ggg_meta['SDSSJ']])
    uni, uni_idx = np.unique(names, return_index=True)
    ggg_meta = ggg_meta[uni_idx]
    nqso = len(ggg_meta)
    #
    meta = Table()
    for key in ['RA', 'DEC', 'zem', 'sig_zem', 'flag_zem']:
        meta[key] = ggg_meta[key]
    meta['STYPE'] = [str('QSO')]*nqso
    # Return
    return meta
'''


def hdf5_adddata(hdf, sname, meta, debug=False, chk_meta_only=False):
    """ Append GGG data to the h5 file

    Parameters
    ----------
    hdf : hdf5 pointer
    sname : str
      Survey name
    meta : Table
    chk_meta_only : bool, optional
      Only check meta file;  will not write

    Returns
    -------

    """
    # Add Survey
    print("Adding {:s} survey to DB".format(sname))
    ggg_grp = hdf.create_group(sname)
    # Load up
    if sname != 'GGG':
        raise IOError("Not expecting this survey..")

    # Build spectra (and parse for meta)
    nspec = len(meta)
    max_npix = 1600  # Just needs to be large enough
    # Init
    data = init_data(max_npix, include_co=False)
    spec_set = hdf[sname].create_dataset('spec', data=data, chunks=True,
                                         maxshape=(None,), compression='gzip')
    spec_set.resize((nspec,))
    Rlist = []
    wvminlist = []
    wvmaxlist = []
    npixlist = []
    speclist = []
    gratinglist = []
    telelist = []
    dateobslist = []
    instrlist = []
    # Loop
    path = os.getenv('RAW_IGMSPEC')+'/GGG/'
    maxpix = 0
    for jj,row in enumerate(meta):
        # Generate full file
        if jj >= nspec//2:
            full_file = path+row['name']+'_R400.fits.gz'
            gratinglist.append('R400')
        else:
            full_file = path+row['name']+'_B600.fits.gz'
            gratinglist.append('B600')
        # Extract
        print("GGG: Reading {:s}".format(full_file))
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
        head = spec.header
        speclist.append(str(fname))
        wvminlist.append(np.min(data['wave'][0][:npix]))
        wvmaxlist.append(np.max(data['wave'][0][:npix]))
        telelist.append(head['OBSERVAT'])
        instrlist.append(head['INSTRUME'])
        tval = Time(head['DATE'], format='isot', out_subfmt='date')
        dateobslist.append(tval.iso)
        npixlist.append(npix)
        if 'R400' in fname:
            Rlist.append(833.)
        else:
            Rlist.append(940.)
        # Only way to set the dataset correctly
        if chk_meta_only:
            continue
        spec_set[jj] = data

    #
    print("Max pix = {:d}".format(maxpix))
    # Add columns
    meta.add_column(Column(speclist, name='SPEC_FILE'))
    meta.add_column(Column(gratinglist, name='GRATING'))
    meta.add_column(Column(telelist, name='TELESCOPE'))
    meta.add_column(Column(instrlist, name='INSTR'))
    meta.add_column(Column(dateobslist, name='DATE-OBS'))
    meta.add_column(Column(npixlist, name='NPIX'))
    meta.add_column(Column(wvminlist, name='WV_MIN'))
    meta.add_column(Column(wvmaxlist, name='WV_MAX'))
    meta.add_column(Column(Rlist, name='R'))
    meta.add_column(Column(np.arange(nspec,dtype=int),name='GROUP_ID'))

    # Add HDLLS meta to hdf5
    if chk_meta(meta):
        if chk_meta_only:
            pdb.set_trace()
        hdf[sname]['meta'] = meta
    else:
        raise ValueError("meta file failed")

    pdb.set_trace()
    # References
    refs = [dict(url='http://adsabs.harvard.edu/abs/2014MNRAS.445.1745W',
                 bib='worseck+14')]
    jrefs = ltu.jsonify(refs)
    hdf[sname]['meta'].attrs['Refs'] = json.dumps(jrefs)
    #
    return


