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

from igmspec.ingest import utils as iiu

igms_path = imp.find_module('igmspec')[1]


def grab_meta():
    """ Grab GGG meta Table
    Returns
    -------

    """
    # This table has units in it!
    ggg_meta = Table.read(os.getenv('RAW_IGMSPEC')+'/GGG/GGG_catalog.fits.gz')
    nqso = len(ggg_meta)
    # Turn off RA/DEC units
    for key in ['RA', 'DEC']:
        ggg_meta[key].unit = None
    #
    # Add zem
    ggg_meta['zem'] = ggg_meta['z_gmos']
    ggg_meta['sig_zem'] = ggg_meta['zerror_gmos']
    ggg_meta['flag_zem'] = [str('GGG')]*nqso
    ggg_meta.add_column(Column([2000.]*nqso, name='EPOCH'))
    #
    return ggg_meta


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
    # Return
    return meta


def hdf5_adddata(hdf, IDs, sname, debug=False, chk_meta_only=False):
    """ Append GGG data to the h5 file

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
    # Add Survey
    print("Adding {:s} survey to DB".format(sname))
    ggg_grp = hdf.create_group(sname)
    # Load up
    meta = grab_meta()
    bmeta = meta_for_build()
    # Checks
    if sname != 'GGG':
        raise IOError("Not expecting this survey..")
    if np.sum(IDs < 0) > 0:
        raise ValueError("Bad ID values")
    # Open Meta tables
    if len(bmeta) != len(IDs):
        raise ValueError("Wrong sized table..")

    # Generate ID array from RA/DEC
    meta_IDs = IDs
    meta.add_column(Column(meta_IDs, name='IGM_ID'))


    # Double up for the two gratings
    meta = vstack([meta,meta])

    # Build spectra (and parse for meta)
    nspec = len(meta)
    max_npix = 1600  # Just needs to be large enough
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
            gratinglist.append('B400')
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
    meta.add_column(Column(np.arange(nspec,dtype=int),name='SURVEY_ID'))

    # Add HDLLS meta to hdf5
    if iiu.chk_meta(meta):
        if chk_meta_only:
            pdb.set_trace()
        hdf[sname]['meta'] = meta
    else:
        raise ValueError("meta file failed")

    # References
    refs = [dict(url='http://adsabs.harvard.edu/abs/2014MNRAS.445.1745W',
                 bib='worseck+14')]
    jrefs = ltu.jsonify(refs)
    hdf[sname]['meta'].attrs['Refs'] = json.dumps(jrefs)
    #
    return


