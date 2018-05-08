""" Module to ingest GGG Survey data

Worseck et al. 2014
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import pdb
import os
import json

from astropy.table import Table, Column, vstack
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u

from linetools.spectra import io as lsio
from linetools import utils as ltu

from specdb.build.utils import chk_meta
from specdb.build.utils import init_data


path = os.getenv('RAW_IGMSPEC')+'/ESI_z6/'

def grab_meta():
    """ Grab GGG meta Table
    Returns
    -------

    """
    # This table has units in it!
    meta = Table.read(os.getenv('RAW_IGMSPEC')+'/ESI_z6/overview_data_igmspec.txt', format='ascii', delimiter='\t')
    nqso = len(meta)
    # Rename
    meta.rename_column('RAdeg', 'RA_GROUP')
    meta.rename_column('DECdeg', 'DEC_GROUP')
    meta.rename_column('z', 'zem_GROUP')
    meta.rename_column('instrument', 'INSTR')
    meta.rename_column('telescope', 'TELESCOPE')
    meta.rename_column('date', 'DATE-OBS')
    #
    # Add zem
    meta['sig_zem'] = 0.
    meta['flag_zem'] = str('ESI_z6')
    meta.add_column(Column([2000.]*nqso, name='EPOCH'))
    #
    meta['STYPE'] = str('QSO')
    meta['DISPERSER'] = str('ECH')
    # Check
    assert chk_meta(meta, chk_cat_only=True)
    # Return
    return meta


def read_spec(row):
    # Filename
    coord = SkyCoord(ra=row['RA_GROUP'], dec=row['DEC_GROUP'], unit='deg')
    filename = 'J{:s}{:s}.txt'.format(coord.ra.to_string(unit=u.hour,sep='',pad=True)[0:4],
                              coord.dec.to_string(sep='',pad=True,alwayssign=True)[0:5])
    # Read
    spec = lsio.readspec(path+filename)
    # Return
    return filename, spec

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
    if sname != 'ESI_z6':
        raise IOError("Not expecting this survey..")

    # Build spectra (and parse for meta)
    nspec = len(meta)
    max_npix = 30000  # Just needs to be large enough
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
    maxpix = 0
    for jj,row in enumerate(meta):
        # Generate full file
        full_file, spec = read_spec(row)
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
        #head = spec.header
        speclist.append(str(fname))
        wvminlist.append(np.min(data['wave'][0][:npix]))
        wvmaxlist.append(np.max(data['wave'][0][:npix]))
        #telelist.append(head['OBSERVAT'])
        #instrlist.append(head['INSTRUME'])
        #tval = Time(head['DATE'], format='isot', out_subfmt='date')
        #dateobslist.append(tval.iso)
        npixlist.append(npix)
        # Only way to set the dataset correctly
        if chk_meta_only:
            continue
        spec_set[jj] = data

    #
    print("Max pix = {:d}".format(maxpix))
    # Add columns
    meta.add_column(Column(speclist, name='SPEC_FILE'))
    #meta.add_column(Column(telelist, name='TELESCOPE'))
    #meta.add_column(Column(instrlist, name='INSTR'))
    #meta.add_column(Column(dateobslist, name='DATE-OBS'))
    meta.add_column(Column(npixlist, name='NPIX'))
    meta.add_column(Column(wvminlist, name='WV_MIN'))
    meta.add_column(Column(wvmaxlist, name='WV_MAX'))
    #meta.add_column(Column(Rlist, name='R'))
    meta.add_column(Column(np.arange(nspec,dtype=int),name='GROUP_ID'))

    # Add HDLLS meta to hdf5
    if chk_meta(meta):
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


def add_ssa(hdf, dset):
    """  Add SSA info to meta dataset

    Parameters
    ----------
    hdf
    dset : str
    """
    from specdb.ssa import default_fields
    Title = '{:s}: Giant Gemini GMOS Survey of z>4 quasars'.format(dset)
    ssa_dict = default_fields(Title, flux='flambda')
    hdf[dset]['meta'].attrs['SSA'] = json.dumps(ltu.jsonify(ssa_dict))
