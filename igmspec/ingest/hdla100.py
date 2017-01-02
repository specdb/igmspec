""" Module to ingest HIRES DLA 100 Survey data

Neeleman et al. 2013
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import pdb
import os, json, glob, imp
import datetime

from astropy.table import Table, Column
from astropy import units as u
from astropy.time import Time

from linetools import utils as ltu
from linetools.spectra import io as lsio

from specdb.build.utils import chk_meta
from specdb.build.utils import init_data
from specdb.build.utils import set_resolution

igms_path = imp.find_module('igmspec')[1]


def grab_meta():
    """ Generates the meta data needed for the IGMSpec build
    Returns
    -------
    meta : Table
    spec_files : list
      List of spec_file names
    """
    # Load DLA
    from pyigm.surveys.dlasurvey import DLASurvey
    hdla100 = DLASurvey.neeleman13_tree()
    # Cut down to unique QSOs
    spec_files = []
    names = []
    ra = []
    dec = []
    coords = hdla100.coord
    cnt = 0
    for coord in coords:
        # Load
        names.append('J{:s}{:s}'.format(coord.ra.to_string(unit=u.hour, sep='', pad=True, precision=2),
                                       coord.dec.to_string(sep='', pad=True, precision=1)))
        # RA/DEC
        ra.append(coord.ra.value)
        dec.append(coord.dec.value)
        # SPEC_FILE
        fname = hdla100._abs_sys[cnt]._datdict['hi res file'].split('/')[-1]
        spec_files.append(fname)
        cnt += 1
    uni, uni_idx = np.unique(names, return_index=True)
    nqso = len(uni_idx)
    #
    meta = Table()
    meta['RA_GROUP'] = np.array(ra)[uni_idx]
    meta['DEC_GROUP'] = np.array(dec)[uni_idx]
    meta['zem_GROUP'] = hdla100.zem[uni_idx]
    meta['sig_zem'] = [0.]*nqso
    meta['flag_zem'] = [str('UNKN')]*nqso
    meta['STYPE'] = [str('QSO')]*nqso
    meta['SPEC_FILE'] = np.array(spec_files)[uni_idx]
    # Check
    assert chk_meta(meta, chk_cat_only=True)
    return meta


def hdf5_adddata(hdf, sname, hdla100_meta, debug=False, chk_meta_only=False,
                 mk_test_file=False):
    """ Append HDLA100 data to the h5 file

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
    # Checks
    if sname != 'HDLA100':
        raise IOError("Not expecting this survey..")

    # Build spectra (and parse for meta)
    #if mk_test_file:
    #    hdla100_full = hdlls_full[0:3]
    max_npix = 192000  # Just needs to be large enough
    data = init_data(max_npix, include_co=False)
    # Init
    spec_set = hdf[sname].create_dataset('spec', data=data, chunks=True,
                                         maxshape=(None,), compression='gzip')
    nspec = len(hdla100_meta)
    spec_set.resize((nspec,))
    Rlist = []
    wvminlist = []
    wvmaxlist = []
    dateobslist = []
    npixlist = []
    gratinglist = []
    # Loop
    for jj,row in enumerate(hdla100_meta):
        kk = jj
        # Extract
        f = os.getenv('RAW_IGMSPEC')+'/HDLA100/'+row['SPEC_FILE']
        spec = lsio.readspec(f)
        # Parse name
        fname = f.split('/')[-1]
        # npix
        head = spec.header
        npix = spec.npix
        if npix > max_npix:
            raise ValueError("Not enough pixels in the data... ({:d})".format(npix))
        # Some fiddling about
        for key in ['wave','flux','sig']:
            data[key] = 0.  # Important to init (for compression too)
        data['flux'][0][:npix] = spec.flux.value
        data['sig'][0][:npix] = spec.sig.value
        data['wave'][0][:npix] = spec.wavelength.value
        # Meta
        wvminlist.append(np.min(data['wave'][0][:npix]))
        wvmaxlist.append(np.max(data['wave'][0][:npix]))
        npixlist.append(npix)
        try:
            Rlist.append(set_resolution(head))
        except ValueError:
            raise ValueError("Header is required for {:s}".format(fname))
        else:
            if '/' in head['DATE-OBS']:
                spl = head['DATE-OBS'].split('/')
                t = Time(datetime.datetime(int(spl[2])+1900, int(spl[1]), int(spl[0])), format='datetime')
            else:
                t = Time(head['DATE-OBS'], format='isot', out_subfmt='date')
        dateobslist.append(t.iso)
        # Grating
        try:
            gratinglist.append(head['XDISPERS'])
        except KeyError:
            try:
                yr = t.value.year
            except AttributeError:
                yr = int(t.value[0:4])
            if yr <= 1997:
                gratinglist.append('RED')
            else:
                pdb.set_trace()
        # Only way to set the dataset correctly
        if chk_meta_only:
            continue
        spec_set[kk] = data

    # Add columns
    nmeta = len(hdla100_meta)
    hdla100_meta.add_column(Column([2000.]*nmeta, name='EPOCH'))
    hdla100_meta.add_column(Column(npixlist, name='NPIX'))
    hdla100_meta.add_column(Column([str(date) for date in dateobslist], name='DATE-OBS'))
    hdla100_meta.add_column(Column(wvminlist, name='WV_MIN'))
    hdla100_meta.add_column(Column(wvmaxlist, name='WV_MAX'))
    hdla100_meta.add_column(Column(Rlist, name='R'))
    hdla100_meta.add_column(Column(np.arange(nmeta,dtype=int),name='GROUP_ID'))
    hdla100_meta.add_column(Column(gratinglist, name='GRATING'))
    hdla100_meta['INSTR'] = ['HIRES']*nspec
    hdla100_meta['TELESCOPE'] = ['Keck-I']*nspec
    #hdla100_meta.rename_column('Z_QSO', 'zem')

    # Add HDLLS meta to hdf5
    if chk_meta(hdla100_meta):
        if chk_meta_only:
            pdb.set_trace()
        hdf[sname]['meta'] = hdla100_meta
    else:
        raise ValueError("meta file failed")
    # References
    refs = [dict(url='http://adsabs.harvard.edu/abs/2013ApJ...769...54N',
                 bib='neeleman+13'),
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
    Title = '{:s}: The Keck/HIRES Survey of 100 Damped Lya Systems'.format(dset)
    ssa_dict = default_fields(Title, flux='normalized')
    hdf[dset]['meta'].attrs['SSA'] = json.dumps(ltu.jsonify(ssa_dict))
