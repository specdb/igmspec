""" Module to ingest HST+FUSE AGN spectra

Cooksey et al. 2010
"""
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import pdb
import warnings
import os, json

from astropy.table import Table, Column
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy import units as u
from astropy.io.fits import Header
from astropy.io import fits
from astropy.time import Time

from linetools.spectra import io as lsio
from linetools import utils as ltu

#from igmspec.ingest import utils as iiu
from specdb.build.utils import chk_meta
from specdb.build.utils import init_data
from specdb import defs

#igms_path = imp.find_module('igmspec')[1]


def grab_meta():
    """ Grab HST/FUSE Cooksey meta Table
    Returns
    -------

    """
    hstc_file = os.getenv('RAW_IGMSPEC')+'/HST_Cooksey/HSTQSO_pre-SM4.lst'
    hstc_meta = Table.read(hstc_file, format='ascii')
    # Cutting those without proper header (for now)
    #badf = ['PKS2005-489lif1a.fits', 'NGC7469lif2a.fits', 'NGC7469sic2b.fits',
    #        'NGC7469lif2b.fits', 'NGC7469sic2a.fits',
    #        'AKN564lif1a.fits', 'AKN564lif1b.fits', 'AKN564lif2a.fits', 'AKN564lif2b.fits',
    #        'AKN564sic1a.fits', 'AKN564sic1b.fits', 'AKN564sic2a.fits', 'AKN564sic2b.fits']
    badf = []
    gdm = np.array([True]*len(hstc_meta))
    for ibadf in badf:
        mt = np.where(hstc_meta['SPEC_FILE'] == ibadf)[0]
        gdm[mt] = False
    for jj, row in enumerate(hstc_meta):  # Skip continua
        if '_c.fits' in row['SPEC_FILE']:
            gdm[jj] = False
        if '_E.fits' in row['SPEC_FILE']:
            gdm[jj] = False
        #if row['INSTR'] == 'GHRS':
        #    gdm[jj] = False
    hstc_meta = hstc_meta[gdm]
    gdf = hstc_meta['INSTR'] == 'FUSE'
    #hstc_meta = hstc_meta[gdf]
    hstc_meta['TELESCOPE'] = 'FUSE'
    hstc_meta[~gdf]['TELESCOPE'] = 'HST'
    #
    hstc_meta.add_column(Column([2000.]*len(hstc_meta), name='EPOCH'))
    hstc_meta['sig_zem'] = 0.
    hstc_meta['flag_zem'] = str('UNKWN')
    hstc_meta['STYPE'] = str('QSO')
    # RENAME
    hstc_meta.rename_column('RA', 'RA_GROUP')
    hstc_meta.rename_column('DEC', 'DEC_GROUP')
    hstc_meta.rename_column('zem', 'zem_GROUP')
    # Check
    assert chk_meta(hstc_meta, chk_cat_only=True)
    return hstc_meta

'''
def meta_for_build():
    """ Generates the meta data needed for the IGMSpec build
    Returns
    -------
    meta : Table
    """
    # Cut down to unique QSOs
    hstc_meta = grab_meta()
    names = hstc_meta['QSO'].data
    uni, uni_idx = np.unique(names, return_index=True)
    hstc_meta = hstc_meta[uni_idx]
    nqso = len(hstc_meta)
    #
    meta = Table()
    meta['RA'] = hstc_meta['RA']
    meta['DEC'] = hstc_meta['DEC']
    meta['zem'] = hstc_meta['zem']
    meta['sig_zem'] = [0.]*nqso
    meta['flag_zem'] = [str('UNKWN')]*nqso
    meta['STYPE'] = [str('QSO')]*nqso
    # Return
    return meta
'''


def hdf5_adddata(hdf, sname, meta, debug=False, chk_meta_only=False,
                 mk_test_file=False):
    """ Append HST/FUSE data to the h5 file

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
    Rdicts = defs.get_res_dicts()
    # Add Survey
    print("Adding {:s} survey to DB".format(sname))
    hstc_grp = hdf.create_group(sname)
    # Checks
    if sname != 'UVpSM4':
        raise IOError("Not expecting this survey..")

    # Build spectra (and parse for meta)
    nspec = len(meta)
    max_npix = 40000  # Just needs to be large enough
    # Init
    data = init_data(max_npix, include_co=True)
    spec_set = hdf[sname].create_dataset('spec', data=data, chunks=True,
                                         maxshape=(None,), compression='gzip')
    spec_set.resize((nspec,))
    Rlist = []
    wvminlist = []
    wvmaxlist = []
    npixlist = []
    gratinglist = []
    datelist = []
    badf = []
    badstis = []
    badghrs = []
    # Loop
    path = os.getenv('RAW_IGMSPEC')+'/HST_Cooksey/'
    maxpix = 0
    for jj,row in enumerate(meta):
        # Generate full file
        full_file = path+'{:s}/{:s}/{:s}'.format(
                         row['QSO'],row['INSTR'],row['SPEC_FILE'])
        # Extract
        if row['INSTR'] == 'FUSE':
            hext = 1
        else:
            hext = 0
        print("HST_Cooksey: Reading {:s}".format(full_file))
        try:
            spec = lsio.readspec(full_file, head_exten=hext, masking='edges')
        except: # BAD HEADER
            hdu = fits.open(full_file)
            head1 = hdu[1].header
            hdu[1].verify('fix')
            tbl = Table(hdu[1].data)
            spec = lsio.readspec(tbl, masking='edges')
            spec.meta['headers'][spec.select] = head1
            # Continuum
            cfile = full_file.replace('.fits', '_c.fits')
            if os.path.isfile(cfile):
                # Watch that mask!
                gdp = ~spec.data['flux'][spec.select].mask
                spec.data['co'][spec.select][gdp] = (fits.open(cfile)[0].data)[gdp]
        # npix
        npix = spec.npix
        if npix > max_npix:
            raise ValueError("Not enough pixels in the data... ({:d})".format(npix))
        else:
            maxpix = max(npix,maxpix)
        # Some fiddling about
        for key in ['wave','flux','sig', 'co']:
            data[key] = 0.  # Important to init (for compression too)
        data['flux'][0][:npix] = spec.flux.value
        data['sig'][0][:npix] = spec.sig.value
        data['wave'][0][:npix] = spec.wavelength.value
        if spec.co_is_set:
            try:
                data['co'][0][:npix] = spec.co.value
            except ValueError:
                pdb.set_trace()
        # Meta
        datet = None
        if row['INSTR'] == 'FUSE':
            if 'HISTORY' in spec.header.keys():
                ncards = len(spec.header['HISTORY'])
                flg_H = True
            else:
                flg_H = False
                hdu = fits.open(full_file)
                head0 = hdu[0].header
                ncards = len(head0)
                # Is this a good one?
                if 'APER_ACT' in head0:
                    pass
                else:  # Need to fight harder for the header
                    # Look for untrim
                    untrim = full_file+'.untrim'
                    if not os.path.isfile(untrim):
                        pdb.set_trace()
                    # Read
                    hduu = fits.open(untrim)
                    if 'PKS2005' in untrim:  # One extra kludge..
                        head0 = hduu[1].header
                        flg_H = True
                        ncards = len(head0['HISTORY'])
                    else:
                        head0 = hduu[0].header
                        ncards = len(head0)
                spec.meta['headers'][spec.select] = head0
            # Read from history
            for ss in range(ncards):
                if flg_H:
                    try:
                        card = Header.fromstring(spec.header['HISTORY'][ss])
                    except:
                        pdb.set_trace()
                    try:
                        ckey = list(card.keys())[0]
                    except IndexError:
                        continue
                    else:
                        card0 = card[0]
                else:
                    ckey, card0 = list(spec.header.keys())[ss], spec.header[ss]
                # Parse
                if ckey == 'APERTURE':
                    aper = card0
                elif ckey == 'DETECTOR':
                    det = card0
                elif ckey == 'APER_ACT': # Extracted aperture
                    ext_ap = card0
                elif ckey == 'DATE': # Extracted aperture
                    datet = card0
            gratinglist.append(ext_ap+det)
        elif row['INSTR'] == 'STIS':
            try:
                datet = spec.header['DATE']
            except KeyError:  # handful of kludged coadds
                if 'HISTORY' not in spec.header.keys():
                    # Grab from the other extension, e.g. PKS0405
                    hdu = fits.open(full_file)
                    head1 = hdu[1].header
                    spec.meta['headers'][0] = head1
                for ihist in spec.header['HISTORY']:
                    if 'TDATEOBS' in ihist:
                        idash = ihist.find('-')
                        datet = ihist[idash-4:idash+6]
                # Grating from name
                i0 = full_file.rfind('_')
                i1 = full_file.rfind('.fits')
                gratinglist.append(full_file[i0+1:i1])
                if datet is None:
                    pdb.set_trace()
            else:
                gratinglist.append(spec.header['OPT_ELEM'])
        elif row['INSTR'] == 'GHRS':
            # Date
            try:
                tmp = spec.header['DATE-OBS']
            except KeyError:
                # Pull header from parallel file
                iM = full_file.find('M_1')
                if iM <= 0:
                    iM = full_file.find('L_1')
                ofile = full_file[:iM+1]+'_F.fits'
                if not os.path.isfile(ofile):
                    if 'NGC4151' in ofile:  # Kludge
                        ofile = ofile.replace('G160M', 'G160Mmd')
                    elif 'PKS2155-304_GHRS_G140L' in ofile:  # Kludge
                        ofile = ofile.replace('G140L', 'G140Llo')
                    elif 'PKS2155-304_GHRS_G160M' in ofile:  # Kludge
                        ofile = ofile.replace('G160M', 'G160Mmd')
                    else:
                        pdb.set_trace()
                hdu = fits.open(ofile)
                head0 = hdu[0].header
                spec.meta['headers'][spec.select] = head0
            # Reformat
            prs = tmp.split('/')
            if prs[2][0] == '9':
                yr = '19'+prs[2]
            else:
                yr = '20'+prs[2]
            datet = yr+'-'+prs[1]+'-{:02d}'.format(int(prs[0]))
            # Grating
            gratinglist.append(spec.header['GRATING'])
        else:
            pdb.set_trace()
        if datet is None:
            try:
                datet = spec.header['DATE-OBS']
            except KeyError:
                print("Missing Header for file: {:s}".format(full_file))
                badf.append(full_file)
                datet = '9999-9-9'
        t = Time(datet, format='isot', out_subfmt='date')  # Fixes to YYYY-MM-DD
        datelist.append(t.iso)
        try:
            Rlist.append(Rdicts[row['INSTR']][gratinglist[-1]])
        except KeyError:
            print(gratinglist[-1])
            pdb.set_trace()
        wvminlist.append(np.min(data['wave'][0][:npix]))
        wvmaxlist.append(np.max(data['wave'][0][:npix]))
        npixlist.append(npix)
        if chk_meta_only:
            continue
        # Only way to set the dataset correctly
        spec_set[jj] = data

    #
    if (len(badstis)) > 0:
        raise ValueError("Somehow have a bad STIS header..")
    if len(badf) > 0:
        print("We still have bad FUSE headers")
        pdb.set_trace()
    if len(badghrs) > 0:
        print("We still have bad GHRS headers")
        pdb.set_trace()
    print("Max pix = {:d}".format(maxpix))
    # Add columns
    meta.add_column(Column(npixlist, name='NPIX'))
    meta.add_column(Column(wvminlist, name='WV_MIN'))
    meta.add_column(Column(Rlist, name='R'))
    meta.add_column(Column(gratinglist, name='DISPERSER'))
    meta.add_column(Column(wvmaxlist, name='WV_MAX'))
    meta.add_column(Column(datelist, name='DATE-OBS'))
    meta.add_column(Column(np.arange(nspec,dtype=int),name='GROUP_ID'))

    # Add HDLLS meta to hdf5
    if chk_meta(meta):
        if chk_meta_only:
            pdb.set_trace()
        hdf[sname]['meta'] = meta
    else:
        pdb.set_trace()
        raise ValueError("meta file failed")
    # References
    refs = [dict(url='http://adsabs.harvard.edu/abs/2010ApJ...708..868C',
                 bib='cooksey10')
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
    Title = '{:s}: HST and FUSE spectra of AGN and Quasars by Cooksey et al. (2010)'.format(dset)
    ssa_dict = default_fields(Title, flux='flambda')
    hdf[dset]['meta'].attrs['SSA'] = json.dumps(ltu.jsonify(ssa_dict))
