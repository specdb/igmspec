#Module for Tables for the igmspec v02 paper
# Imports
from __future__ import print_function, absolute_import, division, unicode_literals


import numpy as np
import os, sys
import json, yaml
import pdb

from astropy.table import Table
#from astropy import units as u


from specdb.specdb import IgmSpec


# Local
#sys.path.append(os.path.abspath("../Analysis/py"))
#import lls_sample as lls_s

#def mktab_all_lines -- LateX table listing all measured transitions 


# Summary table of the NHI models
def mktab_datasets(outfil='tab_datasets.tex'):
    """ Generate a Table summarizing the datasets

    Parameters
    ----------
    outfil

    Returns
    -------

    """
    # Setup [confirm v02 eventually]
    igmsp = IgmSpec()
    groups = igmsp.groups
    groups.sort()

    # Open
    tbfil = open(outfil, 'w')

    # Header
    tbfil.write('\\clearpage\n')
    tbfil.write('\\begin{table}[ht]\n')
    tbfil.write('\\caption{{\\it igmspec} DATA GROUPS \\label{tab:datasets}}\n')
    #tbfil.write('\\tabletypesize{\\tiny}\n')
    tbfil.write('\\begin{tabular}{lcccccc}\n')
    tbfil.write('Group & $N_{\\rm source}^a$ \n')
    tbfil.write('& $N_{\\rm spec}^b$ & $\\lambda_{\\rm min}$ (\\AA) \n')
    tbfil.write('& $\\lambda_{\\rm max}$ (\\AA) & $R^c$  & Flux$^d$\\\\ \n')
    #tbfil.write('& References & Website \n')
    #tbfil.write('} \n')

    tbfil.write('\\hline \n')

    # Looping on systems
    restrict = True
    for survey in groups:
        if survey == 'quasars':
            continue
        # Restrict
        #if survey != 'HD-LLS_DR1':
        #    continue
        if restrict:
            if survey in ['BOSS_DR12', 'SDSS_DR7']:
                continue
        print("Working on survey={:s}".format(survey))
        # Setup
        meta = Table(igmsp.hdf[survey]['meta'].value)
        try:
            ssa = json.loads(igmsp.hdf[survey+'/meta'].attrs['SSA'])
        except KeyError:
            if survey not in ['COS-Halos']:
                pdb.set_trace()
            fluxc = 'MIXED'
        else:
            fluxc = ssa['FluxCalib']

        # Survey
        msurvey = survey.replace('_','\\_')
        tbfil.write(msurvey)

        # N sources
        uniq = np.unique(meta['IGM_ID'])
        tbfil.write('& {:d}'.format(len(uniq)))

        # N spectra
        tbfil.write('& {:d}'.format(len(meta)))

        # Wave min
        sig = igmsp.hdf[survey]['spec']['sig']
        gds = sig > 0.
        gdwv = igmsp.hdf[survey]['spec']['wave'][gds]
        tbfil.write('& {:d}'.format(int(np.round(np.min(gdwv)))))
        # Wave max
        tbfil.write('& {:d}'.format(int(np.round(np.max(gdwv)))))

        # R
        tbfil.write('& {:d}'.format(int(np.round(np.median(meta['R'])))))

        # Flux
        tbfil.write('& {:s}'.format(fluxc))

        # Write
        tbfil.write('\\\\ \n')


    # End
    tbfil.write('\\hline \n')


    #tbfil.write('\\enddata \n')
    #tbfil.write('\\tablecomments{This table is available as a YAML file at ')
    #tbfil.write('http://blah')
    #tbfil.write('} \n')
    tbfil.write('\\multicolumn{6}{l}{{$^a$}{Number of unique sources in the dataset. }} \\\\ \n')
    tbfil.write('\\multicolumn{6}{l}{{$^b$}{Number of unique spectra in the dataset. }} \\\\ \n')
    tbfil.write('\\multicolumn{6}{l}{{$^c$}{Characteristic FWHM resolution of the spectra. }} \\\\ \n')
    tbfil.write('\\multicolumn{6}{l}{{$^d$}{Indicates whether the data are fluxed (absolute or relative) or normalized. The COS-Halos spectra include both fluxed (COS) and normalized (HIRES) spectra.}} \\\\ \n')
    #tbfil.write('\\tablenotetext{a}{Number of positive detections constraining the model.}')
    # End
    tbfil.write('\\end{tabular} \n')
    tbfil.write('\\end{table} \n')

    tbfil.close()



#### ########################## #########################
#### ########################## #########################
#### ########################## #########################

# Command line execution
if __name__ == '__main__':

    flg_tab = 0 
    flg_tab += 2**0  # Datasets
    #flg_tab += 2**1  # Ionization models
    #flg_tab += 2**2  # Edits to COS-Halos

    # NHI fits
    if (flg_tab % 2**1) >= 2**0:
        mktab_datasets()

