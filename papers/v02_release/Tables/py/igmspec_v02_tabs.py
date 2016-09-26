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
    surveys = igmsp.idb.hdf.keys()
    surveys.pop(surveys.index('catalog'))
    surveys.sort()

    # Open
    tbfil = open(outfil, 'w')

    # Header
    tbfil.write('\\clearpage\n')
    tbfil.write('\\begin{deluxetable}{lccccccc}\n')
    tbfil.write('\\rotate\n')
    tbfil.write('\\tablewidth{0pc}\n')
    tbfil.write('\\tablecaption{{\\it igmspec} DATASETS \\label{tab:datasets}}\n')
    #tbfil.write('\\tabletypesize{\\tiny}\n')
    tbfil.write('\\tablehead{\\colhead{Survey} & \\colhead{$N_{\\rm source}^a$} \n')
    tbfil.write('& \\colhead{$N_{\\rm spec}^b$} & \\colhead{$\\lambda_{\\rm min}$}\n')
    tbfil.write('& \\colhead{$\\lambda_{\\rm max}$} & \\colhead{$R^c$} \n')
    tbfil.write('& \\colhead{References} & \\colhead{Website} \n')
    tbfil.write('} \n')

    tbfil.write('\\startdata \n')

    # Looping on systems
    restrict = False
    for survey in surveys:
        if survey == 'quasars':
            continue
        # Restrict
        if restrict:
            if survey == 'BOSS_DR12':
                pdb.set_trace()
            elif survey == 'SDSS_DR7':
                print("SKIPPING SDSS FOR NOW")
                continue
        print("Working on survey={:s}".format(survey))
        # Setup
        meta = Table(igmsp.idb.hdf[survey]['meta'].value)

        # Survey
        survey.replace('_','\\_')
        tbfil.write(survey)

        # N sources
        uniq = np.unique(meta['IGM_ID'])
        tbfil.write('& {:d}'.format(len(uniq)))

        # N spectra
        tbfil.write('& {:d}'.format(len(meta)))

        # Wave min
        sig = igmsp.idb.hdf[survey]['spec']['sig']
        gds = sig > 0.
        gdwv = igmsp.idb.hdf[survey]['spec']['wave'][gds]
        tbfil.write('& {:0.1f}'.format(np.min(gdwv)))
        # Wave max
        tbfil.write('& {:0.1f}'.format(np.max(gdwv)))

        # R
        tbfil.write('& {:d}'.format(int(np.round(np.median(meta['R'])))))

        # Write
        tbfil.write('\\\\ \n')


    # End
    tbfil.write('\\hline \n')


    tbfil.write('\\enddata \n')
    #tbfil.write('\\tablecomments{This table is available as a YAML file at ')
    #tbfil.write('http://blah')
    #tbfil.write('} \n')
    tbfil.write('\\tablenotetext{a}{Number of positive detections constraining the model.}')
    # End
    tbfil.write('\\end{deluxetable} \n')

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

