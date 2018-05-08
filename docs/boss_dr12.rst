.. highlight:: rest

*****************
BOSS DR12 Dataset
*****************

This document describes the BOSS DR12 dataset.

Notebooks
=========

.. toctree::
    :maxdepth: 1

            BOSS_DR12 <BOSS_DR12>

Sources
=======

The igmspec database is intended to include all of the
quasar spectra in the BOSS DR12 data release.  The current
version includes the following catalogs:

============  =============================================== ===========
Catalog       Description                                     Link
============  =============================================== ===========
DR12Q         Official DR12 quasar release of the             `DR12Q <http://data.sdss3.org/datamodel/files/BOSS_QSO/DR12Q/DR12Q.html>`_
              BOSS survey.  This is 297,301 quasars
DR12Q_sup     Supplemental quasars                            `DR12Q_sup <https://data.sdss.org/datamodel/files/BOSS_QSO/DR12Q/DR12Q_sup.html>`_
DR12Q_supbad  Bad spectra identified as quasars               No unique link
============  =============================================== ===========

For DR12Q, I have ignored the ~50 sources with Z_VI>0 but
Z_PCA and Z_PIPE junk.  I noted that at least one of these
doesn't have a proper spectrum.  

Meta Data
=========

See the links provided above to see the full set of meta data
provided with BOSS catalogs.  It is quite extensive.

The *flag_co* column indicates the source of the continuum.  It is a
bitwise flag with 1=Zhu et al.; 2=Lee et al.

Spectra
=======

v1.0 now includes all of the DR12 spectra.  Note that some of these were
downloaded by JXP (those at z>4.8 and those past v_5_7_0), while
the rest were kindly provided by G. Zhu.

Continua
========

There are two sets of continua included in igmspec.  As a default,
we use the continua kindly provided by G. Zhu using the methodology
described in this paper:
`Zhu et al. 2014 <http://adsabs.harvard.edu/abs/2014MNRAS.439.3139Z>`_.
These are primarily useful for analysis *outside* of the Lya forest.
Note that these were only generated for sources with z<4.8 and
for data processed in v_5_7_0.

When available, we have used mean flux regulated continua kindly
provided by K.G. Lee, as described in these papers:
`Lee et al. 2012 <http://adsabs.harvard.edu/abs/2012AJ....143...51L>`_,
`Lee et al. 2013 <http://adsabs.harvard.edu/abs/2013AJ....145...69L>`_.
This continuum is restricted to rest-frame wavelengths of less than 1200A.
