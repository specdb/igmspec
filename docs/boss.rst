.. highlight:: rest

************
BOSS Dataset
************

This document describes the BOSS dataset.

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

==========  =============================================== ===========
Catalog     Description                                     Link
==========  =============================================== ===========
DR12Q       Official DR12 quasar release of the             `DR12Q <http://data.sdss3.org/datamodel/files/BOSS_QSO/DR12Q/DR12Q.html>`_
            BOSS survey.  This is 297,301 quasars
==========  =============================================== ===========

Meta Data
=========

See the links provided above to see the full set of meta data
provided with BOSS catalogs.  It is quite extensive.


Spectra
=======

v1.0 now includes all of the DR12 spectra.

Continua
========

There are two sets of continua included in igmspec.  As a default,
we use the continua kindly provided by G. Zhu using the methodology
described in this paper:
`Zhu et al. 2014 <http://adsabs.harvard.edu/abs/2014MNRAS.439.3139Z>`.
These are primarily useful for analysis *outside* of the Lya forest.

When available, we have used mean flux regulated continua kindly
provided by K.G. Lee, as described in these papers:
`Lee et al. 2012 <http://adsabs.harvard.edu/abs/2012AJ....143...51L>`,
`Lee et al. 2013 <http://adsabs.harvard.edu/abs/2013AJ....145...69L>`.
