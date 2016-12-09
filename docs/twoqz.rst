.. highlight:: rest

***********
2QZ Dataset
***********

This document describes the 2QZ dataset as
published by
`Croom et al. (2004) <http://adsabs.harvard.edu/abs/2004MNRAS.349.1397C>`_.

Notebooks
=========

.. toctree::
    :maxdepth: 1

        2QZ <2QZ>

Sources
=======

The igmspec database is intended to include all of the
quasar spectra in the 2QZ data release.  I have acheived
this by querying the main
`web server <http://www.2dfquasar.org/Spec_Cat/2qzsearch2.html>`_,
restricting to sources with zmin>0.05.

Within the code, I then take any source with "QSO" in
the ID1 or ID2 columns of the meta file.


Meta Data
=========

See `catalog description <http://www.2dfquasar.org/Spec_Cat/catalogue.html>`_
for a full description of the meta data.

Spectra
=======

Spectra are straight from the long format of the
server.  Bad pixels have sig=0.  Note that we have
not included spectra with error arrays having all zeros
or those without a spectrum available via the website.
These do not occur within the source catalog.
