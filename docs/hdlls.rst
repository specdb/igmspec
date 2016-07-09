.. highlight:: rest

***************
HD-LLS Dataaset
***************

This document describes the HD-LLS dataset.

Sources
=======

The high dispersion Lyman Limit System (HD-LLS) sample is a set of
echelle and echellette spectra acquired by
`Prochaska et al. (2015) <http://adsabs.harvard.edu/abs/2015ApJS..221....2P>`_
for the analysis of z~3 LLS.
The quasars are a heterogenous set of sources useful
for such analysis (ie. bright).


Meta Data
=========

The meta data provided with KODIAQ follows from Table 1 of
O'Meara et al. (2015).  These are:

============  ====== =========================================
Key           Type   Description
============  ====== =========================================
pi_date       str    PI and Observation Date
spec_prefix   str    Spectrum prefix
redux_setup   str    Reduction setup
targname      str    Target name set by observer
deckname      str    Name of the HIRES decker used
elaptime      str    Description of the total exposure time
qaflag        str    Description of the data reduction process
kodetime      int    Total exposure time of the co-added spectrum (seconds)
kodwblue      int    Starting wavelength (Ang)
kodwred       int    Ending wavelength (Ang)
kodrelease    int    KODIAQ data release (1=DR1)
============  ====== =========================================


Spectra
=======

In igmspec v1.0 are all of the spectra released in DR1 of
KODIAQ.  These were modified, however, to no longer include
the BZERO header card.
