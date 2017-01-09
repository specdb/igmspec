.. highlight:: rest

**************
KODIAQ Dataset
**************

This document describes the KODIAQ dataset.

Sources
=======

The Keck Observatory Database of Ionized Absorption toward Quasars (KODIAQ)
survey is a public spectra data release of quasars observed with
the HIRES spectrometer on the Keck I telescope.  The first Data Release
(DR1) became available in 2015, as described in
`O'Meara et al. (2015) <http://adsabs.harvard.edu/abs/2015AJ....150..111O>`_.


Meta Data
=========

The additional meta data provided with KODIAQ follows from Table 1 of
O'Meara et al. (2015).  These are:

============  ====== =========================================
Key           Type   Description
============  ====== =========================================
pi_date       str    PI and Observation Date
spec_prefix   str    Spectrum prefix
redux_setup   str    Reduction setup
targname      str    Target name set by observer
deckname      str    Name of the HIRES decker used
elaptime      str    Description of the total exposure time (seconds)
qaflag        str    Description of the data reduction process
kodetime      int    Total exposure time of the co-added spectrum (seconds)
kodwblue      int    Approximate starting wavelength (Ang)
kodwred       int    Approximate ending wavelength (Ang)
kodrelease    int    KODIAQ data release (1=DR1)
============  ====== =========================================


Spectra
=======

In igmspec v02 are all of the spectra released in DR1 of
KODIAQ.  These data are continuum normalized.
