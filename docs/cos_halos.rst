.. highlight:: rest

*****************
COS-Halos Dataset
*****************

This document describes the COS-Halos dataset.

Sources
=======

The COS-Halos Survey
is a set of HST/COS
`Tumlinson et al. (2013) <http://adsabs.harvard.edu/abs/2013ApJ...777...59T>`_
and Keck/HIRES spectra
`Werk et al. (2013) <http://adsabs.harvard.edu/abs/2013ApJS..204...17W>`_
for the analysis of the circumgalactic medium of z~0.2 galaxies.
The quasars are a set of sources from the SDSS with
zem ~ 0.3 to 1. and bright in the FUV.


Meta Data
=========

The additional meta data are from Table 1 and include
the following for the HST/COS rows.

============  ======== =========================================
Key           Type     Description
============  ======== =========================================
m_FUV         float    FUV magnitude from GALEX
t_G130M       float    Exposure time (sec) for the G130M grating
t_G160M       float    Exposure time (sec) for the G160M grating
Visit         str      Visit ID(s)
============  ======== =========================================

Use `pyigm <http://https://github.com/pyigm/pyigm>`_
to access the measurements from the COS-Halos survey.

Spectra
=======

The HST/COS spectra are wavelength and flux calibrated spectra using
CALCOS and then coadded with proprietary code
(see Tumlinson et al. 2013).

The Keck/HIRES spectra are standard HIRedux reductions
(see Werk et al. 2013).
