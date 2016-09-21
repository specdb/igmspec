.. highlight:: rest

***************
ESI_DLA Dataset
***************

This document describes the ESI_DLA dataset.

Sources
=======

The High z ESI DLA sample is a set of Keck/ESI
echellette spectra acquired and published by
`Rafelski et al. (2012) <http://adsabs.harvard.edu/abs/2012ApJ...755...89R>`_
and
`Rafelski et al. (2014) <http://adsabs.harvard.edu/abs/2014ApJ...782L..29R>`_
for the analysis of z>4 DLAs.
The quasars were all drawn from the SDSS.

The data are also available online at
`here <http://www.rafelski.com/data/DLA/hizesi/>`_ along with some basic information.

Meta Data
=========

The meta data provided with ESI_DLA includes some basic information
on the quasars (r,i,z photometry) and a S/N estimate.
The redshifts were taken from SDSS.

============  ======== =========================================
Key           Type     Description
============  ======== =========================================
r_mag         float    Quasar r_mag from SDSS
SN            float    Estimate of S/N per pixel

Name           str        Name
RA                float    Right Ascension (degrees)
DEC              float    Declination (degrees)
Plate             int       SDSS Plate
MJD              int       SDSS Modified Julian Date
FiberID         int       SDSS Fiber ID
zem              float    Quasar emission redshift
r_mag           float    Quasar r_mag from SDSS
i_mag           float    Quasar i_mag from SDSS
z_mag          float    Quasar z_mag from SDSS
ObsDate       str       Date observed with Keck
Exptime        int       Total exposure time
SN                int        Estimate of S/N per pixel
Slit               float     Slit used (0.75" or 0.6")
Reference     str        Data reference paper
DATE-OBS    str        Observed date in 'DATE-OBS' format, YYYY-MM-DD
sig_zem       float     Uncertainty in Quasar emission redshift
flag_zem      str        Where the zem comes from
EPOCH         float     Epoch
IGM_ID        int         Individual IGM ID
SPEC_FILE    str        Filename of spectrum
NPIX            int        Number of pixels
WV_MIN      float     Minimum wavelength
WV_MAX     float     Maximum wavelength
R                 float    Resolution 
SURVEY_ID   int      Survey ID number
TELESCOPE  str       Telescope of observations
INSTR          str       Instrument of observations
GRATING     str       Grating type used
============  ======== =========================================


Spectra
=======

All of the data are fluxed ESI spectra.
