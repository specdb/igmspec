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

The data are also available
`online <http://www.rafelski.com/data/DLA/hizesi/>`_
along with some additional basic information
about the survey.

Meta Data
=========

The meta data provided with ESI_DLA includes some basic information
on the quasars (r,i,z photometry) and a S/N estimate.
The redshifts were taken from SDSS.

============ ======== =========================================
Key          Type     Description
============ ======== =========================================
SN            float   Estimate of S/N per 11 km/s pixel
Name          str     Name
Plate         int     SDSS Plate
MJD           int     SDSS Modified Julian Date
FiberID       int     SDSS Fiber ID
zem           float   Quasar emission redshift
r_mag         float   Quasar r_mag from SDSS
i_mag         float   Quasar i_mag from SDSS
z_mag         float   Quasar z_mag from SDSS
ObsDate       str     Date observed with Keck
Exptime       int     Total exposure time (seconds)
Slit          float   Slit used in arcseconds (0.75" or 0.5")
Reference     str     Data reference paper
============ ======== =========================================


Spectra
=======

All of the data are fluxed ESI spectra.
