.. highlight:: rest

*******************
Datasets of igmspec
*******************

This document describes the datasets of igmspec.

List of Surveys
===============

Below is the list of surveys included in igmspec v02

* :doc:`boss`  [BOSS_DR12]
* :doc:`hdlls` [HD-LLS_DR1]
* :doc:`kodiaq` [KODIAQ_DR1]
* :doc:`sdss`   [SDSS_DR7]
* :doc:`ggg`   [GGG]
* :doc:`hst_z2`   [HST_z2]
* :doc:`twoqz`   [2QZ]
* :doc:`esidla`   [ESI_DLA]
* :doc:`xq100`   [XQ-100]
* :doc:`cos_halos`   [COS-Halos]
* :doc:`cos_dwarfs`   [COS-Dwarfs]
* :doc:`musodla`   [MUSoDLA]
* :doc:`hst_qso`   [HSTQSO]
* :doc:`hdla100`   [HDLA100]

Overview
========
Each dataset included in igmspec is composed of two
components:

1. A Table of meta data
2. A numpy data array containing the spectra

Meta Data
=========

Each survey has its own unique set of meta data describing
the data products.  The following keys are required for
inclusion in igmspec:

==========  ======== ============================================
Key         Type     Description
==========  ======== ============================================
IGM_ID      int      Unique igmspec identifier
zem         float    Emission redshift of background source
RA          float    Right Ascension (deg)
DEC         float    Declination (deg)
EPOCH       float    Coordinate epoch (only 2000 in igmspec v1.0)
DATE-OBS    str      Date observed (YYYY-MM-DD)
R           float    Instrument resolution, :math:`\lambda/\Delta\lambda` (FWHM)
WV_MIN      float    Minimum wavelength of the spectrum
WV_MAX      float    Maximum wavelength of the spectrum
NPIX        int      Number of pixels in the spectrum; may include null values
SURVEY_ID   int      Unique identifier for the survey [not well implemented yet]
SPEC_FILE   str      Spectrum file name
INSTR       str      Instrument file name (see `Instruments and Gratings`_ below for definitions)
GRATING     str      Grating name (see `Instruments and Gratings`_ below for definitions)
TELESCOPE   str      Telescope name (see `Telescopes`_ below for definitions)
==========  ======== ============================================


Instruments and Gratings
------------------------

The complete list of instruments that may be
used in igmspec are provided in the
`specdb <http://specdb.readthedocs.io/en/latest/>`_
documentation.


Telescopes
----------

Similarly, the list of telescopes that may be used
in igmspec are provided in the
`specdb <http://specdb.readthedocs.io/en/latest/>`_
documentation.


Spectral Data
=============

The spectra in igmspec are written as a numpy masked array with
three required columns and one optional:

=============  ======= =============================================
Key            Type    Description
=============  ======= =============================================
wave           float64 Wavelength array; default is Angstroms
flux           float32 Flux array; default is unitless
sig            float32 Error array; same units as flux
co (optional)  float32 Continuum array; same units as flux
=============  ======= =============================================

The software included with specdb read these data into
a XSpectrum1D object from
`linetools <http://linetools.readthedocs.io/en/latest/>`_.
