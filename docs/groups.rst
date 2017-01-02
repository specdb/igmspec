.. highlight:: rest

*****************
Groups of igmspec
*****************

This document describes the data groups of igmspec.

List of Surveys
===============

Below is the list of data groups included in igmspec v02:

* :doc:`boss`  [BOSS_DR12]
* :doc:`cos_dwarfs`   [COS-Dwarfs]
* :doc:`cos_halos`   [COS-Halos]
* :doc:`esidla`   [ESI_DLA]
* :doc:`ggg`   [GGG]
* :doc:`hdla100`   [HDLA100]
* :doc:`hdlls` [HD-LLS_DR1]
* :doc:`hst_z2`   [HST_z2]
* :doc:`hst_qso`   [HSTQSO]
* :doc:`kodiaq` [KODIAQ_DR1]
* :doc:`musodla`   [MUSoDLA]
* :doc:`xq100`   [XQ-100]
* :doc:`sdss`   [SDSS_DR7]
* :doc:`uves_dall`   [UVES_Dall]
* :doc:`uvpsm4`   [UVpSM4]
* :doc:`twoqz`   [2QZ]

Each document provides the survey reference
and additional details on the spectra and
associated meta data.

Overview
========
Each group included in igmspec is composed of two
components:

1. A Table of meta data
2. A numpy data array containing the spectra

Meta Data
=========

Each group has its own unique set of meta data describing
the data products.  The following keys are required for
inclusion in igmspec:

==========  ======== ============================================
Key         Type     Description
==========  ======== ============================================
IGM_ID      int      Unique igmspec identifier
zem_GROUP   float    Emission redshift of background source given by the survey
RA_GROUP    float    Right Ascension (deg) given by the survey
DEC_GROUP   float    Declination (deg) given by the survey
EPOCH       float    Coordinate epoch (only 2000 in igmspec v1.0)
DATE-OBS    str      Date observed (YYYY-MM-DD)
R           float    Instrument resolution, :math:`\lambda/\Delta\lambda` (FWHM)
WV_MIN      float    Minimum wavelength of the spectrum
WV_MAX      float    Maximum wavelength of the spectrum
NPIX        int      Number of pixels in the spectrum; may include null values
GROUP_ID    int      Unique identifier for the group [not well implemented yet]
SPEC_FILE   str      Spectrum file name
INSTR       str      Instrument file name (see `Instruments and Dispersers`_ for definitions)
DISPERSER   str      Disperser name (see `Instruments and Dispersers`_ for definitions)
TELESCOPE   str      Telescope name (see `Telescopes`_ below for definitions)
==========  ======== ============================================

Additional meta data may be provided for
individual surveys.

.. _Instruments and Dispersers:

Instruments and Dispersers
--------------------------

The complete list of instruments and associated
dispersers that may be used in igmspec
are provided in the
`specdb <http://specdb.readthedocs.io/en/latest/>`_
documentation.

.. _Telescopes:

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

The software included with
`specdb <http://specdb.readthedocs.io/en/latest/>`_
read these data into a XSpectrum1D object in the
`linetools <http://linetools.readthedocs.io/en/latest/>`_
software repository.
