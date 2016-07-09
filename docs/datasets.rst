.. highlight:: rest

*******************
Datasets of igmspec
*******************

This document describes the datasets of igmspec.

List of Surveys
===============

Below is the list of surveys included in igmspec v1.0

* :doc:`boss`
* :doc:`hdlls`
* :doc:`kodiaq`

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
DATE-OBS    str      Date observed (see individual survey notes)
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

The instruments used in igmspec are provided in igmspec.defs.instruments.
The following Table summarizes and defines the instruments
used in igmspec v1.0:

==========  ======== ============================================
Instrument  Gratings Description
==========  ======== ============================================
BOSS        BLUE     Blue channel spectrograph
 ..         RED      Red channel spectrograph
 ..         BOTH     Spectrum includes data from both spectrographs
ESI         ECH      Echelette mode on Keck/ESI instrument
GMOS-N      R400     Gemini North GMOS spectrometer
 ..         B600     ..
GMOS-S      R400     Gemini South GMOS spectrometer
 ..         B600     ..
HIRES       BLUE     Blue cross-disperser on HIRES (aka HIRESb)
 ..         RED      Red cross-dispereser on HIRES (aka HIRESr)
 ..         BOTH     Spectrum includes data from both cross-dispersers
MagE        N/A      MagE spectrometer
MIKEb       BLUE     Blue camera of MIKE spectrometer
MIKEr       RED      Red camera of MIKE spectrometer
MIKE        BOTH     Spectrum is a splice of MIKEb and MIKEr data
SDSS        BLUE     Blue channel spectrograph
 ..         RED      Red channel spectrograph
 ..         BOTH     Spectrum includes data from both spectrographs
==========  ======== ============================================

Telescopes
----------

Here are the telescopes currently incorporated in igmspec v1.0:

==============  ====================================================
Telescope       Website
==============  ====================================================
Gemini-N        http://www.gemini.edu
Gemini-S        http://www.gemini.edu
Keck I          http://www.keckobservatory.org/
Keck II         http://www.keckobservatory.org/
Magellan/Clay   http://obs.carnegiescience.edu/Magellan
Magellan/Baade  http://obs.carnegiescience.edu/Magellan
SDSS 2.5-M      https://www.sdss3.org/instruments/telescope.php
==============  ====================================================


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

The software included with igmspec read these data into
a XSpectrum1D object from
`linetools <http://linetools.readthedocs.io/en/latest/>`_.
