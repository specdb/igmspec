
Simple Examples for Using IGMspec (v1.1.1)
==========================================

:download:`Download <examples/Simple_Usage.ipynb>` this notebook.

.. code:: python

    %matplotlib inline

.. code:: python

    # import
    from astropy import units as u
    
    from igmspec import query_catalog as iqcat
    from igmspec import interface_db as igidb
    from igmspec import db_utils as idbu
    from igmspec.igmspec import IgmSpec
    
    from pyigm.surveys.llssurvey import LLSSurvey

Setup Class
-----------

.. code:: python

    igmsp = IgmSpec()


.. parsed-literal::

    Using /raid/IGMSPEC_DB/IGMspec_DB_v01.hdf5 for the catalog file
    Using /raid/IGMSPEC_DB/IGMspec_DB_v01.hdf5 for the DB file
    Available surveys: [u'GGG', u'HD-LLS_DR1', u'KODIAQ_DR1', u'SDSS_DR7']


.. code:: python

    igmsp.qcat




.. parsed-literal::

    <QueryCatalog:  DB_file=/raid/IGMSPEC_DB/IGMspec_DB_v01.hdf5 with 377018 sources
       Loaded surveys are [u'BOSS_DR12', u'GGG', u'HD-LLS_DR1', u'KODIAQ_DR1', u'SDSS_DR7'] 
    >



.. code:: python

    igmsp.idb




.. parsed-literal::

    <InterfaceDB:  DB_file=/raid/IGMSPEC_DB/IGMspec_DB_v01.hdf5 
       Loaded surveys are [u'GGG', u'HD-LLS_DR1', u'KODIAQ_DR1', u'SDSS_DR7'] 
    >



Radial search
-------------

Search around FJ0812+32
~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    ids0812 = igmsp.radial_search(('08:12:40.68','+32:08:09'), 1.*u.arcsec)
    igmsp.show_cat(ids0812)


.. parsed-literal::

    Your search yielded 1 match[es]
    IGM_ID    RA      DEC     zem   flag_survey sig_zem  flag_zem
    ------ -------- -------- ------ ----------- ------- ---------
     50734 123.1695  32.1357  2.698           7   0.000 BOSS_PCA 
    ----------
    Survey key:
        BOSS_DR12: 1
        GGG: 16
        HD-LLS_DR1: 8
        KODIAQ_DR1: 4
        SDSS_DR7: 2


.. code:: python

    # Grab meta
    meta0812 = igmsp.grab_meta('SDSS_DR7', ids0812, show=True)

.. code:: python

    meta0812




.. raw:: html

    &lt;Table length=1&gt;
    <table id="table4581353360">
    <thead><tr><th>zem</th><th>sig_zem</th><th>Z_CONF</th><th>Z_WARN</th><th>PLATE</th><th>MJD</th><th>FIBERID</th><th>FLG_TARG</th><th>RA</th><th>DEC</th><th>PSF_U</th><th>PSF_G</th><th>PSF_R</th><th>PSF_I</th><th>PSF_Z</th><th>PSF_SU</th><th>PSF_SG</th><th>PSF_SR</th><th>PSF_SI</th><th>PSF_SZ</th><th>DATE-OBS</th><th>EPOCH</th><th>R</th><th>IGM_ID</th><th>SPEC_FILE</th><th>NPIX</th><th>WV_MIN</th><th>WV_MAX</th><th>SURVEY_ID</th><th>INSTR</th><th>GRATING</th><th>TELESCOPE</th></tr></thead>
    <thead><tr><th>float32</th><th>float32</th><th>float32</th><th>int16</th><th>int32</th><th>int32</th><th>int32</th><th>int16</th><th>float64</th><th>float64</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>str10</th><th>float64</th><th>float64</th><th>int64</th><th>str28</th><th>int64</th><th>float64</th><th>float64</th><th>int64</th><th>str4</th><th>str4</th><th>str10</th></tr></thead>
    <tr><td>2.704</td><td>0.0015619</td><td>0.0</td><td>0</td><td>861</td><td>52318</td><td>333</td><td>0</td><td>123.170</td><td>32.136</td><td>19.6334</td><td>17.8486</td><td>17.4545</td><td>17.2825</td><td>17.139</td><td>0.02811</td><td>0.017889</td><td>0.021285</td><td>0.018897</td><td>0.028298</td><td>2002-02-13</td><td>2000.0</td><td>2000.0</td><td>50734</td><td>spSpec-52318-0861-333.fit.gz</td><td>3821</td><td>3823.8</td><td>9215.1</td><td>9317</td><td>SDSS</td><td>BOTH</td><td>SDSS 2.5-M</td></tr>
    </table>



.. code:: python

    # Grab spectra
    J0812spec, meta = igmsp.grab_spec('KODIAQ_DR1', ids0812)


.. parsed-literal::

    Staged 1 spectra totalling 0.0032 Gb
    Loaded spectra


.. parsed-literal::

    /Users/xavier/local/Python/linetools/linetools/spectra/xspectrum1d.py:291: UserWarning: Assuming wavelength unit is Angstroms
      warnings.warn("Assuming wavelength unit is Angstroms")


.. code:: python

    J0812spec.plot()



.. image:: Simple_Usage_files/Simple_Usage_13_0.png


Search around J233446.40-090812.3
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    ids2334 = igmsp.radial_search('J233446.40-090812.3', 1.*u.arcsec)
    igmsp.show_cat(ids2334)


.. parsed-literal::

    Your search yielded 1 match[es]
    IGM_ID    RA      DEC     zem   flag_survey sig_zem flag_zem
    ------ -------- -------- ------ ----------- ------- --------
    376530 353.6934  -9.1368  3.317          14   0.001     SDSS
    ----------
    Survey key:
        BOSS_DR12: 1
        GGG: 16
        HD-LLS_DR1: 8
        KODIAQ_DR1: 4
        SDSS_DR7: 2


.. code:: python

    # Grab meta
    meta2334 = igmsp.grab_meta('HD-LLS_DR1', ids2334, show=True)

.. code:: python

    meta2334




.. raw:: html

    &lt;Table length=3&gt;
    <table id="table4583427920">
    <thead><tr><th>Name</th><th>QSO</th><th>RA</th><th>DEC</th><th>zem</th><th>Z_LLS</th><th>logNHI</th><th>sig(logNHI) [2]</th><th>SPEC_FILE</th><th>IGM_ID</th><th>EPOCH</th><th>NPIX</th><th>DATE-OBS</th><th>WV_MIN</th><th>WV_MAX</th><th>R</th><th>SURVEY_ID</th><th>GRATING</th><th>INSTR</th><th>TELESCOPE</th></tr></thead>
    <thead><tr><th>str33</th><th>str19</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>str37</th><th>int64</th><th>float64</th><th>int64</th><th>str10</th><th>float64</th><th>float64</th><th>float64</th><th>int64</th><th>str4</th><th>str5</th><th>str13</th></tr></thead>
    <tr><td>HD-LLS_J233446.40-090812.3_z3.226</td><td>SDSSJ2334-0908</td><td>353.693</td><td>-9.137</td><td>3.317</td><td>3.226</td><td>17.7</td><td>0.1 .. 0.3</td><td>HD-LLS_J233446.40-090812.3_ESI.fits</td><td>376530</td><td>2000.0</td><td>33000</td><td>2002-12-02</td><td>3899.8</td><td>11714.9</td><td>6000.0</td><td>139</td><td>ECH</td><td>ESI</td><td>Keck-II</td></tr>
    <tr><td>HD-LLS_J233446.40-090812.3_z3.226</td><td>SDSSJ2334-0908</td><td>353.693</td><td>-9.137</td><td>3.317</td><td>3.226</td><td>17.7</td><td>0.1 .. 0.3</td><td>HD-LLS_J233446.40-090812.3_HIRES.fits</td><td>376530</td><td>2000.0</td><td>129277</td><td>2007-09-18</td><td>4064.6</td><td>7120.0</td><td>48000.0</td><td>140</td><td>BOTH</td><td>HIRES</td><td>Keck-I</td></tr>
    <tr><td>HD-LLS_J233446.40-090812.3_z3.226</td><td>SDSSJ2334-0908</td><td>353.693</td><td>-9.137</td><td>3.317</td><td>3.226</td><td>17.7</td><td>0.1 .. 0.3</td><td>HD-LLS_J233446.40-090812.3_MAGE.fits</td><td>376530</td><td>2000.0</td><td>16580</td><td>2010-08-13</td><td>3042.1</td><td>10269.6</td><td>5857.14285714</td><td>141</td><td>N/A</td><td>MagE</td><td>Magellan/Clay</td></tr>
    </table>



.. code:: python

    # Grab spectra
    J2334spec, meta_2334 = igmsp.grab_spec('HD-LLS_DR1', ids2334)


.. parsed-literal::

    Staged 3 spectra totalling 0.01008 Gb
    Loaded spectra


.. code:: python

    # Plot the first one (ESI)
    J2334spec.plot(inline=True)



.. image:: Simple_Usage_files/Simple_Usage_19_0.png


--------------

Simple catalog search
---------------------

LLS from SDSS\_DR7 vs. IGMspec
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: python

    sdss_dr7_all = LLSSurvey.load_SDSS_DR7(sample='all')
    sdss_dr7_all


.. parsed-literal::

    SDSS-DR7: Loading LLS file /Users/xavier/local/Python/pyigm/pyigm/data/LLS/SDSS/lls_dr7_stat_LLS.fits.gz
    SDSS-DR7: Loading QSOs file /Users/xavier/local/Python/pyigm/pyigm/data/LLS/SDSS/lls_dr7_qsos_sn2050.fits.gz




.. parsed-literal::

    <IGMSurvey: nsys=1935, type=LLS, ref=SDSS-DR7, nsightlines=3759>



.. code:: python

    # Grab the coord
    lls_coord = sdss_dr7_all.coord
    lls_coord




.. parsed-literal::

    <SkyCoord (ICRS): (ra, dec) in deg
        [(339.61320833, 13.90905556), (160.36441667, 2.65569444),
         (167.247, 3.19108333), ..., (230.02470833, 23.66472222),
         (124.528625, 7.32227778), (150.86345833, 22.97005556)]>



.. code:: python

    # Match
    lls_ids = igmsp.match_coord(lls_coord)


.. parsed-literal::

    Your search yielded 1779 matches


.. code:: python

    # Show cat
    igmsp.show_cat(lls_ids)


.. parsed-literal::

    IGM_ID    RA      DEC     zem   flag_survey sig_zem  flag_zem
    ------ -------- -------- ------ ----------- ------- ---------
       881   0.7514  16.0077  3.662           3   0.002 BOSS_PCA 
       897   0.7640 -10.8641  3.638           3   0.001 BOSS_PCA 
      1680   1.4016  13.9971  3.709           3   0.002 BOSS_PCA 
      2168   1.8785  16.1257  3.519           3   0.002 BOSS_PCA 
      3248   2.8135  14.7672  4.986          19   0.003 BOSS_PCA 
      5192   4.5579  14.4155  4.216           3   0.001 BOSS_PCA 
      5225   4.5863  14.3143  3.938           3   0.004 BOSS_PCA 
      5482   4.8268  15.1032  4.145           3   0.003 BOSS_PCA 
      5617   4.9586  -0.6780  4.319           3   0.010 BOSS_PCA 
      5987   5.3335  15.8571  3.690           3   0.001 BOSS_PCA 
       ...      ...      ...    ...         ...     ...       ...
    374730 326.8571  -8.6430  4.588          18   0.001      SDSS
    374927 329.9334  -8.2762  3.755           2   0.001      SDSS
    375013 331.0089  -8.8288  4.448           2   0.001      SDSS
    375511 338.8385  -8.3576  4.425           2   0.000      SDSS
    375611 340.6793  -9.2622  4.214           2   0.001      SDSS
    375741 342.7194  -8.7667  3.869           2   0.001      SDSS
    375749 342.7879  -8.5271  3.883           2   0.001      SDSS
    375956 345.7561  -9.6585  3.455          14   0.001      SDSS
    376114 347.9044  -8.7360  3.745           2   0.001      SDSS
    376548 353.8939  -8.9943  3.683           2   0.001      SDSS
    376662 355.3705  -9.2696  4.125           2   0.002      SDSS
    Length = 1779 rows
    ----------
    Survey key:
        BOSS_DR12: 1
        GGG: 16
        HD-LLS_DR1: 8
        KODIAQ_DR1: 4
        SDSS_DR7: 2


.. code:: python

    # Grab GGG spec
    LLSspec, LLSmeta = igmsp.grab_spec('GGG', lls_ids)


.. parsed-literal::

    Staged 172 spectra totalling 0.0044032 Gb
    Loaded spectra


.. code:: python

    # Show the first spectrum
    LLSspec.plot(inline=True)



.. image:: Simple_Usage_files/Simple_Usage_28_0.png


.. code:: python

    # Show the 10th spectrum
    LLSspec.select=9
    LLSspec.plot(inline=True)



.. image:: Simple_Usage_files/Simple_Usage_29_0.png


--------------

Time Evolution
--------------

.. code:: python

    ggg_kodiaq = igmsp.cat['flag_survey'] == 18

.. code:: python

    ids_time = igmsp.cat['IGM_ID'][ggg_kodiaq]
    igmsp.show_cat(ids_time)


.. parsed-literal::

    IGM_ID    RA      DEC     zem   flag_survey sig_zem flag_zem
    ------ -------- -------- ------ ----------- ------- --------
    298834  21.2893 -10.7169  4.492          18   0.001     SDSS
    300926  52.8319  -7.6953  4.738          18   0.001     SDSS
    301049  54.6221   0.3656  5.032          18   0.001     SDSS
    331478 170.7229   0.8916  4.551          18   0.001     SDSS
    332948 173.1938  12.1505  5.167          18   0.002     SDSS
    349924 205.1677  28.2245  5.338          18   0.001     SDSS
    350040 205.4228  46.1862  5.023          18   0.001     SDSS
    355109 215.4375  35.2210  4.549          18   0.001     SDSS
    358845 222.7831   2.9377  4.481          18   0.001     SDSS
    366519 238.1793  25.8748  4.666          18   0.000     SDSS
    368827 243.6047  46.6747  5.313          18   0.001     SDSS
    368866 243.6960  20.9842  5.091          18   0.001     SDSS
    369021 244.0921   5.0244  4.872          18   0.001     SDSS
    373638 264.4370  58.4749  4.941          18   0.001     SDSS
    373960 314.3506  -0.5052  4.663          18   0.001     SDSS
    374730 326.8571  -8.6430  4.588          18   0.001     SDSS
    375405 337.1881  -7.9654  5.142          18   0.001     SDSS
    ----------
    Survey key:
        BOSS_DR12: 1
        GGG: 16
        HD-LLS_DR1: 8
        KODIAQ_DR1: 4
        SDSS_DR7: 2


.. code:: python

    meta = igmsp.grab_meta(['GGG','KODIAQ_DR1'], ids_time)

.. code:: python

    spec_time, meta_time = igmsp.grab_spec(['GGG','SDSS_DR7'], ids_time)


.. parsed-literal::

    Staged 34 spectra totalling 0.0008704 Gb
    Loaded spectra
    Staged 17 spectra totalling 0.001088 Gb
    Loaded spectra


.. code:: python

    spec_time




.. parsed-literal::

    [<XSpectrum1D: file=none, nspec=34, select=0, wvmin=4335.6 Angstrom, wvmax=7242.52 Angstrom>,
     <XSpectrum1D: file=none, nspec=17, select=0, wvmin=3800.14 Angstrom, wvmax=9206.62 Angstrom>]



Plot both
~~~~~~~~~

.. code:: python

    spec_time[0].plot(plot_two=spec_time[1],inline=True, scale_two=0.6)



.. image:: Simple_Usage_files/Simple_Usage_38_0.png


--------------

Pairs
-----

QPQ8 like
~~~~~~~~~

Query on separation (angular and redshift)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    ID_fg, ID_bg = igmsp.pairs(30*u.arcsec, 3000.*u.km/u.s)

.. code:: python

    igmsp.show_cat(ID_fg)


.. parsed-literal::

    IGM_ID    RA      DEC     zem   flag_survey sig_zem  flag_zem
    ------ -------- -------- ------ ----------- ------- ---------
       293   0.2461  28.3758  0.985           1  -1.000 BOSS_PCA 
       434   0.3669   8.6782  2.794           1   0.003 BOSS_PCA 
       851   0.7280  21.6517  1.877           1  -1.000 BOSS_PCA 
      1006   0.8555  13.4381  0.628           1  -1.000 BOSS_PCA 
      1638   1.3680  25.7650  2.545           1  -1.000 BOSS_PCA 
      2491   2.1692  17.1640  1.380           1  -1.000 BOSS_PCA 
      2749   2.3791  17.4591  0.665           1  -1.000 BOSS_PCA 
      2831   2.4458  12.0703  2.254           1   0.002 BOSS_PCA 
      2873   2.4857  26.2747  1.997           1   0.001 BOSS_PCA 
      2941   2.5424  32.9975  2.043           1  -1.000 BOSS_PCA 
       ...      ...      ...    ...         ...     ...       ...
    374776 327.5874   0.9728  1.011           2   0.001      SDSS
    374886 329.2297  -0.0669  1.657           2   0.014      SDSS
    374986 330.7025  12.6126  2.073           2   0.002      SDSS
    375051 331.6001  11.5284  0.400           2   0.001      SDSS
    375227 334.5278   0.8732  1.273           2   0.002      SDSS
    375376 336.8406  -1.1105  1.363           2   0.002      SDSS
    375390 336.9536  12.2599  0.978           2   0.001      SDSS
    375725 342.4622  -0.7526  1.356           2   0.002      SDSS
    376169 348.6047  -1.1508  1.640           2   0.002      SDSS
    376362 351.1376  15.4438  0.295           2   0.001      SDSS
    376932 359.7813  13.7919  0.247           2   0.001      SDSS
    Length = 1272 rows
    ----------
    Survey key:
        BOSS_DR12: 1
        GGG: 16
        HD-LLS_DR1: 8
        KODIAQ_DR1: 4
        SDSS_DR7: 2


Check for high dispersion spectrum in b/g QSOs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    highdisp = igmsp.cutid_on_surveys(['HD-LLS_DR1', 'KODIAQ_DR1'], ID_bg)

.. code:: python

    igmsp.show_cat(ID_bg[highdisp])


.. parsed-literal::

    IGM_ID    RA      DEC     zem   flag_survey sig_zem flag_zem
    ------ -------- -------- ------ ----------- ------- --------
    376957  36.4785   0.9144  2.975           4   0.000   SIMBAD
    ----------
    Survey key:
        BOSS_DR12: 1
        GGG: 16
        HD-LLS_DR1: 8
        KODIAQ_DR1: 4
        SDSS_DR7: 2


.. code:: python

    igmsp.show_cat(ID_fg[highdisp])


.. parsed-literal::

    IGM_ID    RA      DEC     zem   flag_survey sig_zem  flag_zem
    ------ -------- -------- ------ ----------- ------- ---------
     35491  36.4786   0.9213  1.770           1   0.000 BOSS_PCA 
    ----------
    Survey key:
        BOSS_DR12: 1
        GGG: 16
        HD-LLS_DR1: 8
        KODIAQ_DR1: 4
        SDSS_DR7: 2


