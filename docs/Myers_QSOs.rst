
Ingesting of Myers' QSO catalog (v1.1)
======================================

:download:`Download <examples/Myers_QSOs.ipynb>` this notebook.

.. code:: python

    # imports
    from igmspec.ingest import myers as imy

Grabbing the catalog
--------------------

http://faraday.uwyo.edu/~admyers/scratch/forjoe/GTR-ADM-QSO-master-wvcv.fits.gz

::

    by JXP on 5 August 2016
    continuing on 13 Aug 2016

Inspect
-------

.. code:: python

    ADM_file = os.getenv('RAW_IGMSPEC')+'/Myers/GTR-ADM-QSO-master-wvcv.fits.gz'
    ADM_qso = Table.read(ADM_file)

.. code:: python

    ADM_qso[0:5]




.. raw:: html

    &lt;Table length=5&gt;
    <table id="table4553456656">
    <thead><tr><th>RA</th><th>DEC</th><th>ZBEST [19]</th><th>SOURCEBIT</th><th>SDSS_UFLG</th><th>AGES_QSO</th><th>AGES_CODE06</th><th>KDE_ZPHOTLO</th><th>KDE_ZPHOTHI</th><th>KDE_ZPHOTPROB</th><th>KDE_LOWZORUVX</th><th>XDQSOZ_PEAKPROB</th><th>XDQSOZ_PEAKFWHM</th><th>XDQSOZ_NPEAKS</th><th>YAPERMAG3</th><th>JAPERMAG3</th><th>HAPERMAG3</th><th>KSAPERMAG3</th><th>YAPERMAG3ERR</th><th>JAPERMAG3ERR</th><th>HAPERMAG3ERR</th><th>KSAPERMAG3ERR</th><th>ZPHOTMINJHK</th><th>ZPHOTBESTJHK</th><th>ZPHOTMAXJHK</th><th>ZPHOTPROBJHK</th></tr></thead>
    <thead><tr><th>float32</th><th>float32</th><th>float32</th><th>int32</th><th>int16</th><th>int16</th><th>int32</th><th>float32</th><th>float32</th><th>float32</th><th>int16</th><th>float32</th><th>float32</th><th>int16</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th><th>float32</th></tr></thead>
    <tr><td>0.027228</td><td>0.515341</td><td>1.82383 .. 0.0</td><td>34305</td><td>0</td><td>0</td><td>0</td><td>1.62</td><td>2.11</td><td>0.719</td><td>1</td><td>0.983058</td><td>0.1924</td><td>1</td><td>-9.99999e+08</td><td>19.1973</td><td>18.8874</td><td>18.3277</td><td>-9.99999e+08</td><td>0.0942081</td><td>0.187276</td><td>0.172265</td><td>1.45</td><td>1.825</td><td>1.95</td><td>0.998998</td></tr>
    <tr><td>0.0339</td><td>0.276301</td><td>1.83638 .. 0.0</td><td>33281</td><td>0</td><td>0</td><td>0</td><td>1.76</td><td>2.14</td><td>0.946</td><td>1</td><td>0.0</td><td>0.0</td><td>0</td><td>-9.99999e+08</td><td>18.3689</td><td>17.716</td><td>17.1701</td><td>-9.99999e+08</td><td>0.0459005</td><td>0.0655821</td><td>0.0615399</td><td>1.7</td><td>1.875</td><td>1.95</td><td>0.996776</td></tr>
    <tr><td>0.038604</td><td>15.2985</td><td>1.19748 .. 0.0</td><td>139777</td><td>1</td><td>0</td><td>0</td><td>1.1</td><td>1.53</td><td>0.912</td><td>1</td><td>0.0</td><td>0.0</td><td>0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr>
    <tr><td>0.039089</td><td>13.9384</td><td>2.2342 .. 2.239</td><td>403073</td><td>1</td><td>0</td><td>0</td><td>0.93</td><td>1.44</td><td>0.634</td><td>1</td><td>0.827387</td><td>0.1976</td><td>3</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr>
    <tr><td>0.039271</td><td>-10.4644</td><td>1.85301 .. 0.0</td><td>132609</td><td>0</td><td>0</td><td>0</td><td>1.45</td><td>2.06</td><td>0.739</td><td>1</td><td>0.939097</td><td>0.3484</td><td>2</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.0</td></tr>
    </table>



.. code:: python

    hdu = fits.open(ADM_file)

.. code:: python

    head0 = hdu[0].header
    head0




.. parsed-literal::

    SIMPLE  =                    T /Dummy Created by MWRFITS v1.8                   
    BITPIX  =                    8 /Dummy primary header created by MWRFITS         
    NAXIS   =                    0 /No data is associated with this header          
    EXTEND  =                    T /Extensions may (will!) be present               



.. code:: python

    head1 = hdu[1].header
    head1['DATE']




.. parsed-literal::

    'Sat Oct 10 07:10:38 2015'



Igmspec
-------

Load
~~~~

.. code:: python

    reload(imy)
    ADM_qso, data = imy.load()

Redshifts
~~~~~~~~~

.. code:: python

    imy.zbest_myers(ADM_qso)

Check
^^^^^

.. code:: python

    myers_binary = [2**0, 2**7, 2**1, 2**2, 2**3, 2**4, 2**5, 2**6, 2**8, 2**11,
                          2**12, 2**13, 2**14, 2**16, 2**17, 2**18]

.. code:: python

    tqz = np.where(ADM_qso['ZEM_SOURCE'] == '2QZ')[0]
    tqz[0]




.. parsed-literal::

    105986



.. code:: python

    ADM_qso[105986]['ZBEST']




.. parsed-literal::

    array([ 0.        ,  1.14339995,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
            0.        ,  0.        ,  1.14300001,  0.        ], dtype=float32)



.. code:: python

    ADM_qso[105986]['SOURCEBIT'] & myers_binary




.. parsed-literal::

    array([     0,      0,      2,      0,      0,      0,      0,      0,
                0,      0,      0,      0,      0,      0, 131072,      0])



.. code:: python

    ADM_qso[105986]['ZEM']




.. parsed-literal::

    1.1434



Cut
~~~

.. code:: python

    # Cut down
    ztrim = (ADM_qso['ZEM'] >= 0.1) & (ADM_qso['ZEM'] <= 7.0)
    coordtrim = (ADM_qso['RA'] >= 0.0) & (ADM_qso['RA'] <= 360.0) & (np.abs(
            ADM_qso['DEC']) <= 90.0)
    keep = ztrim & coordtrim
    #ADM_qso = ADM_qso[keep]

.. code:: python

    np.sum(keep)




.. parsed-literal::

    513996



