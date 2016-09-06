
2QZ [v1.1]
==========

:download:`Download <examples/2QZ.ipynb>` this notebook.

.. code:: python

    # imports

Data
----

::

    Used this website (on Croom recommendation):  http://www.2dfquasar.org/Spec_Cat/2qzsearch2.html
    Restricted to zmin>0.05

Catalog
-------

::

    See file:///Users/xavier/Raw_IGMspec/2dF/2df/html/catalogue.html

.. code:: python

    #catfil = os.getenv('RAW_IGMSPEC')+'/2dF/2df/cat/2QZ_6QZ_pubcat.txt'
    catfil = os.getenv('RAW_IGMSPEC')+'/2dF/2QZ393524355693.out'

.. code:: python

    cat = Table.read(catfil,format='ascii')

.. code:: python

    cat[0:3]




.. raw:: html

    &lt;Table length=3&gt;
    <table id="table4556544912">
    <thead><tr><th>col1</th><th>col2</th><th>col3</th><th>col4</th><th>col5</th><th>col6</th><th>col7</th><th>col8</th><th>col9</th><th>col10</th><th>col11</th><th>col12</th><th>col13</th><th>col14</th><th>col15</th><th>col16</th><th>col17</th><th>col18</th><th>col19</th><th>col20</th><th>col21</th><th>col22</th><th>col23</th><th>col24</th><th>col25</th><th>col26</th><th>col27</th><th>col28</th><th>col29</th><th>col30</th><th>col31</th><th>col32</th><th>col33</th><th>col34</th><th>col35</th><th>col36</th><th>col37</th><th>col38</th><th>col39</th><th>col40</th><th>col41</th><th>col42</th><th>col43</th><th>col44</th><th>col45</th></tr></thead>
    <thead><tr><th>str16</th><th>int64</th><th>int64</th><th>float64</th><th>int64</th><th>int64</th><th>float64</th><th>int64</th><th>str10</th><th>str25</th><th>int64</th><th>int64</th><th>float64</th><th>int64</th><th>int64</th><th>float64</th><th>int64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>int64</th><th>float64</th><th>int64</th><th>str9</th><th>int64</th><th>float64</th><th>int64</th><th>float64</th><th>float64</th><th>int64</th><th>str9</th><th>float64</th><th>float64</th><th>int64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>str20</th><th>str20</th></tr></thead>
    <tr><td>J120310.0-025755</td><td>12</td><td>3</td><td>10.09</td><td>-2</td><td>57</td><td>54.2</td><td>8573</td><td>TQN178_111</td><td>N_178_244</td><td>12</td><td>0</td><td>36.29</td><td>-2</td><td>41</td><td>12.8</td><td>859</td><td>-1051.86</td><td>17295.65</td><td>3.14423174</td><td>-0.04689506</td><td>20.224</td><td>-0.229</td><td>-0.612</td><td>3</td><td>0.9174</td><td>11</td><td>QSO</td><td>20010524</td><td>2441.0</td><td>168</td><td>3.89</td><td>0.9113</td><td>11</td><td>QSO</td><td>20000624.0</td><td>2441.0</td><td>167</td><td>0.41</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.02606</td><td>-</td><td>-</td></tr>
    <tr><td>J120427.0-025447</td><td>12</td><td>4</td><td>27.1</td><td>-2</td><td>54</td><td>46.7</td><td>8601</td><td>TQN178_114</td><td>N_178_179_244</td><td>12</td><td>1</td><td>53.29</td><td>-2</td><td>38</td><td>5.3</td><td>859</td><td>-3200.31</td><td>16946.83</td><td>3.14983133</td><td>-0.04598603</td><td>20.621</td><td>-0.634</td><td>1.112</td><td>3</td><td>0.3276</td><td>11</td><td>QSO</td><td>20010524</td><td>2441.0</td><td>184</td><td>5.97</td><td>0.0</td><td>33</td><td>??</td><td>20000624.0</td><td>2441.0</td><td>190</td><td>0.71</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.02662</td><td>-</td><td>-</td></tr>
    <tr><td>J120434.3-025032</td><td>12</td><td>4</td><td>34.32</td><td>-2</td><td>50</td><td>31.3</td><td>8631</td><td>TQN178_117</td><td>N_178_179_244</td><td>12</td><td>2</td><td>0.51</td><td>-2</td><td>33</td><td>49.9</td><td>859</td><td>-3402.12</td><td>16471.4</td><td>3.15035639</td><td>-0.04474782</td><td>20.798</td><td>-0.319</td><td>0.406</td><td>3</td><td>0.4312</td><td>11</td><td>NELG</td><td>20010524</td><td>2441.0</td><td>190</td><td>3.49</td><td>0.431</td><td>11</td><td>NELG</td><td>20000624.0</td><td>2441.0</td><td>189</td><td>0.26</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.02582</td><td>-</td><td>-</td></tr>
    </table>



Rename columns
~~~~~~~~~~~~~~

.. code:: python

    clms = ['Name', 'RAh00', 'RAm00', 'RAs00', 'DECd00', 'DECm00', 'DECs00',
           'ID','cat_name', 'Sector', 'RAh50', 'RAm50', 'RAs50', 'DECd50', 'DECm50', 'DECs50',
           'UKST', 'XAPM','YAPM','RA50','DEC50','bj','u-b','b-r','Nobs',
           'z1','q1','ID1','date1','fld1','fiber1','SN1',
           'z2','q2','ID2','date2','fld2','fiber2','SN2',
            'zprev','rflux','Xray','EBV','comm1','comm2']

.. code:: python

    for ii in range(1,46):
        cat.rename_column('col{:d}'.format(ii), clms[ii-1])

.. code:: python

    cat[0:3]




.. raw:: html

    &lt;Table length=3&gt;
    <table id="table4748415120">
    <thead><tr><th>Name</th><th>RAh00</th><th>RAm00</th><th>RAs00</th><th>DECd00</th><th>DECm00</th><th>DECs00</th><th>ID</th><th>cat_name</th><th>Sector</th><th>RAh50</th><th>RAm50</th><th>RAs50</th><th>DECd50</th><th>DECm50</th><th>DECs50</th><th>UKST</th><th>XAPM</th><th>YAPM</th><th>RA50</th><th>DEC50</th><th>bj</th><th>u-b</th><th>b-r</th><th>Nobs</th><th>z1</th><th>q1</th><th>ID1</th><th>date1</th><th>fld1</th><th>fiber1</th><th>SN1</th><th>z2</th><th>q2</th><th>ID2</th><th>date2</th><th>fld2</th><th>fiber2</th><th>SN2</th><th>zprev</th><th>rflux</th><th>Xray</th><th>EBV</th><th>comm1</th><th>comm2</th></tr></thead>
    <thead><tr><th>str16</th><th>int64</th><th>int64</th><th>float64</th><th>int64</th><th>int64</th><th>float64</th><th>int64</th><th>str10</th><th>str25</th><th>int64</th><th>int64</th><th>float64</th><th>int64</th><th>int64</th><th>float64</th><th>int64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>int64</th><th>float64</th><th>int64</th><th>str9</th><th>int64</th><th>float64</th><th>int64</th><th>float64</th><th>float64</th><th>int64</th><th>str9</th><th>float64</th><th>float64</th><th>int64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>float64</th><th>str20</th><th>str20</th></tr></thead>
    <tr><td>J120310.0-025755</td><td>12</td><td>3</td><td>10.09</td><td>-2</td><td>57</td><td>54.2</td><td>8573</td><td>TQN178_111</td><td>N_178_244</td><td>12</td><td>0</td><td>36.29</td><td>-2</td><td>41</td><td>12.8</td><td>859</td><td>-1051.86</td><td>17295.65</td><td>3.14423174</td><td>-0.04689506</td><td>20.224</td><td>-0.229</td><td>-0.612</td><td>3</td><td>0.9174</td><td>11</td><td>QSO</td><td>20010524</td><td>2441.0</td><td>168</td><td>3.89</td><td>0.9113</td><td>11</td><td>QSO</td><td>20000624.0</td><td>2441.0</td><td>167</td><td>0.41</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.02606</td><td>-</td><td>-</td></tr>
    <tr><td>J120427.0-025447</td><td>12</td><td>4</td><td>27.1</td><td>-2</td><td>54</td><td>46.7</td><td>8601</td><td>TQN178_114</td><td>N_178_179_244</td><td>12</td><td>1</td><td>53.29</td><td>-2</td><td>38</td><td>5.3</td><td>859</td><td>-3200.31</td><td>16946.83</td><td>3.14983133</td><td>-0.04598603</td><td>20.621</td><td>-0.634</td><td>1.112</td><td>3</td><td>0.3276</td><td>11</td><td>QSO</td><td>20010524</td><td>2441.0</td><td>184</td><td>5.97</td><td>0.0</td><td>33</td><td>??</td><td>20000624.0</td><td>2441.0</td><td>190</td><td>0.71</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.02662</td><td>-</td><td>-</td></tr>
    <tr><td>J120434.3-025032</td><td>12</td><td>4</td><td>34.32</td><td>-2</td><td>50</td><td>31.3</td><td>8631</td><td>TQN178_117</td><td>N_178_179_244</td><td>12</td><td>2</td><td>0.51</td><td>-2</td><td>33</td><td>49.9</td><td>859</td><td>-3402.12</td><td>16471.4</td><td>3.15035639</td><td>-0.04474782</td><td>20.798</td><td>-0.319</td><td>0.406</td><td>3</td><td>0.4312</td><td>11</td><td>NELG</td><td>20010524</td><td>2441.0</td><td>190</td><td>3.49</td><td>0.431</td><td>11</td><td>NELG</td><td>20000624.0</td><td>2441.0</td><td>189</td><td>0.26</td><td>0.0</td><td>0.0</td><td>0.0</td><td>0.02582</td><td>-</td><td>-</td></tr>
    </table>



FITS
----

::

    Long or short?
    Over-wrote CD files with those from the website
    tar -xvf 2QZspec393524355693.tar

.. code:: python

    exfits = os.getenv('RAW_IGMSPEC')+'/2dF/2df/fits/ra12_13/J120310.0-025755a.fits.gz'

.. code:: python

    hdu = fits.open(exfits)
    hdu.info()


.. parsed-literal::

    Filename: /u/xavier/Raw_IGMspec/2dF/2df/fits/ra12_13/J120310.0-025755a.fits.gz
    No.    Name         Type      Cards   Dimensions   Format
    0    PRIMARY     PrimaryHDU      70   (1024,)      int16 (rescales to float32)   
    1    BADPIX      ImageHDU        15   (1024,)      uint8   
    2    VARIANCE    ImageHDU        17   (1024,)      int16 (rescales to float32)   


Resolution
----------

2dF
~~~

::

     arXiv:astro-ph/9804079

.. code:: python

    R = 5000./8.6
    R




.. parsed-literal::

    581.3953488372093



6dF
~~~

::

    Dead link:  https://www.aao.gov.au/ukst/6df.html
    http://adsabs.harvard.edu/abs/2001MNRAS.322L..29C

Time
----

.. code:: python

    from astropy.time import Time

.. code:: python

    t = Time(['2000-06-24'], format='iso',out_subfmt='date')  # Fixes to YYYY-MM-DD
    t.iso




.. parsed-literal::

    array(['2000-06-24'], 
          dtype='|S10')



