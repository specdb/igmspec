
Simple Scripts with IGMspec (v1.2)
==================================

:download:`Download <examples/Simple_Scripts.ipynb>` this notebook.

.. code:: python

    # imports

Downloading the DataBase
------------------------

After installing igmspec, you can grab the latest (or any previous)
version of the DB with the *get\_igmspec* script.

::

    usage: get_igmspec [-h] [-v VERSION]

    Grab the IGMspec DB

    optional arguments:
      -h, --help            show this help message and exit
      -v VERSION, --version VERSION
                            DB version to generate

Examples
~~~~~~~~

v01
^^^

::

    get_igmspec -v v01

--------------

Plot
----

::

    wolverine.local> plot_igmspec -h
    usage: plot_igmspec [-h] [--toler TOLER] [--meta] [--survey SURVEY]
                        [--select SELECT] [--mplot MPLOT]
                        coord

    plot_igmspec script v0.1

    positional arguments:
      coord            Coordinates, e.g. J081240+320808

    optional arguments:
      -h, --help       show this help message and exit
      --toler TOLER    Maximum offset in arcsec [default=5.]
      --meta           Show meta data? [default: True]
      --survey SURVEY  Name of Survey to use
      --select SELECT  Name of Survey to use [default: 0]
      --mplot MPLOT    Use simple matplotlib plot [default: False]

Examples
~~~~~~~~

FJ0812+32
^^^^^^^^^

::

    plot_igmspec J081241+320808 --survey KODIAQ_DR1

J001115.23+144601.8
^^^^^^^^^^^^^^^^^^^

::

    plot_igmspec J001115.23+144601.8

--------------

SDSS/BOSS Database
------------------

::

    profx.ucolick.org> sdss_igmspec -h
    usage: sdss_igmspec [-h] [-s SURVEY] [--select SELECT] [-p] plate fiberid

    sdss_igmspec script v0.1

    positional arguments:
      plate                 Plate
      fiberid               FiberID

    optional arguments:
      -h, --help            show this help message and exit
      -s SURVEY, --survey SURVEY
                            Name of Survey to use (BOSS_DR12 or SDSS_DR7)
      --select SELECT       Index of spectrum to plot (when multiple exist)
      -p, --plot            Plot with lt_xspec

Example
~~~~~~~

::

    UNIX> sdss_igmspec 434 555 -p

