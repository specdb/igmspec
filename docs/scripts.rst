.. highlight:: rest

***************
igmspec Scripts
***************

This file summarizes the igmspec scripts
(used outside Python).  These are installed
within your standard Python script path (e.g.
~/anaconda/bin).

Notebooks
=========

.. toctree::
   :maxdepth: 1

       Simple Scripts <Simple_Scripts>

plot_igmspec
============

Plot a spectrum at the given coordinate.  One can
restrict the survey used and/or select the desired
spectrum from the available list.  By default, the
XSpecGui gui from linetools is called to display
the spectrum.   Here is the help::

   $ plot_igmspec -h
    usage: plot_igmspec [-h] [--tol TOL] [--meta] [--survey SURVEY]
                        [--select SELECT] [--mplot MPLOT]
                        coord

    plot_igmspec script v0.1

    positional arguments:
      coord            Coordinates, e.g. J081240+320808

    optional arguments:
      -h, --help       show this help message and exit
      --tol TOL        Maximum offset in arcsec [default=5.]
      --meta           Show meta data? [default: True]
      --survey SURVEY  Name of Survey to use
      --select SELECT  Name of Survey to use [default: 0]
      --mplot MPLOT    Use simple matplotlib plot [default: False]

