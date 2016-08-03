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
    usage: plot_igmspec [-h] [--tol TOL] [--meta] [-s SURVEY] [--select SELECT]
                        [--mplot MPLOT]
                        coord

    plot_igmspec script v0.2

    positional arguments:
      coord                 Coordinates, e.g. J081240+320808

    optional arguments:
      -h, --help            show this help message and exit
      --tol TOL             Maximum offset in arcsec [default=5.]
      --meta                Show meta data? [default: True]
      -s SURVEY, --survey SURVEY
                            Name of Survey to use
      --select SELECT       Index of spectrum to plot (when multiple exist)
      --mplot MPLOT         Use simple matplotlib plot [default: False]

sdss_igmspec
============

Grab data from the SDSS/BOSS survey.  Current implementation requires
plate/fiber input.  Here is the help::

   $profx.ucolick.org> sdss_igmspec -h
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
