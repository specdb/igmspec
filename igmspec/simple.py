""" Module for simple igmspec routines
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os, glob
import warnings
import pdb


from igmspec import query_catalog as iqcat
from igmspec import cat_utils as icu
from igmspec import interface_db as igidb

from astropy import units as u


