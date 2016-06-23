""" Module for IgmSpec Class
"""
from __future__ import print_function, absolute_import, division, unicode_literals

#import numpy as np
import pdb


from igmspec.query_catalog import QueryCatalog
from igmspec.interface_db import InterfaceDB


class IgmClass(object):
    """ An uber-class for the package
    Ideally one-stop-shopping for most consumers

    Parameters
    ----------

    Attributes
    ----------
    qcat : QueryCatalog
    idb : InterfaceDB
    """

    def __init__(self, **kwargs):
        """
        """
        # Init
        self.qcat = QueryCatalog(**kwargs)
        self.idb = InterfaceDB(**kwargs)
        # Checks
        assert self.idb.db_file == self.qcat.db_file
        assert self.idb.surveys == self.qcat.surveys

    def __getattr__(self, k):
        """ Overload attributes using the underlying classes

        Parameters
        ----------
        k

        Returns
        -------

        """
        # Try Self first
        try:
            return getattr(self, k)
        except AttributeError:
            # DB
            try:
                return getattr(self.idb, k)
            except AttributeError:
                # Try qcat last
                return getattr(self.qcat, k)

    def __repr__(self):
        txt = '<{:s}:  IGM_file={:s} with {:d} sources\n'.format(self.__class__.__name__,
                                            self.db_file, len(self.cat))
        # Surveys
        txt += '   Loaded surveys are {} \n'.format(self.surveys)
        txt += '>'
        return (txt)
