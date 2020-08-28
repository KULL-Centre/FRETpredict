# -*- coding: utf-8 -*-
"""
Rotamer library handling
========================

:mod:`rotamers.library` contains the data (:data:`LIBRARIES`) to load
a rotamer library, represented by a :class:`RotamerLibrary`.

"""

import MDAnalysis

import logging
logger = logging.getLogger("MDAnalysis.app")

import numpy as np
import os.path
import pkg_resources
import yaml

#: Name of the directory in the package that contains the library data.
LIBDIR = "./lib"

def find_file(filename, pkglibdir=LIBDIR):
    """
    Function to find path using MDAnalysis/

    Args:
        filename:
        pkglibdir:

    Returns:

    """
    if os.path.exists(filename):
        return MDAnalysis.lib.util.realpath(filename)
    return pkg_resources.resource_filename(__name__, os.path.join(pkglibdir, filename))


with open(find_file('libraries.yml'), 'r') as yaml_file:
    #: Registry of libraries, indexed by name.
    #: Takes the format = {name: {topology, data, author, license, citation}, ...}
    LIBRARIES = yaml.load(yaml_file, Loader=yaml.FullLoader)

class RotamerLibrary(object):
    """
    Rotamer library

    The library makes available the attributes :attr:`data` and :attr:`weights`.

    Attributes:
        coord (:py:class:`numpy.ndarray`): Array containing the relative coordinates of each rotamer.

        weights (:py:class:`numpy.ndarray`): Array containing the population of each rotamer.

        name (:py:class:`str`): Name of the library.

        lib (:py:class:`dict`): Dictionary containing the file names and meta data for the library :attr:`name`.

    """

    def __init__(self, name):
        """

        Args:
            name (:py:class:`str`): name of the library (must exist in the registry of libraries, :data:`LIBRARIES`)
        """
        self.name = name
        self.lib = {}
        try:
            self.lib.update(LIBRARIES[name])  # make a copy
        except KeyError:
            raise ValueError("No rotamer library with name {0} known: must be one of {1}".format(name,
                                                                                                 list(LIBRARIES.keys())))
        logger.info("Using rotamer library '{0}' by {1[author]}".format(self.name, self.lib))
        logger.info("Please cite: {0[citation]}".format(self.lib))
        # adjust paths
        for k in 'data', 'topology':
            self.lib[k] = find_file(self.lib[k])
        logger.debug("[rotamers] ensemble = {0[data]} with topology = {0[topology]}".format(self.lib))
        logger.debug("[rotamers] populations = {0[data]}".format(self.lib))

        self.top = MDAnalysis.Universe(self.lib['topology'])
        data = np.loadtxt(self.lib['data'], dtype='float32',usecols=(2,3,4,5))
        self.coord = data.reshape((data.shape[0] // len(self.top.atoms),len(self.top.atoms), 4))[:,:,:3].swapaxes(0, 1)
        self.weights = data.reshape((data.shape[0] // len(self.top.atoms),len(self.top.atoms), 4))[:,0,3]

    def __repr__(self):
        """

        Returns:
            (:py:class:`str`): Name of the library in use.
        """
        return "<RotamerLibrary '{0}' by {1} with {2} rotamers>".format(self.name, self.lib['author'],
                                                                        len(self.weights)-2)
