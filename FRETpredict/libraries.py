# -*- coding: utf-8 -*-
"""

Rotamer library handling
========================

:mod:`rotamers.library` contains the data (:data:`LIBRARIES`) to load
a rotamer library, represented by a :class:`RotamerLibrary`.

"""

import MDAnalysis
import logging
import numpy as np
import os.path
import pkg_resources
import yaml

# Create logger
#logger = logging.getLogger("MDAnalysis.app")
#logger.setLevel(logging.DEBUG)
#logger.root.setLevel(logging.DEBUG)

#: Name of the directory in the package that contains the library data.
LIBDIR = "lib"

def find_file(filename, pkglibdir=LIBDIR):

    """

    Function to find path using MDAnalysis

    Parameters
    ==========

        filename: str
            Name of the file to find

        pkglibdir: str
            Name of the directory in the package that contains the library data

    Returns:

        filename path

    """

    # If path already exists, return path rooted at /
    if os.path.exists(filename):
        return MDAnalysis.lib.util.realpath(filename)

    return pkg_resources.resource_filename(__name__, os.path.join(pkglibdir, filename))


# Registry of libraries, indexed by name.
# Takes the format = {name: {topology, data, author, license, citation}, ...}
with open(find_file('libraries.yml'), 'r') as yaml_file:
    LIBRARIES = yaml.load(yaml_file, Loader=yaml.FullLoader)


class RotamerLibrary(object):

    """

    Rotamer library
    ===============

    The library makes available the attributes `data` and `weights`.

    Attributes
    ==========

        coord: numpy.ndarray
            Array containing the relative coordinates of each rotamer.

        weights: numpy.ndarray
            Array containing the population of each rotamer.

        name: str
            Name of the library.

        lib: dict
            Dictionary containing the file names and metadata for the library attribute `name`.

    """

    def __init__(self, name):

        """

        Parameters
        ==========

            name: str
                name of the library (must exist in the registry of libraries, `LIBRARIES`)

        """

        # Keep the name of the library without "cutoff"
        self.name = name.split(' cutoff')[0]

        self.lib = {}

        # Make a copy of the library
        try:
            self.lib.update(LIBRARIES[self.name])
        # No rotamer library found called 'name'
        except KeyError:
            raise ValueError("No rotamer library with name {0} known: must be one of {1}".format(
                name, list(LIBRARIES.keys())))

        # Obtain filename path
        self.lib['filename'] = find_file(self.lib['filename'][:-2] + name.split(' cutoff')[1])

        # Print logging information
        #logger.debug("Using rotamer library '{0}' by {1[author]}".format(self.name, self.lib))
        #logger.debug("Please cite: {0[citation]}".format(self.lib))
        #logger.debug("[rotamers] ensemble {:s} with topology {:s}.pdb".format(self.lib['filename'],self.lib['filename'].split('_cutoff')[0]))
        #logger.debug("[rotamers] populations {:s}".format(self.lib['filename'] + '_weights.txt'))

        # If trajectory is not found
        if not os.path.isfile(self.lib['filename'] + '.dcd'):
            raise ValueError("No trajectory named {0}.dcd".format(self.lib['filename']))

        # If weights are not found
        if not os.path.isfile(self.lib['filename'] + '_weights.txt'):
            raise ValueError("No file named {0}_weights.txt".format(self.lib['filename'] + '_weights.txt'))

        # Rotamers parameters
        self.top = MDAnalysis.Universe(self.lib['filename'].split('_cutoff')[0] + '.pdb')
        self.mu = self.lib['mu']
        self.r = self.lib['r']
        self.positive = self.lib['positive']
        self.negative = self.lib['negative']

        # Obtain chromophore trajectory
        traj = MDAnalysis.Universe(self.lib['filename'].split('_cutoff')[0] + '.pdb',
                                   self.lib['filename'] + '.dcd')

        # TODO: Obtain chromophore trajectory from .XTC or other formats (selectable by user)

        # Extract coordinates from DCD trajectory
        self.coord = traj.trajectory.timeseries(traj.atoms)

        # Extract weights from weight file
        self.weights = np.loadtxt(self.lib['filename'] + '_weights.txt')
        self.weights /= np.sum(self.weights)

    def __repr__(self):

        """

        Returns
        =======

            String with name of the library in use.

        """

        return "<RotamerLibrary '{0}' by {1} with {2} rotamers>".format(self.name, self.lib['author'], len(self.weights) - 2)
