# -*- coding: utf-8 -*-
"""
PREPrediction Class
-------------------

Class to perform Paramagnetic Relaxation Enhancement prediction, employing the Model-Free Solomon-Bloembergen equation.

"""

import os
import numpy as np
import MDAnalysis
from MDAnalysis.coordinates.memory import MemoryReader
from DEERPREdict.lennardjones import vdw, p_Rmin2, eps
import DEERPREdict.libraries as libraries
import logging
import scipy.special as special
from scipy.spatial.distance import cdist

class Operations(object):
    """Calculation of the distance profile between a probe and backbone amide."""

    def __init__(self, protein, **kwargs):
        """
        Args:
            protein (:py:class:`MDAnalysis.core.universe.Universe`): trajectory
        """
        self.protein = protein
        self.libname = kwargs.get('libname', 'MTSSL 175K X1X2')
        self.lib = libraries.RotamerLibrary(self.libname)
        self.temp = kwargs.get('temperature', 300)
        self.z_cutoff = kwargs.get('z_cutoff', 0.05)
        self.ign_H = kwargs.get('ign_H', True)
        self.chains = kwargs.get('chains', [None,None])

    def precalculate_rotamer(self, residue, chain):
        residue_sel = "resid {:d}".format(residue)
        if type(chain) == str:
            residue_sel += " and segid {:s}".format(chain)
        prot_Ca = self.protein.select_atoms('protein and name CA and '+residue_sel)
        prot_Co = self.protein.select_atoms('protein and name C and '+residue_sel)
        prot_N = self.protein.select_atoms('protein and name N and '+residue_sel)
        probe_coords = np.zeros((len(self.lib.top.atoms),1, 3))
        universe = MDAnalysis.Universe(self.lib.top.filename, probe_coords, format=MemoryReader, order='afc')
        return universe, (prot_Ca, prot_Co, prot_N), residue_sel
        
    def rotamer_placement(self, universe, prot_atoms):
        prot_Ca, prot_Co, prot_N = prot_atoms
        offset = prot_Ca.positions.copy()
        Ca_coords = prot_Ca.positions - offset
        Co_coords = prot_Co.positions - offset
        N_coords = prot_N.positions - offset
        x_vector = N_coords - Ca_coords
        x_vector /= np.linalg.norm(x_vector)
        yt_vector = Co_coords - Ca_coords
        yt_vector /= np.linalg.norm(yt_vector)
        z_vector = np.cross(x_vector, yt_vector)
        z_vector /= np.linalg.norm(z_vector)
        y_vector = np.cross(z_vector, x_vector)
        rotation = np.vstack([x_vector, y_vector, z_vector])
        probe_coords = np.tensordot(self.lib.coord,rotation,axes=([2],[0])) + offset
        universe.load_new(probe_coords, format=MemoryReader, order='afc')
        #save aligned rotamers
        #mtssl = universe.select_atoms("all")
        #with MDAnalysis.Writer("mtssl.pdb", mtssl.n_atoms) as W:
        #    for ts in universe.trajectory:
        #        W.write(mtssl)
        return universe

    def lj_calculation(self, fitted_rotamers, residue_sel):
        gas_un = 1.9858775e-3 # CHARMM, in kcal/mol*K
        if self.ign_H:
            proteinNotSite = self.protein.select_atoms("protein and not type H and not ("+residue_sel+")")
            rotamerSel_LJ = fitted_rotamers.select_atoms("not type H and not (name CA or name C or name N or name O)")
        else:
            proteinNotSite = self.protein.select_atoms("protein and not ("+residue_sel+")")
            rotamerSel_LJ = fitted_rotamers.select_atoms("not (name CA or name C or name N or name O)")
            
        eps_rotamer = np.array([eps[probe_atom] for probe_atom in rotamerSel_LJ.types])
        rmin_rotamer = np.array([p_Rmin2[probe_atom] for probe_atom in rotamerSel_LJ.types])*0.5

        eps_protein = np.array([eps[probe_atom] for probe_atom in proteinNotSite.types])
        rmin_protein = np.array([p_Rmin2[probe_atom] for probe_atom in proteinNotSite.types])*0.5
        eps_ij = np.sqrt(np.multiply.outer(eps_rotamer, eps_protein))
        
        rmin_ij = np.add.outer(rmin_rotamer, rmin_protein)
        #Convert atom groups to indices for efficiecy
        proteinNotSite = proteinNotSite.indices
        #Convert indices of protein atoms (constant within each frame) to positions
        proteinNotSite = self.protein.trajectory.ts.positions[proteinNotSite]

        rotamerSel_LJ = rotamerSel_LJ.indices
        lj_energy_pose = np.zeros(len(fitted_rotamers.trajectory))
        for rotamer_counter, rotamer in enumerate(fitted_rotamers.trajectory):
            d = MDAnalysis.lib.distances.distance_array(rotamer.positions[rotamerSel_LJ],proteinNotSite)
            d = np.power(rmin_ij/d,6)
            pair_LJ_energy = eps_ij*(d*d-2.*d)
            lj_energy_pose[rotamer_counter] = pair_LJ_energy.sum()
        return np.exp(-lj_energy_pose/(gas_un*self.temp))  # for new alignment method
        # Slower implementation without for loop
        #rot_coords = fitted_rotamers.trajectory.timeseries(rotamerSel_LJ)
        #d = MDAnalysis.lib.distances.distance_array(rot_coords.reshape(-1,3),proteinNotSite).reshape(rot_coords.shape[0],rot_coords.shape[1],proteinNotSite.shape[0])
        #d = np.power(rmin_ij[:,np.newaxis,:]/d,6)
        #LJ_energy = (eps_ij[:,np.newaxis,:]*(d*d-2.*d)).sum(axis=(0,2))
        #return np.exp(-LJ_energy/(gas_un*self.temp))

    def rotamerWeights(self, rotamersSite, lib_weights_norm, residue_sel):
        # Calculate Boltzmann weights
        boltz = self.lj_calculation(rotamersSite, residue_sel)
        # save external probabilities
        # np.savetxt(self.output_prefix+'-boltz-{:s}.dat'.format(residue_sel.replace(" ", "")),boltz)
        # Set to zero Boltzmann weights that are NaN
        boltz[np.isnan(boltz)] = 0.0

        # Multiply Boltzmann weights by library weights
        boltz = lib_weights_norm * boltz
        return boltz, np.nansum(boltz)
