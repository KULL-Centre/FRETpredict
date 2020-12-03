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
from FRETpredict.lennardjones import lj_parameters
import FRETpredict.libraries as libraries
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
        self.libname_1 = kwargs.get('libname_1', 'Alexa 488')
        self.lib_1 = libraries.RotamerLibrary(self.libname_1)
        self.libname_2 = kwargs.get('libname_2', 'Alexa 594')
        self.lib_2 = libraries.RotamerLibrary(self.libname_2)
        self.temp = kwargs.get('temperature', 300)
        self.z_cutoff = kwargs.get('z_cutoff', 0.05)
        self.ign_H = kwargs.get('ign_H', True)
        self.electrostatic = kwargs.get('electrostatic', False)
        self.chains = kwargs.get('chains', [None,None])
        # scaling factor to modulate probe-protein interactions
        sigma_scaling = kwargs.get('sigma_scaling', 1)
        epsilon_scaling = kwargs.get('epsilon_scaling', 1)
        self.rmin2LJ = {atom:lj_parameters[atom]['p_Rmin2']*sigma_scaling for atom in lj_parameters}
        self.epsLJ = {atom:lj_parameters[atom]['eps']*epsilon_scaling for atom in lj_parameters}

    def precalculate_rotamer(self, residue, chain, lib):
        residue_sel = "resid {:d}".format(residue)
        if type(chain) == str:
            residue_sel += " and segid {:s}".format(chain)
        prot_Ca = self.protein.select_atoms('protein and name CA and '+residue_sel)
        prot_Co = self.protein.select_atoms('protein and name C and '+residue_sel)
        prot_N = self.protein.select_atoms('protein and name N and '+residue_sel)
        probe_coords = np.zeros((len(lib.top.atoms),1, 3))
        universe = MDAnalysis.Universe(lib.top.filename, probe_coords, format=MemoryReader, order='afc')

        # Atom indices and parameters for LJ calculation
        if self.ign_H:
            proteinNotSite = self.protein.select_atoms("protein and not type H and not ("+residue_sel+")")
            rotamerSel_LJ = universe.select_atoms("not type H and not (name CA or name C or name N or name O)")
        else:
            proteinNotSite = self.protein.select_atoms("protein and not ("+residue_sel+")")
            rotamerSel_LJ = universe.select_atoms("not (name CA or name C or name N or name O)")
           
        eps_rotamer = np.array([self.epsLJ[probe_atom] for probe_atom in rotamerSel_LJ.types])
        rmin_rotamer = np.array([self.rmin2LJ[probe_atom] for probe_atom in rotamerSel_LJ.types])
        eps_protein = np.array([self.epsLJ[prot_atom] for prot_atom in proteinNotSite.types])
        rmin_protein = np.array([self.rmin2LJ[prot_atom] for prot_atom in proteinNotSite.types])
        eps_ij = np.sqrt(np.multiply.outer(eps_rotamer, eps_protein))
        rmin_ij = np.add.outer(rmin_rotamer, rmin_protein)

        LJ_data = [proteinNotSite.indices,rotamerSel_LJ.indices,eps_ij,rmin_ij]

        # Charges for Yukawa potential
        if self.electrostatic:
            rot_positive = rotamerSel_LJ.select_atoms('name '+' or name '.join(lib.positive))
            rot_negative = rotamerSel_LJ.select_atoms('name '+' or name '.join(lib.negative))
            charge = lambda x : -1 if x in rot_negative else 0.5 if x in rot_positive else 0
            q_rotamer = np.array([charge(probe_atom) for probe_atom in rotamerSel_LJ])
            prot_positive = proteinNotSite.select_atoms("(name CZ and resname ARG) or (name NZ and resname LYS)").indices
            prot_negative = proteinNotSite.select_atoms("(name CG and resname ASP) or (name CD and resname GLU)").indices
            prot_his = proteinNotSite.select_atoms("(name ND1 and resname HIS) or (name NE2 and resname HIS)").indices
            charge = lambda x : -1 if x in prot_negative else 1 if x in prot_positive else 0.25 if x in prot_his else 0
            q_prot = np.array([charge(prot_atom) for prot_atom in proteinNotSite.indices])
            q_ij = np.multiply.outer(q_rotamer, q_prot)
            LJ_data.append(q_ij)

        return universe, (prot_Ca, prot_Co, prot_N), LJ_data
        
    def rotamer_placement(self, universe, prot_atoms, lib):
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
        print(lib.coord.shape,rotation.shape)
        probe_coords = np.tensordot(lib.coord,rotation,axes=([2],[0])) + offset
        universe.load_new(probe_coords, format=MemoryReader, order='afc')
        #save aligned rotamers
        #mtssl = universe.select_atoms("all")
        #with MDAnalysis.Writer(lib.name+".pdb", mtssl.n_atoms) as W:
        #    for ts in universe.trajectory:
        #        W.write(mtssl)
        return universe

    def lj_calculation(self, fitted_rotamers, LJ_data):
        gas_un = 1.9858775e-3 # CHARMM, in kcal/mol*K
        if self.electrostatic:
            proteinNotSite, rotamerSel_LJ, eps_ij, rmin_ij, q_ij = LJ_data
        else:
            proteinNotSite, rotamerSel_LJ, eps_ij, rmin_ij = LJ_data
        #Convert indices of protein atoms (constant within each frame) to positions
        proteinNotSite = self.protein.trajectory.ts.positions[proteinNotSite]
        lj_energy_pose = np.zeros(len(fitted_rotamers.trajectory))
        dh_energy_pose = np.zeros(len(fitted_rotamers.trajectory))
        for rotamer_counter, rotamer in enumerate(fitted_rotamers.trajectory):
            d = MDAnalysis.lib.distances.distance_array(rotamer.positions[rotamerSel_LJ],proteinNotSite)
            if self.electrostatic:
                cutoff = (q_ij!=0) & (d<20)
                pair_dh_energy = q_ij[cutoff]*7.0/d[cutoff]*np.exp(-d[cutoff]/10.0)
                dh_energy_pose[rotamer_counter] = pair_dh_energy.sum()
            cutoff = d<10
            d = np.power(rmin_ij[cutoff]/d[cutoff],6)
            pair_lj_energy = eps_ij[cutoff]*(d*d-2.*d)
            lj_energy_pose[rotamer_counter] = pair_lj_energy.sum()
        return np.exp(-lj_energy_pose/(gas_un*self.temp)-dh_energy_pose)  # for new alignment method

    def rotamerWeights(self, rotamersSite, lib, LJ_data):
        # Calculate Boltzmann weights
        boltz = self.lj_calculation(rotamersSite, LJ_data)
        # save external probabilities
        #np.savetxt(self.output_prefix+'-boltz-{:s}.dat'.format(residue_sel.replace(" ", "")),boltz)
        # Set to zero Boltzmann weights that are NaN
        boltz[np.isnan(boltz)] = 0.0
        # Multiply Boltzmann weights by library weights
        boltz = lib.weights * boltz
        steric_z = np.sum(boltz)
        return boltz/steric_z, steric_z

    def weightedAvgSDSE(self, values, weights):
        # Calculate the weighted average and standard deviation.
        avg = np.average(values, weights=weights)
        variance = np.average((values-avg)**2, weights=weights)
        return (avg, np.sqrt(variance), np.sqrt(variance/values.size))
