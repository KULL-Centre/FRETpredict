# -*- coding: utf-8 -*-
"""
FRET prediction Class
---------------------

Class to perform Fluorescence Resonance Energy Transfer calculations given two positions for the fluorofores.

"""

# Coordinates and arrays
import numpy as np
import h5py
import MDAnalysis
import math

# Logger
import logging

# Inner imports
from FRETpredict.utils import Operations

class FRETpredict(Operations):
    """Calculation of FRET signal between two chromophores."""

    def __init__(self, protein, residues, record_frames=False, **kwargs):
        """
        Args:
            protein (:py:class:`MDAnalysis.core.universe.Universe`):
            residues (list(:py:class:`str`)):
        :Keywords:
        """
        Operations.__init__(self, protein, **kwargs)
        self.residues = residues
        self.k2array = np.empty(0)
        self.rijarray = np.empty(0)
        # loads all args and kwargs
        self.record_frames = record_frames
        self.ra = -5
        self.re = 20
        self.nr = 501
        self.rax = np.linspace(self.ra, self.re, self.nr)
        self.vari = np.exp(-100*(self.rax)**2)
        self.r0 = kwargs.pop('r0', 5.4)
        if len(residues) != 2:
            raise ValueError("The residue_list must contain exactly 2 "
                             "residue numbers: current value {0}.".format(residues))
            
    def trajectoryAnalysis(self):
        logging.info("Starting rotamer distance analysis of trajectory {:s} with labeled residues "
                     "{:d} and {:d}".format(self.protein.trajectory.filename, self.residues[0], self.residues[1]))
        f = h5py.File(self.output_prefix+'-{:d}-{:d}.hdf5'.format(self.residues[0], self.residues[1]), "w")
        distributions = f.create_dataset("distributions",
                (self.protein.trajectory.n_frames, self.rax.size), fillvalue=0, compression="gzip")
        lib_1_weights_norm = self.lib_1.weights / np.sum(self.lib_1.weights)
        lib_2_weights_norm = self.lib_2.weights / np.sum(self.lib_2.weights)
        rotamer1, prot_atoms1, residue_sel1 = self.precalculate_rotamer(self.residues[0], self.chains[0], self.lib_1)
        rotamer2, prot_atoms2, residue_sel2 = self.precalculate_rotamer(self.residues[1], self.chains[1], self.lib_2)
        # For each trajectory frame, place the probes at the spin-labeled site using rotamer_placement(), calculate
        # Boltzmann weights based on Lennard-Jones interactions and calculate weighted distributions of probe-probe separations
        zarray = np.empty(0) # Array of steric partition functions (sum over Boltzmann weights)
        allk2 = np.empty(0)
        #old_k2array = np.empty(0)
        for frame_ndx, _ in enumerate(self.protein.trajectory):
            # Fit the rotamers onto the protein
            rotamersSite1 = self.rotamer_placement(rotamer1, prot_atoms1, self.lib_1)
            rotamersSite2 = self.rotamer_placement(rotamer2, prot_atoms2, self.lib_2)

            boltz1, z1 = self.rotamerWeights(rotamersSite1, lib_1_weights_norm, residue_sel1)
            boltz2, z2 = self.rotamerWeights(rotamersSite2, lib_2_weights_norm, residue_sel2)
            zarray = np.append(zarray, [z1,z2])
            if (z1 <= self.z_cutoff) or (z2 <= self.z_cutoff):
                 continue
            boltzmann_weights_norm1 = boltz1 / z1
            boltzmann_weights_norm2 = boltz2 / z2
            boltzmann_weights_norm =  boltzmann_weights_norm1.reshape(-1,1) * boltzmann_weights_norm2

            #print('boltz',boltzmann_weights_norm.shape)
            
            rotamer1_c1 = rotamersSite1.select_atoms('name '+self.lib_1.atoms[0])
            rotamer2_c1 = rotamersSite2.select_atoms('name '+self.lib_2.atoms[0])
            rotamer1_c2 = rotamersSite1.select_atoms('name '+self.lib_1.atoms[1])
            rotamer2_c2 = rotamersSite2.select_atoms('name '+self.lib_2.atoms[1])

            c1_rot1_pos = np.array([rotamer1_c1.positions for x in rotamersSite1.trajectory]).squeeze(axis=1)
            c1_rot2_pos = np.array([rotamer2_c1.positions for x in rotamersSite2.trajectory]).squeeze(axis=1)
            c2_rot1_pos = np.array([rotamer1_c2.positions for x in rotamersSite1.trajectory]).squeeze(axis=1)
            c2_rot2_pos = np.array([rotamer2_c2.positions for x in rotamersSite2.trajectory]).squeeze(axis=1)

            """
            c1_rot1_pos = np.array([i for x in rotamersSite1.trajectory for i in rotamer1_c1.positions])
            c1_rot2_pos = np.array([i for x in rotamersSite2.trajectory for i in rotamer2_c1.positions])
            c2_rot1_pos = np.array([i for x in rotamersSite1.trajectory for i in rotamer1_c2.positions])
            c2_rot2_pos = np.array([i for x in rotamersSite2.trajectory for i in rotamer2_c2.positions])
            """

            #print('pos',c1_rot1_pos.shape,c1_rot2_pos.shape)

            mu_1 = c2_rot1_pos - c1_rot1_pos
            mu_1 /= np.linalg.norm(mu_1, axis=1, keepdims=True)

            mu_2 = c1_rot2_pos - c2_rot2_pos
            mu_2 /= np.linalg.norm(mu_2, axis=1, keepdims=True)

            #print('mu',mu_1.shape,mu_2.shape)

            mu_cosine = np.einsum('ik,jk->ij', mu_1, mu_2)

            #print('cosine',mu_cosine.shape,mu_cosine.min(),mu_cosine.max(),mu_cosine.mean())

            #mid_1 = (c2_rot1_pos - c1_rot1_pos)/2

            #mid_2 = (c2_rot2_pos - c1_rot2_pos)/2

            r_12 = c2_rot1_pos[:,np.newaxis,:] - c1_rot2_pos
            r_12 /= np.linalg.norm(r_12, axis=2, keepdims=True)

            r_21 = c1_rot1_pos[:,np.newaxis,:] - c2_rot2_pos
            r_21 /= np.linalg.norm(r_21, axis=2, keepdims=True)

            #print('r12',r_12.shape,'r21',r_21.shape)

            trans_dip_moment12_rot1 = np.einsum('ijk,ik->ij', r_12, mu_1)
            trans_dip_moment21_rot1 = np.einsum('ijk,ik->ij', r_21, mu_1)

            trans_dip_moment_rot1 = (trans_dip_moment12_rot1 + trans_dip_moment21_rot1) / 2.

            trans_dip_moment12_rot2 = np.einsum('ijk,jk->ij', r_12, mu_2)
            trans_dip_moment21_rot2 = np.einsum('ijk,jk->ij', r_21, mu_2)

            trans_dip_moment_rot2 = (trans_dip_moment12_rot2 + trans_dip_moment21_rot2) / 2.

            #print('trans',trans_dip_moment_rot1.shape,trans_dip_moment_rot2.shape)
            
            #print('1',np.mean(trans_dip_moment_rot1))
            #print('2',np.mean(trans_dip_moment_rot2))
            k2 = np.power(mu_cosine - 3*trans_dip_moment_rot1*trans_dip_moment_rot2, 2)
            allk2 = np.append(allk2,k2.flatten())
            self.k2array = np.append(self.k2array,np.sum(k2*boltzmann_weights_norm))

            """
            kappa2_list = []
            weights_list = []

            for index_1, boltz_1 in enumerate(boltzmann_weights_norm1):
                mutual_displacement = np.dot(mu_2, mu_1[index_1].T)

                outer_vect1 = c2_rot1_pos[index_1] - c1_rot2_pos
                outer_vect1 /= np.linalg.norm(outer_vect1, axis=1)[:, np.newaxis]
                outer_vect2 = c1_rot1_pos[index_1] - c2_rot2_pos
                outer_vect2 /= np.linalg.norm(outer_vect2, axis=1)[:, np.newaxis]

                trans_dip_moment1_probe1 = np.dot(outer_vect1, mu_1[index_1].T)
                trans_dip_moment2_probe1 = np.dot(outer_vect2, mu_1[index_1].T)

                for index_2, boltz_2 in enumerate(boltzmann_weights_norm2):  # could we take out this for loop?
                    trans_dip_moment1_probe2 = np.dot(outer_vect1[index_2], mu_2[index_2])
                    trans_dip_moment2_probe2 = np.dot(outer_vect2[index_2], mu_2[index_2])

                    new_k2 = np.power((mutual_displacement[index_2] -
                                       ((3. / 4.) * (trans_dip_moment1_probe1[index_2] +
                                                     trans_dip_moment2_probe1[index_2]) * (
                                            trans_dip_moment1_probe2 +
                                            trans_dip_moment2_probe2))), 2)
                    kappa2_list.append(new_k2)
                    weights_list.append(boltz_1*boltz_2)
            old_k2array = np.append(old_k2array, np.average(kappa2_list, weights=weights_list))
            """

            rotamer1oxigen = rotamersSite1.select_atoms("name O6")
            rotamer2oxigen = rotamersSite2.select_atoms("name O6")

            oxi1_pos = np.array([rotamer1oxigen.positions for x in rotamersSite1.trajectory]).squeeze(axis=1)
            oxi2_pos = np.array([rotamer2oxigen.positions for x in rotamersSite2.trajectory]).squeeze(axis=1)
            # Distances between nitroxide groups
            dists_array = np.linalg.norm(oxi1_pos[:,np.newaxis,:] - oxi2_pos, axis=2) / 10
            self.rijarray = np.append(self.rijarray,np.sum(dists_array*boltzmann_weights_norm))
            dists_array = np.round((self.nr * (dists_array - self.ra)) / (self.re - self.ra)).astype(int).flatten()
            distribution = np.bincount(dists_array, weights=boltzmann_weights_norm.flatten(), minlength=self.rax.size) 
            distributions[frame_ndx] = distribution
        f.close()
        np.savetxt(self.output_prefix+'-Z-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]),zarray.reshape(-1,2))
        np.savetxt(self.output_prefix+'-k2-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]),self.k2array)
        np.savetxt(self.output_prefix+'-allk2-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]),allk2)
        #np.savetxt(self.output_prefix+'-old_k2-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]),old_k2array)

    def save(self,filename):
        f = h5py.File(filename, "r")
        distributions = f.get('distributions')
        if isinstance(self.weights, np.ndarray):
            if self.weights.size != distributions.shape[0]:
                    logging.info('Weights array has size {} whereas the number of frames is {}'.
                            format(self.weights.size, distributions.shape[0]))
                    raise ValueError('Weights array has size {} whereas the number of frames is {}'.
                            format(self.weights.size, distributions.shape[0]))
        elif self.weights == False:
            self.weights = np.ones(distributions.shape[0])
        else:
            logging.info('Weights argument should be a numpy array')
            raise ValueError('Weights argument should be a numpy array')
        distribution = np.nansum(distributions*self.weights.reshape(-1,1), 0)
        frame_inv_distr = np.fft.ifft(distribution) * np.fft.ifft(self.vari)
        smoothed = np.real(np.fft.fft(frame_inv_distr))
        smoothed /= np.trapz(smoothed, self.rax)
        np.savetxt(self.output_prefix + '-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]),
                np.c_[self.rax[100:401], smoothed[200:]],
                   header='distance distribution')
        np.savetxt(self.output_prefix + '-dist-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]),
                np.c_[self.rax, distribution],
                   header='distance distribution')
        f.close()
        e_static = 1. / ( 1 + 2/3./self.k2array * ( np.power(self.rijarray / self.r0, 6) ) )
        e_dynamic = 1. / ( 1 + 2/3./np.sum(self.k2array*self.weights) * (np.power(self.rijarray / self.r0, 6)) )
        np.savetxt(self.output_prefix + '-efficiency-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]),
            np.c_[e_dynamic, e_static, self.weights],header='E_dynamic E_static weights')




    def run(self, **kwargs):
        # Output
        self.output_prefix = kwargs.get('output_prefix', 'res')
        # Input
        self.load_file = kwargs.get('load_file', False)
        # Weights for each frame
        self.weights = kwargs.get('weights', False)
        if self.load_file:
            if os.path.isfile(self.load_file):
                logging.info('Loading pre-computed data from {} - will not load trajectory file.'.format(self.load_file))
            else:
                logging.info('File {} not found!'.format(self.load_file))
                raise FileNotFoundError('File {} not found!'.format(self.load_file))
            self.save(self.load_file)    
        else:
            self.trajectoryAnalysis()
            self.save(self.output_prefix+'-{:d}-{:d}.hdf5'.format(self.residues[0], self.residues[1]))
