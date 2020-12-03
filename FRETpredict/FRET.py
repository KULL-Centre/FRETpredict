# -*- coding: utf-8 -*-
"""
FRET prediction Class
---------------------

Class to perform Fluorescence Resonance Energy Transfer calculations given two positions for the fluorofores.

"""

# Coordinates and arrays
import numpy as np
import pandas as pd
import h5py
import MDAnalysis
import math

# Logger
import logging

# Inner imports
from FRETpredict.utils import Operations

class FRETpredict(Operations):
    """Calculation of FRET signal between two chromophores."""

    def __init__(self, protein, residues, **kwargs):
        """
        Args:
            protein (:py:class:`MDAnalysis.core.universe.Universe`):
            residues (list(:py:class:`str`)):
        :Keywords:
        """
        Operations.__init__(self, protein, **kwargs)
        self.residues = residues
        self.r0 = kwargs.pop('r0', 5.4)
        logging.basicConfig(filename=kwargs.get('log_file', 'log'),level=logging.INFO)
        for i in range(2):
            residue_sel = "resid {:d}".format(self.residues[i])
            if type(self.chains[i]) == str:
                residue_sel += " and segid {:s}".format(self.chains[i])
            logging.info('{:s} = {:s}'.format(residue_sel,self.protein.select_atoms(residue_sel).atoms.resnames[0]))

        dr = 0.05
        self.rmin = -5
        self.rmax = kwargs.get('rmax', 12)*2 - self.rmin
        self.nr = int(round((self.rmax-self.rmin)/dr,0) + 1)
        self.rax = np.linspace(self.rmin, self.rmax, self.nr)

        if len(residues) != 2:
            raise ValueError("The residue_list must contain exactly 2 "
                             "residue numbers: current value {0}.".format(residues))
            
    def trajectoryAnalysis(self):
        logging.info("Starting rotamer distance analysis of trajectory {:s} with labeled residues "
                     "{:d} and {:d}".format(self.protein.trajectory.filename, self.residues[0], self.residues[1]))
        f = h5py.File(self.output_prefix+'-{:d}-{:d}.hdf5'.format(self.residues[0], self.residues[1]), "w")
        distributions = f.create_dataset("distributions",
                (self.protein.trajectory.n_frames, self.rax.size), fillvalue=0, compression="gzip")
        rotamer1, prot_atoms1, LJ_data1 = self.precalculate_rotamer(self.residues[0], self.chains[0], self.lib_1)
        rotamer2, prot_atoms2, LJ_data2 = self.precalculate_rotamer(self.residues[1], self.chains[1], self.lib_2)
        # For each trajectory frame, place the probes at the spin-labeled site using rotamer_placement(), calculate
        # Boltzmann weights based on Lennard-Jones interactions and calculate weighted distributions of probe-probe separations
        zarray = np.empty(0) # Array of steric partition functions (sum over Boltzmann weights)
        k2_avg = np.full(self.protein.trajectory.n_frames, np.nan)
        esta_avg = np.full(self.protein.trajectory.n_frames, np.nan)
        edyn_avg = np.full(self.protein.trajectory.n_frames, np.nan)

        for frame_ndx, _ in enumerate(self.protein.trajectory):
            # Fit the rotamers onto the protein
            rotamersSite1 = self.rotamer_placement(rotamer1, prot_atoms1, self.lib_1)
            rotamersSite2 = self.rotamer_placement(rotamer2, prot_atoms2, self.lib_2)

            boltz1, z1 = self.rotamerWeights(rotamersSite1, self.lib_1, LJ_data1)
            boltz2, z2 = self.rotamerWeights(rotamersSite2, self.lib_2, LJ_data2)
            boltzmann_weights_norm = (boltz1.reshape(-1,1) * boltz2).flatten()
            zarray = np.append(zarray, [z1,z2])
            if (z1 <= self.z_cutoff) or (z2 <= self.z_cutoff):
                 continue
            
            rotamer1_c1 = rotamersSite1.select_atoms('name '+self.lib_1.mu[0])
            rotamer2_c1 = rotamersSite2.select_atoms('name '+self.lib_2.mu[0])
            rotamer1_c2 = rotamersSite1.select_atoms('name '+self.lib_1.mu[1])
            rotamer2_c2 = rotamersSite2.select_atoms('name '+self.lib_2.mu[1])

            c1_rot1_pos = np.array([rotamer1_c1.positions for x in rotamersSite1.trajectory]).squeeze(axis=1)
            c1_rot2_pos = np.array([rotamer2_c1.positions for x in rotamersSite2.trajectory]).squeeze(axis=1)
            c2_rot1_pos = np.array([rotamer1_c2.positions for x in rotamersSite1.trajectory]).squeeze(axis=1)
            c2_rot2_pos = np.array([rotamer2_c2.positions for x in rotamersSite2.trajectory]).squeeze(axis=1)

            mu_1 = c2_rot1_pos - c1_rot1_pos
            mu_1 /= np.linalg.norm(mu_1, axis=1, keepdims=True)

            mu_2 = c1_rot2_pos - c2_rot2_pos
            mu_2 /= np.linalg.norm(mu_2, axis=1, keepdims=True)
            
            mu_cosine = np.einsum('ik,jk->ij', mu_1, mu_2)

            rotamer1centre = rotamersSite1.select_atoms('name '+self.lib_1.r[0])
            rotamer2centre = rotamersSite2.select_atoms('name '+self.lib_2.r[0])

            centre1_pos = np.array([rotamer1centre.positions for x in rotamersSite1.trajectory]).squeeze(axis=1)
            centre2_pos = np.array([rotamer2centre.positions for x in rotamersSite2.trajectory]).squeeze(axis=1)
            
            rvec = centre1_pos[:,np.newaxis,:] - centre2_pos

            rdist = np.linalg.norm(rvec, axis=2, keepdims=True)

            runitvec = rvec / rdist

            rdist = rdist.flatten() / 10

            r_mu_1 = np.einsum('ijk,ik->ij', runitvec, mu_1)
            r_mu_2 = np.einsum('ijk,jk->ij', runitvec, mu_2)
 
            k2 = np.power(mu_cosine - 3*r_mu_1*r_mu_2, 2).flatten()

            r_12 = c2_rot1_pos[:,np.newaxis,:] - c1_rot2_pos
            r_12 /= np.linalg.norm(r_12, axis=2, keepdims=True)

            r_21 = c1_rot1_pos[:,np.newaxis,:] - c2_rot2_pos
            r_21 /= np.linalg.norm(r_21, axis=2, keepdims=True)

            trans_dip_moment12_rot1 = np.einsum('ijk,ik->ij', r_12, mu_1)
            trans_dip_moment21_rot1 = np.einsum('ijk,ik->ij', r_21, mu_1)

            trans_dip_moment_rot1 = (trans_dip_moment12_rot1 + trans_dip_moment21_rot1) / 2.

            trans_dip_moment12_rot2 = np.einsum('ijk,jk->ij', r_12, mu_2)
            trans_dip_moment21_rot2 = np.einsum('ijk,jk->ij', r_21, mu_2)

            trans_dip_moment_rot2 = (trans_dip_moment12_rot2 + trans_dip_moment21_rot2) / 2.

            k2_avg[frame_ndx] = np.dot(k2, boltzmann_weights_norm)

            ratio6 = np.power(rdist / self.r0, 6)

            estatic = 1. / ( 1 + 2/3. * np.divide(ratio6,k2))
            edynamic = 1. / ( 1 + 2/3./k2_avg[frame_ndx] * ratio6 )
 
            esta_avg[frame_ndx] = np.dot( estatic, boltzmann_weights_norm )
            edyn_avg[frame_ndx] = np.dot( edynamic, boltzmann_weights_norm )

            rdist = np.round((self.nr * (rdist - self.rmin)) / (self.rmax - self.rmin)).astype(int).flatten()
            distribution = np.bincount(rdist, weights=boltzmann_weights_norm.flatten(), minlength=self.rax.size)
            distributions[frame_ndx] = distribution
        f.close()
        np.savetxt(self.output_prefix+'-Z-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]),zarray.reshape(-1,2))
        np.savetxt(self.output_prefix+'-k2-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]),k2_avg)
        np.savetxt(self.output_prefix+'-Es-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]),esta_avg)
        np.savetxt(self.output_prefix+'-Ed-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]),edyn_avg)

    def save(self):
        k2 = np.loadtxt(self.output_prefix+'-k2-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]))
        if isinstance(self.weights, np.ndarray):
            if self.weights.size != k2.size:
                    logging.info('Weights array has size {} whereas the number of frames is {}'.
                            format(self.weights.size, k2.size))
                    raise ValueError('Weights array has size {} whereas the number of frames is {}'.
                            format(self.weights.size, k2.size))
        elif self.weights == False:
            self.weights = np.ones(k2.size)
        else:
            logging.info('Weights argument should be a numpy array')
            raise ValueError('Weights argument should be a numpy array')
        f = h5py.File(self.output_prefix+'-{:d}-{:d}.hdf5'.format(self.residues[0], self.residues[1]), "r")
        distributions = f.get('distributions')
        distribution = np.nansum(distributions*self.weights.reshape(-1,1), 0)
        frame_inv_distr = np.fft.ifft(distribution) * np.fft.ifft(np.exp(-0.5*(self.rax/self.stdev)**2))
        smoothed = np.real(np.fft.fft(frame_inv_distr))
        smoothed /= np.trapz(smoothed, self.rax)
        np.savetxt(self.output_prefix + '-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]),
                np.c_[self.rax[100:-100], smoothed[200:]], header='distance distribution')
        f.close()
        estatic = np.loadtxt(self.output_prefix+'-Es-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]))
        edynamic = np.loadtxt(self.output_prefix+'-Ed-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]))
        if k2.size == 1:
            df = pd.Series([k2, estatic, edynamic], index=['k2','Estatic','Edynamic'])
        else:
            df = pd.DataFrame([self.weightedAvgSDSE(k2[np.isfinite(k2)],self.weights[np.isfinite(k2)]), 
                self.weightedAvgSDSE(estatic[np.isfinite(k2)],self.weights[np.isfinite(k2)]), 
                self.weightedAvgSDSE(edynamic[np.isfinite(k2)],self.weights[np.isfinite(k2)])], 
                columns = ['Average', 'SD','SE'], index=['k2','Estatic','Edynamic'])
        df.to_pickle(self.output_prefix + '-data-{:d}-{:d}.pkl'.format(self.residues[0], self.residues[1]))

    def run(self, **kwargs):
        # Output
        self.output_prefix = kwargs.get('output_prefix', 'res')
        # Input
        self.load_file = kwargs.get('load_file', False)
        # Weights for each frame
        self.weights = kwargs.get('weights', False)
        self.stdev = kwargs.get('filter_stdev', 0.02)
        if self.load_file:
            if os.path.isfile(self.load_file):
                logging.info('Loading pre-computed data from {} - will not load trajectory file.'.format(self.load_file))
            else:
                logging.info('File {} not found!'.format(self.load_file))
                raise FileNotFoundError('File {} not found!'.format(self.load_file))
            self.save(self.load_file)    
        else:
            self.trajectoryAnalysis()
            self.save()
