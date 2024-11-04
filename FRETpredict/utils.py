import numpy as np
import MDAnalysis
from MDAnalysis.coordinates.memory import MemoryReader
from .lennardjones import lj_parameters
from .libraries import *


class Operations(object):
    """

    Utility functions for the Calculation of the distance profile between a probe and backbone amide.

        Attributes
        ==========

            protein : MDAnalysis.core.universe.Universe
                Protein MDAnalysis Universe object

            libname_1 : str
                Name of the donor rotamer library

            lib_1 : dict
                Rotamer library of the first chromophore

            libname_2 : str
                Name of the acceptor rotamer library

            lib_2 : dict
                Rotamer library of the second chromophore

            temp : int
                Temperature

            ign_H : Bool
                Consider or ignore H atoms in LJ calculations

            electrostatic : Bool
                Perform or ignore electrostatic calculations

            potential : str
                Select potential for VdW calculations ('lj' or 'gauss')

            chains : str
                Name of the protein chains in which the placement residues are located

            rmin2LJ : float
                Scaled Sigma of the LJ distribution

            epsLJ : float
                Scaled epsilon LJ parameter

        Methods
        =======

            precalculate_rotamer(residue, chain, lib)
                Select placement residue atoms, compute Lennard-Jones and Electrostatic (Debye-Huckel) parameters for
                Chromophore rotamers and Protein.

            rotamer_placement(universe, prot_atoms, lib)
                Place chromophore rotamers on a selected protein residue by rotation and translation, return chromophore
                Universe object with new coordinates.

            lj_calculation(fitted_rotamers, LJ_data)
                Calculate Boltzmann weights for each rotamer pose by summing Electrostatic and potential energy
                contributions.

            rotamerWeights(rotamersSite, lib, LJ_data)
                Calculate Boltzmann distribution for each rotamer library.

            weightedAvgSDSE(values, weights)
                Calculate the weighted average and standard deviation.

            calculate_ws()
                Calculate per-frame weights for reweighting.

            fraction_frames():
                Compute effective fraction of frames contributing to the averages


    """

    def __init__(self, protein, **kwargs):

        """

        Parameters
        ==========

            protein: MDAnalysis.Universe
                MDA Universe of the protein (single structure/trajectory)

        """

        # Protein MDA Universe assignment
        self.protein = protein

        # Chromophores libraries assignments
        # Chromophore 1
        self.libname_1 = kwargs.get('libname_1', 'Alexa 488')
        self.lib_1 = RotamerLibrary(self.libname_1)

        # Chromophore 2
        self.libname_2 = kwargs.get('libname_2', 'Alexa 594')
        self.lib_2 = RotamerLibrary(self.libname_2)

        # Parameters value for Boltzmann factor calculation assignment
        self.temp = kwargs.get('temperature', 300)

        # Calculate only for heavy atoms
        self.ign_H = kwargs.get('ign_H', True)

        # Ignore electrostatic interactions
        self.electrostatic = kwargs.get('electrostatic', False)

        # Potential energy function
        self.potential = kwargs.get('potential', 'lj')

        # Define protein chains
        self.chains = kwargs.get('chains', [None, None])

        # Scaling factor to modulate probe-protein interactions
        # Sigma
        sigma_scaling = kwargs.get('sigma_scaling', 0.5)
        sigma_scaling = sigma_scaling if (self.potential == 'lj') else sigma_scaling / np.power(2, 1. / 6.)
        self.rmin2LJ = {atom: lj_parameters[atom]['p_Rmin2'] * sigma_scaling for atom in lj_parameters}

        # Epsilon
        epsilon_scaling = kwargs.get('epsilon_scaling', 1)
        self.epsLJ = {atom: lj_parameters[atom]['eps'] * epsilon_scaling for atom in lj_parameters}

    def precalculate_rotamer(self, residue, chain, lib):

        """

        Select placement residue atoms, compute Lennard-Jones and Electrostatic (Debye-Huckel) parameters for
        Chromophore rotamers and Protein

        Parameters
        ==========

             residue: int
                Placement residue number

             chain: char
                Chain ID of the placement residue

             lib: dict
                Rotamers library of the chromophore

        Returns
        =======

            universe: MDAnalysis.Universe
                Chromophore Universe

            (prot_Ca, prot_Co, prot_N): MDAnalysis.AtomGroup
                Selection of the placement residue backbone amide atoms

            LJ_data: list
                list containing data for LJ calculations (proteinNotSite.indices, rotamerSel_LJ.indices, eps_ij, rmin_ij, q_ij)

        """

        # Select residue for chromophore placement
        residue_sel = "resid {:d}".format(residue)
        if type(chain) == str:
            residue_sel += " and segid {:s}".format(chain)

        # Select placement residue backbone amide atoms (to be used in rotamer_placement)
        prot_Ca = self.protein.select_atoms('protein and name CA and ' + residue_sel)
        prot_Co = self.protein.select_atoms('protein and name C and ' + residue_sel)
        prot_N = self.protein.select_atoms('protein and name N and ' + residue_sel)

        # Create Chromophore MDAnalysis Universe (to be used in rotamer_placement)
        probe_coords = np.zeros((len(lib.top.atoms), 1, 3))
        universe = MDAnalysis.Universe(lib.top.filename, probe_coords, format=MemoryReader, order='afc')

        # Protein and chromophore atom indices for LJ calculation
        # Without H atoms (heavy atoms only)
        if self.ign_H:
            proteinNotSite = self.protein.select_atoms("protein and not type H and not (" + residue_sel + ")")
            rotamerSel_LJ = universe.select_atoms("not type H and not (name CA or name C or name N or name O)")

        # With H atoms
        else:
            proteinNotSite = self.protein.select_atoms("protein and not (" + residue_sel + ")")
            rotamerSel_LJ = universe.select_atoms("not (name CA or name C or name N or name O)")

        # Parameters for LJ Calculations
        # Chromophore
        eps_rotamer = np.array([self.epsLJ[probe_atom] for probe_atom in rotamerSel_LJ.types])
        rmin_rotamer = np.array([self.rmin2LJ[probe_atom] for probe_atom in rotamerSel_LJ.types])

        # Protein
        eps_protein = np.array([self.epsLJ[prot_atom] for prot_atom in proteinNotSite.types])
        rmin_protein = np.array([self.rmin2LJ[prot_atom] for prot_atom in proteinNotSite.types])

        # Combined parameters
        eps_ij = np.sqrt(np.multiply.outer(eps_rotamer, eps_protein))
        rmin_ij = np.add.outer(rmin_rotamer, rmin_protein)
        LJ_data = [proteinNotSite.indices, rotamerSel_LJ.indices, eps_ij, rmin_ij]

        # Charges for Debye-Huckel calculations (if electrostatic calculations are enabled)
        if self.electrostatic:

            # Chromophore charged atoms selections (if present)
            try:
                rot_positive = rotamerSel_LJ.select_atoms('name ' + ' or name '.join(lib.positive))
            except:
                rot_positive = []

            try:
                rot_negative = rotamerSel_LJ.select_atoms('name ' + ' or name '.join(lib.negative))
            except:
                rot_negative = []

            # Charge assignment
            charge = lambda x: -1 if x in rot_negative else 0.5 if x in rot_positive else 0

            # Chromophore charges array
            q_rotamer = np.array([charge(probe_atom) for probe_atom in rotamerSel_LJ])

            # Charged (+ Histidine) protein atoms selections
            prot_positive = proteinNotSite.select_atoms(
                "(name CZ and resname ARG) or (name NZ and resname LYS)").indices
            prot_negative = proteinNotSite.select_atoms(
                "(name CG and resname ASP) or (name CD and resname GLU)").indices
            prot_his = proteinNotSite.select_atoms("(name ND1 and resname HIS) or (name NE2 and resname HIS)").indices

            # Charge assignment
            charge = lambda x: -1 if x in prot_negative else 1 if x in prot_positive else 0.25 if x in prot_his else 0

            # Protein charges array
            q_prot = np.array([charge(prot_atom) for prot_atom in proteinNotSite.indices])

            # Combined parameters
            q_ij = np.multiply.outer(q_rotamer, q_prot)

            # Concatenate Electrostatic data to LJ data
            LJ_data.append(q_ij)

        return universe, (prot_Ca, prot_Co, prot_N), LJ_data

    def rotamer_placement(self, universe, prot_atoms, lib):

        """

        Place chromophore rotamers on a selected protein residue by rotation and translation, return chromophore
        Universe object with new coordinates, and save PDB with aligned rotamers

        Parameters
        ==========

            universe: MDAnalysis.Universe
                Chromophore MDAnalysis Universe

            prot_atoms: MDAnalysis.AtomGroup
                Selection of the placement residue backbone amide atoms

            lib: dict
                Rotamers library of the chromophore

        Returns
        =======

            universe: MDAnalysis.Universe
                Chromophore rotamer library Universe object with new coordinates

        """

        # Atoms of the chromophore placement residue (from precalculate_rotamer)
        prot_Ca, prot_Co, prot_N = prot_atoms

        # Use the chromophore placement residue Ca atoms positions as offset
        offset = prot_Ca.positions.copy()

        # Set new reference system with origin in placement residue Ca
        # (Translate chromophore placement residue atoms position using the offset)
        Ca_coords = prot_Ca.positions - offset
        Co_coords = prot_Co.positions - offset
        N_coords = prot_N.positions - offset

        # Create unitary x vector (Ca-N bond)
        x_vector = N_coords - Ca_coords
        x_vector /= np.linalg.norm(x_vector)

        # Create unitary y_t vector (Co-Ca bond, to obtain z vector)
        yt_vector = Co_coords - Ca_coords
        yt_vector /= np.linalg.norm(yt_vector)

        # Create unitary z vector (perpendicular to plane formed by Ca-N and Co-Ca bonds)
        z_vector = np.cross(x_vector, yt_vector)
        z_vector /= np.linalg.norm(z_vector)

        # Create unitary y vector (perpendicular to plane formed by z vector and Ca-N bond)
        y_vector = np.cross(z_vector, x_vector)

        # Stack x, y, z vectors vertically to obtain rotation matrix
        rotation = np.vstack([x_vector, y_vector, z_vector])

        # Rotate chromophore rotamers using the rotation matrix, and add the offset to place it on the placement residue
        probe_coords = np.tensordot(lib.coord, rotation, axes=([2], [0])) + offset

        # Load new atoms position into the chromophore MDA Universe
        universe.load_new(probe_coords, format=MemoryReader, order='afc')

        # Save aligned rotamers
        if self.verbose:
            rotamers = universe.select_atoms("all")
            with MDAnalysis.Writer(self.output_prefix+'_'+lib.name.replace(' ','_')+'.pdb', rotamers.n_atoms) as W:
                for ts in universe.trajectory:
                    W.write(rotamers)

        return universe

    def lj_calculation(self, fitted_rotamers, LJ_data):

        """

        Calculate Boltzmann weights for each rotamer pose by summing Electrostatic and potential energy contributions

        Parameters
        ==========

            fitted_rotamers: MDAnalysis.Universe
                Universe object of the Chromophore rotamer library placed on the protein residue

            LJ_data: list
                List containing data for LJ calculations (proteinNotSite.indices, rotamerSel_LJ.indices, eps_ij, rmin_ij, q_ij)

        Return
        ======

            boltzmann_factor: numpy.array
                Array of Boltzmann weights (one for each trajectory frame)

        """
        # Universal gas constant (CHARMM, in kcal/mol*K)
        gas_un = 1.9858775e-3

        # Assign LJ data to different variables (same naming as in precalculate_rotamer)
        if self.electrostatic:
            proteinNotSite, rotamerSel_LJ, eps_ij, rmin_ij, q_ij = LJ_data
        else:
            proteinNotSite, rotamerSel_LJ, eps_ij, rmin_ij = LJ_data

        # Convert indices of protein atoms (constant within each frame) to positions
        proteinNotSite = self.protein.trajectory.ts.positions[proteinNotSite]

        # Allocate arrays for potential_energy and dh_energy
        pot_energy_pose = np.zeros(len(fitted_rotamers.trajectory))
        dh_energy_pose = np.zeros(len(fitted_rotamers.trajectory))

        # Calculate energy (potential+electrostatic) of every rotamer in library
        for rotamer_counter, rotamer in enumerate(fitted_rotamers.trajectory):

            # Calculate array of all distances of rotamer atoms from protein atoms
            d = MDAnalysis.lib.distances.distance_array(rotamer.positions[rotamerSel_LJ], proteinNotSite)

            # Electrostatic energy contribution (Debye-Huckel)
            if self.electrostatic:
                # Define cutoff for the interaction (2 charged atoms within 20 Å)
                cutoff = (q_ij != 0) & (d < 20)

                # Calculate electrostatic energy for every atoms pair
                pair_dh_energy = q_ij[cutoff] * 7.0 / d[cutoff] * np.exp(-d[cutoff] / 10.0)

                # Calculate total electrostatic energy for the rotamer pose (= sum of pairwise energies)
                dh_energy_pose[rotamer_counter] = pair_dh_energy.sum()

            # Potential energy contribution (Lennard-Jones)
            if self.potential == 'lj':
                # Define cutoff for the interaction (2 atoms within 10 Å)
                cutoff = d < 10

                # Calculate the attraction part of the LJ potential for every atoms pair
                d = np.power(rmin_ij[cutoff] / d[cutoff], 6)

                # Calculate the total LJ potential energy for the rotamer pose (= sum of pairwise energies)
                pot_energy_pose[rotamer_counter] = np.sum(eps_ij[cutoff] * (d * d - 2. * d))

            # Potential energy contribution (Gaussian)
            if self.potential == 'gauss':
                # Define cutoff for the interaction (2 atoms within 10 Å)
                cutoff = d < 10

                # Calculate the exponent of Gaussian potential for every atoms pair
                d = np.power(d[cutoff] / rmin_ij[cutoff], 2)

                # Calculate the total Gaussian potential energy (mean = 0) for the rotamer pose
                # (= sum of pairwise energies)
                pot_energy_pose[rotamer_counter] = np.sum(eps_ij[cutoff] * np.exp(-0.5 * d))

        # Calculate Boltzmann weights array (one for each trajectory frame)
        boltzmann_factor = np.exp(-pot_energy_pose / (gas_un * self.temp) - dh_energy_pose)

        return boltzmann_factor

    def rotamerWeights(self, rotamersSite, lib, LJ_data):

        """

        Calculate Boltzmann distribution for each rotamer library

        Parameters
        ==========

            rotamersSite: MDAnalysis.Universe
                Universe object of the Chromophore rotamer library placed on the protein residue (corresponds to
                fitted_rotamers in lj_calculation)

            lib: dict
                Rotamers library of the chromophore

            LJ_data: list
                List containing data for LJ calculations (proteinNotSite.indices, rotamerSel_LJ.indices, eps_ij, rmin_ij, q_ij)

        Return
        ======

            prob: numpy.array
                Boltzmann distribution for the rotamer library

            Z: float
                Boltzmann partition function

        """

        # Calculate Boltzmann weights
        boltz = self.lj_calculation(rotamersSite, LJ_data)

        # Save external probabilities
        # np.savetxt(self.output_prefix+'-boltz-{:s}.dat'.format(residue_sel.replace(" ", "")),boltz)

        # Set to zero Boltzmann weights that are NaN
        boltz[np.isnan(boltz)] = 0.0

        # Multiply Boltzmann weights by library weights (= population of each rotamer)
        boltz = lib.weights * boltz

        # Calculate Boltzmann partition function
        Z = np.sum(boltz)

        # Calculate Boltzmann distribution
        prob = boltz / Z

        return prob, Z

    def weightedAvgSDSE(self, values, weights):

        """

        Calculate the weighted average and standard deviation.

        Parameters
        ==========

            values: numpy.array
                Array of value to be averaged

            weights: numpy.array
                Array of weights to compute weighted average

        Return
        ======

            (avg, std, sem) = (float, float, float)
                Weighted average, Standard deviation and Standard Error of the Mean

        """

        # Calculate weighted average
        avg = np.average(values, weights=weights)

        # Calculate weighted variance
        variance = np.average((values - avg) ** 2, weights=weights)

        # Standard deviation
        std = np.sqrt(variance)

        # Standard error of the mean
        sem = np.sqrt(variance / values.size)

        return (avg, std, sem)

    def calculate_ws(self):

        """ Calculate per-frame weights for reweighting """

        Z = np.genfromtxt(self.output_prefix + '-Z-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]),
                          skip_header=-1,
                          skip_footer=0,
                          delimiter=' ')

        # Check if single-frame structure
        if np.shape(Z) == (2,):

            w_s = np.array([1.0])

        else:

            Z_s = Z[:, 0] * Z[:, 1]
            w_s = Z_s / np.sum(Z_s)

        return w_s

    def fraction_frames(self):

        """ Compute effective fraction of frames contributing to the averages """

        if isinstance(self.weights, np.ndarray):

            if np.array_equal(self.weights, np.ones(len(self.protein.trajectory))):

                w_s = np.genfromtxt(
                    self.output_prefix + '-w_s-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]))

            else:

                w_s = self.weights

        elif not self.weights:

            w_s = np.genfromtxt(
                self.output_prefix + '-w_s-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]))

        ws_correct = w_s[w_s != 0]
        ws_0 = np.zeros(len(ws_correct))
        ws_0[:] = 1 / len(ws_correct)

        S = - np.sum(ws_correct * np.log(ws_correct / ws_0[ws_0 != 0]))

        phi = np.exp(S)

        return phi
