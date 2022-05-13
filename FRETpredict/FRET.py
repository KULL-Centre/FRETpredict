# -*- coding: utf-8 -*-
"""
FRET prediction Class
=====================

Class to perform Fluorescence Resonance Energy Transfer calculations given rotamer libraries for the fluorophores.

"""

# Coordinates and arrays
import os
import numpy as np
import pandas as pd
import h5py
import logging
import MDAnalysis
import re

# Inner imports
from .utils import Operations


class FRETpredict(Operations):
    """

    Calculation of FRET signal between two chromophores.

    Attributes
    ==========

        Operations.chains: str
            Name of the protein chains in which the placement residues are located

        Operations.protein: MDAnalysis.Universe
            Protein MDAnalysis Universe object

        residues: list of str
            Placement residue number

        donor: int
            donor chromophore number

        acceptor: int
            acceptor chromophore number

        r0lib: str
            path to chromophore R0 data

        z_cutoff: float
            Cutoff value for Boltzmann partition function, corresponding to the threshold for states population.

        fixed_R0: bool
            Decide if performing analytical R0 calculation or provide fixed value through the r0 attribute.

        r0 = float
            Forster distance for the chromophores pair

        dr: float
            distance different between different bins in distance distribution calculation

        rmin: int
            Inferior limit for distance distribution calculation

        rmax: int
            Superior limit for distance distribution calculation

        nr: int
            number of bins for distance distribution calculation

        rax: numpy.array
            array of distance values of bins for distance distribution calculation

        output_prefix: str
            Prefix used in saving data to file

        load_file: Bool
            Load data from pre-existing file

        weights: numpy.array
            Boltzmann factors distribution

        stdev: float
            Standard deviation

    Methods
    =======

        compute_chromophore_vectors:
            Calculate chromophores transition dipole moment vector and center-to-center distance distribution for k2
            calculations.

        calculateR0:
            Calculate FRET R0 between the donor and acceptor chromophores.

        trajectoryAnalysis:
            Calculate distribution of <E> (i.e. one <E> for each protein trajectory frame) in static, dynamic1, and
            dynamic2 regimes.

        save:
            Calculate k2 distribution and k2, Static, Dynamic1, Dynamic2 averaging. Save data to file.

        run:
            Run FRET Efficiency calculations by calling trajectoryAnalysis() and save() functions or by loading
            pre-computed data from file.

    """

    def __init__(self, protein, residues, **kwargs):

        """

        Parameters
        ==========

            protein: MDAnalysis.Universe
                Protein MDAnalysis Universe object

            residues: list of str
                Placement residue number

        """

        # Call the constructor of the Operations class
        Operations.__init__(self, protein, **kwargs)

        # Attributes assignments
        self.residues = residues
        self.donor = kwargs.get('donor', 488)
        self.acceptor = kwargs.get('acceptor', 532)
        self.r0lib = kwargs.get('r0lib', 'lib/R0')
        self.z_cutoff = kwargs.get('z_cutoff', 0.05)
        self.fixed_R0 = kwargs.get('fixed_R0', False)
        self.r0 = kwargs.pop('r0', 5.4)
        self.dr = 0.05
        self.rmin = -5
        self.rmax = kwargs.get('rmax', 20) * 2 - self.rmin
        self.nr = int(round((self.rmax - self.rmin) / self.dr, 0) + 1)
        self.rax = np.linspace(self.rmin, self.rmax, self.nr)
        self.output_prefix = kwargs.get('output_prefix', 'res')
        self.load_file = kwargs.get('load_file', False)
        self.weights = kwargs.get('weights', False)
        self.stdev = kwargs.get('filter_stdev', 0.02)

        # Logging set up
        logging.basicConfig(filename=kwargs.get('log_file', 'log'), level=logging.INFO)

        # Write the string for the placement residues selection
        for i in range(2):

            residue_sel = "resid {:d}".format(self.residues[i])

            # If chains are specified, add them to the atom selection
            if type(self.chains[i]) == str:
                residue_sel += " and segid {:s}".format(self.chains[i])

            # Logging information on the selected residues
            logging.info('{:s} = {:s}'.format(residue_sel, self.protein.select_atoms(residue_sel).atoms.resnames[0]))

        # Raise error if the specified placement residues are different than two
        if len(residues) != 2:
            raise ValueError("The residue_list must contain exactly 2 "
                             "residue numbers: current value {0}.".format(residues))

    def compute_chromophore_vectors(self, rotamersSite1, rotamersSite2):

        """

        Calculate chromophores transition dipole moment vector and center-to-center distance distribution for k2
        calculations.

        Parameters
        ==========

            rotamersSite1: MDA.Universe
                donor Universe

            rotamersSite2: MDA.Universe
                acceptor Universe

        Returns
        =======

            mu_cosine: numpy.array of float
                μ_si·μ_sj for k2 calculations

            r_mu_1: numpy.array of float
                R_sij·μ_sj for k2 calculations

            r_mu_2: numpy.array of float
                R_sij·μ_si for k2 calculations

            rdist: numpy.array of float
                chromophore distance distribution

        """

        # Select Chromophore atoms for transition dipole vector calculation
        rotamer1_c1 = rotamersSite1.select_atoms('name ' + self.lib_1.mu[0])
        rotamer2_c1 = rotamersSite2.select_atoms('name ' + self.lib_2.mu[0])
        rotamer1_c2 = rotamersSite1.select_atoms('name ' + self.lib_1.mu[1])
        rotamer2_c2 = rotamersSite2.select_atoms('name ' + self.lib_2.mu[1])

        # Select positions of Chromophore atoms for transition dipole vector
        c1_rot1_pos = np.array([rotamer1_c1.positions for x in rotamersSite1.trajectory]).squeeze(axis=1)
        c1_rot2_pos = np.array([rotamer2_c1.positions for x in rotamersSite2.trajectory]).squeeze(axis=1)
        c2_rot1_pos = np.array([rotamer1_c2.positions for x in rotamersSite1.trajectory]).squeeze(axis=1)
        c2_rot2_pos = np.array([rotamer2_c2.positions for x in rotamersSite2.trajectory]).squeeze(axis=1)

        # Array of Unit vectors for transition dipole of chromophore 1
        mu_1 = c2_rot1_pos - c1_rot1_pos
        mu_1 /= np.linalg.norm(mu_1, axis=1, keepdims=True)

        # Array of Unit vectors for transition dipole of chromophore 2
        mu_2 = c1_rot2_pos - c2_rot2_pos
        mu_2 /= np.linalg.norm(mu_2, axis=1, keepdims=True)

        # Calculate μ_si·μ_sj for k2 calculations
        mu_cosine = np.einsum('ik,jk->ij', mu_1, mu_2)

        # Calculation of the vector uniting donor->acceptor centers (R_sij)
        # Select chromophore centers atoms for each rotamer (specific atoms defined in rotamer library)
        rotamer1centre = rotamersSite1.select_atoms('name ' + self.lib_1.r[0])
        rotamer2centre = rotamersSite2.select_atoms('name ' + self.lib_2.r[0])

        # Obtain coordinates of chromophore centers atom for each rotamer
        centre1_pos = np.array([rotamer1centre.positions for x in rotamersSite1.trajectory]).squeeze(axis=1)
        centre2_pos = np.array([rotamer2centre.positions for x in rotamersSite2.trajectory]).squeeze(axis=1)

        # Calculate R_sij vector and unitary R_sij vector
        rvec = centre1_pos[:, np.newaxis, :] - centre2_pos
        rdist = np.linalg.norm(rvec, axis=2, keepdims=True)
        runitvec = rvec / rdist

        # Convert distance in nm
        rdist = rdist.flatten() / 10

        # Calculate R_sij·μ_sj, and R_sij·μ_si for k2 calculations
        r_mu_1 = np.einsum('ijk,ik->ij', runitvec, mu_1)
        r_mu_2 = np.einsum('ijk,jk->ij', runitvec, mu_2)

        # Different method for calculating transition dipole moments:
        # Mean between C1-C2 vectors
        # r_12 = c2_rot1_pos[:, np.newaxis, :] - c1_rot2_pos
        # r_12 /= np.linalg.norm(r_12, axis=2, keepdims=True)
        # r_21 = c1_rot1_pos[:, np.newaxis, :] - c2_rot2_pos
        # r_21 /= np.linalg.norm(r_21, axis=2, keepdims=True)
        # trans_dip_moment12_rot1 = np.einsum('ijk,ik->ij', r_12, mu_1)
        # trans_dip_moment21_rot1 = np.einsum('ijk,ik->ij', r_21, mu_1)
        # trans_dip_moment_rot1 = (trans_dip_moment12_rot1 + trans_dip_moment21_rot1) / 2.
        # trans_dip_moment12_rot2 = np.einsum('ijk,jk->ij', r_12, mu_2)
        # trans_dip_moment21_rot2 = np.einsum('ijk,jk->ij', r_21, mu_2)
        # trans_dip_moment_rot2 = (trans_dip_moment12_rot2 + trans_dip_moment21_rot2) / 2.

        return mu_cosine, r_mu_1, r_mu_2, rdist

    def calculateR0(self, k2):

        """

        Calculate FRET R0 between the donor and acceptor chromophores.

        Parameters
        ==========

            k2: float
                average orientation factor between two chromophores.

        """

        # Extract donor and acceptor producer and numbers from string
        temp = re.compile("([a-zA-Z]+) ([0-9-a-zA-Z]+)")
        donor_producer = temp.match(self.donor).groups()[0]
        acceptor_producer = temp.match(self.acceptor).groups()[0]

        donor_number = temp.match(self.donor).groups()[1]
        acceptor_number = temp.match(self.acceptor).groups()[1]

        # Read donor spectrum from file and normalize max value to 1
        donor_spectrum = pd.read_csv(f'{self.r0lib}/{donor_producer}{donor_number}.csv')
        donor_spectrum[['Emission', 'Excitation']] = donor_spectrum[['Emission', 'Excitation']] / 100

        # Read acceptor spectrum from file and normalize max value to 1
        acceptor_spectrum = pd.read_csv(f'{self.r0lib}/{acceptor_producer}{acceptor_number}.csv')
        acceptor_spectrum[['Emission', 'Excitation']] = acceptor_spectrum[['Emission', 'Excitation']] / 100

        # Quantum yield and extinction coefficient data for the chromophores
        chromophore_data = pd.read_csv(f'{self.r0lib}/Dyes_extinction_QD.csv', delimiter=',', on_bad_lines='skip',
                                       names=['Type', 'Chromophore', 'Ext_coeff', 'QD'])

        # R0 calculation parameters
        # Initial factor, for R0 expressed in nm
        factor = 0.02108

        # 4th-power of the medium refractive index (water)
        n4 = 1.4 ** 4

        # Quantum yield of the donor in the acceptor absence
        QD = float(chromophore_data['QD'].loc[(chromophore_data['Chromophore'] == donor_number) &
                                              (chromophore_data['Type'] == donor_producer)])

        # Extinction coefficient of the acceptor at its peak absorption value (= max value)
        ext_coeff_max = float(chromophore_data['Ext_coeff'].loc[(chromophore_data['Chromophore'] == acceptor_number) &
                                                                (chromophore_data['Type'] == acceptor_producer)])

        # Extinction coefficient spectrum of the acceptor
        ext_coeff_acceptor = (ext_coeff_max * acceptor_spectrum['Excitation']).fillna(0)

        # Integral of the donor emission spectrum
        donor_spectra_integral = np.trapz(donor_spectrum['Emission'], x=donor_spectrum['Wavelength'])

        # Overlap integral between donor-acceptor (normalized by the donor emission spectrum)
        J = np.trapz(donor_spectrum['Emission'] * ext_coeff_acceptor * donor_spectrum['Wavelength'] ** 4,
                     x=donor_spectrum['Wavelength']) / donor_spectra_integral

        # Forster R0 distance between the donor-acceptor pair, in nm
        self.r0 = factor * np.power(k2 * QD / n4 * J, 1 / 6)

    def trajectoryAnalysis(self):

        """

        Calculate distribution of <E> (i.e. one <E> for each protein trajectory frame) in static and dynamic regimes

        """

        # Print info
        logging.info("Starting rotamer distance analysis of trajectory {:s} with labeled residues "
                     "{:d} and {:d}".format(self.protein.trajectory.filename, self.residues[0], self.residues[1]))

        # Create H5PY file
        f = h5py.File(self.output_prefix + '-{:d}-{:d}.hdf5'.format(self.residues[0], self.residues[1]), "w")

        # Initialize a H5PY dataset named "distributions", with shape = (n_frames, rax.size))
        distributions = f.create_dataset("distributions", (self.protein.trajectory.n_frames, self.rax.size),
                                         fillvalue=0, compression="gzip")

        # Select placement residue atoms, compute Lennard-Jones and Electrostatic (Debye-Huckel) parameters for
        # Chromophore rotamers and Protein
        rotamer1, prot_atoms1, LJ_data1 = self.precalculate_rotamer(self.residues[0], self.chains[0], self.lib_1)
        rotamer2, prot_atoms2, LJ_data2 = self.precalculate_rotamer(self.residues[1], self.chains[1], self.lib_2)

        # For each trajectory frame of the protein, place the probes at the spin-labeled site using rotamer_placement(),
        # calculate Boltzmann weights based on Lennard-Jones interactions and calculate weighted distributions of
        # probe-probe separations

        # Variable allocation
        zarray = np.empty(0)  # Array of Boltzmann partition functions (sum over Boltzmann weights)
        k2_avg = np.full(self.protein.trajectory.n_frames, np.nan)
        esta_avg = np.full(self.protein.trajectory.n_frames, np.nan)
        edyn1_avg = np.full(self.protein.trajectory.n_frames, np.nan)
        edyn2_avg = np.full(self.protein.trajectory.n_frames, np.nan)
        allk2 = np.empty(0)
        allZ = np.empty(0)

        for frame_ndx, _ in enumerate(self.protein.trajectory):

            print(f'\nFrame {frame_ndx + 1}/{len(self.protein.trajectory)}')

            # Fit the rotamers onto the protein
            # Each protein structure (i.e. trajectory frame) has m conformations for chromophore 1 and l conformations
            # for chromophore 2
            rotamersSite1 = self.rotamer_placement(rotamer1, prot_atoms1, self.lib_1)
            rotamersSite2 = self.rotamer_placement(rotamer2, prot_atoms2, self.lib_2)

            # Calculate Boltzmann weights based on Lennard-Jones interactions
            boltz1, z1 = self.rotamerWeights(rotamersSite1, self.lib_1, LJ_data1)
            boltz2, z2 = self.rotamerWeights(rotamersSite2, self.lib_2, LJ_data2)

            # Calculate normalized Boltzmann weights (combined probability of the two dyes for each frame)
            boltzmann_weights_norm = (boltz1.reshape(-1, 1) * boltz2).flatten()

            # Create array of partition function values
            allZ = np.append(allZ, boltzmann_weights_norm.flatten())

            # Append Boltzmann partition functions of the two dyes for the frame to the array
            zarray = np.append(zarray, [z1, z2])

            # Compare Boltzmann partition function with cutoff
            if (z1 <= self.z_cutoff) or (z2 <= self.z_cutoff):
                # Warning for Z < Z_cutoff
                print('\nZ < Z_cutoff')

                # If Z value < cutoff then create an empty array for k2 values with same dimension as Z array, to save
                allk2 = np.zeros_like(allZ, dtype=float)

                # Skip to next iteration
                continue

            # Calculate vector factors for k2 calculations
            mu_cosine, r_mu_1, r_mu_2, rdist = self.compute_chromophore_vectors(rotamersSite1, rotamersSite2)

            # Calculate Orientation factor k2 for each donor-acceptor rotamer combination
            k2 = np.power(mu_cosine - 3 * r_mu_1 * r_mu_2, 2).flatten()

            # Create array of k2 values
            allk2 = np.append(allk2, k2.flatten())

            # Orientation factor for dynamic averaging is weighted by rotamers probability (Boltzmann weights),
            # to account for reduced weight of high energy conformations

            # Calculate Average Orientation factor <k2>
            k2_avg[frame_ndx] = np.dot(k2, boltzmann_weights_norm)

            # Calculate R0 for the donor-acceptor pair, based on computed k2, if R0 calculations are enabled
            if self.fixed_R0 == False:

                self.calculateR0(k2_avg[frame_ndx])

                # If dyes pair has no spectral overlap, skip iteration
                if self.r0 == 0:

                    print('\nDyes pair R0 = 0!')
                    continue

            # Calculate (r/r0)^6 factor for FRET efficiency calculations
            ratio6 = np.power(rdist / self.r0, 6)

            # Calculation of average FRET efficiencies distributions
            # Static Regime
            estatic = 1. / (1 + 2 / 3. * np.divide(ratio6, k2))
            esta_avg[frame_ndx] = np.dot(estatic, boltzmann_weights_norm)

            # Dynamic1 Regime
            edynamic = 1. / (1 + 2 / 3. / k2_avg[frame_ndx] * ratio6)
            edyn1_avg[frame_ndx] = np.dot(edynamic, boltzmann_weights_norm)

            # Dynamic2 Regime
            A_avg = np.dot(3. / 2. * k2 / ratio6, boltzmann_weights_norm)
            edyn2_avg[frame_ndx] = A_avg / (A_avg + 1)

            # Calculate normalized rdist between 0 and 1
            rdist = np.round((self.nr * (rdist - self.rmin)) / (self.rmax - self.rmin)).astype(int).flatten()

            # Calculate weighted distribution of center-to-center distances between the two chromophores
            distribution = np.bincount(rdist, weights=boltzmann_weights_norm.flatten(), minlength=self.rax.size)

            # Write the calculated distribution for each frame in the H5PY dataset
            try:

                distributions[frame_ndx] = distribution

            except:

                continue

        # Close H5PY file
        f.close()

        # Save distributions to file
        # Save k2 weighted distribution
        np.savetxt(self.output_prefix + '-Pk2-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]),
                   np.c_[np.arange(0.0, 4.02, .04)[:-1] + .02,
                         np.histogram(allk2, bins=np.arange(0.0, 4.02, .04), density=True, weights=allZ)[0]])

        # Save Partition function distribution
        np.savetxt(self.output_prefix + '-Z-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]),
                   zarray.reshape(-1, 2))

        # Save <k2> distribution
        np.savetxt(self.output_prefix + '-k2-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]), k2_avg)

        # Save <E>_static distribution
        np.savetxt(self.output_prefix + '-Es-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]), esta_avg)

        # Save <E>_dynamic1 distribution
        np.savetxt(self.output_prefix + '-Ed1-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]), edyn1_avg)

        # Save <E>_dynamic2 distribution
        np.savetxt(self.output_prefix + '-Ed2-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]), edyn2_avg)

    def save(self):

        """

        Calculate k2 distribution and k2, Static, Dynamic1, Dynamic2 averaging. Save data to file.

        """

        # Load <k2> distribution from file
        k2 = np.loadtxt(self.output_prefix + '-k2-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]))

        # Check if weights is an array
        if isinstance(self.weights, np.ndarray):

            # Check if every k2 instance is associated with one weight
            if self.weights.size != k2.size:
                # Warning log
                logging.info('Weights array has size {} whereas the number of frames is {}'.
                             format(self.weights.size, k2.size))

                # Raise error
                raise ValueError('Weights array has size {} whereas the number of frames is {}'.
                                 format(self.weights.size, k2.size))

        # If we want the average to be unweighted
        elif self.weights == False:

            # Associate every k2 instance with a unitary weight
            self.weights = np.ones(k2.size)

        # Other options will raise errors
        else:

            # Warning log
            logging.info('Weights argument should be a numpy array')

            # Raise error
            raise ValueError('Weights argument should be a numpy array')

        # Read data from H5PY file
        f = h5py.File(self.output_prefix + '-{:d}-{:d}.hdf5'.format(self.residues[0], self.residues[1]), "r")

        # Read distribution data from H5PY file
        distributions = f.get('distributions')

        # Calculate the sum of the weighted distribution
        distribution = np.nansum(distributions * self.weights.reshape(-1, 1), 0)

        # Calculate smoothed distance distribution
        frame_inv_distr = np.fft.ifft(distribution) * np.fft.ifft(np.exp(-0.5 * (self.rax / self.stdev) ** 2))
        smoothed = np.real(np.fft.fft(frame_inv_distr))
        smoothed /= np.trapz(smoothed, self.rax)

        # Save distance distribution to file
        np.savetxt(self.output_prefix + '-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]),
                   np.c_[self.rax[100:-100], smoothed[200:]], header='distance distribution')

        # Close H5PY file
        f.close()

        # Read FRET Efficiency data from file
        estatic = np.loadtxt(self.output_prefix + '-Es-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]))
        edynamic1 = np.loadtxt(self.output_prefix + '-Ed1-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]))
        edynamic2 = np.loadtxt(self.output_prefix + '-Ed2-{:d}-{:d}.dat'.format(self.residues[0], self.residues[1]))

        # k2, Static, Dynamic1, and Dynamic2 averaging
        # Create DataFrame if only one value of k2 is present (only 1 protein structure/frame provided)
        if k2.size == 1:

            df = pd.Series([k2, estatic, edynamic1, edynamic2], index=['k2', 'Estatic', 'Edynamic1', 'Edynamic2'])

        # Create DataFrame of weighted averaged values
        else:

            df = pd.DataFrame([self.weightedAvgSDSE(k2[np.isfinite(k2)], self.weights[np.isfinite(k2)]),
                               self.weightedAvgSDSE(estatic[np.isfinite(k2)], self.weights[np.isfinite(k2)]),
                               self.weightedAvgSDSE(edynamic1[np.isfinite(k2)], self.weights[np.isfinite(k2)]),
                               self.weightedAvgSDSE(edynamic2[np.isfinite(k2)], self.weights[np.isfinite(k2)])],
                              columns=['Average', 'SD', 'SE'], index=['k2', 'Estatic', 'Edynamic1', 'Edynamic2'])

        # Save DataFrame in pickle format
        df.to_pickle(self.output_prefix + '-data-{:d}-{:d}.pkl'.format(self.residues[0], self.residues[1]))

    def run(self, **kwargs):

        """

        Run FRET Efficiency calculations by calling trajectoryAnalysis() and save() functions or by loading
        pre-computed data from file.

        **kwargs:
        ========

            data_filepath: str
                Path of the data file with pre-computed data

        """

        data_filepath = kwargs.get('data_filepath', '')

        # If a file with pre-computed data is already present
        if self.load_file:

            # If the file path is correct
            if os.path.isfile(data_filepath):

                # Info log
                logging.info(
                    'Loading pre-computed data from {} - will not load trajectory file.'.format(data_filepath))

            # If the file path is not correct
            else:

                # Warning log
                logging.info('File {} not found!'.format(data_filepath))

                # Raise error
                raise FileNotFoundError('File {} not found!'.format(data_filepath))

            # Calculate k2 distribution and k2, Static, Dynamic1, Dynamic2 averaging. Save data to file.
            self.save()

        # If a file with pre-computed data isn't already present
        else:

            # Calculate distribution of <E> (i.e. one <E> for each protein trajectory frame) in static, dynamic1, and
            # dynamic2 regimes.
            self.trajectoryAnalysis()

            # Calculate k2 distribution and k2, Static, Dynamic1, Dynamic2 averaging. Save data to file.
            self.save()

            print('\nDone.')
