import MDAnalysis
import numpy as np
import pandas as pd
import mdtraj as md
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
from scipy.signal import find_peaks
from scipy.cluster.vq import kmeans2


class RotamerClusters(object):
    """

    Calculation of a rotamer library starting from a dye+linker trajectory.

    Attributes
    ==========

        libpath: str
            Path to folder where dye+linker simulation data is stored.

        path: str
            Path to folder where rotamer library data will be saved.

        cutoff: list of int
            List of cluster population used to filter conformations.

        dye: str
            Name of the dye+linker combination, as written in the DataFrame.

        df: pandas.DataFrame
            index: dye
            columns: 'indices' -> indices of the atoms forming the dihedrals
                     'peaks' -> peak dihedral values


    Methods
    =======

        calcDihe:
            Calculate dihedral angles on the dye+linker trajectory.

        genPeaks:
            Compute linker dihedral peaks.

        genClusters:
            1- Generate combinations of dihedral angles from the peaks (cluster centers C1).
            2- K-means Clustering
                - Assign each trajectory frame to the cluster center C1 of least square deviation.
                - Calculate the average over the dihedral angles that were assigned to the same cluster center.
                    This results in a set of new centers (C2).
                - Assign each trajectory frame to the cluster center C2 of least square deviation
            3- Find the trajectory frame that best represents the cluster center.

        filterCluster:
            1- Filter the cluster centers C2 based on a cutoff on N (cluster population).
                This results in a different number of centers (C3).
            2- Reassign the discarded frames to the remaining C3 cluster center of least square deviation.

        genRotLib:
            Translate + Rotate C3 cluster centers conformations, and write data to file.

        plotClustHist:
            Plot Dihedral distribution and peaks with cluster centers C3 dihedrals.

        plotClustPolar:
            Plot Dihedral distribution and peaks with cluster centers C3 dihedrals on a polar plot.

        run:
            Run all the calculations to generate a rotamer library from a dye+linker trajectory.

    """

    def __init__(self, **kwargs):

        self.libpath = kwargs.get('libpath', 'lib/')
        self.path = kwargs.get('path', 'lib/genLIB2/')
        self.cutoff = kwargs.get('cutoff', [50])
        self.dye = kwargs.get('dye', 'A48_C1R')
        self.df = kwargs.get('df', pd.DataFrame({'indices': [], 'peaks': []}).T)
        self.traj_extension = kwargs.get('traj_extension', 'xtc')

    def calcDihe(self):

        """ Calculate dihedral angles on the dye+linker trajectory. """

        # Load dye+linker trajectory
        if self.traj_extension == 'xtc':
            t = md.load_xtc(self.libpath + self.dye + '/traj.xtc', self.libpath + self.dye + '/conf_ed.gro')

        if self.traj_extension == 'dcd':
            t = md.load_dcd(self.libpath + self.dye + '/traj.dcd', self.libpath + self.dye + '/conf_ed.gro')

        elif self.traj_extension == 'pdb':
            t = md.load_pdb(self.libpath + self.dye + '/traj.pdb')

        # Calculate dihedrals from atom indices
        dihe = md.compute_dihedrals(t, self.df.loc[self.dye, 'indices'], periodic=True, opt=True) / np.pi * 180

        # Save data to text file
        np.savetxt(self.path + 'dihedrals/' + self.dye + '.txt', dihe)

    def genPeaks(self):

        """ Compute linker dihedral peaks. """

        # Load computed dihedrals from file
        dihe = np.loadtxt(self.path + 'dihedrals/' + self.dye + '.txt')

        peaks = []

        fig, axes = plt.subplots(np.shape(dihe)[1], figsize=(50, 20), sharex=True, sharey=True)

        # Iterate for every subplot (= dihedral angle)
        for i, ax in enumerate(axes.flatten()):

            # Calculate dihedral angle distribution
            h, b = np.histogram(dihe[:, i], bins=np.arange(-180, 181, 1), density=True)

            # Rebinning
            bins = b[:-1] + (b[1] - b[0]) / 2

            # Find dihedral distribution peaks
            p, properties = find_peaks(h, prominence=.0005, width=2, distance=60)

            # Plot the line with the distribution
            ax.plot(bins, h, linewidth=10)

            if h.argmax() not in p:
                p = np.append(p, h.argmax())

            # Print bin number corresponding to peak values, and corresponding dihedral angle
            print(f'Dihedral {i + 1}')
            print(f'Peak indices {p} \n peak dihedral angles {bins[p]}')

            # List of peak dihedral angles
            peaks.append(bins[p])

        # Do not display the plot
        plt.close(fig)

        # Append peak dihedral angles to the dye+linker dataframe

        self.df.loc[self.dye, 'peaks'] = peaks  # np.array(peaks)

    def genClusters(self):

        """

        1- Generate combinations of dihedral angles from the peaks (cluster centers C1).

        2- K-means Clustering
            - Assign each trajectory frame to the cluster center C1 of least square deviation.
            - Calculate the average over the dihedral angles that were assigned to the same cluster center.
                This results in a set of new centers (C2).
            - Assign each trajectory frame to the cluster center C2 of least square deviation

        3- Find the trajectory frame that best represents the cluster center.

        """

        # Load dihedral data from file
        dihedrals = np.loadtxt(self.path + 'dihedrals/' + self.dye + '.txt').astype(np.float16)

        # -1-
        # Generate all possible combinations for dihedral angle peaks (cluster centers C1)
        peaks = np.array(list(itertools.product(*self.df.loc[self.dye, 'peaks']))).astype(np.float16)

        # Print number of combinations and number of dihedral angles
        print(f'Peak combinations (C1): {peaks.shape[0]}')

        # -2-
        # Cluster centers C2
        cluster_dict = {}

        # k-means clustering, using peaks as the initial number of clusters, assigns every
        # trajectory frame to a centroid.
        # labels.shape = (num_frames, 1)
        centroids, labels = kmeans2(dihedrals.astype(float), peaks.astype(float), minit='matrix')

        # Number of cluster centers C2
        print(f'C2 Cluster centers: {len(np.unique(labels))}')

        # -3-
        # Iterate on trajectory frame
        for frame_index, centroid_id in enumerate(labels):

            # Sum of square distance values of the rotamer conformation to its centroid
            sum_sqdist = np.sum((dihedrals[frame_index] - centroids[centroid_id]) ** 2)

            # Add C2 centroid to the dictionary
            # If not already present
            if centroid_id not in cluster_dict.keys():

                # 'frame': Frame index of the closest trajectory conformation to the centroid
                # 'square_dist': Least sum of square distance between rotamers and centroid
                # 'dihe': Dihedral angles of the closest rotamer to the centroid
                # 'N': Number of frames belonging to the centroid (centroid population)
                cluster_dict[centroid_id] = {'frame': frame_index,
                                             'square_dist': sum_sqdist,
                                             'dihe': dihedrals[frame_index],
                                             'N': 1}
            # If already present
            else:

                # If square deviation is lower that that already present
                if sum_sqdist < cluster_dict[centroid_id]['square_dist']:
                    cluster_dict[centroid_id]['frame'] = frame_index
                    cluster_dict[centroid_id]['square_dist'] = sum_sqdist
                    cluster_dict[centroid_id]['dihe'] = dihedrals[frame_index]
                    cluster_dict[centroid_id]['N'] += 1

                else:
                    cluster_dict[centroid_id]['N'] += 1

        # Save C2 data to pickle file
        pd.DataFrame(cluster_dict).T.to_pickle(self.path + 'clusters_{:s}_1step.pkl'.format(self.dye))

    def filterCluster(self, cutoff):

        """

        1- Filter the cluster centers C2 based on a cutoff on N (cluster population).
            This results in a different number of centers (C3).

        2- Reassign the discarded frames to the remaining C3 cluster center of least square deviation.

        Parameters
        ==========

            cutoff: int
                cluster population used to filter conformations.

        """

        # Read C2 data from pickle file
        clusters = pd.read_pickle(self.path + '/clusters_{:s}_1step.pkl'.format(self.dye))

        # Print cutoff
        print(f'cutoff: {cutoff}')

        # -1-
        # Discard sparsely populated C3 cluster centers, based on cutoff
        # Get frames of discarded cluster centers
        discarded_frames = clusters[clusters.N < cutoff].frame.values.astype(int)

        # Get C1 peak index of discarded cluster centers
        discarded_peaks = clusters[clusters.N < cutoff].index

        # Get population of discarded cluster centers
        frequency = clusters[clusters.N < cutoff].N.values

        # Remove sparsely populated cluster centers
        clusters.drop(discarded_peaks, inplace=True)

        # Print number of clusters with population > cutoff
        print(f'C3 Final clusters: {len(clusters.index)}')

        # -2-
        # Load dihedral data from file
        dihedrals = np.loadtxt(self.path + 'dihedrals/' + self.dye + '.txt')

        # Dihedral angle values of the remaining cluster centers
        sel_dihe = np.array(clusters.dihe.tolist())

        # If the number of filtered clusters is zero
        if len(sel_dihe) == 0:

            # Cannot create a rotamer library with zero clusters
            print(f'\nNo cluster has population > {cutoff}, so the total number of clusters is 0!')

            raise ValueError

        # If there are no discarded frames
        elif len(discarded_frames) == 0:

            # Save cluster centers to pickle file
            clusters.to_pickle(self.path + 'clusters_{:s}_{:d}_cutoff.pkl'.format(self.dye, cutoff))

        # If there are discarded frames to reassign
        else:

            # Iterate on every frame associated with the discarded cluster centers C2
            for frame_index, frame_dihe in enumerate(dihedrals[discarded_frames]):

                # Compute square distance of each frame dihedral angle with cluster center C3
                sqdist = (frame_dihe - sel_dihe) ** 2

                # Compute sum of square distances of every dihedral angle for every peak combination
                sumleastsq = np.sum(sqdist, axis=1)

                # Assign every frame of discarded cluster center to closest cluster center C3
                clusters.iloc[sumleastsq.argmin()]['N'] += frequency[frame_index]

                # Save filtered cluster centers to pickle file
                clusters.to_pickle(self.path + 'clusters_{:s}_{:d}_cutoff.pkl'.format(self.dye, cutoff))

    def genRotLib(self, cutoff):

        """

        Translate + Rotate C3 cluster centers conformations, and write data to file.

        Parameters
        ==========

            cutoff: int
                cluster population used to filter conformations.

        """

        # Read C3 clusters data
        clusters = pd.read_pickle(self.path + 'clusters_{:s}_{:d}_cutoff.pkl'.format(self.dye, cutoff))

        # Create Universe for the dye+linker trajectory
        if self.traj_extension == 'xtc':
            u = MDAnalysis.Universe(self.libpath + self.dye + '/conf_ed.gro', self.libpath + self.dye + '/traj.xtc')

        if self.traj_extension == 'dcd':
            u = MDAnalysis.Universe(self.libpath + self.dye + '/conf_ed.gro', self.libpath + self.dye + '/traj.dcd')

        elif self.traj_extension == 'pdb':
            u = MDAnalysis.Universe(self.libpath + self.dye + '/traj.pdb')

        # Select dye+linker atoms
        chromophore = u.select_atoms('all and not (resname ACE or resname NHE)')

        # Write PDB for first frame
        chromophore.write(self.path + 'rot_lib_{:s}.pdb'.format(self.dye))

        # Select Ca, N, and C atoms of the linker, to position on the target protein residue
        Ca_pos = u.select_atoms('all and not (resname ACE or resname NHE) and name CA')
        N_pos = u.select_atoms('all and not (resname ACE or resname NHE) and name N')
        C_pos = u.select_atoms('all and not (resname ACE or resname NHE) and name C')

        new_coords = np.empty(0)

        # Iterate on the trajectory frames corresponding to the C3 cluster centers
        for _ in u.trajectory[clusters.frame.values]:
            # Use Ca as reference frame origin, and translate other atoms
            offset = Ca_pos.positions.copy()

            Ca_coords = Ca_pos.positions - offset
            N_coords = N_pos.positions - offset
            C_coords = C_pos.positions - offset
            chromophore_coords = chromophore.positions - offset

            # Create unitary x vector (Ca-N bond)
            x_vector = N_coords - Ca_coords
            x_vector /= np.linalg.norm(x_vector)

            # Create unitary y_t vector (Co-Ca bond, to obtain z vector)
            yt_vector = C_coords - Ca_coords
            yt_vector /= np.linalg.norm(yt_vector)

            # Create unitary z vector (perpendicular to plane formed by Ca-N and Co-Ca bonds)
            z_vector = np.cross(x_vector, yt_vector)
            z_vector /= np.linalg.norm(z_vector)

            # Create unitary y vector (perpendicular to plane formed by z vector and Ca-N bond)
            y_vector = np.cross(z_vector, x_vector)

            # Stack x, y, z vectors vertically to obtain rotation matrix
            rotation = np.array((x_vector, y_vector, z_vector)).T

            # Rotate chromophore rotamers using the rotation matrix
            chromophore_coords = np.dot(chromophore_coords, rotation.reshape(3, 3))
            chromophore_coords = chromophore_coords.reshape((len(chromophore_coords), 3))

            # Add new cluster center C3 coordinates to the array
            new_coords = np.append(new_coords, chromophore_coords)

        # Reshape and rescale the coordinates (MDAnalysis length unit is Å, MDTraj length unit is nm)
        new_coords = new_coords.reshape((len(clusters.frame.values), len(chromophore), 3)) / 10

        # Write cluster data (index, atom number, coordinates, cluster population) for each
        # atom to file
        output_file = open(self.path + 'rot_lib_matrix_{:s}_{:d}.txt'.format(self.dye, cutoff), 'w')

        for index, conformer in enumerate(new_coords):

            for atom_index, atom in enumerate(conformer):
                output_file.write(
                    '{0:>3} {1:>3} {2[0]:> 10.6f} {2[1]:> 10.6f} {2[2]:> 10.6f} {3:>5}\n'.format(index + 1,
                                                                                                 chromophore[
                                                                                                     atom_index].id -
                                                                                                 chromophore[1].id,
                                                                                                 atom,
                                                                                                 clusters.N.values[
                                                                                                     index]))

        output_file.close()

        # Load topology from 1st frame of dye+linker trajectory
        top = md.load(self.path + 'rot_lib_{:s}.pdb'.format(self.dye)).top

        # Number of frames of the trajectory
        n_frames = new_coords.shape[0]

        # Create trajectory for the new coordinates of the cluster centers
        traj = md.Trajectory(new_coords, top, unitcell_lengths=[[10, 10, 10]] * n_frames,
                             unitcell_angles=[[90, 90, 90]] * n_frames)

        # Save cluster centers trajectory to file
        traj.save_dcd(self.path + '{:s}_cutoff{:d}.dcd'.format(self.dye, cutoff))

        # Save first frame of the trajectory as PDB
        traj[0].save_pdb(self.path + '{:s}.pdb'.format(self.dye))

        # Save cluster populations as weights for FRET Efficiency calculations
        np.savetxt(self.path + '{:s}_cutoff{:d}_weights.txt'.format(self.dye, cutoff), clusters.N.values)

    def plotClustHist(self, dye, cutoff):

        """

        Plot Dihedral distribution and peaks with cluster centers C3 dihedrals.

        Parameters
        ==========

            dye: str
                Name of the dye+linker combination, as written in the DataFrame.

            cutoff: int
                cluster population used to filter conformations.

        """

        self.dye = dye

        # Read dihedral data from file
        dihe = np.loadtxt(self.path + 'dihedrals/' + self.dye + '.txt')
        num_dihedrals = np.shape(dihe)[1]

        # Read C3 filtered data from pickle file
        clusters_cutoff = pd.read_pickle(self.path + 'clusters_{:s}_{:d}_cutoff.pkl'.format(self.dye, cutoff))

        # Plot
        sns.set_style('darkgrid')

        fig, axes = plt.subplots(nrows=np.round(num_dihedrals / 3).astype(int), ncols=3,
                                 sharex=True, sharey=True, figsize=(9, 6))

        for i, ax in enumerate(axes.flatten()):

            if i == num_dihedrals:
                ax.set_visible(False)
                continue

            # Dihedral histogram
            h, b = np.histogram(dihe[:, i], bins=np.arange(-180, 181, 2), density=True)

            bins = b[:-1] + (b[1] - b[0]) / 2

            # Vertical lines corresponding to the cluster center dihedrals
            ax.vlines(np.array(clusters_cutoff.dihe.tolist())[:, i], ymin=0, ymax=h.max(), color='r', lw=0.5, ls=':')

            # Vertical lines corresponding to the dihedral peaks
            ax.vlines(self.df.loc[self.dye, 'peaks'][i], ymin=0, ymax=1, color='k')

            ax.plot(bins, h, lw=2)

            ax.set_ylim(0, h.max() + 0.01)

            ax.set_title("$\chi_" + '{' + f'{i + 1}' + '}$')

        # Set labels and titles
        for i in list(range(2, num_dihedrals, 3)) + list(range(0, num_dihedrals, 3)):
            axes.flatten()[i].set_ylabel(r'$P(\theta)$')

        for i in range(0, num_dihedrals):
            axes.flatten()[i].set_xlabel(r'$\theta$ / deg')

        for i in range(2, num_dihedrals, 3):
            axes.flatten()[i].yaxis.set_ticks_position('right')
            axes.flatten()[i].yaxis.set_label_position("right")

        fig.suptitle(self.dye.replace('_', ' ') + ' cutoff {:d}, {:d} rotamers'.format(cutoff,
                                                                                       len(clusters_cutoff.index)))

        plt.savefig(self.path + '/hist_{:s}_{:d}.pdf'.format(self.dye, cutoff))

        plt.tight_layout()

        plt.show()

    def plotClustPolar(self, dye, cutoff):

        """

        Plot Dihedral distribution and peaks with cluster centers C3 dihedrals on a polar plot.

        Parameters
        ==========

            dye: str
                Name of the dye+linker combination, as written in the DataFrame.

            cutoff: int
                cluster population used to filter conformations.

        """

        self.dye = dye

        # Read dihedral data from file
        dihe = np.loadtxt(self.path + 'dihedrals/' + self.dye + '.txt')
        num_dihedrals = np.shape(dihe)[1]

        # Read C3 filtered data from pickle file
        clusters_cutoff = pd.read_pickle(self.path + 'clusters_{:s}_{:d}_cutoff.pkl'.format(self.dye, cutoff))

        # Plot
        sns.set_style('darkgrid')

        fig, axes = plt.subplots(np.round(num_dihedrals / 3).astype(int), 3,
                                 sharex=False, sharey=False, subplot_kw=dict(polar=True), figsize=(9, 7))

        for i, ax in enumerate(axes.flatten()):

            if i == num_dihedrals:
                ax.set_visible(False)
                continue

            # Dihedral histogram
            h, b = np.histogram(dihe[:, i], bins=np.arange(-180, 181, 2), density=True)

            bins = b[:-1] + (b[1] - b[0]) / 2

            # Vertical lines corresponding to dihedral peaks
            ax.vlines((self.df.loc[self.dye, 'peaks'][i] + 180) / 180 * np.pi, ymin=0, ymax=h.max(), color='k')

            # Vertical lines corresponding to cluster centers C3
            ax.vlines((np.array(clusters_cutoff.dihe.tolist())[:, i] + 180) / 180 * np.pi, ymin=0, ymax=h.max(),
                      color='r', lw=0.5, ls=':')

            # Plot dihedral distribution in polar graph
            ax.plot((bins + 180) / 180 * np.pi, h, lw=2)

            # Plot settings
            ax.set_xticks(np.arange(0, 2 * np.pi, np.pi / 2))
            ax.set_title("$\chi_" + '{' + f'{i + 1}' + '}$')
            ax.set_yticks([])
            ax.grid(False)

        Nrotamers = len(clusters_cutoff.index)

        fig.suptitle(self.dye.replace('_', ' ') + ' cutoff {:d}, {:d} rotamers'.format(cutoff, Nrotamers))

        fig.tight_layout()

        fig.savefig(self.path + '/{:s}_{:d}.pdf'.format(self.dye, cutoff))

        plt.show()

    def run(self):

        """ Run all the calculations to generate a rotamer library from a dye+linker trajectory. """

        for dye_name in self.df.index:

            # Display dye+linker under analysis
            self.dye = dye_name
            print(f'Dye + Linker: {self.dye}')

            # Calculate dihedral angles for all the chromophore+linkers
            print("\nCalculating dihedral angles...")
            self.calcDihe()

            # Compute dihedral distributions and peaks
            print("\nComputing dihedral distributions and peaks...")
            self.genPeaks()

            # Generate C2 cluster centers
            print("\nGenerating cluster centers...")
            self.genClusters()

            # Filter C2 cluster centers for population, generate cluster centers C3,
            # save generated rotamer libraries.
            for _, co in enumerate(self.cutoff):

                try:

                    print(f"\nFiltering cluster centers...")
                    self.filterCluster(co)

                    print("\nGenerating rotamer library...")
                    self.genRotLib(co)

                    print('Done.\n')

                except ValueError:

                    print(f'\nCould not generate rotamer library for cutoff {co}')
