# -*- coding: utf-8 -*-
"""
FRETPrediction Class
--------------------

Class to perform Fluorescence Resonance Energy Transfer calculations given two positions for the fluorofores.

"""

import os

# Coordinates and arrays
import numpy as np
import MDAnalysis
import math
import MDAnalysis.analysis.distances as mda_dist

# Logger
import logging

# Inner imports
from DEERpredict.utils import Operations
import DEERpredict.libraries

logger = logging.getLogger("MDAnalysis.app")

class FRETpredict(Operations):
    """Calculation of FRET signal between two chromophores."""

    def __init__(self, protein_structure, residues, record_frames=False, **kwargs):
        """

        Args:
            protein_structure (:py:class:`MDAnalysis.core.universe.Universe`):
            residues (list(:py:class:`str`)):
        :Keywords:
        """

        Operations.__init__(self, protein_structure, **kwargs)
        self.libname_1 = kwargs.pop('libname_1', 'Alexa 488 150cutoff')
        self.libname_2 = kwargs.pop('libname_2', 'Alexa 594 150cutoff')
        self.lib_1 = DEERpredict.libraries.RotamerLibrary(self.libname_1)
        self.lib_2 = DEERpredict.libraries.RotamerLibrary(self.libname_2)
        self.k2_calc_static = kwargs.pop('k2_calc_static', False)
        self.k2_calc_dyn = kwargs.pop('k2_calc_dynamic', False)
        # loads all args and kwargs
        self.record_frames = record_frames
        self.r0 = kwargs.pop('r0', 5.4)
        #

        if self.chains:
            output_file, ext = os.path.splitext(kwargs.pop('output_file', 'distance_profile'))
            ext = ext or ".dat"
            self.output_file = "{0}-{1[0]}{2[0]}-{1[1]}{2[1]}{3}".format(output_file, residues, self.chains, ext)
            ext = ".png"
            self.output_plot = "{0}-{1[0]}{2[0]}-{1[1]}{2[1]}{3}".format(output_file, residues, self.chains, ext)
        else:
            output_file, ext = os.path.splitext(kwargs.pop('output_file', 'distance_profile'))
            ext = ext or ".dat"
            self.output_file = "{0}-{1[0]}-{1[1]}{2}".format(output_file, residues, ext)
            ext = ".png"
            self.output_plot = "{0}-{1[0]}-{1[1]}{2}".format(output_file, residues, ext)

        # Binning
        rax = np.linspace(0, 1, 100)
        distributions = np.zeros((self.replicas + 1, rax.size))

        if len(residues) != 2:
            raise ValueError("The residue_list must contain exactly 2 residue numbers: current "
                             "value {0}.".format(residues))

        logger.info("Starting rotamer distance analysis of trajectory "
                    "{0}...".format(protein_structure.trajectory.filename))
        logger.info("Rotamer library for position 1 = '{0}'".format(self.lib_1.name))
        logger.info("Rotamer library for position 2 = '{0}'".format(self.lib_2.name))
        logger.debug("Results will be written to {0}.".format(self.output_file))

        # Pre-processing weights
        lib1_norm = self.lib_1.weights / np.sum(self.lib_1.weights)
        lib2_norm = self.lib_2.weights / np.sum(self.lib_2.weights)

        self.frames_pre_replica = int((protein_structure.trajectory.n_frames - self.discard_frames) / self.replicas)
        progressmeter = MDAnalysis.lib.log.ProgressMeter(math.ceil((self.stop_frame-self.start_frame)/self.jump_frame),

                                                         interval=1)
        if self.k2_calc_dyn:
            self.kappa2 = []
            weights_list = []
            for protein in protein_structure.trajectory[
                           self.start_frame:self.stop_frame:self.jump_frame]:  # discard first on gromacs xtc
                if self.chains:
                    rotamersSite1 = self.rotamer_placement(self.lib_1.data,
                                                           protein_structure,
                                                           residues[0],
                                                           self.chains[0],
                                                           probe_library=self.lib_1)
                    rotamersSite2 = self.rotamer_placement(self.lib_2.data,
                                                           protein_structure,
                                                           residues[1],
                                                           self.chains[1],
                                                           probe_library=self.lib_2)

                else:
                    rotamersSite1 = self.rotamer_placement(self.lib_1.data,
                                                           protein_structure,
                                                           residues[0],
                                                           probe_library=self.lib_1)
                    rotamersSite2 = self.rotamer_placement(self.lib_2.data,
                                                           protein_structure,
                                                           residues[1],
                                                           probe_library=self.lib_2)
                boltz1 = self.lj_calculation(rotamersSite1, protein_structure, residues[0], fret=True)
                boltz1 = np.multiply(lib1_norm, boltz1)
                z_1 = np.nansum(boltz1)
                if z_1 == 0:
                    continue
                boltzman_weights_norm1 = boltz1 / z_1

                boltz2 = self.lj_calculation(rotamersSite2, protein_structure, residues[1], fret=True)
                boltz2 = np.multiply(lib2_norm, boltz2)
                z_2 = np.nansum(boltz2)
                if z_2 == 0:
                    continue
                boltzman_weights_norm2 = boltz2 / z_2

                # define the atoms to measure the distances between
                if self.libname_1.split()[1] == '488':
                    rotamer1_c11 = rotamersSite1.select_atoms("name C11")
                    rotamer2_c11 = rotamersSite2.select_atoms("name C17")
                    rotamer1_c12 = rotamersSite1.select_atoms("name C12")
                    rotamer2_c12 = rotamersSite2.select_atoms("name C18")
                else:
                    rotamer1_c11 = rotamersSite1.select_atoms("name C17")
                    rotamer2_c11 = rotamersSite2.select_atoms("name C11")
                    rotamer1_c12 = rotamersSite1.select_atoms("name C18")
                    rotamer2_c12 = rotamersSite2.select_atoms("name C12")
                c11_probe1_pos = np.array([i for x in rotamersSite1.trajectory for i in rotamer1_c11.positions])
                c11_probe2_pos = np.array([i for x in rotamersSite2.trajectory for i in rotamer2_c11.positions])

                c12_probe1_pos = np.array([i for x in rotamersSite1.trajectory for i in rotamer1_c12.positions])
                c12_probe2_pos = np.array([i for x in rotamersSite2.trajectory for i in rotamer2_c12.positions])

                inner_probe1_vect = c12_probe1_pos - c11_probe1_pos
                inner_probe1_vect /= np.linalg.norm(inner_probe1_vect, axis=1)[:, np.newaxis]

                inner_probe2_vect = c11_probe2_pos - c12_probe2_pos
                inner_probe2_vect /= np.linalg.norm(inner_probe2_vect, axis=1)[:, np.newaxis]

                for position1_index, position1 in enumerate(lib1_norm):
                        mutual_displacement = np.dot(inner_probe2_vect, inner_probe1_vect[position1_index].T)

                        outer_vect1 = c12_probe1_pos[position1_index] - c11_probe2_pos
                        outer_vect1 /= np.linalg.norm(outer_vect1, axis=1)[:, np.newaxis]
                        outer_vect2 = c11_probe1_pos[position1_index] - c12_probe2_pos
                        outer_vect2 /= np.linalg.norm(outer_vect2, axis=1)[:, np.newaxis]

                        trans_dip_moment1_probe1 = np.dot(outer_vect1, inner_probe1_vect[position1_index].T)
                        trans_dip_moment2_probe1 = np.dot(outer_vect2, inner_probe1_vect[position1_index].T)

                        for pos2_index, element in enumerate(lib2_norm):  # could we take out this for loop?
                            trans_dip_moment1_probe2 = np.dot(outer_vect1[pos2_index], inner_probe2_vect[pos2_index])
                            trans_dip_moment2_probe2 = np.dot(outer_vect2[pos2_index], inner_probe2_vect[pos2_index])

                            new_k2 = np.power((mutual_displacement[pos2_index] -
                                               ((3. / 4.) * (trans_dip_moment1_probe1[pos2_index] +
                                                             trans_dip_moment2_probe1[pos2_index]) * (
                                                    trans_dip_moment1_probe2 +
                                                    trans_dip_moment2_probe2))), 2)
                            self.kappa2.append(new_k2)
                            weights_list.append(boltzman_weights_norm1[position1_index] *
                                                boltzman_weights_norm2[pos2_index])
        if self.k2_calc_dyn:
            # self.kappa2 = np.average(self.kappa2)
            self.kappa2 = np.average(self.kappa2, weights=weights_list)

            # kappa2 = np.average(kappa2)
            logger.info('Dynamic K2 calculated: {:.3f}'.format(self.kappa2))
        if self.k2_calc_static:
            kappa2_hold = []
            weights_list = []
        for protein in protein_structure.trajectory[self.start_frame:self.stop_frame:self.jump_frame]:  # discard first on gromacs xtc
            progressmeter.echo(int(protein.frame/self.jump_frame))
            self.current_replica = ((protein.frame - self.discard_frames) // self.frames_pre_replica) # FIXME: not correct
            if self.chains:
                rotamersSite1 = self.rotamer_placement(self.lib_1.data,
                                                       protein_structure,
                                                       residues[0],
                                                       self.chains[0],
                                                       probe_library=self.lib_1)
                rotamersSite2 = self.rotamer_placement(self.lib_2.data,
                                                       protein_structure,
                                                       residues[1],
                                                       self.chains[1],
                                                       probe_library=self.lib_2)

            else:
                rotamersSite1 = self.rotamer_placement(self.lib_1.data,
                                                       protein_structure,
                                                       residues[0],
                                                       probe_library=self.lib_1)
                rotamersSite2 = self.rotamer_placement(self.lib_2.data,
                                                       protein_structure,
                                                       residues[1],
                                                       probe_library=self.lib_2)

            boltz1 = self.lj_calculation(rotamersSite1, protein_structure, residues[0], fret=True)
            # print 'boltz 1'
            # print boltz1
            boltz1 = np.multiply(lib1_norm, boltz1)
            z_1 = np.nansum(boltz1)
            # print z_1
            if z_1 == 0:
                print('no rotamer could be placed for site 1')
                continue
            boltzman_weights_norm1 = boltz1 / z_1

            boltz2 = self.lj_calculation(rotamersSite2, protein_structure, residues[1], fret=True)
            boltz2 = np.multiply(lib2_norm, boltz2)
            z_2 = np.nansum(boltz2)
            # print z_2
            if z_2 == 0:
                print('no rotamer could be placed for site 2')
                continue
            boltzman_weights_norm2 = boltz2 / z_2

            # define the atoms to measure the distances between
            frame_distributions = np.zeros((rax.size))
            rotamer1oxigen = rotamersSite1.select_atoms("name O6")
            rotamer2oxigen = rotamersSite2.select_atoms("name O6")

            size = len(rotamersSite2.trajectory)
            oxi1_pos = np.array([rotamer1oxigen.positions for x in rotamersSite1.trajectory])
            oxi2_pos = np.array([i for x in rotamersSite2.trajectory for i in rotamer2oxigen.positions])
            if self.k2_calc_dyn:
                # self.kappa2 = np.full(len(oxi2_pos), self.kappa2)  # average k2 calculated previously
                pass
            else:
                self.kappa2 = np.full(len(oxi2_pos), 2 / 3)
            if self.k2_calc_static:
                if self.libname_1.split()[1] == '488':
                    rotamer1_c11 = rotamersSite1.select_atoms("name C11")
                    rotamer2_c11 = rotamersSite2.select_atoms("name C17")
                    rotamer1_c12 = rotamersSite1.select_atoms("name C12")
                    rotamer2_c12 = rotamersSite2.select_atoms("name C18")
                else:
                    rotamer1_c11 = rotamersSite1.select_atoms("name C17")
                    rotamer2_c11 = rotamersSite2.select_atoms("name C11")
                    rotamer1_c12 = rotamersSite1.select_atoms("name C18")
                    rotamer2_c12 = rotamersSite2.select_atoms("name C12")

                c11_probe1_pos = np.array([i for x in rotamersSite1.trajectory for i in rotamer1_c11.positions])
                c11_probe2_pos = np.array([i for x in rotamersSite2.trajectory for i in rotamer2_c11.positions])

                c12_probe1_pos = np.array([i for x in rotamersSite1.trajectory for i in rotamer1_c12.positions])
                c12_probe2_pos = np.array([i for x in rotamersSite2.trajectory for i in rotamer2_c12.positions])

                inner_probe1_vect = c12_probe1_pos - c11_probe1_pos
                inner_probe1_vect /= np.linalg.norm(inner_probe1_vect, axis=1)[:, np.newaxis]

                inner_probe2_vect = c11_probe2_pos - c12_probe2_pos
                inner_probe2_vect /= np.linalg.norm(inner_probe2_vect, axis=1)[:, np.newaxis]

            dists_array = np.zeros((1, size), dtype=np.float64)
            p_of_k = []

            for position1_index, position1 in enumerate(oxi1_pos):
                mda_dist.distance_array(position1, oxi2_pos, result=dists_array, backend="OpenMP")
                if self.k2_calc_static:
                    self.kappa2 = np.full(len(oxi2_pos), 2 / 3)
                    mutual_displacement = np.dot(inner_probe2_vect, inner_probe1_vect[position1_index].T)

                    outer_vect1 = c12_probe1_pos[position1_index] - c11_probe2_pos
                    outer_vect1 /= np.linalg.norm(outer_vect1, axis=1)[:, np.newaxis]
                    outer_vect2 = c11_probe1_pos[position1_index] - c12_probe2_pos
                    outer_vect2 /= np.linalg.norm(outer_vect2, axis=1)[:, np.newaxis]

                    trans_dip_moment1_probe1 = np.dot(outer_vect1, inner_probe1_vect[position1_index].T)
                    trans_dip_moment2_probe1 = np.dot(outer_vect2, inner_probe1_vect[position1_index].T)

                    for pos2_index, element in enumerate(dists_array[0]):  # could we take out this for loop?
                        trans_dip_moment1_probe2 = np.dot(outer_vect1[pos2_index], inner_probe2_vect[pos2_index])
                        trans_dip_moment2_probe2 = np.dot(outer_vect2[pos2_index], inner_probe2_vect[pos2_index])

                        new_k2 = np.power((mutual_displacement[pos2_index] -
                                           ((3. / 4.) * (trans_dip_moment1_probe1[pos2_index] +
                                                         trans_dip_moment2_probe1[pos2_index]) * (
                                                trans_dip_moment1_probe2 +
                                                trans_dip_moment2_probe2))), 2)
                        self.kappa2[pos2_index] = new_k2
                        kappa2_hold.append(new_k2)
                        weights_list.append((boltzman_weights_norm1[position1_index] *
                                                     boltzman_weights_norm2[pos2_index]))
                dists_array /= 10.0
                #print 'debugger dists_array'
                #print dists_array
                efficiency_array = self.fret_efficiency(dists_array, k2=self.kappa2, r0=self.r0)
                #print 'efficiency_array'
                #print efficiency_array

                dists_array = np.round(100 * efficiency_array)-1  # correct?
                #print 'debugger dists_array'
                #if np.isnan(dists_array).any():
                #    print 'nan found in dists_array'


                for pos2_index, element in enumerate(dists_array[0]):  # could we take out this for loop?
                    element = int(element)
                    frame_distributions[element] += (boltzman_weights_norm1[position1_index] *
                                                     boltzman_weights_norm2[pos2_index])
            # print 'b weights debug 1, then 2'
            # print boltzman_weights_norm1
            # print boltzman_weights_norm2
            distributions[0] += frame_distributions
        if self.k2_calc_static:
            logger.info('Static K2 calculated: {:.3f}'.format(np.average(kappa2_hold, weights=weights_list)))
            # for i in np.linspace(0, 4, num=40):
            #     if i < 1:
            #         p_of_k.append((1 / (2 * np.sqrt(3 * i))) * np.log(2 + np.sqrt(3)))
            #     else:
            #         p_of_k.append((1 / (2 * np.sqrt(3 * i))) * np.log(((2 + np.sqrt(3)) / (np.sqrt(i) + np.sqrt(i - 1)))))
            # print np.average(kappa2_hold, weights=weights_list)
            # sns.set(style="ticks")
            # # plt.hist(kappa2_hold, bins=40, normed=True)
            # first = np.max(np.histogram(kappa2_hold, normed=True, weights=weights_list, bins=40)[0])
            # print first
            # plt.hist(kappa2_hold, bins=40, weights=weights_list, normed=True)
            # plt.plot(np.linspace(0, 4, num=40), p_of_k, color='red')
            # plt.axvline(np.average(kappa2_hold, weights=weights_list), color='green')
            # plt.text(x=np.average(kappa2_hold, weights=weights_list)+0.1, y=0.9*first,
            #          s=r'$\left\langle\kappa^{{2}}\right\rangle$ = {0:.3f}'.format(np.average(kappa2_hold,
            #                                                                               weights=weights_list)))
            # plt.title('Cutoff: {}'.format(self.output_file.split('_')[-1][:-15]))
            # plt.xlabel('$\kappa^{2}$')
            # plt.ylabel('$p(\kappa^{2})$')
            # sns.despine()
            # plt.savefig('p_of_k_{0}cutoff.png'.format(self.output_file.split('_')[-1][:-15]), dpi=600)
            # plt.close()
            self.kappa2 = np.average(kappa2_hold, weights=weights_list)


        # print 'debugger kappa2'
        # print self.kappa2
        #
        # print 'debugger rax'
        # print rax
        #
        # print 'debugger distr'
        # print distributions
        try:
            self.efficiency = np.average(rax, weights=~np.isnan(distributions[0]))
            print('hello')
            logger.info('FRET efficiency: {:.3f}'.format(self.efficiency))
        except ZeroDivisionError:
            logger.error('Unable to place rotamers, try a bigger rotamer library.')
            return

        if self.replicas == 1:
            self.plot_and_save_fret(rax, residues, distributions)
            # distributions = distributions[0] / (math.ceil((self.stop_frame-self.start_frame)/self.jump_frame))
            # plt.plot(rax, distributions)
            # plt.xlabel("FRET efficiency")
            # plt.ylabel("Probability density")
            # plt.savefig(self.output_plot, dpi=300)
            # plt.close()
            #
            # with open(self.output_file, 'w') as OUTPUT:
            #     for index, value in enumerate(rax):
            #         OUTPUT.write('{0:>7.4f} {1:>7.4f}\n'.format(value, distributions[index]))
            # logger.info("Distance distribution for residues {0[0]} - {0[1]} "
            #             "was written to {1}".format(residues, self.output_file))
            # logger.info("Distance distribution for residues {0[0]} - {0[1]} "
            #             "was plotted to {1}".format(residues, self.output_plot))

        # else:
        #     inv_distr = np.fft.ifft(distributions[0]) * np.fft.ifft(vari)
        #     distributions_global = np.real(np.fft.fft(inv_distr))
        #     distributions_global = distributions_global / np.sum(distributions_global)
        #     plt.plot(rax[100:401], distributions_global[200:501])
        #     plt.xlim([1, 5])
        #     plt.ylim([-np.max(distributions_global[200:501]) / 20,
        #               np.max(distributions_global[200:501]) + np.max(distributions_global[200:501]) / 20])
        #     plt.xlabel(r"Spin-label distance $d$ ($\AA$)")
        #     plt.ylabel("Probability density")
        #     plt.savefig(self.output_plot, dpi=300)
        #     plt.close()
        #
        #     with open(self.output_file, 'w') as OUTPUT:
        #         for index, distance in enumerate(rax[100:401]):
        #             OUTPUT.write('{0:>7.4f} {1:>7.4f}\n'.format(distance, distributions_global[index + 200]))
        #     logger.info("Distance distribution for residues {0[0]} - {0[1]} "
        #                 "was written to {1}".format(residues, self.output_file))
        #     logger.info("Distance distribution for residues {0[0]} - {0[1]} "
        #                 "was plotted to {1}".format(residues, self.output_plot))
        #
        #     for replica in range(1, self.replicas + 1):
        #         inv_distr = np.fft.ifft(distributions[replica]) * np.fft.ifft(vari)
        #         distributions_replica = np.real(np.fft.fft(inv_distr))
        #         distributions_replica = distributions_replica / np.sum(distributions_replica)
        #         plt.plot(rax[100:401], distributions_replica[200:501])
        #         plt.xlim([1, 5])
        #         plt.ylim([-np.max(distributions_replica[200:501]) / 20,
        #                   np.max(distributions_replica[200:501]) + np.max(distributions_replica[200:501]) / 20])
        #         plt.xlabel(r"Spin-label distance $d$ ($\AA$)")
        #         plt.ylabel("Probability density")
        #         if self.chains:
        #             ext = ".png"
        #             self.output_plot = "{0}-{1[0]}{2[0]}-{1[1]}{2[1]}_replica{3}{4}".format(output_file,
        #                                                                                     residues,
        #                                                                                     self.chains,
        #                                                                                     replica,
        #                                                                                     ext)
        #             ext = '.dat'
        #             self.output_file = "{0}-{1[0]}{2[0]}-{1[1]}{2[1]}_replica{3}{4}".format(output_file,
        #                                                                                     residues,
        #                                                                                     self.chains,
        #                                                                                     replica,
        #                                                                                     ext)
        #         else:
        #             ext = ".png"
        #             self.output_plot = "{0}-{1[0]}-{1[1]}_replica{3}{4}".format(output_file,
        #                                                                         residues,
        #                                                                         replica,
        #                                                                         ext)
        #             ext = '.dat'
        #             self.output_file = "{0}-{1[0]}-{1[1]}_replica{3}{4}".format(output_file,
        #                                                                         residues,
        #                                                                         replica,
        #                                                                         ext)
        #         plt.savefig(self.output_plot, dpi=300)
        #         plt.close()
        #         with open(self.output_file, 'w') as OUTPUT:
        #             for index, distance in enumerate(rax[100:401]):
        #                 OUTPUT.write('{0:>7.4f} {1:>7.4f}\n'.format(distance, distributions_replica[index + 200]))
        #
        #         logger.info("Distance distribution for residues {0[0]} - {0[1]} and replica {1} "
        #                     "was written to {2}".format(residues, replica, self.output_file))
        #         logger.info("Distance distribution for residues {0[0]} - {0[1]} and replica {1} "
        #                     "was plotted to {2}".format(residues, replica, self.output_plot))
    def results(self):
        if self.k2_calc_static is False and self.k2_calc_dyn is False:
            self.kappa2 = self.kappa2[0]
        return self.libname_1, self.libname_2, self.efficiency, self.kappa2
