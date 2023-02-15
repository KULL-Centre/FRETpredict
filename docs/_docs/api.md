# FRETpredict API Reference

FRETpredict has four main classes. Two of them are dedicated to FRET efficiency calculations and the others to Rotamer Library generation.

## FRET efficiency calculation

### `FRETpredict` 
Performs predictions of FRET efficiencies. It can be found in the `FRET.py` file. <br>

Methods:<br> 
- `compute_chromophore_vectors`: calculates chromophores transition dipole moment vector and center-to-center distance distribution for k2
            calculations, 
- `calculateR0`: calculates FRET R0 between the donor and acceptor chromophores, 
- `trajectoryAnalysis`: calculates <E> distribution (i.e. one <E> for each protein trajectory frame) in static, dynamic1, and dynamic2 regimes, 
- `save`: calculates k2 distribution and k2, Static, Dynamic1, Dynamic2 averaging, with or without reweighting, and save data to file, 
- `reweight`: Alias for reweigthing calculations, calls save() function with weights for reweighting,
- `run`: runs FRET efficiency calculations by calling trajectoryAnalysis() and saving data to file.

### `Operations` 
It is the base class containing attributes and methods inherited and used by the calculation classes. It can be found in the `utils.py` file.<br> 

Methods:<br>
- `precalculate_rotamer`: selects placement residue atoms, compute Lennard-Jones and Electrostatic (Debye-Huckel) parameters for chromophore rotamers and protein.
   
The Lennard-Jones parameters of the CHARMM36 force field in the `lennardjones.py` file are used to calculate the external energy contribution to the Boltzmann weight of each conformer.

~~~python 
FRETpredict.lennardjones.lj_parameters = {
    'C': {
        'vdw': 1.70,
        'p_q': 0,
        'p_Rmin2': 2.275,
        'eps': -0.020
    }, 
    ...
}
~~~
   
- `rotamer_placement`: places chromophore rotamers on a selected protein residue by rotation and translation, return chromophore Universe object with new coordinates,
- `lj_calculation`: calculates Boltzmann weights for each rotamer pose by summing Electrostatic and potential energy contributions,
- `rotamerWeights`: calculates Boltzmann distribution for each rotamer library,
- `weightedAvgSDSE`: calculates the weighted average and standard deviation,
- `calculate_ws`: calculates per-frame weights for reweighting,
- `fraction_frames`: computes effective fraction of frames contributing to the averages.
   
 __[`Tutorial_FRETpredict_pp11`](https://github.com/Monte95/FRETpredict/blob/62ee39e82e82691a237da8e927d686378aff5fb1/tests/tutorials/Tutorial_FRETpredict_pp11.ipynb)__ is a Jupyter Notebook with simple tutorials on how to compute the FRET efficiency and use the FRETpredict functions on the trajectory of a Poliproline 11 (pp11) system.

## Rotamer Library generation
   
### `RotamerClusters`
Creates a rotamer library starting from a dye + linker (FRET probe) trajectory. It can be found in the `rotamer_libraries.py` file.<br>

Methods:<br>
- `calcDihe`: calculates dihedral angles on the dye+linker trajectory,
- `genPeaks`: computes linker dihedral peaks,
- `genClusters`:<br>
   1. generate combinations of dihedral angles from the peaks (cluster centers C1).
   2. K-means Clustering
      - Assign each trajectory frame to the cluster center C1 of least square deviation.
      - Calculate the average over the dihedral angles that were assigned to the same cluster center. This results in a set of new centers (C2).
      - Assign each trajectory frame to the cluster center C2 of least square deviation
   3. Find the trajectory frame that best represents the cluster center.
- `filterCluster`:<br>
   1. Filter the cluster centers C2 based on a cutoff on N (cluster population). This results in a different number of centers (C3).
   2. Reassign the discarded frames to the remaining C3 cluster center of least square deviation.
- `genRotLib`: translates + rotates C3 cluster centers conformations, and write data to file, 
- `plotClustHist`: plots dihedral distribution and peaks with cluster centers C3 dihedrals,
- `plotClustPolar`: plots dihedral distribution and peaks with cluster centers C3 dihedrals on a polar plot,
- `run`: runs all the calculations to generate a rotamer library from a dye+linker trajectory.

### `RotamerLibrary`

Makes available the rotamer library attributes. Rotamers libraries consist of three different files:
 - a structure file with a protein conformation (PDB format), 
 - a trajectory file with all the rotamer conformations (DCD format),
 - a text file for the populations (i.e., weight) of each chromophore conformation. These files are included in the `FRETpredict/lib` folder.
   
 Specific fluorophore parameters are specified in `FRETpredict/lib/libraries.yml`.
   
  ~~~python
   'AlexaFluor 488 C1R': {'author': 'D Montepietra, G Tesei, JM Martins, MBA Kunze, RB Best, K Lindorff-Larsen',
                          'citation': 'https://doi.org/10.1101/2023.01.27.525885',
                          'filename': 'A48_C1R_cutoff30',
                          'licence': 'GPLv3',
                          'mu': ['C2', 'C13 and resname A48'],
                          'negative': ['S1 and resname A48',
                          'S2 and resname A48',
                          'C20 and resname A48'],
                          'positive': ['N1 and resname A48', 'N2 and resname A48'],
                          'r': ['C7 and resname A48']},
  ~~~

__[`Tutorial_generate_new_rotamer_libraries.ipynb`](https://github.com/Monte95/FRETpredict/blob/eef8bf0d219109ada605e943ecc4b1aa9dde86df/tests/tutorials/Tutorial_generate_new_rotamer_libraries.ipynb)__ is Jupyter Notebook tutorial on how to create and add new rotamer libraries.


