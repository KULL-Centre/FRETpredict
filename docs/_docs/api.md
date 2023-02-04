# FRETpredict API Reference

FRETpredict has two main classes:
- `FRETpredict` performs predictions of FRET efficiencies. <br>
   Functions: `compute_chromophore_vectors`, `calculateR0`, `trajectoryAnalysis`, `run`, and `save`.
- `Operations` is the base class containing attributes and methods inherited and used by the calculation classes. <br> 
   Functions: `precalculate_rotamer`, `rotamer_placement`, `lj_calculation`, `rotamerWeights`, and `weightedAvgSDSE`.

## RotamerLibrary class

`FRETpredict.libraries.LIBRARIES`: Rotamers libraries consist of a PDB file, a DCD files and a text file for the weights. These files are included in the `FRETpredict/lib` folder.

## Lennard-Jones parameters

`FRETpredict.lennardjones`: Lennard-Jones parameters of the CHARMM36 force field used to calculate the external 
energy contribution to the Boltzmann weight of each conformer.

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
