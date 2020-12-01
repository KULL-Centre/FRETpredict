# -*- coding: utf-8 -*-
"""

Parameters for calculation of VDW interaction between in-range atoms

"""

#: CHARMM36 FF parameters for Lennard-Jones calculation
lj_parameters = {
    'C': {
        'p_Rmin2': 2.0245,
        'eps': 0.0639
    },
    'N': {
        'p_Rmin2': 1.8929,
        'eps': 0.1543
    },
    'O': {
        'p_Rmin2': 1.6930,
        'eps': 0.1264
    },
    'S': {
        'p_Rmin2': 2.1000,
        'eps': 0.4700
    },
    'H': {
        'p_Rmin2': 0.9836,
        'eps': 0.0347
    }
}

#Create simpler dictionaries for optimization:
p_Rmin2 = {atom:lj_parameters[atom]['p_Rmin2'] for atom in lj_parameters}
eps = {atom:lj_parameters[atom]['eps'] for atom in lj_parameters}

