# -*- coding: utf-8 -*-
"""

Parameters for calculation of VDW interaction between in-range atoms

"""

#: CHARMM36 FF parameters for Lennard-Jones calculation
lj_parameters = {
    'C': {
        'vdw': 1.70,
        'p_q': 0,
        'p_Rmin2': 2.275,
        'eps': -0.020
    },
    'N': {
        'vdw': 1.55,
        'p_q': 0,
        'p_Rmin2': 1.850,
        'eps': -0.2000
    # 'N': {
    #     'vdw': 2.0,
    #     'p_q': 0,
    #     'p_Rmin2': 2.950,
    #     'eps': -0.2000
    },
    'O': {
        'vdw': 1.52,
        'p_q': 0,
        'p_Rmin2': 1.700,
        'eps': -0.120
    },
    'S': {
        'vdw': 1.80,
        'p_q': 0,
        'p_Rmin2': 2.000,
        'eps': -0.450
    },
    'H': {
        'vdw': 0.500,
        'p_q': 0,
        'p_Rmin2': 0.5,
        'eps': -0.030
    }
}

#Create simpler dictionaries for optimization:
vdw = {atom:lj_parameters[atom]['vdw'] for atom in lj_parameters}
p_Rmin2 = {atom:lj_parameters[atom]['p_Rmin2'] for atom in lj_parameters}
eps = {atom:lj_parameters[atom]['eps'] for atom in lj_parameters}

