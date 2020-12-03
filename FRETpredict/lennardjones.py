# -*- coding: utf-8 -*-
"""

Parameters for calculation of VDW interaction between in-range atoms

"""

#: CHARMM36 FF parameters for Lennard-Jones calculation
lj_parameters = {
    'C': {
        'p_Rmin2': 2.02446316,
        'eps': -0.06394724
    },
    'N': {
        'p_Rmin2': 1.89285714,
        'eps': -0.15428571
    },
    'O': {
        'p_Rmin2': 1.693,
        'eps': -0.12642017
    },
    'S': {
        'p_Rmin2': 2.1,
        'eps': -0.47
    },
    'H': {
        'p_Rmin2': 0.98357778,
        'eps': -0.03466645
    }
}

