import pytest
import MDAnalysis
import pandas as pd
import numpy as np
import os
from FRETpredict import FRETpredict

if not os.path.isdir('tests/test_systems/Hsp90/output'):
    os.mkdir('tests/test_systems/Hsp90/output')

# Test with single frame, explicit R0 calculation
def test_Hsp90_R0_calculation():

    # Import protein structure
    u = MDAnalysis.Universe('tests/test_systems/Hsp90/openHsp90.pdb')

    # Instantiate class object
    FRET = FRETpredict(protein=u, residues=[452, 637], temperature=293,
                   chains=['A', 'B'],
                   donor='AlexaFluor 594', acceptor='AlexaFluor 568',
                   electrostatic=True,
                   libname_1='AlexaFluor 594 C1R cutoff30',
                   libname_2='AlexaFluor 568 C1R cutoff30',
                   output_prefix='tests/test_systems/Hsp90/output/E30',
                   verbose=True)

    eff = []

    # Run FRETpredict calculations
    FRET.run()
    eff.append(pd.read_pickle('tests/test_systems/Hsp90/output/E30-data-452-637.pkl')['Estatic'])

    # Reweight FRETpredict calculations
    FRET.reweight(user_weights=np.asarray([0.5]))
    eff.append(pd.read_pickle('tests/test_systems/Hsp90/output/E30-data-452-637.pkl')['Estatic'])

    # Read computed data from file
    print(eff)
    assert np.allclose(eff,0.46766,rtol=0,atol=1E-05)

# Test with single frame, fixed R0 value
def test_Hsp90_fixed_R0():

    # Import protein structure
    u = MDAnalysis.Universe('tests/test_systems/Hsp90/openHsp90.pdb')

    # Instantiate class object
    FRET_fixedR0 = FRETpredict(protein=u, residues=[452, 637], temperature=293,
                           chains=['A', 'B'],
                           fixed_R0=True, r0=5.5, electrostatic=True,
                           libname_1='AlexaFluor 594 C1R cutoff30',
                           libname_2='AlexaFluor 568 C1R cutoff30',
                           output_prefix='tests/test_systems/Hsp90/output/E30_fixedR0',
                           verbose=True)

    eff = []

    # Run FRETpredict calculations
    FRET_fixedR0.run()
    eff.append(pd.read_pickle('tests/test_systems/Hsp90/output/E30_fixedR0-data-452-637.pkl')['Estatic'])

    # Reweight FRETpredict calculations
    FRET_fixedR0.reweight(user_weights=np.asarray([0.5]))
    eff.append(pd.read_pickle('tests/test_systems/Hsp90/output/E30_fixedR0-data-452-637.pkl')['Estatic'])

    # Read computed data from file
    print(eff)
    assert np.allclose(eff,0.53814,rtol=0,atol=1E-05)
