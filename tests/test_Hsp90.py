import pytest
import MDAnalysis
import pandas as pd
import numpy as np
from FRETpredict import FRETpredict

# Test with single frame, explicit R0 calculation
def test_Hsp90_R0_calculation():
    
    # Import protein structure
    u = MDAnalysis.Universe('tests/test_systems/Hsp90/openHsp90.pdb')
    
    # Instantiate class object
    FRET = FRETpredict(protein=u, residues=[452, 637], temperature=293, 
                   chains=['A', 'B'], 
                   donor='AlexaFluor 594', acceptor='AlexaFluor 568', 
                   electrostatic=True, r0lib='../FRETpredict/lib/R0/',
                   libname_1='AlexaFluor 594 C1R cutoff30',
                   libname_2='AlexaFluor 568 C1R cutoff30', 
                   output_prefix='data/E30')
    
    # Run FRETpredict calculations
    FRET.run()
#     "Cannot calculate FRET Efficiency with explicit R0 calculation"
    
    # Read computed data from file
    assert pd.read_pickle('data/E30-data-452-637.pkl')['Estatic'] != np.nan,  \
    "Could not read FRET data from file (explicit R0 calculations)"

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
                           output_prefix='data/E30_fixedR0')
    
    # Run FRETpredict calculations
    FRET_fixedR0.run()
#     "Cannot calculate FRET Efficiency with fixed R0"
    
    # Read computed data from file
    assert pd.read_pickle('data/E30_fixedR0-data-452-637.pkl')['Estatic'] != np.nan,  \
    "Could not read FRET data from file (fixed R0)"

    
    
    
 
