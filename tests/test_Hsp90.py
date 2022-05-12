import pytest
import MDAnalysis
import numpy as np
from FRETpredict import FRETpredict

# Test with single frame, explicit R0 calculation
def test_Hsp90_R0_calculation():
    
    # Import protein structure
    u = MDAnalysis.Universe('../test_systems/Hsp90/openHsp90.pdb')
    
    # Instantiate class object
    FRET = FRETpredict(protein=u, residues=[452, 637], temperature=293, 
                   chains=['A', 'B'], 
                   donor='AlexaFluor 594', acceptor='AlexaFluor 568', 
                   electrostatic=True, r0lib='../lib/R0/',
                   libname_1=f'AlexaFluor 594 C1R cutoff30',
                   libname_2=f'AlexaFluor 568 C1R cutoff30', 
                   output_prefix=f'data/E30')
    
    # Run FRETpredict calculations
    assert FRET.run(), "Cannot calculate FRET Efficiency with explicit R0 calculation"
    
    # Read computed data from file
    assert pd.read_pickle(r'data/E30-data-594-568.pkl'), 
    "Could not read FRET data from file (explicit R0 calculations)"

# Test with single frame, fixed R0 value
def test_Hsp90_fixed_R0():
    
    # Import protein structure
    u = MDAnalysis.Universe('../test_systems/Hsp90/openHsp90.pdb')
    
    # Instantiate class object
    FRET_fixedR0 = FRETpredict(protein=u, residues=[452, 637], temperature=293, 
                           chains=['A', 'B'], 
                           fixed_R0=True, r0=5.5, electrostatic=True,
                           libname_1=f'AlexaFluor 594 C1R cutoff30',
                           libname_2=f'AlexaFluor 568 C1R cutoff30', 
                           output_prefix=f'data/E30_fixedR0')
    
    # Run FRETpredict calculations
    assert FRET.run(), "Cannot calculate FRET Efficiency with fixed R0"
    
    # Read computed data from file
    assert pd.read_pickle(r'data/E30_fixedR0-data-594-568.pkl'), 
    "Could not read FRET data from file (fixed R0)"

# Test with trajectory, explicit R0 calculations
def test_Hsp90_trajectory_R0_calculation():
    
    # Import protein trajectory
    u = MDAnalysis.Universe('../test_systems/Hsp90/conf.pdb', 
                            '../test_systems/Hsp90/Hsp90_open_all.xtc')
    
    # Instantiate class object
    FRET_traj = FRETpredict(protein=u, residues=[452, 637], temperature=293, 
                        chains=['A', 'B'], 
                        donor='AlexaFluor 594', acceptor='AlexaFluor 568', 
                        electrostatic=True, r0lib='../lib/R0/',
                        libname_1=f'AlexaFluor 594 C1R cutoff30',
                        libname_2=f'AlexaFluor 568 C1R cutoff30', 
                        output_prefix=f'data/E30_traj')
    
    # Run FRETpredict calculations
    assert FRET.run(), "Cannot calculate FRET Efficiency for trajectory with explicit R0 calculation"
    
    # Read computed data from file
    assert pd.read_pickle(r'data/E30_traj-data-594-568.pkl'), 
    "Could not read FRET data from file (trajectory with explicit R0 calculations)"

# Test with trajectory, fixed R0 value
def test_Hsp90_trajectory_fixed_R0():
    
    # Import protein trajectory
    u = MDAnalysis.Universe('../test_systems/Hsp90/conf.pdb', 
                            '../test_systems/Hsp90/Hsp90_open_all.xtc')
    
    # Instantiate class object
    FRET_traj = FRETpredict(protein=u, residues=[452, 637], temperature=293, 
                        chains=['A', 'B'], 
                        fixed_R0=True, r0=5.5, electrostatic=True,
                        libname_1=f'AlexaFluor 594 C1R cutoff30',
                        libname_2=f'AlexaFluor 568 C1R cutoff30', 
                        output_prefix=f'data/E30_traj_fixedR0')
    
    # Run FRETpredict calculations
    assert FRET.run(), "Cannot calculate FRET Efficiency for trajectory with fixed R0"
    
    # Read computed data from file
    assert pd.read_pickle(r'data/E30_traj_fixedR0-data-594-568.pkl'), 
    "Could not read FRET data from file (trajectory with fixed R0)"
    
    
    
 