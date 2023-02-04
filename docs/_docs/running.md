# Running Calculations

## FRETpredict

The class for running the PRE calculations is `FRETpredict`.

Here is how to calculate FRET efficiencies from a single protein structure using the rotamer-library approach

~~~ python
import MDAnalysis
from FRETpredict import FRETpredict

# Create a MDAnalysis.Universe object for the protein structure
u = MDAnalysis.Universe('tests/test_systems/Hsp90/openHsp90.pdb')

# Create instance of the FRETpredict class
FRET = FRETpredict(protein=u, residues=[452, 637], chains=['A', 'B'], temperature=293, 
                   fixed_R0=True, r0=6.3, electrostatic=True,
                   libname_1='AlexaFluor 594 C1R cutoff30',
                   libname_2='AlexaFluor 568 C1R cutoff30', 
                   output_prefix='E30_594_568')

# Run FRET efficiency calculations
FRET.run()
~~~

The program generates the pickle file named `E30_594_568-data-452-637.pkl` which can be loaded as a pandas DataFrame and contains average, standard deviation, and standard error of the orientation factor (`k2`) and the FRET efficiencies calculated in the Static (`Estatic`), Dynamic (`Edynamic1`), and Dynamic+ (`Edynamic2`) regimes.

### Reweighting

Here is how to calculate FRET efficiencies from a protein trajectory using the rotamer-library approach. In this example, we reweight the trajectory by per-frame weights, which are passed as a numpy array to the `weights` attribute

~~~ python
import MDAnalysis
from FRETpredict import FRETpredict

# Create a MDAnalysis.Universe object for the protein trajectory.
u_0M = MDAnalysis.Universe('actr.gro', 'actr_urea0.xtc')

# Load per-frame weights (from BME reweighting or enhanced sampling technique)
weights = np.load('weights.npy')

# Create instance of the FRETpredict class
FRET = FRETpredict(protein=u_0M, residues=[3, 61],
                    fixed_R0=True, r0=5.40,
                    electrostatic=True,
                    libname_1='AlexaFluor 488 C1R cutoff20',
                    libname_2='AlexaFluor 594 C1R cutoff20',
                    weights=weights)

# Run FRET efficiency calculations
FRET.run()
~~~

