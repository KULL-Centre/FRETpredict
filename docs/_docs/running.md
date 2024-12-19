# Running FRETpredict calculations

The class for running the FRET efficiency calculations is `FRETpredict`.

Here is how to calculate FRET efficiencies from protein trajectory using the Rotamer Library Approach (RLA).

~~~ python
import MDAnalysis
from FRETpredict import FRETpredict

# Create a MDAnalysis.Universe object for the protein structure.
u = MDAnalysis.Universe('tests/test_systems/pp11/pp11.pdb', 'tests/test_systems/pp11/pp11.xtc')

# Create instance of the FRETpredict class
FRET = FRETpredict(protein=u, residues=[0, 12], chains=['A', 'A'], temperature=298, 
                   donor='AlexaFluor 488', acceptor='AlexaFluor 594', electrostatic=True,
                   libname_1=f'AlexaFluor 488 C1R cutoff10',
                   libname_2=f'AlexaFluor 594 C1R cutoff10',  
                   output_prefix='E_pp11')

# Run FRET efficiency calculations.
FRET.run()
~~~

The program generates the pickle file named `E_pp11-data-0-12.pkl` which can be loaded as a pandas DataFrame and contains average, standard deviation, and standard error of the orientation factor (`k2`) and the FRET efficiencies calculated in the Static (`Estatic`), Dynamic (`Edynamic1`), and Dynamic+ (`Edynamic2`) regimes.

### Reweighting

Here is how to reweight the trajectory by per-frame dye-protein weights. 
~~~ python
FRET.reweight(boltzmann_weights=True, reweight_output_prefix='E_pp11_reweighted')
~~~
Per-frame dye-protein weights can be combined with user-provided weights (numpy array) from previous calculations (e.g., Enhanced sampling simulations), which are passed as a numpy array to the `user_weights` attribute.
~~~ python
FRET.reweight(boltzmann_weights=True, user_weights=user_weights_pp11, reweight_output_prefix='E_pp11_reweighted')
~~~
User-provided weights alone can be used to reweight the ensemble. In this case Boltzmann weights are not used and only frames with steric clashes between dye and protein are discarded (default behaviour).
~~~ python
FRET.reweight(boltzmann_weights=False, user_weights=user_weights_pp11, reweight_output_prefix='E_pp11_reweighted')
~~~

Each of these lines generates a pickle file named `E_pp11_reweighted-data-0-12.pkl` with all the reweighted parameters, which can be loaded as a pandas DataFrame.

In this other example, we can directly obtain the reweighted FRET efficiencies combining the dye-protein weights with custom pre-computed weights, which are passed as a numpy array to the `user_weights` attribute in the `FRETpredict` class.

~~~ python
user_weights_pp11 = np.load('user_weights_pp11.npy')

FRET = FRETpredict(protein=u, residues=[0, 12], chains=['A', 'A'], temperature=298, 
                   donor='AlexaFluor 488', acceptor='AlexaFluor 594', electrostatic=True,
                   libname_1=f'AlexaFluor 488 C1R cutoff10',
                   libname_2=f'AlexaFluor 594 C1R cutoff10',  
                   boltzmann_weights=True,
                   user_weights=user_weights_pp11,
                   output_prefix='E_pp11_reweighted')
~~~
