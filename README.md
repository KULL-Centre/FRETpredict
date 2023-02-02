FRETpredict
===========

Overview
--------

A package for FRET Efficiency prediction of protein structures and trajectories, based on the Rotamer Library Approach (RLA).

Installation
------------

To install FRETpredict, use the [PyPI package](https://pypi.org/project/FRETpredict):

```bash
  pip install FRETpredict
```

or clone the repo and install locally:

```bash
  git clone https://github.com/KULL-Centre/FRETpredict.git
  cd FRETpredict

  pip install -e . 
```

The software requires Python 3.6+.

Testing
-------

```bash
  pip install pytest

  python -m pytest
```

Code Example
------------

```python

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

```
To add reweighting calculations

```python

FRET.save(reweight=True, reweight_prefix='E_pp11_reweighted')

```


Tutorials
---------

- __[Tutorial_FRETpredict_Hsp90](https://github.com/Monte95/FRETpredict/blob/eef8bf0d219109ada605e943ecc4b1aa9dde86df/tests/tutorials/Tutorial_FRETpredict_Hsp90.ipynb)__ : Jupyter Notebook with simple tutorials on how to use the code on the Hsp90 system.

- __[Generate new rotamer libraries](https://github.com/Monte95/FRETpredict/blob/eef8bf0d219109ada605e943ecc4b1aa9dde86df/tests/tutorials/Tutorial_generate_new_rotamer_libraries.ipynb)__ : Jupyter Notebook on how to create and add new rotamer libraries.


Structure
---------
```
FRETpredict/
├─ FRETpredict/
│  ├─ lib/
│  │  ├─ R0/
│  ├─ FRET.py
│  ├─ lennardjones.py
│  ├─ libraries.py
│  ├─ R0_calculation.py
│  ├─ rotamer_libraries.py
│  ├─ utils.py
│  ├─ __init__.py
├─ tests/
|  ├─ test_Hsp90.py
│  ├─ test_systems/
│  │  ├─ Hsp90/
│  ├─ tutorials/
│  │  ├─ genLIB/
│  │  ├─ test/
│  │  ├─ Tutorial_FRETpredict_Hsp90.ipynb
│  │  ├─ Tutorial_generate_new_rotamer_libraries.ipynb
├─ LICENSE
├─ README.md
├─ setup.py
```

Contributors
-------------

[Daniele Montepietra (@Monte95)](https://github.com/Monte95)

[Giulio Tesei (@gitesei)](https://github.com/gitesei)

[João M Martins (@joaommartins)](https://github.com/joaommartins)

[Micha BA Kunze (@mbakunze)](https://github.com/mbakunze)

[Robert Best (@bestlab)](https://github.com/bestlab)

[Kresten Lindorff-Larsen (@lindorff-larsen)](https://github.com/lindorff-larsen)
