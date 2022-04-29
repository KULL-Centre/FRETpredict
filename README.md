FRETpredict
===========

Overview
--------

A package for FRET Efficiency prediction from protein structures and trajectories, based on the Rotamer Library Approach (RLA).

Installation
------------

To install FRETpredict, use the [PyPI package](https://pypi.org/project/FRETpredict):

```bash
  pip install FRETpredict
```

or clone the repo:

```bash
  git clone https://github.com/KULL-Centre/FRETpredict.git
  cd FRETpredict

  pip install -e . 
```

The software requires Python 3.6+.

Code Example
------------

```python

import MDAnalysis
from FRET import FRETpredict

# Create a MDAnalysis.Universe object for the protein structure.
u = MDAnalysis.Universe('test_systems/Hsp90/openHsp90.pdb')

# Create instance of the FRETpredict class
FRET = FRETpredict(protein=u, residues=[452, 637], chains=['A', 'B'], temperature=293, 
                   donor=594, acceptor=568, electrostatic=True,
                   libname_1='Alexa 594 cutoff10',
                   libname_2='Alexa 568 cutoff10', 
                   output_prefix='test/E10_594_568')

# Run FRET efficiency calculations.
FRET.run()

```

Tutorial
--------

Jupyter Notebook with simple tutorials on how to use the code on the Hsp90 system: https://github.com/Monte95/FRETpredict/blob/e3311a7af03c045000e4ac69928307d7aca0d684/FRETpredict/tutorial_FRETpredict.ipynb

Documentation
-------------

Testing
-------

```bash
  git clone https://github.com/KULL-Centre/FRETpredict.git
  cd FRETpredict

  python -m pytest
```
Contributors
-------------

[João M Martins (@joaommartins)](https://github.com/joaommartins)

[Micha BA Kunze (@mbakunze)](https://github.com/mbakunze)

[Ramon Crehuet (@rcrehuet)](https://github.com/rcrehuet)

[Giulio Tesei (@gitesei)](https://github.com/gitesei)

[Daniele Montepietra (@Monte95)](https://github.com/Monte95)

[Kresten Lindorff-Larsen (@lindorff-larsen)](https://github.com/lindorff-larsen)
