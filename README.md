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

import MDAnalysis as md
from FRET import FRETpredict

res1 = 452
chain_res1 = 'A'
chromophore_1 = 532

res2 = 637
chain_res2 = 'B'
chromophore_2 = 647

u = MDAnalysis.Universe('test_systems/Hsp90/openHsp90.pdb')

cutoff = 10

sigma = 1
epsilon = 1

FRET = FRETpredict(protein=u, residues=[res1, res2], temperature=293, chains=[chain_res1, chain_res2], 
                   donor=chromophore_1, acceptor=chromophore_2, 
                   sigma_scaling=sigma, epsilon_scaling=epsilon, electrostatic=True,
                   libname_1='Alexa {} cutoff{:d}'.format(chromophore_1, cutoff),
                   libname_2='Alexa {} cutoff{:d}'.format(chromophore_2, cutoff), 
                   output_prefix='prova/E{:d}_{}_{}'.format(cutoff, sigma, epsilon))

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
