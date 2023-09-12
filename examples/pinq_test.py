from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.nwchem import NWChem
from ase.io import write

h2 = Atoms('H2', positions=[[0, 0, 0], [0.0, 0, 0.7]])

import pinq

energy = pinq.run(h2)

print("Energy ", energy);

