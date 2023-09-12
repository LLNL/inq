from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.nwchem import NWChem
from ase.io import write

print("*** Starting pinq test ***")

h2 = Atoms('H2', positions=[[0, 0, 0], [0.0, 0, 0.7]], cell = [3.0, 3.0, 4.0])

print(h2.get_cell()[:])

import pinq

energy = pinq.run(h2)

print("Energy ", energy);
