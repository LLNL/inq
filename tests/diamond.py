from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.nwchem import NWChem
from ase.io import write
import pinq

alat =  3.567095

atoms = Atoms('C2', positions = [[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]])
atoms.set_cell([[0.0, alat/2.0, alat/2.0],
                [alat/2.0, 0.0, alat/2.0],
                [alat/2.0, alat/2.0, 0.0]], scale_atoms=True)

energy = pinq.run(atoms)

print("Energy = ", energy);

assert abs(energy - -10.949196617732) < 3.0e-5
