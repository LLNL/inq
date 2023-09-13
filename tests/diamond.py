import importlib.util
ase_spec = importlib.util.find_spec("ase")
if(ase_spec is None):
  print("Cannot find ASE, this test will be skipped")
  exit()

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

atoms.calc = pinq.calculator(ecut = 70.0, extra_bands = 1)
atoms.calc.calculate(atoms)
energy = atoms.get_potential_energy()

print("Energy = ", energy);

assert abs(energy - 2.0*-10.949196617732) < 3.0e-5
