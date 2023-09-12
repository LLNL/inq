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

atoms = Atoms('H2', positions=[[0, 0, 0], [0.0, 0, 0.7]], cell = [3.0, 3.0, 4.0])

calc = pinq.calculator()
atoms.calc = calc
energy = atoms.get_potential_energy()

print("Energy = ", energy);

assert abs(energy - -2.4473819221017314) < 3.0e-5
