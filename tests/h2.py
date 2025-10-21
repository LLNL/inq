import importlib.util
ase_spec = importlib.util.find_spec("ase")
if(ase_spec is None):
  print("Cannot find ASE, this test will be skipped")
  exit()

from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.nwchem import NWChem
from ase.io import write
import _pinq

atoms = Atoms('H2', positions=[[0, 0, 0], [0.0, 0, 0.7]], cell = [3.0, 3.0, 4.0])
atoms.set_pbc(True)

atoms.calc = _pinq.calculator(ecut = 80.0, xc = 'LDA')
atoms.calc.calculate(atoms)
energy = atoms.get_potential_energy()


print("Energy = ", energy);

assert abs(energy - -32.636969831061656) < 3.0e-5

forces = atoms.get_forces()

assert abs(forces[0, 0]) < 3.0e-5
assert abs(forces[0, 1]) < 3.0e-5
assert abs(forces[0, 2] - -4.71251309e+00) < 3.0e-5
assert abs(forces[1, 0]) < 3.0e-5
assert abs(forces[1, 1]) < 3.0e-5
assert abs(forces[1, 2] -  4.71251139e+00) < 3.0e-5

#Run a second time to check it works

atoms.calc.calculate(atoms)
energy = atoms.get_potential_energy()

assert abs(energy - -32.63697189005673) < 3.0e-5

#Test the scf_step function
array = atoms.calc.get_density()
array.fill(0.0)

for x in range(5):
  density = atoms.calc.scf_step(atoms, array);
  print("energy step ", x, " = ", atoms.get_potential_energy());

assert abs(atoms.get_potential_energy() - -35.79364282124202) < 3.0e-5
  
