import numpy as np
from ase.calculators.calculator import Calculator, all_changes
import _pinq

class PinqCalculator(Calculator):
    implemented_properties = ['energy', 'forces']

    def __init__(self, pcalc=None, **kwargs):
        Calculator.__init__(self, **kwargs)
        self.kwargs =kwargs
        if pcalc is None:
            self.pcalc = _pinq.calculator(**self.kwargs)

    def calculate(self, atoms=None, properties=['energy'], system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)

        if len(system_changes) > 0 :
            # if 'positions' or 'cell' in system_changes :
            self.pcalc.calculate(self.atoms)
            self.results['energy'] = self.pcalc.get_potential_energy(self.atoms)

        if 'forces' in properties:
            forces = self.pcalc.get_forces(self.atoms)
            self.results['forces'] = np.array(forces).reshape((-1,3))

        if 'stress' in properties:
            print("!WARN: stress not implemented yet.")
            stress_voigt = np.zeros(6)
            self.results['stress'] = stress_voigt
