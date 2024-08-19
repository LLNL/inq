
import pinq

pinq.clear()

dist = 2.01

pinq.cell.cubic(4.0, "angstrom")
pinq.cell.status()

pinq.ions.insert("Na", [0.0, 0.0, -dist/2.0], "angstrom")
pinq.ions.insert("Na", [0.0, 0.0,  dist/2.0], "angstrom")

pinq.species.file("Na", pinq.util.test_data() + "/Na.pz-n-vbc.UPF");

pinq.electrons.cutoff(24.0, "Hartree")
pinq.electrons.extra_electrons(-1)
pinq.ground_state.tolerance(1e-8)

pinq.run.ground_state()

assert pinq.util.match(pinq.results.ground_state.energy.total(),            -0.643875190035, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.kinetic(),           0.018878577118, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.eigenvalues(),      -0.089771561959, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.hartree(),           0.000453780325, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.external(),          0.076259198549, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.non_local(),         0.002713603553, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.xc(),               -0.376876179126, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.nvxc(),             -0.188530501831, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.exact_exchange(),    0.000000000000, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.ion(),              -0.365304170454, 3e-5)
  
assert pinq.util.match(pinq.results.ground_state.forces()[0][0],  7.26217835896163880644e-11, 3e-5)
assert pinq.util.match(pinq.results.ground_state.forces()[0][1], -8.66073022870166591849e-11, 3e-5)
assert pinq.util.match(pinq.results.ground_state.forces()[0][2],  1.46550554264721099619e-03, 3e-5)
assert pinq.util.match(pinq.results.ground_state.forces()[1][0],  8.64681672462082104618e-11, 3e-5)
assert pinq.util.match(pinq.results.ground_state.forces()[1][1], -5.46191861868080879087e-11, 3e-5)
assert pinq.util.match(pinq.results.ground_state.forces()[1][2], -1.46550559849585658283e-03, 3e-5)

pinq.real_time.num_steps(3000)
pinq.real_time.time_step(0.075, "atu")
pinq.real_time.observables.dipole()
pinq.perturbations.kick([0, 0, 0.01])

pinq.run.real_time()

energy = pinq.results.real_time.total_energy()

assert pinq.util.match(energy[   0], -0.6438225831700184, 3e-5)
assert pinq.util.match(energy[  10], -0.6438225831586412, 3e-5)
assert pinq.util.match(energy[ 300], -0.6438225835346540, 3e-5)
assert pinq.util.match(energy[ 444], -0.6438225837767503, 3e-5)
assert pinq.util.match(energy[1000], -0.6438225844307714, 3e-5)
assert pinq.util.match(energy[1663], -0.6438225852107715, 3e-5)
assert pinq.util.match(energy[2000], -0.6438225854501384, 3e-5)
assert pinq.util.match(energy[2101], -0.6438225856304741, 3e-5)
assert pinq.util.match(energy[2748], -0.6438225861809898, 3e-5)
assert pinq.util.match(energy[2937], -0.6438225863718710, 3e-5)
assert pinq.util.match(energy[3000], -0.6438225864372554, 3e-5)

