
import pinq

pinq.clear()

dist = 2.01

pinq.cell.cubic(4.0, "angstrom")
pinq.cell.status()

pinq.ions.insert("Na", [0.0, 0.0, -dist/2.0], "angstrom")
pinq.ions.insert("Na", [0.0, 0.0,  dist/2.0], "angstrom")

pinq.species.file("Na", pinq.util.test_data() + "/Na.pz-n-vbc.UPF");

pinq.electrons.cutoff(30.0, "Hartree")
pinq.electrons.extra_electrons(-1)
pinq.ground_state.tolerance(1e-8)

pinq.run.ground_state()

assert pinq.util.match(pinq.results.ground_state.energy.total(),            -0.643838412931, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.kinetic(),           0.018892092880, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.eigenvalues(),      -0.089748441881, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.hartree(),           0.000453719133, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.external(),          0.076259382274, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.non_local(),         0.002717510413, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.xc(),               -0.376856947176, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.nvxc(),             -0.188524865714, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.exact_exchange(),    0.000000000000, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.ion(),              -0.365304170454, 3e-5)
  
assert pinq.util.match(pinq.results.ground_state.forces()[0][0],  1.39972194211486813973e-09, 3e-5)
assert pinq.util.match(pinq.results.ground_state.forces()[0][1], -1.64167367761154538002e-11, 3e-5)
assert pinq.util.match(pinq.results.ground_state.forces()[0][2],  1.49641828458383328859e-03, 3e-5)
assert pinq.util.match(pinq.results.ground_state.forces()[1][0],  1.39264459119190319341e-09, 3e-5)
assert pinq.util.match(pinq.results.ground_state.forces()[1][1], -5.34464806097899905402e-10, 3e-5)
assert pinq.util.match(pinq.results.ground_state.forces()[1][2], -1.49641826162808617637e-03, 3e-5)
