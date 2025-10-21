
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

assert pinq.util.match(pinq.results.ground_state.energy.total(),            -0.643814287043, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.kinetic(),           0.018925591948, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.eigenvalues(),      -0.089705264331, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.hartree(),           0.000453291561, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.external(),          0.076260504480, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.non_local(),         0.002721827596, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.xc(),               -0.376871332174, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.nvxc(),             -0.188519771478, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.exact_exchange(),    0.000000000000, 3e-5)
assert pinq.util.match(pinq.results.ground_state.energy.ion(),              -0.365304170454, 3e-5)
  
assert pinq.util.match(pinq.results.ground_state.forces()[0],  [3.74786824263793630366e-11, -1.05853471242772412887e-10,  1.31008348141167883100e-03], 3e-5)
assert pinq.util.match(pinq.results.ground_state.forces()[1],  [6.54160820496192356550e-11, -6.38554194489870691095e-11, -1.31008351907754708524e-03], 3e-5)

pinq.real_time.num_steps(3000)
pinq.real_time.time_step(0.0666, "atu")
pinq.real_time.observables.dipole()
pinq.perturbations.kick([0, 0, 0.01])

pinq.run.real_time()

energy = pinq.results.real_time.total_energy()

print(energy[   0])
print(energy[  10])
print(energy[ 300])
print(energy[ 444])
print(energy[1000])
print(energy[1663])
print(energy[2000])
print(energy[2101])
print(energy[2748])    
print(energy[2937])
print(energy[3000])

assert pinq.util.match(energy[   0], -0.643761678270877 , 3e-5)
assert pinq.util.match(energy[  10], -0.643761678261219 , 3e-5)
assert pinq.util.match(energy[ 300], -0.6437616784360484, 3e-5)
assert pinq.util.match(energy[ 444], -0.6437616785328906, 3e-5)
assert pinq.util.match(energy[1000], -0.6437616789277596, 3e-5)
assert pinq.util.match(energy[1663], -0.643761679451264 , 3e-5)
assert pinq.util.match(energy[2000], -0.6437616795980003, 3e-5)
assert pinq.util.match(energy[2101], -0.6437616796968362, 3e-5)
assert pinq.util.match(energy[2748], -0.6437616800638288, 3e-5)
assert pinq.util.match(energy[2937], -0.6437616801719617, 3e-5)
assert pinq.util.match(energy[3000], -0.643761680176714 , 3e-5)

dipole = pinq.results.real_time.dipole()

print(dipole[   0])
print(dipole[  10])
print(dipole[ 300])
print(dipole[ 444])
print(dipole[1000])
print(dipole[1663])
print(dipole[2000])
print(dipole[2101])
print(dipole[2748])
print(dipole[2937])
print(dipole[3000])

assert pinq.util.match(dipole[   0], [-0.17853907, -0.17853954, -0.22782616], 3e-5)
assert pinq.util.match(dipole[  10], [-0.178539  , -0.17853941, -0.22839654], 3e-5)
assert pinq.util.match(dipole[ 300], [-0.1784962 , -0.17849663, -0.22849769], 3e-5)
assert pinq.util.match(dipole[ 444], [-0.17851319, -0.17851274, -0.22856397], 3e-5)
assert pinq.util.match(dipole[1000], [-0.1785875 , -0.17858712, -0.22776555], 3e-5)
assert pinq.util.match(dipole[1663], [-0.17841327, -0.17841309, -0.22773484], 3e-5)
assert pinq.util.match(dipole[2000], [-0.1786417 , -0.1786419 , -0.2278392 ], 3e-5)
assert pinq.util.match(dipole[2101], [-0.17843316, -0.1784334 , -0.22758917], 3e-5)
assert pinq.util.match(dipole[2748], [-0.17855025, -0.17855   , -0.2274164 ], 3e-5)
assert pinq.util.match(dipole[2937], [-0.17853585, -0.17853627, -0.22748399], 3e-5)
assert pinq.util.match(dipole[3000], [-0.17857928, -0.17857941, -0.22809489], 3e-5)











