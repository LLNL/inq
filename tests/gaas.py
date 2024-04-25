import pinq

alat = 565.315

pinq.clear()

pinq.cell.lattice([0.0, 1.0/2.0, 1.0/2.0], [1.0/2.0, 0.0, 1.0/2.0], [1.0/2.0, 1.0/2.0, 0.0], alat, "pm")
pinq.cell.status()

pinq.ions.insert_fractional("Ga", [0.00, 0.00, 0.00])
pinq.ions.insert_fractional("As", [0.25, 0.25, 0.25])
pinq.ions.insert_fractional("As", [0.25, 0.25, 0.25])
pinq.ions.remove(1);
pinq.ions.status()

pinq.electrons.cutoff(30.0, "Hartree")
pinq.electrons.extra_states(2)
pinq.electrons.temperature(300.0, "Kelvin")
pinq.electrons.status()

pinq.kpoints.shifted_grid(2, 2, 2)

pinq.ground_state.tolerance(1e-9)
pinq.ground_state.status()

pinq.run.ground_state()

pinq.results.ground_state.status()
pinq.results.ground_state.energy.status()
