import pinq

alat = 565.315

pinq.clear()

pinq.cell.lattice([0.0, 1.0/2.0, 1.0/2.0], [1.0/2.0, 0.0, 1.0/2.0], [1.0/2.0, 1.0/2.0, 0.0], alat, "pm")
pinq.cell.show()

pinq.ions.insert_fractional("Ga", [0.00, 0.00, 0.00])
pinq.ions.insert_fractional("As", [0.25, 0.25, 0.25])
pinq.ions.show()

pinq.electrons.cutoff(30.0, "Hartree")
pinq.electrons.extra_states(2)
pinq.electrons.temperature(300.0, "Kelvin")
pinq.electrons.show()

pinq.run.ground_state()
