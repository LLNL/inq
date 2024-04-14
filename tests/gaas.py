import pinq

alat = 565.315

pinq.clear()
pinq.cell.lattice([0.0, 1.0/2.0, 1.0/2.0], [1.0/2.0, 0.0, 1.0/2.0], [1.0/2.0, 1.0/2.0, 0.0], alat, "pm")
pinq.cell.show()

