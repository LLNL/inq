#!/bin/sh

set -e #make the script fail if a command fails
set -x #output commands to the terminal

inq clear

inq cell orthorhombic 10.0 10.0 12.0 bohr

#inq species N file `inq util test-data`/N_diagonal.upf.gz
inq species N file `inq util test-data`/N_non_diagonal.upf.gz

inq ions insert N 0 0 -1.1 bohr
inq ions insert N 0 0  1.1 bohr

inq electrons cutoff 40.0 Hartree
inq ground-state tolerance 1e-9

inq run ground-state

inq util match `inq results ground-state energy total`          -20.763133982527  1e-5
inq util match `inq results ground-state energy kinetic`         13.230584158096  1e-5
inq util match `inq results ground-state energy eigenvalues`     -5.238706955641  1e-5
inq util match `inq results ground-state energy hartree`         14.564582838069  1e-5
inq util match `inq results ground-state energy external`       -39.563429064622  1e-5
inq util match `inq results ground-state energy non-local`       -1.621656527689  1e-5
inq util match `inq results ground-state energy xc`              -5.855780921531  1e-5
inq util match `inq results ground-state energy nvxc`            -6.413371197564  1e-5
inq util match `inq results ground-state energy exact-exchange`   0.000000000000  1e-5
inq util match `inq results ground-state energy ion`             -1.517434464849  1e-5

inq util match `inq results ground-state forces 0` -3.39800762999163621530e-08 -4.21756503390703733819e-08  1.42438586921724785750e-01  3e-5
inq util match `inq results ground-state forces 1` -5.79904200238517265052e-08  1.35365826593040522793e-08 -1.42438160496473564809e-01  3e-5
