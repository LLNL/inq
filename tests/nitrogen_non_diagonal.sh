#!/bin/sh

set -e #make the script fail if a command fails
set -x #output commands to the terminal

inq clear

inq cell orthorhombic 10.0 10.0 12.0 bohr

inq species N file `inq util test-data`/N_diagonal.upf.gz

inq ions insert N 0 0 -1.1 bohr
inq ions insert N 0 0  1.1 bohr

inq electrons cutoff 40.0 Hartree
inq ground-state tolerance 1e-9

inq run ground-state

printf "*****************************\n\n  Checking diagonal results\n\n*****************************\n"
inq util match `inq results ground-state energy total`          -20.763133982527  1e-5
inq util match `inq results ground-state energy non-local`       -1.621656527689  1e-5
inq util match `inq results ground-state energy kinetic`         13.230584158096  1e-5
inq util match `inq results ground-state energy eigenvalues`     -5.238706955641  1e-5
inq util match `inq results ground-state energy hartree`         14.564582838069  1e-5
inq util match `inq results ground-state energy external`       -39.563429064622  1e-5
inq util match `inq results ground-state energy xc`              -5.855780921531  1e-5
inq util match `inq results ground-state energy nvxc`            -6.413371197564  1e-5
inq util match `inq results ground-state energy exact-exchange`   0.000000000000  1e-5
inq util match `inq results ground-state energy ion`             -1.517434464849  1e-5

inq util match `inq results ground-state forces 0` -3.39800762999163621530e-08 -4.21756503390703733819e-08  1.42438586921724785750e-01  3e-5
inq util match `inq results ground-state forces 1` -5.79904200238517265052e-08  1.35365826593040522793e-08 -1.42438160496473564809e-01  3e-5

inq species N file `inq util test-data`/N_non_diagonal.upf.gz

inq run ground-state

printf "*********************************\n\n  Checking non-diagonal results\n\n*********************************\n"
inq util match `inq results ground-state energy total`           -20.763189284257 1e-5
inq util match `inq results ground-state energy non-local`        -1.622238914430 1e-5
inq util match `inq results ground-state energy kinetic`          13.231118947204 1e-5
inq util match `inq results ground-state energy eigenvalues`      -5.238721752365 1e-5
inq util match `inq results ground-state energy hartree`          14.564651460870 1e-5
inq util match `inq results ground-state energy external`        -39.563500607215 1e-5
inq util match `inq results ground-state energy xc`               -5.855785705837 1e-5
inq util match `inq results ground-state energy nvxc`             -6.413404099665 1e-5
inq util match `inq results ground-state energy exact-exchange`    0.000000000000 1e-5
inq util match `inq results ground-state energy ion`              -1.517434464849 1e-5

inq util match `inq results ground-state forces 0` -3.72931746878889079665e-11 -9.81254878920072358131e-11  1.42687266557671499356e-01  3e-4
inq util match `inq results ground-state forces 1` -4.37555277175191418715e-11  1.03612913819364776335e-10 -1.42687264686675835401e-01  3e-4

