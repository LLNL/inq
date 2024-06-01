#!/bin/sh

set -e #make the script fail if a command fails
set -x #output commands to the terminal

inq clear

inq ions file `inq util test-data`/bn.poscar

inq electrons cutoff 35 hartree
inq electrons extra states 3

inq kpoints shifted grid 2 2 2

inq theory pbe

inq ground-state tolerance 1e-8

inq run ground-state

inq util match `inq results ground-state energy total`           -13.415121196063 3e-5
inq util match `inq results ground-state energy kinetic`           9.561233691561 3e-5
inq util match `inq results ground-state energy eigenvalues`      -0.946326812486 3e-5
inq util match `inq results ground-state energy hartree`           1.767819938212 3e-5
inq util match `inq results ground-state energy external`         -7.924084911288 3e-5
inq util match `inq results ground-state energy non-local`        -1.120000193499 3e-5
inq util match `inq results ground-state energy xc`               -4.410027302009 3e-5
inq util match `inq results ground-state energy nvxc`             -4.999115275683 3e-5
inq util match `inq results ground-state energy exact-exchange`    0.000000000000 3e-5
inq util match `inq results ground-state energy ion`             -11.290062419039 3e-5
