#!/bin/sh

set -e #make the script fail if a command fails
set -x #output commands to the terminal

inq clear
inq cell orthorhombic 10.0 10.0 12.0 bohr
inq ions insert N 0 0 -1.1 bohr
inq ions insert N 0 0  1.1 bohr

inq electrons cutoff 40.0 Hartree
inq ground-state tolerance 1e-9
inq theory lda

inq run ground-state

inq util match `inq results ground-state energy total`          -20.642638443188  1e-5
inq util match `inq results ground-state energy kinetic`         13.163479179497  1e-5
inq util match `inq results ground-state energy eigenvalues`     -5.266545171994  1e-5
inq util match `inq results ground-state energy hartree`         14.494239375977  1e-5
inq util match `inq results ground-state energy external`       -39.444396032274  1e-5
inq util match `inq results ground-state energy non-local`       -1.620128440167  1e-5
inq util match `inq results ground-state energy xc`              -5.718398061370  1e-5
inq util match `inq results ground-state energy nvxc`            -6.353978631002  1e-5
inq util match `inq results ground-state energy exact-exchange`   0.000000000000  1e-5
inq util match `inq results ground-state energy ion`             -1.517434464849  1e-5

inq util match `inq results ground-state forces 0` -0.00000000250910 -0.00000000185763  0.14573092924241  3e-5
inq util match `inq results ground-state forces 1`  0.00000000341961  0.00000000029088 -0.14573050416125  3e-5
#extra check to verify the indiviual component query
inq util match `inq results ground-state forces 1 z` -0.14573050416125  3e-5

inq real-time time-step 0.025 atu
inq real-time num-steps 10
inq real-time ions ehrenfest

inq run real-time


