#!/bin/sh

set -e #make the script fail if a command fails
set -x #output commands to the terminal

inq clear

inq ions file `inq util test-data`/H.cif

inq electrons cutoff 30 hartree
inq theory pbe
inq ground-state tolerance 1e-8

inq run ground-state

inq util match `inq energy total`           -2.355128752518 3e-5
inq util match `inq energy kinetic`          2.125707369840 3e-5
inq util match `inq energy eigenvalues`     -1.511390321649 3e-5
inq util match `inq energy hartree`          0.903321603071 3e-5
inq util match `inq energy external`        -3.319351558502 3e-5
inq util match `inq energy non-local`       -0.389801056728 3e-5
inq util match `inq energy xc`              -1.351311210649 3e-5
inq util match `inq energy nvxc`            -1.734588282400 3e-5
inq util match `inq energy exact-exchange`   0.000000000000 3e-5
inq util match `inq energy ion`             -0.323693899550 3e-5

#TODO: check the eigenvalues to be -0.401182208666 and -0.354512952159

