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

inq util match `inq energy total`           -13.415882369125 3e-5
inq util match `inq energy kinetic`           9.561946756363 3e-5
inq util match `inq energy eigenvalues`      -0.946744600841 3e-5
inq util match `inq energy hartree`           1.768191600618 3e-5
inq util match `inq energy external`         -7.925687586347 3e-5
inq util match `inq energy non-local`        -1.120110161748 3e-5
inq util match `inq energy xc`               -4.410160558972 3e-5
inq util match `inq energy nvxc`             -4.999276810345 3e-5
inq util match `inq energy exact-exchange`    0.000000000000 3e-5
inq util match `inq energy ion`             -11.290062419039 3e-5

