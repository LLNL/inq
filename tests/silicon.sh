#!/bin/sh

set -e #make the script fail if a command fails
set -x #output commands to the terminal

inq clear
inq cell cubic 10.18 bohr

inq ions insert fractional Si  0.0  0.0  0.0 
inq ions insert fractional Si  0.25 0.25 0.25
inq ions insert fractional Si  0.5  0.5  0.0 
inq ions insert fractional Si  0.75 0.75 0.25
inq ions insert fractional Si  0.5  0.0  0.5 
inq ions insert fractional Si  0.75 0.25 0.75
inq ions insert fractional Si  0.0  0.5  0.5 
inq ions insert fractional Si  0.25 0.75 0.75

inq electrons cutoff 25 hartree
inq theory non-interacting
inq kpoints shifted grid 2 1 1
inq ground-state tolerance 1e-9

inq run ground-state

inq util match `inq energy total`            -23.834202378912 3e-5
inq util match `inq energy kinetic`           14.428062931586 3e-5
inq util match `inq energy eigenvalues`        7.649418116188 3e-5
inq util match `inq energy hartree`            0.000000000000 3e-5
inq util match `inq energy external`         -12.019317709315 3e-5
inq util match `inq energy non-local`          5.240672893917 3e-5
inq util match `inq energy xc`                 0.000000000000 3e-5
inq util match `inq energy nvxc`               0.000000000000 3e-5
inq util match `inq energy exact-exchange`     0.000000000000 3e-5
inq util match `inq energy ion`              -31.483620495100 3e-5
