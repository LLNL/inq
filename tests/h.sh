#!/bin/sh

set -e #make the script fail if a command fails
set -x #output commands to the terminal

inq clear
inq cell cubic 9.0 bohr finite
inq ions insert H 0.0 0.0 0.0 bohr
inq electrons cutoff 25 hartree
inq electrons spin polarized
inq ground-state tolerance 1e-8
inq run ground-state

inq util match `inq results ground-state energy total`           -0.503522815653 1e-5
inq util match `inq results ground-state energy kinetic`          0.450527223482 1e-5
inq util match `inq results ground-state energy eigenvalues`     -0.292974894898 1e-5
inq util match `inq results ground-state energy hartree`          0.292837199614 1e-5
inq util match `inq results ground-state energy external`        -0.887386400176 1e-5
inq util match `inq results ground-state energy non-local`       -0.064476885703 1e-5
inq util match `inq results ground-state energy xc`              -0.295023952870 1e-5
inq util match `inq results ground-state energy nvxc`            -0.377313231730 1e-5
inq util match `inq results ground-state energy exact-exchange`   0.000000000000 1e-5
inq util match `inq results ground-state energy ion`              0.000000000000 1e-5
