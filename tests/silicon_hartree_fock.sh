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
inq electrons extra-states 4

inq ground-state tolerance 1e-4
inq run ground-state

inq theory Hartree-Fock
inq ground-state tolerance 1e-8
inq run ground-state

inq util match `inq result energy total`          -30.503810445836 3e-5
inq util match `inq result energy kinetic`         13.263456976312 3e-5
inq util match `inq result energy eigenvalues`     -6.184068316670 3e-5
inq util match `inq result energy hartree`          2.508149586725 3e-5
inq util match `inq result energy external`        -9.273071664516 3e-5
inq util match `inq result energy non-local`        4.153303103402 3e-5
inq util match `inq result energy xc`               0.000000000000 3e-5
inq util match `inq result energy nvxc`             0.000000000000 3e-5
inq util match `inq result energy exact-exchange`  -9.672027952659 3e-5
inq util match `inq result energy ion`            -31.483620495100 3e-5
