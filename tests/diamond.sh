#!/bin/sh

set -e #make the script fail if a command fails
set -x #output commands to the terminal

inq clear
inq cell  0.0 0.5 0.5  0.5 0.0 0.5  0.5 0.5 0.0 scale 3.567095 A

inq ions insert fractional C 0.00 0.00 0.00
inq ions insert fractional C 0.25 0.25 0.25

inq electrons cutoff 35.0 Ha
inq electrons extra-states 3

inq kpoints gamma
inq kpoints insert -0.5 -0.5 -0.5 0.0

inq ground-state tolerance 1e-8

inq run ground state

inq util match `inq result energy total`            -10.949196120641 3e-5
inq util match `inq result energy kinetic`           11.411454824474 3e-5
inq util match `inq result energy eigenvalues`        0.795475249124 3e-5
inq util match `inq result energy hartree`            1.473883894765 3e-5
inq util match `inq result energy external`          -7.034444373077 3e-5
inq util match `inq result energy non-local`         -1.496762228839 3e-5
inq util match `inq result energy xc`                -4.568604109362 3e-5
inq util match `inq result energy nvxc`              -5.032540762964 3e-5
inq util match `inq result energy exact-exchange`     0.000000000000 3e-5
inq util match `inq result energy ion`              -10.734724128603 3e-5

