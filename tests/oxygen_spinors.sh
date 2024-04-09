#!/bin/sh

set -e #make the script fail if a command fails
set -x #output commands to the terminal

distance=121.0 #pm

inq clear
inq cell cubic 10.0 bohr finite
inq ions insert O -$distance/2 0 0 pm
inq ions insert O  $distance/2 0 0 pm

inq electrons extra states 4
inq electrons spacing 0.43 bohr
inq electrons spin polarized
inq electrons temperature 1000 K

inq ground-state mixing 0.2
inq ground-state tolerance 1e-8

inq electrons spin non collinear
inq theory non interacting

inq run ground-state

inq util match `inq results ground-state energy total`           -72.587780837035 1e-5
inq util match `inq results ground-state energy kinetic`          35.942275641284 1e-5
inq util match `inq results ground-state energy eigenvalues`     -88.331896202714 1e-5
inq util match `inq results ground-state energy hartree`           0.000000000000 1e-5
inq util match `inq results ground-state energy external`       -116.686239321533 1e-5
inq util match `inq results ground-state energy non-local`        -7.587932522466 1e-5
inq util match `inq results ground-state energy xc`                0.000000000000 1e-5
inq util match `inq results ground-state energy nvxc`              0.000000000000 1e-5
inq util match `inq results ground-state energy exact-exchange`    0.000000000000 1e-5
inq util match `inq results ground-state energy ion`              15.744115365679 1e-5

inq util match  `inq results ground-state magnetization`  0.000000000000 0.000000000000 0.000000001822  1e-5
