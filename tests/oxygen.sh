#!/bin/sh

set -e #make the script fail if a command fails
set -x #output commands to the terminal

distance=121.0 #pm

inq clear
inq cell cubic 10.0 bohr finite
inq ions insert O -$distance/2 0 0 pm
inq ions insert O  $distance/2 0 0 pm

inq electrons extra states 2
inq electrons spacing 0.43 bohr
inq electrons spin polarized
inq electrons temperature 1000 K

inq ground-state mixing 0.2
inq ground-state tolerance 1e-8

inq run ground-state

inq util match `inq result energy total`          -32.885878268619 1e-5
inq util match `inq result energy kinetic`         20.663840528309 1e-5
inq util match `inq result energy eigenvalues`     -7.271346087686 1e-5
inq util match `inq result energy hartree`         42.107413409441 1e-5
inq util match `inq result energy external`       -98.967748293591 1e-5
inq util match `inq result energy non-local`       -4.258411059856 1e-5
inq util match `inq result energy xc`              -8.175088218600 1e-5
inq util match `inq result energy nvxc`            -8.923854081430 1e-5
inq util match `inq result energy exact-exchange`   0.000000000000 1e-5
inq util match `inq result energy ion`             15.744115365679 1e-5

inq util match `inq result magnetization x` 0.000000000000 1e-5
inq util match `inq result magnetization y` 0.000000000000 1e-5
inq util match `inq result magnetization z` 1.999999000000 1e-5

inq electrons spin non collinear
inq theory non interacting

inq run ground-state

