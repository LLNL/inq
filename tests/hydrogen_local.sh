#!/bin/sh

set -e #make the script fail if a command fails
set -x #output commands to the terminal

inq clear

inq cell cubic 15.0 bohr finite

inq species add Hloc element H
inq species Hloc file `inq util test-data`/H.blyp-vbc.UPF
inq ions insert Hloc 150.0 -30.0 0.0 bohr

inq electrons cutoff 40 hartree
inq theory non interacting
inq ground-state tolerance 1e-8
inq ground-state max-steps 400

inq run gs

inq util match `inq results ground-state energy total`           -0.499022722967 3e-5
inq util match `inq results ground-state energy kinetic`          0.488268502855 3e-5
inq util match `inq results ground-state energy eigenvalues`     -0.499022722967 3e-5
inq util match `inq results ground-state energy hartree`          0.000000000000 3e-5
inq util match `inq results ground-state energy external`        -0.987291225822 3e-5
inq util match `inq results ground-state energy non-local`        0.000000000000 3e-5
inq util match `inq results ground-state energy xc`               0.000000000000 3e-5
inq util match `inq results ground-state energy nvxc`             0.000000000000 3e-5
inq util match `inq results ground-state energy exact-exchange`   0.000000000000 3e-5
inq util match `inq results ground-state energy ion`              0.000000000000 3e-5

inq theory lda
inq run gs

inq util match `inq results ground-state energy total`           -0.445160072256 3e-5
inq util match `inq results ground-state energy kinetic`          0.414315464604 3e-5
inq util match `inq results ground-state energy eigenvalues`     -0.234029035766 3e-5
inq util match `inq results ground-state energy hartree`          0.281309132025 3e-5
inq util match `inq results ground-state energy external`        -0.909266382097 3e-5
inq util match `inq results ground-state energy non-local`        0.000000000000 3e-5
inq util match `inq results ground-state energy xc`              -0.231518286787 3e-5
inq util match `inq results ground-state energy nvxc`            -0.301696382322 3e-5
inq util match `inq results ground-state energy exact-exchange`   0.000000000000 3e-5
inq util match `inq results ground-state energy ion`              0.000000000000 3e-5

inq theory b3lyp
inq run gs

inq util match `inq results ground-state energy total`          -0.447429725736  3e-5
inq util match `inq results ground-state energy kinetic`         0.421657099440  7e-5
inq util match `inq results ground-state energy eigenvalues`    -0.247459996333  3e-5
inq util match `inq results ground-state energy hartree`         0.282608756325  3e-5
inq util match `inq results ground-state energy external`       -0.916549935273  7e-5
inq util match `inq results ground-state energy non-local`       0.000000000000  3e-5
inq util match `inq results ground-state energy xc`             -0.206884771385  3e-5
inq util match `inq results ground-state energy nvxc`           -0.261262923463  3e-5
inq util match `inq results ground-state energy exact-exchange` -0.028260874843  3e-5
inq util match `inq results ground-state energy ion`             0.000000000000  3e-5

