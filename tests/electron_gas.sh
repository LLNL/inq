#!/bin/sh

set -e #make the script fail if a command fails
set -x #output commands to the terminal

inq clear
inq cell  0.0 0.5 0.5  0.5 0.0 0.5  0.5 0.5 0.0 scale 10.0 bohr

inq electrons extra states 2
inq electrons extra electrons 18.0
inq electrons cutoff 30 hartree
inq electrons spin unpolarized
inq electrons temperature 300 K

inq theory lda

inq ground-state tolerance 1e-9

inq run ground-state

inq util match `inq energy total`            3.023858111742 1e-5
inq util match `inq energy kinetic`          9.474820236862 1e-5
inq util match `inq energy eigenvalues`      1.054657739286 1e-5
inq util match `inq energy hartree`          0.000000000711 1e-5
inq util match `inq energy external`         0.000000000000 1e-5
inq util match `inq energy non-local`        0.000000000000 1e-5
inq util match `inq energy xc`              -6.450962125831 1e-5
inq util match `inq energy nvxc`            -8.420162498999 1e-5
inq util match `inq energy exact-exchange`   0.000000000000 1e-5
inq util match `inq energy ion`              0.000000000000 1e-5
