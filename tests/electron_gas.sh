#!/bin/sh

set -e #make the script fail if a command fails
set -x #output commands to the terminal

#RUN 1

inq clear
inq cell cubic 10.0 bohr periodic

inq electrons extra states 2
inq electrons extra electrons 14.0
inq electrons cutoff 30 hartree
inq electrons spin unpolarized
inq electrons temperature 300 K

inq kpoints grid 1 1 3

inq theory lda

inq ground-state tolerance 1e-9

inq run ground-state

inq util match `inq result energy total`           -0.567967486012 3e-5
inq util match `inq result energy kinetic`          2.485678153650 3e-5
inq util match `inq result energy eigenvalues`     -1.488506069666 3e-5
inq util match `inq result energy hartree`          0.000000110256 3e-5
inq util match `inq result energy external`         0.000000000000 3e-5
inq util match `inq result energy non-local`        0.000000000000 3e-5
inq util match `inq result energy xc`              -3.053645749918 3e-5
inq util match `inq result energy nvxc`            -3.974184443829 3e-5
inq util match `inq result energy exact-exchange`   0.000000000000 3e-5
inq util match `inq result energy ion`              0.000000000000 3e-5

# RUN 2

inq clear

#This is the same as before but the lattice is rotated 45 degrees along the z axis
a=10
inq cell  "$a/sqrt(2.0)" $a/2 $a/2  "-$a/sqrt(2)" $a/2 $a/2  0 "-$a/sqrt(2.0)" "$a/sqrt(2.0)"  bohr

inq electrons extra states 2
inq electrons extra electrons 14.0
inq electrons cutoff 30 hartree
inq electrons spin unpolarized
inq electrons temperature 300 K

inq kpoints grid 1 1 3

inq theory lda

inq ground-state tolerance 1e-9

inq run ground-state

inq util match `inq result energy total`           -0.567967360088 1e-5
inq util match `inq result energy kinetic`          2.485678164227 1e-5
inq util match `inq result energy eigenvalues`     -1.488505578638 1e-5
inq util match `inq result energy hartree`          0.000000587891 1e-5
inq util match `inq result energy external`         0.000000000000 1e-5
inq util match `inq result energy non-local`        0.000000000000 1e-5
inq util match `inq result energy xc`              -3.053646112207 1e-5
inq util match `inq result energy nvxc`            -3.974184918648 1e-5
inq util match `inq result energy exact-exchange`   0.000000000000 1e-5
inq util match `inq result energy ion`              0.000000000000 1e-5

# RUN 3

inq clear
inq cell  0 1/2 1/2  1/2 0 1/2  1/2 1/2 0 scale 10.0 bohr

inq electrons extra states 2
inq electrons extra electrons 18.0
inq electrons cutoff 30 hartree
inq electrons spin unpolarized
inq electrons temperature 300 K

inq theory lda

inq ground-state tolerance 1e-9

inq run ground-state

inq util match `inq result energy total`            3.023858111742 1e-5
inq util match `inq result energy kinetic`          9.474820236862 1e-5
inq util match `inq result energy eigenvalues`      1.054657739286 1e-5
inq util match `inq result energy hartree`          0.000000000711 1e-5
inq util match `inq result energy external`         0.000000000000 1e-5
inq util match `inq result energy non-local`        0.000000000000 1e-5
inq util match `inq result energy xc`              -6.450962125831 1e-5
inq util match `inq result energy nvxc`            -8.420162498999 1e-5
inq util match `inq result energy exact-exchange`   0.000000000000 1e-5
inq util match `inq result energy ion`              0.000000000000 1e-5
