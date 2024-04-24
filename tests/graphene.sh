#!/bin/sh

set -e #make the script fail if a command fails
set -x #output commands to the terminal

dcc=1.42
a=`inq util calc "sqrt(3)*$dcc"`

inq clear
inq cell   $a 0 0  -$a/2 "(sqrt(3)/2)*$a" 0   0 0 5.2917721 Angstrom 2D

inq ions insert H 0 0    0 Angstrom #unnecessary, but just to test remove
inq ions insert C 0 0    0 Angstrom
inq ions insert C 0 $dcc 0 Angstrom
inq ions remove 0 #test remove

inq electrons extra states 2
inq electrons spacing $a/15 Angstrom

inq theory functional XC_GGA_X_RPBE XC_GGA_C_PBE

inq ground-state tolerance 1e-8

inq run ground-state

inq util match `inq results ground-state energy total`          -11.805691437798 3e-5
inq util match `inq results ground-state energy kinetic`          9.569622557902 3e-5
inq util match `inq results ground-state energy eigenvalues`     -3.508150981796 3e-5
inq util match `inq results ground-state energy hartree`        -11.111524722499 3e-5
inq util match `inq results ground-state energy external`        15.082403832279 3e-5
inq util match `inq results ground-state energy non-local`       -1.140029519851 3e-5
inq util match `inq results ground-state energy xc`              -4.392447120595 3e-5
inq util match `inq results ground-state energy nvxc`            -4.797098407128 3e-5
inq util match `inq results ground-state energy exact-exchange`   0.000000000000 3e-5
inq util match `inq results ground-state energy ion`            -19.813716465033 3e-5
