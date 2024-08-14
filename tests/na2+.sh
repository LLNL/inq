#!/bin/sh

set -e #make the script fail if a command fails
set -x #output commands to the terminal

inq clear
inq cell cubic 4 angstrom

dist=2.01

inq ions insert Na  0.0  0.0 -$dist/2 angstrom
inq ions insert Na  0.0  0.0  $dist/2 angstrom

inq species Na file `inq util test-data`/Na.pz-n-vbc.UPF 

inq electrons cutoff 30 hartree
inq electrons extra electrons -1
inq ground-state tolerance 1e-8

inq run ground-state

inq util match `inq results ground-state energy total`            -0.643838412931  3e-5
inq util match `inq results ground-state energy kinetic`           0.018892092880  3e-5
inq util match `inq results ground-state energy eigenvalues`      -0.089748441881  3e-5
inq util match `inq results ground-state energy hartree`           0.000453719133  3e-5
inq util match `inq results ground-state energy external`          0.076259382274  3e-5
inq util match `inq results ground-state energy non-local`         0.002717510413  3e-5
inq util match `inq results ground-state energy xc`               -0.376856947176  3e-5
inq util match `inq results ground-state energy nvxc`             -0.188524865714  3e-5
inq util match `inq results ground-state energy exact-exchange`    0.000000000000  3e-5
inq util match `inq results ground-state energy ion`              -0.365304170454  3e-5
  
inq util match `inq results ground-state forces 0`                 1.39972194211486813973e-09 -1.64167367761154538002e-11 1.49641828458383328859e-03 3e-5
inq util match `inq results ground-state forces 1`                 1.39264459119190319341e-09 -5.34464806097899905402e-10 -1.49641826162808617637e-03 3e-5
