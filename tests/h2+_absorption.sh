#!/bin/sh

set -e #make the script fail if a command fails
set -x #output commands to the terminal

inq clear
inq cell orthorhombic 8 8 10 bohr finite

inq ions insert H 0.0 0.0 -1.0 bohr
inq ions insert H 0.0 0.0  1.0 bohr

inq electrons extra-electrons -1
inq electrons spacing 0.45 bohr
inq ground-state tolerance 1e-8

inq theory non-interacting

inq run ground-state

inq util match `inq results ground-state energy total`           -0.601409734044 1e-5
inq util match `inq results ground-state energy kinetic`          0.576555647167 1e-5
inq util match `inq results ground-state energy eigenvalues`     -1.101409734044 1e-5
inq util match `inq results ground-state energy hartree`          0.000000000000 1e-5
inq util match `inq results ground-state energy external`        -1.591616828069 1e-5
inq util match `inq results ground-state energy non-local`       -0.086348553142 1e-5
inq util match `inq results ground-state energy xc`               0.000000000000 1e-5
inq util match `inq results ground-state energy nvxc`             0.000000000000 1e-5
inq util match `inq results ground-state energy exact-exchange`   0.000000000000 1e-5
inq util match `inq results ground-state energy ion`              0.500000000000 1e-5

inq rt num-steps 2000
inq rt time-step 0.075 atu
inq rt observables dipole
inq perturbation kick 0 0 0.01

inq run rt

inq util match `inq results rt total-steps`               2000    1e-16
inq util match `inq results rt total-time`                150.00  1e-12

inq util match `inq results rt total-energy    0`          -6.01357152790818916266e-01 3e-5
inq util match `inq results rt total-energy  250`          -6.01357647700651587463e-01 3e-5
inq util match `inq results rt total-energy  500`          -6.01357809869941584147e-01 3e-5
inq util match `inq results rt total-energy  750`          -6.01357900582282778323e-01 3e-5
inq util match `inq results rt total-energy 1000`          -6.01357961479969005403e-01 3e-5
inq util match `inq results rt total-energy 1250`          -6.01358006598493255446e-01 3e-5
inq util match `inq results rt total-energy 1500`          -6.01358042261849523591e-01 3e-5
inq util match `inq results rt total-energy 1750`          -6.01358071666526439181e-01 3e-5
inq util match `inq results rt total-energy 2000`          -6.01358096593417656983e-01 3e-5

#check the individual component query works too
inq util match `inq results rt dipole    0 x`              -1.54859253809264759034e-03 3e-5
inq util match `inq results rt dipole    0 y`              -1.54854628737975110950e-03 3e-5

inq util match `inq results rt dipole    0`                -1.54859253809264759034e-03	-1.54854628737975110950e-03	-3.49869451722191438740e-04 3e-5
inq util match `inq results rt dipole  250`                -1.55123568829897552461e-03	-1.55115793549341239241e-03	 2.06165612349959205540e-02 3e-5
inq util match `inq results rt dipole  500`                -1.54923816113996517738e-03	-1.54919118476384279001e-03	-1.39623342683573240658e-02 3e-5
inq util match `inq results rt dipole  750`                -1.55113772430512039668e-03	-1.55115884097505193039e-03	-1.24870703321999683422e-02 3e-5
inq util match `inq results rt dipole 1000`                -1.54800224387945778957e-03	-1.54807469329935692114e-03	 2.11928890496754437911e-02 3e-5
inq util match `inq results rt dipole 1250`                -1.55209413244008836638e-03	-1.55216002509228840137e-03	-2.20788482860621122425e-03 3e-5
inq util match `inq results rt dipole 1500`                -1.54748415599520895009e-03	-1.54749095708308758573e-03	-2.07210877929851905455e-02 3e-5
inq util match `inq results rt dipole 1750`                -1.55103034708906919147e-03	-1.55097244580162832266e-03  1.47500374842999444625e-02 3e-5
inq util match `inq results rt dipole 2000`                -1.55059504965341212260e-03	-1.55051876215510841747e-03  1.01723392737756209575e-02 3e-5
