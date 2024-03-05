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

inq util match `inq results rt total-energy    0`          -6.01357152790813587195e-01 3e-5
inq util match `inq results rt total-energy  250`          -6.01357647700651587463e-01 3e-5
inq util match `inq results rt total-energy  500`          -6.01357809869950465931e-01 3e-5
inq util match `inq results rt total-energy  750`          -6.01357900582293214420e-01 3e-5
inq util match `inq results rt total-energy 1000`          -6.01357961479974778563e-01 3e-5
inq util match `inq results rt total-energy 1250`          -6.01358006598494143624e-01 3e-5
inq util match `inq results rt total-energy 1500`          -6.01358042261855518795e-01 3e-5
inq util match `inq results rt total-energy 1750`          -6.01358071666531324162e-01 3e-5
inq util match `inq results rt total-energy 2000`          -6.01358096593422097875e-01 3e-5

inq util match `inq results rt dipole    0 x`              -1.54859253834443923344e-03 3e-5
inq util match `inq results rt dipole  250 x`              -1.55123568896079291370e-03 3e-5
inq util match `inq results rt dipole  500 x`              -1.54923816168807686362e-03 3e-5
inq util match `inq results rt dipole  750 x`              -1.55113772430569762592e-03 3e-5
inq util match `inq results rt dipole 1000 x`              -1.54800224333070143672e-03 3e-5
inq util match `inq results rt dipole 1250 x`              -1.55209413178612428731e-03 3e-5
inq util match `inq results rt dipole 1500 x`              -1.54748415575500157353e-03 3e-5
inq util match `inq results rt dipole 1750 x`              -1.55103034746264861407e-03 3e-5
inq util match `inq results rt dipole 2000 x`              -1.55059505033224389638e-03 3e-5

inq util match `inq results rt dipole    0 y`              -1.54854628723687407985e-03 3e-5
inq util match `inq results rt dipole  250 y`              -1.55115793548575640719e-03 3e-5
inq util match `inq results rt dipole  500 y`              -1.54919118490723763462e-03 3e-5
inq util match `inq results rt dipole  750 y`              -1.55115884115211862503e-03 3e-5
inq util match `inq results rt dipole 1000 y`              -1.54807469337626285068e-03 3e-5
inq util match `inq results rt dipole 1250 y`              -1.55216002500586881460e-03 3e-5
inq util match `inq results rt dipole 1500 y`              -1.54749095690195968660e-03 3e-5
inq util match `inq results rt dipole 1750 y`              -1.55097244567385791555e-03 3e-5
inq util match `inq results rt dipole 2000 y`              -1.55051876218004008011e-03 3e-5

inq util match `inq results rt dipole    0 z`              -3.49869452401776624608e-04 3e-5
inq util match `inq results rt dipole  250 z`               2.06165612343310565568e-02 3e-5
inq util match `inq results rt dipole  500 z`              -1.39623342672730247488e-02 3e-5
inq util match `inq results rt dipole  750 z`              -1.24870703322320798084e-02 3e-5
inq util match `inq results rt dipole 1000 z`               2.11928890486084778311e-02 3e-5
inq util match `inq results rt dipole 1250 z`              -2.20788482789152813723e-03 3e-5
inq util match `inq results rt dipole 1500 z`              -2.07210877923769513920e-02 3e-5
inq util match `inq results rt dipole 1750 z`               1.47500374831711651008e-02 3e-5
inq util match `inq results rt dipole 2000 z`               1.01723392738944738017e-02 3e-5
