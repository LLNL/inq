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

inq util match `inq results ground-state energy total`            -0.602091471326 1e-5
inq util match `inq results ground-state energy kinetic`           0.578156580803 1e-5
inq util match `inq results ground-state energy eigenvalues`      -1.102091471326 1e-5
inq util match `inq results ground-state energy hartree`           0.000000000000 1e-5
inq util match `inq results ground-state energy external`         -1.592916754227 1e-5
inq util match `inq results ground-state energy non-local`        -0.087331297902 1e-5
inq util match `inq results ground-state energy xc`                0.000000000000 1e-5
inq util match `inq results ground-state energy nvxc`              0.000000000000 1e-5
inq util match `inq results ground-state energy exact-exchange`    0.000000000000 1e-5
inq util match `inq results ground-state energy ion`               0.500000000000 1e-5

inq rt num-steps 2000
inq rt time-step 0.07 atu
inq rt observables dipole
inq perturbation kick 0 0 0.01

inq run rt

inq util match `inq results rt total-steps`               2000    1e-16
inq util match `inq results rt total-time`                140.00  1e-12

inq util match `inq results rt total-energy    0`          -6.02038713767435340607e-01 3e-5
inq util match `inq results rt total-energy  250`          -6.02039400615363429203e-01 3e-5
inq util match `inq results rt total-energy  500`          -6.02039563777437947145e-01 3e-5
inq util match `inq results rt total-energy  750`          -6.02039650569539008806e-01 3e-5
inq util match `inq results rt total-energy 1000`          -6.02039708902089576448e-01 3e-5
inq util match `inq results rt total-energy 1250`          -6.02039752399125771554e-01 3e-5
inq util match `inq results rt total-energy 1500`          -6.02039786801750076428e-01 3e-5
inq util match `inq results rt total-energy 1750`          -6.02039815142913492529e-01 3e-5
inq util match `inq results rt total-energy 2000`          -6.02039839201153115233e-01 3e-5

#check the individual component query works too
inq util match `inq results rt dipole    0 x`              -1.54124808281475814305e-03 3e-5
inq util match `inq results rt dipole    0 y`              -1.54124808014932277930e-03 3e-5

inq util match `inq results rt dipole    0`                -1.54124808281475814305e-03	-1.54124808014932277930e-03	-3.18539677983979620680e-04 3e-5
inq util match `inq results rt dipole  250`                -1.54615234532997779422e-03	-1.54615234914674175640e-03	2.14187233787422916242e-02 3e-5
inq util match `inq results rt dipole  500`                -1.54334128450166500543e-03	-1.54334128304983372870e-03	8.97097561028314602338e-03 3e-5
inq util match `inq results rt dipole  750`                -1.54194820177244642373e-03	-1.54194820262243642399e-03	-1.79777352895105942132e-02 3e-5
inq util match `inq results rt dipole 1000`                -1.54130017749138726116e-03	-1.54130016731222983101e-03	-1.72305596896911057681e-02 3e-5
inq util match `inq results rt dipole 1250`                -1.54237594471584344054e-03	-1.54237595285780125162e-03	1.00418523973392712773e-02 3e-5
inq util match `inq results rt dipole 1500`                -1.53780963116460864942e-03	-1.53780963394002006053e-03	2.11680482416481442753e-02 3e-5
inq util match `inq results rt dipole 1750`                -1.54273831826737956030e-03	-1.54273831975235533112e-03	-1.54542217063192597812e-03 3e-5
inq util match `inq results rt dipole 2000`                -1.54542886241270673552e-03	-1.54542887026533934218e-03	-2.22503596592770538920e-02 3e-5
