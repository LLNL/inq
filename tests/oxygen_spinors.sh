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

inq electrons spin non collinear
inq theory non interacting

inq run ground-state

