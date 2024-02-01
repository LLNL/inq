#!/bin/sh

set -e #make the script fail if a command fails
set -x #output commands to the terminal

inq clear
inq cell cubic 9.0 bohr finite
inq ions add H 0.0 0.0 0.0 bohr
inq electrons cutoff 25 hartree
inq electrons polarized
inq run ground_state
