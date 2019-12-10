# INQ

[Xavier Andrade](mailto:xavier@llnl.gov) (LLNL), [Alfredo A. Correa](mailto:correaa@llnl.gov) (LLNL)


INQ is a set of library components designed to build a parallel, portable code for electronic structure calculations, including ground state and excited states.
It concentrates on algorithms as a portability hardware layer, generic types and flexibility of usage.

## Components

* Electronic states, complex and real fields
* Solvers
* Diagonalization
* CPU/GPU generic algorithms
* Parallel distribution MPI/GPU/threads

## Basic installation

```
sudo apt install libgsl-dev libblas-dev liblapack-dev libfftw3-dev
```

Instructions for compiling with Cuda

```bash
git clone ...
cd inq
autoreconf -i
cd your_build
export CXX="/usr/lib/cuda/bin/nvcc -x cu"
export CXXFLAGS="-D_DISABLE_CUDA_SLOW -O3 -std=c++14 --expt-relaxed-constexpr --compiler-options -std=c++14,-Wall,-Wfatal-errors"
export CXXLD=/usr/lib/cuda/bin/nvcc
../../inq/configure -prefix=$HOME --enable-cuda --with-cuda-prefix
```

This instructions might be incomplete, to see how to have a basic install in a standard distribution see [`.gitlab-ci.yml`](https://gitlab.com/npnq/inq/blob/master/.gitlab-ci.yml).