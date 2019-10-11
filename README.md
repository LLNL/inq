# inq

This is not a PW code

## Basic installation

```
sudo apt install libgsl-dev libblas-dev liblapack-dev libfftw3-dev
```

```bash
git clone ...
cd inq
autoreconf -i
cd your_build
export CXX="/usr/lib/cuda/bin/nvcc -x cu"
export CXXFLAGS="-D_DISABLE_CUDA_SLOW -O3 -std=c++14 --expt-relaxed-constexpr --compiler-options -std=c++14,-Wall,-Wfatal-errors"
export CXXLD=/usr/lib/cuda/bin/nvcc
../../inq/configure -prefix-$HOME --enable-cuda --with-cuda-prefix
```