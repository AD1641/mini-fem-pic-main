# mini-pic-sycl

## Introduction
Prototype implementation of the PIC (Particle In Cell) algorithm on unstructured grids.
The particle push (`MoveParticles` in `src/mini-pic.cpp`) runs as a SYCL kernel.

## Building

### macOS (AdaptiveCpp)
```
brew install cmake adaptivecpp
cmake -S . -B build
cmake --build build
```
On a Mac the SYCL kernel runs on the CPU cores. Apple GPUs can't be used as they don't
support double precision.

### Linux cluster (Intel oneAPI)
```
module load intel
module load OpenBLAS
cmake -S . -B build -DCMAKE_CXX_COMPILER=icpx
cmake --build build
```

## Running
```
cd examples/one_stream/coarse
../../../build/mini-pic system.param
```
The FE solver uses dense matrices, so the `fine` and `x_fine` examples need far more
memory than a laptop has.

## Acknowledgements
Based on the `fem-pic` application by Lubos Brieda for Advanced PIC 2015 Lesson 8
  See https://www.particleincell.com/2015/fem-pic/ for more information
