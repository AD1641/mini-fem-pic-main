# mini-pic
module load intel

module load OpenBLAS


cd mini-pic-fem-main

mkdir build

cd build

cmake -DCMAKE_CXX_COMPILER=icpx ..

make


cd mini-pic-fem-main/examples/one_stream/coarse

../../../build/mini-pic system.param
## Introduction
Prototype implementation of the PIC (Particle In Cell) algorithm on unstructured grids

## Acknowledgements
Based on the `fem-pic` application by Lubos Brieda for Advanced PIC 2015 Lesson 8
  See https://www.particleincell.com/2015/fem-pic/ for more information
# mini-pic-sycl


rm -rf build; mkdir build; cd build; cmake -DCMAKE_CXX_COMPILER=icpx ..; make