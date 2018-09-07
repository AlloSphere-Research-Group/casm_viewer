#!/bin/bash

#allolib/run.sh vdv_group/simulator.cpp

mkdir -p build
cd build
mkdir -p Release
cd Release
cmake -DCMAKE_BUILD_TYPE="Release" ../.. 
cmake --build . --config Release --target casm_viewer_run

# if [ $? == 0 ]; then
#     cd ../bin
#     ./casm_viewer
# fi
