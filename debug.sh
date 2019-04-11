#!/bin/bash

#allolib/run.sh vdv_group/simulator.cpp

mkdir -p build
cd build
mkdir -p Debug
cd Debug
cmake -DCMAKE_BUILD_TYPE="Debug" ../.. 
cmake --build . --config Debug --target casm_viewer_run_debug

# if [ $? == 0 ]; then
#     cd ../bin
#     ./casm_viewer
# fi
