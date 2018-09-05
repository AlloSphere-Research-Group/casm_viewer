#!/bin/bash

#allolib/run.sh vdv_group/simulator.cpp

mkdir -p build
cd build
cmake -DCMAKE_GENERATOR_PLATFORM=x64 ..
cmake --build .

if [ $? == 0 ]; then
    cd ../bin
    ./casm_viewer
fi
