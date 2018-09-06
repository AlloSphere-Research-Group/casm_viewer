#!/bin/bash

#allolib/run.sh vdv_group/simulator.cpp

mkdir -p build
cd build
cmake ..
cmake --build .

if [ $? == 0 ]; then
    cd ../bin
    ./casm_viewer
fi
