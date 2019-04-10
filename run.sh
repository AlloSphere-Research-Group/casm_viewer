#!/bin/bash

#allolib/run.sh vdv_group/simulator.cpp

CMAKE_BINARY=cmake

if [ -d "C:\Program Files (x86)\Microsoft Visual Studio\2019" ]; then
  GENERATOR="Visual Studio 15 2017 Win64"
  CMAKE_BINARY="C:/Program Files (x86)/Microsoft Visual Studio/2017/Community/Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin/cmake.exe"
  echo Building for Visual Studio 16 2019 Win64
elif [ -d "C:\Program Files (x86)\Microsoft Visual Studio\2017" ]; then
  GENERATOR="-G 'Visual Studio 15 2017 Win64'"
  CMAKE_BINARY="C:/Program Files (x86)/Microsoft Visual Studio/2017/Community/Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin/cmake.exe"
  echo Building for Visual Studio 15 2017 Win64
fi

mkdir -p build
cd build
mkdir -p Release
cd Release

set -x
"${CMAKE_BINARY}" -G "${GENERATOR}" -DCMAKE_GENERATOR_INSTANCE="C:\Program Files (x86)\Microsoft Visual Studio\2017\Community" -DCMAKE_BUILD_TYPE="Release" ../..
"${CMAKE_BINARY}" --build . --config Release --target casm_viewer_run

# if [ $? == 0 ]; then
#     cd ../bin
#     ./casm_viewer
# fi
