#!/bin/bash

#allolib/run.sh vdv_group/simulator.cpp

CMAKE_BINARY=cmake
GENERATOR="Unix Makefiles"


mkdir -p build
cd build
mkdir -p Release
cd Release


if [ -d "C:\Program Files (x86)\Microsoft Visual Studio\2017" ]; then
    GENERATOR='Visual Studio 15 2017 Win64'
    CMAKE_BINARY="C:/Program Files (x86)/Microsoft Visual Studio/2017/Community/Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin/cmake.exe"
    if [ ! -f "${CMAKE_BINARY}" ]; then
      CMAKE_BINARY="C:/Program Files (x86)/Microsoft Visual Studio/2017/BuildTools/Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin/cmake.exe"
    fi
    if [ ! -f "${CMAKE_BINARY}" ]; then
      echo Trying to use cmake on PATH as Visual Studio Cmake not found
      CMAKE_BINARY="cmake.exe"
    fi

    echo Building for Visual Studio 15 2017 Win64

    set -x
    "${CMAKE_BINARY}" -G "${GENERATOR}" -DCMAKE_GENERATOR_INSTANCE="C:\Program Files (x86)\Microsoft Visual Studio\2017\Community" -DCMAKE_BUILD_TYPE="Release" ../..
else
    set -x
    "${CMAKE_BINARY}" -G "${GENERATOR}" -DCMAKE_BUILD_TYPE="Release" ../..


fi

"${CMAKE_BINARY}" --build . --config Release --target casm_viewer_run

# if [ $? == 0 ]; then
#     cd ../bin
#     ./casm_viewer
# fi
