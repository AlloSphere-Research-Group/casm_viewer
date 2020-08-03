#!/bin/bash

CMAKE_BINARY=cmake
GENERATOR="Unix Makefiles"

mkdir -p build
cd build
mkdir -p Release
cd Release


if [ -d "C:\Program Files (x86)\Microsoft Visual Studio\2017" ]; then
  GENERATOR='Visual Studio 15 2017 Win64'
  
  VS_CMAKE_BINARY="C:/Program Files (x86)/Microsoft Visual Studio/2017/BuildTools/Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin/cmake.exe"
  GENERATOR="Visual Studio 15 2017 Win64"
  echo Tryng VS 2017 build.
  if [ -f "${VS_CMAKE_BINARY}" ]; then
    echo Using VS 2017 CMaake
    CMAKE_BINARY=$VS_CMAKE_BINARY
  fi

  set -x
  "${CMAKE_BINARY}" -G "${GENERATOR}" -DCMAKE_GENERATOR_INSTANCE="C:\Program Files (x86)\Microsoft Visual Studio\2017\Community" -DCMAKE_BUILD_TYPE="Release" ../..
  set +x
elif [ -d "C:\Program Files (x86)\Microsoft Visual Studio\2019" ]; then

  VS_CMAKE_BINARY="C:/Program Files (x86)/Microsoft Visual Studio/2019/Community/Common7/IDE/CommonExtensions/Microsoft/CMake/CMake/bin/cmake.exe"
  GENERATOR="Visual Studio 16 2019"
  echo Tryng VS 2019 build.
    if [ ! -f "${CMAKE_BINARY}" ]; then
    echo Using VS 2019 CMaake
    CMAKE_BINARY=$VS_CMAKE_BINARY
  fi
  set -x
  "${CMAKE_BINARY}" -G "${GENERATOR}" -DTINC_BUILD_EXAMPLES=ON -DCMAKE_BUILD_TYPE="Release" ../..
  set +x
else
  set -x
  "${CMAKE_BINARY}" -G "${GENERATOR}" -DTINC_BUILD_EXAMPLES=ON -DCMAKE_BUILD_TYPE="Release" ../..
  set +x
fi

set -x
"${CMAKE_BINARY}" --build . --config Release --target casm_viewer -j7
set +x

if [ $? == 0 ]; then
    cd ../../bin
    ./casm_viewer
fi
