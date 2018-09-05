#allolib/run.sh vdv_group/simulator.cpp

mkdir -p build
cd build
cmake -DCMAKE_GENERATOR_PLATFORM=x64 ..
cmake --build .

cd ../bin
./casm_viewer
