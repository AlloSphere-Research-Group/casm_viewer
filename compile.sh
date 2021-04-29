mkdir -p build/release
cd build/release
cmake -DCMAKE_BUILD_TYPE=Release -Wno-deprecated ../..
# cmake -G Xcode -DCMAKE_BUILD_TYPE=Release -Wno-deprecated ../..
cd ../..
cmake --build build/release -j7
cp res/macos/Info.plist bin/casm_viewer.app/Contents/
mkdir bin/casm_viewer.app/Contents/Resources
cp res/macos/myrioi_icon.icns bin/casm_viewer.app/Contents/Resources/
cp res/macos/run bin/casm_viewer.app/Contents/MacOS/
cp build/release/external/tinc/external/allolib/external/rtmidi/RtMidiConfig.cmake build/release