# Keyboard Commands

    1-9    - Change View
    [,]    - Scale layer view
    -, =   - prev/next layer in perspective view
    O, P   - Temperature prev/next
    L, ;   - Chempot A prev/next
    . , /  - Chempot B prev/next
    Q, W   - Preset prev/next
    B, V   - Perspective rotate step up/down
    A, S   - Start/Stop sequence
    Shift-A, (Shift-)S - Start/stop record
    SPACE - Toggle GUI visibility

# Configuration
The *casm_viewer.toml* file is used to configure the application. There you can set system specific paths and dataset paths and settings.

    datasets = ["cooling", "heating_12"] - The list of datasets to explore
    pythonScriptsPath - The path to the python scripts folder for this project
    pythonBinary - The python binary to be sued for scripts

View options:

    sliceBackground - Background color for slice (RGBA)
    fontSize = 32.000000000000000
    background - General bakcgroun color (RGBA)

# Downloading

(Only follow these building and downloading instructions if not using the casm_viewer repo. Updating and build should be done with the scripts provided in that repo)

You need allolib:

    git clone https://github.com/AlloSphere-Research-Group/allolib.git
    cd allolib
    git submodule init
    git submodule update

The VdV repos are at nonce. You need special access for this. Let mantaraya36@gmail.com know. Clone them inside the allolib repo root folder.

    git clone arg@nonce.mat.ucsb.edu:vdv_group.git
    git clone arg@nonce.mat.ucsb.edu:vdv_group_python.git

# Building


## Linux and OS X
 Clone this repo inside an allolib repo

Remember to initialize and update the submodules within allolib.

Then type : 

     ./run.sh vdv_group/simulator.cpp

This should build and run the application.

## Windows
### Visual Studio 2017
Open the CMakeLists.txt file found in the windows folder in Visual Studio. Set build to 64 bits. On the cmake menu select 'Build All' and cross your fingers...
