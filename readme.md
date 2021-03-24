# CASM Viewer

Interactive visualization of CASM output.

![CASM Viewer-Grand Canonical](images/GC.PNG)

![CASM Viewer-KMC](images/CASM_viewer-KMC-cumulative-path.PNG)

# Keyboard Commands

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
    pythonScriptsPath - The path to the python scripts folder for this project
    pythonBinary - The python binary to be sued for scripts

View options:

    sliceBackground - Background color for slice (RGBA)
    fontSize = 32.000000000000000
    background - General bakcgroun color (RGBA)

# Building

Make sure you run:

    git submodule update --init --recursive
    
To get all the needed submodules.

You can then use the ```run.sh``` script inside a bash shell (on Windows, Linux or macOS) or open the CMakeLists.txt file in an IDE that supports cmake.

# Running

You will need python3, numpy and netcdf4 installed and accessible. You can set the python 
to be used in the Help tab.

## Windows
For anaconda on Windows use:

    c:/Users/YourName/anaconda3/Scripts/conda.exe run -n base python

There is a bug in anaconda that will prevent it from running if the path to conda
has spaces. In this case, you can install another version of python like the
one provided in the Microsoft Store and install numpy in a shell through pip.
