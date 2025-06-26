# MorphCut
## Usage
### Step 1 - Install dependencies
- CGAL (ver 6.0) (Use the modified version in the uploaded cgal.zip)
- Qt6 (For visualization in runtime)
- Eigen
- RapidJSON
- NLopt
- SDF (Use the modified version in the uploaded sdf.zip)
- OpenMP
### Step2 - Build from source
type the commands in the terminal:
cd code & mkdir build & cd build & cmake .. & make
### Step 3 - Run the compiled executable
now you are in the build directory, type in the terminal:
./bcd path-to-sample.off
To test a sample, users can unzip the data.zip, and then keep the folder structure as it is, and then pass the sample.off in the folder to the bcd executable. A configuration file customizing the parameters is named sample.off_config.json in the same folder.
e.g., ./bcd /test_MorphCut/data/sample1/sample1.off
### Step 4 - Inspect the results
Results are saved in the same directory of the sample.off
## Tested platforms
We mainly tested MorphCut on Linux and MacOS:
- MacOS, Sequoia 15.2
- Windows Subsystem for Linux (WSL)
- Ubuntu 24
- Linux-based High performance facilities
