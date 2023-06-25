# Averaged Minimum Distances
Averaged Minimum Distances (AMD) encode crystal structures to distance vectors. More details can be found in our paper published in MATCH Communications in Mathematical and in Computer Chemistry Journal: http://dx.doi.org/10.46793/match.87-3.529W

The project has been compiled and run only on Windows x64.

For a more recent version written in Python, refer to the github link of the first author https://github.com/dwiddo/average-minimum-distance.

## Installation
1. Install Visual Studio (e.g. Community) choosing the following `single components`:
    - `.NET SDK`
    - `.NET 7.0 Runtime`
    - `MSVC v142 - C++ Build tools`
    - `MSVC v143 - C++ Build tools`
    - `C++/CLI for Build Tools v142`
    - `Windows 10 SDK (10.2.20348.0)`
    - `CMake C++ Tools for Windows`
    - `Git for Windows`

2. Download the repository

3. The software below must be installed through Windows installers:

    - [boost_1_69_0-msvc-14.1-64.exe](https://sourceforge.net/projects/boost/files/boost-binaries/1.69.0/)

      Note: Install in `AMD\External\boost_1_69_0` folder.

3. Add the following environmental variables as User:
    - `VS_DIR`: Path to Visual Studio with all folders (e.g. `Microsoft Visual Studio\2022\Community`)

5. Run the following in a `Command Prompt` (no PowerShell) to install it:
    ```
      "%VS_DIR%\VC\Auxiliary\Build\vcvarsall.bat" x86_amd64
      cd AMD/windows-setup
      dotnet build
      dotnet run
    ```
6. Restart the command prompt and run the executable: `amd.exe`.

## Preprocessing the dataset

Before running AMD on a dataset, all input crystals must have the full motif explicitely listed in each file.
First you need to preprocess the dataset with the python script below.

Install the required python libraries:

- [CCDC library](https://www.ccdc.cam.ac.uk/solutions/csd-core/components/csd-python-api/) (needed subscription to install the entire software)
- [gemmi](https://gemmi.readthedocs.io/en/latest/install.html)

Then run the script to output the full motif cif files: 

- python.exe Script\update_cif_to_fullmotif.py [Input folder] [Output folder]

## Run AMDs computations

Finally, the dataset with the full motif can be used to run the AMD computation.
Running the executable 'amd.exe' without parameters outputs the usage as follows: 

Required options:

- -inputdir [Input Folder with (full motif) CIF files]
- -outputdir [Output Folder for Results]

Optional:

- -mindistances [Number of Averaged Minimum Distances, default: 200]

```
amd.exe -inputdir "path\to\cif_folder" -outputdir "path\to\output_dir"
```