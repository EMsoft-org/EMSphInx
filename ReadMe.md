![EMSphInx Logo](icons/sphinx.png)

# EMSphInx
*EMSphInx* is a collection of modules and programs that can be used to index a variety of diffraction patterns, ranging from EBSD to TKD, ECP, and x-ray Laue patterns. EMSphInx is currently a public beta; please report bugs to the issue tracker. If you use these programs for research, consider citing the corresponding papers:

* [EBSD Indexing](https://doi.org/10.1016/j.ultramic.2019.112841)
* [Pseudo-symmetry Prediction](https://doi.org/10.1107/S1600576719011233)

## Financial Support 
The *EMSphInx* code was developed with support from an ONR Vannevar Bush Faculty Fellowship grant, N00014-­16-­1-­2821. The central indexing algorithm is covered by a provisional patent application.

## Master Pattern Database
A [database of master patterns](https://github.com/EMsoft-org/SHTdatabase) in the new SHT file format is now available. New EMSphInx releases (0.2+) use the same file version as the database (1.1) but EMSphInx 0.1 uses an older format and is incompatable.

## Build Instructions
Nightly builds will be available soon for a variety of operating systems. Binaries are also available as [assets for tagged releases](https://github.com/EMsoft-org/EMSphInx/releases).

*EMSphInx* requires [CMake 3.14 or higher](https://www.cmake.org/download) to build. All dependencies are downloaded and compiled as part of the build process by default. The easiest way to build a non-default version of *EMSphInx* is with the cmake gui or ccmake. If you are restricted to the command line and only need the default configuration you can build with the following sequence:

Download the source

> git clone https://github.com/EMsoft-org/EMSphInx

Create a build directory and move into it

> mkdir EMSphInxBuild

> cd EMSphInxBuild

Run cmake and build, if you would like to build the GUIs you can optionally set the GUI CMake flag (e.g. -DEMSPHINX_BUILD_GUIS=ON)

> cmake ../EMSphInx

> make -j

FFTW can compile SIMD instructions on some platforms even if they are not available on the current hardware. If you encounter illegal instructions at runtime try compiling with SIMD disabled (-DEMSPHINX_FFTW_SIMD=OFF). AVX2 instructions are disabled by default but can be enabled with EMSPHINX_FFTW_AVX2.

## Documentation

Detailed documentation is available for the EMSphInxEBSD GUI program on [readthedocs](https://emsphinx.readthedocs.io/en/latest/). The functionality of the commandline programs is mostly equivlanet to the GUI and instruction are printed by running the program with no arguments.

## New features in 0.2
- Bug fixes
- GUI improvements

## What's coming in future versions
- Additional diffraction modalities
- Python bindings
- A GUI for pseudo-symmetry prediction (currently available through the command line program 'MasterXCorr')

Feel free to request additional features using the repo's [issue tracker](https://github.com/EMsoft-org/EMSphInx/issues) (please be sure to tag your issue with the 'enhancement flag')

## License ##

*EMSphInx* source files are distributed under GNU General Public License v2.0 (GPL2), see the license.txt file for details.

*EMSphInx* also includes several files from BSD licensed (3-clause) projects (please refer to the individual files for details):

- include/miniz/miniz.c
