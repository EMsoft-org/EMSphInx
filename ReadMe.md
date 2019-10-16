![EMSphInx Logo](icons/sphinx.png)

# EMSphInx
*EMSphInx* is a collection of modules and programs that can be used to index a variety of diffraction patterns, ranging from EBSD to TKD, ECP, and x-ray Laue patterns. EMSphInx is currently a public beta; please report bugs to the issue tracker. If you use these programs for research, consider citing the corresponding papers:

* [EBSD Indexing](https://doi.org/10.1016/j.ultramic.2019.112841)
* [Pseudo-symmetry Prediction](https://doi.org/10.1107/S1600576719011233)

## Financial Support 
The *EMSphInx* code was developed with support from an ONR Vannevar Bush Faculty Fellowship grant, N00014-­16-­1-­2821. The central indexing algorithm is covered by a provisional patent application.

## Build Instructions
Nightly builds will be available soon for a variety of operating systems

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

## Utility Program Overview

1. IndexEBSD

   index EBSD patterns on the command line with a namelist file

2. MasterXcorr

   compute spherical cross correlation between 2 spherical master patterns

3. ShtWisdom

   build FFTW wisdom on new systems (reduces initialization time on first execution of other programs)

4. mp2sht

   convert from EMsoft EBSD master patterns to the new SHT file format used by the indexing programs

5. EMSphInxEBSD (only if EMSPHINX_BUILD_GUIS=ON)

   graphical user interface to build namelist files for IndexEBSD and/or index patterns directly

Help files for these programs are available as wiki pages on [github.com:EMsoft-org/EMSphInx/wiki]() or in the documentation folder of this repository.

## New features in 0.1
- Public Beta

## What's coming in future versions
- Additional diffraction modalities
- Python bindings

Feel free to request additional features using the repo's [issue tracker](https://github.com/EMsoft-org/EMSphInx/issues) (please be sure to tag your issue with the 'enhancement flag')

## License ##

*EMSphInx* source files are distributed under GNU General Public License v2.0 (GPL2), see the license.txt file for details.

*EMSphInx* also includes several files from BSD licensed (3-clause) projects (please refer to the individual files for details):

- include/miniz/miniz.c
