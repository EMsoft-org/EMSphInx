
EMSphInx Documentation
==================================================

EMSphInx is a collection of modules and programs that can be used to index a variety of diffraction patterns, ranging from EBSD to TKD, ECP, and x-ray Laue patterns. EMSphInx is currently a public beta; please report bugs to the issue tracker. If you use these programs for research, consider citing the corresponding papers:

    * `EBSD Indexing <https://doi.org/10.1016/j.ultramic.2019.112841>`_
    * `Pseudo-symmetry Prediction <https://doi.org/10.1107/S1600576719011233>`_


Financial Support
------------------------------------
EMSphInx was developed with support from an ONR Vannevar Bush Faculty Fellowship grant, N00014-­16-­1-­2821. The central indexing algorithm is covered by a provisional patent application.

Installation
------------

Nightly builds will be available soon for a variety of operating systems. Binaries are also available as assets for `tagged releases <https://github.com/EMsoft-org/EMSphInx/releases>`_
.

EMSphInx requires `CMake 3.14 or higher <https://www.cmake.org/download>`_ to build. All dependencies are downloaded and compiled as part of the build process by default. The easiest way to build a non-default version of EMSphInx is with the cmake gui or ccmake. If you are restricted to the command line and only need the default configuration you can build with the following sequence:

1. Download the source, create a build directory, and move into it

.. code-block:: console

   $ git clone https://github.com/EMsoft-org/EMSphInx
   $ mkdir EMSphInxBuild
   $ cd EMSphInxBuild

2. Run cmake and build


Run cmake and build. If you would like to build the GUIs you can optionally set the GUI CMake flag (EMSPHINX_BUILD_GUIS).

.. code-block:: console

   $ cmake ../EMSphInx
   $ make -j

FFTW can compile SIMD instructions on some platforms even if they are not available on the current hardware. If you encounter illegal instructions at runtime try compiling with SIMD disabled (EMSPHINX_FFTW_SIMD). AVX2 instructions are disabled by default but can be enabled (EMSPHINX_FFTW_AVX2).

For example to build the GUI with no SIMD:

.. code-block:: console

   $ cmake ../EMSphInx -DEMSPHINX_BUILD_GUIS=ON -DEMSPHINX_FFTW_SIMD=OFF
   $ make -j


Contribute
----------

- Issue Tracker: https://github.com/EMsoft-org/EMSphInx/issues
- Source Code: https://github.com/EMsoft-org/EMSphInx


License
------------

EMSphInx is distributed under a GPL2 license with several external dependencies. Please refer to the license page for details.

Table of Contents
----------

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   self
   emsphinxebsd
   license


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

