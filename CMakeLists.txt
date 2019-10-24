# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#                                                                     #
# Copyright (c) 2019, De Graef Group, Carnegie Mellon University      #
# All rights reserved.                                                #
#                                                                     #
# Author: William C. Lenthe                                           #
#                                                                     #
# This package is free software; you can redistribute it and/or       #
# modify it under the terms of the GNU General Public License as      #
# published by the Free Software Foundation; either version 2 of the  #
# License, or (at your option) any later version.                     #
#                                                                     #
# This program is distributed in the hope that it will be useful,     #
# but WITHOUT ANY WARRANTY; without even the implied warranty of      #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       #
# GNU General Public License for more details.                        #
#                                                                     #
# You should have received a copy of the GNU General Public License   #
# along with this program; if not, check the Free Software Foundation #
# website: <https://www.gnu.org/licenses/old-licenses/gpl-2.0.html>   #
#                                                                     #
#                                                                     #
# Interested in a commercial license? Contact:                        #
#                                                                     #
# Center for Technology Transfer and Enterprise Creation              #
# 4615 Forbes Avenue, Suite 302                                       #
# Pittsburgh, PA 15213                                                #
#                                                                     #
# phone. : 412.268.7393                                               #
# email  : innovation@cmu.edu                                         #
# website: https://www.cmu.edu/cttec/                                 #
#                                                                     #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

project(EMSphInx)
# cmake_minimum_required(VERSION 3.12) # 3.12 for the parallel option for cmake --build
cmake_minimum_required(VERSION 3.14) # for FetchContent_MakeAvailable (could fall back to 3.12 with alternate pattern but it is now deprecated)

set(CMAKE_CXX_STANDARD 11) # apparently you still need to explicitly request 2011 conformance in 2019...

################################
#       package options        #
################################

option(EMSPHINX_BUILD_SHARED "should shared or static libraries be preferred"                                           OFF)
option(EMSPHINX_FFTW_F       "enable float templates by linking against fftw_f"                                         OFF)
option(EMSPHINX_FFTW_D       "enable double templates by linking against fftw"                                          ON )
option(EMSPHINX_FFTW_L       "enable long double templates by linking against fftw_l"                                   OFF)
option(EMSPHINX_BUILD_FFTW   "fetch and build FFTW instead of using an existing system installation"                    ON )
option(EMSPHINX_FFTW_SIMD    "should SSE/SSE2/AVX instructions be enabled for fftw build (faster run, slower compiler)" ON )
option(EMSPHINX_FFTW_AVX2    "should AVX2 instructions be enabled for fftw build (not supported on all processors)    " OFF)
option(EMSPHINX_BUILD_HDF5   "fetch and build HDF5 instead of using an existing system installation"               ON )
option(EMSPHINX_BUILD_TESTS  "should test programs (/test/*) be built"                                             ON )
option(EMSPHINX_BUILD_GUIS   "should the GUI programs be built (experimental)"                                     OFF)
if(EMSPHINX_BUILD_GUIS)
	option(EMSPHINX_BUILD_wxWidgets "should wxWidgets be built from source" ON)
endif()

################################
#     multicore compiling      #
################################

include(ProcessorCount)
ProcessorCount(NCORES)
if(NOT NCORES EQUAL 0)
	set(CMAKE_BUILD_PARALLEL_LEVEL NCORES)
endif()

################################
# compiler specific adjustment #
################################

# visual studio
if(MSVC)
	add_definitions(-D_SCL_SECURE_NO_WARNINGS) # these warnings are annoying and not critical
	add_definitions(-D_CRT_SECURE_NO_WARNINGS) # these warnings are annoying and not critical
	# add_compile_options($<$<CXX_COMPILER_ID:MSVC>:/MP>) # use multi core compiling
	add_compile_options(/MP) # the generator expression version breaks wxWidgets
endif()

# gcc
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-pthread" )# std::thread requires pthreads on gcc
endif()

# non clang osx
# if(APPLE AND NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang") # only for gcc, not clang
	# set(CMAKE_INSTALL_RPATH_USE_LINK_PATH OFF) # there seems to be an rpath issue on apple + GCC
# endif()
set(MACOSX_RPATH ON) # still haven't fully resolved this...

################################
#   dynamic library copying    #
################################

# get prefix/extension for built libraries
if(EMSPHINX_BUILD_SHARED)
	set(LIB_PRE ${CMAKE_SHARED_LIBRARY_PREFIX})
	set(LIB_EXT ${CMAKE_SHARED_LIBRARY_SUFFIX})
else()
	set(LIB_PRE ${CMAKE_STATIC_LIBRARY_PREFIX})
	set(LIB_EXT ${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

# determine name of files to link against (windows makes .dll for runtime but also keeps .lib for linking) and copy into build directory
if(MSVC)
	set(LNK_EXT ${CMAKE_STATIC_LIBRARY_SUFFIX}) # link against the .lib
	set(CPY_DIR bin) # location to copy shared libraries from
else()
	set(LNK_EXT ${LIB_EXT}) # link against static or dynamic as appropriate
	set(CPY_DIR lib) # location to copy shared libraries from
endif()

################################
#        git versioning        #
################################

find_package(Git REQUIRED)
execute_process( # get branch name
	COMMAND ${GIT_EXECUTABLE} rev-parse --abbrev-ref HEAD
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	OUTPUT_VARIABLE EMSPHINX_GIT_BRANCH
	OUTPUT_STRIP_TRAILING_WHITESPACE
)
execute_process( # get abbreviated hash
	COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
	WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
	OUTPUT_VARIABLE EMSPHINX_GIT_HASH
	OUTPUT_STRIP_TRAILING_WHITESPACE
)

if("${EMSPHINX_GIT_BRANCH}" STREQUAL "")
	message(FATAL_ERROR "failed to determine git branch)")
endif()
if("${EMSPHINX_GIT_HASH}" STREQUAL "")
	message(FATAL_ERROR "failed to determine git hash)")
endif()

add_definitions("-DEMSPHINX_GIT_HASH=${EMSPHINX_GIT_HASH}")
add_definitions("-DEMSPHINX_GIT_BRANCH=${EMSPHINX_GIT_BRANCH}")

################################
#         dependencies         #
################################

# fetch SHT file format content
# this could be done with a git submodule but the user would have to call:
#   git submodule update
#   git clone --recursive
# which is a bit annoying
include(FetchContent)
FetchContent_Declare(
	SHTfile
	GIT_REPOSITORY "https://github.com/EMsoft-org/SHTfile"
	GIT_TAG        "e49ad6b"
	# GIT_TAG        "v3.1.2" # just get a recent version for now
	# GIT_PROGRESS   TRUE # its currently only 1 file, we probably don't need to print out progress
)
# set(FETCHCONTENT_QUIET NO) # again only 1 file
FetchContent_MakeAvailable(SHTfile)
FetchContent_GetProperties(SHTfile BINARY_DIR SHTfile_BINARY_DIR)
include_directories(${SHTfile_BINARY_DIR})

# compiled dependencies
include(${CMAKE_CURRENT_LIST_DIR}/depends/FFTW.cmake) # build or find existing FFTW
include(${CMAKE_CURRENT_LIST_DIR}/depends/HDF5.cmake) # build or find existing HDF5

if(EMSPHINX_BUILD_GUIS)
	include(${CMAKE_CURRENT_LIST_DIR}/depends/wxWidgets.cmake) # build or find existing HDF5
endif()

# png support, for now I've included these files directly due to download issues on some versions of cmake (if libcurl doesn't support https)
set(BuildMiniZ OFF)
if(BuildMiniZ)
	set(MINIZ_VERSION "2.0.8")
	ExternalProject_Add(
		miniz PREFIX "miniz"
		URL https://github.com/richgel999/miniz/releases/download/${MINIZ_VERSION}/miniz-${MINIZ_VERSION}.zip
		URL_HASH MD5=0692c3f080267b24419ab67c1eefc881 # make sure the download wasn't corrupted (/ someone hasn't injected a different version)
		CONFIGURE_COMMAND "" ${CMAKE_COMMAND} -E copy_directory ${CMAKE_BINARY_DIR}/miniz/src/miniz/ ${CMAKE_SOURCE_DIR}/include/miniz # copy everything (source + license are most critical)
		BUILD_COMMAND "" INSTALL_COMMAND ""
	)
endif(BuildMiniZ)

################################
#     add actual components    #
################################

# include headers
include_directories(${CMAKE_CURRENT_LIST_DIR}/include)
include_directories(${CMAKE_CURRENT_LIST_DIR}/icons)

# add test programs if needed
if(EMSPHINX_BUILD_TESTS)
	add_subdirectory(${CMAKE_CURRENT_LIST_DIR}/test ${CMAKE_BINARY_DIR}/test) # don't put in /cpp/test
endif()

# add actual programs
include(${CMAKE_CURRENT_LIST_DIR}/programs/CMakeLists.txt) # include instead of add_subdirectory so we don't get a subfolder in the build
