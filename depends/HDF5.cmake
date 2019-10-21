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

################################################
#                                              #
#                     HDF5                     #
#                                              #
################################################

set(HDF5_LIBRARIES "") # store all hdf5 libraries to link against here

if(EMSPHINX_BUILD_HDF5)
	set(HDF5_VERS "hdf5-1_8_20") # version of hdf5 to build, 1_10_x requires cmake 3.10 or newer
	set(HDF5_URL "https://bitbucket.hdfgroup.org/scm/hdffv/hdf5.git") # git repo to fetch hdf5 source from
	set(HDF5_OPTIONS -DBUILD_TESTING=OFF -DHDF5_BUILD_EXAMPLES=OFF -DHDF5_BUILD_TOOLS=OFF -DHDF5_DISABLE_COMPILER_WARNINGS=ON -DHDF5_BUILD_HL_LIB=OFF) # build options (arguments for ccmake)
	set(HDF5_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/hdf5) # location to build

	if(EMSPHINX_BUILD_SHARED)
		set(HDF5_OPTIONS ${HDF5_OPTIONS} -DBUILD_SHARED_LIBS=ON ) # build options (arguments for ccmake)
	else()
		set(HDF5_OPTIONS ${HDF5_OPTIONS} -DBUILD_SHARED_LIBS=OFF) # build options (arguments for ccmake)
	endif()

	# now set up hdf5 build
	include(ExternalProject)
	ExternalProject_add(hdf5 PREFIX hdf5 GIT_REPOSITORY ${HDF5_URL} GIT_TAG ${HDF5_VERS} GIT_SHALLOW TRUE
		# UPDATE_DISCONNECTED 1 # this should keep the cmake from trying to repull from the repo every time make is run but doesn't...
		UPDATE_COMMAND "" # this does keep cmake from redoing everythin on every make but means that if you change HDF5_VERS it won't automatically pull and build the new one
		CMAKE_ARGS ${HDF5_OPTIONS} -DCMAKE_INSTALL_PREFIX=${HDF5_BUILD_DIR}/install
		BUILD_COMMAND ${CMAKE_COMMAND} --build . --parallel ${NCORES}
	)
	if(EMSPHINX_BUILD_SHARED) # copy shared library to binary directory if needed
		ExternalProject_Add_Step(hdf5 CopyLibs # execute copy command as part of build
			COMMAND ${CMAKE_COMMAND} -E copy_directory ${HDF5_BUILD_DIR}/install/lib/ ${CMAKE_CURRENT_BINARY_DIR} # hdf5 seems to use /lib/ even when /lib64/ is preferred
			COMMENT "copying HDF5 binarys to build folder"
			DEPENDEES install # don't copy from install dir until files have been copied into it
		)
	endif()

	# hdf5 names librarys [lib]hdf5_x[-static].(a/lib) or [lib]hdf5_x[-shared].(so/dll)
	if(MSVC)
		set(HDF5_PRE  "lib") # windows always uses libhdf5(_x).lib
		set(HDF5_POST ""   ) # windows doesn't use a library suffix for hdf5
		set(HDF5_DEB  "_D" ) # windows appends debug library names with _D 
	else()
		set(HDF5_PRE ${LIB_PRE})
		if(EMSPHINX_BUILD_SHARED)
			set(HDF5_POST "-shared")
		else()
			set(HDF5_POST "-static")
		endif()
		set(HDF5_DEB "")
	endif()

	# finally accumulate hdf5 libraries and headers
	include_directories(${HDF5_BUILD_DIR}/install/include)
	set(HDF5_NAMES hdf5_cpp hdf5) # hdf5_hl hdf5_hl_cpp
	# set(HDF5_NAMES hdf5 hdf5_hl hdf5_cpp hdf5_hl_cpp)
	foreach(HDF5_LIB ${HDF5_NAMES})
		list(APPEND HDF5_LIBRARIES ${HDF5_BUILD_DIR}/install/lib/${HDF5_PRE}${HDF5_LIB}${HDF5_POST}$<$<CONFIG:Debug>:${HDF5_DEB}>${LNK_EXT}) # hdf5 seems to use /lib/ even when /lib64/ is preferred
	endforeach()

	if(UNIX AND NOT APPLE)
        list(APPEND HDF5_LIBRARIES ${CMAKE_DL_LIBS}) # this isn't needed for (but doesn't disrupt) debian, is required for ubuntu
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ldl")
	endif()

else(EMSPHINX_BUILD_HDF5)
	find_package(HDF5 REQUIRED COMPONENTS CXX PATHS "${HDF5_ROOT}") # search for HDF5
	if(HDF5_FOUND) # hdf5 was found automatically
		include_directories(${HDF5_INCLUDE_DIR}) # include hdf5 headers

		#should be able to just use find package but it seems to be broken so I'll manually find libraries instead
		if("${HDF5_LIBRARIES}" STREQUAL "") # did find_package actually get HDF5_LIBRARIES?
			get_filename_component(HDF5_ROOT "${HDF5_INCLUDE_DIR}" DIRECTORY CACHE PATH) # parent of include directory
			find_library(HDF5_C_LIBRARIES hdf5-shared PATH "${HDF5_ROOT}/lib") # prefer shared libraries
			find_library(HDF5_CXX_LIBRARIES hdf5_cpp-shared PATH "${HDF5_ROOT}/lib") # prefer shared libraries
			if(NOT ${HDF5_C_LIBRARIES})
				find_library(HDF5_C_LIBRARIES hdf5-static PATH "${HDF5_ROOT}/lib") # fall back to static libraries
			endif()
			if(NOT ${HDF5_CXX_LIBRARIES})
				find_library(HDF5_CXX_LIBRARIES hdf5_cpp-static PATH "${HDF5_ROOT}/lib") # fall back to static libraries
			endif()
			set(HDF5_LIBRARIES ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES})
		endif()
	else(HDF5_FOUND) # hdf5 wasn't found automatically
		message(FATAL_ERROR "This project requires HDF5 configured with CMake, try setting the HDF5_ROOT variable")
	endif(HDF5_FOUND)
endif(EMSPHINX_BUILD_HDF5)
