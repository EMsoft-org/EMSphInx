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
#                     FFTW                     #
#                                              #
################################################

# apply switches for which versions of fftw to use
if(NOT (${EMSPHINX_FFTW_F} OR ${EMSPHINX_FFTW_D} OR ${EMSPHINX_FFTW_L}) )
	message(FATAL_ERROR "must link against at least one of fftw_f, fftw, or fftw_l (set EMSPHINX_FFTW_{F,D, or L} to ON")
endif()

# define preprocessor macros to define function that use float/double/long double versions of fftw
if(${EMSPHINX_FFTW_F})
	add_definitions(-DEM_USE_F) # if cmake is 3.11+ should use add_compile_definitions(EM_USE_F)
endif()
if(${EMSPHINX_FFTW_D})
	add_definitions(-DEM_USE_D) # if cmake is 3.11+ should use add_compile_definitions(EM_USE_D)
endif()
if(${EMSPHINX_FFTW_L})
	add_definitions(-DEM_USE_L) # if cmake is 3.11+ should use add_compile_definitions(EM_USE_L)
endif()

set(FFTW_LIBRARIES "") # store all fftw libraries to link against here
set(FFTW_DEPENDS "") # store all fftw dependencies here

if(EMSPHINX_BUILD_FFTW) # download + build fftw
	# version of fftw to build, any version 3+ should be fine
	set(FFTW_VER "3.3.7")
	
	# name to use for projects (folder name for build + internal cmake variable name)
	set(FFTW_NAME_F fftwf)
	set(FFTW_NAME_D fftwd)
	set(FFTW_NAME_L fftwl)

	# URL to fetch FFTW source from
	set(FFTW_URL "http://www.fftw.org/fftw-${FFTW_VER}.tar.gz")

	# name of file to save fetched source as
	set(FFTW_SAVE "fftw-${FFTW_VER}.tar.gz")
	
	# common ccmake args for fftw versions
	set(FFTW_OPTIONS -DBUILD_TESTS=OFF -DCMAKE_INSTALL_LIBDIR=lib) # force install to /lib/ (fftw defaults to /lib64/ on some systems)
	if(EMSPHINX_BUILD_SHARED)
		set(FFTW_OPTIONS ${FFTW_OPTIONS} -DBUILD_SHARED_LIBS=ON ) # build options (arguments for ccmake)
	else()
		set(FFTW_OPTIONS ${FFTW_OPTIONS} -DBUILD_SHARED_LIBS=OFF) # build options (arguments for ccmake)
	endif()

	# simd ccmake args for fftw versions (only float + double)
	if(EMSPHINX_FFTW_SIMD)
		if(EMSPHINX_FFTW_AVX2)
			set(FFTW_SIMD -DENABLE_SSE=ON -DENABLE_SSE2=ON -DENABLE_AVX=ON -DENABLE_AVX2=ON ) # flags to enable SIMD instructions for fftw
		else()
			set(FFTW_SIMD -DENABLE_SSE=ON -DENABLE_SSE2=ON -DENABLE_AVX=ON -DENABLE_AVX2=OFF) # flags to enable SIMD instructions for fftw
		endif()
	else()
		set(FFTW_SIMD "") # flags to enable SIMD instructions for fftw
	endif()

	# location to install builds
	set(FFTW_BUILD_DIR_F ${CMAKE_CURRENT_BINARY_DIR}/${FFTW_NAME_F})
	set(FFTW_BUILD_DIR_D ${CMAKE_CURRENT_BINARY_DIR}/${FFTW_NAME_D})
	set(FFTW_BUILD_DIR_L ${CMAKE_CURRENT_BINARY_DIR}/${FFTW_NAME_L})

	# now set up each version of fftw required
	include(ExternalProject)

	# build float fftw if needed
	if(${EMSPHINX_FFTW_F})
		ExternalProject_add(${FFTW_NAME_F} PREFIX ${FFTW_NAME_F} URL ${FFTW_URL} DOWNLOAD_NAME ${FFTW_SAVE}
			CMAKE_ARGS ${FFTW_OPTIONS} ${FFTW_SIMD} -DCMAKE_INSTALL_PREFIX=${FFTW_BUILD_DIR_F}/install -DENABLE_FLOAT=ON # agruments for ccmake
			BUILD_COMMAND ${CMAKE_COMMAND} --build . --parallel ${NCORES}
		)
		if(EMSPHINX_BUILD_SHARED) # copy shared library to binary directory if needed
			ExternalProject_Add_Step(${FFTW_NAME_F} CopyLibs # execute copy command as part of build
				COMMAND ${CMAKE_COMMAND} -E copy_directory ${FFTW_BUILD_DIR_F}/install/${CPY_DIR}/ ${CMAKE_CURRENT_BINARY_DIR}
				COMMENT "copying fftwf binarys to build folder"
				DEPENDEES install # don't copy from install dir until files have been copied into it
			)
		endif()
	endif(${EMSPHINX_FFTW_F})

	# build double fftw if needed
	if(${EMSPHINX_FFTW_D})
		ExternalProject_add(${FFTW_NAME_D} PREFIX ${FFTW_NAME_D} URL ${FFTW_URL} DOWNLOAD_NAME ${FFTW_SAVE}
			CMAKE_ARGS ${FFTW_OPTIONS} ${FFTW_SIMD} -DCMAKE_INSTALL_PREFIX=${FFTW_BUILD_DIR_D}/install # agruments for ccmake
			BUILD_COMMAND ${CMAKE_COMMAND} --build . --parallel ${NCORES}
		)
		if(EMSPHINX_BUILD_SHARED) # copy shared library to binary directory if needed
			ExternalProject_Add_Step(${FFTW_NAME_D} CopyLibs # execute copy command as part of build
				COMMAND ${CMAKE_COMMAND} -E copy_directory ${FFTW_BUILD_DIR_D}/install/${CPY_DIR}/ ${CMAKE_CURRENT_BINARY_DIR}
				COMMENT "copying fftw binarys to build folder"
				DEPENDEES install # don't copy from install dir until files have been copied into it
			)
		endif()
	endif(${EMSPHINX_FFTW_D})

	# build long double fftw if needed
	if(${EMSPHINX_FFTW_L})
		ExternalProject_add(${FFTW_NAME_L} PREFIX ${FFTW_NAME_L} URL ${FFTW_URL} DOWNLOAD_NAME ${FFTW_SAVE}
			CMAKE_ARGS ${FFTW_OPTIONS} -DCMAKE_INSTALL_PREFIX=${FFTW_BUILD_DIR_L}/install -DENABLE_LONG_DOUBLE=ON # agruments for ccmake
			BUILD_COMMAND ${CMAKE_COMMAND} --build . --parallel ${NCORES}
		)
		if(EMSPHINX_BUILD_SHARED) # copy shared library to binary directory if needed
			ExternalProject_Add_Step(${FFTW_NAME_L} CopyLibs # execute copy command as part of build
				COMMAND ${CMAKE_COMMAND} -E copy_directory ${FFTW_BUILD_DIR_L}/install/${CPY_DIR}/ ${CMAKE_CURRENT_BINARY_DIR}
				COMMENT "copying fftwl binarys to build folder"
				DEPENDEES install # don't copy from install dir until files have been copied into it
			)
		endif()
	endif(${EMSPHINX_FFTW_L})

	# include fftw3.h once
	if(${EMSPHINX_FFTW_F})
		include_directories(${FFTW_BUILD_DIR_F}/install/include)
	elseif(${EMSPHINX_FFTW_D})
		include_directories(${FFTW_BUILD_DIR_D}/install/include)
	elseif(${EMSPHINX_FFTW_L})
		include_directories(${FFTW_BUILD_DIR_L}/install/include)
	endif()

	# build list of fftw libraries / projects
	if(${EMSPHINX_FFTW_F})
		list(APPEND FFTW_LIBRARIES ${FFTW_BUILD_DIR_F}/install/lib/${LIB_PRE}fftw3f${LNK_EXT}) # add built library to list for linking
		list(APPEND FFTW_DEPENDS fftwf) # accumulate dependencies
	endif()
	if(${EMSPHINX_FFTW_D})
		list(APPEND FFTW_LIBRARIES ${FFTW_BUILD_DIR_D}/install/lib/${LIB_PRE}fftw3${LNK_EXT} ) # add built library to list for linking
		list(APPEND FFTW_DEPENDS fftwd) # accumulate dependencies
	endif()
	if(${EMSPHINX_FFTW_L})
		list(APPEND FFTW_LIBRARIES ${FFTW_BUILD_DIR_L}/install/lib/${LIB_PRE}fftw3l${LNK_EXT}) # add built library to list for linking
		list(APPEND FFTW_DEPENDS fftwl) # accumulate dependencies
	endif()
else(EMSPHINX_BUILD_FFTW) # use existing fftw builds
	# find float libraries if needed
	if(${EMSPHINX_FFTW_F})
		find_library(FFTW_LIBRARY_F fftw3f "fftwfloat library")
		list(APPEND FFTW_LIBRARIES ${FFTW_LIBRARY_F})
	endif()

	# find double libraries if needed
	if(${EMSPHINX_FFTW_D})
		find_library(FFTW_LIBRARY_D fftw3 "fftw double library")
		list(APPEND FFTW_LIBRARIES ${FFTW_LIBRARY_D})
	endif()

	# find long double libraries if needed
	if(${EMSPHINX_FFTW_L})
		find_library(FFTW_LIBRARY_L fftw3l "fftw long double library")
		list(APPEND FFTW_LIBRARIES ${FFTW_LIBRARY_L})
	endif()

	# find the header
	find_file(FFTW_HEADER fftw.h "fftw header")
	get_filename_component(FFTW_INCLUDE ${FFTW_HEADER} DIRECTORY)
	include_directories(${FFTW_INCLUDE})
endif(EMSPHINX_BUILD_FFTW)
