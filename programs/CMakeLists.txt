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

################################
#        add executables       #
################################

add_executable(IndexEBSD     ${CMAKE_CURRENT_LIST_DIR}/index_ebsd.cpp    )
add_executable(MasterXcorr   ${CMAKE_CURRENT_LIST_DIR}/master_xcorr.cpp  )
add_executable(PatternRepack ${CMAKE_CURRENT_LIST_DIR}/pattern_repack.cpp)
add_executable(EBSPDims      ${CMAKE_CURRENT_LIST_DIR}/ebsp_dims.cpp     )
add_executable(mp2sht        ${CMAKE_CURRENT_LIST_DIR}/mp2sht.cpp        )
add_executable(sht2png       ${CMAKE_CURRENT_LIST_DIR}/sht2png.cpp       )
add_executable(ShtWisdom     ${CMAKE_CURRENT_LIST_DIR}/sht_wisdom.cpp    )

install(TARGETS IndexEBSD MasterXcorr mp2sht sht2png ShtWisdom COMPONENT clt DESTINATION bin) # command line tools

if(${EMSPHINX_BUILD_GUIS})

	set(MACOSX_BUNDLE_ICON_FILE sphinx.icns) # add icon to Info.plist
	set(SPHINX_ICON ${CMAKE_SOURCE_DIR}/icons/sphinx.icns) # location of icon file in source
	set_source_files_properties(${SPHINX_ICON} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources") # copy the icon file to the Resources folder of the bundle

	if(WIN32)
		set(SPHINX_ICON ${CMAKE_SOURCE_DIR}/icons/sphinx.rc) # location of icon file in source
	endif()

	add_executable(EMSphInxEBSD WIN32 MACOSX_BUNDLE ${CMAKE_CURRENT_LIST_DIR}/ebsd_wizard.cpp ${SPHINX_ICON})
	set_target_properties(EMSphInxEBSD PROPERTIES
        MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_LIST_DIR}/MacOSXBundleInfo.plist.in)
	install(TARGETS EMSphInxEBSD COMPONENT gui BUNDLE DESTINATION .) # command line tools
endif()

################################
#         dependencies         #
################################

if(${EMSPHINX_BUILD_FFTW})
	add_dependencies(IndexEBSD    ${FFTW_DEPENDS})
	add_dependencies(MasterXcorr  ${FFTW_DEPENDS})
	add_dependencies(mp2sht       ${FFTW_DEPENDS})
	add_dependencies(sht2png      ${FFTW_DEPENDS})
	add_dependencies(ShtWisdom    ${FFTW_DEPENDS})
endif()
if(${EMSPHINX_BUILD_HDF5})
	add_dependencies(IndexEBSD     hdf5)
	add_dependencies(MasterXcorr   hdf5)
	add_dependencies(PatternRepack hdf5)
	add_dependencies(EBSPDims      hdf5)
	add_dependencies(mp2sht        hdf5)
	add_dependencies(sht2png       hdf5)
endif()
if(${EMSPHINX_BUILD_GUIS})
	if(${EMSPHINX_BUILD_FFTW})
		add_dependencies(EMSphInxEBSD ${FFTW_DEPENDS})
	endif()
	if(${EMSPHINX_BUILD_HDF5})
		add_dependencies(EMSphInxEBSD hdf5)
	endif()
endif()

################################
#           linking            #
################################

# link
target_link_libraries(IndexEBSD        ${FFTW_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(MasterXcorr      ${FFTW_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(PatternRepack                      ${HDF5_LIBRARIES})
target_link_libraries(EBSPDims                           ${HDF5_LIBRARIES})
target_link_libraries(mp2sht           ${FFTW_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(sht2png          ${FFTW_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(ShtWisdom        ${FFTW_LIBRARIES}                  )

if(${EMSPHINX_BUILD_GUIS})
	target_link_libraries(EMSphInxEBSD ${FFTW_LIBRARIES} ${HDF5_LIBRARIES} core base propgrid html)
endif()
