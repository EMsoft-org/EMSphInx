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

# high level tests
add_executable(TestDict   ${CMAKE_CURRENT_LIST_DIR}/dict.cpp      )
add_executable(TestDiag   ${CMAKE_CURRENT_LIST_DIR}/diagram.cpp   )

# tests for include/sht/*
add_executable(TestWigner ${CMAKE_CURRENT_LIST_DIR}/sht/wigner.cpp     )
add_executable(TestSquare ${CMAKE_CURRENT_LIST_DIR}/sht/square_sht.cpp )
add_executable(TestXCorr  ${CMAKE_CURRENT_LIST_DIR}/sht/sht_xcorr.cpp  )

# tests for include/util/*
add_executable(TestBase64 ${CMAKE_CURRENT_LIST_DIR}/util/base64.cpp    )
add_executable(TestLinAlg ${CMAKE_CURRENT_LIST_DIR}/util/linalg.cpp    )
add_executable(TestThread ${CMAKE_CURRENT_LIST_DIR}/util/threadpool.cpp)
add_executable(TestTimer  ${CMAKE_CURRENT_LIST_DIR}/util/timer.cpp     )
add_executable(TestColor  ${CMAKE_CURRENT_LIST_DIR}/util/colorspace.cpp)
add_executable(TestNML    ${CMAKE_CURRENT_LIST_DIR}/util/nml.cpp       )

# tests for include/xtal/*
add_executable(TestRot    ${CMAKE_CURRENT_LIST_DIR}/xtal/rotations.cpp )
add_executable(TestQuat   ${CMAKE_CURRENT_LIST_DIR}/xtal/quaternion.cpp)
add_executable(TestSym    ${CMAKE_CURRENT_LIST_DIR}/xtal/symmetry.cpp  )
add_executable(TestPos    ${CMAKE_CURRENT_LIST_DIR}/xtal/position.cpp  )
add_executable(TestHM     ${CMAKE_CURRENT_LIST_DIR}/xtal/hm.cpp        )

################################
#         dependencies         #
################################

if(${EMSPHINX_BUILD_FFTW})
	# high level tests
	add_dependencies(TestDict   ${FFTW_DEPENDS})
	add_dependencies(TestDiag   ${FFTW_DEPENDS})

	# tests for include/sht/*
	add_dependencies(TestSquare ${FFTW_DEPENDS})
	add_dependencies(TestXCorr  ${FFTW_DEPENDS})
endif()
if(${EMSPHINX_BUILD_HDF5})
	# high level tests
	add_dependencies(TestDict  hdf5)
	add_dependencies(TestDiag  hdf5)

	# tests for include/sht/*
	add_dependencies(TestXCorr hdf5)

	# tests for include/xtal/*
	add_dependencies(TestNML   hdf5)
endif()

if(${BuildMiniZ})
	add_dependencies(TestDiag miniz)
endif()

################################
#           linking            #
################################

# high level tests
target_link_libraries(TestDict   ${FFTW_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(TestDiag   ${FFTW_LIBRARIES} ${HDF5_LIBRARIES})

# include/sht/*
target_link_libraries(TestSquare ${FFTW_LIBRARIES})
target_link_libraries(TestXCorr  ${FFTW_LIBRARIES} ${HDF5_LIBRARIES})

# include/util/*
target_link_libraries(TestNML                      ${HDF5_LIBRARIES})

