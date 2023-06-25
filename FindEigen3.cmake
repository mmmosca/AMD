#[[*
 * Permission is granted to copy, distribute and/or modify the documents
 * in this directory and its subdirectories unless otherwise stated under
 * the terms of the GNU Free Documentation License, Version 1.1 or any later version 
 * published by the Free Software Foundation; with no Invariant Sections, 
 * no Front-Cover Texts and no Back-Cover Texts. A copy of the license 
 * is available at the website of the GNU Project.
 * The programs and code snippets in this directory and its subdirectories
 * are free software; you can redistribute them and/or modify it under the 
 * terms of the GNU General Public License as published by the Free Software 
 * Foundation; either version 2 of the License, or (at your option) any later
 * version. This code is distributed in the hope that it will be useful, 
 * but WITHOUT ANY WARRANTY; without even the implied warranty of 
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * Author Marco M. Mosca, email: marcomichele.mosca@gmail.com
*]]

IF(NOT DEFINED ENV{EIGEN3_ROOT})
	message(STATUS "Using local package for Eigen3 library...")
	set(PACKAGE_DIR_EIGEN "${CMAKE_CURRENT_SOURCE_DIR}/External/eigen-3.3.7")	
ELSE()
	message(STATUS "Using system environmental variable package for Eigen3 library...")
	set(PACKAGE_DIR_EIGEN "$ENV{EIGEN3_ROOT}")

ENDIF()

set(EIGEN3_INCLUDE_DIR "${PACKAGE_DIR_EIGEN}")
set(EIGEN3_SRC_DIR "${PACKAGE_DIR_EIGEN}/Eigen/src")

IF( (NOT EXISTS ${EIGEN3_INCLUDE_DIR}) OR (NOT EXISTS ${EIGEN3_SRC_DIR}) )
	message(SEND_ERROR "EIGEN3 directories do not exist, Please assure they are inside the package folder")
	return()
ELSE()
	message(STATUS "Eigen library detected at: ${PACKAGE_DIR_EIGEN}")
ENDIF()

include_directories(${EIGEN3_INCLUDE_DIR})

set (Eigen3_FOUND 1)