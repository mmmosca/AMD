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
IF(NOT DEFINED ENV{GEMMI_ROOT})
	message(STATUS "Using local package for Gemmi library...")
	set(PACKAGE_DIR_GEMMI "${CMAKE_CURRENT_SOURCE_DIR}/External/gemmi-0.3.3")	
ELSE()
	message(STATUS "Using system environmental variable package for Gemmi library...")
	set(PACKAGE_DIR_GEMMI "$ENV{GEMMI_ROOT}")

ENDIF()

set(GEMMI_INCLUDE_DIR "${PACKAGE_DIR_GEMMI}/include")
set(3RD_INCLUDE_DIR "${PACKAGE_DIR_GEMMI}/third_party" )
set(GEMMI_SRC_DIR "${PACKAGE_DIR_GEMMI}/src")

IF( (NOT EXISTS ${GEMMI_INCLUDE_DIR}) OR (NOT EXISTS ${GEMMI_SRC_DIR}) )
	message(SEND_ERROR "GEMMI directories do not exist, Please assure they are inside the package folder")
	return()
ELSE()
	message(STATUS "Gemmi library detected at: ${PACKAGE_DIR_GEMMI}")
ENDIF()

include_directories(${GEMMI_INCLUDE_DIR} ${3RD_INCLUDE_DIR} )

set (Gemmi_FOUND 1)