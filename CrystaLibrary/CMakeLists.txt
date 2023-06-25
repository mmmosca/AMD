#[[
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
cmake_minimum_required (VERSION 3.8)
project(CrystaLibrary)

set( CURR_DIR "${CMAKE_CURRENT_SOURCE_DIR}" )
set( INCLUDE_DIR "${CURR_DIR}/include" )
set( SRC_DIR "${CURR_DIR}/src" )

# Aggiungere l'origine all'eseguibile di questo progetto.

IF(WIN32)
	# Add 'lib' at the beginning of the file
	# set(CMAKE_SHARED_LIBRARY_PREFIX "lib")

	# Export symbols in a .lib library file
	set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
ENDIF()

file(GLOB CFILES ${SRC_DIR}/*.cpp)
file(GLOB HFILES ${INCLUDE_DIR}/*.h)

if (SHARED)
	add_library(${PROJECT_NAME} SHARED ${CFILES} ${HFILES})
else()
	add_library(${PROJECT_NAME} STATIC ${CFILES} ${HFILES})
ENDIF()

target_include_directories(${PROJECT_NAME} PUBLIC ${INCLUDE_DIR} )
target_link_libraries(${PROJECT_NAME} ${BABEL_LIBRARIES})

set_target_properties(${PROJECT_NAME} PROPERTIES 
	LIBRARY_OUTPUT_DIRECTORY "${WORKDIR}"
	RUNTIME_OUTPUT_DIRECTORY "${WORKDIR}"
	ARCHIVE_OUTPUT_DIRECTORY "${WORKDIR}"
)

install(TARGETS ${PROJECT_NAME} DESTINATION lib)
install(FILES ${HFILES} DESTINATION include)
