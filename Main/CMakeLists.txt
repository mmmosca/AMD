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
project(main)

add_executable(amd compute_amds.cpp)

add_dependencies(amd CrystaLibrary)

target_link_libraries(amd CrystaLibrary)

set_target_properties(amd PROPERTIES 
	LIBRARY_OUTPUT_DIRECTORY "${WORKDIR}"
	RUNTIME_OUTPUT_DIRECTORY "${WORKDIR}"
	ARCHIVE_OUTPUT_DIRECTORY "${WORKDIR}")

install(TARGETS amd DESTINATION bin)