/*
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
*/
#ifndef _STRINGHANDLER_H
#define _STRINGHANDLER_H

#include <string>
#include <vector>
#include <set>
#include <boost/algorithm/string.hpp>

std::vector<std::string> splitSingleAtomQuery(std::string allquery, std::string delimiter);
std::set<std::string> splitSingleAtomQuery_Set(std::string allquery, std::string delimiter);
std::vector<std::set<std::string>> splitAtomQuery_Set(std::string allquery);
std::set<std::string> splitAtoms(std::string allquery);
#endif //! _STRINGHANDLER_H