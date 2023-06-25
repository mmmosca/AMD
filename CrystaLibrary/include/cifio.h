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
#ifndef _CIF_IO_H
#define _CIF_IO_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <istream>
#include <fstream>
#include <string>
#include <sstream>

#include "gemmi/cif.hpp"
#include "gemmi/to_cif.hpp"
#include "gemmi/symmetry.hpp"

class CIFHandler {
private:

	const static std::string cifexample;

public:

	const static std::vector<std::string> cell_parameters_tags;

	static int findValidBlockIndex(gemmi::cif::Document doc);

	/*
	*	Copy the CIF Document 
	*	param1: CIF file as Document class
	*	return: gemmi::gemmi::cif::Document class 
	*/
	static gemmi::cif::Document* copyCIFrom(gemmi::cif::Document doc);
	
	template<typename T>
	static std::vector<std::string> getElementStrings(std::vector<T> v) {
		std::stringstream sout;
		sout << std::setprecision(25);
		std::vector<std::string> strings(v.size());
		for (int i = 0; i < v.size(); ++i) {
			sout << v[i];
			strings[i] = sout.str();
			sout.str(std::string());
			sout.clear();
		}
		return strings;
	}

	/*
	*	Change a pair of a CIF Document
	*	param1: Document to modify
	*	param2: vector of 6 parameters in this order { length_a, length_b, length_c, alpha, beta, gamma }
	*	return: gemmi::gemmi::cif::Document class with parameters
	*/
	static void modifyorAddCIFPairs(gemmi::cif::Document* doc, const std::vector<std::string> parameters_tags, std::vector<std::string> &parameters);

	static void addCIFPairs(gemmi::cif::Document* doc, const std::vector<std::string> parameters_tags, std::vector<std::string> &parameters);

	/*
	*	Change or add a loop of a CIF Document
	*	param1: Document to modify
	*	param2: Loop prefix. Ex: "_atom_site_"
	*	param3: Loop tags. Ex: "label" or "type_symbol"
	*	param4: Columns of new table
	*	return: gemmi::gemmi::cif::Document class with parameters
	*/
	static bool modifyOrAddCIFLoop(gemmi::cif::Document* doc, const std::string loop_prefix, const std::vector<std::string> loop_tags, const std::vector<std::vector<std::string>> loop_rows);

	/*
	*	Perturb the lattice about the quantity specified in parameters_offset
	*	param1: Document to perturb
	*	param2: vector of 6 offsets for present parameters in this order { length_a, length_b, length_c, alpha, beta, gamma }
	*	param3: Characteristic vector of parameters_offset used to mark the elements to use 
	*	return: true if it has been perturbed, false if some pair in the Document has not been found (file corrupted)
	*/
	static bool perturbTheLattice(gemmi::cif::Document* doc, std::vector<double> &parameters_offset, std::vector<unsigned short> characteristic_vector = { 1,1,1,1,1,1 });

	/*
	*	Create a simple CIF Document with cell parameters provided as input
	*	param1: vector of 6 parameters in this order { length of a, length of b, length of c, alpha, beta, gamma }
	*	return: gemmi::gemmi::cif::Document Simple class with parameters
	*/
	static gemmi::cif::Document createEmptyCIF(std::string filename);

	static std::vector<std::string> getCIFPairs(gemmi::cif::Document doc, const std::vector<std::string> parameters_tags);

	static std::vector<double> getCIFPairsValues(gemmi::cif::Document doc, const std::vector<std::string> parameters_tags);
	
	static std::vector<std::vector<double>> getCIFAtomCoordinates(gemmi::cif::Document doc);

	/*
	*	Create the skeleton of a CIF Document from an existing CIF Document class (only with tags)
	*	param1: CIF file as Document class
	*	return: gemmi::gemmi::cif::Document class with only tags
	*/
	static gemmi::cif::Document* createEmptyCIFDocumentFrom(gemmi::cif::Document doc, std::string filename);

};

#endif // !_CIF_IO_H
