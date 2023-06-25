
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
#include <stringhandler.h>


std::vector<std::string> splitSingleAtomQuery(std::string allquery, std::string delimiter) {
	std::vector<std::string> query_atoms;
	boost::algorithm::split(query_atoms, allquery, boost::is_any_of(delimiter), boost::token_compress_on);

	for (int i = 0; i < query_atoms.size(); ++i) {
		boost::trim(query_atoms[i]);
	}

	return query_atoms;
}

std::set<std::string> splitSingleAtomQuery_Set(std::string allquery, std::string delimiter) {
	std::vector<std::string> query_atoms;
	std::set<std::string> atoms_set;
	boost::algorithm::split(query_atoms, allquery, boost::is_any_of(delimiter), boost::token_compress_on);

	for (int i = 0; i < query_atoms.size(); ++i) {
		boost::trim(query_atoms[i]);
		atoms_set.insert(query_atoms[i]);
	}

	return atoms_set;
}

std::vector<std::set<std::string>> splitAtomQuery_Set(std::string allquery) {
	std::vector<std::set<std::string>> left_right_hands;
	std::vector<std::string> queries;

	boost::algorithm::split(queries, allquery, boost::is_any_of(","), boost::token_compress_on);

	std::set<std::string> left_hand, right_hand;
	for (std::string& q : queries) {
		std::vector<std::string> query_atoms;
		boost::trim(q);
		boost::algorithm::split(query_atoms, q, boost::is_any_of("-"), boost::token_compress_on);

		boost::trim(query_atoms[0]);
		boost::trim(query_atoms[1]);
		left_hand.insert(query_atoms[0]);
		right_hand.insert(query_atoms[1]);
	}

	left_right_hands.push_back(left_hand);
	left_right_hands.push_back(right_hand);
	return left_right_hands;
}

std::set<std::string> splitAtoms(std::string allquery) {
	std::set<std::string> atom_map;
	std::vector<std::string> queries;

	boost::algorithm::split(queries, allquery, boost::is_any_of(","), boost::token_compress_on);

	for (std::string& atom : queries) {
		boost::trim(atom);
		atom_map.insert(atom);
	}

	return atom_map;
}
