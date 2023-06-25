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
#include <cifio.h>

/**************************\
|*******  PUBLIC  *********|
\**************************/
const std::string CIFHandler::cifexample = "cifexample.cif";

const std::vector<std::string> CIFHandler::cell_parameters_tags = { "_cell_length_a", "_cell_length_b", "_cell_length_c",
							"_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma" };

int CIFHandler::findValidBlockIndex(gemmi::cif::Document doc) {
	for (int i = 0; i < doc.blocks.size(); ++i) {
		// Check if it is the right block
		if (doc.blocks[i].find("_atom_site_", { "label", "fract_x", "fract_y" , "fract_z" }).get_loop()) {
			return i;
		}
	}
	return -1;
}

gemmi::cif::Document* CIFHandler::copyCIFrom(gemmi::cif::Document doc) {

	gemmi::cif::Document* new_doc = new gemmi::cif::Document();
	new_doc->source = doc.source;
	std::vector<std::string> curr_row;
	for (gemmi::cif::Block &block : doc.blocks) {
		gemmi::cif::Block* new_block = &(new_doc->add_new_block(block.name));

		for (gemmi::cif::Item& item : block.items) {
			if (item.type == gemmi::cif::ItemType::Pair) {
				new_block->set_pair(item.pair.at(0), item.pair.at(1));
			}
			if (item.type == gemmi::cif::ItemType::Loop) {
				gemmi::cif::Loop* new_loop;
				int tags_number = item.loop.tags.size(), row_number = item.loop.values.size() / item.loop.tags.size();
				new_loop = &(new_block->init_loop(std::string(), item.loop.tags));
				
				// This part copies values from a doc to the new one
				for (int i = 0; i < row_number; ++i) {
					for (int j = 0; j < tags_number; ++j) {
						curr_row.push_back(item.loop.values[i*tags_number+j]);
					}
					new_loop->add_row(curr_row);
					curr_row.clear();
				}				
			}
		}
	}
	return new_doc;
}

void CIFHandler::modifyorAddCIFPairs(gemmi::cif::Document* doc, const std::vector<std::string> parameters_tags, std::vector<std::string> &parameters) {

	if (parameters.size() != 6) {
		std::cerr << "Error in modifyCIF: The vector must have 6 parameters, 3 lengths and 3 angles" << std::endl;
	}
	for (gemmi::cif::Block &block : doc->blocks) {
		if (block.find_pair("_cell_length_a")) {
			for (int i = 0; i < parameters_tags.size(); ++i) {
				block.set_pair(parameters_tags[i], parameters[i]);
			}
		}
	}	
}

void CIFHandler::addCIFPairs(gemmi::cif::Document* doc, const std::vector<std::string> parameters_tags, std::vector<std::string> &parameters) {

	if (parameters.size() != 6) {
		std::cerr << "Error in modifyCIF: The vector must have 6 parameters, 3 lengths and 3 angles" << std::endl;
	}
	for (gemmi::cif::Block &block : doc->blocks) {
		for (int i = 0; i < parameters_tags.size(); ++i) {
			block.set_pair(parameters_tags[i], parameters[i]);
		}
	}
}

bool CIFHandler::modifyOrAddCIFLoop(gemmi::cif::Document* doc, const std::string loop_prefix, const std::vector<std::string> loop_tags, const std::vector<std::vector<std::string>> loop_rows) {
	bool found = false;
	for (gemmi::cif::Block& block : doc->blocks) {
		// Check if it is the right block
		if (block.find_pair("_cell_length_a")) {
			found = true;
			gemmi::cif::Loop& loop = block.init_loop(loop_prefix, loop_tags);
			for (std::vector<std::string> s : loop_rows) {
				loop.add_row(s);
			}			
		}
	}
	if (!found) {
		std::cerr << "No valid crystal block has been found" << std::endl;
	}
	return found;
}

bool CIFHandler::perturbTheLattice(gemmi::cif::Document* doc, std::vector<double> &parameters_offset, std::vector<unsigned short> characteristic_vector) {
	if ( (parameters_offset.size() != 6 ) && ( characteristic_vector.size() != 6 ) ) {
		std::cerr << "Error in perturbTheLattice: The vectors must have 6 parameters" << std::endl;
		return false;
	}
	bool found = false;
	for (gemmi::cif::Block &block : doc->blocks) {
		for (unsigned short i = 0; i < CIFHandler::cell_parameters_tags.size(); ++i) {
			const gemmi::cif::Pair* p = block.find_pair(CIFHandler::cell_parameters_tags[i]);
			if (p != nullptr) {
				if (characteristic_vector[i] == 1) {
					block.set_pair(CIFHandler::cell_parameters_tags[i], std::to_string( stod( p->at(1) ) + parameters_offset[i]));
				}
				found = true;
			}
			else {
				// some pair with tags is missing
				found = false;
				break;
			}
		}
	}
	if (!found) {
		std::cerr << "No pairs have been found" << std::endl;
	}
	return found;
}

gemmi::cif::Document CIFHandler::createEmptyCIF(std::string filename) {

	gemmi::cif::Document new_doc = gemmi::cif::Document();
	
	new_doc.add_new_block(filename);

	return new_doc;

}

std::vector<std::string> CIFHandler::getCIFPairs(gemmi::cif::Document doc, const std::vector<std::string> parameters_tags) {

	std::vector<std::string> values;
	bool found = false;
	for (gemmi::cif::Block &block : doc.blocks) {
		for (int i = 0; i < parameters_tags.size(); ++i) {
			const gemmi::cif::Pair* p = block.find_pair(parameters_tags[i]);
			if (p != nullptr) {
				values.push_back(p->at(1));
				found = true;
			}
			else {
				// some pair with tags is missing
				found = false;
				break;
			}
		}
	}
	if (!found) {
		std::cerr << "Some pairs have not been found. File corrupted!" << std::endl;
		values.clear(); 
	}
	return values;
}

std::vector<double> CIFHandler::getCIFPairsValues(gemmi::cif::Document doc, const std::vector<std::string> parameters_tags) {

	std::vector<double> values;
	bool found = false;
	for (gemmi::cif::Block &block : doc.blocks) {
		for (int i = 0; i < parameters_tags.size(); ++i) {
			const gemmi::cif::Pair* p = block.find_pair(parameters_tags[i]);
			if (p != nullptr) {
				values.push_back(std::stod(p->at(1)));
				found = true;
			}
			else {
				// some pair with tags is missing
				found = false;
				break;
			}
		}
	}
	if (!found) {
		std::cerr << "Some pairs have not been found. File corrupted!" << std::endl;
		values.clear();
	}
	return values;
}

std::vector<std::vector<double>> CIFHandler::getCIFAtomCoordinates(gemmi::cif::Document doc) {
	std::vector<std::vector<double>> coordinates;
	bool found = false;
	std::string loop_tag("_atom_site_label"), strings;
	int x, y, z;
	for (gemmi::cif::Block &block : doc.blocks) {
		gemmi::cif::Column col = block.find_loop(loop_tag);
		gemmi::cif::Loop* loop = col.get_loop();
		if (loop == nullptr) {
			found = false;
			continue;
		}
		found = true;
		x = loop->find_tag(std::string("_atom_site_fract_x"));
		y = loop->find_tag(std::string("_atom_site_fract_y"));
		z = loop->find_tag(std::string("_atom_site_fract_z"));
		//std::cout << "INDEXES: " << x << " " << y << " " << z << " " << std::endl;
		// number of atoms
		int n_atoms = loop->length();
		// number of strings in the atom line
		int n_values = loop->width();
		//std::cout << length << " " << width << std::endl;

		// Iterate the line
		for (int i = 0; i < n_atoms; ++i) {
			std::vector<double> coord(3);
			coord[0] = std::stod(loop->values[x]);
			coord[1] = std::stod(loop->values[y]);
			coord[2] = std::stod(loop->values[z]);
			//for (auto &c : coord) std::cout << std::to_string(c) << " "; std::cout << std::endl;
			coordinates.push_back(coord);
			x += n_values; y += n_values; z += n_values;
		}
	}

	if (!found) {
		std::cerr << "Some coordinates have not been found. File corrupted!" << std::endl;
		coordinates.clear();
	}
	return coordinates;
}

gemmi::cif::Document* CIFHandler::createEmptyCIFDocumentFrom(gemmi::cif::Document doc, std::string filename) {

	gemmi::cif::Document* new_doc = new gemmi::cif::Document();
	new_doc->source = filename;
	std::vector<std::string> curr_row;
	int spec_block_number = 0;
	for (gemmi::cif::Block &block : doc.blocks) {

		gemmi::cif::Block* new_block = (block.name.compare("global") == 0) ? &(new_doc->add_new_block(block.name)) : &(new_doc->add_new_block(std::string("specifics_") + std::to_string(spec_block_number++)));

		for (gemmi::cif::Item& item : block.items) {
			if (item.type == gemmi::cif::ItemType::Pair) {
				new_block->set_pair(item.pair.at(0), std::string("\t?"));
			}
			if (item.type == gemmi::cif::ItemType::Loop) {
				gemmi::cif::Loop* new_loop;
				int tags_number = item.loop.tags.size();
				new_loop = &(new_block->init_loop(std::string(), item.loop.tags));
				// At least one add_row must be used to print the tags
				for (int j = 0; j < tags_number; ++j) {
					curr_row.push_back(std::string("?"));
				}
				new_loop->add_row(curr_row);
				curr_row.clear();
			}
		}
	}
	return new_doc;
}
