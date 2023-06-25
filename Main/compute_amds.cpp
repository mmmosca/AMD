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

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <ctime>
#include <crystalstructure.h>
#include <amds.h>
#include <boost/filesystem.hpp>
#include "gemmi/dirwalk.hpp"
#include <cmd.h>
#include <boost/asio/thread_pool.hpp>
#include <boost/asio.hpp>
#include <boost/bind.hpp>



#define PRINTUSAGE() printf("|--- USAGE ---|:\nRequired options:\n"\
							"\t-inputdir [Input Folder with CIF files]\n"\
							"\t-outputdir [Output Folder for Results]\n"\
							"Optional:\n"\
							"\t-mindistances [Number of Averaged Minimum Distances, default: 200]\n"\
							"\t-threads [t]: Run AMDs computation on t threads (by default t=1)\n");

namespace boost_dir = boost::filesystem;

struct AMDThreadCommons {
	std::vector<std::string> atompairs_toselect_map;
	bool INPUTDIR_OPT = false, OUTPUTDIR_OPT = false, POINTS_OPT = false,
		MINDIST_OPT = false, EXTEND_OPT = false;
	int n_mindists, extend = 3, threads = 1;
	std::string output_dir;
};

struct AMDThreadIO {
	std::string inputfile_path, filename, outputfile_string;
	bool reduced;
};

void AMDs_Thread(AMDThreadCommons& thread_commons, AMDThreadIO& thread_io) {
	//std::cout << thread_io.inputfile_path << std::endl;
	gemmi::cif::Document d = gemmi::cif::read_file(thread_io.inputfile_path);
	int block_index;
	if ((block_index = CIFHandler::findValidBlockIndex(d)) == -1) {
		thread_io.outputfile_string = std::string("File corrupted: " + thread_io.inputfile_path);
		return;
		// Terminate thread.
	}
	UndirectedGraph compound_graph;
	// map the Atom label in the atomic structure index
	std::map<std::string, int> mapAtomLabelToIndex;

	gemmi::SmallStructure crystal_structure = gemmi::make_small_structure_from_block(d.blocks[block_index]);
	gemmi::UnitCell unit_cell = crystal_structure.cell;
	std::vector<double> cell_parameters{ unit_cell.a, unit_cell.b ,unit_cell.c, unit_cell.alpha, unit_cell.beta, unit_cell.gamma };

	Motif motif;
	std::vector<std::vector<double>> cell_motif;
	// Filter filename
	thread_io.filename = std::string(thread_io.inputfile_path);
	thread_io.filename.erase(0, thread_io.filename.find_last_of("\\/") + 1);
	thread_io.filename.erase(thread_io.filename.find_last_of('.'), thread_io.filename.length());

	if (thread_commons.threads == 1) {
		std::cout << '\r' << "Processing file: " << thread_io.filename;
	}
	motif.setCellParameters(cell_parameters);

	for (auto& site : crystal_structure.sites) {
		std::vector<double> coordinates(3);
		coordinates[0] = site.fract.x;
		coordinates[1] = site.fract.y;
		coordinates[2] = site.fract.z;
		cell_motif.push_back(coordinates);
	}
	motif.setCellMotif(cell_motif);
	motif.setTypedUnitcellMotif(crystal_structure.sites);
	
	// Update to cartesian coordinates
	for (int i = 0; i < crystal_structure.sites.size(); ++i) {
		gemmi::Position cartesian_position = unit_cell.orthogonalize(crystal_structure.sites[i].fract);
		motif.cell_motif_cart[i] = Eigen::Vector3d(cartesian_position.x, cartesian_position.y, cartesian_position.z);
	}

	std::vector<gemmi::Restraints::Bond> bonds = getBondsFromAtomCoordinates(motif.typed_cell_motif, motif.cell_motif_cart);
	compound_graph = makeUndirectedGraphFromAtomicStructure(crystal_structure, bonds, mapAtomLabelToIndex);
	std::vector<std::vector<int>> molecules_indexes = getMoleculesFromUndirectedGraphCompound(compound_graph);
	motif.setAtomIndexesInMolecules(molecules_indexes);
	motif.generateMolecularCentres(Motif::MolecularCentreMode::Geometric);

	motif.updateMotifToReducedMode();

	AMDs amds(motif, thread_commons.n_mindists);

	std::string path_to_output_file;
	if (thread_commons.POINTS_OPT) {	
		path_to_output_file = thread_commons.output_dir + std::string("/") + thread_io.filename + std::string("_extendedmotif.csv");
		motif.writeExtendedMotifPointsToCSVFormatFile(path_to_output_file.c_str());
		path_to_output_file = thread_commons.output_dir + std::string("/") + thread_io.filename + std::string("_unitcellmotif.csv");
		motif.writeUnitCellMotifPointsToCSVFormatFile(path_to_output_file.c_str());
	}
	
	if (thread_commons.MINDIST_OPT) {
		// Averaged Minimum Distances
		Eigen::VectorXd amds_vector;
		Eigen::MatrixXd PDD;
		Eigen::Vector3d v(0,0,0);
		std::stringstream sout;

		sout << std::setprecision(25);
		amds_vector = amds.getAMDsVector();
		sout << thread_io.filename;
		for (int i = 0; i < thread_commons.n_mindists; ++i) {
			sout << "," << amds_vector[i];
		}
		sout << std::endl;

		thread_io.outputfile_string = sout.str();
	}
}


int main(int argc, char* argv[]) {
	std::string input_dir, output_dir;
	AMDThreadCommons t_c;
	CommandLine cmd;
	bool THREADS_OPT = false;
	t_c.n_mindists = 200;
	t_c.MINDIST_OPT = true;
	t_c.extend = 3;
	t_c.EXTEND_OPT = true;
	char* w;
	while ( ( w = cmd.mygetoptW(argc, argv, "inputdir:|outputdir:|mindistances:|threads:|") ) != NULL) {
		if (strcmp(w, "inputdir") == 0) {
			input_dir.assign(cmd.myoptarg);
			t_c.INPUTDIR_OPT = true;
			continue;
		}
		if (strcmp(w, "outputdir") == 0) {
			output_dir.assign(cmd.myoptarg);
			t_c.OUTPUTDIR_OPT = true;
			continue;
		}
		if (strcmp(w, "threads") == 0) {
			t_c.threads = stoi(cmd.myoptarg);
			THREADS_OPT = true;
			continue;
		}
		if (strcmp(w, "mindistances") == 0) {
			t_c.n_mindists = stoi(cmd.myoptarg);
			t_c.MINDIST_OPT = true;
			continue;
		}
		
	}
	
	if (!(t_c.INPUTDIR_OPT && t_c.OUTPUTDIR_OPT) || !t_c.MINDIST_OPT) {
		PRINTUSAGE();
		exit(EXIT_FAILURE);
	}
	t_c.output_dir = output_dir;

	std::cout << "Please wait..." << std::endl;
	
	auto start = std::chrono::system_clock::now();

	boost::asio::thread_pool t_pool(t_c.threads);
	std::vector<AMDThreadIO> thread_results;

	for (std::string& cif_filepath : gemmi::CifWalk(input_dir)) {
		AMDThreadIO t_io;
		t_io.inputfile_path = std::string(cif_filepath);
		thread_results.push_back(t_io);
		
	}
	for (AMDThreadIO& t_io : thread_results) {
		boost::asio::post(t_pool, boost::bind(&AMDs_Thread, boost::ref(t_c), boost::ref(t_io)));
	}
	t_pool.join();

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	int elapsed_s = elapsed_seconds.count();
	int hours = floor(elapsed_s / 3600), minutes = floor(fmod(elapsed_s,3600) / 60), seconds = fmod(elapsed_s, 60);
	std::cout << "Elapsed time: " << hours << ':' << minutes << ':' << seconds << "s\n";

	/***********************\
	|*** POST-PROCESSING ***|
	\***********************/
	start = std::chrono::system_clock::now();

	std::string	distancetuples_file_name(output_dir + "/AMDs_" + std::to_string(t_c.n_mindists) + ".csv");
	std::ofstream out_disttuples;

	out_disttuples.open(distancetuples_file_name);
	out_disttuples << "ID";
	for (int i = 1; i <= t_c.n_mindists; ++i) {
		out_disttuples << ",AMD_" << i;
	}
	out_disttuples << std::endl;

	std::cout << "Writing all crystals AMDs on file..." << std::endl;
	for (int i = 0; i < thread_results.size(); ++i) {
		out_disttuples << thread_results[i].outputfile_string;
	}
	out_disttuples.close();
	

	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	elapsed_s = elapsed_seconds.count();
	hours = floor(elapsed_s / 3600);
	minutes = floor(fmod(elapsed_s, 3600) / 60);
	seconds = fmod(elapsed_s, 60);
	std::cout << "Elapsed time: " << hours << ':' << minutes << ':' << seconds << "s\n";

	std::cout << "Completed!" << std::endl;
}