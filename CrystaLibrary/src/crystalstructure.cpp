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
#include "crystalstructure.h"

/**************************\
|*******  UNITCELL  *******|
\**************************/

double myUnitCell::getVolume() {
	return abs((Eigen::Vector3d(this->v_a).cross(Eigen::Vector3d(this->v_b))).dot(this->v_c));
}


/**************************\
|*******  LATTICE  ********|
\**************************/

std::vector<double> Lattice::getCellParameters() {
	return this->cell_parameters;
}

void Lattice::setCellParameters(std::vector<double> cell_params) {
	this->cell_parameters.assign(cell_params.begin(), cell_params.end());
}

void Lattice::updateCellParameters() {
	Eigen::Matrix3d axis;

	assert(this->cell_parameters.size() == 6);
	axis = getTransformationMatrixFromFractionalToCartesian(this->cell_parameters);
	this->v_a = axis.col(0);
	this->v_b = axis.col(1);
	this->v_c = axis.col(2);
}

bool Lattice::updateReducedCellParameters() {
	Eigen::Matrix3d fixedX_axis, reduced_transform, reduced_axis;
	bool total_reduced = false;
	
	assert(this->cell_parameters.size() == 6);
#ifdef DEBUG
	std::cout << "\tLattice:" << std::endl;
#endif
	reduced_axis = reduceUnitCell(this->cell_parameters, reduced_transform, total_reduced);
	// The following conversion is needed to fix A as X axis
	fixedX_axis = getTransformationMatrixFromFractionalToCartesian(this->cell_parameters);
	this->v_a = fixedX_axis.col(0);
	this->v_b = fixedX_axis.col(1);
	this->v_c = fixedX_axis.col(2);

	return total_reduced;
}

bool Lattice::updateParametersToReducedMode() {
	this->clearCombinations();
	return this->updateReducedCellParameters();

}

void Lattice::updateParametersToOrginalMode() {
	this->clearCombinations();
	this->updateCellParameters();
}

void Lattice::spanTheLattice(int n, bool positive) {
	//assert(n > 1);
	this->clearCombinations();
	this->addGenerator(this->v_a);
	this->addGenerator(this->v_b);
	this->addGenerator(this->v_c);
	if (positive) {
		this->updateCoefficientsForPositiveDirection(n);
	}
	else {
		// n = i * 2 + 1
		// i = (n-1)/2 = floor(n/2)
		// For n = 2, I want coefficients [-1,0,1]
		// n = 3, I want coefficients [-1,0,1,2]
		// n = 4, I want coefficients [-2,-1,0,1,2]
		// n = 5, [-2,-1,0,1,2,3]
		int k = n / 2;
		if (n % 2 == 0) {
			this->updateCoefficientsForEveryDirection(k);
		}
		else {
			//std::cout << -k << " " << k+1 << std::endl;
			this->updateCoefficientsFromRange(-k,k+1);
		}
	}
	this->updateCoefficientCombinations();
	this->updateGeneratorCombinations();
}

void Lattice::clearTheLattice() {
	this->clearCombinations();
	this->cell_parameters.clear();
	this->v_a.setZero();
	this->v_b.setZero();
	this->v_c.setZero();
}

void Lattice::print_lattice_info() {
	this->print_info();
	std::cout << "Cell parameters details:" << std::endl;
	for (auto &p : this->cell_parameters) {
		std::cout << to_string(p) << std::endl;
	}
	std::cout << "Cartesian vectors details (by row):" << std::endl;

	std::cout << this->v_a.transpose() << std::endl;
	std::cout << this->v_b.transpose() << std::endl;
	std::cout << this->v_c.transpose() << std::endl;
}


/**************************\
|********  MOTIF  *********|
\**************************/

void Motif::setCellMotif(std::vector<std::vector<double>> cell_particles) {
	assert(cell_particles.size() > 0);
	std::vector<double> v = cell_particles[0];
	int size = v.size();
	for (int p_i = 0; p_i < cell_particles.size(); p_i++) {
		Eigen::VectorXd c(size);
		for (int i = 0; i < size; ++i) {
			c[i] = cell_particles[p_i][i];
		}
		this->cell_motif.push_back(c);
	}

	this->removeatomtype_mask = std::vector<bool>(this->cell_motif.size(), false);
	this->cell_motif_cart = std::vector<Eigen::VectorXd>(this->cell_motif.size());
}

void Motif::setTypedUnitcellMotif(std::vector<gemmi::SmallStructure::Site> atom_list) {
	assert(atom_list.size() > 0);
	this->typed_cell_motif.assign(atom_list.begin(), atom_list.end());
}

void Motif::setAtomIndexesInMolecules(std::vector<std::vector<int>> atom_indexes) {
	assert(atom_indexes.size() > 0);
	this->atomindexes_molecules.assign(atom_indexes.begin(), atom_indexes.end());
}

void Motif::generateMolecularCentres(MolecularCentreMode mode) {
	if (!(this->atomindexes_molecules).empty()) {
		for (std::vector<int> mol : this->atomindexes_molecules) {
			switch (mode) {
				case (MolecularCentreMode::Geometric) : {
					molecular_centres.push_back(getGeometricCentre(this->cell_motif, mol));
					break;
				}
				case (MolecularCentreMode::Mass) : {
					molecular_centres.push_back(getCentreOfMass(this->cell_motif, this->typed_cell_motif, mol));
					break;
				}
			}			
		}
	}
	else {
		std::cerr << "No Molecule detected!" << std::endl;
	}
}

void Motif::addMolecularCentresToMotif() {
	// Add molecular centres to the motif
	Eigen::MatrixXd fracToCart = getTransformationMatrixFromFractionalToCartesian(this->cell_parameters);
	for (int i = 0; i < molecular_centres.size(); ++i) {
		this->cell_motif.push_back(molecular_centres[i]);
		this->cell_motif_cart.push_back(fracToCart * molecular_centres[i]);
		gemmi::SmallStructure::Site molcentre_site;
		// occ, u_iso, charge and element are already set
		molcentre_site.fract = gemmi::Fractional(molecular_centres[i][0], molecular_centres[i][1], molecular_centres[i][2]);
		molcentre_site.label = std::string("CM" + to_string(i+1));
		molcentre_site.type_symbol = "CM";
		this->typed_cell_motif.push_back(molcentre_site);
		this->removeatomtype_mask.push_back(false);
	}
}

void Motif::unSelectAtomTypes(std::set<std::string> atomtypes_toremove) {
	// if no specific atom to select return and all atoms will be considered (all true) 
	if (atomtypes_toremove.empty()) {
		// Initialize the mask if not initialized before
		this->removeatomtype_mask = std::vector<bool>(this->cell_motif.size(), false);
	}
	// else mark as false the elements to unselect
	else {
		// Initialize the mask if not initialized before
		this->removeatomtype_mask = std::vector<bool>(this->cell_motif.size(), false);
		for (int i = 0; i < this->cell_motif.size(); ++i) {
			if (atomtypes_toremove.find(this->typed_cell_motif[i].type_symbol) != atomtypes_toremove.end()) {
				this->removeatomtype_mask[i] = true; 
			}
		}
	}
	/*
	std::cout << std::endl;
	std::cout << "Motif number: " << this->cell_motif.size() << std::endl;
	for (int i = 0; i < this->cell_motif.size(); ++i) {
		std::cout << this->atomtype_mask[i] << "-" << this->typed_cell_motif[i].type_symbol << " ";
	}
	std::cout << std::endl;
	*/
}

void Motif::loadSymmetryOperations() {

}

void Motif::updateToFullCrystalMotif() {

}

std::vector<double> Motif::getCellParameters() {
	return this->cell_parameters;
}

void Motif::setCellParameters(std::vector<double> cell_params) {
	this->cell_parameters.assign(cell_params.begin(), cell_params.end());
}

void Motif::updateCellMotif() {
	Eigen::Matrix3d axis;

	assert(this->cell_parameters.size() == 6);
	axis = getTransformationMatrixFromFractionalToCartesian(this->cell_parameters);
	this->basis_matrix = axis;
	this->v_a = axis.col(0);
	this->v_b = axis.col(1);
	this->v_c = axis.col(2);
	// Compute cartesian values for motif points (only original cell) from fractional coordinates
	if (!(this->cell_motif.empty())) {
		for (int i = 0; i < this->cell_motif.size(); ++i) {
			// Convert the new fractional coordinates to the Cartesian ones
			this->cell_motif_cart[i] = axis * this->cell_motif[i];
		}
	}
}

bool Motif::updateReducedCellMotif() {
	std::vector<std::string> v;
	Eigen::Matrix3d old_axis, fixedX_axis, reduced_transform, reduced_axis;
	bool total_reduced = false;

	assert(this->cell_parameters.size() == 6);
	reduced_axis = reduceUnitCell(this->cell_parameters, reduced_transform, total_reduced);

	if (total_reduced) {
		// The following conversion is needed to fix A as X axis
		fixedX_axis = getTransformationMatrixFromFractionalToCartesian(this->cell_parameters);
		this->basis_matrix = fixedX_axis;
		this->v_a = fixedX_axis.col(0);
		this->v_b = fixedX_axis.col(1);
		this->v_c = fixedX_axis.col(2);

		// Apply reduction to every molecule
		// Transform atoms of the molecule and process its centre
		if (!(this->molecular_centres.empty())) {
			//std::cout << "Number of molecule centres: " << this->molecular_centres.size() << "\n";
			//std::cout << "Number of atoms: " << this->cell_motif.size() << "\n";
			for (int i = 0; i < this->molecular_centres.size(); ++i) {
				Eigen::VectorXd mol_centre = this->molecular_centres[i];
				Eigen::VectorXd mol_centre_transformed(3);
				// Transform the centre by reduced axis
				mol_centre_transformed = (mol_centre.transpose() * reduced_transform).transpose();
				// Move the transformed vector inside the unit cell
				Eigen::VectorXd mol_centre_final(mol_centre_transformed);
				processNewFractionalCoordinates(mol_centre, mol_centre_final);
				//std::cout << "\n Orig Mol centre: " << mol_centre.transpose() << " -- Transformed centre: " << mol_centre_transformed.transpose() << " -- Fixed Centre: " << mol_centre_final.transpose() << std::endl;
				//std::cout << "Original Centre of Mass: " << centre_of_mass.transpose() << " -- Transformed: " << centre_of_mass_transformed.transpose() << " -- Final: " << centre_of_mass_final.transpose() << std::endl;
				// Get the translation vector to move the entire molecule
				Eigen::VectorXd translation_vector(mol_centre_final - mol_centre_transformed);
				// transform mol centre
				this->molecular_centres[i] = mol_centre_final;
				for (int atom_ind : atomindexes_molecules[i]) {
					Eigen::VectorXd new_c, translated_c;
					new_c = (this->cell_motif[atom_ind].transpose() * reduced_transform).transpose();
					translated_c = new_c + translation_vector;
					this->cell_motif[atom_ind] = translated_c;
					this->cell_motif_cart[atom_ind] = fixedX_axis * translated_c;
				}
			}
			//std::cout << "Number of atoms (Cart): " << this->cell_motif_cart.size() << "\n";

		}
		else {
			//std::cerr << "No molecular centre detected!" << std::endl;
			if (!(this->cell_motif.empty())) {
				for (int i = 0; i < this->cell_motif.size(); ++i) {
					Eigen::VectorXd new_c(this->cell_motif[i].size());
					// fractional -> new fractional -> cartesian
					new_c = (this->cell_motif[i].transpose() * reduced_transform).transpose();
					//std::cout << "New fractional: " << c.transpose() << " -> " << new_c.transpose() << " -> ";
					// Fix fractional coordinate that go outside the unit cell			
					processNewFractionalCoordinates(this->cell_motif[i], new_c);
					//std::cout << new_c.transpose() << std::endl;
					this->cell_motif[i] = new_c;
					// Convert the new fractional coordinates to the Cartesian ones
					this->cell_motif_cart[i] = fixedX_axis * new_c;
				}
			}
		}
	}
	else {
		this->updateCellMotif();
	}
	return total_reduced;
}

bool Motif::updateMotifToReducedMode() {
	this->clearCombinations();
	return this->updateReducedCellMotif();
}

void Motif::updateMotifToOrginalMode() {
	this->clearCombinations();
	this->updateCellMotif();
}

void Motif::repeatTheMotif(int n, bool positive) {
	//assert(n > 1);
	this->clearCombinations();
	this->addGenerator(this->v_a);
	this->addGenerator(this->v_b);
	this->addGenerator(this->v_c);
	if (positive) {
		this->updateCoefficientsForPositiveDirection(n - 1);
	}
	else {
		// To generate motif point consider lattice point - 1 layer
		// otherwise motif goes outside
		// n = i * 2 + 1 or n = i * 2
		// For n = 2, I want coefficients [-1,0]
		// n = 3, I want coefficients [-1,0,1]
		// n = 4, I want coefficients [-2,-1,0,1]
		// n = 5, [-2,-1,0,1,2]
		int k = n / 2;
		if (n % 2 == 0) {
			this->updateCoefficientsFromRange(-k, k - 1);
		}
		else {
			this->updateCoefficientsForEveryDirection(k);
		}
	}
	this->updateCoefficientCombinations();
	this->updateGeneratorCombinationsForPoints(this->cell_motif_cart);
}

void Motif::repeatTheMotif_Fract(int n, bool positive) {
	//assert(n > 1);
	this->clearCombinations();
	this->addGenerator(Eigen::Vector3d(1, 0, 0));
	this->addGenerator(Eigen::Vector3d(0, 1, 0));
	this->addGenerator(Eigen::Vector3d(0, 0, 1));
	if (positive) {
		this->updateCoefficientsForPositiveDirection(n - 1);
	}
	else {
		// To generate motif point consider lattice point - 1 layer
		// otherwise motif goes outside
		// n = i * 2 + 1 or n = i * 2
		// For n = 2, I want coefficients [-1,0]
		// n = 3, I want coefficients [-1,0,1]
		// n = 4, I want coefficients [-2,-1,0,1]
		// n = 5, [-2,-1,0,1,2]
		int k = n / 2;
		if (n % 2 == 0) {
			this->updateCoefficientsFromRange(-k, k - 1);
		}
		else {
			this->updateCoefficientsForEveryDirection(k);
		}
	}
	this->updateCoefficientCombinations();
	this->updateGeneratorCombinationsForPoints(this->cell_motif);
	for (int i = 0; i < this->getPointCloud().size(); ++i) {
		Eigen::Vector3d cart = this->basis_matrix * this->getPointCloud()[i];
		this->getPointCloud()[i][0] = cart[0];
		this->getPointCloud()[i][1] = cart[1];
		this->getPointCloud()[i][2] = cart[2];
	}
}

void Motif::print_motif_info() {
	this->print_info();
	std::cout << "Cell parameters details:" << std::endl;
	for (auto& p : this->cell_parameters) {
		std::cout << to_string(p) << std::endl;
	}
	std::cout << "Cartesian vectors details (by row):" << std::endl;

	std::cout << this->v_a.transpose() << std::endl;
	std::cout << this->v_b.transpose() << std::endl;
	std::cout << this->v_c.transpose() << std::endl;
}

int Motif::getCellIndexByCellShift(Eigen::VectorXd shift) {
	int cell_index = -1;
	// Retrieve the shift index (cell index)
	for (int i = 0; i < this->c_combination_count; ++i) {
		if ((this->C_combinations[i][0] == shift[0]) && (this->C_combinations[i][1] == shift[1]) && (this->C_combinations[i][2] == shift[2])) {
			cell_index = i;
			break;
		}
	}
	return cell_index;
}

int Motif::getStartingExtensionFactor(int n_mindists, int selected_atoms) {
	// Start with an extension factor that includes enough distances as requested (upper bound of distances)
	int start_extension_factor = MAX(std::ceil(std::cbrt(n_mindists / (double)(selected_atoms - 1))), 3);
	// If the starting extension is not enough increment by 1
	if ((int)std::pow(start_extension_factor, 3) * selected_atoms - 1 < n_mindists) {
		++start_extension_factor;
	}
	// Assure it is odd to have a central unit cell
	if (start_extension_factor % 2 == 0)
		++start_extension_factor;

	return start_extension_factor;
}

std::vector<Eigen::VectorXd> Motif::getMotifByCellShift(Eigen::VectorXd shift) {
	int cell_index = this->getCellIndexByCellShift(shift), start, size;

	size = this->cell_motif.size();
	start = cell_index * size;

	return std::vector<Eigen::VectorXd>(this->cell_motif_cart.begin() + start, this->cell_motif_cart.begin() + start + size);
}

void Motif::getNeighboursDataFromPoint(int point_index, std::multimap<double, int>& neighbour_distances, int n_mindists, int starting_neighbour) {
	// Iterate all the repeated motif 
	double curr_distance;
	int relative_index, motif_number = this->cell_motif.size(), k_distances = neighbour_distances.size();
	for (int ng = starting_neighbour; ng < this->getPointCloud().size(); ++ng) {
		relative_index = ng % motif_number;
		if (this->removeatomtype_mask[relative_index]) continue;
		if (point_index == ng) continue;
		curr_distance = getEuclideanDistance(this->getPointCloud()[point_index], this->getPointCloud()[ng]);
		neighbour_distances.insert(std::pair<double, int>(curr_distance, ng));
		++k_distances;
		if (k_distances > n_mindists) {
			neighbour_distances.erase(--neighbour_distances.end());
		}
	}
}

template<class T>
bool areVectorMatricesEqual(std::vector<std::vector<T>> &m1, std::vector<std::vector<T>> &m2) {
	for (int i = 0; i < m1.size(); ++i) {
		for (int j = 0; j < m1[i].size(); ++j) {
			if (m1[i][j] != m2[i][j]) return false;
		}
	}
	return true;
}

int Motif::autoRepeatTheMotifByType(int n_mindists, std::vector<std::string> atom_types) {
	int motif_number = this->cell_motif.size(), type1_index, selected_motif_number = 0, cell_index = -1, start_point_index, neighbours_count = n_mindists + 1, n_extension = 0, start_extension_factor;
	Eigen::VectorXd shift(3);
	// Check the neighbours matrix and extend the point cloud
	std::vector<std::vector<size_t>> *previous_neighboursIds_matrix = nullptr;
	bool sameNeighbours = false, found_type1 = false;
	// Masks for original and modified selected elements
	std::vector<Eigen::VectorXd> motif_point_cloud_eigen;
	std::map<std::string, std::pair<std::string,int>> map_neighbourstype_count;

	// Get all atom queries
	for (std::string &a_type : atom_types) {
		std::vector<std::string> query_atoms = splitSingleAtomQuery(a_type, "-");
		map_neighbourstype_count.insert(std::make_pair(query_atoms[1], std::make_pair(query_atoms[0],0)));
	}

	// Set the number of neighbours of every type (among user-defined types) to find a good number of distances
	for (int i = 0; i < motif_number; ++i) {
		auto it = map_neighbourstype_count.find(this->typed_cell_motif[i].type_symbol);
		if (it != map_neighbourstype_count.end()) {
			(it->second).second++;
		}
	}
	// Get the minimum number of neighbours
	std::pair<std::string, std::pair<std::string, int>> min_count_pair("X",std::make_pair("X", -1));
	for (auto &e : map_neighbourstype_count) {
		if (min_count_pair.second.second == -1)
			min_count_pair = e;
		else {
			if (e.second < min_count_pair.second) {
				min_count_pair = e;
			}
		}
	}
	// Get the number of selected atoms to compute the starting extension factor
	selected_motif_number = min_count_pair.second.second;
	// Choose an atom from which we start computing distances inside the unit cell (atom related to the minimum neighbour)
	for (int i = 0; i < motif_number; ++i) {		
		if (this->typed_cell_motif[i].type_symbol == min_count_pair.second.first) {
			type1_index = i;
			break;
		}
	}
	// Start with an extension factor that includes enough distances as requested (upper bound of distances)
	start_extension_factor = MAX(std::ceil(std::cbrt(n_mindists / (double)(selected_motif_number - 1))), 3);
	// If the starting extension is not enough increment by 1
	if ((int)std::pow(start_extension_factor, 3) * selected_motif_number - 1 < n_mindists) {
		++start_extension_factor;
	}
	// Assure it is odd to have a central unit cell
	if (start_extension_factor % 2 == 0)
		++start_extension_factor;

	//std::cout << "Extension: " << start_extension_factor << std::endl;
	//std::cout << "motif number: " << motif_number << std::endl;
	//std::cout << "n_mindists: " << n_mindists << std::endl;

	this->repeatTheMotif(start_extension_factor);

	shift << 0, 0, 0;
	// Look for the shift index (cell index) just once
	for (int i = 0; i < this->c_combination_count; ++i) {
		if ((this->C_combinations[i][0] == shift[0]) && (this->C_combinations[i][1] == shift[1]) && (this->C_combinations[i][2] == shift[2])) {
			cell_index = i;
			break;
		}
	}
	start_point_index = cell_index * motif_number;
	int selected_atoms = 0, relative_site_index;

	std::vector<double> query(3);
	query[0] = this->G_combinations[start_point_index + type1_index][0];
	query[1] = this->G_combinations[start_point_index + type1_index][1];
	query[2] = this->G_combinations[start_point_index + type1_index][2];
	
	while (!sameNeighbours) {
		motif_point_cloud_eigen.clear();
		// Add type1 atom to check distances
		motif_point_cloud_eigen.push_back(this->G_combinations[start_point_index + type1_index]);
		// Add neighbours to point cloud
		for (int j = 0; j < this->g_combination_count; ++j) {
			// Point added before, skip it
			if (j == start_point_index + type1_index) continue;
			relative_site_index = j % motif_number;
			if (this->typed_cell_motif[relative_site_index].type_symbol == min_count_pair.second.first) {
				motif_point_cloud_eigen.push_back(this->G_combinations[j]);
			}
		}
		// populate the KD-tree
		KDTreeVectorOfVectorsAdaptor< std::vector<Eigen::VectorXd>, double > mat_index(3, motif_point_cloud_eigen, 10);
		std::vector<std::vector<size_t>> neighboursIds_matrix;
		mat_index.index->buildIndex();

		std::vector<size_t> ret_indexes(neighbours_count);
		std::vector<double> out_dists_sqr(neighbours_count);
		nanoflann::KNNResultSet<double> resultSet(neighbours_count);
		resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
		mat_index.index->findNeighbors(resultSet, &query[0], nanoflann::SearchParams(10));
		// Skip the first neighbour as it is the query point itself
		neighboursIds_matrix.push_back(std::vector<size_t>(ret_indexes.begin() + 1, ret_indexes.end()));
		//std::cout << "Search: " << std::to_string(i) << " completed" << std::endl;

		// If exist a previous matrix check if they are equal
		if (previous_neighboursIds_matrix != nullptr) {
			sameNeighbours = areVectorMatricesEqual(*previous_neighboursIds_matrix, neighboursIds_matrix);
			if (!sameNeighbours) {
				// Extend by 1 layer (for every direction) at each iteration: Index of the central unit cell does NOT change
				this->extendBasisVectorsCoefficientsAndLinearCombinations(1);
				n_extension += 2;
			}
			else {
				n_extension -= 2;
			}
		}
		else {
			this->extendBasisVectorsCoefficientsAndLinearCombinations(1);
			n_extension += 2;
		}
		previous_neighboursIds_matrix = &neighboursIds_matrix;
	}
	//std::cout << "Final Extension: " << start_extension_factor+n_extension << std::endl;
	return start_extension_factor + n_extension;
}

Eigen::VectorXd Motif::getAveragedMinimumDistancesByType(Eigen::VectorXd shift, int n_mindists, std::string type1, std::string type2, Eigen::MatrixXd* neighbours_distances) {

	int motif_number = this->cell_motif.size();
	// Masks for original and modified selected elements
	std::vector<Eigen::VectorXd> motif_point_cloud_eigen;

	// Allocate neighbour distances matrix
	if (neighbours_distances != nullptr) {
		neighbours_distances->resize(motif_number, n_mindists);
	}

	int cell_index = -1;
	// Retrieve the shift index (cell index)
	for (int i = 0; i < this->c_combination_count; ++i) {
		if ((this->C_combinations[i][0] == shift[0]) && (this->C_combinations[i][1] == shift[1]) && (this->C_combinations[i][2] == shift[2])) {
			cell_index = i;
			break;
		}
	}
	/*
	std::cout << std::endl;
	std::cout << "Motif number: " << motif_number << std::endl;
	std::cout << "Linear combinations count: " << this->g_combination_count << std::endl;
	std::cout << "Selected Atoms: " << selected_atoms << std::endl;
	std::cout << "New motif size: " << motif_point_cloud_eigen.size() << std::endl;
	for (int i = 0; i < motif_number; ++i) {
		std::cout << element_mask[i] << "-" << this->typed_cell_motif[i].type_symbol << " ";
	}
	std::cout << std::endl;
	*/

	Eigen::VectorXd amds_vector(n_mindists);
	//std::cout << "amds_vector vector lenght: " << amds_vector.size() << std::endl;
	amds_vector.setZero();
	//std::cout << "amds_vector zeroed" << std::endl;

	int start = cell_index * motif_number, selected_atoms = 0, relative_site_index;
	double curr_dist;
	// Iterate all the repeated motif and search for the k neighbours
	for (int i = 0; i < motif_number; ++i) {
		// Clear the point cloud
		motif_point_cloud_eigen.clear();

		//Add the point in the unit cell to the point cloud
		relative_site_index = (start + i) % motif_number;
		if (this->typed_cell_motif[relative_site_index].type_symbol == type1) {
			motif_point_cloud_eigen.push_back(this->G_combinations[start + i]);
		}
		else {
			// not selected type
			continue;
		}

		for (int j = 0; j < this->g_combination_count; ++j) {
			// Point added before, skip it
			if (j == start + i) continue;
			relative_site_index = j % motif_number;
			if (this->typed_cell_motif[relative_site_index].type_symbol == type2) {
				motif_point_cloud_eigen.push_back(this->G_combinations[j]);
			}
		}
		//std::cout << "Search: " << std::to_string(i) << " started" << std::endl;
		++selected_atoms;

		KDTreeVectorOfVectorsAdaptor< std::vector<Eigen::VectorXd>, double > mat_index(3, motif_point_cloud_eigen, 10);
		mat_index.index->buildIndex();

		std::vector<size_t> ret_indexes(n_mindists + 1);
		std::vector<double> out_dists_sqr(n_mindists + 1);
		nanoflann::KNNResultSet<double> resultSet(n_mindists + 1);

		std::vector<double> query(3);
		query[0] = this->G_combinations[start + i][0];
		query[1] = this->G_combinations[start + i][1];
		query[2] = this->G_combinations[start + i][2];
		resultSet.init(&ret_indexes[0], &out_dists_sqr[0]);
		mat_index.index->findNeighbors(resultSet, &query[0], nanoflann::SearchParams(10));

		// n_mindist+1 because the search tree include the point itself (dist 0)

		//std::cout << "Search class initialized" << std::endl;
		for (size_t j = 1; j < n_mindists + 1; j++) {
			curr_dist = std::sqrt(out_dists_sqr[j]);
			//std::cout << curr_dist << "\t";
			if (neighbours_distances != nullptr) {
				(*neighbours_distances)(i, j - 1) = curr_dist;
			}
			amds_vector[j - 1] += curr_dist;
		}
	}
	amds_vector /= selected_atoms;
	return amds_vector;

}

void Motif::writeExtendedMotifPointsToCSVFormatFile(const char* filename) {
	int motif_number = this->cell_motif.size();
	std::ofstream out_points;
	out_points.open(filename);
	out_points << "X,Y,Z" << std::endl;
	int relative_site_index;
	for (int i = 0; i < this->g_combination_count; ++i) {
		relative_site_index = i % motif_number;
		if (!this->removeatomtype_mask[relative_site_index]) {
			for (int j = 0; j < 3; ++j) {
				if (j == 2) {
					out_points << this->G_combinations[i][j];
				}
				else {
					out_points << this->G_combinations[i][j] << ",";
				}
			}
			out_points << std::endl;
		}
	}
	out_points.close();
}

void Motif::writeUnitCellMotifPointsToCSVFormatFile(const char* filename) {
	int motif_number = this->cell_motif_cart.size();
	std::ofstream out_points;
	out_points.open(filename);
	out_points << "X,Y,Z" << std::endl;
	int relative_site_index;
	for (int i = 0; i < motif_number; ++i) {
		if (!this->removeatomtype_mask[i]) {
			for (int j = 0; j < 3; ++j) {
				if (j == 2) {
					out_points << this->cell_motif_cart[i][j];
				}
				else {
					out_points << this->cell_motif_cart[i][j] << ",";
				}
			}
			out_points << std::endl;
		}
	}
	out_points.close();
}

void Motif::writeFractUnitCellMotifPointsToCSVFormatFile(const char* filename) {
	int motif_number = this->cell_motif.size();
	std::ofstream out_points;
	out_points.open(filename);
	out_points << "X,Y,Z" << std::endl;
	int relative_site_index;
	for (int i = 0; i < motif_number; ++i) {
		if (!this->removeatomtype_mask[i]) {
			for (int j = 0; j < 3; ++j) {
				if (j == 2) {
					out_points << this->cell_motif[i][j];
				}
				else {
					out_points << this->cell_motif[i][j] << ",";
				}
			}
			out_points << std::endl;
		}
	}
	out_points.close();
}
