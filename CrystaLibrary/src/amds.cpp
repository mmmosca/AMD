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
#include <amds.h>

int AMDs::getStartingEstimatedExtensionFactor(int n_mindists, int selected_atoms) {
	// Start with an extension factor that includes enough distances as requested (upper bound of distances)
	int start_extension_factor;
	
	if (selected_atoms == 1) {
		// Fix for 1 motif point unit cell
		start_extension_factor = MAX(std::ceil(std::cbrt(n_mindists / static_cast<double>(selected_atoms))), 3);
	}
	else {
		start_extension_factor = MAX(std::ceil(std::cbrt(n_mindists / static_cast<double>(selected_atoms - 1))), 3);
	}
	// If the starting extension is not enough increment by 1
	if ((int)std::pow(start_extension_factor, 3) * selected_atoms - 1 < n_mindists) {
		++start_extension_factor;
	}
	// Assure it is odd to have a central unit cell
	if (start_extension_factor % 2 == 0)
		++start_extension_factor;

	return start_extension_factor;
}

template<typename T>
bool areMatricesEqual(std::vector<std::vector<T>>& m1, std::vector<std::vector<T>>& m2) {
	for (int i = 0; i < m1.size(); ++i) {
		for (int j = 0; j < m1[i].size(); ++j) {
			if (m1[i][j] != m2[i][j]) {
				return false;
			}
		}
	}
	return true;
}

void AMDs::updateNeighbourhoodData_MultiMap(Motif& motif, int n_mindists) {
	Eigen::Vector3d shift(0, 0, 0);
	int motif_number = motif.cell_motif.size(), selected_atoms = 0;

	// Get the number of selected atoms to compute the starting extension factor
	
	for (int i = 0; i < motif_number; ++i) {
		if (!motif.removeatomtype_mask[i]) {
			++selected_atoms;
		}
	}
	// Start with an extension factor that includes enough distances as requested (upper bound of distances)
	int start_extension_factor = this->getStartingEstimatedExtensionFactor(n_mindists, selected_atoms);
	motif.repeatTheMotif(start_extension_factor);

	// Retrieve the shift index (cell index)
	int cell_index = motif.getCellIndexByCellShift(shift), starting_neighbour = 0;
	int start_central = cell_index * motif_number, neighbour_counter;
	bool sameNeighbours = false;
	std::vector<std::vector<int>> previous_neighboursIds_matrix;
	std::vector<std::multimap<double, int>> neighbourhooddata_matrix(motif_number);

	// Extend the motif point cloud if neighbour matrices are different
	while (!sameNeighbours) {

		std::vector<std::vector<int>> neighboursIds_matrix;
		// Iterate all the repeated motif and search for the k neighbours
		for (int i = 0; i < motif_number; ++i) {
			std::vector<int> ith_neighbour_indeces;
			if (motif.removeatomtype_mask[i]) continue;

			neighbour_counter = 0;
			motif.getNeighboursDataFromPoint(start_central + i, neighbourhooddata_matrix[i], n_mindists, starting_neighbour);
			for (auto it = neighbourhooddata_matrix[i].begin(); it != neighbourhooddata_matrix[i].end() ; ++it) {
				ith_neighbour_indeces.push_back(it->second);
				++neighbour_counter;
			}
			neighboursIds_matrix.push_back(ith_neighbour_indeces);
		}
		starting_neighbour = motif.getPointCloud().size();

		// If exist a previous matrix check if they are equal
		if (!previous_neighboursIds_matrix.empty()) {
			sameNeighbours = areMatricesEqual(previous_neighboursIds_matrix, neighboursIds_matrix);
			if (!sameNeighbours) {
				// Extend by 1 layer (for every direction) at each iteration: Index of the central unit cell does NOT change
				motif.extendBasisVectorsCoefficientsAndLinearCombinations(1, motif.cell_motif_cart);
			}
		}
		else {
			// No previous matrix (first loop)
			motif.extendBasisVectorsCoefficientsAndLinearCombinations(1, motif.cell_motif_cart);
		}
		previous_neighboursIds_matrix = neighboursIds_matrix;
	}

	// Update PDD and AMDs with distances
	for (int mp = 0; mp < neighbourhooddata_matrix.size(); ++mp) {
		if (!neighbourhooddata_matrix[mp].empty()) {
			int ng = 0;
			for (auto& m : neighbourhooddata_matrix[mp]) {
				this->PDD_matrix(mp, ng) = m.first;
				this->AMDs_vector[ng] += m.first;
				++ng;
			}
		}
	}
	this->AMDs_vector /= selected_atoms;
}

void AMDs::updateNeighbourhoodData_KDTree(Motif& motif, int n_mindists) {
	Eigen::Vector3d shift(0, 0, 0);
	int motif_number = motif.cell_motif.size(), selected_atoms = 0;
	std::vector<Eigen::VectorXd> motif_point_cloud_eigen;

	for (int i = 0; i < motif_number; ++i) {
		if (!motif.removeatomtype_mask[i]) {
			++selected_atoms;
		}
	}
	// Start with an extension factor that includes enough distances as requested (upper bound of distances)
	int start_extension_factor = this->getStartingEstimatedExtensionFactor(n_mindists, selected_atoms);
	motif.repeatTheMotif(start_extension_factor);
	// motif.repeatTheMotif_Fract(start_extension_factor);
	// Retrieve the shift index (cell index)
	int cell_index = motif.getCellIndexByCellShift(shift), starting_neighbour = 0;
	int start_central = cell_index * motif_number, neighbour_counter;
	bool sameNeighbours = false;
	std::vector<std::vector<size_t>> previous_neighboursIds_matrix;
	
	//std::cout << "Motif number in the unit cell: " << motif.cell_motif_cart.size() << std::endl;

	// Extend the motif point cloud if neighbour matrices are different
	std::vector<std::vector<double>> previous_sqr_neighbourdistances_matrix;
	std::vector<std::vector<double>> sqr_neighbourdistances_matrix(motif_number);
	while (!sameNeighbours) {
		std::vector<std::vector<size_t>> neighboursIds_matrix;

		// Iterate all the repeated motif and search for the k neighbours
		int relative_site_index, total_points = motif.getPointCloud().size();
		//std::cout << "Starting neighbour: " << starting_neighbour << std::endl;
		//std::cout << "Number of points: " << total_points << std::endl;
		for (int i = starting_neighbour; i < total_points; ++i) {
			relative_site_index = i % motif_number;
			if (!motif.removeatomtype_mask[relative_site_index]) {
				motif_point_cloud_eigen.push_back(motif.getPointCloud()[i]);
			}
		}
		starting_neighbour = total_points;

		// populate the KD-tree
		KDTreeVectorOfVectorsAdaptor< std::vector<Eigen::VectorXd>, double > kdtree(3, motif_point_cloud_eigen, 10);
		kdtree.index->buildIndex();

		for (int i = 0; i < motif_number; ++i) {
			if (motif.removeatomtype_mask[i]) continue;
			std::vector<size_t> ith_neighbour_indeces(n_mindists+1);
			sqr_neighbourdistances_matrix[i] = std::vector<double>(n_mindists + 1);
			nanoflann::KNNResultSet<double> resultSet(n_mindists + 1);

			std::vector<double> query(3);
			query[0] = motif.getPointCloud()[start_central + i][0];
			query[1] = motif.getPointCloud()[start_central + i][1];
			query[2] = motif.getPointCloud()[start_central + i][2];

			resultSet.init(&ith_neighbour_indeces[0], &sqr_neighbourdistances_matrix[i][0]);
			kdtree.index->findNeighbors(resultSet, &query[0], nanoflann::SearchParams(10));
			// Skip the first neighbour as it is the query point itself
			neighboursIds_matrix.push_back(std::vector<size_t>(ith_neighbour_indeces.begin() + 1, ith_neighbour_indeces.end()));
		}
		//std::cout << "PREVIOUS pointer: " << &previous_neighboursIds_matrix << std::endl;
		//std::cout << "NEXT pointer: " << &neighboursIds_matrix << std::endl;
		
		// If exist a previous matrix check if they are equal
		//if (!previous_neighboursIds_matrix.empty()) {
		if (!previous_sqr_neighbourdistances_matrix.empty()) {
			//sameNeighbours = areMatricesEqual(previous_neighboursIds_matrix, neighboursIds_matrix);
			sameNeighbours = areMatricesEqual(previous_sqr_neighbourdistances_matrix, sqr_neighbourdistances_matrix);
			//std::cout << "SAME_NEIGHBOURS: " << sameNeighbours << std::endl;
			if (!sameNeighbours) {
				// Extend by 1 layer (for every direction) at each iteration: Index of the central unit cell does NOT change
				motif.extendBasisVectorsCoefficientsAndLinearCombinations(1, motif.cell_motif_cart);
				/*
				motif.extendBasisVectorsCoefficientsAndLinearCombinations(1, motif.cell_motif);
				for (int i = starting_neighbour; i < total_points; ++i) {
					Eigen::Vector3d cart = motif.basis_matrix * motif.getPointCloud()[i];
					motif.getPointCloud()[i][0] = cart[0];
					motif.getPointCloud()[i][1] = cart[1];
					motif.getPointCloud()[i][2] = cart[2];
				}
				*/
			}
			else {
				//Point Cloud is correctly extended -> Update data structure
				double current_distance;
				for (int i = 0; i < sqr_neighbourdistances_matrix.size(); i++) {
					// Start from 1 j=1 to skip the 0 distances
					for (int j = 1; j < sqr_neighbourdistances_matrix[i].size(); j++) {
						current_distance = std::sqrt(sqr_neighbourdistances_matrix[i][j]);
						this->PDD_matrix(i, j - 1) = current_distance;
						this->AMDs_vector[j - 1] += current_distance;
					}
				}
				this->AMDs_vector /= selected_atoms;
			}
		}
		else {
			// No previous matrix (first loop)
			//std::cout << "FIRST LOOP EXTENSION" << std::endl;
			motif.extendBasisVectorsCoefficientsAndLinearCombinations(1, motif.cell_motif_cart);
			/*
			motif.extendBasisVectorsCoefficientsAndLinearCombinations(1, motif.cell_motif);
			for (int i = starting_neighbour; i < total_points; ++i) {
				Eigen::Vector3d cart = motif.basis_matrix * motif.getPointCloud()[i];
				motif.getPointCloud()[i][0] = cart[0];
				motif.getPointCloud()[i][1] = cart[1];
				motif.getPointCloud()[i][2] = cart[2];
			}
			*/
		}
		previous_neighboursIds_matrix = neighboursIds_matrix;
		previous_sqr_neighbourdistances_matrix = sqr_neighbourdistances_matrix;
	}
}

AMDs::AMDs(Motif &motif, int n_mindists)
{
	this->PDD_matrix = Eigen::MatrixXd(motif.cell_motif.size(), n_mindists);
	this->PDD_matrix.setZero();
	this->AMDs_vector = Eigen::VectorXd(n_mindists);
	this->AMDs_vector.setZero();

	//this->updateNeighbourhoodData_MultiMap(motif, n_mindists);
	this->updateNeighbourhoodData_KDTree(motif, n_mindists);
}

Eigen::VectorXd AMDs::getAMDsVector() {
	return this->AMDs_vector;
}

Eigen::MatrixXd AMDs::getPDDMatrix() {
	return this->PDD_matrix;
}