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
#ifndef _AMDS_H
#define _AMDS_H

#include <crystalstructure.h>

class AMDs {
private:
	Eigen::VectorXd AMDs_vector;
	Eigen::MatrixXd PDD_matrix;

	int getStartingEstimatedExtensionFactor(int n_mindists, int selected_atoms);
	void updateNeighbourhoodData_MultiMap(Motif& motif, int n_mindists);
	void updateNeighbourhoodData_KDTree(Motif& motif, int n_mindists);
public:
	AMDs(Motif& motif, int n_mindists);
	Eigen::VectorXd getAMDsVector();
	Eigen::MatrixXd getPDDMatrix();
};

#endif // !_AMDS_H