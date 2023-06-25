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
#ifndef _STRUCTURE_H
#define _STRUCTURE_H

//#define DEBUG

#include <stdio.h>
#include <Eigen/Dense>
#include <set>
#include <map>
#include <iterator>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <thread> // Standard library threads

#include <stringhandler.h>
#include <cifio.h>
#include <chem.h>
#include <geom.h>
#include <unitcellreduction.h>
#include <periodicpointcloud.h>

#include <KDTreeVectorOfVectorsAdaptor.h>

class myUnitCell {
public:
	std::vector<double> cell_parameters;
	Eigen::Vector3d v_a, v_b, v_c;
	Eigen::Matrix3d fractToCart;
		
	double getVolume();
};


class Lattice : public PeriodicPointCloud, public myUnitCell {
private:	
	//	Update the cell parameters (lengths, angles and cartesian vectors) from Cif Document this->doc
	void updateCellParameters();
		
	//	Update and reduce the cell parameters (lengths, angles and cartesian vectors) from Cif Document this->doc
	bool updateReducedCellParameters();

public:
	Lattice() : PeriodicPointCloud(3) {};
	std::vector<double> getCellParameters();

	void setCellParameters(std::vector<double> cell_params);

	//	Update the Lattice to Niggli reduced mode
	bool updateParametersToReducedMode();

	/*
	*	Update the Lattice to Original mode, only cell vectors are updated.
	*	Use spanTheLattice to recalculate points
	*/
	void updateParametersToOrginalMode();

	/*
	*	Generate points in the Lattice around the origin by n copies in all three dimensions
	*	param n: number of unit cells for each direction
	*	param positive: If true unit cell will be extended towards positive directions
	*/
	void spanTheLattice(int n, bool positive = false);

	// Clear the instance of the Lattice
	void clearTheLattice();

	//Print info
	void print_lattice_info();
};

class Motif : public PeriodicPointCloud, public myUnitCell {
private:
	//	Update the cell parameters (lengths, angles and cartesian vectors) from Cif Document this->doc
	void updateCellMotif();

	//	Update and reduce the cell parameters (lengths, angles and cartesian vectors) from Cif Document this->doc
	bool updateReducedCellMotif();

public:
	Motif() : PeriodicPointCloud(3) {};
	// Transformation Matrix to Cartesian
	Eigen::MatrixXd basis_matrix;
	// Typed Fractional coordinates of the motif
	std::vector<gemmi::SmallStructure::Site> typed_cell_motif;
	// Fractional coordinates of the motif
	std::vector<Eigen::VectorXd> cell_motif;
	// Cartesian coordinates of the motif 
	std::vector<Eigen::VectorXd> cell_motif_cart;
	// Atoms indeces grouped by molecules
	std::vector<std::vector<int>> atomindexes_molecules;
	// Fractional Molecular Centres coordinates (Mass or Geometric)
	std::vector<Eigen::VectorXd> molecular_centres;
	// Symmetry operations strings
	std::vector<std::string> symmetry_operations;
	// MAsk for atom type selection
	std::vector<bool> removeatomtype_mask;

	std::vector<double> getCellParameters();

	void setCellParameters(std::vector<double> cell_params);

	/* Load motif points in the data structure and update atomtype_mask*/
	void setCellMotif(std::vector<std::vector<double>> cell_particles);

	void setTypedUnitcellMotif(std::vector<gemmi::SmallStructure::Site> atom_list);

	void setAtomIndexesInMolecules(std::vector<std::vector<int>> atom_indexes);

	/*
	*	Molecular Centre mode
	*/
	enum class MolecularCentreMode {Mass, Geometric};
	
	/*
	*	Param 1 : Molecular Centre Mode (Geometric or Mass)
	*	Generate molecular centres and populate 'molecular_centres' vector
	*/
	void generateMolecularCentres(MolecularCentreMode mode);

	void addMolecularCentresToMotif();

	void loadSymmetryOperations();

	void updateToFullCrystalMotif();

	void unSelectAtomTypes(std::set<std::string> atomtypes_toremove);

	//	Update the Motif to Niggli reduced mode
	bool updateMotifToReducedMode();

	/*
	*	Update the Motif to Original mode, only cell vectors are updated.
	*	Use repeatTheMotif to recalculate points
	*/
	void updateMotifToOrginalMode();

	/*
	*	Generate Cartesian points in n Unit Cell in all three dimensions
	*	param n: number of Unit cells with the showed motif
	*	param positive: If true the points components will be only positives
	*/
	void repeatTheMotif(int n, bool positive = false);

	void repeatTheMotif_Fract(int n, bool positive = false);
	
	void print_motif_info();

	int getCellIndexByCellShift(Eigen::VectorXd shift);

	int getStartingExtensionFactor(int n_mindists, int selected_atoms);

	std::vector<Eigen::VectorXd> getMotifByCellShift(Eigen::VectorXd shift);

	void getNeighboursDataFromPoint(int point_index, std::multimap<double, int>& neighbour_distances, int n_mindists, int starting_neighbour);

	/*
	*	Estimate the Extension Factor and expand the unit cell domain
	*
	*	param1 - int n_mindists:	Number of neghbours to consider for auto expansion
	*	param2 - std::vector<std::string> atom_types:	Atom queries vector (e.g "C-O")
	*	return - int:				Estimated Extension Factor
	*/
	int Motif::autoRepeatTheMotifByType(int n_mindists, std::vector<std::string> atom_types);

	Eigen::VectorXd getAveragedMinimumDistancesByType(Eigen::VectorXd shift, int n_mindists, std::string type1, std::string type2, Eigen::MatrixXd* neighbours_distances = nullptr);

	void writeExtendedMotifPointsToCSVFormatFile(const char* filename);
	void writeUnitCellMotifPointsToCSVFormatFile(const char* filename);
	void writeFractUnitCellMotifPointsToCSVFormatFile(const char* filename);

};
#endif // !_STRUCTURE_H