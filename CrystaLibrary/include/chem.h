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
#ifndef _CHEM_H
#define _CHEM_H

#include <utility>
#include <iostream>
#include <string>

#include <Eigen/Dense>

#include "gemmi/chemcomp.hpp"
#include "gemmi/smcif.hpp"

#include "boost/graph/graph_traits.hpp"
#include "boost/graph/adjacency_list.hpp"
#include <boost/graph/depth_first_search.hpp>

#include <openbabelsupport.h>

enum class MolecularCentreMode_old { Mass, Geometric };

struct AtomVertex_old {
	int id;
	std::string label;
	gemmi::Element element = gemmi::El::X;
};

struct BondEdge_old {
	gemmi::BondType type;
};

using UndirectedGraph = boost::adjacency_list<boost::vecS, boost::vecS,
	boost::undirectedS,
	AtomVertex_old, BondEdge_old>;

class AtomVertexVisitor_old : public boost::default_dfs_visitor
{
public:
	AtomVertexVisitor_old() : current_molecule(new std::vector<int>()) {}

	void discover_vertex(int v, const UndirectedGraph& g)
	{
		current_molecule->push_back(g[v].id);
		// else keep going with searching
		return;
	}

	void clearCurrentMolecule() {
		current_molecule->clear();
	}

	std::vector<int> getMoleculeAtomIndexes() const {
		return *current_molecule;
	}

private:
	boost::shared_ptr<std::vector<int>> current_molecule;
};


std::vector<gemmi::Restraints::Bond> getBondsFromBlock(gemmi::cif::Block& block);

std::vector<gemmi::Restraints::Bond> getBondsFromAtomCoordinates(std::vector<gemmi::SmallStructure::Site> typed_unitcellmotif, std::vector<Eigen::VectorXd> motif_points);

UndirectedGraph makeUndirectedGraphFromAtomicStructure(gemmi::SmallStructure& crystal_structure, std::vector<gemmi::Restraints::Bond>& bonds, std::map<std::string, int>& mapAtomLabelToIndex);

UndirectedGraph buildCompoundGraphFromAtomicStructure(std::vector<gemmi::SmallStructure::Site>& atom_sites, std::vector<gemmi::Restraints::Bond>& bonds, std::map<std::string, int>& mapAtomLabelToIndex);

std::vector<std::vector<int>> getMoleculesFromUndirectedGraphCompound(UndirectedGraph graph_compound);

Eigen::Vector3d getCentreOfMass(std::vector<Eigen::VectorXd> atoms, std::vector<gemmi::SmallStructure::Site> atom_sites, std::vector<int> atomindexes_molecule);
Eigen::Vector3d getCentreOfMass(std::vector<Eigen::Vector3d> atoms, std::vector<gemmi::SmallStructure::Site> atom_sites, std::vector<int> atomindexes_molecule);

Eigen::Vector3d getGeometricCentre(std::vector<Eigen::VectorXd> atoms, std::vector<int> atomindexes_molecule);
Eigen::Vector3d getGeometricCentre(std::vector<Eigen::Vector3d> atoms, std::vector<int> atomindexes_molecule);

#endif //_CHEM_H