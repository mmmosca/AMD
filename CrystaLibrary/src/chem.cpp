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
#include "chem.h"

/*****************************************\
|**********  CHEMICAL STRUCTURE  *********|
\*****************************************/

std::vector<gemmi::Restraints::Bond> getBondsFromBlock(gemmi::cif::Block& block) {
	std::vector<gemmi::Restraints::Bond> bonds;

	for (auto row : block.find("_geom_bond_",
		{ "atom_site_label_1", "atom_site_label_2",              // 0, 1
		 "distance" })) {  // 2
		double dist = row.has(2) ? gemmi::cif::as_number(row[2]) : NAN;
		gemmi::Restraints::Bond b;
		b.id1 = { 1, row.str(0) };
		b.id2 = { 1, row.str(1) };
		b.value = dist;
		b.type = gemmi::BondType::Unspec;
		bonds.push_back(b);
	}

	return bonds;
}

std::vector<gemmi::Restraints::Bond> getBondsFromAtomCoordinates(std::vector<gemmi::SmallStructure::Site> typed_unitcellmotif, std::vector<Eigen::VectorXd> motif_points) {
	
	assert(typed_unitcellmotif.size() == motif_points.size());
	std::vector<gemmi::Restraints::Bond> bonds;
	OpenBabel::OBMol compound;
	
	int atom_count = 0;
	
	compound.SetDimension(3);
	compound.BeginModify();
	for (int i = 0; i < typed_unitcellmotif.size(); ++i)
	{
		OpenBabel::OBPairData *label = new OpenBabel::OBPairData;
		OpenBabel::OBAtom *ob_atom = compound.NewAtom();
		ob_atom->SetVector(motif_points[i][0], motif_points[i][1], motif_points[i][2]);
		ob_atom->SetType(typed_unitcellmotif[i].type_symbol);
		ob_atom->SetAtomicNum(typed_unitcellmotif[i].element.atomic_number());
		label->SetAttribute("AtomLabel");
		label->SetValue(typed_unitcellmotif[i].label);
		ob_atom->SetData(label);
		++atom_count;
	}
	
	//std::cout << "Atom count: " << atom_count << std::endl;
	compound.ConnectTheDots();
	//compound.PerceiveBondOrders();
	
	/*
	OpenBabel::OBForceField *ff = dynamic_cast<OpenBabel::OBForceField*>(OpenBabel::OBPlugin::GetPlugin("forcefields", "MMFF94s"));
	if (!ff) {
		std::cout << "Could not find forcefield." << std::endl;
	}
	else {
		std::cout << "Found forcefield: " << ff->Description() << std::endl;
	}
	ff->SetLogFile(&std::cerr);
	ff->SetLogLevel(OBFF_LOGLVL_LOW);
	OpenBabel::OBChargeModel* cm = dynamic_cast<OpenBabel::OBChargeModel*>(OpenBabel::OBPlugin::GetPlugin("charges", "gasteiger"));
	if (!cm) {
		std::cout << "ERROR : Charge model error" << std::endl;
	}

	bool report = cm->ComputeCharges(compound);
	compound.EndModify();

	if (!report) {
		std::cout << "Warning: Charge Model: missing parameter for some atoms" << std::endl;
	}
	if (!(ff->Setup(compound))) {
		std::cout << "Setup molecule with forcefield Failed!" << std::endl;
	}
	ff->EnableCutOff(true);
	ff->SetVDWCutOff(10.0);
	ff->SetElectrostaticCutOff(20.0);
	ff->SetUpdateFrequency(10);
	ff->SteepestDescent(100);
	*/

	int bond_count = 0;
	for (OpenBabel::OBMolBondIter b(compound) ; b ; ++b) {
		OpenBabel::OBAtom* ob_atom1 = b->GetEndAtom();
		OpenBabel::OBAtom* ob_atom2 = b->GetBeginAtom();
		gemmi::Restraints::Bond b_gemmi;
		b_gemmi.id1.atom = ob_atom1->GetData("AtomLabel")->GetValue();
		b_gemmi.id2.atom = ob_atom2->GetData("AtomLabel")->GetValue();
		b_gemmi.value = b->GetLength();
		b_gemmi.type = gemmi::BondType::Unspec;
		bonds.push_back(b_gemmi);
		++bond_count;
	}
	//std::cout << "Bond count: " << bond_count << std::endl;
	return bonds;
}

UndirectedGraph makeUndirectedGraphFromAtomicStructure(gemmi::SmallStructure& crystal_structure, std::vector<gemmi::Restraints::Bond>& bonds, std::map<std::string, int>& mapAtomLabelToIndex) {
	int size = crystal_structure.sites.size();
	UndirectedGraph g(size);
	// Atoms in the graph have the same indexing of the original array where points are (crystal_structure.sites[i])
	for (int i = 0; i < size; ++i) {
		g[i].label = crystal_structure.sites[i].label;
		g[i].element = crystal_structure.sites[i].element;
		g[i].id = i;
		//std::cout << "Inizialize graph vertex: ID: " << g[i].id << " - Label: " << g[i].label << " - Element: " << g[i].element.name() << std::endl;
		mapAtomLabelToIndex.insert(std::make_pair(crystal_structure.sites[i].label, i));
	}
	for (const gemmi::Restraints::Bond& bond : bonds) {
		int n1 = mapAtomLabelToIndex[bond.id1.atom];
		int n2 = mapAtomLabelToIndex[bond.id2.atom];
		boost::add_edge(n1, n2, BondEdge_old{ bond.type }, g);
	}
	return g;
}

UndirectedGraph buildCompoundGraphFromAtomicStructure(std::vector<gemmi::SmallStructure::Site>& atom_sites, std::vector<gemmi::Restraints::Bond>& bonds, std::map<std::string, int>& mapAtomLabelToIndex) {
	int size = atom_sites.size();
	UndirectedGraph g(size);
	// Atoms in the graph have the same indexing of the original array where points are (atom_sites[i])
	for (int i = 0; i < size; ++i) {
		g[i].label = atom_sites[i].label;
		g[i].element = atom_sites[i].element;
		g[i].id = i;
		//std::cout << "Inizialize graph vertex: ID: " << g[i].id << " - Label: " << g[i].label << " - Element: " << g[i].element.name() << std::endl;
		mapAtomLabelToIndex.insert(std::make_pair(atom_sites[i].label, i));
	}
	for (const gemmi::Restraints::Bond& bond : bonds) {
		int n1 = mapAtomLabelToIndex[bond.id1.atom];
		int n2 = mapAtomLabelToIndex[bond.id2.atom];
		boost::add_edge(n1, n2, BondEdge_old{ bond.type }, g);
	}
	return g;
}

std::vector<std::vector<int>> getMoleculesFromUndirectedGraphCompound(UndirectedGraph graph_compound) {
	std::vector<std::vector<int>> molecules_indexes;
	AtomVertexVisitor_old vis;
	
	auto indexmap = boost::get(boost::vertex_index, graph_compound);
	auto colormap = boost::make_vector_property_map<boost::default_color_type>(indexmap);
	
	UndirectedGraph::vertex_iterator it, itEnd;

	for (boost::tie(it, itEnd) = boost::vertices(graph_compound); it != itEnd; ++it) {
		boost::put(colormap, *it, boost::default_color_type::white_color);
		vis.initialize_vertex(*it, graph_compound);
	}

	for (boost::tie(it, itEnd) = boost::vertices(graph_compound); it != itEnd; it++)
	{
		//std::cout << "DFSearch at: " << *it << " - Atom: " << graph_compound[*it].label << " - color: " << colormap[*it] << std::endl;
		if (colormap[*it] == boost::default_color_type::white_color) {
			//std::cout << "Depth first visit starts at: " << *it << " - Atom: " << graph_compound[*it].label << " - color: " << colormap[*it] << std::endl;
			boost::depth_first_visit(graph_compound, *it, vis, colormap);
			std::vector<int> mol_indexes = vis.getMoleculeAtomIndexes();
			molecules_indexes.push_back(mol_indexes);
			vis.clearCurrentMolecule();
		}
	}
	
	return molecules_indexes;
}

Eigen::Vector3d getCentreOfMass(std::vector<Eigen::VectorXd> atoms, std::vector<gemmi::SmallStructure::Site> atom_sites, std::vector<int> atomindexes_molecule) {
	Eigen::Vector3d centre_of_mass;
	double molecularmass = 0;
	centre_of_mass.setZero();
	for (int ind : atomindexes_molecule) {
		centre_of_mass[0] += atom_sites[ind].element.weight() * atoms[ind][0];
		centre_of_mass[1] += atom_sites[ind].element.weight() * atoms[ind][1];
		centre_of_mass[2] += atom_sites[ind].element.weight() * atoms[ind][2];
		molecularmass += atom_sites[ind].element.weight();
	}
	centre_of_mass /= molecularmass;
	return centre_of_mass;
}

Eigen::Vector3d getCentreOfMass(std::vector<Eigen::Vector3d> atoms, std::vector<gemmi::SmallStructure::Site> atom_sites, std::vector<int> atomindexes_molecule) {
	Eigen::Vector3d centre_of_mass;
	double molecularmass = 0;
	centre_of_mass.setZero();
	for (int ind : atomindexes_molecule) {
		centre_of_mass[0] += atom_sites[ind].element.weight() * atoms[ind][0];
		centre_of_mass[1] += atom_sites[ind].element.weight() * atoms[ind][1];
		centre_of_mass[2] += atom_sites[ind].element.weight() * atoms[ind][2];
		molecularmass += atom_sites[ind].element.weight();
	}
	centre_of_mass /= molecularmass;
	return centre_of_mass;
}

Eigen::Vector3d getGeometricCentre(std::vector<Eigen::VectorXd> atoms, std::vector<int> atomindexes_molecule) {
	Eigen::VectorXd geometric_centre(3);
	geometric_centre.setZero();
	for (int ind : atomindexes_molecule) {
		geometric_centre += atoms[ind];
	}
	return geometric_centre / (double)atomindexes_molecule.size();
}

Eigen::Vector3d getGeometricCentre(std::vector<Eigen::Vector3d> atoms, std::vector<int> atomindexes_molecule) {
	Eigen::VectorXd geometric_centre(3);
	geometric_centre.setZero();
	for (int ind : atomindexes_molecule) {
		geometric_centre += atoms[ind];
	}
	return geometric_centre / (double)atomindexes_molecule.size();
}