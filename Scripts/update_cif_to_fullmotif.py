"""
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
"""

import sys
import os

# Check parameters
if (len(sys.argv) != 3):
    print("Must have 2 parameters: input dir and an output dir!")
    sys.exit()
inputdir = sys.argv[1]
outputdir = sys.argv[2]
# Check they are directories
if not(os.path.isdir(inputdir) & os.path.isdir(outputdir)) :
    print("Parameters are not directories!")
    sys.exit()

import ccdc
from ccdc import search, io, molecule
from ccdc.io import MoleculeReader, CrystalReader, EntryReader, CrystalWriter
import numpy as np
import gemmi
from gemmi import cif

# sys.argv[1] Required to have '/' in the path
files = list(gemmi.CifWalk(inputdir))

for file_path in files:
    file_path_split = file_path.split("/")
    file_name = file_path_split[len(file_path_split)-1]
    doc = cif.read(file_path)
    block = doc[0]
    for b in doc:
        #print b.name
        if (b.find_loop('_atom_site_') != None): 
            block = b

    loop_coords = block.init_loop('_atom_site_', ['label', 'type_symbol', 'fract_x', 'fract_y', 'fract_z'])
    #loop_bonds = block.init_loop('_geom_bond_', ['atom_site_label_1', 'atom_site_label_2', 'distance'])

    print("Processed file --> " + file_name)
    crystal_reader = CrystalReader(file_path)
    crystal = crystal_reader[0]
    crystal.assign_bonds()

    print([len(comp.atoms) for comp in crystal.molecule.components])

    fullcrystal_mols = []
    print(list(crystal.symmetry_operators))
    for symmop in list(crystal.symmetry_operators):
        fullcrystal_mols.append(crystal.symmetric_molecule(symmop, force=True))

    for mol in fullcrystal_mols:
        for a in mol.atoms:
            loop_coords.add_row([a.label, a.atomic_symbol, str(a.fractional_coordinates.x), str(a.fractional_coordinates.y), str(a.fractional_coordinates.z) ])

    print(fullcrystal_mols)
    print("Number of Asymmetric Molecules: ", len(fullcrystal_mols))
