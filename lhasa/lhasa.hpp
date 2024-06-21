/* lhasa/lhasa.hpp
 * 
 * Copyright 2023 by Global Phasing Ltd.
 * Author: Jakub Smulski
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifndef LHASA_HPP
#define LHASA_HPP
#include <rdkit/GraphMol/RWMol.h>
#include <string>
#include <memory>
#include "../layla/ligand_editor_canvas.hpp"
#include <emscripten/val.h>


namespace lhasa {

std::unique_ptr<RDKit::RWMol> rdkit_mol_from_smiles(std::string smiles);
std::unique_ptr<RDKit::RWMol> rdkit_mol_from_pickle(std::string pickle_string);

std::string rdkit_mol_to_smiles(RDKit::ROMol& mol);

unsigned int append_from_smiles(CootLigandEditorCanvas& canvas, std::string smiles);
unsigned int append_from_pickle_base64(CootLigandEditorCanvas& canvas, std::string pickle_string);

std::unique_ptr<coot::ligand_editor_canvas::ActiveTool> make_active_tool(emscripten::val t);
coot::ligand_editor_canvas::ElementInsertion element_insertion_from_symbol(std::string sym);

}


#endif