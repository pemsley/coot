#ifndef LHASA_HPP
#define LHASA_HPP
#include <rdkit/GraphMol/RWMol.h>
#include <string>
#include <memory>
#include "../layla/ligand_editor_canvas.hpp"
#include <emscripten/val.h>

namespace lhasa {

std::unique_ptr<RDKit::RWMol> rdkit_mol_from_smiles(std::string smiles);
void append_from_smiles(CootLigandEditorCanvas& canvas, std::string smiles);
std::string rdkit_mol_to_smiles(RDKit::ROMol& mol);
std::unique_ptr<coot::ligand_editor_canvas::ActiveTool> make_active_tool(emscripten::val t);
coot::ligand_editor_canvas::ElementInsertion element_insertion_from_symbol(std::string sym);

}


#endif