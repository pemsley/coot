/* lhasa/lhasa.cpp
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
#include "lhasa.hpp"
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>
#include <rdkit/GraphMol/MolPickler.h>
#include "glog_replacement.hpp"

#ifdef LHASA_STANDALONE
#define TINYGLTF_IMPLEMENTATION
#include "../coot-utils/tiny_gltf.h"
#else
// This must not be here because it will cause a linker error (in Moorhen)
// #define TINYGLTF_IMPLEMENTATION
// #include "../coot-utils/tiny_gltf.h"
// ..and so we have to do a forward declaration
namespace tinygltf {
    std::string base64_encode(unsigned char const *bytes_to_encode, unsigned int in_len);
    std::string base64_decode(std::string const &encoded_string);
}
#endif

std::unique_ptr<RDKit::RWMol> lhasa::rdkit_mol_from_smiles(std::string smiles) {
    std::unique_ptr<RDKit::RWMol> ret(RDKit::SmilesToMol(smiles));
    return ret;
}

std::unique_ptr<RDKit::RWMol> lhasa::rdkit_mol_from_pickle(std::string pickle_string) {
    std::unique_ptr<RDKit::RWMol> ret = std::make_unique<RDKit::RWMol>();
    RDKit::MolPickler::molFromPickle(pickle_string, ret.get());
    return ret;
}

std::string lhasa::rdkit_mol_to_smiles(RDKit::ROMol& mol) {
    auto ret = RDKit::MolToSmiles(mol);
    return ret;
}


void lhasa::append_from_smiles(CootLigandEditorCanvas& canvas, std::string smiles) {
    canvas.append_molecule(rdkit_mol_from_smiles(smiles));
}

void lhasa::append_from_pickle_base64(CootLigandEditorCanvas& canvas, std::string base64_pickle_string) {
    std::string pickle_string = tinygltf::base64_decode(base64_pickle_string);
    auto appendee = rdkit_mol_from_pickle(pickle_string);
    auto smiles = rdkit_mol_to_smiles(*appendee.get());
    g_info("Smiles from pickle: %s -> %s", base64_pickle_string.c_str(), smiles.c_str());
    canvas.append_molecule(std::move(appendee));
}

std::unique_ptr<coot::ligand_editor_canvas::ActiveTool> lhasa::make_active_tool(emscripten::val tool) {
    using namespace coot::ligand_editor_canvas;

    // Just returns 'object'
    //std::string type_name = tool.typeOf().as<std::string>();

    std::string classname = tool["__proto__"]["constructor"]["name"].as<std::string>();
    // std::cout<< classname << '\n';

    if(classname == "ElementInsertion") {
        return std::make_unique<ActiveTool>(tool.as<ElementInsertion>());
    }
    if(classname == "BondModifier") {
        return std::make_unique<ActiveTool>(tool.as<BondModifier>());
    }
    if(classname == "TransformTool") {
        return std::make_unique<ActiveTool>(tool.as<TransformTool>());
    }
    if(classname == "StructureInsertion") {
        return std::make_unique<ActiveTool>(tool.as<StructureInsertion>());
    }
    if(classname == "FlipTool") {
        return std::make_unique<ActiveTool>(tool.as<FlipTool>());
    }
    if(classname == "DeleteTool") {
        return std::make_unique<ActiveTool>(tool.as<DeleteTool>());
    }
    if(classname == "ChargeModifier") {
        return std::make_unique<ActiveTool>(tool.as<ChargeModifier>());
    }
    if(classname == "GeometryModifier") {
        return std::make_unique<ActiveTool>(tool.as<GeometryModifier>());
    }
    if(classname == "FormatTool") {
        return std::make_unique<ActiveTool>(tool.as<FormatTool>());
    }
    if(classname == "RemoveHydrogensTool") {
        return std::make_unique<ActiveTool>(tool.as<RemoveHydrogensTool>());
    }

    g_critical("%s does not correspond to any known tool type. Returning empty ActiveTool.", classname.c_str());
    return std::make_unique<ActiveTool>();

}

coot::ligand_editor_canvas::ElementInsertion lhasa::element_insertion_from_symbol(std::string sym) {
    return coot::ligand_editor_canvas::ElementInsertion(sym.c_str());
}