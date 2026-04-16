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
#include <rdkit/GraphMol/FileParsers/FileWriters.h>
#include <rdkit/GraphMol/FileParsers/FileParsers.h>
#include <rdkit/RDGeneral/RDConfig.h>
#ifdef RDK_BUILD_INCHI_SUPPORT
#include <rdkit/GraphMol/inchi.h>
#endif
#include <rdkit/GraphMol/chemdraw.h>
#include "glog_replacement.hpp"
#include "../utils/base64-encode-decode.hh"



std::unique_ptr<RDKit::RWMol> lhasa::rdkit_mol_from_pickle(std::string pickle_string) {
    std::unique_ptr<RDKit::RWMol> ret = std::make_unique<RDKit::RWMol>();
    RDKit::MolPickler::molFromPickle(pickle_string, ret.get());
    return ret;
}

std::string lhasa::rdkit_mol_to_pickle_base64(const RDKit::ROMol& mol) {
    std::string pickle_string;
    unsigned int pickleFlags = RDKit::PicklerOps::AtomProps | RDKit::PicklerOps::BondProps | RDKit::PicklerOps::MolProps | RDKit::PicklerOps::CoordsAsDouble | RDKit::PicklerOps::AllProps;
    RDKit::MolPickler::pickleMol(mol, pickle_string, pickleFlags);
    return moorhen_base64::base64_encode((const unsigned char*) pickle_string.c_str(), pickle_string.size());
}

std::string lhasa::export_mol_to_pickle_base64(CootLigandEditorCanvas& canvas, unsigned int molecule_idx) {
    const auto& mol = canvas.get_rdkit_molecule(molecule_idx);
    return rdkit_mol_to_pickle_base64(mol);
}

unsigned int lhasa::append_from_pickle_base64(CootLigandEditorCanvas& canvas, std::string base64_pickle_string) {
    std::string pickle_string = moorhen_base64::base64_decode(base64_pickle_string);
    auto appendee = rdkit_mol_from_pickle(pickle_string);
    auto smiles = rdkit_mol_to_smiles(*appendee.get());
    g_info("Smiles from pickle: %s -> %s", base64_pickle_string.c_str(), smiles.c_str());
    return canvas.append_molecule(std::move(appendee));
}

std::string lhasa::rdkit_mol_to_smiles(const RDKit::ROMol& mol) {
    auto ret = RDKit::MolToSmiles(mol, true);
    return ret;
}

std::unique_ptr<RDKit::RWMol> lhasa::rdkit_mol_from_smiles(std::string smiles) {
    std::unique_ptr<RDKit::RWMol> ret(RDKit::SmilesToMol(smiles, 0, false));
    return ret;
}

unsigned int lhasa::append_from_smiles(CootLigandEditorCanvas& canvas, std::string smiles) {
    return canvas.append_molecule(rdkit_mol_from_smiles(smiles));
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


std::string lhasa::export_mol(CootLigandEditorCanvas& canvas, unsigned int molecule_idx, lhasa::CheminformaticsFileFormat format) {
    const auto& mol = canvas.get_rdkit_molecule(molecule_idx);
    switch(format) {
        case CheminformaticsFileFormat::Molfile: {
            return RDKit::MolToMolBlock(mol);
        }
        case CheminformaticsFileFormat::SDF: {
            return RDKit::MolToMolBlock(mol) + "$$$$\n";
        }
        case CheminformaticsFileFormat::InChI: {
#ifdef RDK_BUILD_INCHI_SUPPORT
            RDKit::ExtraInchiReturnValues rv;
            return RDKit::MolToInchi(mol, rv);
#else
            throw std::runtime_error("RDKit was built without InChI support");
#endif
        }
        case CheminformaticsFileFormat::CDXML: {
            return RDKit::v2::MolToChemDrawBlock(mol);
        }
        default: {
            throw std::runtime_error("Unknown file format");
        }
    }
}

unsigned int lhasa::append_from_import(CootLigandEditorCanvas& canvas, std::string data, lhasa::CheminformaticsFileFormat format) {
    std::shared_ptr<RDKit::RWMol> mol;
    switch(format) {
        case CheminformaticsFileFormat::Molfile:
        case CheminformaticsFileFormat::SDF: {
            mol = RDKit::v2::FileParsers::MolFromMolBlock(data);
            break;
        }
        case CheminformaticsFileFormat::InChI: {
#ifdef RDK_BUILD_INCHI_SUPPORT
            RDKit::ExtraInchiReturnValues rv;
            mol.reset(RDKit::InchiToMol(data, rv));
#else
            throw std::runtime_error("RDKit was built without InChI support");
#endif
            break;
        }
        case CheminformaticsFileFormat::CDXML: {
            auto mols = RDKit::v2::MolsFromChemDrawBlock(data);
            if (mols.empty()) {
                throw std::runtime_error("No molecules found in CDXML data");
            }
            if (mols.size() > 1) {
                g_warning("Multiple molecules found in CDXML data. Only the first one will be imported.");
            }
            /// Hmmm, maybe we should import all the molecules, not just the first one?
            mol = std::move(mols.front());
            break;
        }
        default: {
            throw std::runtime_error("Unknown file format");
        }
    }
    if (!mol) {
        throw std::runtime_error("Failed to parse molecule from input data");
    }
    return canvas.append_molecule(std::move(mol));
}