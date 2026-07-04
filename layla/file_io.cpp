/* layla/file_io.cpp
 *
 * Copyright 2026 by Global Phasing Ltd.
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
#include "file_io.hpp"
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>
#include <rdkit/GraphMol/FileParsers/FileParsers.h>
#include <rdkit/GraphMol/FileParsers/FileWriters.h>
#include <rdkit/RDGeneral/RDConfig.h>
#ifdef RDK_BUILD_INCHI_SUPPORT
#include <rdkit/GraphMol/inchi.h>
#endif
#include <rdkit/GraphMol/chemdraw.h>
#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace coot::layla::io {

std::unique_ptr<RDKit::RWMol> mol_from_string(const std::string& data, CheminformaticsFileFormat format) {
    std::unique_ptr<RDKit::RWMol> mol;
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
            // Only the first molecule is imported (matches Lhasa's behaviour).
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
    return mol;
}

std::string mol_to_string(const RDKit::ROMol& mol, CheminformaticsFileFormat format) {
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

std::unique_ptr<RDKit::RWMol> rdkit_mol_from_smiles(const std::string& smiles) {
    std::unique_ptr<RDKit::RWMol> ret(RDKit::SmilesToMol(smiles, 0, false));
    return ret;
}

std::string rdkit_mol_to_smiles(const RDKit::ROMol& mol) {
    return RDKit::MolToSmiles(mol, true);
}

std::optional<CheminformaticsFileFormat> format_from_extension(const std::string& path) {
    auto dot = path.find_last_of('.');
    if (dot == std::string::npos) {
        return std::nullopt;
    }
    std::string ext = path.substr(dot + 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c){ return std::tolower(c); });
    if (ext == "mol")   return CheminformaticsFileFormat::Molfile;
    if (ext == "mdf")   return CheminformaticsFileFormat::Molfile;
    if (ext == "sdf")   return CheminformaticsFileFormat::SDF;
    if (ext == "inchi") return CheminformaticsFileFormat::InChI;
    if (ext == "cdxml") return CheminformaticsFileFormat::CDXML;
    return std::nullopt;
}

std::unique_ptr<RDKit::RWMol> mol_from_file(const std::string& path,
                                            std::optional<CheminformaticsFileFormat> format) {
    if (!format) {
        format = format_from_extension(path);
        if (!format) {
            throw std::runtime_error("Could not determine file format from extension: " + path);
        }
    }
    std::ifstream in(path, std::ios::binary);
    if (!in) {
        throw std::runtime_error("Could not open file for reading: " + path);
    }
    std::ostringstream buffer;
    buffer << in.rdbuf();
    return mol_from_string(buffer.str(), *format);
}

void mol_to_file(const RDKit::ROMol& mol, const std::string& path,
                 std::optional<CheminformaticsFileFormat> format) {
    if (!format) {
        format = format_from_extension(path);
        if (!format) {
            throw std::runtime_error("Could not determine file format from extension: " + path);
        }
    }
    std::string data = mol_to_string(mol, *format);
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
        throw std::runtime_error("Could not open file for writing: " + path);
    }
    out << data;
    if (!out) {
        throw std::runtime_error("Failed while writing file: " + path);
    }
}

} // namespace coot::layla::io
