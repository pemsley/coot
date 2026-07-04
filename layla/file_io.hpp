/* layla/file_io.hpp
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

#ifndef LAYLA_FILE_IO_HPP
#define LAYLA_FILE_IO_HPP

#include <rdkit/GraphMol/RWMol.h>
#include <memory>
#include <optional>
#include <string>

// Cheminformatics file I/O shared by Layla (native GTK) and Lhasa (WASM).
// Deliberately free of GTK, emscripten and the ligand-editor canvas: it only
// converts between RDKit molecules and strings/files. Each front-end handles
// its own canvas append/get and its own way of reading/writing bytes.
namespace coot::layla::io {

enum class CheminformaticsFileFormat {
    Molfile,
    SDF,
    InChI,
    CDXML
};

// Molecule <-> in-memory string in the given format. Both throw std::exception
// on parse/serialization failure (mol_from_string never returns nullptr).
std::unique_ptr<RDKit::RWMol> mol_from_string(const std::string& data, CheminformaticsFileFormat format);
std::string mol_to_string(const RDKit::ROMol& mol, CheminformaticsFileFormat format);

std::unique_ptr<RDKit::RWMol> rdkit_mol_from_smiles(const std::string& smiles);
std::string rdkit_mol_to_smiles(const RDKit::ROMol& mol);

// Native file helpers (unused by the WASM front-end). If `format` is nullopt the
// format is inferred from the path's extension. Both throw on failure.
std::optional<CheminformaticsFileFormat> format_from_extension(const std::string& path);
std::unique_ptr<RDKit::RWMol> mol_from_file(const std::string& path,
                                            std::optional<CheminformaticsFileFormat> format = std::nullopt);
void mol_to_file(const RDKit::ROMol& mol, const std::string& path,
                 std::optional<CheminformaticsFileFormat> format = std::nullopt);

} // namespace coot::layla::io

#endif
