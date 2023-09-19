/* layla/generators.hpp
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

#ifndef LAYLA_GENERATORS_HPP
#define LAYLA_GENERATORS_HPP
#include <string>
#include <optional>
#include <variant>
#include <vector>
#include <gio/gio.h>
#include "notifier.hpp"

namespace coot::layla {

struct AcedrgOptions {
    bool p;
    bool z;
};

struct Grade2Options {

};

struct GeneratorRequest {
    enum class InputFormat: unsigned char {
        SMILES,
        MolFile
    } input_format;
    enum class Generator: unsigned char {
        Acedrg,
        Grade2
    } generator;

    std::string monomer_id;
    std::string molecule_smiles;
    std::optional<std::string> executable_path;
    std::variant<Grade2Options, AcedrgOptions> generator_settings;

    std::string get_input_filename() const;
    std::string get_output_filename() const;
    std::vector<std::string> build_commandline() const;
};

inline GCancellable* global_generator_request_task_cancellable;

GCancellable* run_generator_request(GeneratorRequest request, CootLaylaNotifier* notifier);

} // namespace coot::layla

#endif // LAYLA_GENERATORS_HPP