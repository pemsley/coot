/*
 * validation-graphs/validation-graphs.cc
 *
 * Copyright 2023 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */
#include "validation-graphs.hh"

std::string coot::validation_graph_type_to_human_name(coot::validation_graph_type graph_type) {
    switch (graph_type) {
        case validation_graph_type::density_fit: {
            return "Density fit";
            break;
        }
        case validation_graph_type::geometry: {
            return "Geometry analysis";
            break;
        }
        case validation_graph_type::ncs: {
            return "NCS";
            break;
        }
        case validation_graph_type::omega: {
            return "Omega";
            break;
        }
        case validation_graph_type::rama: {
            return "Ramachandran";
            break;
        }
        case validation_graph_type::rota: {
            return "Rotamer analysis";
            break;
        }
        case validation_graph_type::temp_factor: {
            return "Temp factor";
            break;
        }
        default: {
            return "unknown";
            break;
        }
    }
}
