/*
 * validation-graphs/validation-graphs.hh
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

#ifndef VALIDATION_GRAPHS_HH
#define VALIDATION_GRAPHS_HH

#include <gtk/gtk.h>
#include <type_traits>
#include <string>

namespace coot {

    enum class validation_graph_type : unsigned char {
        rota,
        temp_factor,
        density_fit,
        density_correlation,
        rama,
        omega,
        geometry,
        ncs
    };
    /// For using as unsigned char like so:
    /// `static_cast<validation_graph_type_repr_t>(enum_instance)`
    typedef std::underlying_type<validation_graph_type>::type validation_graph_type_repr_t;
    std::string validation_graph_type_to_human_name(validation_graph_type graph_type);
} // namespace coot

#endif // VALIDATION_GRAPHS_HH

