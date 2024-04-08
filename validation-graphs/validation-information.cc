/*
 * validation-graphs/validation-information.cc
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
#include "validation-information.hh"


unsigned int coot::validation_information_t::get_index_for_chain(const std::string &chain_id) {

   for (unsigned int i=0; i<cviv.size(); i++) {
    if (chain_id == cviv[i].chain_id)
        return i;
    }
    chain_validation_information_t cvi(chain_id);
    cviv.push_back(cvi);
    return cviv.size() -1;
}

bool coot::should_hang_down(coot::graph_data_type type) {

    using ty = coot::graph_data_type;
    switch (type) {
        // case ty::Correlation: {
        //     return true;
        // }
        default: {
            return false;
        }
    }
}

bool is_probability_plot(coot::graph_data_type type) {
    using ty = coot::graph_data_type;
    switch (type) {
        case ty::LogProbability:
        case ty::Probability: {
            return true;
        }
        default: {
            return false;
        }
    }
}

// void coot::validation_information_t::add_residue_valiation_informtion(
//     const residue_validation_information_t &rvi, 
//     const std::string &chain_id) {

//     unsigned int idx = get_index_for_chain(chain_id);
//     cviv[idx].add_residue_validation_information(rvi);
// }
