/*
 * src/cc-interface-alignment.hh
 *
 * Copyright 2018 by Medical Research Council
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

#include <string>


/*! \file
  \brief Coot Scripting Interface - Alignment utilities
*/

//!  \brief associate an alignment with a chain in a model molecule
//!
//! The pir_alignment is a string (with newlines)
//!
void associate_pir_alignment(int imol, std::string chain_id, std::string pir_alignment);

//!  \brief associate an alignment in a file with a chain in a model molecule
//!
void associate_pir_alignment_from_file(int imol, std::string chain_id, std::string pir_alignment_file_name);

//!  \brief apply the mutations of the associated alignment
//!
void apply_pir_alignment(int imol, std::string chain_id);

