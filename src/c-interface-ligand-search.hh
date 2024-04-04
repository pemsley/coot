/*
 * src/c-interface-ligand-search.hh
 *
 * Copyright 2016 by Medical Research Council
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

#ifndef C_INTERFACE_LIGAND_SEARCH_HH
#define C_INTERFACE_LIGAND_SEARCH_HH

#include "ligand/wligand.hh"
#include "c-interface-ligands-widgets.hh"

std::vector<int> execute_ligand_search_internal(coot::wligand *wlig_p);
ligand_wiggly_ligand_data_t ligand_search_install_wiggly_ligands();

#endif // C_INTERFACE_LIGAND_SEARCH_HH

