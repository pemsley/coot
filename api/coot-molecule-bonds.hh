/*
 * api/coot-molecule-bonds.hh
 * 
 * Copyright 2020 by Medical Research Council
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
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifndef COOT_MOLECULE_BONDS_HH
#define COOT_MOLECULE_BONDS_HH

#include "coot-utils/simple-mesh.hh"
#include "coords/graphical-bonds-container.hh"

void make_graphical_bonds_cis_peptides(coot::simple_mesh_t &m, const graphical_bonds_container &gbc);


#endif // COOT_MOLECULE_BONDS_HH
