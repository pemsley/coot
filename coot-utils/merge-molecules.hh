/*
 * coot-utils/merge-molecules.hh
 *
 * Copyright 2021 by Medical Research Council
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
#ifndef MERGE_MOLECULES_HH
#define MERGE_MOLECULES_HH

#include <vector>
#include <mmdb2/mmdb_manager.h>

namespace coot {
   // only copy the first chain of the mol_otheres
   void merge_molecules(mmdb::Manager *mol, std::vector<mmdb::Manager *> mol_others);
}



#endif // MERGE_MOLECULES_HH
