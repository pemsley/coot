/*
 * src/animated-ligand.hh
 *
 * Copyright 2015 by Medical Research Council
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

#include "pli/flev-annotations.hh"
#include "gl-bits.hh"

namespace coot {

   class animated_ligand_interactions_t : public pli::fle_ligand_bond_t {
   public:
      animated_ligand_interactions_t(const fle_ligand_bond_t &lb) :
         fle_ligand_bond_t(lb) { }
      void draw(mmdb::Manager *mol,
                const gl_context_info_t &gl_info,
                const long &start_time) const;
   };

   // encapsulate this into molecule_class_info_t?
   //
   // trivial helper class to get specs and distance for atoms when a
   // link is made.
   //
   class dict_link_info_t {
      bool check_for_order_switch(mmdb::Residue *residue_ref,
                                  mmdb::Residue *residue_new,
                                  const std::string &link_type,
                                  const protein_geometry &geom) const;
   public:
      // this can throw a std::runtime_error
      dict_link_info_t (mmdb::Residue *residue_ref, mmdb::Residue *residue_new,
                        const std::string &link_type, const protein_geometry &geom);
      atom_spec_t spec_ref;
      atom_spec_t spec_new;
      double dist;
   };
}
