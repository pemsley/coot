/* src/graphics-ligand-view.hh
 * 
 * Copyright 2011 by The University of Oxford
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
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include <iostream>
#include "lidia-core/lig-build.hh"
#include "lidia-core/lbg-molfile.hh"
#include "geometry/protein-geometry.hh"
#include "coot-render.hh"  // for colour_t

class graphics_ligand_atom : public lig_build::atom_t {
   void bitmap_text(const std::string &s) const;
public:
   graphics_ligand_atom(lig_build::pos_t pos_in, std::string ele_in, int formal_charge_in) :
      lig_build::atom_t(pos_in, ele_in, formal_charge_in) {
      font_colour = "";
   }
   std::string font_colour;

   void make_text_item(const lig_build::atom_id_info_t &atom_id_info_in,
		       const coot::colour_t &fc) const;
   coot::colour_t get_colour(bool against_a_dark_background) const;
   
};

class graphics_ligand_bond : public lig_build::bond_t {
public:
   graphics_ligand_bond(int first, int second, lig_build::bond_t::bond_type_t type) :
      lig_build::bond_t(first, second, type) {}
   void gl_bond(const lig_build::pos_t &pos_1, const lig_build::pos_t &pos_2,
		bool shorten_first, bool shorten_second, lig_build::bond_t::bond_type_t bt);
   void gl_bond_double_aromatic_bond(const lig_build::pos_t &pos_1, const lig_build::pos_t &pos_2,
				     bool shorten_first, bool shorten_second);
   void gl_bond_double_bond(const lig_build::pos_t &pos_1, const lig_build::pos_t &pos_2, bool shorten_first, bool shorten_second);
};

// for graphics ligand view (bottom left)
// 
class graphics_ligand_molecule : public lig_build::molecule_t<graphics_ligand_atom,
							      graphics_ligand_bond> {

   GLuint display_list_tag;
   void init_from_molfile_molecule(const lig_build::molfile_molecule_t &mol,
				   bool dark_background_flag = true);
   std::pair<bool, double> scale_correction; // push down to base class at some stage
                                             // (wmolecule has one of these too).
   void gl_bonds(bool against_a_dark_background);

public:
   graphics_ligand_molecule() {
      display_list_tag = 0;
      imol = -1; // unset
   } 
   ~graphics_ligand_molecule();
   // some OpenGL stuff.
   void render();
   void generate_display_list(bool against_a_dark_background);
   bool setup_from(int imol, mmdb::Residue *r,
		   const std::string &alt_conf,
		   coot::protein_geometry *geom_p, bool against_a_dark_background);
   int imol;
};
