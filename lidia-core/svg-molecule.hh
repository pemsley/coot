/*
 * lidia-core/svg-molecule.hh
 *
 * Copyright 2022 by Medical Research Council
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
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 */

#ifndef SVG_MOLECULE_HH
#define SVG_MOLECULE_HH

#include "lig-build.hh"
#include "lbg-shared.hh"
#include "use-rdkit.hh"
#include "svg-container.hh"

class svg_atom_t : public lig_build::atom_t {
   // use self element to set the colour
   std::string colour;
public:
   svg_atom_t(const lig_build::pos_t &pos_in, const std::string &ele_in, int formal_charge_in) :
      lig_build::atom_t(pos_in, ele_in, formal_charge_in) { set_colour(); solvent_accessibility = -1; }

   std::vector<coot::bash_distance_t> bash_distances;
   double solvent_accessibility;
   double get_solvent_accessibility() const { return solvent_accessibility; }
   void add_solvent_accessibility(double sa) {
      solvent_accessibility = sa;
   }

   void set_colour(bool against_a_dark_background=false);
   std::string make_text_item(const lig_build::atom_id_info_t &atom_id_info_in,
                              const lig_build::pos_t &centre,
                              double scale, double median_bond_length) const;
};

class svg_bond_t : public lig_build::bond_t {
   std::string draw_sheared_or_darted_wedge_bond(const lig_build::pos_t &pos_1,
                                                 const lig_build::pos_t &pos_2,
                                                 const std::string &bond_colour,
                                                 const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
                                                 const lig_build::pos_t &centre,
                                                 double scale) const;
   std::string make_bond_line_string(const lig_build::pos_t &p1, const lig_build::pos_t &p2,
                                     const std::string &bond_colour) const;
   std::string make_dashed_bond_line_string(const lig_build::pos_t &p1, const lig_build::pos_t &p2,
                                            const std::string &bond_colour) const;
public:
   svg_bond_t(int first, int second, lig_build::bond_t::bond_type_t type) :
      lig_build::bond_t(first, second, type) {}

   std::string draw_bond(const svg_atom_t &at_1, const svg_atom_t &at_2,
                  bool at_1_in_ring_flag, bool at_2_in_ring_flag,
                  lig_build::bond_t::bond_type_t bt,
                  const std::string &bond_colour,
                  bool shorten_first, bool shorten_second,
                  const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
                  const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
                  const lig_build::pos_t &centre, double scale);
   std::string draw_double_in_ring_bond(const lig_build::pos_t &pos_1,
                                        const lig_build::pos_t &pos_2,
                                        const std::string &bond_colour,
                                        bool shorten_first,
                                        bool shorten_second,
                                        const lig_build::pos_t &centre,
                                        double scale, bool dashed_inner=false);
   std::string draw_double_bond(const lig_build::atom_t &at_1,
                                const lig_build::atom_t &at_2,
                                const std::string &bond_colour,
                                bool shorten_first, bool shorten_second,
                                const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_first_atom,
                                const std::vector<std::pair<lig_build::atom_t, lig_build::bond_t> > &other_connections_to_second_atom,
                                const lig_build::pos_t &centre, double scale);
};

class svg_molecule_t : public lig_build::molecule_t<svg_atom_t, svg_bond_t> {
public:
   svg_molecule_t() { median_bond_length_ = 1.0; }
   void import_rdkit_mol(RDKit::ROMol *mol, int iconf);
   std::string render_to_svg_string(double scale_factor, bool dark_background_flag);
   svg_container_t make_svg(double scale_factor, bool dark_background_flag);
   double median_bond_length_;
   double get_scale() const;

   static lig_build::pos_t mol_coords_to_svg_coords(const lig_build::pos_t &pos_1,
                                                    const lig_build::pos_t &centre,
                                                    double scale);
};


#endif // SVG_MOLECULE_HH
