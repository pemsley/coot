// -*-c++-*-
/* coords/Bond_lines.h
 *
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2016 by Medical Research Council
 * Author: Paul Emsley
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


// OK, the thinking here is that molecule_class_info_t functions
// create Bond_lines_containers via various constructors and perhaps
// other manipulations.  These then get converted to a
// graphical_bonds_container, and the drawing routine knows how to
// iterate through that to make the right coloured lines - that
// routine uses display_bonds(), which operates on a
// graphical_bonds_container bonds_box (the drawing function also
// takes a bond_width)
//
//

#ifndef BOND_LINES_H
#define BOND_LINES_H

#include <vector>
#include <string>

#include "geometry/protein-geometry.hh"
#include "Cartesian.hh"
#include "phenix-geo.hh"

#include "coot-utils/coot-rama.hh" // for ramachandran scoring of intermediate atoms
#include "ramachandran-container.hh"
#include "rotamer-container.hh"
#include "ligand/rotamer.hh"
#include "coot-utils/coot-coord-utils.hh" // is this needed?

#include "geometry/bonded-quad.hh"
#include "mmdb-crystal.hh"

namespace coot {

   static std::string b_factor_bonds_scale_handle_name;

   enum coords_bond_colour_t { COLOUR_BY_CHAIN=0,
      COLOUR_BY_CHAIN_C_ONLY=20,
      COLOUR_BY_CHAIN_GOODSELL=21,
      COLOUR_BY_ATOM_TYPE=1,
      COLOUR_BY_SEC_STRUCT=2,
      DISULFIDE_COLOUR=3,
      COLOUR_BY_MOLECULE=4,
      COLOUR_BY_RAINBOW=5,
      COLOUR_BY_OCCUPANCY=6,
      COLOUR_BY_B_FACTOR=7,
      COLOUR_BY_USER_DEFINED_COLOURS=8,
      COLOUR_BY_HYDROPHOBIC_SIDE_CHAIN=9 };

   enum hydrophobic_side_chain_t {
                                  HYDROPHOBIC_TYPE_MAIN_CHAIN,
                                  HYDROPHOBIC_TYPE_HYDROPHOBIC,
                                  HYDROPHOBIC_TYPE_HYDROPHILIC };

   hydrophobic_side_chain_t get_type(mmdb::Residue *residue_p);

   class my_atom_colour_map_t {
   public:
      my_atom_colour_map_t() {
         atom_colour_map.resize(50, "---");
      }
      std::vector<std::string> atom_colour_map;

      unsigned int index_for_chain(const std::string &chain_id);
      // These colours ranges need to be echoed in the GL bond drawing
      // routine.
      int index_for_rainbow(float wheel_colour) {
         return int(30.0*wheel_colour);
      }
      int index_for_occupancy(float wheel_colour) {
         return int(5.0*wheel_colour);
      }
      int index_for_b_factor(float fraction) {
         const unsigned n_b_factor_colours = 48;
         return static_cast<int>(static_cast<float>(n_b_factor_colours) * fraction); // 48 so that we don't overlap with CA colours
      }

      void fill_chain_id_map(const atom_selection_container_t &SelAtom);

   };

   class model_bond_atom_info_t {
      std::vector<mmdb::PAtom> hydrogen_atoms_;
      std::vector<mmdb::PAtom> non_hydrogen_atoms_;
   public:
      mmdb::PPAtom     Hydrogen_atoms() const;
      mmdb::PPAtom non_Hydrogen_atoms() const;
      int n_H() const { return hydrogen_atoms_.size(); }
      int n_non_H() const { return non_hydrogen_atoms_.size(); }
      void add_atom(mmdb::Atom *atom) {
         std::string element = atom->element;
         if (element == " H" || element == " D") {
            hydrogen_atoms_.push_back(atom);
         } else {
            non_hydrogen_atoms_.push_back(atom);
         }
      }
   };
}

#include "graphics-line.hh"

#include "graphical-bonds-container.hh"

// Bond_lines is a container class, containing a colour index
// and a vector of pairs of coot::Cartesians that should be drawn
// _in_ that colour.
//
class Bond_lines {
   int colour;
   std::vector<graphics_line_t> points;

 public:
   explicit Bond_lines(const graphics_line_t &pts);
   Bond_lines() { colour = 0; }
   explicit Bond_lines(int col) { colour = col; }

   void add_bond(const coot::CartesianPair &p,
                 graphics_line_t::cylinder_class_t cc,
                 bool begin_end_cap,
                 bool end_end_cap,
                 int model_number_in,
                 int atom_index_1, int atom_index_2);
   int size() const;

   // return the coordinates of the start and finish points of the i'th bond.
   //
   const coot::Cartesian &GetStart(unsigned int i) const;
   const coot::Cartesian &GetFinish(unsigned int i) const;
   const graphics_line_t &operator[](unsigned int i) const;
   void update(mmdb::Atom **atom_selection, int n_atoms);
};

//
enum symm_keys {NO_SYMMETRY_BONDS};

class Bond_lines_container {

   bool verbose_reporting;
   bool do_disulfide_bonds_flag;
   bool do_bonds_to_hydrogens;
   int udd_has_ca_handle;
   float b_factor_scale;
   bool for_GL_solid_model_rendering;
   bool do_sticks_for_waters;
   int n_atoms_in_atom_selection; // for fast not-in-no-bonds-to-these-atoms check

   // we rely on SelAtom.atom_selection being properly constucted to
   // contain all atoms
   //
   // if model_number is 0, display all models.
   //
   void construct_from_asc(const atom_selection_container_t &SelAtom,
                           int imol,
                           float min_dist, float max_dist,
                           int atom_colour_type,
                           short int is_from_symmetry_flag,
                           bool draw_missing_loops_flag,
                           int model_number,
                           bool do_rama_markup=false,
                           bool do_rota_markup=false);

   // PDBv3 FIXME
   bool is_hydrogen(const std::string &ele) const {
      if (ele == " H" || ele == " D")
         return true;
      else
         return false;
   }
   bool is_deuterium(const std::string &ele) const {
      if (ele == " D")
         return true;
      else
         return false;
   }

   void construct_from_atom_selection(const atom_selection_container_t &asc,
                                      const mmdb::PPAtom atom_selection_1,
                                      int n_selected_atoms_1,
                                      const mmdb::PPAtom atom_selection_2,
                                      int n_selected_atoms_2,
                                      int imol,
                                      float min_dist, float max_dist,
                                      int atom_colour_type,
                                      bool are_different_atom_selections,
                                      bool have_udd_atoms,
                                      int udd_handle);
   void atom_selection_missing_loops(const atom_selection_container_t &asc,
                                     int udd_atom_index_handle,
                                     int udd_fixed_during_refinement_handle);

   void construct_from_model_links(mmdb::Model *model, int udd_atom_index_handle,
                                   int udd_user_defined_atom_colour_index_handle,
                                   int atom_colour_type);
   // which wraps...
   void add_link_bond(mmdb::Model *model_p, int udd_atom_index_handle, int udd_user_defined_atom_colour_index_handle,
                      int atom_colour_type, mmdb::Link *link);
   void add_link_bond(mmdb::Model *model_p, int udd_atom_index_handle, int udd_user_defined_atom_colour_index_handle,
                      int atom_colour_type, mmdb::LinkR *linkr);

   template<class T> void add_link_bond_templ(mmdb::Model *model_p, int udd_atom_index_handle,
                                              int udd_user_defined_atom_colour_index_handle,
                                              int atom_colour_type, T *link);

   // now wit optional arg.  If atom_colour_type is set, then use/fill
   // it to get colour indices from chainids.
   void handle_MET_or_MSE_case (mmdb::PAtom mse_atom, int udd_handle,
                                int udd_handle_for_atom_index,
                                int udd_user_defined_atom_colour_index_handle,
                                int atom_colour_type,
                                coot::my_atom_colour_map_t *atom_colour_map = 0);
   void handle_long_bonded_atom(mmdb::PAtom atom,
                                int udd_handle_bond,
                                int udd_atom_index_handle,
                                int udd_user_defined_atom_colour_index_handle,
                                int atom_colour_type,
                                coot::my_atom_colour_map_t *atom_colour_map = 0);

   // void check_atom_limits(atom_selection_container_t SelAtoms) const;

   void write(std::string) const;

   mmdb::PPAtom trans_sel(atom_selection_container_t AtomSel,
                          const std::pair<symm_trans_t, Cell_Translation> &symm_trans) const;

   void add_zero_occ_spots(const atom_selection_container_t &SelAtom);
   void add_deuterium_spots(const atom_selection_container_t &SelAtom);
   void add_ramachandran_goodness_spots(const atom_selection_container_t &SelAtom);
   void add_rotamer_goodness_markup(const atom_selection_container_t &SelAtom);
   void add_atom_centres(int imol,
                         const atom_selection_container_t &SelAtom,
                         int atom_colour_type,
                         coot::my_atom_colour_map_t *atom_colour_map = 0);

   int add_ligand_bonds(const atom_selection_container_t &SelAtom, int imol,
                        mmdb::PPAtom ligand_atoms_selection, int n_ligand_atoms);
   std::vector<rotamer_markup_container_t> dodecs;
   // filled by return value of
   std::vector<rotamer_markup_container_t> get_rotamer_dodecs(const atom_selection_container_t &SelAtom) const;
   // which uses:

   // partially fill dodecs
   static void
   add_rotamer_markups(const std::vector<unsigned int> &indices,
                       const std::vector<std::pair<mmdb::Residue *, mmdb::Atom *> > &residues,
                       coot::rotamer_probability_tables *rpt,
                       std::vector<rotamer_markup_container_t> *dodecs);
   // and
   static
   rotamer_markup_container_t get_rotamer_probability(const std::pair<mmdb::Residue *, mmdb::Atom *> &ra,
                                                      coot::rotamer_probability_tables *rpt);


   void add_cis_peptide_markup(const atom_selection_container_t &SelAtom, int model_number);

   // no longer use this - it's too crude
   bool draw_these_residue_contacts(mmdb::Residue *this_residue, mmdb::Residue *env_residue,
                                    coot::protein_geometry *protein_geom);

   // we want to filter out atom contacts along the main chain.  Previously we did that by
   // checking that the residues were not next to each other (above) - but I want to see contacts
   // between bases in DNA, so now, filter out distances based on atom names (and residue numbering)
   //
   bool draw_these_atom_contacts(mmdb::Atom *this_residue, mmdb::Atom *env_residue,
                                 coot::protein_geometry *protein_geom); // protein_geom is modifiable

   // abstract this out of construct_from_atom_selection for cleanliness.
   //
   void mark_atoms_as_bonded(mmdb::Atom *atom_p_1, mmdb::Atom *atom_p_2, bool have_udd_atoms, int udd_handle, bool done_bond_udd_handle) const;


   void add_half_bonds(const coot::Cartesian &atom_1,
                       const coot::Cartesian &atom_2,
                       mmdb::Atom *at_1,
                       mmdb::Atom *at_2,
                       graphics_line_t::cylinder_class_t cc,
                       int model_number,
                       int atom_index_1,
                       int atom_index_2,
                       int atom_colour_type,
                       int udd_user_defined_atom_colour_index_handle,
                       coot::my_atom_colour_map_t *atom_colour_map_p,
                       bool add_begin_end_cap,
                       bool add_end_end_cap);

   // double and delocalized bonds (default (no optional arg) is double).
   // We pass udd_atom_index_handle because we need the atom index (not residue atom index) for
   // using no_bonds_to_these_atoms
   void add_double_bond(int imol, int imodel, int iat_1, int iat_2, mmdb::PPAtom atoms, int n_atoms,
                        int atom_colour_type, coot::my_atom_colour_map_t *atom_colour_map_p,
                        int udd_atom_index_handle,
                        int udd_user_defined_atom_colour_index_handle,
                        const std::vector<coot::dict_bond_restraint_t> &bond_restraints,
                        bool is_deloc=false);
   // used by above, can throw an exception
   clipper::Coord_orth get_neighb_normal(int imol, int iat_1, int iat_2, mmdb::PPAtom atoms, int n_atoms,
                                         bool also_2nd_order_neighbs=0) const;
   void add_triple_bond(int imol, int imodel, int iat_1, int iat_2, mmdb::PPAtom atoms, int n_atoms,
                        int atom_colour_type,
                        coot::my_atom_colour_map_t *atom_colour_map_p,
                        int udd_atom_index_handle,
                        int udd_user_defined_atom_colour_index_handle,
                        const std::vector<coot::dict_bond_restraint_t> &bond_restraints);


   bool have_dictionary;
   bool use_deuteranomaly_mode;
   const coot::protein_geometry *geom;
   enum { NOT_HALF_BOND, HALF_BOND_FIRST_ATOM, HALF_BOND_SECOND_ATOM };

 protected:
   std::vector<Bond_lines> bonds;
   std::vector<coot::Cartesian>  zero_occ_spots;
   std::vector<coot::Cartesian>  bad_CA_CA_dist_spots;
   std::vector<coot::Cartesian>  deuterium_spots;
   std::vector<std::pair<coot::Cartesian, coot::util::phi_psi_t> >  ramachandran_goodness_spots;
   std::vector<graphical_bonds_atom_info_t>  atom_centres;
   std::vector<int>        atom_centres_colour;
   void addBond(int colour, const coot::Cartesian &first, const coot::Cartesian &second,
                graphics_line_t::cylinder_class_t cc,
                int model_number,
                int atom_index_1,
                int atom_index_2,
                bool add_begin_end_cap = true,
                bool add_end_end_cap = true);
   void addBondtoHydrogen(const coot::Cartesian &first, const coot::Cartesian &second);
   // void add_deloc_bond_lines(int colour, const coot::Cartesian &first, const coot::Cartesian &second,
   // int deloc_half_flag);
   void add_dashed_bond(int col,
                        const coot::Cartesian &start,
                        const coot::Cartesian &end,
                        int half_bond_type_flag,
                        graphics_line_t::cylinder_class_t cc,
                        int model_number,
                        int atom_index_1, int atom_index_2);
   void addAtom(int colour, const coot::Cartesian &pos);
   int atom_colour(mmdb::Atom *at, int bond_colour_type,
                   int udd_user_defined_atom_colour_index_handle,
                   coot::my_atom_colour_map_t *atom_colour_map = 0);
   void bonds_size_colour_check(int icol) {
      int bonds_size = bonds.size();
      if (icol >= bonds_size)
         bonds.resize(icol+1);
   }

   // return the UDD handle
   int set_rainbow_colours(mmdb::Manager *mol);
   int set_b_factor_colours(mmdb::Manager *mol);

   void do_colour_by_chain_bonds_carbons_only(const atom_selection_container_t &asc,
                                              int imol,
                                              bool draw_missing_loops_flag,
                                              int atom_colour_type, // C-only or goodsell
                                              int draw_hydrogens_flag);
   void bond_by_distance(const atom_selection_container_t &asc, int imol, std::vector<mmdb::Residue *> &residues,
                         bool have_udd_atoms, int udd_found_bond_handle);
   void do_colour_by_chain_bonds_carbons_only_internals(int imol, int imodel,
                                                        int chain_idx,
                                                        mmdb::Atom *at1, mmdb::Atom *at2,
                                                        int iat_1, int iat_2,
                                                        std::vector<std::pair<bool, mmdb::Residue *> > *het_residues_p,
                                                        const std::string &element_1,
                                                        const std::string &element_2,
                                                        const coot::Cartesian &atom_1,
                                                        const coot::Cartesian &atom_2,
                                                        int atom_colour_type,
                                                        int uddHnd,
                                                        int udd_user_defined_atom_colour_index_handle);
   void do_colour_by_chain_bonds_internals_goodsell_mode(int imol, int imodel,
                                                         int chain_idx,
                                                         mmdb::Atom *at1, mmdb::Atom *at2,
                                                         int iat_1, int iat_2,
                                                         std::vector<std::pair<bool, mmdb::Residue *> > *het_residues_p,
                                                         const std::string &element_1,
                                                         const std::string &element_2,
                                                         const coot::Cartesian &atom_1,
                                                         const coot::Cartesian &atom_2,
                                                         int uddHnd,
                                                         int udd_user_defined_atom_colour_index_handle);

   void do_colour_by_ncs_related_chains_atoms_only(const atom_selection_container_t &asc,
                                                   int imol,
                                                   std::vector<std::vector<mmdb::Chain *> > ncs_related_chains,
                                                   bool change_c_only_flag, bool goodsell_mode);

   void do_colour_by_dictionary_and_by_chain_bonds(const atom_selection_container_t &asc,
                                                   int imol,
                                                   int draw_hydrogens_flag,
                                                   bool draw_missing_loops_flag,
                                                   short int change_c_only_flag,
                                                   bool do_goodsell_colour_mode,
                                                   bool do_rota_markup);

   void add_residue_monomer_bonds(const std::map<std::string, std::vector<mmdb::Residue *> > &residue_monomer_map,
                                  int imol, int model_number,
                                  int atom_colour_type,
                                  coot::my_atom_colour_map_t *atom_colour_map_p,
                                  int udd_atom_index_handle,
                                  int udd_bond_handle,
                                  int udd_user_defined_atom_colour_index_handle,
                                  int draw_hydrogens_flag,
                                  bool do_goodsell_colour_mode);

   void do_colour_by_dictionary_and_by_chain_bonds_carbons_only(const atom_selection_container_t &asc,
                                                                int imol,
                                                                int draw_hydrogens_flag,
                                                                bool draw_missing_loops_flag,
                                                                bool do_goodsell_colour_mode,
                                                                bool do_rota_markup);
   // and the bonds between the above monomers
   void add_polymer_bonds(const atom_selection_container_t &asc,
                          int atom_colour_type,
                          coot::my_atom_colour_map_t *atom_colour_map_p,
                          int draw_hydrogens_flag,
                          bool do_goodsell_colour_mode);
   void add_peptide_bonds(const atom_selection_container_t &asc,
                          int atom_colour_type,
                          coot::my_atom_colour_map_t *atom_colour_map_p,
                          int draw_hydrogens_flag,
                          bool do_goodsell_colour_mode);
   void add_phosphodiester_bonds(const atom_selection_container_t &asc,
                                 int atom_colour_type,
                                 coot::my_atom_colour_map_t *atom_colour_map_p,
                                 int draw_hydrogens_flag,
                                 bool do_goodsell_colour_mode);
   void add_carbohydrate_bonds(const atom_selection_container_t &asc, // oh. Tricky.
                               int atom_colour_type,
                               coot::my_atom_colour_map_t *atom_colour_map_p,
                               int draw_hydrogens_flag,
                               bool do_goodsell_colour_mode);
   void add_polymer_bonds_generic(const atom_selection_container_t &asc,
                                  int atom_colour_type,
                                  coot::my_atom_colour_map_t *atom_colour_map_p,
                                  int draw_hydrogens_flag,
                                  const std::string &res_1_atom_name, // in "res1"
                                  const std::string &res_2_atom_name, // in "res2"
                                  bool allow_het_group_link_bond,
                                  bool do_goodsell_colour_mode);
   void add_SS_bonds(const atom_selection_container_t &asc,
                     int atom_colour_type,
                     int draw_hydrogens_flag,
                     bool do_goodsell_colour_mode);
   void add_link_bonds(const atom_selection_container_t &asc,
                       int atom_colour_type,
                       int draw_hydrogens_flag,
                       bool do_goodsell_colour_mode);

   // the atoms have been added in order 0 is bonded to 1, 1 is bonded to 2, 2 is bonded to 3 etc.
   // and there is a double bond between 0 and 1, 2 and 3, and 4 to 5. Or maybe we could explicitly
   // add that to the the ring_atoms data.
   void draw_6_membered_ring(const std::string &residue_name,
                             const std::vector<mmdb::Atom *> &ring_atoms, int imodel,
                             int atom_colour_type, coot::my_atom_colour_map_t *atom_colour_map_p,
                             int udd_atom_index_handle,
                             int udd_user_defined_atom_colour_index_handle);
   // this calls the above function
   void draw_phenyl_ring_outer(mmdb::Residue *residue_p, int model_number,
                               int atom_colour_type, coot::my_atom_colour_map_t *atom_colour_map_p,
                               int udd_atom_index_handle, int udd_user_defined_atom_colour_index_handle);
   void draw_trp_rings(const std::vector<mmdb::Atom *> &ring_atoms, int imodel,
                       int atom_colour_type, coot::my_atom_colour_map_t *atom_colour_map_p,
                       int udd_atom_index_handle,
                       int udd_user_defined_atom_colour_index_handle);
   void draw_GA_rings(const std::vector<mmdb::Atom *> &ring_atoms, int imodel,
                      int atom_colour_type, coot::my_atom_colour_map_t *atom_colour_map_p,
                      int udd_atom_index_handle,
                      int udd_user_defined_atom_colour_index_handle);
   void draw_trp_ring_outer(mmdb::Residue *residue_p, int model_number,
                            int atom_colour_type, coot::my_atom_colour_map_t *atom_colour_map_p,
                            int udd_atom_index_handle,
                            int udd_user_defined_atom_colour_index_handle);
   void draw_GA_rings_outer(mmdb::Residue *residue_p, int model_number,
                            int atom_colour_type, coot::my_atom_colour_map_t *atom_colour_map_p,
                            int udd_atom_index_handle,
                            int udd_user_defined_atom_colour_index_handle);
   void draw_CUT_ring(mmdb::Residue *residue_p, int model_number,
                      int atom_colour_type, coot::my_atom_colour_map_t *atom_colour_map_p,
                      int udd_atom_index_handle, int udd_user_defined_atom_colour_index_handle);
   void draw_het_group_rings(mmdb::Residue *residue_p,
                             const std::vector<bonded_quad_atom_names> &bonded_quad_atom_names,
                             int model_number, int atom_colour_type,
                             coot::my_atom_colour_map_t *atom_colour_map,
                             int udd_atom_index_handle,
                             int udd_user_defined_atom_colour_index_handle);
   void draw_bonded_quad_atoms_rings(const std::vector<bonded_quad_atoms> &ring_atoms,
                                     int imodel, int atom_colour_type,
                                     coot::my_atom_colour_map_t *atom_colour_map_p,
                                     int udd_atom_index_handle,
                                     int udd_user_defined_atom_colour_index_handle);

   std::vector<std::pair<std::string, std::string> >
   get_aromatic_bonds(const coot::dictionary_residue_restraints_t &restraints) const;

   void try_set_b_factor_scale(mmdb::Manager *mol);
   graphical_bonds_container make_graphical_bonds_with_thinning_flag(bool thinning_flag) const;
   void add_bonds_het_residues(const std::vector<std::pair<bool, mmdb::Residue *> > &het_residues,
                               const atom_selection_container_t &sel_atoms,
                               int imol, int atom_colour_t, short int have_udd_atoms,
                               int udd_found_bond_handle, int udd_atom_index_handle,
                               int udd_user_defined_atom_colour_index_handle);
   void het_residue_aromatic_rings(mmdb::Residue *res, const coot::dictionary_residue_restraints_t &restraints,
                                   int udd_atom_index_handle, int col);
   // pass a list of atom name that are part of the aromatic ring system.
   void add_aromatic_ring_bond_lines(const std::vector<std::string> &ring_atom_names, mmdb::Residue *res,
                                     int udd_atom_index_handle, int col);
   bool invert_deloc_bond_displacement_vector(const clipper::Coord_orth &vect,
                                              int iat_1, int iat_2, mmdb::PPAtom residue_atoms, int n_atoms,
                                              const std::vector<coot::dict_bond_restraint_t> &bond_restraints) const;

   // add to het_residues maybe
   bool add_bond_by_dictionary_maybe(int imol, mmdb::Atom *atom_p_1,
                                     mmdb::Atom *atom_p_2,
                                     std::vector<std::pair<bool, mmdb::Residue *> > *het_residues);


   void do_colour_by_hydrophobic_side_chains(const atom_selection_container_t &asc,
                                             int imol,
                                             bool draw_missing_loops_flag,
                                             int draw_hydrogens_flag);

   std::vector<coot::util::cis_peptide_quad_info_t> cis_peptide_quads;

   // for user defined colours:
   //
   // return a colour index, and -1 on failure
   //
   int get_user_defined_col_index(mmdb::Atom *at, int udd_handle) const;

   void init();


public:

   enum bond_representation_type { COLOUR_BY_REGULAR_ATOM_MODE=601,
                                   COLOUR_BY_OCCUPANCY=602,
                                   COLOUR_BY_B_FACTOR=603,
                                   COLOUR_BY_USER_DEFINED_COLOURS=604
   };

   // Constructor A
   //
   // getting caught out with Bond_lines_container dependencies?
   // We need:  mmdb-extras.h which needs mmdb-manager.h and <string>
   // Bond_lines_container(const atom_selection_container_t &asc);
   explicit Bond_lines_container(const atom_selection_container_t &asc,
                                 int imol,
                                 int include_disulphides=0,
                                 int include_hydrogens=1,
                                 bool do_rama_markup=false,
                                 bool do_rota_markup=false,
                                 coot::rotamer_probability_tables *tables_p=0);

   // Constructor B
   //
   // the constructor for bond by dictionary - should use this most of the time.
   // geom_in can be null if you don't have it.
   //
   // if model_number is 0, display all models. If it is not 0 then
   // display only the given model_number (if possible, of course).
   //
   // 20220508-PE This is the constructor used by make_moving_atoms_graphics_object()
   //
   Bond_lines_container(const atom_selection_container_t &asc,
                        int imol,
                        const std::set<int> &no_bonds_to_these_atoms,
                        const coot::protein_geometry *geom_in,
                        int include_disulphides,
                        int include_hydrogens,
                        bool draw_missing_loops_flag,
                        int model_number,
                        std::string dummy,
                        bool do_rama_markup=false,
                        bool do_rota_markup=false,
                        bool do_sticks_for_waters=true,
                        coot::rotamer_probability_tables *tables_p=0);

   // Constructor C
   //
   Bond_lines_container(atom_selection_container_t asc, int imol,
                        float max_dist);

   // Constructor D
   //
   Bond_lines_container(atom_selection_container_t asc, int imol,
                        float min_dist, float max_dist);

   // Constructor E
   //
   // The constructor for ball and stick, this constructor implies that
   // for_GL_solid_model_rendering is set.
   //
   // geom_in can be null.
   Bond_lines_container (atom_selection_container_t asc,
                         int imol,
                         const coot::protein_geometry *geom);

   // Constructor F
   //
   Bond_lines_container(atom_selection_container_t SelAtom, coot::Cartesian point,
                        float symm_distance,
                        std::vector<symm_trans_t> symm_trans); // const & FIXME

   // Constructor G
   //
   // This finds bonds between a residue and the protein (in SelAtom).
   // It is used for the environment bonds box.
   //
   Bond_lines_container(const atom_selection_container_t &SelAtom,
                        mmdb::PPAtom residue_atoms,
                        int n_residue_atoms,
                        coot::protein_geometry *protein_geom, // modifiable, currently
                        bool residue_is_water_flag,
                        bool draw_env_distances_to_hydrogens_flag,
                        float min_dist,
                        float max_dist);

   // Constructor H
   //
   // same as above, except this does symmetry contacts too.
   //
   Bond_lines_container(const atom_selection_container_t &SelAtom,
                        mmdb::PPAtom residue_atoms,
                        int n_residue_atoms,
                        float min_dist,
                        float max_dist,
                        bool draw_env_distances_to_hydrogens_flag,
                        short int do_symmetry);

   // Constructor I
   //
   // This is the one for user-defined, occupancy and B-factor representation
   //
   Bond_lines_container (const atom_selection_container_t &SelAtom, int imol,
                         const coot::protein_geometry *protein_geom,
                         bond_representation_type br_type);

   // Constructor J
   //
   explicit Bond_lines_container(int col);

   // Constructor K
   //
   explicit Bond_lines_container(symm_keys key);

   // Constructor L
   //
   // Used by various CA-mode bonds
   //
   explicit Bond_lines_container(coot::protein_geometry *protein_geom,
                                 std::string dummy_flag_for_CA_mode,  // 20240210-PE this makes this constructor different to the next one
                                                                      // now that we want no_bonds_to_these_atoms to be used in CA mode
                                 const std::set<int> &no_bonds_to_these_atoms_in,
                                 bool do_bonds_to_hydrogens_in=true) :
      no_bonds_to_these_atoms(no_bonds_to_these_atoms_in) {

      init();
      verbose_reporting = false;
      do_bonds_to_hydrogens = do_bonds_to_hydrogens_in;
      do_disulfide_bonds_flag = true;
      b_factor_scale = 1.0;
      have_dictionary = false;
      geom = protein_geom;
      udd_has_ca_handle = -1;
      if (protein_geom)
         have_dictionary = true;
      for_GL_solid_model_rendering = 0;
      if (bonds.size() == 0) {
         for (int i=0; i<13; i++) {  // 13 colors now in bond_colours
            Bond_lines a(i);
            bonds.push_back(a);
         }
      }
   }

   // Constructor M
   //
   // Used by make_colour_by_chain_bonds() - and others in the future?
   // Used by get_bonds_mesh_for_selection_instanced() in coot_molecule_bonds_instanced.cc
   //
   Bond_lines_container(coot::protein_geometry *protein_geom,
                        const std::set<int> &no_bonds_to_these_atoms_in,
                        bool do_bonds_to_hydrogens_in=true) : no_bonds_to_these_atoms(no_bonds_to_these_atoms_in) {
      verbose_reporting = false;
      do_bonds_to_hydrogens = do_bonds_to_hydrogens_in;
      do_disulfide_bonds_flag = true;
      b_factor_scale = 1.0;
      have_dictionary = false;
      geom = protein_geom;
      init();
      geom = protein_geom;
      udd_has_ca_handle = -1;
      if (protein_geom)
         have_dictionary = true;
      for_GL_solid_model_rendering = 0;
      if (bonds.size() == 0) {
         for (int i=0; i<13; i++) {  // 13 colors now in bond_colours
            Bond_lines a(i);
            bonds.push_back(a);
         }
      }
   }

   // Constructor N
   //
   // Phenix Geo
   //
   Bond_lines_container(mmdb::Manager *mol, const coot::phenix_geo::phenix_geometry &pg);

   // Constructor O
   //
   // initial constructor, added to by  addSymmetry_vector_symms from update_symmetry()
   //
   Bond_lines_container() {
      verbose_reporting = false;
      geom = 0;
      do_disulfide_bonds_flag = true;
      do_bonds_to_hydrogens = 1;  // added 20070629
      b_factor_scale = 1.0;
      have_dictionary = 0;
      for_GL_solid_model_rendering = 0;
      udd_has_ca_handle = -1;
      do_sticks_for_waters = true;
      init();
      if (bonds.size() == 0) {
         for (int i=0; i<13; i++) { // 13 colors now in bond_colours
            Bond_lines a(i);
            bonds.push_back(a);
         }
      }
   }
   // used by above:
   void stars_for_unbonded_atoms(mmdb::Manager *mol, int UddHandle);
   void set_udd_unbonded(mmdb::Manager *mol, int UddHandle);
   std::set<int> no_bonds_to_these_atoms; // exclude these atoms because they are part of moving atoms

   coot::rotamer_probability_tables *rotamer_probability_tables_p;
   void add_rotamer_tables(coot::rotamer_probability_tables *tables_p) {
      rotamer_probability_tables_p = tables_p;
   }

   // arguments as above.
   //
   // FYI: there is only one element to symm_trans, the is called from
   // the addSymmetry_vector_symms wrapper
   graphical_bonds_container
   addSymmetry(const atom_selection_container_t &SelAtom, int imol,
                  coot::Cartesian point,
                  float symm_distance,
                  const std::vector<std::pair<symm_trans_t, Cell_Translation> > &symm_trans,
                  short int symmetry_as_ca_flag,
                  short int symmetry_whole_chain_flag);
   graphical_bonds_container
   addSymmetry_with_mmdb(const atom_selection_container_t &SelAtom, int imol,
                            coot::Cartesian point,
                            float symm_distance,
                            const std::vector<std::pair<symm_trans_t, Cell_Translation> > &symm_trans,
                            short int symmetry_as_ca_flag);

   // FYI: there is only one element to symm_trans, the is called from
   // the addSymmetry_vector_symms wrapper
   graphical_bonds_container
      addSymmetry_calphas(const atom_selection_container_t &SelAtom,
                          const coot::Cartesian &point,
                          float symm_distance,
                          const std::vector<std::pair<symm_trans_t, Cell_Translation> > &symm_trans);

   // FYI: there is only one element to symm_trans, the is called from
   // the addSymmetry_vector_symms wrapper
   graphical_bonds_container
   addSymmetry_whole_chain(const atom_selection_container_t &SelAtom, int imol,
                             const coot::Cartesian &point,
                             float symm_distance,
                             const std::vector <std::pair<symm_trans_t, Cell_Translation> > &symm_trans);

   // symmetry with colour-by-symmetry-operator: If
   // symmetry_whole_chain_flag is set, then if any part of a symmetry
   // related molecule is found, then select the whole chain.
   //
   std::vector<std::pair<graphical_bonds_container, std::pair<symm_trans_t, Cell_Translation> > >
   addSymmetry_vector_symms(const atom_selection_container_t &SelAtom,
                            int imol,
                            coot::Cartesian point,
                            float symm_distance,
                            const std::vector<std::pair<symm_trans_t, Cell_Translation> > &symm_trans,
                            short int symmetry_as_ca_flag,
                            short int symmetry_whole_chain_flag,
                            short int draw_hydrogens_flag,
                            bool do_intermolecular_symmetry_bonds);

   std::vector<std::pair<graphical_bonds_container, symm_trans_t> >
   add_NCS(const atom_selection_container_t &SelAtom,
           int imol,
           coot::Cartesian point,
           float symm_distance,
           std::vector<std::pair<coot::coot_mat44, symm_trans_t> > &strict_ncs_mats,
           short int symmetry_as_ca_flag,
           short int symmetry_whole_chain_flag);

   graphical_bonds_container
   add_NCS_molecule(const atom_selection_container_t &SelAtom,
                    int imol,
                    const coot::Cartesian &point,
                    float symm_distance,
                    const std::pair<coot::coot_mat44, symm_trans_t> &strict_ncs_mat,
                    short int symmetry_as_ca_flag,
                    short int symmetry_whole_chain_flag);

   graphical_bonds_container
      add_NCS_molecule_calphas(const atom_selection_container_t &SelAtom,
                               const coot::Cartesian &point,
                               float symm_distance,
                               const std::pair<coot::coot_mat44, symm_trans_t> &strict_ncs_mat);

   graphical_bonds_container
   add_NCS_molecule_whole_chain(const atom_selection_container_t &SelAtom,
                                int imol,
                                const coot::Cartesian &point,
                                float symm_distance,
                                const std::pair<coot::coot_mat44, symm_trans_t> &strict_ncs_mat);

   graphical_bonds_container make_graphical_bonds() const;
   graphical_bonds_container make_graphical_bonds_no_thinning() const;

   graphical_bonds_container make_graphical_bonds(const ramachandrans_container_t &rc,
                                                  bool do_ramachandran_markup,
                                                  bool do_rotamer_markup) const;


   // debugging function
   void check() const;

   void no_symmetry_bonds();


   // void make_graphical_symmetry_bonds() const;
   graphical_bonds_container make_graphical_symmetry_bonds() const;

   void check_graphical_bonds() const;
   void check_static() const;
   void do_disulphide_bonds(atom_selection_container_t, int imodel);
   // which calls either
   void do_disulphide_bonds_by_header(atom_selection_container_t SelAtom, int imodel);
   // or
   void do_disulphide_bonds_by_distance(atom_selection_container_t SelAtom, int imodel);

   // This can get called for intermediate atoms before the restraints have been made
   // (and the FixedDuringRefinement UDD is set), so that loops can flash on
   // then off (when FixedDuringRefinement UDDs *are* set).
   // Oh well, a brief flash is better than permanently on during refinement.
   void do_Ca_bonds(atom_selection_container_t SelAtom,
                    float min_dist, float max_dist, bool draw_missing_loops_flag);
   // make bonds/lies dependent on residue order in molecule - no neighbour search needed. Don't show HETATMs
   coot::my_atom_colour_map_t do_Ca_or_P_bonds_internal(atom_selection_container_t SelAtom,
                                                        const char *backbone_atom_id,
                                                        coot::my_atom_colour_map_t acm,
                                                        float min_dist, float max_dist,
                                                        bool draw_missing_loops_flag,
                                                        int bond_colour_type);
   coot::my_atom_colour_map_t do_Ca_or_P_bonds_internal_old(atom_selection_container_t SelAtom,
                                                        const char *backbone_atom_id,
                                                        coot::my_atom_colour_map_t acm,
                                                        float min_dist, float max_dist,
                                                        int bond_colour_type);
   void do_Ca_loop(int imod, int ires, int nres,
                   mmdb::Chain *chain_p, mmdb::Residue *residue_prev, mmdb::Residue *residue_this,
                   int udd_atom_index_handle,
                   int udd_fixed_during_refinement_handle);

   void do_Ca_plus_ligands_bonds(atom_selection_container_t SelAtom,
                                 int imol,
                                 coot::protein_geometry *pg,
                                 float min_dist, float max_dist,
                                 bool draw_missing_loops_flag,
                                 bool do_bonds_to_hydrogens_in);
   void do_Ca_plus_ligands_bonds(atom_selection_container_t SelAtom,
                                 int imol,
                                 coot::protein_geometry *pg,
                                 float min_dist, float max_dist,
                                 bool draw_missing_loops_flag,
                                 int atom_colour_type,
                                 bool do_bonds_to_hydrogens_in);
   void do_Ca_plus_ligands_and_sidechains_bonds(atom_selection_container_t SelAtom,
                                                int imol,
                                                coot::protein_geometry *pg,
                                                float min_dist_ca, float max_dist_ca,
                                                float min_dist, float max_dist,
                                                bool draw_missing_loops_flag,
                                                bool do_bonds_to_hydrogens_in);
   void do_Ca_plus_ligands_and_sidechains_bonds(atom_selection_container_t SelAtom,
                                                int imol,
                                                coot::protein_geometry *pg,
                                                float min_dist_ca, float max_dist_ca,
                                                float min_dist, float max_dist,
                                                bool draw_missing_loops_flag,
                                                int atom_colour_type,
                                                bool do_bonds_to_hydrogens_in);
   void do_colour_by_chain_bonds(const atom_selection_container_t &asc,
                                 bool use_asc_atom_selection_flag,
                                 int imol,
                                 int draw_hydrogens_flag,
                                 bool draw_missing_loops_flag,
                                 short int change_c_only_flag,
                                 bool do_goodsell_colour_mode,
                                 bool do_ramachandran_markup); // 20221011-PE we want bond by dict *and* rota dodecs!
   void do_colour_by_molecule_bonds(const atom_selection_container_t &asc,
                                    int imol,
                                    int draw_hydrogens_flag);
   void do_normal_bonds_no_water(const atom_selection_container_t &asc,
                                 int imol,
                                 float min_dist, float max_dist);
   void do_colour_sec_struct_bonds(const atom_selection_container_t &asc,
                                   int imol,
                                   float min_dist, float max_dist);
   void do_Ca_plus_ligands_colour_sec_struct_bonds(const atom_selection_container_t &asc,
                                                   int imol,
                                                   coot::protein_geometry *pg,
                                                   float min_dist, float max_dist,
                                                   bool draw_missing_loops_flag,
                                                   bool do_bonds_to_hydrogens_in);

   void do_colour_by_ncs_related_chain_bonds(const atom_selection_container_t &asc,
                                             int imol,
                                             std::vector<std::vector<mmdb::Chain *> > ncs_related_chains,
                                             int draw_mode, // just_atoms_or_bonds_and_atoms mode
                                             bool change_c_only_flag, bool goodsell_mode);

   atom_selection_container_t
      ContactSel(mmdb::PPAtom trans_sel, mmdb::Contact *contact, int ncontacts) const;
   void do_symmetry_Ca_bonds(atom_selection_container_t SelAtom,
                             symm_trans_t symm_trans);

   void set_verbose_reporting(short int i) { verbose_reporting = i;}

   std::vector<coot::torus_description_t> rings; // for OpenGL rendering of aromaticity bonding.
   bool have_rings() const {
      return ! rings.empty();
   }

   class symmetry_atom_bond {
   public:
      // bond between at_1 and symmetry-related copy of at_2
      mmdb::Atom *at_1;
      mmdb::Atom *at_2;
      symm_trans_t st;
      Cell_Translation ct;
      symmetry_atom_bond(mmdb::Atom *at_1_in, mmdb::Atom *at_2_in,
                         const symm_trans_t &st_in,
                         const Cell_Translation &ct_in) : st(st_in), ct(ct_in) {
         at_1 = at_1_in;
         at_2 = at_2_in;
      }
      // modify m
      int GetTMatrix(mmdb::Manager *mol, mmdb::mat44 *m) const {
         return mol->GetTMatrix(*m, st.isym(), st.x(), st.y(), st.z());
      }
   };

   std::vector<symmetry_atom_bond> find_intermolecular_symmetry(const atom_selection_container_t &SelAtom) const;

   graphical_bonds_container
   intermolecular_symmetry_graphical_bonds(mmdb::Manager *mol,
                                           const std::vector <symmetry_atom_bond> &sabv,
                                           const std::pair<symm_trans_t, Cell_Translation> &symm_trans);

   void set_use_deuteranomaly_mode() { use_deuteranomaly_mode = true; }

   // 20251014-PE lightweight update the bonds without recalculating the bonding. make_graphical_bonds() will
   // still be needed for now.
   // Note that some cylinder (in kekulized mode, say) will not be terminating at atom positions.
   // Using this update() will make them do so.
   //
   // Note to self: I need to update the atoms too.
   //
   void update(mmdb::Atom **atom_selection, int n_atoms);

};

#endif /* BOND_LINES_H */
