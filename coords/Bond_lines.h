// -*-c++-*-
/* coords/Bond_lines.h
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2016 by Medical Research Council
 * Author: Paul Emsley
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
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
#include "Cartesian.h"
#include "phenix-geo.hh"

#include "coot-utils/coot-rama.hh" // for ramachandran scoring of intermediate atoms
#include "ramachandran-container.hh"
#include "rotamer-container.hh"
#include "ligand/rotamer.hh"
#include "coot-utils/coot-coord-utils.hh" // is this needed?

#include "geometry/bonded-quad.hh"

namespace coot {

   static std::string b_factor_bonds_scale_handle_name;

   enum bond_colour_t { COLOUR_BY_CHAIN=0,
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
      int index_for_b_factor(float wheel_colour) {
         return int(30.0*wheel_colour);
      }
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

   // For OpenGL solid model, it is much better to draw a torus than a
   // set of short sticks (particularly for the ring representing
   // aromaticity).  So now (20100831 Bond_lines_container contains a
   // number of torus descriptions).
   //
   class torus_description_t {
   public:
      double inner_radius;
      double outer_radius;
      int n_sides;
      int n_rings;
      clipper::Coord_orth centre;
      clipper::Coord_orth normal;
      torus_description_t(const clipper::Coord_orth &pt,
			  const clipper::Coord_orth &normal_in,
			  double ir1, double ir2, int n1, int n2) {
	 inner_radius = ir1;
	 outer_radius = ir2;
	 n_sides = n1;
	 n_rings = n2;
	 centre = pt;
	 normal = normal_in;
      }
   };
}

class graphics_line_t {
public:
   enum cylinder_class_t { UNK, SINGLE, DOUBLE, TRIPLE }; // so that double bonds can be drawn thinner
                                                          // than single bonds (likewise triple)
   cylinder_class_t cylinder_class;
   coot::CartesianPair positions;
   bool has_begin_cap;
   bool has_end_cap;
   // int residue_index;
   // restore this when finished
   // mmdb::Residue *residue_p; // the residue for the bond (maybe there should be 2 residues_ps? because
                             // sometimes there will be 2 residues for the same graphics_line_t.
                             // Hmm.
   int model_number; // -1 is unset
   int atom_index_1;
   int atom_index_2;
#if 0
   // default single bond constructor
   graphics_line_t(const coot::CartesianPair &p, bool b, bool e) {
      positions = p;
      has_begin_cap = b;
      has_end_cap = e;
      cylinder_class = SINGLE;
      // residue_index = -1; // unset
      //residue_p = 0;
      atom_index_1 = -1;
      atom_index_2 = -1;
      model_number = -1;
   }
#endif
   // we want atom indices now, not just the residue
   // graphics_line_t(const coot::CartesianPair &p, cylinder_class_t cc, bool b, bool e, mmdb::Residue *residue_p_in);

   graphics_line_t(const coot::CartesianPair &p, cylinder_class_t cc, bool b, bool e,
		   int model_no_in,
		   int atom_index_1_in, int atom_index_2_in) {
      positions = p;
      has_begin_cap = b;
      has_end_cap = e;
      cylinder_class = cc;
      atom_index_1 = atom_index_1_in;
      atom_index_2 = atom_index_2_in;
      model_number = model_no_in;
   }
   graphics_line_t() { }
};
 
// A poor man's vector. Contains the set of lines for each element (well,
// basically colour) type.
// 
template<class T> class graphical_bonds_lines_list {

 public:
   int num_lines;
   T *pair_list;
   bool thin_lines_flag;

   graphical_bonds_lines_list() {
      pair_list = NULL;
      thin_lines_flag = 0;
   }
};


class graphical_bonds_atom_info_t {
public:
   bool is_hydrogen_atom;
   bool is_water; // don't not display this in sticks-only mode - or the water
                  // will disappear (needs rebonding otherwise - not just a display flag)
   float radius_scale; // Waters (and perhaps metals) should have big radii, so that
                       // they are easier to see.
   coot::Cartesian position;
   mmdb::Atom *atom_p; // this should be a shared pointer I think.
                       // we don't want to be looking at this pointer
                       // if some other part of the code has deleted the atom.
   int model_number; // -1 is unset
   int atom_index;
   graphical_bonds_atom_info_t(const coot::Cartesian &pos, int atom_index_in, bool is_hydrogen_atom_in) {
      model_number = -1;
      position = pos;
      is_hydrogen_atom = is_hydrogen_atom_in;
      is_water = false;
      atom_index = atom_index_in;
      atom_p = 0;
      radius_scale = 1.0;
   }
   graphical_bonds_atom_info_t() {
      model_number = -1;
      is_hydrogen_atom = false;
      atom_index = -1; // unset
      radius_scale = 1.0;
      atom_p = 0;
      is_water = false;
   }
   // this is a bit of a weird construction
   bool radius_for_atom_should_be_big(mmdb::Atom *atom_p) const {

      // 20190822-PE: you might like to add other tests here.
      // 20200608-PE: I did!
      mmdb::Residue *r = atom_p->GetResidue();
      if (r) {
         std::string res_name = r->GetResName();
         if (res_name == "HOH")
            return true;
         if (res_name == "CA")
            return true;
         if (res_name == "MG")
            return true;
         if (res_name == "IOD")
            return true;
         if (res_name == "CL")
            return true;
         if (res_name == "NA")
            return true;
         if (res_name == "K")
            return true;
      }
      return false;
   }
};

template<class T> class graphical_bonds_points_list {

public:
   unsigned int num_points;
   unsigned int current_count;

   // use a is-H-atom-flag for first
   T *points;
   
   graphical_bonds_points_list() {
      current_count = 0;
      num_points = 0;
      points = NULL;
   }

   graphical_bonds_points_list(unsigned int size) {
      current_count = 0;
      num_points = size;
      points = new T[size];
   }

   void add_point(const T &pt) {
      points[current_count] = pt;
      current_count++;
   } 
};

class graphical_bonds_cis_peptide_markup {

public:

   bool is_pre_pro_cis_peptide;
   bool is_twisted; // twisted trans
   int model_number; // -1 is unset
   coot::Cartesian pt_ca_1; 
   coot::Cartesian pt_c_1;
   coot::Cartesian pt_n_2;
   coot::Cartesian pt_ca_2;
   coot::atom_index_quad atom_index_quad;
   graphical_bonds_cis_peptide_markup(const coot::Cartesian &pt_ca_1_in,
				      const coot::Cartesian &pt_c_1_in,
				      const coot::Cartesian &pt_n_2_in,
				      const coot::Cartesian &pt_ca_2_in,
				      bool is_pre_pro_cis_peptide_in,
				      bool is_twisted_in,
				      int model_number_in) {
      pt_ca_1 = pt_ca_1_in;
      pt_c_1  = pt_c_1_in;
      pt_n_2  = pt_n_2_in;
      pt_ca_2 = pt_ca_2_in;
      is_pre_pro_cis_peptide = is_pre_pro_cis_peptide_in;
      is_twisted = is_twisted_in;
      model_number = model_number_in;
   }

   void add_atom_index_quad(const coot::atom_index_quad &iq) {
      atom_index_quad = iq;
   }

   graphical_bonds_cis_peptide_markup() {
      is_pre_pro_cis_peptide = false;
   } 
};

// Uses graphical_bonds_lines_list
// 
// a graphical_bonds_container is used by draw_molecule() and are
// created from a Bond_lines_container (which uses a vector).
// 
class graphical_bonds_container { 

 public:
   
   int num_colours;
   graphical_bonds_lines_list<graphics_line_t> *bonds_; 

   int symmetry_has_been_created;
   graphical_bonds_lines_list<graphics_line_t> *symmetry_bonds_;

   // if the distance between CAs in a missing loop is longer than is possible given
   // the residue number difference, then we want to mark up that line with
   // dots along the line joining the residues.  This should work similarly with P-P
   // for nucleic acid chains - but I won't change the function name.
   coot::Cartesian *bad_CA_CA_dist_spots_ptr;
   coot::Cartesian *zero_occ_spots_ptr;
   coot::Cartesian *deuterium_spots_ptr;
   std::pair<coot::Cartesian, float> *ramachandran_goodness_spots_ptr;
   int n_zero_occ_spots;
   int n_bad_CA_CA_dist_spots;
   int n_deuterium_spots;
   int n_ramachandran_goodness_spots;

   // first is is-H-atom-flag
   graphical_bonds_atom_info_t *atom_centres_;
   int n_atom_centres_;
   int *atom_centres_colour_;
   std::vector<coot::torus_description_t> rings;
   int n_consolidated_atom_centres;
   graphical_bonds_points_list<graphical_bonds_atom_info_t> *consolidated_atom_centres;
   int n_cis_peptide_markups;
   graphical_bonds_cis_peptide_markup *cis_peptide_markups;
   int n_rotamer_markups;
   rotamer_markup_container_t *rotamer_markups;
   
   graphical_bonds_container() { 
      num_colours = 0; 
      bonds_ = NULL;
      symmetry_has_been_created = 0; 
      symmetry_bonds_ = NULL; 
      zero_occ_spots_ptr = NULL;
      bad_CA_CA_dist_spots_ptr = NULL;
      n_bad_CA_CA_dist_spots = 0;
      n_zero_occ_spots = 0;
      deuterium_spots_ptr = NULL;
      n_deuterium_spots = 0;
      atom_centres_colour_ = NULL;
      atom_centres_ = NULL; 
      n_atom_centres_ = 0;
      n_ramachandran_goodness_spots = 0;
      ramachandran_goodness_spots_ptr = NULL;
      consolidated_atom_centres = NULL;
      n_consolidated_atom_centres = 0;
      n_cis_peptide_markups = 0;
      cis_peptide_markups = NULL;
      n_rotamer_markups = 0;
      rotamer_markups = NULL;
   }

   void clear_up() {

      if (bonds_)
	 for (int icol=0; icol<num_colours; icol++)
	    delete [] bonds_[icol].pair_list;
      if (symmetry_bonds_)
	 for (int icol=0; icol<num_colours; icol++)
	    delete [] symmetry_bonds_[icol].pair_list;

      delete [] bonds_;  // null testing part of delete
      delete [] symmetry_bonds_; 
      delete [] atom_centres_;
      delete [] atom_centres_colour_;
      bonds_ = NULL;
      symmetry_bonds_ = NULL;
      atom_centres_ = NULL;
      atom_centres_colour_ = NULL;
      if (n_zero_occ_spots) 
	 delete [] zero_occ_spots_ptr;
      if (n_deuterium_spots)
	 delete [] deuterium_spots_ptr;
      if (n_ramachandran_goodness_spots)
	 delete [] ramachandran_goodness_spots_ptr;
      n_zero_occ_spots = 0;
      n_deuterium_spots = 0;
      n_ramachandran_goodness_spots = 0;
      n_atom_centres_ = 0;
      if (consolidated_atom_centres) {
	 for (int i=0; i<n_consolidated_atom_centres; i++)
	    delete [] consolidated_atom_centres[i].points;
	 delete [] consolidated_atom_centres;
	 consolidated_atom_centres = NULL;
      }
      delete [] cis_peptide_markups;
      cis_peptide_markups = NULL;
      n_rotamer_markups = 0;
      delete [] rotamer_markups;
      rotamer_markups = NULL;
   }

   graphical_bonds_container(const std::vector<graphics_line_t> &a) { 

      std::cout << "constructing a graphical_bonds_container from a vector " 
		<< "of size " << a.size() << std::endl;

      num_colours = 1;
      
      bonds_ = new graphical_bonds_lines_list<graphics_line_t>[1]; // only 1 graphical_bonds_lines_list needed
      bonds_[0].pair_list = new graphics_line_t[(a.size())];
      bonds_[0].num_lines = a.size();

      // copy over
      for(int i=0; i<bonds_[0].num_lines; i++)
	 bonds_[0].pair_list[i] = a[i];

      symmetry_bonds_ = NULL; 
      symmetry_has_been_created = 0; 
      zero_occ_spots_ptr = NULL;
      n_zero_occ_spots = 0;
      deuterium_spots_ptr = NULL;
      n_deuterium_spots = 0;
      atom_centres_colour_ = NULL;
      atom_centres_ = NULL; 
      n_atom_centres_ = 0;
      n_ramachandran_goodness_spots = 0;
      ramachandran_goodness_spots_ptr = NULL;
      consolidated_atom_centres = NULL;
      n_consolidated_atom_centres = 0;
      n_cis_peptide_markups = 0;
      cis_peptide_markups = NULL;
   }

   // This function used by skeletonization (and nothing else)
   void add_colour(const  std::vector<graphics_line_t> &a);
      
   void add_zero_occ_spots(const std::vector<coot::Cartesian> &spots);
   void add_bad_CA_CA_dist_spots(const std::vector<coot::Cartesian> &spots);
   void add_deuterium_spots(const std::vector<coot::Cartesian> &spots);
   void add_ramachandran_goodness_spots(const std::vector<std::pair<coot::Cartesian,
					coot::util::phi_psi_t> > &spots,
					const ramachandrans_container_t &rc);
   void add_rotamer_goodness_markup(const std::vector<rotamer_markup_container_t> &ric);

   void add_atom_centres(const std::vector<graphical_bonds_atom_info_t> &centres,
			 const std::vector<int> &colours);
   bool have_rings() const { return rings.size(); }
   bool empty() const { return (bonds_ == NULL); }
   void add_cis_peptide_markup(const std::vector<coot::util::cis_peptide_quad_info_t> &cis_peptide_quads);
};



// Bond_lines is a container class, containing a colour index
// and a vector of pairs of coot::Cartesians that should be drawn
// _in_ that colour.
//
class Bond_lines { 
   int colour;
   std::vector<graphics_line_t> points;

 public:
   Bond_lines(const graphics_line_t &pts);
   Bond_lines(); 
   Bond_lines(int col);

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
};

// 
enum symm_keys {NO_SYMMETRY_BONDS}; 

class Bond_lines_container { 

   enum { NO_BOND,
	  BONDED_WITH_STANDARD_ATOM_BOND,
	  BONDED_WITH_BOND_TO_HYDROGEN, 
	  BONDED_WITH_HETATM_BOND /* by dictionary */ };
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

   void construct_from_model_links(mmdb::Model *model, int udd_atom_index_handle, int atom_colour_type);
   // which wraps...
   void add_link_bond(mmdb::Model *model_p, int udd_atom_index_handle, int atom_colour_type, mmdb::Link *link);
   void add_link_bond(mmdb::Model *model_p, int udd_atom_index_handle, int atom_colour_type, mmdb::LinkR *linkr);

   template<class T> void add_link_bond_templ(mmdb::Model *model_p, int udd_atom_index_handle, int atom_colour_type, T *link);

   // now wit optional arg.  If atom_colour_type is set, then use/fill
   // it to get colour indices from chainids.
   void handle_MET_or_MSE_case (mmdb::PAtom mse_atom, int udd_handle,
				int udd_handle_for_atom_index, int atom_colour_type,
				coot::my_atom_colour_map_t *atom_colour_map = 0);
   void handle_long_bonded_atom(mmdb::PAtom atom,
				int udd_handle_bond,
				int udd_atom_index_handle,
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
   void add_atom_centres(const atom_selection_container_t &SelAtom, int atom_colour_type,
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
		       int model_number,
		       int atom_index_1,
		       int atom_index_2,
		       int atom_colour_type,
                       coot::my_atom_colour_map_t *atom_colour_map_p);

   // double and delocalized bonds (default (no optional arg) is double).
   // We pass udd_atom_index_handle because we need the atom index (not residue atom index) for
   // using no_bonds_to_these_atoms
   void add_double_bond(int imol, int imodel, int iat_1, int iat_2, mmdb::PPAtom atoms, int n_atoms,
                        int atom_colour_type, coot::my_atom_colour_map_t *atom_colour_map_p,
			int udd_atom_index_handle,
			const std::vector<coot::dict_bond_restraint_t> &bond_restraints,
			bool is_deloc=false);
   // used by above, can throw an exception
   clipper::Coord_orth get_neighb_normal(int imol, int iat_1, int iat_2, mmdb::PPAtom atoms, int n_atoms, 
	 				 bool also_2nd_order_neighbs=0) const;
   void add_triple_bond(int imol, int imodel, int iat_1, int iat_2, mmdb::PPAtom atoms, int n_atoms,
                        int atom_colour_type,
                        coot::my_atom_colour_map_t *atom_colour_map_p,
			int udd_atom_index_handle,
			const std::vector<coot::dict_bond_restraint_t> &bond_restraints);


   bool have_dictionary;
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
		bool add_begin_end_cap = false,
		bool add_end_end_cap = false);
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
   int atom_colour(mmdb::Atom *at, int bond_colour_type, coot::my_atom_colour_map_t *atom_colour_map = 0);
   void bonds_size_colour_check(int icol) {
      int bonds_size = bonds.size();
      if (icol >= bonds_size) 
	 bonds.resize(icol+1);
   }

   // return the UDD handle
   int set_rainbow_colours(mmdb::Manager *mol);
   void do_colour_by_chain_bonds_carbons_only(const atom_selection_container_t &asc,
					      int imol,
                                              bool draw_missing_loops_flag,
					      int atom_colour_type, // C-only or goodsell
					      int draw_hydrogens_flag);
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
							int uddHnd);
   void do_colour_by_chain_bonds_internals_goodsell_mode(int imol, int imodel,
							int chain_idx,
							mmdb::Atom *at1, mmdb::Atom *at2,
							int iat_1, int iat_2,
							std::vector<std::pair<bool, mmdb::Residue *> > *het_residues_p,
							const std::string &element_1,
							const std::string &element_2,
							const coot::Cartesian &atom_1,
							const coot::Cartesian &atom_2,
							int uddHnd);

   void do_colour_by_dictionary_and_by_chain_bonds(const atom_selection_container_t &asc,
                                                   int imol,
                                                   int draw_hydrogens_flag,
                                                   bool draw_missing_loops_flag,
                                                   short int change_c_only_flag,
                                                   bool do_goodsell_colour_mode);

   void add_residue_monomer_bonds(const std::map<std::string, std::vector<mmdb::Residue *> > &residue_monomer_map,
                                  int imol, int model_number,
                                  int atom_colour_type,
                                  coot::my_atom_colour_map_t *atom_colour_map_p,
                                  int udd_atom_index_handle,
                                  int udd_bond_handle,
                                  int draw_hydrogens_flag,
                                  bool do_goodsell_colour_mode);

   void do_colour_by_dictionary_and_by_chain_bonds_carbons_only(const atom_selection_container_t &asc,
                                                                int imol,
                                                                int draw_hydrogens_flag,
                                                                bool draw_missing_loops_flag,
                                                                bool do_goodsell_colour_mode);
   // and the bonds between the above monomers
   void add_polymer_bonds(const atom_selection_container_t &asc,
                          int atom_colour_type,
                          coot::my_atom_colour_map_t *atom_colour_map_p,
                          int draw_hydrogens_flag,
                          bool draw_missing_loops_flag,
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
   void draw_phenyl_ring(const std::vector<mmdb::Atom *> &ring_atoms, int imodel,
                         int atom_colour_type, coot::my_atom_colour_map_t *atom_colour_map_p,
                         int udd_atom_index_handle);
   // this calls the above function
   void draw_phenyl_ring_outer(mmdb::Residue *residue_p, int model_number,
                               int atom_colour_type, coot::my_atom_colour_map_t *atom_colour_map_p,
                               int udd_atom_index_handle);
   void draw_trp_rings(const std::vector<mmdb::Atom *> &ring_atoms, int imodel,
                       int atom_colour_type, coot::my_atom_colour_map_t *atom_colour_map_p,
                       int udd_atom_index_handle);
   void draw_trp_ring_outer(mmdb::Residue *residue_p, int model_number,
                            int atom_colour_type, coot::my_atom_colour_map_t *atom_colour_map_p,
                            int udd_atom_index_handle);
   void draw_GA_rings_outer(mmdb::Residue *residue_p, int model_number,
                            int atom_colour_type, coot::my_atom_colour_map_t *atom_colour_map_p,
                            int udd_atom_index_handle);
   void draw_CUT_ring_outer(mmdb::Residue *residue_p, int model_number,
                            int atom_colour_type, coot::my_atom_colour_map_t *atom_colour_map_p,
                            int udd_atom_index_handle);
   void draw_het_group_rings(mmdb::Residue *residue_p,
                             const std::vector<bonded_quad_atom_names> &bonded_quad_atom_names,
                             int model_number, int atom_colour_type,
                             coot::my_atom_colour_map_t *atom_colour_map,
                             int udd_atom_index_handle);
   void draw_bonded_quad_atoms_rings(const std::vector<bonded_quad_atoms> &ring_atoms,
                             int imodel, int atom_colour_type,
                             coot::my_atom_colour_map_t *atom_colour_map_p,
                             int udd_atom_index_handle);
   
   std::vector<std::pair<std::string, std::string> >
   get_aromatic_bonds(const coot::dictionary_residue_restraints_t &restraints) const;

   void try_set_b_factor_scale(mmdb::Manager *mol);
   graphical_bonds_container make_graphical_bonds_with_thinning_flag(bool thinning_flag) const;
   void add_bonds_het_residues(const std::vector<std::pair<bool, mmdb::Residue *> > &het_residues,
                               int imol, int atom_colour_t, short int have_udd_atoms,
                               int udd_found_bond_handle, int udd_atom_index_handle);
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

   // we can put other things here
   void init() {
      rotamer_probability_tables_p = NULL;
      do_sticks_for_waters = false;
      n_atoms_in_atom_selection = 0;
   }
   

public:
   enum bond_representation_type { COLOUR_BY_OCCUPANCY, COLOUR_BY_B_FACTOR, COLOUR_BY_USER_DEFINED_COLOURS}; 

   // getting caught out with Bond_lines_container dependencies?
   // We need:  mmdb-extras.h which needs mmdb-manager.h and <string>
   // Bond_lines_container(const atom_selection_container_t &asc);
   Bond_lines_container(const atom_selection_container_t &asc,
			int imol,
			int include_disulphides=0,
			int include_hydrogens=1,
			bool do_rama_markup=false,
			bool do_rota_markup=false,
			coot::rotamer_probability_tables *tables_p=0);

   // the constructor for bond by dictionary - should use this most of the time.
   // geom_in can be null if you don't have it.
   //
   // if model_number is 0, display all models. If it is not 0 then
   // display only the given model_number (if possible, of course).
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

   Bond_lines_container(atom_selection_container_t, int imol, float max_dist);

   Bond_lines_container(atom_selection_container_t asc,
			int imol,
			float min_dist, float max_dist);

   // The constructor for ball and stick, this constructor implies that
   // for_GL_solid_model_rendering is set.
   // 
   // geom_in can be null.
   Bond_lines_container (atom_selection_container_t asc,
			 int imol,
			 const coot::protein_geometry *geom);

   Bond_lines_container(atom_selection_container_t SelAtom, coot::Cartesian point,
			float symm_distance,
			std::vector<symm_trans_t> symm_trans); // const & FIXME

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

   // same as above, except this does symmetry contacts too.
   // 
   Bond_lines_container(const atom_selection_container_t &SelAtom,
			mmdb::PPAtom residue_atoms,
			int n_residue_atoms,
			float min_dist,
			float max_dist, 
			bool draw_env_distances_to_hydrogens_flag,
			short int do_symmetry);

   // This is the one for occupancy and B-factor representation
   // 
   Bond_lines_container (const atom_selection_container_t &SelAtom,
			 int imol, bond_representation_type by_occ);

   Bond_lines_container(int col);
   Bond_lines_container(symm_keys key);


   // Used by various CA-mode bonds
   //
   Bond_lines_container(coot::protein_geometry *protein_geom, bool do_bonds_to_hydrogens_in=true) {
      do_bonds_to_hydrogens = do_bonds_to_hydrogens_in;
      b_factor_scale = 1.0;
      have_dictionary = false;
      geom = protein_geom;
      init();
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

   // Used by make_colour_by_chain_bonds() - and others in the future?
   //
   Bond_lines_container(coot::protein_geometry *protein_geom,
			const std::set<int> &no_bonds_to_these_atoms_in,
			bool do_bonds_to_hydrogens_in=true) : no_bonds_to_these_atoms(no_bonds_to_these_atoms_in) {
      do_bonds_to_hydrogens = do_bonds_to_hydrogens_in;
      b_factor_scale = 1.0;
      have_dictionary = false;
      geom = protein_geom;
      init();
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

   // Phenix Geo
   // 
   Bond_lines_container(mmdb::Manager *mol, const coot::phenix_geo_bonds &gb);

   // initial constructor, added to by  addSymmetry_vector_symms from update_symmetry()
   // 
   Bond_lines_container() {
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
				 bool do_goodsell_colour_mode);
   void do_colour_by_molecule_bonds(const atom_selection_container_t &asc,
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

   atom_selection_container_t
      ContactSel(mmdb::PPAtom trans_sel, mmdb::Contact *contact, int ncontacts) const;
   void do_symmetry_Ca_bonds(atom_selection_container_t SelAtom,
			     symm_trans_t symm_trans);

   void set_verbose_reporting(short int i) { verbose_reporting = i;}

   std::vector<coot::torus_description_t> rings; // for OpenGL rendering of aromaticity bonding.
   bool have_rings() const {
      return rings.size();
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
			 const Cell_Translation &ct_in) {
	 at_1 = at_1_in;
	 at_2 = at_2_in;
	 ct = ct_in;
	 st = st_in;
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


};

#endif /* BOND_LINES_H */
