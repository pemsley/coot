// -*-c++-*-
/* coords/Bond_lines.h
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
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

#include "protein-geometry.hh"

namespace coot { 

   static std::string b_factor_bonds_scale_handle_name;
   
   enum bond_colour_t { COLOUR_BY_CHAIN=0, COLOUR_BY_ATOM_TYPE=1, 
		       COLOUR_BY_SEC_STRUCT=2, DISULFIDE_COLOUR=3,
                       COLOUR_BY_MOLECULE=4, COLOUR_BY_RAINBOW=5, 
		       COLOUR_BY_OCCUPANCY=6, COLOUR_BY_B_FACTOR=7};

  class my_atom_colour_map_t { 

    public:
    std::vector<std::string> atom_colour_map;
    int index_for_chain(const std::string &chain) { 
      int isize = atom_colour_map.size();
      for(int i=0; i<isize; i++) { 
	if (atom_colour_map[i] == chain) {
	  return i;
	}
      }
      atom_colour_map.push_back(chain);
      return isize;
    }
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
      std::vector<PCAtom> hydrogen_atoms_;
      std::vector<PCAtom> non_hydrogen_atoms_;
   public:
      PPCAtom     Hydrogen_atoms() const;
      PPCAtom non_Hydrogen_atoms() const;
      int n_H() const { return hydrogen_atoms_.size(); }
      int n_non_H() const { return non_hydrogen_atoms_.size(); }
      void add_atom(CAtom *atom) {
	 std::string element = atom->element;
	 if (element == " H" || element == " D") {
	    hydrogen_atoms_.push_back(atom);
	 } else {
	    non_hydrogen_atoms_.push_back(atom);
	 }
      }
   };
   
}
 
// A poor man's vector.  For use when we can't use vectors
// 
class Lines_list { 

 public:   
   // contain a number of elements
   int num_lines; 
   coot::CartesianPair *pair_list;
   bool thin_lines_flag;
   
   Lines_list() { 
      pair_list = NULL;
      thin_lines_flag = 0;
   }
   
};

// Uses Lines_list
// 
// a graphical_bonds_container is used by draw_molecule() and are
// created from a Bond_lines_container (which uses a vector).
// 
class graphical_bonds_container { 

 public:
   
   int num_colours;
   Lines_list *bonds_; 

   int symmetry_has_been_created;
   Lines_list *symmetry_bonds_;

   coot::Cartesian *zero_occ_spot;
   int n_zero_occ_spot;

   coot::Cartesian *atom_centres_;
   int n_atom_centres_;
   int *atom_centres_colour_;

   graphical_bonds_container() { 
      num_colours = 0; 
      bonds_ = NULL;
      symmetry_has_been_created = 0; 
      symmetry_bonds_ = NULL; 
      zero_occ_spot = NULL;
      n_zero_occ_spot = 0;
      atom_centres_colour_ = NULL;
      atom_centres_ = NULL; 
      n_atom_centres_ = 0;
   }

   void clear_up() { 

      if (bonds_)
	 for (int icol=0; icol<num_colours; icol++)
	    delete [] bonds_[icol].pair_list;
      if (symmetry_bonds_)
	 for (int icol=0; icol<num_colours; icol++)
	    delete [] symmetry_bonds_[icol].pair_list;

//       if (bonds_)
// 	 std::cout << " clearing bonds " << bonds_ << std::endl;
//       if (symmetry_bonds_)
// 	 std::cout << " clearing symmetry bonds " << symmetry_bonds_ << std::endl;
      
      delete [] bonds_;  // null testing part of delete
      delete [] symmetry_bonds_; 
      delete [] atom_centres_;
      delete [] atom_centres_colour_;
      bonds_ = NULL;
      symmetry_bonds_ = NULL;
      atom_centres_ = NULL;
      atom_centres_colour_ = NULL;
      if (n_zero_occ_spot) 
	 delete [] zero_occ_spot;
      n_zero_occ_spot = 0;
      n_atom_centres_ = 0;
   }

   graphical_bonds_container(std::vector<coot::CartesianPair> &a) { 

      std::cout << "constructing a graphical_bonds_container from a vector " 
		<< "of size " << a.size() << std::endl;

      num_colours = 1;
      
      bonds_ = new Lines_list[1]; // only 1 Lines_list needed
      bonds_[0].pair_list = new coot::CartesianPair[(a.size())];

      bonds_[0].num_lines = a.size();

      // copy over
      for(int i=0; i<bonds_[0].num_lines; i++) { 
	bonds_[0].pair_list[i] = a[i]; 
      }
      symmetry_bonds_ = NULL; 
      symmetry_has_been_created = 0; 
      zero_occ_spot = NULL;
      n_zero_occ_spot = 0;
      atom_centres_colour_ = NULL;
      atom_centres_ = NULL; 
      n_atom_centres_ = 0;
   }
      
   void add_colour( std::vector<coot::CartesianPair> &a ) {
 /*       cout << "filling a graphical_bonds_container from a vector "  */
/* 	   << "of size " << a.size() << endl; */

     Lines_list *new_bonds_ = new Lines_list[num_colours+1];
     if ( bonds_ != NULL ) {
       for (int i = 0; i < num_colours; i++ ) new_bonds_[i] = bonds_[i];
       delete[] bonds_;
     }
     bonds_ = new_bonds_;
     bonds_[num_colours].pair_list = new coot::CartesianPair[(a.size())];
     bonds_[num_colours].num_lines = a.size();

     // copy over
     for(unsigned int i=0; i<a.size(); i++) { 
       bonds_[num_colours].pair_list[i] = a[i]; 
     }
     num_colours++;

     symmetry_bonds_ = NULL;
     symmetry_has_been_created = 0; 
   }

   void add_zero_occ_spots(const std::vector<coot::Cartesian> &spots);
   void add_atom_centre(const coot::Cartesian &pos, int colour_index);
   void add_atom_centres(const std::vector<coot::Cartesian> &centres,
			 const std::vector<int> &colours);
};


// Bond_lines is a container class, containing a colour index
// and a vector of pairs of coot::Cartesians that should be drawn
// _in_ that colour.
//
class Bond_lines { 
   int colour;
   std::vector<coot::CartesianPair> points;

 public:
   Bond_lines(coot::CartesianPair pts);
   Bond_lines(); 
   Bond_lines(int col);

   void add_bond(coot::CartesianPair);
   int size() const; 

   // return the coordinates of the start and finish points of the i'th bond.
   //
   coot::Cartesian GetStart(int i) const;
   coot::Cartesian GetFinish(int i) const;
   coot::CartesianPair operator[](int i) const;
};

// 
enum symm_keys {NO_SYMMETRY_BONDS}; 

class Bond_lines_container { 

   short int verbose_reporting;
   short int do_disulfide_bonds_flag;
   short int do_bonds_to_hydrogens;
   int udd_has_ca_handle;
   float b_factor_scale; 

   // we rely on SelAtom.atom_selection being properly constucted to
   // contain all atoms
   //
   void construct_from_asc(const atom_selection_container_t &SelAtom, 
			   float min_dist, float max_dist, 
			   int atom_colour_type, 
			   short int is_from_symmetry_flag);

   void construct_from_atom_selection(const atom_selection_container_t &asc,
				      const PPCAtom atom_selection_1,
				      int n_selected_atoms_1,
				      const PPCAtom atom_selection_2,
				      int n_selected_atoms_2,
				      float min_dist, float max_dist,
				      int atom_colour_type,
				      short int are_different_atom_selections,
				      short int have_udd_atoms,
				      int udd_handle);
   void construct_from_model_links(CModel *model, int atom_colour_type);

   void handle_MET_or_MSE_case (PCAtom mse_atom, int udd_handle, int atom_colour_type);
   void handle_long_bonded_atom(PCAtom     atom, int udd_handle, int atom_colour_type);

   // void check_atom_limits(atom_selection_container_t SelAtoms) const;
   
   void write(std::string) const;

   PPCAtom trans_sel(atom_selection_container_t AtomSel, 
		     const std::pair<symm_trans_t, Cell_Translation> &symm_trans) const;

   void add_zero_occ_spots(const atom_selection_container_t &SelAtom);
   void add_atom_centres(const atom_selection_container_t &SelAtom, int atom_colour_type);
   int add_ligand_bonds(const atom_selection_container_t &SelAtom, 
			PPCAtom ligand_atoms_selection,
			int n_ligand_atoms);

   bool draw_these_residue_contacts(CResidue *this_residue, CResidue *env_residue,
				    coot::protein_geometry *protein_geom);

                           
 protected:
   std::vector<Bond_lines> bonds; 
   std::vector<coot::Cartesian>  zero_occ_spot;
   std::vector<coot::Cartesian>  atom_centres;
   std::vector<int>        atom_centres_colour;
   void addBond(int colour, const coot::Cartesian &first, const coot::Cartesian &second);
   void addBondtoHydrogen(const coot::Cartesian &first, const coot::Cartesian &second);
   void add_dashed_bond(int col,
			const coot::Cartesian &start,
			const coot::Cartesian &end,
			bool is_half_bond_flag);
   void addAtom(int colour, const coot::Cartesian &pos);
   int atom_colour(CAtom *at, int bond_colour_type);
   void bonds_size_colour_check(int icol) {
      int bonds_size = bonds.size();
      if (icol >= bonds_size) 
	 bonds.resize(icol+1);
   }

   // return the UDD handle
   int set_rainbow_colours(int selHnd_ca, CMMDBManager *mol);
   void do_colour_by_chain_bonds_change_only(const atom_selection_container_t &asc,
					     int draw_hydrogens_flag);

   void try_set_b_factor_scale(CMMDBManager *mol);

public:
   enum bond_representation_type { COLOUR_BY_OCCUPANCY, COLOUR_BY_B_FACTOR}; 

   // getting caught out with Bond_lines_container dependencies?
   // We need:  mmdb-extras.h which needs mmdb-manager.h and <string>
   // Bond_lines_container(const atom_selection_container_t &asc);

   Bond_lines_container(const atom_selection_container_t &asc, 
			int include_disulphides=0,
			int include_hydrogens=1);
   Bond_lines_container(atom_selection_container_t, float max_dist);
   Bond_lines_container(atom_selection_container_t, 
			float min_dist, float max_dist);
   Bond_lines_container(atom_selection_container_t SelAtom, coot::Cartesian point,
			float symm_distance,
			std::vector<symm_trans_t> symm_trans); // const & FIXME

   // This finds bonds between a residue and the protein (in SelAtom).
   // It is used for the environment bonds box.
   // 
   Bond_lines_container(const atom_selection_container_t &SelAtom,
			PPCAtom residue_atoms,
			int n_residue_atoms,
			coot::protein_geometry *protein_geom, // modifiable, currently
			short int residue_is_water_flag,
			float min_dist,
			float max_dist);

   // same as above, except this does symmetry contacts too.
   // 
   Bond_lines_container(const atom_selection_container_t &SelAtom,
			PPCAtom residue_atoms,
			int n_residue_atoms,
			float min_dist,
			float max_dist, 
			short int do_symmetry);

   // This is the one for occupancy and B-factor representation
   // 
   Bond_lines_container (const atom_selection_container_t &SelAtom,
			 bond_representation_type by_occ);

   Bond_lines_container(int col);
   Bond_lines_container(symm_keys key);

   // initial constructor, added to by  addSymmetry_vector_symms from update_symmetry()
   Bond_lines_container() {
      do_bonds_to_hydrogens = 1;  // added 20070629
      b_factor_scale = 1.0;
      if (bonds.size() == 0) { 
	 for (int i=0; i<10; i++) { 
	    Bond_lines a(i);
	    bonds.push_back(a);
	 }
      }
   }

   // arguments as above.
   // 
   // FYI: there is only one element to symm_trans, the is called from
   // them addSymmetry_vector_symms wrapper
   graphical_bonds_container
      addSymmetry(const atom_selection_container_t &SelAtom,
		  coot::Cartesian point,
		  float symm_distance,
		  const std::vector<std::pair<symm_trans_t, Cell_Translation> > &symm_trans,
		  short int symmetry_as_ca_flag,
		  short int symmetry_whole_chain_flag);
   graphical_bonds_container
      addSymmetry_with_mmdb(const atom_selection_container_t &SelAtom,
			    coot::Cartesian point,
			    float symm_distance,
			    const std::vector<std::pair<symm_trans_t, Cell_Translation> > &symm_trans,
			    short int symmetry_as_ca_flag); 

   // FYI: there is only one element to symm_trans, the is called from
   // them addSymmetry_vector_symms wrapper
   graphical_bonds_container
      addSymmetry_calphas(const atom_selection_container_t &SelAtom,
			  const coot::Cartesian &point,
			  float symm_distance,
			  const std::vector<std::pair<symm_trans_t, Cell_Translation> > &symm_trans);

   // FYI: there is only one element to symm_trans, the is called from
   // them addSymmetry_vector_symms wrapper
   graphical_bonds_container
     addSymmetry_whole_chain(const atom_selection_container_t &SelAtom,
			     const coot::Cartesian &point,
			     float symm_distance, 
			     const std::vector <std::pair<symm_trans_t, Cell_Translation> > &symm_trans);

   // symmetry with colour-by-symmetry-operator: If
   // symmetry_whole_chain_flag is set, then if any part of a symmetry
   // related molecule is found, then select the whole chain.
   // 
   std::vector<std::pair<graphical_bonds_container, std::pair<symm_trans_t, Cell_Translation> > >
      addSymmetry_vector_symms(const atom_selection_container_t &SelAtom,
			       coot::Cartesian point,
			       float symm_distance,
			       const std::vector<std::pair<symm_trans_t, Cell_Translation> > &symm_trans,
			       short int symmetry_as_ca_flag,
			       short int symmetry_whole_chain_flag,
			       short int draw_hydrogens_flag);

   std::vector<std::pair<graphical_bonds_container, symm_trans_t> >
      add_NCS(const atom_selection_container_t &SelAtom,
	      coot::Cartesian point,
	      float symm_distance,
	      std::vector<std::pair<coot::coot_mat44, symm_trans_t> > &strict_ncs_mats,
	      short int symmetry_as_ca_flag,
	      short int symmetry_whole_chain_flag);

   graphical_bonds_container
      add_NCS_molecule(const atom_selection_container_t &SelAtom,
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
				   const coot::Cartesian &point,
				   float symm_distance,
				   const std::pair<coot::coot_mat44, symm_trans_t> &strict_ncs_mat);

   graphical_bonds_container make_graphical_bonds() const;

     
   // debugging function
   void check() const; 

   void no_symmetry_bonds();


   // void make_graphical_symmetry_bonds() const;
   graphical_bonds_container make_graphical_symmetry_bonds() const; 

   void check_graphical_bonds() const; 
   void check_static() const; 
   void do_disulphide_bonds(atom_selection_container_t, int imodel);
   void do_Ca_bonds(atom_selection_container_t SelAtom, 
		    float min_dist, float max_dist); 
   coot::my_atom_colour_map_t do_Ca_or_P_bonds_internal(atom_selection_container_t SelAtom,
							const char *backbone_atom_id,
							coot::my_atom_colour_map_t acm,
							float min_dist, float max_dist,
							int bond_colour_type); 
   void do_Ca_plus_ligands_bonds(atom_selection_container_t SelAtom, 
				 float min_dist, float max_dist);
   void do_Ca_plus_ligands_bonds(atom_selection_container_t SelAtom, 
				 float min_dist, float max_dist,
				 int atom_colour_type); 
   void do_colour_by_chain_bonds(const atom_selection_container_t &asc,
				 int draw_hydrogens_flag,
				 short int change_c_only_flag);
   void do_colour_by_molecule_bonds(const atom_selection_container_t &asc,
				    int draw_hydrogens_flag);
   void do_normal_bonds_no_water(const atom_selection_container_t &asc,
				 float min_dist, float max_dist);
   void do_colour_sec_struct_bonds(const atom_selection_container_t &asc,
				   float min_dist, float max_dist); 
   void do_Ca_plus_ligands_colour_sec_struct_bonds(const atom_selection_container_t &asc,
						   float min_dist, float max_dist); 

   atom_selection_container_t
      ContactSel(PPCAtom trans_sel, PSContact contact, int ncontacts) const;
   void do_symmetry_Ca_bonds(atom_selection_container_t SelAtom,
			     symm_trans_t symm_trans);

   void set_verbose_reporting(short int i) { verbose_reporting = i;}

};

#endif /* BOND_LINES_H */
