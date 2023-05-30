/* coords/mmdb-crystal.h
 * 
 * Copyright 2006 by The University of York
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

// -*-c++-*-

#ifndef MMDB_CRYSTAL
#define MMDB_CRYSTAL

#include "Cartesian.h" // uncommented so that r899 compiles for Bill. Strange.
#include "coot-utils/atom-selection-container.hh"

namespace coot { 
   class coot_v4 {
   public:
      std::vector<float> v4;
   };

   class coot_mat44 {
   public:
      std::vector<coot_v4> m;
      coot_mat44() {
	 m.resize(4);
	 for (int i=0; i<4; i++)
	    m[i].v4.resize(4);
      }
      bool is_close_to_unit_matrix() const {
	 float sum_dist = 0.0;
	 for (int i=0; i<4; i++) {
	    for (int j=0; j<4; j++) {
	       if (i==j) {
		  if (i==3)
		     sum_dist += fabsf(m[i].v4[j]);
		  else
		     sum_dist += fabsf(m[i].v4[j]-1);
	       } else{
		  sum_dist += fabsf(m[i].v4[j]);
	       }
	    }
	 }
	 return (sum_dist < 0.001);
      }
   };

   class trans_selection_t {
   public:
      Cartesian front;
      Cartesian back;
      Cartesian left;
      Cartesian right;
      Cartesian top;
      Cartesian bottom;
      bool point_is_in_box(const coot::Cartesian &point) const;
   };
}

std::ostream & operator<<(std::ostream &s, const coot::coot_mat44 &m);

//! A class for symmetry operators
class symm_trans_t { 

   int symm_no, x_shift_, y_shift_, z_shift_;

 public:
   //! symmetry as symm
   std::string symm_as_string;
   symm_trans_t(int n, int x, int y, int z) {
      symm_no = n; x_shift_ = x; y_shift_ = y; z_shift_ = z; mmdb::Mat4Init(mat); };
   explicit symm_trans_t(int idx) { symm_no = idx; x_shift_ = 0; y_shift_ = 0; z_shift_ = 0; }
   symm_trans_t() {};
   mmdb::mat44 mat;

   friend std::ostream & operator<<(std::ostream &s, const symm_trans_t &st);

   //! the operator index (from SYMINFO)
   int isym() const { return symm_no;};
   //! the x-shift
   int x()    const { return x_shift_;};
   //! the y-shift
   int y()    const { return y_shift_;};
   //! the z-shift
   int z()    const { return z_shift_;};
   void add_shift(int xs, int ys, int zs) {
      x_shift_ += xs;
      y_shift_ += ys;
      z_shift_ += zs;
   } 

   //! @return true if this symmetry operator is the identity matrix
   bool is_identity();

   //! symmetry as string
   std::string str(bool expanded_flag) const;

   //! fill mat using mol
   void as_mat44(mmdb::mat44 *mat, mmdb::Manager *mol);

   //! fill m
   void fill_mat(mmdb::Manager *mol) {
      as_mat44(&mat, mol);
   }

};

//! class for the cell translation for generation/application of symmery-related molecules
//!
class Cell_Translation { 

 public: 

   //! unit cell steps
   int us, vs, ws; 
   
   //! constructor
   Cell_Translation() {us=0; vs=0; ws=0;}
   //! constructor
   Cell_Translation(int a, int b, int c);
   friend std::ostream& operator<<(std::ostream &s, Cell_Translation ct);
   //! inversion
   Cell_Translation inv() const { return Cell_Translation(-us, -vs, -ws); }
};

std::string to_string(const std::pair<symm_trans_t, Cell_Translation> &sts);



class molecule_extents_t { 

   // coordinates of the most limiting atoms in the faces. 
   //
   coot::Cartesian front, back, left, right, top, bottom, centre;
   // front, back, minimum and maximum in z;
   // left, right, minimum and maximum in x;
   // bottom, top  minimum and maximum in y;

   mmdb::PPAtom extents_selection;
   float expansion_size_;
   // Grrr.. we cant have a function that returns an mmdb symmetry matrix.
   // So modify it in place.
   void shift_matrix(mmdb::Manager *mol,
		     mmdb::mat44 my_matt,
		     int x_shift, int y_shift, int z_shift,
		     mmdb::mat44 new_matrix) const;

   Cell_Translation atom_sel_cell_trans;  // The reverse
					  // transformation to bring
					  // atom sel current position
					  // to close to the origin.
					  // So, it's the
					  // origin->current_position
					  // cell translation. Set by
					  // the constructor.

 public:

   // expansion_size is typically the symmetry radius
   molecule_extents_t(atom_selection_container_t, float expansion_size);
   ~molecule_extents_t();
   coot::Cartesian get_front() const; 
   coot::Cartesian get_back()  const; 
   coot::Cartesian get_left() const; 
   coot::Cartesian get_right()  const; 
   coot::Cartesian get_top() const; 
   coot::Cartesian get_bottom()  const; 

   Cell_Translation 
      coord_to_unit_cell_translations(coot::Cartesian point,
				      atom_selection_container_t AtomSel) const; 

   // Return vector size 0 when there is no symmetry
   // (GetNumberOfSymOps returns 0)
   // 
   std::vector<std::pair<symm_trans_t, Cell_Translation> >
   which_boxes(coot::Cartesian point,
	       atom_selection_container_t AtomSel,
	       int shift_search_size = 1) const;
   
   std::vector<std::pair<int, symm_trans_t> >
   which_strict_ncs(const coot::Cartesian &centre_pt,
		    atom_selection_container_t &AtomSel,
		    const std::vector<coot::coot_mat44> &strict_ncs_mats,
		    const Cell_Translation &c_t) const;

   // new style

   // use extents to fill transsel, use cryst from mol (not coords of mol)
   coot::trans_selection_t trans_sel_o(mmdb::Manager *mol, const symm_trans_t &symm_trans) const;
   mmdb::PPAtom trans_sel(mmdb::Cryst *my_cryst, symm_trans_t symm_trans) const;
   mmdb::PPAtom trans_sel(mmdb::Manager *mol, const symm_trans_t &symm_trans) const;
   mmdb::PPAtom trans_sel(mmdb::Manager *mol, mmdb::mat44 my_mat,
		     int x_shift, int y_shift, int z_shift) const;


   bool point_is_in_box(const coot::Cartesian &point, mmdb::PPAtom TransSel) const;

   friend std::ostream& operator<<(std::ostream &s, const molecule_extents_t &e);
};



// Let's use symmetry in a properly OO manner, rather than the fortran 
// way that I dislike.  Hide all that here.
// 
class SymmMatrix {

   double mat[4][4];

 public:

   SymmMatrix(); // creates identity matrix.
   SymmMatrix(double** in_mat); // creates from a mmdb::mat44.

   double** getMat() const;
   // float[4][4] testing() const;  // oh.  C++ does not allow us to return
                                    // a mmdb::mat44 - grumble.

   void add_unit_shift(int x, int y, int z);
   friend std::ostream& operator<<(std::ostream&, SymmMatrix);

};


// return an atom selection that has had the symm_trans
// applied to it.
//
mmdb::PPAtom translated_atoms(atom_selection_container_t AtomSel,
			 symm_trans_t symm_trans);

coot::Cartesian translate_atom(atom_selection_container_t AtomSel, 
			 int i, 
			 symm_trans_t symm_trans);

coot::Cartesian translate_atom_with_pre_shift(atom_selection_container_t AtomSel, 
					      int i, 
					      const std::pair<symm_trans_t, Cell_Translation> &symm_trans);

// Tinker with asc (actually, internally, the mmdb::CMMDBCryst of asc)
// 
// Return 1 on success, 0 on failure.
int set_mmdb_cell_and_symm(atom_selection_container_t asc, 
			   std::pair<std::vector<float>, std::string> cell_spgr);
			   

#endif // MMDB_CRYSTAL
