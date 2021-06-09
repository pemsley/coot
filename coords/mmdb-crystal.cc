
/* coords/mmdb-crystal.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 The University of York
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1335, USA
 */

// mmdb related functions that also (used to) depend on clipper.  We
// put them here so that mmdb and mmdb-extras are not exposed to
// clipper complexities.
//
#include <string.h>  // for strlen
#include <string>
#include <vector>
#include <math.h>

#include <sys/types.h>  // stating
#include <sys/stat.h>

#if !defined _MSC_VER
#include <unistd.h>
#else
#define PKGDATADIR "C:/coot/share"
#define snprintf _snprintf
#endif // _MSC_VER

#include "Cartesian.h"

#include <mmdb2/mmdb_manager.h>
#include "mmdb-extras.h"
#include "mmdb.h"
#include "mmdb-crystal.h"
#include "utils/coot-utils.hh" // for Upper

#include "mini-mol/mini-mol.hh"

// Note that with expansion_size small or 0, symmetry operators are
// missed.  We only consider symmetry operators that put the point in
// the tranformed box, BUT there could well be symmetry operators that
// don't have the point in the box but DO have atoms that are 20
// Angstroms away (if the symmetry radius was 21 say).
// 
// The error is most clearly seen with symmetry radius at 100.  Move
// by a few angstroms to the next residue (say) and a whole new other
// set of atoms from other symmetry operators appear.
//
// So if we pass the symmetry radius as the expansion_size, that
// should fix it - which is what we do.
// 
molecule_extents_t::molecule_extents_t(atom_selection_container_t selection,
				       float expansion_size) {

   float atom_x, atom_y, atom_z;
   float max_x, max_y, max_z, min_x, min_y, min_z;
   expansion_size_ = expansion_size;

   max_x = -99999999.9;
   max_y = -99999999.9;
   max_z = -99999999.9;
   
   min_x = 99999999.9;
   min_y = 99999999.9;
   min_z = 99999999.9;
   
   if (selection.n_selected_atoms > 0 ) {

      // we need to reset these, lets not rely on what they used to be.

      for (int i=0; i< selection.n_selected_atoms; i++) {

	 atom_x = selection.atom_selection[i]->x;
	 atom_y = selection.atom_selection[i]->y;
	 atom_z = selection.atom_selection[i]->z;

	 // if there is only one atom, it will be all the limits,
	 // so we don't use the else.
	 //
	 if (atom_x > max_x) {
	    max_x = atom_x;
	    right = coot::Cartesian(atom_x, atom_y, atom_z);
	 }

	 if (atom_x < min_x) {
	    min_x = atom_x; 
	    left = coot::Cartesian(atom_x, atom_y, atom_z);
	 } 

	 if (atom_y > max_y) { 
	    max_y = atom_y;
	    top = coot::Cartesian(atom_x, atom_y, atom_z);
	 } 

	 if (atom_y < min_y) {
	    min_y = atom_y;
	    bottom = coot::Cartesian(atom_x, atom_y, atom_z);
	 } 

	 if (atom_z > max_z) { // back is at max_z;
	    max_z = atom_z;
	    back = coot::Cartesian(atom_x, atom_y, atom_z);
	 } 

	 if (atom_z < min_z) {  // front is at min_z;
	    min_z = atom_z;
	    front = coot::Cartesian(atom_x, atom_y, atom_z);
	 }
      }
   }

   float mid_x, mid_y, mid_z;
   mid_x = (left.get_x()+right.get_x())*0.5;
   mid_y = (bottom.get_y()+top.get_y())*0.5;
   mid_z = (front.get_z()+back.get_z())*0.5;

   // Adjust the extents so that they are on the midpoints of the other axes.
   //
   left   = coot::Cartesian( left.get_x() - expansion_size, mid_y, mid_z);
   right  = coot::Cartesian(right.get_x() + expansion_size, mid_y, mid_z);
   front  = coot::Cartesian(mid_x, mid_y, front.get_z() - expansion_size);
   back   = coot::Cartesian(mid_x, mid_y,  back.get_z() + expansion_size);
   bottom = coot::Cartesian(mid_x, bottom.get_y() - expansion_size, mid_z);
   top    = coot::Cartesian(mid_x, top.get_y()    + expansion_size, mid_z);

//    std::cout << "DEBUG:: left:   " << left   << std::endl;
//    std::cout << "DEBUG:: right:  " << right  << std::endl;
//    std::cout << "DEBUG:: front:  " << front  << std::endl;
//    std::cout << "DEBUG:: back  : " << back   << std::endl;
//    std::cout << "DEBUG:: bottom: " << bottom << std::endl;
//    std::cout << "DEBUG:: top:    " << top    << std::endl;

   // now make the centre for the above coordinates
   // just for reference
   centre = front + back + left + right + top + bottom;
   centre *= .16666666;

   // cout << "centre at: " << centre << endl;
   
   extents_selection = new mmdb::PAtom[6];

   extents_selection[0] = new mmdb::Atom;
   extents_selection[0]->SetCoordinates(front.get_x(), front.get_y(),
					front.get_z(), 1.0 ,99.9);

   extents_selection[1] = new mmdb::Atom; // back is at max_z;
   extents_selection[1]->SetCoordinates(back.get_x(), back.get_y(),
					back.get_z(), 1.0 ,99.9);

   extents_selection[2] = new mmdb::Atom;
   extents_selection[2]->SetCoordinates(left.get_x(), left.get_y(),
					left.get_z(), 1.0 ,99.9);

   extents_selection[3] = new mmdb::Atom;
   extents_selection[3]->SetCoordinates(right.get_x(), right.get_y(),
					right.get_z(), 1.0 ,99.9);

   extents_selection[4] = new mmdb::Atom;
   extents_selection[4]->SetCoordinates(bottom.get_x(), bottom.get_y(),
					bottom.get_z(), 1.0 ,99.9);

   extents_selection[5] = new mmdb::Atom;
   extents_selection[5]->SetCoordinates(top.get_x(), top.get_y(),
					top.get_z(), 1.0 ,99.9);


   coot::minimol::residue res(1, "EXT");
   for (int i=0; i<6; i++) {
      coot::minimol::atom at(" CA ", " C",
			     extents_selection[i]->x,
			     extents_selection[i]->y,
			     extents_selection[i]->z, "", 10.0, 1.0);
      res.addatom(at);
   }

   atom_sel_cell_trans = coord_to_unit_cell_translations(centre, selection);
}

molecule_extents_t::~molecule_extents_t() {

   for (int i=0; i<6; i++) {
      delete extents_selection[i];
   }
   delete [] extents_selection;
} 

coot::Cartesian
molecule_extents_t::get_front() const {
   return front;
}

coot::Cartesian
molecule_extents_t::get_back() const {
   // Jojo
   return back;
} 

coot::Cartesian
molecule_extents_t::get_left() const {
   return left;
}

coot::Cartesian
molecule_extents_t::get_right() const {
   return right;
} 
coot::Cartesian
molecule_extents_t::get_top() const {
   return top;
}

coot::Cartesian
molecule_extents_t::get_bottom() const {
   return bottom;
} 


// So we have a coordinate somewhere in space and the extents of
// the protein.  We want to find the atoms and the bonds of those
// atoms around that point.
//
// We what to find the unit cell that that coordinate is in, so that
// we know which my_matt sets to look over. 
// 
// So we convert the real-space coordinates to fractional and take
// the whole parts to find the unit cell.  This gives us the central
// point, around which we will do +/-x,yz searching using my_matt
// (and all the symmetry of each).
//  That is 6x27xNsymop checks.
//
// We need to know the symmetry of the molecule.  Read it from the molecule.
// Question: if the sg is in the pdb file, does mmdb find it?
//
// Yes it does: MMDBManager->GetSpaceGroup() and we can use that to look up
// symmetry elements (cf filter3.cc). 
// 
// If it is not there, then presume that it is the symmetry of the
// last map that was read in.
// 
//
//



std::ostream& operator<<(std::ostream &s, molecule_extents_t e) {

   s << "front:  " << e.front << std::endl;
   s << "back :  " << e.back  << std::endl;
   s << "left :  " << e.left  << std::endl;
   s << "right:  " << e.right << std::endl;
   s << "top  :  " << e.top   << std::endl;
   s << "bottom: " << e.bottom  << std::endl;

   return s;
} 


//  Let's try again.
//
// For each symm, and 3 trans:
//    Ask:  Is the point in the box?
//       if yes, return that box.
//
// Note that this is probably the slow way of doing it.  As Eugene noted, it
// will be faster to apply the translations/rotations to the search point
// and check against the fixed atoms - and then apply the reverse transformation
// to the atoms at the end.
//
// However, to do that, we first need to know what answer we should be finding
// from that method, hence this method.
// 
// Return vector size 0 when there is no symmetry (GetNumberOfSymOps returns 0)
//
// shift_search_size is an optional argument, default is 1.
// 
std::vector<std::pair<symm_trans_t, Cell_Translation> >
molecule_extents_t::which_boxes(coot::Cartesian point,
				atom_selection_container_t AtomSel,
				int shift_search_size) const {

   std::vector<std::pair<symm_trans_t, Cell_Translation> > symm_trans;
   int n = AtomSel.mol->GetNumberOfSymOps();
   // mmdb::PPAtom trans_selection; old style
   coot::trans_selection_t trans_selection_obj;
   bool in;

   if (n > 0) { 

      // First we move point as close to the origin as we can.  Remember
      // the shift, we'll be needing it.
      //
      // Then generate all the sym ops around that point.
      // Then apply the reverse unit cell(s) shift
      // Then ask if point is in box.
      //    if it is, push that symm_trans on the the returned
      //       symm_trans vector.

      // OK, so what now is the fractional position of the point?
      mmdb::realtype u, v, w, u_cs, v_cs, w_cs;
      int point_unit_cell[3];
      AtomSel.mol->Orth2Frac(point.x(), point.y(), point.z(), u, v, w);
      point_unit_cell[0] = int(u+0.5);
      point_unit_cell[1] = int(v+0.5);
      point_unit_cell[2] = int(w+0.5);
      if (u<0)
	 point_unit_cell[0] -= 1;
      if (v<0)
	 point_unit_cell[1] -= 1;
      if (w<0)
	 point_unit_cell[2] -= 1;

      // ishift is the translation search, translation of the unit
      // cell for each symmetry operator.  For a given search radius
      // (expansion_size_), we want more ishifts for small unit cells.
      // First find the minium cell dimension and compare that to the
      // search radius.
      int ishift=1;
      mmdb::realtype a[6];
      mmdb::realtype vol;
      mmdb::realtype min_cell = 9999999.9; // A
      int orthcode;
      AtomSel.mol->GetCell(a[0], a[1], a[2], a[3], a[4], a[5], vol, orthcode);
      for (int i=0; i<3; i++) 
	 if (a[i]< min_cell)
	    min_cell = a[i];
      ishift = shift_search_size + int((expansion_size_)/min_cell);
//       std::cout << "old ishift " << ishift << " "
// 		<< expansion_size_ << "/" << min_cell << std::endl;
      // a shift of 2 finds the missing symmetry molecule
//       ishift = 2 + int((expansion_size_)/min_cell);
//       std::cout << "new ishift " << ishift << std::endl;

      // Now how about the centre of the atom selection:
      //
      // the AtomSel, after all is the thing to which we apply the
      // symm_trans, we don't want that swinging about all over the
      // place.
      //
      // So, add to the constructor the setting of
      // atom_selection_to_origin shifts, which tries to get the
      // centre of the atom selection as close to the origin as
      // possible, so that when we call trans_sel(), that applies
      // these atom_selection_to_origin shifts before the symop is
      // applied, then re-applies it backwards to the 6 points of a
      // transsel;

      // OK, so point_unit_cell is an integer number of cell shifts
      // (1.6, -1.2, 0.8) -> (2, -1, 1)

      for (int ii=0; ii<n; ii++) { // number of symm ops
	 for (int x_shift= -ishift; x_shift<=ishift; x_shift++) {
	    for (int y_shift= -ishift; y_shift<=ishift; y_shift++) { 
	       for (int z_shift= -ishift; z_shift<=ishift; z_shift++) {

		  // Let's concern ourselves only with extents.
		  // The atom_sel is not important here.
		  // 
		  // we want to generate symm trans of the extents
		  // that have been shifted close to the origin.
		  //
		  // shift extents close to origin,
		  symm_trans_t s_t(ii, x_shift, y_shift, z_shift);
		  
		  // Now add the transformation to get from close
		  // to the origin to point
		  s_t.add_shift(point_unit_cell[0],
				point_unit_cell[1],
				point_unit_cell[2]);
		  
		  // We only pass mol to get to the cryst object (sigh)
		  // We want to move extents close to origin, then
		  // multiply by s_t
		  trans_selection_obj = trans_sel_o(AtomSel.mol, s_t);
		  in = trans_selection_obj.point_is_in_box(point);
		  
		  if (in == 1) {
 		     if ( ! ((ii == 0) &&
 			     (atom_sel_cell_trans.us == x_shift + point_unit_cell[0]) &&
 			     (atom_sel_cell_trans.vs == y_shift + point_unit_cell[1]) &&
 			     (atom_sel_cell_trans.ws == z_shift + point_unit_cell[2]) ) ) {
			//		     if (1) {
			
			s_t.symm_as_string = AtomSel.mol->GetSymOp(ii);
			// coords needs a atom_unit_cell_shift (back
			// to origin) applied to them before s_t is applied
			std::pair<symm_trans_t, Cell_Translation> p(s_t, atom_sel_cell_trans);
			symm_trans.push_back(p);
// 		     } else {
// 			s_t.symm_as_string = AtomSel.mol->GetSymOp(ii);
// 			std::cout << "DEBUG:: rejecting " << s_t << " "
// 				  << "atom_sel_cell_trans " << atom_sel_cell_trans << " "
// 				  << "point_unit_cell" << " " << point_unit_cell[0] << " "
// 				  << point_unit_cell[1] << " "
// 				  << point_unit_cell[2] << " "
// 				  << std::endl;
		     } 
		  }
	       }
	    }
	 }
      }
   }
   return symm_trans;
}


std::vector<std::pair<int, symm_trans_t> >
molecule_extents_t::which_strict_ncs(const coot::Cartesian &centre_pt,
				     atom_selection_container_t &AtomSel,
				     const std::vector<coot::coot_mat44> &strict_ncs_matrices,
				     const Cell_Translation &c_t) const {

   std::vector<std::pair<int, symm_trans_t> > r;
   int n = strict_ncs_matrices.size();
   // mmdb::mat44 m[n]; GNU code
   mmdb::mat44 *m = new mmdb::mat44[n]; // give it back at end;

   for (unsigned int imat=0; imat<strict_ncs_matrices.size(); imat++) { 
      m[imat][0][0] = strict_ncs_matrices[imat].m[0].v4[0];
      m[imat][0][1] = strict_ncs_matrices[imat].m[0].v4[1];
      m[imat][0][2] = strict_ncs_matrices[imat].m[0].v4[2];
      m[imat][1][0] = strict_ncs_matrices[imat].m[1].v4[0];
      m[imat][1][1] = strict_ncs_matrices[imat].m[1].v4[1];
      m[imat][1][2] = strict_ncs_matrices[imat].m[1].v4[2];
      m[imat][2][0] = strict_ncs_matrices[imat].m[2].v4[0];
      m[imat][2][1] = strict_ncs_matrices[imat].m[2].v4[1];
      m[imat][2][2] = strict_ncs_matrices[imat].m[2].v4[2];
      m[imat][0][3] = strict_ncs_matrices[imat].m[0].v4[3];  // t
      m[imat][1][3] = strict_ncs_matrices[imat].m[1].v4[3];  // t
      m[imat][2][3] = strict_ncs_matrices[imat].m[2].v4[3];  // t
   }

   mmdb::Atom atom;
   mmdb::Atom trans_atom;
   mmdb::Atom tmp_atom;
   atom.SetCoordinates(centre_pt.x(), centre_pt.y(), centre_pt.z(), 1.0, 10.0);
   mmdb::realtype diff_x, diff_y, diff_z, u, v, w, u_cs, v_cs, w_cs;
   symm_trans_t  symm_trans_this;

   for (int ii=0; ii<n; ii++) {
      trans_atom.Copy(&atom); // atom not modified.

      trans_atom.Transform(m[ii]);

      // OK, so what now is the fractional difference between atom and trans_atom?
      //
      diff_x = centre_pt.x() - trans_atom.x;
      diff_y = centre_pt.y() - trans_atom.y;
      diff_z = centre_pt.z() - trans_atom.z;

      AtomSel.mol->Orth2Frac(diff_x, diff_y, diff_z, u, v, w); // fill u, v, w.
      u_cs = rint(u);
      v_cs = rint(v);
      w_cs = rint(w);

      float min_dist = 99999999999.9;
      float b, dist;
      mmdb::mat44 shifted_mat;
      
      for (int x_shift= -1; x_shift<= 1; x_shift++) { 
	 for (int y_shift= -1; y_shift<= 1; y_shift++) { 
	    for (int z_shift= -1; z_shift<= 1; z_shift++) {
	       
	       if ( ! (x_shift == 0 && y_shift == 0 && z_shift == 0 && ii == 0) ) { 
		  tmp_atom.Copy(&trans_atom); // not modify trans_atom
		  
		  // Add x_shift, y_shift, z_shift to my_matt
		  // We need cryst info to do this, so pass the mol.
		  shift_matrix(AtomSel.mol, m[ii], x_shift, y_shift, z_shift, shifted_mat);
		  tmp_atom.Transform(shifted_mat);
		  b = 0.0; 
		  dist = tmp_atom.x - centre_pt.x();
		  b += dist*dist;
		  dist = tmp_atom.y - centre_pt.y();
		  b += dist*dist;
		  dist = tmp_atom.z - centre_pt.z();
		  b += dist*dist;
		  if (b < min_dist) {
		     min_dist = b;
		     symm_trans_this = symm_trans_t(ii, x_shift, y_shift, z_shift);
		     symm_trans_this.symm_as_string = AtomSel.mol->GetSymOp(ii);
		  }
	       }
	    }
	 }
      }

      // so now we have sym_trans_this, which corresponds to the
      // closest approach of the NCS molecules to the centre point.
      int px_shift = symm_trans_this.x();
      int py_shift = symm_trans_this.y();
      int pz_shift = symm_trans_this.z();
      mmdb::PPAtom trans_selection;
      int in;
      
      for (int x_shift=px_shift-1; x_shift<=px_shift+1; x_shift++) { 
	 for (int y_shift=py_shift-1; y_shift<=py_shift+1; y_shift++) { 
	    for (int z_shift=pz_shift-1; z_shift<=pz_shift+1; z_shift++) { 
	       
	       if ( ! (x_shift == 0 && y_shift == 0 && z_shift == 0 && ii == 0) ) { 

		  // return transformed the 6 atom limits:
		  trans_selection = trans_sel(AtomSel.mol, m[ii], x_shift, y_shift, z_shift);
		  in = point_is_in_box(centre_pt, trans_selection);
		  if (in == 1) {
// 		     std::cout << "shifts: " << x_shift << " " << y_shift << " "
// 			       << z_shift << " is IN box\n";
		     symm_trans_t s_t(ii, x_shift, y_shift, z_shift);
		     r.push_back(std::pair<int, symm_trans_t> (ii, s_t));
// 		  } else {
// 		     std::cout << "shifts: " << x_shift << " " << y_shift << " "
// 			       << z_shift << " is OUT of box\n";
		  }
		  for (int it=0; it<6; it++)  
		     delete trans_selection[it];
		  delete [] trans_selection;
	       }
	    }
	 }
      }
   }
   delete [] m;
   return r;
}

// Is this nonsense really necessary?
void
molecule_extents_t::shift_matrix(mmdb::Manager *mol,
				 mmdb::mat44 my_matt,
				 int x_shift, int y_shift, int z_shift,
				 mmdb::mat44 new_matrix) const {
   mmdb::mat44 amat;
   mol->GetTMatrix(amat, 0, x_shift, y_shift, z_shift);
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++)
	 new_matrix[i][j] = my_matt[i][j];
   
   for (int i=0; i<3; i++) 
	 new_matrix[i][3] = amat[i][3];

}



// 
Cell_Translation
molecule_extents_t::coord_to_unit_cell_translations(coot::Cartesian point,
						    atom_selection_container_t AtomSel) const {

   // fractionalize point and return a set of integers for iv, iv, iw.
   //
   // we don't need to use clipper, we can use:
   //
   // mmdb::Manager::Orth2Frac Orthogonal-to-fractional transformation
   // of coordinates and mmdb::Manager::Frac2Orth.
   //
   // mind-bogglingly ugly code necessary to use an mmdb::CMMDBCryst (we can't just use
   // the simple object because the destructor will crush it).
   // 
   // Pmmdb::CMMDBCryst my_cryst_p = (mmdb::CMMDBCryst *) &(AtomSel.mol->get_cell());

   mmdb::realtype u, v, w;

   // 20031105 Well, we don't need Pmmdb::CMMDBCryst my_cryst_p, we can use
   // the function at the mmdb manager level
   // AtomSel.mol->Orth2Frac(point.get_x(), point.get_y(), point.get_z(),
   //                        u, v, w);

   // 
   //
   AtomSel.mol->Orth2Frac(point.get_x(), point.get_y(), point.get_z(),
			  u, v, w);

   // doub ud = u
   int iu = int (rint(u));
   int iv = int (rint(v));
   int iw = int (rint(w));  // 0.6 -> 1 and -0.6 -> -1
   
   return Cell_Translation (iu, iv, iw);

} 

// Using a const Object & conventionally means that you would normally
// pass a reference to the object, const means that we do not modify
// the Object.
//
// I want to use const mmdb::CMMDBCryst &my_cryst, because the argument is
// not modified, and the only reason I am passing a reference is
// because I can't pass the object, because MMDB crushes it with its
// (IMHO) bogus constructors/destructors (grrr).  But I can't use
// const because GetTMatrix() is not a const function.  Here is a bit
// of my screaming at this concept: ARRRRRRRGGGHHH!
//
// I tried getting round this by passing the mmdb::mat44, but that gives problems
// because a mmdb::mat44 is not a class it is a typedef mmdb::realtype[4][4] and you
// can't pass arguments and create functions that return arguments of
// that type.  Arrrrgh!
//
// I am also angry because it just takes *so long* to track down problems
// with mmdb.
//
// And here's another thing: look at atom1.Copy(mmdb::PAtom atom2).  Now,
// does that copy inside to out (atom2 to atom1) or the other way
// round (atom 1 to atom2)?  Impossible to tell without const.  It
// would have been so much easier to say atom2 = atom1.Copy(); but no,
// in mmdb functions don't return objects - grrrr.
//
//
// 20031105: Actually, a lot of this irritation is misplaced.  You
// don't need to have access to Pmmdb::CMMDBCryst, you can do it through the
// MMDBManager top-level.  Not elegant, obviously, but not worth the
// rant.
// 
mmdb::PPAtom 
molecule_extents_t::trans_sel(mmdb::Cryst *my_cryst,
			      symm_trans_t symm_trans) const {

   mmdb::Atom atom;
   mmdb::PPAtom trans_selection = new mmdb::PAtom[6];
   mmdb::mat44 my_matt;
   
   // Modify my_matt so that it is a coordinate transformation
   // matrix.
   //
   my_cryst->GetTMatrix(my_matt, symm_trans.isym(), symm_trans.x(),
			symm_trans.y(), symm_trans.z());

   for (int ii=0; ii<6; ii++) {
      trans_selection[ii] = new mmdb::Atom;

      trans_selection[ii]->SetCoordinates(extents_selection[ii]->x,
					  extents_selection[ii]->y,
					  extents_selection[ii]->z,
					  1.0, 99.9);
      trans_selection[ii]->Transform(my_matt);
   }
   return trans_selection;
}

mmdb::PPAtom
molecule_extents_t::trans_sel(mmdb::Manager *mol,
			      const symm_trans_t &symm_trans) const {

   mmdb::Atom atom;
   mmdb::PPAtom trans_selection = new mmdb::PAtom[6];
   mmdb::mat44 my_matt;
   
   // Modify my_matt so that it is a coordinate transformation
   // matrix.
   //
   mol->GetTMatrix(my_matt, symm_trans.isym(), symm_trans.x(),
		   symm_trans.y(), symm_trans.z());

   for (int ii=0; ii<6; ii++) {
      trans_selection[ii] = new mmdb::Atom;

      trans_selection[ii]->SetCoordinates(extents_selection[ii]->x,
					  extents_selection[ii]->y,
					  extents_selection[ii]->z,
					  1.0, 99.9);
      trans_selection[ii]->Transform(my_matt);
   }
   return trans_selection;
   

}

// Return the 6 modified coordinates of the extents in a mmdb::PPAtom.  The
// coordinates are transformed by the my_mat matrix and are shifted by
// the cell shifts.
// 
mmdb::PPAtom
molecule_extents_t::trans_sel(mmdb::Manager *mol, mmdb::mat44 my_mat,
			      int x_shift, int y_shift, int z_shift) const {

   mmdb::Atom atom;
   mmdb::PPAtom trans_selection = new mmdb::PAtom[6];
   
   // Modify my_matt so that it is a coordinate transformation
   // matrix.
   //

   mmdb::mat44 amat;
   mol->GetTMatrix(amat, 0, x_shift, y_shift, z_shift);
   for (int i=0; i<3; i++) 
      for (int j=0; j<3; j++)
	 amat[i][j] = my_mat[i][j];

   
//     std::cout << "amat is: " << std::endl
// 	<< amat[0][0] << " "  << amat[0][1] << " "
// 	<< amat[0][2] << " "  << amat[0][3] << " "  << std::endl
// 	<< amat[1][0] << " "  << amat[1][1] << " "
// 	<< amat[1][2] << " "  << amat[1][3] << " "  << std::endl
// 	<< amat[2][0] << " "  << amat[2][1] << " "
// 	<< amat[2][2] << " "  << amat[2][3] << " "  << std::endl
// 	<< amat[3][0] << " "  << amat[3][1] << " "
// 	<< amat[3][2] << " "  << amat[3][3] << " "  << std::endl; 

// a sensible matrix:
// 0 -1 0 -280.056 
// 1  0 0      0
// 0  0 1    -22.573 
// 0  0 0      1 

   for (int ii=0; ii<6; ii++) {
      trans_selection[ii] = new mmdb::Atom;

      trans_selection[ii]->SetCoordinates(extents_selection[ii]->x,
					  extents_selection[ii]->y,
					  extents_selection[ii]->z,
					  1.0, 99.9);
      trans_selection[ii]->Transform(amat);
   }
   return trans_selection;
}

// use extents to fill transsel
//
// We only pass mol to get to the cryst object (sigh)
// We want to move extents close to origin, then
// multiply by s_t
//
coot::trans_selection_t
molecule_extents_t::trans_sel_o(mmdb::Manager *mol, const symm_trans_t &symm_trans) const {

   coot::trans_selection_t t;
   mmdb::Atom atom;
   mmdb::mat44 my_matt;
   mmdb::mat44 to_origin_matt;
   
   mol->GetTMatrix(my_matt, symm_trans.isym(), symm_trans.x(),
		   symm_trans.y(), symm_trans.z());

   mol->GetTMatrix(to_origin_matt, 0,
		   -atom_sel_cell_trans.us,
		   -atom_sel_cell_trans.vs,
		   -atom_sel_cell_trans.ws);
		   
   atom.Copy(extents_selection[0]);
   atom.Transform(to_origin_matt);
   atom.Transform(my_matt);
   t.front = coot::Cartesian(atom.x, atom.y, atom.z);
   
   atom.Copy(extents_selection[1]);
   atom.Transform(to_origin_matt);
   atom.Transform(my_matt);
   t.back =  coot::Cartesian(atom.x, atom.y, atom.z);

   atom.Copy(extents_selection[2]);
   atom.Transform(to_origin_matt);
   atom.Transform(my_matt);
   t.left =  coot::Cartesian(atom.x, atom.y, atom.z);
   
   atom.Copy(extents_selection[3]);
   atom.Transform(to_origin_matt);
   atom.Transform(my_matt);
   t.right = coot::Cartesian(atom.x, atom.y, atom.z);
   
   atom.Copy(extents_selection[4]);
   atom.Transform(to_origin_matt);
   atom.Transform(my_matt);
   t.bottom = coot::Cartesian(atom.x, atom.y, atom.z);

   atom.Copy(extents_selection[5]);
   atom.Transform(to_origin_matt);
   atom.Transform(my_matt);
   t.top =    coot::Cartesian(atom.x, atom.y, atom.z);
   
   return t;
} 


// Not used?
// 
bool
molecule_extents_t::point_is_in_box(const coot::Cartesian &point, mmdb::PPAtom TransSel) const { 

   // front back left right bottom top
   //      z         x            y
   // 
   coot::Cartesian  front(TransSel[0]->x, TransSel[0]->y, TransSel[0]->z);
   coot::Cartesian   back(TransSel[1]->x, TransSel[1]->y, TransSel[1]->z);
   coot::Cartesian   left(TransSel[2]->x, TransSel[2]->y, TransSel[2]->z);
   coot::Cartesian  right(TransSel[3]->x, TransSel[3]->y, TransSel[3]->z);
   coot::Cartesian bottom(TransSel[4]->x, TransSel[4]->y, TransSel[4]->z);
   coot::Cartesian    top(TransSel[5]->x, TransSel[5]->y, TransSel[5]->z);

   coot::Cartesian back_to_front = front - back;
   coot::Cartesian left_to_right = right - left;
   coot::Cartesian bottom_to_top = top - bottom;
   
   coot::Cartesian back_to_point   = point - back;
   coot::Cartesian left_to_point   = point - left;
   coot::Cartesian bottom_to_point = point - bottom;

   coot::Cartesian front_to_point = point - front;
   coot::Cartesian right_to_point = point - right;
   coot::Cartesian top_to_point   = point - top;

   if (coot::dot_product(back_to_front, back_to_point) >= 0.0) {
      if (coot::dot_product(left_to_right, left_to_point) >= 0.0) { 
	 if (coot::dot_product(bottom_to_top, bottom_to_point) >=0.0) {

	    //
	    if (coot::dot_product(back_to_front, front_to_point) <= 0.0) {
	       if (coot::dot_product(left_to_right, right_to_point) <= 0.0) {
		  if (coot::dot_product(bottom_to_top, top_to_point) <= 0.0) {

// 		     cout << "hit: for point: " << point << endl;
// 		     cout << "front:  " << front << endl;
// 		     cout << "back :  " << back  << endl;
// 		     cout << "left :  " << left  << endl;
// 		     cout << "right:  " << right << endl;
// 		     cout << "top  :  " << top   << endl;
// 		     cout << "bottom: " << bottom  << endl;

		     return 1;
		  }
	       }
	    }
	 }
      }
   }
   return 0;
}



bool
coot::trans_selection_t::point_is_in_box(const coot::Cartesian &point) const {

   coot::Cartesian back_to_front = front - back;
   coot::Cartesian left_to_right = right - left;
   coot::Cartesian bottom_to_top = top - bottom;
   
   coot::Cartesian back_to_point   = point - back;
   coot::Cartesian left_to_point   = point - left;
   coot::Cartesian bottom_to_point = point - bottom;

   coot::Cartesian front_to_point = point - front;
   coot::Cartesian right_to_point = point - right;
   coot::Cartesian top_to_point   = point - top;

   if (coot::dot_product(back_to_front, back_to_point) >= 0.0) {
      if (coot::dot_product(left_to_right, left_to_point) >= 0.0) { 
	 if (coot::dot_product(bottom_to_top, bottom_to_point) >=0.0) {

	    if (coot::dot_product(back_to_front, front_to_point) <= 0.0) {
	       if (coot::dot_product(left_to_right, right_to_point) <= 0.0) {
 		  if (coot::dot_product(bottom_to_top, top_to_point) <= 0.0) {
 		     return 1;
 		  }
 	       }
 	    }
	 }
      }
   }
   return 0;
}


std::ostream & operator<<(std::ostream &s, const symm_trans_t &t) {
   s << "symm: " << t.symm_as_string << " (op-idx: "  << t.isym() << ") trans: "
     << t.x_shift_ << " " << t.y_shift_ << " " << t.z_shift_;
   return s;
}


Cell_Translation::Cell_Translation(int a, int b, int c) {
   //
   us = a;
   vs = b;
   ws = c;
}

std::ostream& operator<<(std::ostream &s, Cell_Translation ct) {

   s << "Cell Trans: (" << ct.us << " " << ct.vs << " " << ct.ws << ")";
   return s;

} 


double **
SymmMatrix::getMat() const  {
   return (double**) (mat);
}

// creates from a mmdb::mat44.
SymmMatrix::SymmMatrix(double** in_mat) {

   for (int ii=0; ii<4; ii++) { 
      for (int jj=0; jj<4; jj++) {
	 mat[ii][jj] = in_mat[ii][jj];
      }
   }
}

void
SymmMatrix::add_unit_shift(int x, int y, int z) {

   mat[0][3] += x; 
   mat[1][3] += y; 
   mat[2][3] += z; 
} 

std::ostream& operator<<(std::ostream& s, SymmMatrix m) {

   s << m.mat[0][0] << " "  << m.mat[0][1] << " "
     << m.mat[0][2] << " "  << m.mat[0][3] << " "  << std::endl
     << m.mat[1][0] << " "  << m.mat[1][1] << " "
     << m.mat[1][2] << " "  << m.mat[1][3] << " "  << std::endl
     << m.mat[2][0] << " "  << m.mat[2][1] << " "
     << m.mat[2][2] << " "  << m.mat[2][3] << " "  << std::endl
     << m.mat[3][0] << " "  << m.mat[3][1] << " "
     << m.mat[3][2] << " "  << m.mat[3][3] << " "  << std::endl;
   return s;
}
   
bool
symm_trans_t::is_identity() {

   if ( (symm_no == 0) && (x_shift_ == 0) &&
	(y_shift_ == 0) && (z_shift_ == 0)) {
      return 1;
   } else {
      return 0;
   }
}

//
std::string
symm_trans_t::str(short int expanded_flag) const {

   //
   std::string b; 
   if (expanded_flag) {
      b = coot::util::Upper(symm_as_string);
   } else {
      b = " #s ";
      b += coot::util::int_to_string(symm_no+1);
   }
   b += " + (";
   b += coot::util::int_to_string(x());
   b += " ";
   b += coot::util::int_to_string(y());
   b += " ";
   b += coot::util::int_to_string(z());
   b += ")";
   return b;
}

std::string
to_string(const std::pair<symm_trans_t, Cell_Translation> &sts) {

   std::string b;
   b = coot::util::Upper(sts.first.symm_as_string);
   b += " + (";
   b += coot::util::int_to_string(sts.first.x());
   b += " ";
   b += coot::util::int_to_string(sts.first.y());
   b += " ";
   b += coot::util::int_to_string(sts.first.z());
   b += ") ";
   b += "& {";
   b += " ";
   b += coot::util::int_to_string(sts.second.us);
   b += " ";
   b += coot::util::int_to_string(sts.second.vs);
   b += " ";
   b += coot::util::int_to_string(sts.second.ws);
   b += "}";
   return b;
} 


// return an atom selection that has had the symm_trans
// applied to it.
//
mmdb::PPAtom
translated_atoms(atom_selection_container_t AtomSel,
		symm_trans_t symm_trans) {

   mmdb::mat44 my_matt;
   int err = AtomSel.mol->GetTMatrix(my_matt, symm_trans.isym(), symm_trans.x(),
				      symm_trans.y(), symm_trans.z());
   
   if (err != 0) {
      std::cout << "!!!!!!!!!!!!!! something BAD with mmdb::CMMDBCryst.GetTMatrix"
	   << std::endl;
   }
      
   mmdb::PPAtom trans_selection = new mmdb::PAtom[AtomSel.n_selected_atoms];
   for (int ii=0; ii<AtomSel.n_selected_atoms; ii++) {

      trans_selection[ii] = new mmdb::Atom;
      trans_selection[ii]->Copy(AtomSel.atom_selection[ii]);
      trans_selection[ii]->Transform(my_matt);
      trans_selection[ii]->SetResidue(    AtomSel.atom_selection[ii]->GetResidue());
      //if (ii == 10) {
      //	 cout << (trans_selection[ii]) << " vs. " << (AtomSel.atom_selection[ii]) << std::endl;
      //}

   }
   return trans_selection;
}

coot::Cartesian translate_atom(atom_selection_container_t AtomSel, int ii,
			 symm_trans_t symm_trans) {
   //
   mmdb::mat44 my_matt;
   
   int err = AtomSel.mol->GetTMatrix(my_matt,
				     symm_trans.isym(),
				     symm_trans.x(),
				     symm_trans.y(),
				     symm_trans.z());
   
   if (err != 0) {
      std::cout << "!!!!!!! something BAD with mmdb::CMMDBCryst.GetTMatrix in "
	   << "coot::Cartesian translate_atom(..)" 
	   << std::endl;
   }

   mmdb::PAtom trans_atom = new mmdb::Atom;

   trans_atom->Copy(AtomSel.atom_selection[ii]);
   trans_atom->Transform(my_matt);

   coot::Cartesian c(trans_atom->x, trans_atom->y, trans_atom->z);
   delete trans_atom;
   return c;

}

coot::Cartesian
translate_atom_with_pre_shift(atom_selection_container_t AtomSel, int ii,
			      const std::pair<symm_trans_t, Cell_Translation> &symm_trans) {
   //
   mmdb::mat44 my_matt;
   mmdb::mat44 pre_shift_matt;
   
   int err = AtomSel.mol->GetTMatrix(my_matt,
				     symm_trans.first.isym(),
				     symm_trans.first.x(),
				     symm_trans.first.y(),
				     symm_trans.first.z());
   
   int err2 = AtomSel.mol->GetTMatrix(pre_shift_matt, 0,
				      -symm_trans.second.us,
				      -symm_trans.second.vs,
				      -symm_trans.second.ws);
   
   if (err != 0) {
      std::cout << "!!!!!!! something BAD with mmdb::CMMDBCryst.GetTMatrix in "
		<< "coot::Cartesian translate_atom(..)" 
		<< std::endl;
   }

   mmdb::Atom trans_atom;

   trans_atom.Copy(AtomSel.atom_selection[ii]);
   trans_atom.Transform(pre_shift_matt);
   trans_atom.Transform(my_matt);

   coot::Cartesian c(trans_atom.x, trans_atom.y, trans_atom.z);
   return c;
}


int set_mmdb_cell_and_symm(atom_selection_container_t asc, 
			   std::pair<std::vector<float>, std::string> cell_spgr) { 

   int istat = 0; 
   if (cell_spgr.first.size() == 6) { 
      std::vector<float> a = cell_spgr.first; // short name
      asc.mol->SetCell(a[0], a[1], a[2], a[3], a[4], a[5]);
      asc.mol->SetSpaceGroup((char *)cell_spgr.second.c_str());
      std::cout << "successfully set cell and symmetry" << std::endl;
      istat = 1;
   } else { 
      std::cout << "WARNING:: failure to set cell on this molecule" << std::endl;
   } 
   return istat;
}



atom_selection_container_t read_standard_residues() {

   std::string standard_env_dir = "COOT_STANDARD_RESIDUES";
   atom_selection_container_t standard_residues_asc;
   
   const char *filename = getenv(standard_env_dir.c_str());
   if (! filename) {

      std::string standard_file_name = PKGDATADIR;
      standard_file_name += "/";
      standard_file_name += "standard-residues.pdb";

      struct stat buf;
      int status = stat(standard_file_name.c_str(), &buf);  
      if (status != 0) { // standard-residues file was not found in
			 // default location either...
	 std::cout << "WARNING: environment variable for standard residues ";
	 std::cout << standard_env_dir << "\n";
	 std::cout << "         is not set.";
	 std::cout << " Mutations will not be possible\n";
	 // mark as not read then:
	 standard_residues_asc.read_success = 0;
	 // std::cout << "DEBUG:: standard_residues_asc marked as
	 // empty" << std::endl;
      } else { 
	 // stat success:
	 standard_residues_asc = get_atom_selection(standard_file_name, true, false, false);
      }
   } else { 
      standard_residues_asc = get_atom_selection(filename, true, false, false);
   }

   return standard_residues_asc;
}

