/* skeleton/BuildCas.cc
 * 
 * Copyright 2001, 2002, 2003 by The University of York
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

#include <string.h>
#include "clipper/mmdb/clipper_mmdb.h" // for clipper mmdb stuff (obviously)
                                       // convert_to_atoms_internal().

#include "compat/coot-sysdep.h"
#include "BuildCas.h"  // contains Cartesian.h needed by graphical_skel.h

#include "graphical_skel.h"


// #include <GL/glut.h> // for timings. 

// MOVEME:
// 
// Not really the place for this I suspect.  
// 
coot::Cartesian
average_Cartesians(std::vector<coot::Cartesian> c) {

   if (c.size() == 0) {

      std::cout << "WARNING: averaging zero Cartesian, returning default"
	   << std::endl;
      return coot::Cartesian();

   } else {
      
      float sum_x = 0;
      float sum_y = 0;
      float sum_z = 0;

      for (unsigned int ii=0; ii< c.size(); ii++) { 
	 sum_x += c[ii].x(); 
	 sum_y += c[ii].y(); 
	 sum_z += c[ii].z();
      }
      
      return coot::Cartesian(sum_x/c.size(), sum_y/c.size(), sum_z/c.size());
   }
   
}

void 
BuildCas::set_density_map_and_cut(clipper::Xmap<float> &map_in, float map_cut) { 

   d_map_p = &map_in; 
   map_cut_off = map_cut; 
} 


void 
BuildCas::setup_internal() {

   // angle_info.setup_angle_torsion_table(5,5);
   n_fitted_in_current_segment = 0; 
   i_current_build = 0; 
   i_max_build = 0; // i_max_build is not the maximum number of
                    // builds, but the number of the top build in the
                    // build array.
   
   grid_dependent_distance_param = 0.3; // should be set depending on
   // grid_sampling();
   
   segment_map_filled = 0; // not filled
   branch_point_have_been_expanded_flag = 0; 
   
   // 
   treenodemap_is_filled = 0; // not filled.
   
   // 
   expansion_centre_is_set = 0;   // not set. 

   
}




// diff_residue_flag tells us that we should put each atom into a separate residue, 
// rather than bunging 10 atoms into a residue.
// 
atom_selection_container_t
BuildCas::convert_to_atoms_internal(clipper::Spacegroup spg,
				    clipper::Cell cell, 
				    const std::vector<coot::Cartesian> &c,
				    short int diff_residue_flag,
				    std::string molecule_name) const
{
   // we need a MyMMDBManager (not just an mmdb::Manager) because it
   // become part of an atom_selection_container_t. 
   // 
   mmdb::Manager *MMDBManager;

   //   Make routine initializations
   //
   mmdb::InitMatType();

   // MyMMDBManager = new Mymmdb::Manager();
   // MMDBManager = new Mymmdb::Manager;

   // There were problems converting between clipper space group
   // and mmdb space group (in the commented out code, in the case
   // of "R 3 2" - which mmdb did not recognise.
   //
   // So, we use a clipper::MMDB which is a wrapper for mmdb::Manager
   // (and other things). 
   //
   // We create a "new" object so that it (the mmdb::Manager) does not
   // get destroyed, but note that we cannot now delete the clipper::MMDB
   // (which is not ideal).
   //
   // Consider making the clipper::MMDB*clmmdb a class member.  Then
   // we can have access to it and delete it.
   // 
   
//    clipper::MMDB* clmmdb = new clipper::MMDB( spg, cell );
//    MMDBManager = reinterpret_cast<Mymmdb::Manager*>(clmmdb->pcmmdbmanager());

   clipper::MMDBManager *clmmdb;
   clmmdb = new clipper::MMDBManager;
   clmmdb->set_spacegroup(spg);
   clmmdb->set_cell(cell);

   MMDBManager = clmmdb;

//    cout << "step 0: " << endl; 
//    cout << "DEBUG: There are "
// 	<< MMDBManager->get_cell_p()->GetNumberOfSymOps()
// 	<< " sym ops" << endl; 

//    for (int isym = 0; isym < MMDBManager->get_cell_p()->GetNumberOfSymOps(); isym++) {
//       mmdb::mat44 my_matt;
//       int err2 = MMDBManager->get_cell_p()->GetTMatrix(my_matt, isym, 0, 0, 0); 
      
//       if (err2 != 0) {
// 	 cout << "!! something BAD with MMDB mmdb::CMMDBCryst.GetTMatrix in convert_to_atoms_internal"
// 	      << endl;
//       } else {
// 	 cout << "DEBUG: symop " << isym << " clipper::MMDB seems OK..." << endl;
//       }
//    }



   atom_selection_container_t atom_selection_container;

   mmdb::Chain* chain_p = new mmdb::Chain;
   chain_p->SetChainID ( "X" );  // new chain ID

//    // start of with space for the first residue (actually, the only residue)
//    // 
//    mmdb::Residue* res_p = new mmdb::Residue;
//    res_p->seqNum = 1;   // set the residue number to 1;

//    strcpy(res_p->name, "SKELETONBONE"); 

   int iadd; // mmdb error flag on atom addition


//    cout << "step 1: " << endl; 
//    cout << "DEBUG: There are "
// 	<< MMDBManager->get_cell_p()->GetNumberOfSymOps()
// 	<< " sym ops" << endl; 

//    for (int isym = 0; isym < MMDBManager->get_cell_p()->GetNumberOfSymOps(); isym++) {
//       mmdb::mat44 my_matt;
//       int err2 = MMDBManager->get_cell_p()->GetTMatrix(my_matt, isym, 0, 0, 0); 
      
//       if (err2 != 0) {
// 	 cout << "!! something BAD with MMDB mmdb::CMMDBCryst.GetTMatrix in convert_to_atoms_internal"
// 	      << endl;
//       } else {
// 	 cout << "DEBUG: symop " << isym << " clipper::MMDB seems OK..." << endl;
//       }
//    }

   int i_atom_loop_count = 0; 
   int i_res_add = 0; 
   mmdb::Residue* res_p = 0;

   std::cout << "we were passed " << c.size() << " atoms to convert " << std::endl; 
   for (unsigned int ii=0; ii<c.size(); ii++) { 
      // for (int ii=0; ii<c.size(); ii++) {  // or 230 for testing

      if (i_atom_loop_count == 0 || diff_residue_flag == 1) { 
	    res_p = new mmdb::Residue;
	    res_p->seqNum = 1 + i_res_add;   // set the residue number to 1 + ... 
	    strcpy(res_p->name, molecule_name.c_str()); 
	    chain_p->AddResidue ( res_p );
      }
   
      mmdb::PAtom atom = new mmdb::Atom;

      atom->SetCoordinates(c[ii].x(), c[ii].y(), c[ii].z(), 1.0, 99);
      atom->SetResidue(res_p); 
      
      atom->SetAtomName    ( " CA " );  // it has to be a PDB name!
      atom->SetElementName ( " C"   );
      
      iadd = res_p->AddAtom(atom);
      
      if ( iadd < 0 ) {
	 std::cout << "Atom  Addition error: " << iadd << std::endl;
      }

      i_atom_loop_count++; 
      if (diff_residue_flag == 1) { 
	 i_res_add++;
      } else {
	 if (i_atom_loop_count == 10) { // 10 atoms per residue in mmdb; 
	    i_res_add++; 
	    i_atom_loop_count = 0; 
	 }
      }
   }



   mmdb::Model* model_p = new mmdb::Model;
   model_p->AddChain(chain_p);

   MMDBManager->AddModel(model_p);
   
   MMDBManager->PDBCleanup ( mmdb::PDBCLEAN_SERIAL | mmdb::PDBCLEAN_INDEX );

   /*
   MMDBManager->SetCell(l1.cell().descr().a(),
			l1.cell().descr().b(),
			l1.cell().descr().c(),
			l1.cell().descr().alpha()*RADTODEG,
			l1.cell().descr().beta() *RADTODEG,
			l1.cell().descr().gamma()*RADTODEG, 1);
   
  // cout << "In convert_to_atoms O: trying to set spacegroup: " 
  //	<< l1.spacegroup().spgr_name().c_str() << endl;

	// int ierr = MMDBManager->SetSpaceGroup((char *)l1.spacegroup().spgr_name().c_str());

   if (ierr == 0) {
      cout << "ierr = 0 (OK) for setting space group: " <<
	 l1.spacegroup().spgr_name() << endl;
   } else { 
      cout << endl <<  "!! !! !! WARNING !! !! !! error setting space group: " 
	   << l1.spacegroup().spgr_name() << endl << endl;
   }
   */

   mmdb::PPAtom SelAtom;
   
   int selHnd = MMDBManager->NewSelection();
   int nSelAtoms;
   MMDBManager->SelectAtoms(selHnd, 0,"*",mmdb::ANY_RES,"*",mmdb::ANY_RES,
			    "*","*",  // EndInsertionCode, RNames
			    "*","*",  // ANames, Elements
			    "*" );    // Alternate locations.


   MMDBManager->GetSelIndex(selHnd, SelAtom, nSelAtoms);
   std::cout << nSelAtoms << " atoms selected." << std::endl;

   // give it a mol (mmdbmanager), the number of atom and the selection.

   atom_selection_container.atom_selection = SelAtom;
   atom_selection_container.n_selected_atoms = nSelAtoms;
   atom_selection_container.mol = MMDBManager;

   // debug:  test the symmetry operators
   //

//    cout << "In convert_to_atoms: there are " 
// 	<< atom_selection_container.mol->GetNumberOfSymOps() 
// 	<< " sym ops" << endl;




   
   // now test that clipper casting:
   //
//    cout << "DEBUG: There are "
// 	<< atom_selection_container.mol->get_cell_p()->GetNumberOfSymOps()
// 	<< " sym ops" << endl; 
//    mmdb::mat44 my_matt;
//    for (int isym = 0; isym < atom_selection_container.mol->get_cell_p()->GetNumberOfSymOps(); isym++) {
//       int err1 = atom_selection_container.mol->get_cell_p()->GetTMatrix(my_matt, isym, 0, 0, 0); 
//       int err2 = MMDBManager->get_cell_p()->GetTMatrix(my_matt, isym, 0, 0, 0); 
      
//       if (err1 != 0) {
// 	 cout << "!! something BAD with asc.mol mmdb::CMMDBCryst.GetTMatrix in convert_to_atoms_internal"
// 	      << endl;
//       } else {
// 	 cout << "DEBUG: symop " << isym << "       asc.mol seems OK..." << endl;
//       }
//       if (err2 != 0) {
// 	 cout << "!! something BAD with MMDB mmdb::CMMDBCryst.GetTMatrix in convert_to_atoms_internal"
// 	      << endl;
//       } else {
// 	 cout << "DEBUG: symop " << isym << " clipper::MMDB seems OK..." << endl;
//       }
//    }



   
   return atom_selection_container; 
}

// Kevin would use a template and in instatiation here, of course (and
// so should I :-)
//
atom_selection_container_t
BuildCas::convert_to_atoms(const clipper::Xmap<float> &l1,
			   const std::vector<coot::Cartesian> &c,
			   std::string molecule_name) const
{

   return convert_to_atoms_internal(l1.spacegroup(), l1.cell(), c, 1, molecule_name); 

}

atom_selection_container_t
BuildCas::convert_to_atoms(const clipper::Xmap<int> &l1,
			   const std::vector<coot::Cartesian> &c, 
			   std::string molecule_name) const
{

   return convert_to_atoms_internal(l1.spacegroup(), l1.cell(), c, 1, molecule_name); 

}

atom_selection_container_t
BuildCas::convert_to_atoms(const std::vector<coot::Cartesian> &c, std::string molecule_name) const { 

   return convert_to_atoms_internal(d_map_p->spacegroup(), d_map_p->cell(), c, 1, molecule_name);

} 


// we want a list of atoms that are "radius" Angstroems from the
// current_point. We want symmetry bones_points too, but there is no
// symmetry marking - we just want a list of atoms, well, points
// really (Cartesians, that is). 
// 
// This is really slow.  We should use the extents architecture to speed it up.
//
std::vector <coot::Cartesian>
BuildCas::point_list_by_symmetry(atom_selection_container_t AtomSel,
				 const std::vector<clipper::Coord_grid> &grids,
				 coot::Cartesian current_point, float radius,
				 short int use_grids) { // default arg

// note that we are no longer const because we construct big_ball_grid

   // std::vector <Cartesian_and_Grid> big_ball_l; // of points, _l local, no shadowing
   std::vector <coot::Cartesian> big_ball_l; // of points, _l local, no shadowing
   mmdb::mat44 my_matt;
   mmdb::Contact *contact;
   int ncontacts;
   
   if (AtomSel.n_selected_atoms > 0) { 
      
      mmdb::PAtom point_atom_p = new mmdb::Atom;
      point_atom_p->SetCoordinates(current_point.get_x(),
				   current_point.get_y(),
				   current_point.get_z(), 1.0, 99.9);
      
      // mmdb::CMMDBCryst *cryst_p =  (mmdb::CMMDBCryst *) &AtomSel.mol->get_cell();

      std::cout << "DEBUG: There are " << AtomSel.mol->GetNumberOfSymOps() << " sym ops" << std::endl; 
      std::cout << "symmetry expanding about " << current_point << std::endl; 

      for (int ix = -1; ix < 2; ix++) { 
	 for (int iy = -1; iy < 2; iy++) { 
	    for (int iz = -1; iz < 2; iz++) {

	       for (int isym = 0; isym < AtomSel.mol->GetNumberOfSymOps(); isym++) { 

		  int err = AtomSel.mol->GetTMatrix(my_matt, isym, ix, iy, iz); 

		  if (err != 0)
		     std::cout << "!! something BAD with mmdb::CMMDBCryst.GetTMatrix"
			  << std::endl;
		  
		  // cout << "Here: " << ix << " " << iy << " " << iz << " " << isym << endl; 

		  mmdb::PPAtom trans_selection = new mmdb::PAtom[AtomSel.n_selected_atoms];
		  for (int ii=0; ii<AtomSel.n_selected_atoms; ii++) {
		     
		     trans_selection[ii] = new mmdb::Atom;
		     trans_selection[ii]->Copy(AtomSel.atom_selection[ii]);
		     trans_selection[ii]->Transform(my_matt);
		  }
		  

		  // so now trans_selection are available for contact seeking
		  
		  contact = NULL;
		  AtomSel.mol->SeekContacts(point_atom_p, trans_selection,
					    AtomSel.n_selected_atoms,
					    0.0, radius,
					    0,  // seqDist
					    contact, ncontacts);
		  
		  if (ncontacts > 0 ) {

		     //cout << "got some contacts for  " 
		     //  << ix << " " << iy << " " << iz << " " << "symm: " << isym << endl;
		     
		     for (int ii=0; ii<ncontacts; ii++) {
			
			// consider pushing back contact[ii].id2, with 
			// std::vector <mmdb::PAtom> big_ball;
			// 
			coot::Cartesian a(trans_selection [ contact[ii].id2 ]->x,
				    trans_selection [ contact[ii].id2 ]->y,
				    trans_selection [ contact[ii].id2 ]->z); 

// 			coot::Cartesian_and_Grid cag(a, grids[ii]); 
// 			big_ball_l.push_back(cag);
			big_ball_l.push_back(a); 

			if (use_grids)
			   big_ball_grid.push_back(grids[ contact[ii].id2 ]); 
		     }
		     
		     // tidy up
		     delete [] contact;
		  }

		  // now tidy up trans_selection:
		  // 

		  for (int ii=0; ii<AtomSel.n_selected_atoms; ii++) {
		     delete trans_selection[ii]; 
		  } 
		  delete [] trans_selection; 
	       }
	    }
	 }
      }
   }
   return big_ball_l; 
}


// #include <algorithm> // needed?  (for std::map ?)
#include <map>
// see http://www.parashift.com/c++-faq-lite/containers-and-templates.html

// We passed a std::vector of bones points that are within a certain
// distance of centre_point. We want a set of std::vectors of clusters.
// How do we determine that a new position is (or is not) a member of
// an already-defined cluster?  We take the dot product of this std::vector
// (the difference of the bones point in question and centre_point)
// and members of other clusters (one at a time (obviously)) and find
// the cosine of the angle between them.  If the cosine of the angle
// (corresponding to some arbitary splitting angle (here we choose 20
// degrees)) is less than the cuttoff, then we have a new cluster.
// 
std::vector<std::vector<coot::Cartesian_and_Grid> > 
BuildCas::cluster_bones_points(std::vector<coot::Cartesian_and_Grid> bones_points, 
			       coot::Cartesian centre_point) const { 
  
   std::vector<cv_list_t> distinct_cv_list; // [cv: current_vector]. A
				       // cv_list_t is (geometrical)
				       // vector and atom index. 
   float cos_30 = 0.866; // as every schoolboy knows.
   int index;
   coot::Cartesian current_point;
   coot::Cartesian current_vector;
   short int not_distinct_flag;
   int matching_index; 

   // usage, cf: std::map<std::string, int, std::less<std::string> >  age;
   // 
   // 
   // Setup the associative array: The array index is atom index and
   // the referenced value (the thing in the pot) is the index of a
   // "matching" atom (i.e. an atom (bones point) that is not the far
   // (angular wise) from the indexing atom.
   // 
   // We use an associative array because there are (typically) 3 or 4
   // pots and if we indexed on the "matching" atom, we would need an
   // array of however many bones points we have... (10's of 1000s
   // perhaps).
   // 
   std::map<int, int> pot;

   std::cout << "cluster_bones_points will cluster " 
	<< bones_points.size() << " bones points" << std::endl; 

   std::vector<std::vector<coot::Cartesian_and_Grid> > cluster_vec;

   for (unsigned int ii=0; ii< bones_points.size(); ii++) { 
    
      // 
      current_point = bones_points[ii].pos;

      current_vector = current_point - centre_point; 

      not_distinct_flag = FALSE; // start by presuming a new distinct vector

      for (unsigned int i_d_cv=0; i_d_cv < distinct_cv_list.size(); i_d_cv++) {

	 // Find the dot product of this test vector (current_vector)
	 // and each distinct (geomtric) vector:
	 // 
	 float abcos_theta = coot::dot_product(distinct_cv_list[i_d_cv].geometrical_vector, 
					 current_vector); 

	 float a_length = distinct_cv_list[i_d_cv].geometrical_vector.amplitude();
	 float b_length =  current_vector.amplitude(); 

	 if ( a_length > 0.0 && b_length > 0.0) { 
	
	    float cos_theta = abcos_theta / (a_length * b_length); 
	
	    if (cos_theta > cos_30) {  // i.e. within the angle (seen similar before)
	       not_distinct_flag = TRUE; 
	       // note: atom_index indexes the bones_points vector.
	       matching_index = distinct_cv_list[i_d_cv].atom_index; 
	       break; 
	    }
	 } else { 
	    std::cout << "Error: bogus lengths in cluster_bones_points" << std::endl; 
	 } 
      }

      if (not_distinct_flag) { 

	 // this has already been seen (or something like it more likely)
	 // 
	 // so add this bones_point to the list of bones point of
	 // which distinct_cv_list[ii] was a member.
	 // (add it to the cluster_vec): 

	 // index is the cluster vector number and pot converts 
	 // from mmdb atom index to cluster vector number. 

	 index = pot[matching_index]; // associative array

	 cluster_vec[index].push_back(bones_points[ii]); 

      } else {

	 // it *was* distinct.  Add another element to (the thing that
	 // gets returned) cluster_vec:

	 // (but first) add a corresponding cv_list_t to distinct_cv_list: 
	 // 
	 distinct_cv_list.push_back(cv_list_t(current_vector, ii));

	 pot[ii] = cluster_vec.size(); 

	 // make room for the appropriate (new, different) bones_points:
	 // 
	 cluster_vec.resize(cluster_vec.size() + 1);
	 
	 // and add it also to the cluster vec:
	 //
	 // by filling the end slot: 
	 cluster_vec[cluster_vec.size() - 1].push_back(bones_points[ii]);

      }
   }
   
   // cluster_vec diagnostics
   for (unsigned int ii=0; ii < cluster_vec.size(); ii++) { 
      std::cout << "cluster " << ii << " has size " 
	   << cluster_vec[ii].size() << std::endl;
   }
   
   return cluster_vec; 
}


// cluster_centres.  Given a list of list of points that are close,
// return a list of central points, along with a corresponding grid. 
// 
std::vector<coot::Cartesian_and_Grid>
BuildCas::cluster_centres(std::vector<std::vector<coot::Cartesian_and_Grid> > cluster_vec) const { 

   std::vector<coot::Cartesian_and_Grid> cluster_centres; 

   // argh!  I want 'map'! (the scheme function). C++ hoops. Grrr.

   for (unsigned int ii=0; ii < cluster_vec.size(); ii++) { 

      // create a temporary vector<Cartesian> 
      std::vector<coot::Cartesian> v; 
      for (unsigned int j=0; j<cluster_vec[ii].size(); j++) {
	 v.push_back(cluster_vec[ii][j].pos);
      }

      coot::Cartesian_and_Grid av_pos(average_Cartesians(v), cluster_vec[ii][0].near_grid_point);
      
      cluster_centres.push_back(av_pos); 
   }
   return cluster_centres; 
}

#include "clipper/core/map_interp.h"

// pick the point from the branch points with the maximum density
//
score_and_cart
BuildCas::fit_first_in_segment(const clipper::Xmap<float> &map) {

   // 
   std::cout << "-----> starting initial fitting...." << std::endl; 
   std::cout << "searching " << branch_points_symm_expanded.size() 
	<< " coordinates (for branch points (symm expanded) " << std::endl; 

   score_and_cart sc = peak_search_simple(); 

   build.resize(i_current_build +1);
   build[i_current_build].resize(0); 
   build[i_current_build].push_back(sc); 
   n_fitted_in_current_segment++; 

   if ( ! (sc.score > 0.0)) { 
      // something bad.
      std::cout << "BADNESS in fit_first_in_segment" << std::endl; 
   } 

   std::cout << "Putting first atom at: " << sc.pos << std::endl; 
   
   std::cout << "------> done initial fitting...." << std::endl; 

   return sc; 
}

// return the maximum of the angstoms/grid-point  
float 
BuildCas::maximum_gridding(const clipper::Xmap<float> &map) const {
   
   float a_gridding = map.cell().a()/map.grid_sampling().nu(); 
   float b_gridding = map.cell().b()/map.grid_sampling().nv(); 
   float c_gridding = map.cell().c()/map.grid_sampling().nw(); 

   float max = 0; 
   
   if ( a_gridding > max) 
      max = a_gridding; 
   if ( b_gridding > max) 
      max = b_gridding; 
   if ( c_gridding > max) 
      max = c_gridding; 

   return max; 
} 

// the same except implicitly use the segement_map.
// 
float 
BuildCas::maximum_gridding() const {
   
   float a_gridding = segment_map.cell().a()/segment_map.grid_sampling().nu(); 
   float b_gridding = segment_map.cell().b()/segment_map.grid_sampling().nv(); 
   float c_gridding = segment_map.cell().c()/segment_map.grid_sampling().nw(); 

   float max = 0; 
   
   if ( a_gridding > max) 
      max = a_gridding; 
   if ( b_gridding > max) 
      max = b_gridding; 
   if ( c_gridding > max) 
      max = c_gridding; 

   return max; 
} 



score_and_cart
BuildCas::peak_search_simple() const { 

   score_and_cart score_and_cartesian; 

   coot::Cartesian current_max;
   float current_max_density_val = 0.0; 
   float dv; 

   
   for(unsigned int ii=0; ii<branch_points_symm_expanded.size(); ii++) {

      float score = 
	 prebuilt_exclusion_score(branch_points_symm_expanded[ii]); 
      if (score > 0) {

	 // we need to find the density at the coordinates
	 // 
	 clipper::Coord_orth c_o(branch_points_symm_expanded[ii].x(), 
				 branch_points_symm_expanded[ii].y(), 
				 branch_points_symm_expanded[ii].z()); 

	 clipper::Coord_frac c_f = c_o.coord_frac(d_map_p->cell()); 

	 clipper::Interp_linear::interp(*d_map_p, 
					c_f.coord_map(d_map_p->grid_sampling()), dv); 
      

	 if (dv > current_max_density_val) { 
	 
	    current_max_density_val = dv; 
	    current_max = branch_points_symm_expanded[ii];
	 }
      }
   }

   score_and_cartesian.score = current_max_density_val; 
   score_and_cartesian.pos = current_max; 

   return score_and_cartesian; 

} 

// Fitting the second point, we only have a distance constraint.  
// 
// So, we make a grid, with point at the centre and find the maxium 
// of some function f(density_value, deviation_from_ideal_length_score). 
// 
// Remember to include distance to nearest branch point at some stage.
// 
// Also include extra score for being on the same segment (use segment_map). 
//
// Also include extra score for not being too close to atoms already built.
// 
score_and_cart
BuildCas::peak_search_distance(coot::Cartesian previous_atom, 
			       coot::Cartesian point) const {

   float grid_size = 0.75; 
   float grid_step = 0.25; 

   // 
   float best_score = 0,  score; 
   // float d; // distance (of point to trial Ca position); 

   clipper::Coord_orth c_o_point(point.x(), point.y(), point.z());
   clipper::Coord_frac c_f_point = c_o_point.coord_frac(d_map_p->cell()); 
   clipper::Coord_grid c_g_point = c_f_point.coord_grid(d_map_p->grid_sampling()); 

   clipper::Coord_orth c_o_prev(previous_atom.x(), 
				previous_atom.y(), 
				previous_atom.z());
   clipper::Coord_frac c_f_prev = c_o_prev.coord_frac(d_map_p->cell()); 
   clipper::Coord_grid c_g_prev = c_f_prev.coord_grid(d_map_p->grid_sampling()); 

   
   float segment_score_val = segment_score(c_g_point, c_g_prev); 

   float mpds = mid_points_density_score(previous_atom, point); 


   score_and_cart score_and_cartesian; 

   // we only need to calculate prebuilt_exclusion_score for the centre point
   // 
   float prebuilt_exclusion_score_val = prebuilt_exclusion_score(point);

   float save_deviation_from_ideal_length_score = 0.0; 
   float save_branch_point_proximity_score = 0.0; 

   float save_prebuilt_exclusion_score = 0.0;
   float save_segment_score = 0.0;

   // float sc;
   scores myscores; 
      
   // 
   for ( float xa = -grid_size; xa <= grid_size; xa += grid_step) { 
      for ( float ya = -grid_size; ya <= grid_size; ya += grid_step) { 
	 for ( float za = -grid_size; za <= grid_size; za += grid_step) { 

	    coot::Cartesian trial_point(point.x()+xa, point.y()+ya, point.z()+za); 

	    // consider removing d_map_p from the call.
	    myscores = non_angle_micro_point_score(previous_atom, trial_point); 

	    score = myscores.score*segment_score_val*prebuilt_exclusion_score_val*mpds; 

	    if (score > best_score) { 

	       best_score = score; 
	       
	       score_and_cartesian.score = score; 
	       score_and_cartesian.pos = trial_point; 
	       score_and_cartesian.near_grid_point = c_g_point; 

	       // debugging
	       save_deviation_from_ideal_length_score = 
	         myscores.deviation_from_ideal_length_score_val; 
	       save_branch_point_proximity_score =
	         myscores.branch_point_proximity_score_val;

	       // invariant of xa, ya, za: 
	       // 
	       save_prebuilt_exclusion_score = prebuilt_exclusion_score_val; 
	       save_segment_score = segment_score_val; 
	    }
	 }
      }
   }

   // debugging
   std::cout << "      deviation_from_ideal_length_score(d)      " 
	<< save_deviation_from_ideal_length_score << std::endl; 
   std::cout << "      branch_point_proximity_score(trial_point) " 
	<< save_branch_point_proximity_score<< std::endl; 
   std::cout << "      prebuilt_exclusion_score_val              " 
	<< save_prebuilt_exclusion_score << std::endl; 
   std::cout << "      segment_score_val                         " 
	<< save_segment_score << std::endl; 
   std::cout << "      midpoint density score                    " 
	<< mpds << std::endl; 
   std::cout << "      TOTAL SCORE    ----->                     " 
	<< score_and_cartesian.score << std::endl; 

   std::cout << "peak_search_distance filled near_grid_point: " 
	<< score_and_cartesian.near_grid_point.format() << std::endl; 

   return score_and_cartesian; 
}

// c.f. peak_search_distance()
// 
// Note it is quite possible that the return .score is negative (since there
// was no plausible place to put a grid checked point
// 
score_and_cart
BuildCas::peak_search_distance_theta_2(const TreeNode *node) //trial_centre_point
				       const
{

   // return this:
   score_and_cart score_and_cartesian; 
   score_and_cartesian.score = 0;

   // protection
   if (! node) return score_and_cartesian;
   if (! node->parent) return score_and_cartesian;
   if (! node->parent->parent) return score_and_cartesian;

   // recover the old variables:
   coot::Cartesian point =                      node->pos; 
   coot::Cartesian ith_point =  node->parent->parent->pos; 
   coot::Cartesian ith_plus_one_point = node->parent->pos;

   // to match with peak_search_distance(): 
   // 
   coot::Cartesian previous_atom = ith_plus_one_point; 

   float grid_size = 0.75;
   float grid_step = 0.25;

   // we only need to calculate prebuilt_exclusion_score for the centre point
   // 

   // int t0 = glutGet(GLUT_ELAPSED_TIME); 
   float prebuilt_exclusion_score_val = prebuilt_exclusion_score(point);
   // int t1 = glutGet(GLUT_ELAPSED_TIME); 
   // cout << "INFO: prebuilt_exclusion_score: elapsed time: " << t1-t0 << "ms" << endl; 

   if (prebuilt_exclusion_score_val > 0) { 

      // 
      float best_score = 0,  score; 
      // float d; // distance (of point to trial Ca position); 

//       clipper::Coord_orth c_o_point(point.x(), point.y(), point.z());
//       clipper::Coord_frac c_f_point = c_o_point.coord_frac(d_map_p->cell()); 
//       clipper::Coord_grid c_g_point = c_f_point.coord_grid(d_map_p->grid_sampling()); 

      clipper::Coord_grid c_g_point = node->near_grid_point; 

//       clipper::Coord_orth c_o_prev(previous_atom.x(), 
// 				   previous_atom.y(), 
// 				   previous_atom.z());
//       clipper::Coord_frac c_f_prev = c_o_prev.coord_frac(d_map_p->cell()); 
//       clipper::Coord_grid c_g_prev = c_f_prev.coord_grid(d_map_p->grid_sampling()); 

      clipper::Coord_grid c_g_prev = node->parent->near_grid_point; 

      // set to some obviously bogus values, so we can see when the don't get 
      // overwritten: 
      float save_deviation_from_ideal_length_score = -999.9; 
      float save_branch_point_proximity_score      = -999.9; 

      float save_prebuilt_exclusion_score = -999.9; 
      float save_segment_score            = -999.9; 
      float save_theta_2_score            = -999.9; 
      float save_m_p_score_score          = -999.9; 
      float save_density_score            = -999.9; 

      float segment_score_val = segment_score(c_g_point, c_g_prev); 

      // mid point of the ca pseudo bond
      // (This is more or less the same for each micropoint, so lets only 
      // calculate it for each trial (cluster centre) point
      //
      float mpds = mid_points_density_score(ith_plus_one_point, point); 

      // float sc;
      float theta_2_sc; 
      scores m_p_scores; 
      
      // 
      // int t0 = glutGet(GLUT_ELAPSED_TIME); 

      for ( float xa = -grid_size; xa <= grid_size; xa += grid_step) { 
	 for ( float ya = -grid_size; ya <= grid_size; ya += grid_step) { 
	    for ( float za = -grid_size; za <= grid_size; za += grid_step) { 

	       coot::Cartesian trial_point(point.x()+xa, point.y()+ya, point.z()+za); 

	       // remove map from the call


	       m_p_scores = non_angle_micro_point_score(previous_atom, 
							trial_point); 
	       theta_2_sc = theta_2_score(ith_point, 
					  ith_plus_one_point, 
					  trial_point);
	       score = m_p_scores.score
		  * segment_score_val
		  * prebuilt_exclusion_score_val
		  * theta_2_sc
		  * mpds; 
	    
	       if (score > best_score) { 

		  best_score = score; 
	       
		  score_and_cartesian.score = score; 
		  score_and_cartesian.pos = trial_point;
		  score_and_cartesian.near_grid_point = c_g_point; 

		  // debug:
		  // The need for this went away when I stopped getting
		  // negative scores (for theta_2)
		  // 
		  //cout << "peak_search_distance_theta_2:" << endl; 
		  //cout << "updating returnable pos to " 
		  //    << score_and_cartesian.pos << " from: " << endl
		  //	    << "trial_point: " << trial_point << endl
		  //	    << "      point: " << point << endl;
	      
		  // debugging
		  save_deviation_from_ideal_length_score = 
		     m_p_scores.deviation_from_ideal_length_score_val; 
		  save_branch_point_proximity_score =
		     m_p_scores.branch_point_proximity_score_val;
		  save_theta_2_score = theta_2_sc; 
		  save_m_p_score_score = m_p_scores.score; 
		  save_density_score   = m_p_scores.density_score_val; 

		  // invariant of xa, ya, za: 
		  // 
		  save_prebuilt_exclusion_score = prebuilt_exclusion_score_val; 
		  save_segment_score = segment_score_val;  


	       }
	    }
	 }
      }

      // int t1 = glutGet(GLUT_ELAPSED_TIME); 
      // cout << "peak_search_distance_theta_2 inner time: " << t1-t0 << "ms" << endl; 



      // debugging
      std::cout << "      micro_point_scores.score                  " 
	   << save_m_p_score_score << std::endl; 
      std::cout << "      density_score                             " 
	   << save_density_score << std::endl; 
      std::cout << "      deviation_from_ideal_length_score(d)      " 
	   << save_deviation_from_ideal_length_score << std::endl; 
      std::cout << "      branch_point_proximity_score(trial_point) " 
	   << save_branch_point_proximity_score<< std::endl; 
      std::cout << "      prebuilt_exclusion_score_val              " 
	   << save_prebuilt_exclusion_score << std::endl; 
      std::cout << "      segment_score_val                         " 
	   << save_segment_score << std::endl; 
      std::cout << "      theta_2_score_val                         " 
	   << save_theta_2_score << std::endl; 
      std::cout << "      midpoint density score                    " 
	   << mpds << std::endl; 
      std::cout << "      TOTAL SCORE    ----->                     " 
	   << best_score << std::endl; 

   }
      
   return score_and_cartesian; 
	    
}

// c.f. peak_search_distance()
// 
// Note it is quite possible that the return .score is negative (since there
// was no plausible place to put a grid checked point
// 
score_and_cart
BuildCas::peak_search_distance_angle_torsion(const TreeNode *node) //trial_centre_point
				             const 
{

   // return this:
   score_and_cart score_and_cartesian; 
   score_and_cartesian.score = 0; 

   // recover the old variables:
   coot::Cartesian point =                      node->pos; 
   // coot::Cartesian ith_point =  node->parent->parent->pos;  // not used
   coot::Cartesian ith_plus_one_point = node->parent->pos;

   // to match with peak_search_distance(): 
   // 
   coot::Cartesian previous_atom = ith_plus_one_point; 

//    float grid_size = 0.75; // FIXME
//    float grid_step = 0.25; 

   float grid_size = 0.75; // FIXME
   float grid_step = 0.25; 

   // we only need to calculate prebuilt_exclusion_score for the centre point
   // 

   // int t0 = glutGet(GLUT_ELAPSED_TIME); 
   float prebuilt_exclusion_score_val = prebuilt_exclusion_score(point);
   // int t1 = glutGet(GLUT_ELAPSED_TIME); 
   // cout << "INFO: prebuilt_exclusion_score: elapsed time: " << t1-t0 << "ms" << endl; 

   if (prebuilt_exclusion_score_val > 0) { 

      // float dv; // density value

      // 
      float best_score = 0,  score; 
      // float d; // distance (of point to trial Ca position); 

//       clipper::Coord_orth c_o_point(point.x(), point.y(), point.z());
//       clipper::Coord_frac c_f_point = c_o_point.coord_frac(d_map_p->cell()); 
//       clipper::Coord_grid c_g_point = c_f_point.coord_grid(d_map_p->grid_sampling()); 

      clipper::Coord_grid c_g_point = node->near_grid_point;

//       clipper::Coord_orth c_o_prev(previous_atom.x(), 
// 				   previous_atom.y(), 
// 				   previous_atom.z());
//       clipper::Coord_frac c_f_prev = c_o_prev.coord_frac(d_map_p->cell()); 
//       clipper::Coord_grid c_g_prev = c_f_prev.coord_grid(d_map_p->grid_sampling()); 

      clipper::Coord_grid c_g_prev = node->parent->near_grid_point; 

      // set to some obviously bogus values, so we can see when the don't get 
      // overwritten: 
      float save_deviation_from_ideal_length_score = -999.9; 
      float save_branch_point_proximity_score      = -999.9; 

      float save_prebuilt_exclusion_score = -999.9; 
      float save_segment_score            = -999.9; 
      float save_angle_torsion_score      = -999.9; 
      float save_m_p_score_score          = -999.9; 
      float save_density_score            = -999.9; 

      float segment_score_val = segment_score(c_g_point, c_g_prev); 

      // mid point of the ca pseudo bond
      // (This is more or less the same for each micropoint, so lets only 
      // calculate it for each trial (cluster centre) point
      //
      float mpds = mid_points_density_score(ith_plus_one_point, point); 

      // float sc;
      float angle_torsion_sc; 
      scores m_p_scores; 
      TreeNode trial_node; 
      // 
      // int t0 = glutGet(GLUT_ELAPSED_TIME); 

      for ( float xa = -grid_size; xa <= grid_size; xa += grid_step) { 
	 for ( float ya = -grid_size; ya <= grid_size; ya += grid_step) { 
	    for ( float za = -grid_size; za <= grid_size; za += grid_step) { 

	       coot::Cartesian trial_point(point.x()+xa, point.y()+ya, point.z()+za); 

	       // remove map from the call


	       m_p_scores = non_angle_micro_point_score(previous_atom, 
							trial_point); 

	       trial_node.setup(node->parent, trial_point); 
	       angle_torsion_sc = angle_torsion_score(&trial_node); 
							 
	       score = m_p_scores.score
		  * segment_score_val
		  * prebuilt_exclusion_score_val
		  * angle_torsion_sc
		  * mpds; 
	    
	       if (score > best_score) { 

		  best_score = score; 
	       
		  score_and_cartesian.score = score; 
		  score_and_cartesian.pos = trial_point;
		  score_and_cartesian.near_grid_point = c_g_point; 

		  // debugging
		  save_deviation_from_ideal_length_score = 
		     m_p_scores.deviation_from_ideal_length_score_val; 
		  save_branch_point_proximity_score =
		     m_p_scores.branch_point_proximity_score_val;
		  save_angle_torsion_score = angle_torsion_sc; 
		  save_m_p_score_score = m_p_scores.score; 
		  save_density_score   = m_p_scores.density_score_val; 

		  // invariant of xa, ya, za: 
		  // 
		  save_prebuilt_exclusion_score = prebuilt_exclusion_score_val; 
		  save_segment_score = segment_score_val;  


	       }
	    }
	 }
      }


      // debugging
      std::cout << "      m_p_scores.score                          " 
	   << save_m_p_score_score << std::endl; 
      std::cout << "      density_score                             " 
	   << save_density_score << std::endl; 
      std::cout << "      deviation_from_ideal_length_score(d)      " 
	   << save_deviation_from_ideal_length_score << std::endl; 
      std::cout << "      branch_point_proximity_score(trial_point) " 
	   << save_branch_point_proximity_score<< std::endl; 
      std::cout << "      prebuilt_exclusion_score_val              " 
	   << save_prebuilt_exclusion_score << std::endl; 
      std::cout << "      segment_score_val                         " 
	   << save_segment_score << std::endl; 
      std::cout << "      angle_torsion_score_val                   " 
	   << save_angle_torsion_score << std::endl; 
      std::cout << "      midpoint density score                    " 
	   << mpds << std::endl; 
      std::cout << "      TOTAL SCORE    ----->                     " 
	   << best_score << std::endl; 

   }

   if (score_and_cartesian.score == 0.0) { 
      std::cout << "returning ZERO score and garbage near_grid_point" << std::endl; 
   } 
   return score_and_cartesian; 
	    
}

// Change the indexing of the atoms to match make-angles.awk
// this is theta_2 and we are trialling positions of i+3. 
// 
// This code took some thinking about!
// 
float 
BuildCas::theta_2_score(coot::Cartesian ith_plus_1_point, 
			coot::Cartesian ith_plus_2_point, 
			coot::Cartesian ith_plus_3_point) const { 

   // note: the construction of these vectors must be the same
   // as the one that we pulled out the angles for the database
   // i.e. make-angles.awk
   // 
   // 
   coot::Cartesian diff2 = ith_plus_1_point - ith_plus_2_point; 
   coot::Cartesian diff3 = ith_plus_3_point - ith_plus_2_point; 

   float cos_theta_2 = 
      coot::dot_product(diff2,diff3)/(diff2.amplitude()*diff3.amplitude()); 

   float theta_2 = acos(cos_theta_2)*RADTODEG; 

      
   float theta_2_sc = angle_info.theta_2_score(theta_2); 

//    cout << "        theta_2 (in degs): " << theta_2 << " score: " 
// 	<< theta_2_sc << endl; 

   return theta_2_sc*0.001; 

}


// we will make up a Treenode tree from a previously read in pdb file
// 
// and add debugging code to angle_torsion_score so that it writes out
// the angles and torsions.  We will use a helical fragment from fffear.
// 
void 
BuildCas::check_angle_torsion(atom_selection_container_t asc) const {

   TreeNode *node = NULL; 

   // yep, it's the same backwards...
   // for (int i=asc.n_selected_atoms-1; i>=0; i--) { 

   for (int i=0; i< asc.n_selected_atoms; i++) { 

      if (std::string(asc.atom_selection[i]->name) == " CA " ) { 

	 coot::Cartesian pos(asc.atom_selection[i]->x,
			     asc.atom_selection[i]->y,
			     asc.atom_selection[i]->z); 

	 std::cout << "Got a CA at " << pos << std::endl; 

	 TreeNode *new_node = new TreeNode; 
	 new_node->setup(node, pos);

	 node = new_node;  // setup for next round
      }
   }

   // node is now the top of the tree
   // 

   float ats; 

   while (1) { 

      if (node != NULL) { 
	 if (node->parent != NULL) { 
	    if (node->parent->parent != NULL) { 
	       if (node->parent->parent->parent != NULL) { 
		  
		  ats = angle_torsion_score(node); 

		  std::cout << "angle_torsion_score: " << ats << std::endl; 

		  node = (TreeNode *) node->parent;  // cast away constness.

	       } else { 
		  break;
	       } 
	    } else { 
	       break;
	    } 
	 } else { 
	    break;
	 } 
      } else { 
	 break;
      }
   }
   delete node;
} 

// much like the above:
float 
BuildCas::angle_torsion_score(const TreeNode *node) const { 

   //                        / i+3      // we are scoring the position of i+3
   //                       /           // (which is node)
   //                      /
   //     i+1             /
   //     ---------------- i+2
   //    /
   //   /
   //  /
   //  i
   //

   coot::Cartesian ith_plus_3_point = node->pos; 
   coot::Cartesian ith_plus_2_point = node->parent->pos;
   coot::Cartesian ith_plus_1_point = node->parent->parent->pos;
   coot::Cartesian ith_point        = node->parent->parent->parent->pos; 

   coot::Cartesian v1 = ith_plus_1_point - ith_point; 
   coot::Cartesian v2 = ith_plus_2_point - ith_plus_1_point; 
   coot::Cartesian v3 = ith_plus_3_point - ith_plus_2_point; 
   coot::Cartesian v3b= ith_plus_2_point - ith_plus_3_point; // backwards

   coot::Cartesian v_prod_1 = coot::cross_product(v1, v2); 
   coot::Cartesian v_prod_2 = coot::cross_product(v2, v3); 

   float dot_prod  = coot::dot_product(v_prod_1, v_prod_2); 
   
   float cos_tor = dot_prod/(v_prod_1.amplitude()*v_prod_2.amplitude()); 

   // sign
   coot::Cartesian cross_1_2 = coot::cross_product(v_prod_1,v_prod_2); 
   
   float sign = coot::dot_product(cross_1_2, v2); 

   float tor = sign < 0 ? acos(cos_tor)*RADTODEG : -acos(cos_tor)*RADTODEG; 

   float cos_theta_2 = coot::dot_product(v2, v3b)/(v2.amplitude()*v3.amplitude()); 
   float theta_2 = acos(cos_theta_2)*RADTODEG; 

   float table_val = angle_info.prob_angle_torsion(theta_2, tor); 

   // debug
//     cout << "in angle_torsion_score: angle=" << theta_2 
//  	<< " and torsion=" << tor << ", "; 

   // the results we get from this are: 
   // for theoretical helix: 
   // 
   // angle=89.4643 and torsion=-50.4069, angle_torsion_score: 11000
   
   return table_val; 

} 



// The bigger the length deviation (length), the lower the score.
// Good score: 2.5: bad score 0.5
float 
BuildCas::deviation_from_ideal_length_score(float length) const { 

   return 1/(grid_dependent_distance_param + fabs(length - 3.8)); 
} 


// segment_score returns small if the closest skeleton point to the previous one 
// is on a different segment and big if it is on the same segment.
// 
float 
BuildCas::segment_score(const clipper::Coord_grid &c_g_point, 
			const clipper::Coord_grid &c_g_previous_atom) const { 


   int i_previous_atom_segment = segment_map.get_data(c_g_previous_atom); 

   // The segment value of the (centre) trial point: 
   // 

   int i_trial_segment = segment_map.get_data(c_g_point); 

   float score; 

   if (i_trial_segment  == 0) { 
      if (i_previous_atom_segment == 0 ) {
	 score = 0.101;  // ooops.  Why were we not on a skeleton point?
      } else { 
	 score = 0.102; // previous atom OK, but target is not
      } 

   } else { 

      // trial point was on a skeleton:
      
      if (i_previous_atom_segment == 0 ) {
	 score = 0.103; // but previous was not
      } else { 
	 if (i_previous_atom_segment == i_trial_segment) { 
	    
	    std::cout << "depth search testing " << c_g_point.format() 
		 << " and " << c_g_previous_atom.format() << std::endl; 

	    if (treenodemap.get_data(c_g_previous_atom).neighbs.size() == 0) {
	       std::cout << "woops! no neighbours for depth search start " << c_g_previous_atom.format() << std::endl; 
	       if (treenodemap.get_data(treenodemap.get_data(c_g_previous_atom).near_grid_point).neighbs.size() == 0) {
		  std::cout << "woops! but constistantly bad in depth search start " 
		       << treenodemap.get_data(c_g_previous_atom).near_grid_point.format() << std::endl; 
	       } else { 
		  std::cout << "whoooo! inconsistantly bad!  in depth search start" 
		       << treenodemap.get_data(c_g_previous_atom).near_grid_point.format() << std::endl; 
	       } 
	    }
	    
	    if (depth_search_skeleton(c_g_previous_atom, c_g_point)) {
	       score = 100;
	    } else { 
		
	       // was in the same skeleton, but not reachable:
	       score = 10; 
	    }
	 } else { 
	    
	    // we were for sure in defined different segments

	    std::cout << "certain different segment: " << i_previous_atom_segment
		 << " and " << i_trial_segment << std::endl;
	    
	    score = 10; 
	 } 
      }
   }

   return score; 

}


void 
BuildCas::transfer_segment_map(clipper::Xmap<int> *skel_p) const { 

   
   clipper::Xmap_base::Map_reference_index ix;

   for (ix = skel_p->first(); !ix.last(); ix.next() ) {
      (*skel_p)[ix] = segment_map[ix]; 
   }
} 

// debugging function
void
BuildCas::show_segment_map() { 

   std::vector<coot::Cartesian> a; 

} 


score_and_cart
BuildCas::fit_forth_in_segment(const clipper::Xmap<float> &map) {
   
   score_and_cart score_and_cartesian; 

   return score_and_cartesian; 

}

score_and_cart
BuildCas::fit_next_in_segment(const clipper::Xmap<float> &map) {

   score_and_cart score_and_cartesian; 

   return score_and_cartesian; 

}


// Given a vector of cartesians, return a vector of cartesians that
// are "close to" a target point (start_point).
//
// If were were using mmdb (i.e. mmdb::PAtoms, this function may be faster -
// consider recoding if this function is slow).
//
// We naively check all points in the big_ball.
// 
std::vector<coot::Cartesian_and_Grid>
BuildCas::select_by_distance(coot::Cartesian start_point,
			     float near_point, float far_point) const { // 3.7A +/- a bit

   std::vector<coot::Cartesian_and_Grid> result;
   coot::Cartesian dist_vec; 

   for(unsigned int ii=0; ii<big_ball.size(); ii++) {

      dist_vec = start_point - big_ball[ii];

      if (dist_vec.amplitude() < far_point) { 
	 if (dist_vec.amplitude() > near_point) {

	    coot::Cartesian_and_Grid cag(big_ball[ii],big_ball_grid[ii]); 

	    result.push_back(cag);
	 }
      }
   }
   return result; 
}

// We are passed a map and an atom_selection_container asc which contains the
// bones atoms in the ASU - all of them, not just the branch points
//
// We return big_ball, a set of coordinates that are (typically)
// within a 100A sphere of the origin.
//
// Note to self, we don't need to return a value here.  We can use
// void and avoid the convert_to_atoms here (historically, we returned
// a vaule so that I could display them on the graphics - and we can
// use asc_to_graphics here if we want to do that).
// 
atom_selection_container_t 
BuildCas::build_big_ball(const clipper::Xmap<float> &map, 
			 atom_selection_container_t asc, 
			 const std::vector<clipper::Coord_grid> &grids) {

   // FIXME: use instead the molecule centre finding routine used in ncs-search.
   //        (uses Bond_lines_ext)
   // 
   // coot::Cartesian centre_point(30.0,30.0,30.0); // calculate this from a ncs-search.
   coot::Cartesian centre_point(0.0,0.0,0.0); // calculate this from a ncs-search.

   expansion_centre = centre_point;  // set the class variable
   expansion_centre_is_set = 1;      // it is set. 

   big_ball = point_list_by_symmetry(asc, grids, centre_point, 50.0); 

   // We want to create big_ball_grid, which corresponds to big_ball, 
   // i.e, big_ball[i] is a Cartestian and big_ball_grid[i] is the grid 
   // point that corresponds to that Cartestian. 
   // 
   // Actually, that is done in point_list_by_symmetry(). 
   // 
   // now big_ball_grid has been filled.
   

   if (big_ball.size() == 0 ) { 
      // strange
      std::cout << "There are no (symmetry expanded) skeleton (not just branch) points" << std::endl;
   } 
   // 
   std::cout << "DEBUG: build_big_ball call convert_to_atoms..." << std::endl; 
   atom_selection_container_t asc_out = convert_to_atoms(map, big_ball, "BIG BALL"); 
   std::cout << "DEBUG: build_big_ball call convert_to_atoms done!" << std::endl; 
   
   return asc_out; 
}

// So we have a point somewhere is space (middle_mol).  The idea is to
// apply symmetry and unit cell translations to it until it is as
// close to the target_point as posible and then return that point.
// 
coot::Cartesian
BuildCas::move_by_symmetry(coot::Cartesian middle_mol, 
			   coot::Cartesian target_point,
			   mmdb::Cryst *cryst_p) const {

   float current_min_dist = (middle_mol - target_point).amplitude(); 
   float test_dist; 
   mmdb::mat44 my_matt;
   mmdb::PAtom atom = new mmdb::Atom;
   mmdb::PAtom trans_atom = new mmdb::Atom; 

   // we will use mmdb atoms to have the symmetry applied to it.
   // 
   atom->SetCoordinates(middle_mol.x(), middle_mol.y(), middle_mol.z(), 1.0, 99);

   // 
   std::cout << "atom from middle_mol: " << atom->x << " " 
	<< atom->y << " " << atom->z << std::endl;

   short int moved_it_flag = TRUE; 

   while (moved_it_flag) { 

      moved_it_flag = FALSE; 

      for (int ix = -1; ix < 2; ix++) { 
	 for (int iy = -1; iy < 2; iy++) { 
	    for (int iz = -1; iz < 2; iz++) {
	    
	       for (int isym = 0; isym < cryst_p->GetNumberOfSymOps(); isym++) { 
	       
		  cryst_p->GetTMatrix(my_matt, isym, ix, iy, iz);

		  trans_atom->Copy(atom);
		  trans_atom->Transform(my_matt);

		  coot::Cartesian cart_atom(trans_atom->x, 
				      trans_atom->y, 
				      trans_atom->z); 

		  std::cout << "testing atom at: " << cart_atom << std::endl; 

		  test_dist = (cart_atom - target_point).amplitude(); 

		  if (test_dist < current_min_dist) { 
		  
		     atom->Copy(trans_atom);
		     current_min_dist = test_dist; 
		     moved_it_flag = TRUE;

		  }
	       }
	    }
	 }
      }
   }

   coot::Cartesian a(atom->x, atom->y, atom->z); 
   delete atom;  
   delete trans_atom;
   
   return a; 

}


#include "../src/xmap-utils.h" // temporary
// convert to atoms
// 
// atom_selection_container_t 
asc_and_grids
BuildCas::all_skel_pts_in_asu(const clipper::Xmap<float> &map,
			      const clipper::Xmap<int>   &l1,
			      float cut_off) const
{
   //
   std::vector<coot::Cartesian> vc;
   std::vector<clipper::Coord_grid> vec_grids; 
   coot::Cartesian t; 

   // cout << "in all_skel_pts_in_asu map_density_distribution of l1:" << endl; 
   // map_density_distribution(l1, 1); 


   clipper::Xmap_base::Map_reference_index ix;

   for (ix = l1.first(); !ix.last(); ix.next() ) {

      if (l1[ix] > 0) {

	 if (map[ix] > cut_off) { 

	    clipper::Coord_orth c_o =
	       ix.coord().coord_frac(l1.grid_sampling()).coord_orth(l1.cell());
	    
	    vc.push_back(coot::Cartesian(c_o.x(), c_o.y(), c_o.z()));
	    vec_grids.push_back(ix.coord()); 
	 }
      }
   }

   std::cout << "all_skel_pts_in_asu calls convert_to_atoms...." << std::endl; 

   asc_and_grids asg(convert_to_atoms(map, vc, "ASU SKEL PTS"), vec_grids);
   return asg; 
}

asc_and_grids
BuildCas::toplevel_skel_pts_in_asu() const { 

   return all_skel_pts_in_asu(*d_map_p, segment_map, map_cut_off); 

} 


// Manipulate skel so that each strand/segment (set of connected
// skeleton points) is marked with a different integer.
//
// I believe that we move over the ASU boundaries with impunity,
// clipper Does The Right Thing (unlike Quanta). 
// 
// 
// Return the number of segments. 
// 
// Note to self: Show the segments in graphics, each segment a
// different colour (so that we know that we have done it right).
//
// Another note to self: I am worrying about why in the test, atom 3
// or so is built onto another segment despite this being massively
// downweighted in comparison to a same segment target.  I realise
// that segment map contains all the bone side chains from the raw
// clipper skeleton and that is the map that is used in segment_score
// (which also does the skeleton depth search (to see if the target is
// reachable in 10 or so grid steps)).
// 
// So, we first want to get rid of side chains in the segment map and
// also use the map cutoff which until now (21Oct2002) was not being
// used.
// 
int 
BuildCas::count_and_mark_segments(const clipper::Xmap<int>   &skel,
				  const clipper::Xmap<float> &map,
				  float cut_off) { 


   clipper::Xmap_base::Map_reference_index ix;

   // redimension segment_map
   // 
   segment_map.init(d_map_p->spacegroup(), 
		d_map_p->cell(), 
		d_map_p->grid_sampling()); 
   segment_map_filled = 1; // filled

   int toplevel = 0; 
   clipper::Xmap_base::Map_reference_index iy;
   for ( iy = skel.first(); !iy.last(); iy.next() )
      toplevel = clipper::Util::max( toplevel, skel[iy] );

   int n_toplevel = 0; 
   for ( iy = skel.first(); !iy.last(); iy.next() ) { 
      if (skel[iy] == toplevel) 
	 n_toplevel++; 
   }
   std::cout << "DEBUG: There were " << n_toplevel << " toplevel (" 
	<< toplevel << ") skeleton points"  << std::endl; 

   // First mark everything that is non-zero as -1 we will re-invert
   // later on in this routine.  (recall that we may have marked the
   // skeleton according to depth (i.e. tips have lower score than
   // mainchain core).
   // 
   // zero it:
   // 
   for (ix = segment_map.first(); !ix.last(); ix.next() )
      segment_map[ix] = 0;
   // 
   // add values:
   //
   for (ix = segment_map.first(); !ix.last(); ix.next() )
      if ( skel[ix] == toplevel)
	 segment_map[ix] = - toplevel; 

   int i_segment_number = 0; 
   clipper::Skeleton_basic::Neighbours neighb( segment_map );

   for (ix = segment_map.first(); !ix.last(); ix.next() ) {	

      if ( segment_map[ix] < 0 ) { // unmarked, so far
	 
	 i_segment_number++;
	 trace_along(ix.coord(), neighb, i_segment_number, 
		     toplevel, cut_off); // recursive function
	 
	 // debugging 
	 if ( ! (segment_map[ix] == i_segment_number) ) { 
	    std::cout << "ERROR ERROR ERROR ERROR" << std::endl; 
	    std::cout << "ERROR ERROR ERROR ERROR: in segment_map segment number" << std::endl; 
	 } 
      }
   }
   return i_segment_number; 
}


// FIXME:  is there a way not to use set_data? and instead us [] indexing.
void
BuildCas::trace_along(const clipper::Coord_grid &c_g_start, 
		      const clipper::Skeleton_basic::Neighbours &neighb,
		      int i_segment_number,
		      int i_max_level, 
		      float cut_off) {

   // question: does neighbs contain the central point?
   // answer: no.
   // 
   segment_map.set_data(c_g_start, i_segment_number); 
   clipper::Coord_grid c_g; 

   for(int i=0; i< neighb.size(); i++) {

      c_g = c_g_start + neighb[i]; 

      if (segment_map.get_data(c_g) == - i_max_level ) { 

	 // if (d_map_p->get_data(c_g) > cut_off) { 

	 segment_map.set_data(c_g, i_segment_number);  // ? necessary?
	    
	 trace_along (c_g, neighb, i_segment_number, i_max_level, cut_off); 
	    // }
      }
   }

}

// theta_2 is the second opening angle in a torsion
// distance is the length between Atom3 and Atom4. 
// 
// We return the position of atom4.
// 
// We pass theta_2 and torsion in radians. 
// 
coot::Cartesian
BuildCas::position_by_torsion(float theta_2, float torsion, float dist) const { 

   coot::Cartesian Atom_1 = build[i_current_build][n_fitted_in_current_segment-3].pos;
   coot::Cartesian Atom_2 = build[i_current_build][n_fitted_in_current_segment-2].pos;  
   coot::Cartesian Atom_3 = build[i_current_build][n_fitted_in_current_segment-1].pos; 

   return position_by_torsion(Atom_1, Atom_2, Atom_3, theta_2, torsion, dist); 

}

coot::Cartesian
BuildCas::position_by_torsion(coot::Cartesian Atom_1, coot::Cartesian Atom_2, coot::Cartesian Atom_3,
			      float theta_2, float torsion, float dist) const { 

   coot::Cartesian a1a2 = Atom_2 - Atom_1; 
   coot::Cartesian a2a3 = Atom_3 - Atom_2; 

   coot::Cartesian z_r_normal = a2a3; z_r_normal.normalize(); 
   coot::Cartesian y_r_normal = coot::cross_product(a1a2, a2a3); y_r_normal.normalize(); 
   coot::Cartesian x_r_normal = coot::cross_product(y_r_normal, z_r_normal); x_r_normal.normalize(); 

   float z_r = dist*sin(theta_2 - PI_BY_2); 

   float l = dist*cos(theta_2 - PI_BY_2); // consider when theta_2 is less than PI_BY_2
                                          // do we need a fabs here? 
                                          // Checked.  It seems not.
   
   float x_r = l*cos(torsion); 
   float y_r = l*sin(torsion); 

   coot::Cartesian x_r_vec = x_r_normal.by_scalar(x_r); 
   coot::Cartesian y_r_vec = y_r_normal.by_scalar(y_r); 
   coot::Cartesian z_r_vec = z_r_normal.by_scalar(z_r); 

   coot::Cartesian sum_bits = x_r_vec + y_r_vec + z_r_vec; 

   coot::Cartesian Atom_4 = Atom_3 + sum_bits; 

   return Atom_4; 
}


// Search through l1 looking for branch points.  A branch point 
// has 3 neighbours - except when this point is part of a smalltriangle
// (a smalltriangle is a cube with 3 connected face diagonals).
// 
// If we get a smalltriangle, then make the branch point be the middle
// of the vertices of the triangle.
// 
// And symmetry expand them! (but don't pass pass back the symmetry expanded ones)
// 
std::vector<coot::Cartesian>
BuildCas::find_branch_points(const clipper::Xmap<float> &map,
			     const clipper::Xmap<int>   &l1,
			     float cut_off)
{
   //
   std::vector<coot::Cartesian> vc;
   std::vector<coot::Cartesian> vc_small_tri; // we keep small triangles separate
                                   // and combine these vectors at the
                                   // end. 

   std::vector<coot::Cartesian> vc_not_new_tr; 
   int n_skel_neighbs;

   // a couple of counters, just for debugging really.
   // 
   int n_simple_branch = 0;
   int n_triangle_branch = 0; 
   int n_not_new_small_tri = 0; 

   // Face diagonals must be great than an edge (1) but less than
   // a body diagonal (sqrt(3)). 
   // 
   // limits between 1 and sqrt(2)
   //        between sqrt(2) and sqrt(3).
   // 
   // edge_neighb finds edges and face diagonals.
   // 
   
   clipper::Skeleton_basic::Neighbours neighb( map );
   clipper::Skeleton_basic::Neighbours face_neighb(map, 1.2, 2.5);
   clipper::Skeleton_basic::Neighbours edge_neighb(map, 0.8, 2.5);

   std::cout << "There are " << neighb.size() << " neighbs and "
	     << face_neighb.size() << " face_neighbs " << std::endl;
   
   std::cout << "The grid sampling is " 
	     << d_map_p->grid_sampling().format() << std::endl;

   std::cout << "The cell is :" 
	     << d_map_p->cell().a() << " " 
	     << d_map_p->cell().b() << " " 
	     << d_map_p->cell().c() << std::endl; 

   std::cout << "Grid points every: " 
	<< d_map_p->cell().a()/d_map_p->grid_sampling().nu() << " "
	<< d_map_p->cell().b()/d_map_p->grid_sampling().nv() << " "
	<< d_map_p->cell().c()/d_map_p->grid_sampling().nw() << " "
	     << "A" << std::endl; 
   std::cout << "find_branch_points is using map cut_off: " << cut_off << std::endl; 


   clipper::Xmap_base::Map_reference_index ix;
   clipper::Coord_grid c_g; 

   small_triangle_thing sm_tr_res; 

   std::vector<clipper::Coord_grid> vec_coord_grid_running; 
   
   for (ix = l1.first(); !ix.last(); ix.next() ) {
      
      n_skel_neighbs = 0; 
      if (l1[ix] > 0) {
	
	 for (int in=0; in < neighb.size(); in++ ) {

	    c_g = ix.coord() + neighb[in];
	    
	    if (l1.get_data(c_g) > 0) {		       
	       if (d_map_p->get_data(c_g) > cut_off) {

		  n_skel_neighbs++; 
	       }
	    }
	 }

	 if (n_skel_neighbs >= 3) {

// 	    if (isSmallTriangle(l1,map, cut_off,
// 				face_neighb, edge_neighb, ix.coord())) {

// This finds all the branch points. 
//
// 	    clipper::Coord_orth c_o = 
// 	       ix.coord().coord_frac(l1.grid_sampling()).coord_orth(l1.cell());
// 	    vc_small_tri.push_back(Cartesian(c_o.x(), c_o.y(), c_o.z())); 

	       
	    sm_tr_res =
	       isSmallTriangle_new(l1,map, cut_off,
				   face_neighb, edge_neighb, ix.coord(),
				   vc_small_tri); 
	       
	    if (sm_tr_res.is_small_triangle) { 
	       
	       // a branch point that was a small triangle.
	       
	       if (sm_tr_res.is_new_small_triangle) { 
		  
		  // which we haven't found before:
		  
		  vc_small_tri.push_back( sm_tr_res.pos);
		  n_triangle_branch++;
		  
	       } else { 
		  
		  n_not_new_small_tri++; 

		  clipper::Coord_orth c_o = 
		     ix.coord().coord_frac(l1.grid_sampling()).coord_orth(l1.cell());
		  vc_not_new_tr.push_back(coot::Cartesian(c_o.x(), c_o.y(), c_o.z())); 

	       }

	    } else { 

	       // a branch point that was NOT a small triangle
	       // i.e. is conventional-looking.
	       
	       clipper::Coord_orth c_o =
		  ix.coord().coord_frac(l1.grid_sampling()).coord_orth(l1.cell());
	       
	       vc.push_back(coot::Cartesian(c_o.x(), c_o.y(), c_o.z()));
	       
	       n_simple_branch++;
	    }
	 }
      }
   }

   std::cout << "DEBUG: there were " << std::endl << "     " 
	<< n_triangle_branch << " (new) triangle branches" << std::endl << "     "
	<< n_not_new_small_tri << " triangle branches that were previously "
	<< "considered " << std::endl << "     " << n_simple_branch   
	<< " simple branches " << std::endl;

   // copy vc to the data member branch_points
   // 
   // branch_points.resize() and branch_points[i] = vc[i] may be faster.
   // 
   for (unsigned int i=0; i< vc.size(); i++)
      branch_points.push_back(vc[i]);
   
   // asc_to_graphics(convert_to_atoms(vc_not_new_tr), "Not New Tris", ATOM_BONDS, 0.1);

   // asc_to_graphics(convert_to_atoms(vc), "Simple Branch Points", ATOM_BONDS, 0.1);

   // copy the small triangles 
   // 
   for (unsigned int i=0; i< vc_small_tri.size(); i++)
      branch_points.push_back(vc_small_tri[i]);

   // which uses the class variable branch_points, of course.
   symmetry_expand_branch_points(); 

   return vc; 
}



// Test face diagonal equilateral triangles as well as edge-edge-face-diagonal
// 
short int
BuildCas::isSmallTriangle(const clipper::Xmap<int> &l1,
			  const clipper::Xmap<float> &map,
			  float cut_off, 
			  const clipper::Skeleton_basic::Neighbours &fd_neighb,
			  const clipper::Skeleton_basic::Neighbours &edge_neighb,
			  const clipper::Coord_grid &pos) const {

   //
   int n_stn = 0; // small triangle neighbours. 

   for (int in=0; in<fd_neighb.size(); in++) {
      //

      clipper::Coord_grid c_g = pos + fd_neighb[in];

      if (l1.get_data(c_g) > 0 ) {
	 if (d_map_p->get_data(c_g) > cut_off) {
	    n_stn++;
	 }
      } 
   }

   if (n_stn > 3) {
      std::cout << "n_stn: " << n_stn << " at " << pos.format() << std::endl;
   }

   if (n_stn > 3) {
      return 1;
   } else { 
      return 0; 
   }

   // 
    
}

// Set is_small_triangle if it is so.
// set is_new_small_trangle if it was not considered before
// 
// set pos only if this is a new small triangle; 
// 
small_triangle_thing
BuildCas::isSmallTriangle_new(const clipper::Xmap<int> &l1,
			      const clipper::Xmap<float> &map,
			      float cut_off, 
			      const clipper::Skeleton_basic::Neighbours &fd_neighb,
			      const clipper::Skeleton_basic::Neighbours &edge_neighb,
			      const clipper::Coord_grid &pos, 
			      const std::vector<coot::Cartesian> &small_tri_running_list) const {

   // We already know (from find_branch_points) that pos is a branch
   // point before this function is called (indeed, this function is
   // called because of that fact).
   // 
   // So we want to know: are the neighbours of this point (pos)
   // neighbours of each other? (as Kevin suggests).
   // 
   // First find the neighbours that have more than 2 neighbours and
   // put them in a vector when we find them.  The run through the
   // vector looking for neighbours of the element of the vector.  If
   // one of those neighbours is in the vector list, then we have a
   // small triangle.
   // 

   // return this:
   small_triangle_thing smt;
   smt.is_small_triangle = 0; 

   std::vector<clipper::Coord_grid> vcg; 
   clipper::Coord_grid c_g, c_g_n;
   int i_n_neighbours;

   for (int in=0; in < fd_neighb.size(); in++ ) {

      c_g = pos + fd_neighb[in];

      if (l1.get_data(c_g) > 0 ) { 

	 if (d_map_p->get_data(c_g) > cut_off) {

	    // we are at a neighbour of pos that is a skeleton point

	    // now count the neighbours of *that*

	    i_n_neighbours = 0; 
	    
	    for (int i_n_n = 0; i_n_n < fd_neighb.size(); i_n_n++ ) {

	       c_g_n = c_g + fd_neighb[i_n_n]; 

	       if (l1.get_data(c_g_n) > 0) { 
		  
		  if (d_map_p->get_data(c_g) > cut_off) {

		     // the neighbour's neighbour was a skeleton point
		     i_n_neighbours++; 
		  }
	       }
	    }

	    if (i_n_neighbours > 2) { 

	       vcg.push_back(c_g); 

	    }
	 } 
      }
   }

   
   int isize = vcg.size(); 
   // cout << "vcg has " << isize << " elements" << endl ;

   short int i_seen_before_flag = 0; 
   clipper::Coord_orth c_o_1; 
   clipper::Coord_orth c_o_2; 
   clipper::Coord_orth c_o_3; 

   std::vector<coot::Cartesian> triangle_points; 
   coot::Cartesian tri_centre; 

   if (isize > 1) { 
      for (int i=0; i<isize; i++) { 
	 for (int j=i+1; j<isize; j++) { 

	    // does vcg[i] have a neigbour which is vcg[j]? 
	    // 
	    for (int in=0; in < fd_neighb.size(); in++ ) {

	       c_g = vcg[i] + fd_neighb[in];

	       // c_g == vcg[j] ? 

	       if (c_g == vcg[j]) { 

		  // OK, we have a small triangle.
		  // Now we need to check if we have set the centre point
		  // before.

		  smt.is_small_triangle = 1; 
		  
		  // First, lets find the centre point of the small triangle:
		  // 
		  
		  c_o_1 = 
		     pos.coord_frac(l1.grid_sampling()).coord_orth(l1.cell());
		  c_o_2 = 
		     vcg[j].coord_frac(l1.grid_sampling()).coord_orth(l1.cell());
		  c_o_3 = 
		     vcg[i].coord_frac(l1.grid_sampling()).coord_orth(l1.cell());

		  triangle_points.push_back(coot::Cartesian(c_o_1.x(), c_o_1.y(), c_o_1.z())); 
		  triangle_points.push_back(coot::Cartesian(c_o_2.x(), c_o_2.y(), c_o_2.z())); 
		  triangle_points.push_back(coot::Cartesian(c_o_3.x(), c_o_3.y(), c_o_3.z())); 

		  tri_centre = average_Cartesians(triangle_points); 

// 		  cout << "debug: averaging: " << endl 
// 		       << "     " << c_o_1.format() << endl 
// 		       << "     " << c_o_2.format() << endl 
// 		       << "     " << c_o_3.format() << endl 
// 		       << "gives" << endl
// 		       << "     " << tri_centre << endl; 
		  

		  i_seen_before_flag = 0; 

		  for(unsigned int k=0; k<small_tri_running_list.size(); k++) {
		     
		     if (fabs(tri_centre.x() - small_tri_running_list[k].x()) < 0.01) { 
			if (fabs(tri_centre.y() - small_tri_running_list[k].y()) < 0.01) { 
			   if (fabs(tri_centre.z() - small_tri_running_list[k].z()) < 0.01) { 
			      
			      i_seen_before_flag = 1; 
// 			      cout << "found a match " << tri_centre << " and " 
// 				   << small_tri_running_list[k] << endl; 
				 
			      break; 
			   }
			}
		     }
		  }

		  if (i_seen_before_flag == 0) { 

		     // i.e. was new
		     // lets set the return value:
		     
		     smt.is_new_small_triangle = 1; 
		     smt.pos = tri_centre; 
		     

		  } else { 

		     // it *was* seen before
		     // 
		     // Do nothing. (set a flag to do nothing)

		     smt.is_new_small_triangle = 0; 
		     
		  }
	       }
	    }
	 }
      }
   }

   return smt; 

}

void 
BuildCas::depth_search_skeleton_testing_2() { 

   make_tree_node_map(); 

   // this ca_1 is close to where we force the first atom.  It is the 
   // nearest branch point.
   // 
   // and ca_2 is what gets built. It is bang on a skeleton grid
   // point. They are for sure connected.  result must be positive if
   // we are to continue (using test map, obviously)
   // 
   clipper::Coord_orth ca_1(-0.698166,-8.4516,-13.8573); 
   clipper::Coord_orth ca_2(-4.189,-6.30823,-13.2358); 

   clipper::Coord_frac ca_f_1 = ca_1.coord_frac(d_map_p->cell());
   clipper::Coord_frac ca_f_2 = ca_2.coord_frac(d_map_p->cell());

   clipper::Coord_grid ca_g_1 = ca_f_1.coord_grid(d_map_p->grid_sampling());
   clipper::Coord_grid ca_g_2 = ca_f_2.coord_grid(d_map_p->grid_sampling());

   std::cout << "depth_search_skeleton for        " << ca_g_1.format() 
	<< " to " << ca_g_2.format() << std::endl; 

   std::cout << "depth_search_skeleton (unit) for " 
	<< ca_g_1.unit(d_map_p->grid_sampling()).format() << " to " 
	<< ca_g_2.unit(d_map_p->grid_sampling()).format() << std::endl; 

   std::cout << "initially: we get grid values: " << segment_map.get_data(ca_g_1)
	<< " and " << segment_map.get_data(ca_g_2) << std::endl;

   std::cout << "initially here are the neighbours of " << ca_g_1.format()
	<< std::endl; 

   for (unsigned int i=0; i<treenodemap.get_data(ca_g_1).neighbs.size(); i++) {
      std::cout << "      " << treenodemap.get_data(ca_g_1).neighbs[i].format() << std::endl; 
   }
   std::cout << "initially here are the neighbours of " << ca_g_2.format()
	<< std::endl; 

   for (unsigned int i=0; i<treenodemap.get_data(ca_g_2).neighbs.size(); i++) {
      std::cout << "      " << treenodemap.get_data(ca_g_2).neighbs[i].format() << std::endl; 
   }

   short int i = depth_search_skeleton(ca_g_1, ca_g_2); 
   
   std::cout << " result: " << i << std::endl; 
   
}

void 
BuildCas::depth_search_skeleton_testing() { 

   make_tree_node_map(); 

   clipper::Coord_grid start; 
   clipper::Coord_grid target; 

   clipper::Xmap_base::Map_reference_index ix;

   int count; 
   short int isearch; 

   // note i=0 is bad.
   // (because start and target do not get set). 
   for (int i=1; i<=200; i++) { 

      count = 0; 

      for (ix = segment_map.first(); !ix.last(); ix.next() ) {

	 if (segment_map[ix] > 0) { 

	    if ((*d_map_p)[ix] > map_cut_off) { 

	       count++; 
	       if (count == i) 
		  start = ix.coord(); 
	       
	       if (count == 2*i) { 
		  target = ix.coord(); 
		  break; 
	       }
	    }
	 }
      }

      isearch = depth_search_skeleton(start, target); 

      std::cout << "result of that: testing " << start.format() 
	   << " to  " << target.format() << " is "; 
      std::cout << isearch << std::endl << std::endl; 
   } 
} 

short int
BuildCas::depth_search_skeleton(const clipper::Coord_grid &start, 
				const clipper::Coord_grid &target) const { 


   clipper::Coord_grid prev_prev;  // dummy for first usage
   clipper::Coord_grid previous;   //  ditto.

   if (treenodemap.get_data(start).neighbs.size() == 0) {

      std::cout << "woops! no neighbours for depth search start " << start.format() << std::endl; 
   }


   short int i = depth_search_skeleton_internal(start, previous, prev_prev, target, 10, 0); 

   return i; 

}

// we currently presume that previous and target and prev_prev (except
// for dummy first usage) are real skeleton points.
// 
// We should test segment_map.get_data(previous) if this is not so. 
// 
short int
BuildCas::depth_search_skeleton_internal(const clipper::Coord_grid &current,
					 const clipper::Coord_grid &previous,
					 const clipper::Coord_grid &prev_prev,
					 const clipper::Coord_grid &target, 
					 int depth, int length) const { 

//    cout << "current is " << current.format() 
// 	<< " and depth is " << depth 
// 	<< " and length is " << length << endl; 

   int j; 
   int imatch = 0; 

   if (depth == 0) { 

      return is_same_gridpoint(current, target); 
      
   } else { 

      if ( is_same_gridpoint(current,target)) { 
	 
	 std::cout << "!!! A depth_search_skeleton_internal hit at " << current.format() << std::endl;
	 return 1; 
	  
      } else { 
 
	 j = 0; 

	 //
	 if (treenodemap.get_data(current).neighbs.size() == 0) {

	    std::cout << "woops! no neighbours for " << current.format() << std::endl; 
	 }

	 // 

	 clipper::Coord_frac cf = current.coord_frac(d_map_p->grid_sampling()); 
	 clipper::Coord_orth co = cf.coord_orth(d_map_p->cell()); 
	 
// 	 cout << "testing the " << treenodemap.get_data(current).neighbs.size() 
// 	      << " neighbs of " << current.format() << " at position " << co.format() << endl; 
	 for(unsigned int i=0; i<treenodemap.get_data(current).neighbs.size(); i++) {

	    if (length > 2) {
	       
	       // 
	       // not a small triangle
	       // going round in circles (or triangles) gets us nowhere
	       
	       if ( ! (treenodemap.get_data(current).neighbs[i] == prev_prev)) {
		  
		  // forward and back, not interesting, kill it
		  
		  if (! (treenodemap.get_data(current).neighbs[i] == previous)) { 
		     
		     
		     j = depth_search_skeleton_internal(treenodemap.get_data(current).neighbs[i], 
							current,
							previous, 
							target, 
							depth - 1, 
							length + 1);
		  }
	       }

	    } else {
	       
	       // short length
	       // 
	       if (length > 1) { 
		  // forward and back, not interesting, kill it
		  
		  if (! (treenodemap.get_data(current).neighbs[i] == previous)) {

		     j = depth_search_skeleton_internal(treenodemap.get_data(current).neighbs[i], 
							current,
							previous, 
							target, 
							depth - 1, 
							length + 1); 
		  }
	       } else { 

		  // we do this only at the start

		  j = depth_search_skeleton_internal(treenodemap.get_data(current).neighbs[i], 
						     current,
						     previous, // dummy first time
						     target, 
						     depth - 1, 
						     length + 1);
	       }
	    }

	    // 
	    imatch += j;

	    if (imatch) 
	       break; 

	 }
	 return imatch; 
      }
   }
}


coot::Cartesian
BuildCas::SmallTriangle_to_branch_point(const clipper::Xmap<int> &l1,
					const clipper::Skeleton_basic::Neighbours &fd_neighb,
					const clipper::Coord_grid &pos) const {
   coot::Cartesian a; 
   std::vector<coot::Cartesian> vc;

   clipper::Coord_orth pos_o =
      pos.coord_frac(l1.grid_sampling()).coord_orth(l1.cell());
   
   vc.push_back(coot::Cartesian(pos_o.x(), pos_o.y(), pos_o.z()));
   
   for (int in = 0; in < fd_neighb.size(); in++ ) {
      //
      clipper::Coord_grid c_g = pos + fd_neighb[in];
      
      if (l1.get_data(c_g) > 0 ) {
	 
	 clipper::Coord_orth c_o = c_g.coord_frac(l1.grid_sampling()).coord_orth(l1.cell());
	 
	 vc.push_back(coot::Cartesian(c_o.x(), c_o.y(), c_o.z()));
	 
      }
   }
   
//    cout << "vc is of size: " << vc.size() << endl; 
//    cout << "averaging: " << endl; 
//    for (int i=0; i<vc.size(); i++)
//       cout << "    " << vc[i] << endl;
      
   return average_Cartesians(vc);
   // return vc[0]; // "pos" as a Cartesian
   
}


// we presume that segment_map has been inited by now.
// 
// i.e. count_and_mark_segments has been called already.
// 
atom_selection_container_t 
BuildCas::grown_Cas() const { 

   std::cout << "sample grown vectors: " << std::endl; 
   for (int ii=1; ii<10; ii++)
      std::cout << "grown_ca vector " << ii << " " 
	   << build[i_current_build][ii].pos <<  " " 
	   << build[i_current_build][ii].near_grid_point.format() << std::endl;

   if (! segment_map_filled) 
      std::cout << " !!!!! WARNING !!!!!" << " garbage grown atoms " << std::endl; 

   std::vector<coot::Cartesian> built_atom_vector;

   // We start building at residue one these days.
   // 
   for (unsigned int i=1; i< build[i_current_build].size(); i++)
      built_atom_vector.push_back(build[i_current_build][i].pos); 

      
   return convert_to_atoms_internal(segment_map.spacegroup(), 
				    segment_map.cell(), 
				    built_atom_vector,
				    1,
				    "GROWN CAS"); // put each atom in a different residue
}

void 
BuildCas::ca_grow_recursive() { 

   int n_trials = 5; 

   build.resize(n_trials); 

   score_and_cart score_and_cartesian; 

   // The first point is special.  We need to build up from several
   // different starting points and choose the best (of course).
   // 
   std::cout << "--------------------------------------------------" << std::endl; 
   std::cout << "         Creating Skeleton Tree Node Map          " << std::endl;
   std::cout << "--------------------------------------------------" << std::endl; 

   make_tree_node_map(); 
   std::cout << "done" << std::endl; 

   std::cout << "--------------------------------------------------" << std::endl; 
   std::cout << "               Finding first point" << std::endl;
   std::cout << "--------------------------------------------------" << std::endl; 

   // score_and_cartesian = build_first_recursive(); 
   score_and_cartesian = build_first_cheat(); 

   build[i_current_build].resize(100); 
   build[0][1] = score_and_cartesian; 

   // put these back (they were manipulated in build_first_recursive(); 
   i_max_build = 1; // maximum number of builds, not the top of the array.
   i_current_build = 0; 
   
   TreeNode first_node(score_and_cartesian); // top of the tree node

   std::cout << "--------------------------------------------------" << std::endl; 
   std::cout << "               Done first point" << std::endl;
   std::cout << "--------------------------------------------------" << std::endl; 




   TreeNode *previous_node = &first_node; 

   for (int ires=2; ires<=40; ires++) {

      std::cout << "=-=-=-=-=-=-=-=-=- fitting residue " << ires 
	   << " =-=-=-=-=-=-=-=-=-" << std::endl;
      
      score_and_cartesian = recursive_build(previous_node, ires,1);

      if (score_and_cartesian.score <= 0) { 
	 std::cout << "bad score: " << score_and_cartesian.score << std::endl
	      << "Building terminated" << std::endl; 
	 break; 

      } else { 

	 std::cout << "=-=-=-=-=-=-=-=-=- fitted a point =-=-=-=-=-=-=-=-=-" 
	      << std::endl
	      << score_and_cartesian.pos 
	      << " with score " << score_and_cartesian.score << std::endl; 
	 std::cout << "saving with near_grid_point " 
	      << score_and_cartesian.near_grid_point.format()
	      << std::endl; 

	 if (treenodemap.get_data(score_and_cartesian.near_grid_point).neighbs.size() == 0) { 
	    std::cout << "woops! in ca_grow_recursive: fitted point near_grid_point no neighbs" << std::endl; 
	 }

	 build[i_current_build][ires] = score_and_cartesian;
      
	 TreeNode *current_node = new TreeNode; // where does this get deleted?
	 current_node->setup(previous_node, score_and_cartesian); 

	 // and now for next round:
	 previous_node = current_node; 
      }
   }
}

score_and_cart
BuildCas::build_first_cheat() { 

   // CA /1/chainid="A"/1194/ASP/

   score_and_cart sc; 


   coot::Cartesian a(-0.698166,-8.4516,-13.8573); 
   
   clipper::Coord_orth ca_1(-0.698166,-8.4516,-13.8573); 
   clipper::Coord_frac ca_f_1 = ca_1.coord_frac(d_map_p->cell());
   clipper::Coord_grid ca_g_1 = ca_f_1.coord_grid(d_map_p->grid_sampling());

   sc.near_grid_point = ca_g_1; 
   sc.pos = a; 
   sc.score = 1; 

   // note we must make sure that sc.near_grid_point is not zero.
   // 
   std::cout << "build_first_cheat: gives grid " << sc.near_grid_point.format() 
	<< " and segment map value: " << segment_map.get_data(sc.near_grid_point) << std::endl; 


   if (segment_map.get_data(sc.near_grid_point) <= 0 ) { 

      std::cout << "ERROR ERROR ERROR ERROR ERROR : unexpected zero grid" << std::endl;
      std::cout << "ERROR ERROR ERROR ERROR ERROR " << std::endl;

   } 


   // we need to push back a score_and_cart
   // build[0].push_back(sc); 

   return sc; 
}

// Here is where we deal with the special case that is the first point. 
// We use the i_current_build and i_max_build. 
// 
// We need to check several starting points ("recursive_build"ing
// them) and then return the best one.
//
// We use prebuilt_exclusion_score to move down the list of good
// branch_points_symm_expanded - but that is done in
// peak_search_simple()
// 
// 
score_and_cart
BuildCas::build_first_recursive() { 

   score_and_cart best; 

   best.score = 0; 

   int n_trials = 5; 

   // build[][] is used in prebuilt_exclusion_score:
   build.resize(n_trials); 
   for (int i=0; i<n_trials; i++)
      std::cout << "DEBUG: build[" << i << "] size is: " << build[i].size() << std::endl; 

   std::vector<score_and_cart> sc_vec(n_trials); 
   std::vector<score_and_cart> peak_search_raw(n_trials); 

   for (int i=0; i<n_trials; i++) {

      i_current_build = i; 
      i_max_build = i+1; 
      std::cout << "DEBUG: i_max_build is " << i_max_build << std::endl; 

      peak_search_raw[i] = peak_search_simple(); 
      
      TreeNode first_node(peak_search_raw[i].pos); 
      // put it in the build array so that it excludes itself when
      // peak_search_simple() finds the next peak. 
      build[i_current_build].push_back(peak_search_raw[i]); 
      sc_vec[i] = recursive_build(&first_node, 2,1); 
   }

   std::cout << "-------------- --------------- " << std::endl; 
   std::cout << "Here are the intial builds: " << std::endl; 
   std::cout << "-------------- --------------- " << std::endl; 
   for (int i=0; i<n_trials; i++)
      std::cout << sc_vec[i].pos << " with raw pos: " << peak_search_raw[i].pos
	   << " with score " << sc_vec[i].score << std::endl;



   // I would use map here:
   // 
   for (int i=0; i<n_trials; i++) {
      
      if (sc_vec[i].score > best.score) { 
	 best = sc_vec[i]; 
      }
   }
   
   std::cout << "-------------- --------------- " << std::endl; 
   std::cout << "build_first_recursive selecting:" << std::endl; 
   std::cout << "    " << best.pos << " with score: " << best.score << std::endl; 
   std::cout << "-------------- --------------- " << std::endl; 
   
   return best; 

} 

// compare with build_big_ball.  We ought to combine at some stage.
// 
void
BuildCas::symmetry_expand_branch_points() {

   if (segment_map_filled == 1) { 
   
      atom_selection_container_t asc = 
	 convert_to_atoms_internal(segment_map.spacegroup(), 
				   segment_map.cell(), 
				   branch_points,
				   1, // all atoms in different residues
				   "symmetry branch points"); 

      if (! expansion_centre_is_set) { 
	 std::cout << "ERROR ERROR! Need to set expansion centre first!" << std::endl; 
      } 

      std::vector<clipper::Coord_grid> dummy_vec; // hmm... 
      
      branch_points_symm_expanded = 
	 // point_list_by_symmetry(asc, expansion_centre,  50.0);

	 // 0 mean don't use the dummy_vec. 
	 // 
	 point_list_by_symmetry(asc, dummy_vec, expansion_centre, 50.0, 0); 
      branch_point_have_been_expanded_flag = 1; 
      
   } else { 
      
      std::cout << "Error - need to fill segment map first "
	   << "(symmetry_expand_branch_points)" 
	   << std::endl; 
   }
}

//  We need to expand the branch points like we expanded the
//  skeleton points, and thus make the data member
//  branch_points_symm_expanded.
// 
float 
BuildCas::branch_point_proximity_score(coot::Cartesian trial_point) const { 

   float min_dist = 9999999.9; ; 
   float dist; 

   // there should be 1000 or 1000s.
   // 
   if (branch_point_have_been_expanded_flag == 0) { 
      std::cout << "Error - branch_points need symmetry expanding first" << std::endl; 
   }

   if ( ! (branch_points_symm_expanded.size() > 1)) { 
      std::cout << "!!! WARNING !!! branch_points_symm_expanded.size() is " 
	   << branch_points_symm_expanded.size() << std::endl; 
   } 

   for (unsigned int i=0; i< branch_points_symm_expanded.size(); i++) {

      if ( fabs(branch_points_symm_expanded[i].x() - trial_point.x()) < 4) {
	 if (fabs(branch_points_symm_expanded[i].y() - trial_point.y()) < 4) {
	    if (fabs(branch_points_symm_expanded[i].z() - trial_point.z()) < 4) {

	       dist = (branch_points_symm_expanded[i] - trial_point).amplitude();

	       if (dist < min_dist) 
		  min_dist = dist; 
	    }
	 }
      }
   }

   // I think that this should be a Fermi-Dirac
   // 

   // float score = pow((min_dist + 0.1),-2); 

   float score = 1/(min_dist + 0.3); 

//    if (min_dist < 2) { 
//       score = 10.0; 
//    } else { 
//       score = 1.0; 
//    } 

   return score; 
}



float
BuildCas::prebuilt_exclusion_score(coot::Cartesian trial_point) const { 

   float min_dist = 9999999.9;
   float dist; 

   float dist_crit = 2.0; // 
   short int ibreakflag = 0; 

   // recall that i_current_build starts at 0 (i.e. 0 has a
   // built/building segment in it). 

   for (int i=0; i<i_max_build; i++) { 

      // std::cout << "DEBUG: prebuilt_exclusion_score: i=" << i << std::endl; 
      
      for (unsigned int n=0; n< build[i].size(); n++) { 

	 // cout << "DEBUG: prebuilt_exclusion_score: n=" << n << endl; 
	 // cout << "build[" << i << "][" << n <<"] is " << build[i][n] << endl; 
	 dist = (build[i][n].pos - trial_point).amplitude(); 

	 if (dist < min_dist) 
	    min_dist = dist; 
	 
	 if (min_dist < dist_crit) { 
	    ibreakflag = 1; 
	    break; 
	 }
      }
      if (ibreakflag)
	 break; 
   }

   if (min_dist > 9999999.0) {
      std::cout << "!!!! WARNING !!!! prebuilt atoms not found "
	   << "in prebuilt_exclusion_score. " << std::endl; 
      min_dist = 9.9; 
   }

   float score; 

   if (min_dist < dist_crit) { 
      score = 0;
   } else { 
      score = 1; 
   }

   //cout << "prebuilt_exclusion_score: " << score << endl; 

   return score; 
} 


// For now, while we are debugging, we well return a "scores". Later,
// we can optimize (perhaps) by returning only a float.
// (We use a scores so that we can debug the scoring upstream).
// 
scores
BuildCas::non_angle_micro_point_score(coot::Cartesian previous_atom, 
				      coot::Cartesian trial_point) const { 

   scores sc; 
   float dv; 

   // segment_map must be filled by now:
   // 
   if (segment_map_filled == 0)
      std::cout << "Error: must fill segment map before "
	   << "non_angle_micro_point_score" << std::endl; 
   
   dv = density_at_point(trial_point); 
   
//    Cartesian bit(0.01, 0.01, 0.01); 
//    int t0 = glutGet(GLUT_ELAPSED_TIME); 
//    for (int i=0; i<1000; i++)
//       dv = density_at_point(trial_point+bit.by_scalar(i)); 
//    int t1 = glutGet(GLUT_ELAPSED_TIME); 
//    cout << "density at point   time: " << t1-t0 << "ms" << endl; 

//    t0 = glutGet(GLUT_ELAPSED_TIME); 
//    for (int i=0; i<1000; i++)
//       branch_point_proximity_score(trial_point);
//    t1 = glutGet(GLUT_ELAPSED_TIME); 
//    cout << "branch point prox. time: " << t1-t0 << "ms" << endl; 

   
   float d = (previous_atom - trial_point).amplitude(); 

   sc.deviation_from_ideal_length_score_val = 
      deviation_from_ideal_length_score(d); 

   sc.branch_point_proximity_score_val =
      branch_point_proximity_score(trial_point); 

   double e = 2.71828;
   sc.density_score_val  = pow(e, (double) 2*dv); // guess

   sc.score = sc.density_score_val
      * sc.deviation_from_ideal_length_score_val
      * sc.branch_point_proximity_score_val;

   return sc; 

}

//
float
BuildCas::interconnectedness(int ntips) const { 

   float v = 0.0; 

   if (branch_points.size() == 0) { 

      std::cout << "interconnectedness: must have branch_points first" << std::endl; 

   } else { 

      if (ntips == 0) { 

	 std::cout << "interconnectedness: must have some non-zero number of tips" << std::endl; 

      } else { 

	 if (! segment_map_filled) { 

	    std::cout << "interconnectedness: must fill the segment_map first" << std::endl; 

	 } else { 

	    int n_skel_pts = 0; 

	    clipper::Xmap_base::Map_reference_index ix;
	 
	    for (ix = segment_map.first(); !ix.last(); ix.next() )
	       if (segment_map[ix] > 0)
		  n_skel_pts++; 

	    std::cout << "interconnectedness: " << std::endl
		 << "    number of branch points: " << branch_points.size() << std::endl
		 << "    number of tips (passed): " << ntips << std::endl
		 << "    number of (segment) skeletoned points:  " << n_skel_pts << std::endl; 
		  
	    v = (float (branch_points.size() - ntips)/float(n_skel_pts)); 
	 }
      }
   }

   return v; 
}


//
float 
BuildCas::mid_point_density_score(coot::Cartesian prev, 
				  coot::Cartesian trial) const { 

   coot::Cartesian mid_point = prev.mid_point(trial); 

   float dv = density_at_point(mid_point); 

   // cout << "DEBUG: mid_point density: " << dv << endl; 

   float score = exp(2.0*(dv-map_cut_off)-2.3);  // -2.3 is /10. 

   return score; 
} 


float 
BuildCas::mid_points_density_score(coot::Cartesian prev, 
				   coot::Cartesian trial) const { 

   // 
   std::vector<coot::Cartesian> mid_points = prev.third_points(trial); 

//    cout << "debug: midpoints of " 
// 	<< "      " << prev << " and " << endl
// 	<< "      " << trial << " are " << endl
// 	<< "      " << mid_points[0] << " and " << endl
// 	<< "      " << mid_points[1] << endl; 

   std::cout << "debug: density values: " << std::endl
	<< "               " << density_at_point(prev)  << std::endl 
	<< "               " << density_at_point(trial) << std::endl 
	<< "               " << density_at_point(mid_points[0]) << std::endl 
	<< "               " << density_at_point(mid_points[1]) << std::endl ; 
      

   float score = 
      exp(5.0*(density_at_point(mid_points[0])-map_cut_off)) * 
      exp(5.0*(density_at_point(mid_points[1])-map_cut_off)); 

   return score; 
}

float
BuildCas::density_at_point(coot::Cartesian trial_point) const { 

   float dv; 

   clipper::Coord_orth c_o(trial_point.x(), trial_point.y(), trial_point.z()); 
   
   clipper::Coord_frac c_f = c_o.coord_frac(segment_map.cell()); 
   clipper::Coord_map  c_m = c_f.coord_map (segment_map.grid_sampling()); 
   clipper::Interp_linear::interp(*d_map_p, c_m, dv); 

   return dv; 

} 


// 
std::vector<coot::Cartesian_and_Grid>
BuildCas::fitting_targets(const TreeNode *node, float max_gridding) const { 

   // we return this:
   //
   std::vector<coot::Cartesian_and_Grid> cluster_centres_vec; 

   // we should use 3.7 +/- bit; where bit is maximum gridding (an
   // xmap-utils feature really) (say 1.3 for a 3A map) multiplied by
   // some wiggle-room factor (1.3, say) - that accounts for errors in
   // the previous position.

   float wiggle_factor = 1.0; 

   coot::Cartesian ires_pos = node->pos; 

   // 
   std::vector<coot::Cartesian_and_Grid> points_within_distance = 
      select_by_distance(ires_pos, 
			 3.7-wiggle_factor*max_gridding, 
			 3.7+wiggle_factor*max_gridding); // 3.7A +/- a bit
   // create a vector of vectors of points that are clustered together
   // on a skeleton.
   //
   std::vector<std::vector<coot::Cartesian_and_Grid> > cluster_vec = 
      cluster_bones_points(points_within_distance, ires_pos); 

   if (cluster_vec.size() == 0) { 
      // cannot happen
      // 
      std::cout << "INFO:: There are no clusters.  This cannot happen" << std::endl; 
      
   } else { 

      cluster_centres_vec = cluster_centres(cluster_vec); 

      std::cout << "INFO:: There were " << cluster_centres_vec.size() 
	   << " cluster center vectors in fitting_targets" << std::endl; 

      if (cluster_centres_vec.size() > 10) { 
	 //
	 std::cout << "WARNING: (fitting targets) " 
	      << "strange (clustered) vectors:" << std::endl; 
	 // 
	 // 
	 for(unsigned int i=0; i < cluster_centres_vec.size(); i++) { 
	    std::cout << "cluster center " << i << " at " 
		 << cluster_centres_vec[i].pos << std::endl;
	    for(unsigned int j=0; j<cluster_vec[i].size(); j++) { 
	       std::cout << "     " << cluster_vec[i][j].pos << std::endl; 
	    }
	 } 
      }
   }

   for(unsigned int i=0; i < cluster_centres_vec.size(); i++) { 
      if (treenodemap.get_data(cluster_centres_vec[i].near_grid_point).neighbs.size() == 0) { 
	 std::cout << "woops! in fitting targets element " <<  i << " near_grid_point no neighbs: "
	      << cluster_centres_vec[i].near_grid_point.format() << std::endl; 
      }
   }

   return cluster_centres_vec; 

} 


// The return value needs to contain a near_grid_point too now.  
// Which means that all the peak_search_*s should have it.
// 
score_and_cart
BuildCas::peak_search_wrapper(const TreeNode *node, int ith_res, int depth) { 

   std::cout << "DEBUG: ith_res: " << ith_res 
	<< " node pos " << node->pos << " depth " << depth << std::endl; 

   score_and_cart score_and_cartesian; 

   
   switch (ith_res) {

   case 0:
      std::cout << "SHOULD NEVER HAPPEN! BADNESS! ith_res is 0 " << std::endl; 
      break; 

   case 1:
      std::cout << "SHOULD NEVER HAPPEN! BADNESS!  ith_res is 1 " << std::endl; 
      break; 

   case 2:
      score_and_cartesian = 
	 peak_search_distance(node->parent->pos, node->pos); 
      break; 

   case 3:
      
      if (node->parent->parent == NULL) { 
	 std::cout << "ERROR null node->parent->parent.  " 
	      << "This should not happen" << std::endl; 
      } 
      score_and_cartesian = 
	 peak_search_distance_theta_2(node); 
      break; 

   default:
      score_and_cartesian = 
	 peak_search_distance_angle_torsion(node); 
      break; 

   } 
   if (score_and_cartesian.score > 0) { 
      if (treenodemap.get_data(score_and_cartesian.near_grid_point).neighbs.size() == 0) { 
	 std::cout << "woops! in peak_search_wrapper: fitted point near_grid_point no neighbs: "
	      << score_and_cartesian.near_grid_point.format() << std::endl; 
      }
   }

   return score_and_cartesian; 
}

// Return a score_and_cartesian
//
// When we get called, we have build upto the (ith_res -1)th residue.
// node contains the position of the (ith_res -1)th residue. 
// 
// ith_res is the Ca to be build this round. Therefore (ith_res -1) is
// the previous one (from which we make the new_targets (cluster
// centre vector).
// 
// depth should be greater than 0 in the outer call.
// 
// This is also called "dynamic evalution" in the literature...
// 
// ... which then goes on to describing MiniMax and the need for good
// pruning (I'm not sure we need Minimax here, but good pruning sounds
// worthwhile).  
// 
// get rid of depth in the call to peak_search_wrapper when this works.
// 
// We need to return a coord_grid in near_grid_point in the return
// value.   i.e. peak_search_wrapper should have it.
// 
score_and_cart
BuildCas::recursive_build(const TreeNode *node, int ith_res, int depth) {

   score_and_cart score_and_cartesian; 
   float max_gridding = maximum_gridding(); // take me outside for more speed
                                            // consider a class data member.

   // beat this:
   score_and_cartesian.score = 0; 


   std::vector<coot::Cartesian> fourth_set;       // debugging (temporary).
   std::vector<coot::Cartesian> fourth_set_grids; //  ditto. 

   if (depth == 0) {
      score_and_cartesian.score = 1;
      return score_and_cartesian;
   } else {

      std::vector<coot::Cartesian_and_Grid> new_targets = 
	 fitting_targets(node, max_gridding);
      
      // now make those fitting targets nodes:
      // 

      // TreeNode new_node(new_targets.size()); 
      std::vector<TreeNode> new_node(new_targets.size()); 

      for (unsigned int i=0; i< new_targets.size(); i++) { 

	 // The grid point of the new_targets[i] is actually at the
	 // same place as the Cartesian - we could bring it with us
	 // from fitting_targets for extra speed, but for now we will
	 // go back using clipper.  The resulting Coord_grid should
	 // *always* be on a skeleton point.
	 //
	 // See the FIXME below.  We want the same c_g.

// 	 clipper::Coord_orth c_o(new_targets[i].pos.x(), 
// 				 new_targets[i].pos.y(), 
// 				 new_targets[i].pos.z());
// 	 clipper::Coord_frac c_f = c_o.coord_frac(d_map_p->cell()); 
//       clipper::Coord_grid c_g = c_f.coord_grid(d_map_p->grid_sampling()); 

	 clipper::Coord_grid c_g = new_targets[i].near_grid_point; 

// 	 // debugging
	 if (ith_res == 30) { 

	    clipper::Coord_frac cfg = 
	       new_targets[i].near_grid_point.coord_frac(d_map_p->grid_sampling());
	    clipper::Coord_orth cfo = cfg.coord_orth(d_map_p->cell()); 

	    fourth_set.push_back(new_targets[i].pos); 
	    fourth_set_grids.push_back(coot::Cartesian(cfo.x(), cfo.y(), cfo.z()));
	    
	    std::cout << "DEBUG: fourth " << new_targets[i].near_grid_point.format()
		 << " to " << coot::Cartesian(cfo.x(), cfo.y(), cfo.z()) 
		 << " with pos " << new_targets[i].pos << std::endl;
	    std::cout << "grid closeness: " 
		 << treenodemap.get_data(new_targets[i].near_grid_point).near_grid_point.format()
		 << " "
		 << treenodemap.get_data(node->near_grid_point).near_grid_point.format()
		 << std::endl; 
	 }

	 // These are fine
// 	 cout << "new_targets[" << i << "].near_grid_point is " 
// 	      << new_targets[i].near_grid_point.format() << endl; 

	 if (! (segment_map.get_data(c_g) > 0) ) { 

	    std::cout << std::endl
		 << "ERROR! ERROR! ERROR! recursive_build c_g "
		 << c_g.format() << " is not a skeleton point. " 
		 << "value: " << segment_map.get_data(c_g) << std::endl << std::endl; 
	 } 

	 new_node[i].setup(node, new_targets[i].pos, new_targets[i].near_grid_point); 

      }


      for (unsigned int i=0; i< new_node.size(); i++) { 
	 if (treenodemap.get_data(new_node[i].near_grid_point).neighbs.size() == 0) { 
	    std::cout << "woops! in recursive_build: new node messed up: " 
		 << " for " << new_node[i].near_grid_point.format() << std::endl; 
	    if (treenodemap.get_data(new_targets[i].near_grid_point).neighbs.size() == 0) { 
	       std::cout << "woops! in recursive_build: new targets messed up too: " 
		    << " for " << new_targets[i].near_grid_point.format() << std::endl; 
	    } else { 
	       std::cout << "woops! in recursive_build: new targets *NOT* messed up!!!: " 
		    << " for " << new_targets[i].near_grid_point.format() << std::endl; 
	    }
	 }
      }


      // convert cluster_centres_vec to a vector of micropoint 
      // refined vector<Cartesian>
      // 
      
      std::vector<score_and_cart> score_and_cart_vec(new_targets.size()); 
      score_and_cart tmp_sc; 
      
      for (unsigned int i=0; i< new_targets.size(); i++) { 

	 // dynamic score the new_targets: 
	 // 
	 // and copy over the near_grid_point from a TreeNode to a scores.
	 // 

	 tmp_sc = peak_search_wrapper(&new_node[i], ith_res, depth); 
	 
	 std::cout << "post peak_search_wrapper: " << tmp_sc.near_grid_point.format() << std::endl; 
	 
	 score_and_cart_vec[i].score = 
	    tmp_sc.score *
	    recursive_build(&new_node[i], ith_res +1, depth-1).score; 
	 
	 score_and_cart_vec[i].pos = tmp_sc.pos;
	 score_and_cart_vec[i].near_grid_point = tmp_sc.near_grid_point; 

      }


      // pick and set the best
      // 
      for (unsigned int i=0; i< new_targets.size(); i++) { 
	 if (score_and_cart_vec[i].score > score_and_cartesian.score) { 
	    score_and_cartesian = score_and_cart_vec[i]; 

	    // redundant - FIXME
	    score_and_cartesian.near_grid_point = score_and_cart_vec[i].near_grid_point; 
	 }
      }
   }

   if (treenodemap.get_data(score_and_cartesian.near_grid_point).neighbs.size() == 0) { 
      std::cout << "woops! in recursive_build: returning near_grid_point no neighbs" 
	   << " for " << score_and_cartesian.near_grid_point.format() << std::endl; 
   }


   return score_and_cartesian; 

}

void
BuildCas::export_coordinates(atom_selection_container_t asc, 
			     std::string filename) const { 

   int err = asc.mol->WritePDBASCII(filename.c_str()); 
   
   if (err) { 
      std::cout << "There was an error in writing " << filename << std::endl; 
   } 

} 


void
BuildCas::make_tree_node_map() { 

   // Chomp up the lovely memory! Yum!
   // 
   treenodemap.init(d_map_p->spacegroup(), 
		d_map_p->cell(), 
		d_map_p->grid_sampling()); 

   clipper::Coord_grid c_g; 


   clipper::Skeleton_basic::Neighbours skel_neighbs(segment_map); 
      
   clipper::Xmap_base::Map_reference_index ix;

   for (ix = segment_map.first(); !ix.last(); ix.next() ) {

      if (segment_map[ix] > 0) { 

	 if ((*d_map_p)[ix] > map_cut_off) { 
	    
	    // we are at a skeleton point
	    //
	    SkeletonTreeNode stn; 
	    
	    for(int i=0; i<skel_neighbs.size(); i++) {

	       c_g = ix.coord() + skel_neighbs[i]; 

// 	       std::cout << "adding " << skel_neighbs[i].format() << " to "
// 		    << ix.coord().format() << " gives " << c_g.format() << endl;
	       
	       if (segment_map.get_data(c_g) > 0 ) { 

		  // map_cut_off is set at set_density_map_and_cut time.
		  // (i.e. constuction of a BuildCas time)
		  // 
		  if ( d_map_p->get_data(c_g) > map_cut_off ) { 

		     // OK, so this node has a neighbour:
		     // 
		     stn.neighbs.push_back(c_g); 
		  }
	       }
	    }
	    stn.near_grid_point = ix.coord();  // Strange but true!
	                                       // 
	                                       // We do this because "out of cell" reference
	                                       // (e.g.  uvw = (  -1, -12, -19)) will get wrapped 
	                                       // to some (hidden) value.  To get the wrapped
	                                       // value (i.e the grid), we look it up here. 
	                                       // Cunning (if it works). 
	    treenodemap[ix] = stn; 
	 }
      }
   }
 
   treenodemap_is_filled = 1; // they think it's all filled... it is now.

}
