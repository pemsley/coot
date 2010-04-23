/* coords/Bond_lines.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2009 by the University of Oxford
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
 
#include <string.h> // for strcmp
#include <stdlib.h> // for labs
#include <string>
#include <fstream>
#include <vector>

#include "Cartesian.h"
#include "mmdb_manager.h"
#include "mmdb-extras.h"
#include "mmdb.h"
#include "mmdb-crystal.h"  // should be merged with extras

#include "Bond_lines.h"
#include "coot-coord-utils.hh"


static std::string b_factor_bonds_scale_handle_name = "B-factor-bonds-scale";

Bond_lines::Bond_lines(coot::CartesianPair pair) {

   points.push_back(pair);

}

// We arrange things like this because the other constructor now uses
// construct_from_asc() too.
//
// This is a tiny bit clumsy having so many constructors.
// Heyho - historical cruft.
// 
Bond_lines_container::Bond_lines_container(const atom_selection_container_t &SelAtom,
					   int do_disulphide_bonds_in,   // default argument
					   int do_bonds_to_hydrogens_in) // default argument
{
   // how to initiaize bonds to 0.
   //
   do_disulfide_bonds_flag = do_disulphide_bonds_in;
   do_bonds_to_hydrogens = do_bonds_to_hydrogens_in;

//    std::cout << "====================== Blc constructor "
// 	     << do_bonds_to_hydrogens << " =======" << std::endl;
   
   b_factor_scale = 1.0;
   // 1.7 will not catch MET bonds (1.791 and 1.803) nor MSE bonds (1.95)
   // but SO4 bonds (1.46 are fine).
   // They should have special case, handle_MET_or_MSE_case
   // However, for VNP thingy, S1 has bonds to carbons of 1.67 1.77.  Baah.
   construct_from_asc(SelAtom, 0.01, 1.64, coot::COLOUR_BY_ATOM_TYPE, 0);
   verbose_reporting = 0;
   udd_has_ca_handle = -1;
}

Bond_lines_container::Bond_lines_container(atom_selection_container_t SelAtom,
					   float max_dist) {

   verbose_reporting = 0;
   do_disulfide_bonds_flag = 1;
   udd_has_ca_handle = -1;
   do_bonds_to_hydrogens = 1;
   b_factor_scale = 1.0;
   construct_from_asc(SelAtom, 0.01, max_dist, coot::COLOUR_BY_ATOM_TYPE, 0); 
}

Bond_lines_container::Bond_lines_container(atom_selection_container_t SelAtom,
					   float min_dist, float max_dist) {

   verbose_reporting = 0;
   do_disulfide_bonds_flag = 1;
   udd_has_ca_handle = -1;
   do_bonds_to_hydrogens = 1;
   b_factor_scale = 1.0;
   construct_from_asc(SelAtom, min_dist, max_dist, coot::COLOUR_BY_ATOM_TYPE, 0); 
}

// This is the one for occupancy and B-factor representation
// 
Bond_lines_container::Bond_lines_container (const atom_selection_container_t &SelAtom,
					    Bond_lines_container::bond_representation_type by_occ) {

   verbose_reporting = 0;
   do_disulfide_bonds_flag = 1;
   udd_has_ca_handle = -1;
   do_bonds_to_hydrogens = 1;
   b_factor_scale = 1.0;
   float max_dist = 1.64;
   if (by_occ == Bond_lines_container::COLOUR_BY_OCCUPANCY) {
      construct_from_asc(SelAtom, 0.01, max_dist, coot::COLOUR_BY_OCCUPANCY, 0); 
   } else {
      if (by_occ == Bond_lines_container::COLOUR_BY_B_FACTOR) {
	 try_set_b_factor_scale(SelAtom.mol);
	 construct_from_asc(SelAtom, 0.01, max_dist, coot::COLOUR_BY_B_FACTOR, 0);
      }
   }
}

void
Bond_lines_container::try_set_b_factor_scale(CMMDBManager *mol) {

   int udd_b_factor_handle =  mol->GetUDDHandle(UDR_HIERARCHY,
						(char *) coot::b_factor_bonds_scale_handle_name.c_str());
   // std::cout << "debug:: Got b factor udd handle: " << udd_b_factor_handle << std::endl;
   if (udd_b_factor_handle > 0) {
      realtype scale;
      if (mol->GetUDData(udd_b_factor_handle, scale) == UDDATA_Ok) {
	 b_factor_scale = scale;
      }
   }
}

void
Bond_lines_container::construct_from_atom_selection(const atom_selection_container_t &asc,
						    const PPCAtom atom_selection_1,
						    int n_selected_atoms_1,
						    const PPCAtom atom_selection_2,
						    int n_selected_atoms_2,
						    float min_dist, float max_dist,
						    int atom_colour_type,
						    short int are_different_atom_selections,
						    short int have_udd_atoms,
						    int udd_handle) {

   PSContact contact = NULL;
   int ncontacts;
   long i_contact_group = 1;

   // matrix stuff
   mat44 my_matt;
   CSymOps symm;

   // update my_matt;  You can't do this if you haven't set the space group.
   // 
   // symm.GetTMatrix(my_matt, 0);

   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;
      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;


   // cout << "my_matt is: " << my_matt << endl; // Argh! doesnt work.
   //
//    cout << "my_matt is: " << endl
// 	<< my_matt[0][0] << " "  << my_matt[0][1] << " "
// 	<< my_matt[0][2] << " "  << my_matt[0][3] << " "  << endl
// 	<< my_matt[1][0] << " "  << my_matt[1][1] << " "
// 	<< my_matt[1][2] << " "  << my_matt[1][3] << " "  << endl
// 	<< my_matt[2][0] << " "  << my_matt[2][1] << " "
// 	<< my_matt[2][2] << " "  << my_matt[2][3] << " "  << endl
// 	<< my_matt[3][0] << " "  << my_matt[3][1] << " "
// 	<< my_matt[3][2] << " "  << my_matt[3][3] << " "  << endl; 

//    std::cout << "Here are the atoms for which we will seek contacts\n";
//    for (int ii=0; ii<n_selected_atoms_1; ii++) { 
//       std::cout << atom_selection_1[ii] << std::endl;
//    } 

   // Note, it has happened a couple of times now, when we get a crash
   // in mmdb's MakeBricks (or something like that) from here, that's
   // because we are passing an atom that has a nan for a coordinate.

   if (0) { 
      std::cout << "Seeking contact: selection 1 " << std::endl;
      for (int ii=0; ii<n_selected_atoms_1; ii++)
	 std::cout << "   " << ii << " " << atom_selection_1[ii] << " :"
		   << atom_selection_1[ii]->isTer() << ":" << std::endl;
      std::cout << "Seeking contact: selection 2 " << std::endl;
      for (int ii=0; ii<n_selected_atoms_2; ii++)
	 std::cout << "   " << ii << " " << atom_selection_2[ii] << " :"
		   << atom_selection_2[ii]->isTer() << ":" << std::endl;
   }

   asc.mol->SeekContacts(atom_selection_1, n_selected_atoms_1,
			 atom_selection_2, n_selected_atoms_2,
			 min_dist, max_dist, // min, max distances
			 0,        // seqDist 0 -> in same res also
			 contact, ncontacts,
			 0, &my_matt, i_contact_group);

   if (0) {
      for (int i=0; i< ncontacts; i++) {
	 CAtom *atom_p_1 = atom_selection_1[ contact[i].id1 ];
	 CAtom *atom_p_2 = atom_selection_2[ contact[i].id2 ];
	 std::cout << "raw contact " << atom_p_1 << " to " << atom_p_2 << std::endl;
      }
   }

   //     if (verbose_reporting)
//    std::cout << "in construct_from_atom_selection found " << ncontacts << " contacts from "
//  	     << n_selected_atoms_1 << " vs " << n_selected_atoms_2
//  	     << " selected atoms. " <<  std::endl;

   std::string element_1, element_2;
   int col; // atom colour

   if (ncontacts > 0) {
      
      for (int i=0; i< ncontacts; i++) {

 	 if (are_different_atom_selections ||
 	      (contact[i].id2 > contact[i].id1) ) {

	    // 	 if (1) {

// 	    found_contact[contact[i].id1] = 1;  // true
// 	    found_contact[contact[i].id2] = 1;

	    CAtom *atom_p_1 = atom_selection_1[ contact[i].id1 ];
	    CAtom *atom_p_2 = atom_selection_2[ contact[i].id2 ];

	    std::string chain_id1(atom_p_1->GetChainID());
	    std::string chain_id2(atom_p_2->GetChainID());

	    std::string aloc_1(atom_p_1->altLoc);
	    std::string aloc_2(atom_p_2->altLoc);

	    element_1 = atom_p_1->element;
	    element_2 = atom_p_2->element;

	    coot::Cartesian atom_1(atom_p_1->x, atom_p_1->y, atom_p_1->z);
	    coot::Cartesian atom_2(atom_p_2->x, atom_p_2->y, atom_p_2->z);

	    if (chain_id1 == chain_id2) {

	       // alternate location test
	       // 
	       // 
	       // 	    short int v =  (aloc_1=="") || (aloc_2=="") || (aloc_1==aloc_2);
	       // 	    std::cout << "alt loc test result is " << v << "  comparing :" 
	       // 		      << aloc_1 << ": and :"<< aloc_2 << ":" << std::endl;
	       if ( (aloc_1=="") || (aloc_2=="") || (aloc_1==aloc_2) ) {

		  int res_1 = atom_p_1->GetSeqNum();
		  int res_2 = atom_p_2->GetSeqNum();

		  // this +/- 1 residue test
		  if (labs(res_1 - res_2) < 2 ||
		      labs(atom_p_1->residue->index - atom_p_2->residue->index) < 2) {

		     //  		  std::cout << "Adding bond " << atom_selection_1[ contact[i].id1 ]
		     //  			    << " to "
		     //  			    << atom_selection_2[ contact[i].id2 ] << std::endl;

		     if (atom_selection_1[ contact[i].id1 ]->GetModel() ==
			 atom_selection_2[ contact[i].id2 ]->GetModel()) {

			if (have_udd_atoms) {
 			   if (! ((!strcmp(atom_selection_1[ contact[i].id1 ]->element, " S")) ||
				  (!strcmp(atom_selection_1[ contact[i].id1 ]->element, "SE")) ||
				  (!strcmp(atom_selection_1[ contact[i].id1 ]->element, "CL")) ||
				  (!strcmp(atom_selection_1[ contact[i].id1 ]->element, "BR")) ||
				  (!strcmp(atom_selection_1[ contact[i].id1 ]->element, "Cl")) ||
				  (!strcmp(atom_selection_1[ contact[i].id1 ]->element, "Br")) ||
				  (!strcmp(atom_selection_1[ contact[i].id1 ]->element, " P")))) { 
			      atom_selection_1[ contact[i].id1 ]->PutUDData(udd_handle, 1);
// 			      std::cout << "marking udd handle "
// 					<< atom_selection_1[ contact[i].id1 ]
// 					<< std::endl;
			   } 

 			   if (! ((!strcmp(atom_selection_2[ contact[i].id2 ]->element, " S")) ||
				  (!strcmp(atom_selection_2[ contact[i].id2 ]->element, "SE")) ||
				  (!strcmp(atom_selection_1[ contact[i].id1 ]->element, "CL")) ||
				  (!strcmp(atom_selection_1[ contact[i].id1 ]->element, "BR")) ||
				  (!strcmp(atom_selection_1[ contact[i].id1 ]->element, "Cl")) ||
				  (!strcmp(atom_selection_1[ contact[i].id1 ]->element, "Br")) ||
				  (!strcmp(atom_selection_2[ contact[i].id2 ]->element, " P")))) { 
			      atom_selection_2[ contact[i].id2 ]->PutUDData(udd_handle, 1);
// 			      std::cout << "marking udd handle "
// 					<< atom_selection_2[ contact[i].id2 ]
// 					<< std::endl;
			   }
			}

			if (element_1 != element_2) {
		  
			   // Bonded to different atom elements.
			   //

			   if ((element_1 != " H") && (element_2 != " H")) {

			      add_half_bonds(atom_1, atom_2,
					     atom_selection_1[contact[i].id1],
					     atom_selection_2[contact[i].id2],
					     atom_colour_type);
			   } else {
			      
			      // Bonds to hydrogens are one colour -
			      // HYDROGEN_GREY_BOND, not half-bonds.
			      // 
			      // Except hydrogens on waters are
			      // treated differently to other
			      // hydrogens (if they are not then we
			      // don't get to see the oxygen).
			      
			      std::string resname_1 = atom_p_1->GetResName();
			      std::string resname_2 = atom_p_2->GetResName();
			      if (resname_1 == "HOH" || resname_2 == "HOH") {
				 add_half_bonds(atom_1, atom_2,
						atom_selection_1[contact[i].id1],
						atom_selection_2[contact[i].id2],
						atom_colour_type);
			      } else { 
				 
				 addBond(HYDROGEN_GREY_BOND, atom_1, atom_2);
			      }
			} 
		  
			} else {
		  
			   // Bonded to same an atom of the same element.
			   //
			   col = atom_colour(atom_selection_1[ contact[i].id1 ], atom_colour_type);
			   addBond(col, atom_1, atom_2);
			}
		     }
		  }
	       }
	    }
	 } // contact atom is higher up the list check.
// 	 else {
// 	    std::cout << "debug:: ignoring contact " << i << std::endl;
// 	 }
      } // i over ncontacts
      
      delete [] contact;
   }
}


void
Bond_lines_container::add_half_bonds(const coot::Cartesian &atom_1,
				     const coot::Cartesian &atom_2,
				     CAtom *at_1,
				     CAtom *at_2,
				     int atom_colour_type) {
   
   coot::Cartesian bond_mid_point = atom_1.mid_point(atom_2);
   int col = atom_colour(at_1, atom_colour_type);
   addBond(col, atom_1, bond_mid_point);
   
   col = atom_colour(at_2, atom_colour_type);
   addBond(col, bond_mid_point, atom_2);
			      
}


void
Bond_lines_container::construct_from_model_links(CModel *model_p,
						 int atom_colour_type) {

   // Interestingly, when we add a LINK to a PDB file, if there are
   // less residues in a chain than is specified in a LINK line, mmdb
   // expands the residue list in a chain with NULL residues!
   
   int n_links = model_p->GetNumberOfLinks();
   if (n_links > 0) { 
      for (int i_link=1; i_link<=n_links; i_link++) {
	 PCLink link = model_p->GetLink(i_link);
	 PCAtom atom_1 = NULL;
	 PCAtom atom_2 = NULL;
	 int n_chains = model_p->GetNumberOfChains();
	 for (int ich=0; ich<n_chains; ich++) {
	    CChain *chain_p = model_p->GetChain(ich);
	    if (chain_p) { 
	       if (std::string(chain_p->GetChainID()) ==
		   std::string(link->chainID1)) {
		  int n_residues = model_p->GetNumberOfResidues();
		  for (int i_res=0; i_res<n_residues; i_res++) {
		     CResidue *res_p = chain_p->GetResidue(i_res);
		     if (res_p) {
			if (res_p->GetSeqNum() == link->seqNum1) {
			   if (std::string(res_p->GetInsCode()) ==
			       std::string(link->insCode1)) {
			      int n_atoms = res_p->GetNumberOfAtoms();
			      for (int iat=0; iat<n_atoms; iat++) {
				 CAtom *at = res_p->GetAtom(iat);
				 if (! at->isTer()) { 
				    if (std::string(at->name) ==
					std::string(link->atName1)) {
				       if (std::string(at->altLoc) ==
					   std::string(link->aloc1)) {
					  atom_1 = at;
					  break;
				       }
				    }
				 }
				 if (atom_1) break;
			      }
			   }
			}
		     } // null residue test
		     if (atom_1) break;
		  }
	       }
	    } // chain_p test
	    if (atom_1) break;
	 }

	 if (atom_1) {
	    for (int ich=0; ich<n_chains; ich++) {
	       CChain *chain_p = model_p->GetChain(ich);
	       if (chain_p) { 
		  if (std::string(chain_p->GetChainID()) ==
		      std::string(link->chainID2)) {
		     int n_residues = model_p->GetNumberOfResidues();
		     for (int i_res=0; i_res<n_residues; i_res++) {
			CResidue *res_p = chain_p->GetResidue(i_res);
			if (res_p) { 
			   if (res_p->GetSeqNum() == link->seqNum2) {
			      if (std::string(res_p->GetInsCode()) ==
				  std::string(link->insCode2)) {
				 int n_atoms = res_p->GetNumberOfAtoms();
				 for (int iat=0; iat<n_atoms; iat++) {
				    CAtom *at = res_p->GetAtom(iat);
				    if (! at->isTer()) { 
				       if (std::string(at->name) ==
					   std::string(link->atName2)) {
					  if (std::string(at->altLoc) ==
					      std::string(link->aloc2)) {
					     atom_2 = at;
					     break;
					  }
				       }
				    }
				    if (atom_2) break;
				 }
			      }
			   }
			} // res_p test
			if (atom_2) break;
		     }
		  }
		  if (atom_2) break;
	       } // chain_p test
	    }
	 } 

	 // OK, make the link bond then!
	 if (atom_1 && atom_2) {
	    coot::Cartesian pos_1(atom_1->x, atom_1->y, atom_1->z);
	    coot::Cartesian pos_2(atom_2->x, atom_2->y, atom_2->z);
	    std::string ele_1 = atom_1->element;
	    std::string ele_2 = atom_2->element;
	    if (ele_1 == ele_2) {
	       int col = atom_colour(atom_1, atom_colour_type);
	       add_dashed_bond(col, pos_1, pos_2, 0);
	    } else {
	       coot::Cartesian bond_mid_point = pos_1.mid_point(pos_2);
	       int col = atom_colour(atom_1, atom_colour_type);
	       add_dashed_bond(col, pos_1, bond_mid_point, 1);
	       col = atom_colour(atom_2, atom_colour_type);
	       add_dashed_bond(col, bond_mid_point, pos_2, 1);
	    } 
	 }
      }
   }
}


PPCAtom
coot::model_bond_atom_info_t::Hydrogen_atoms() const {

   PPCAtom H_atoms = new PCAtom[hydrogen_atoms_.size()];
   for (int i=0; i<hydrogen_atoms_.size(); i++) {
      H_atoms[i] = hydrogen_atoms_[i];
   }
   return H_atoms;
}

PPCAtom
coot::model_bond_atom_info_t::non_Hydrogen_atoms() const {

   PPCAtom non_H_atoms = new PCAtom[non_hydrogen_atoms_.size()];
   for (int i=0; i<non_hydrogen_atoms_.size(); i++) {
      non_H_atoms[i] = non_hydrogen_atoms_[i];
   }
   return non_H_atoms;
}

// we rely on SelAtom.atom_selection being properly constucted to
// contain all atoms
void
Bond_lines_container::construct_from_asc(const atom_selection_container_t &SelAtom, 
					 float min_dist, float max_dist,
					 int atom_colour_type,
					 short int is_from_symmetry_flag) {

   // initialize each colour in the Bond_lines_container
   //
   if (bonds.size() == 0) { 
      for (int i=0; i<10; i++) { 
	 Bond_lines a(i);
	 bonds.push_back(a);
      }
   }
   float star_size = 0.4;

   // initialize the hydrogen bonding flag:
   //
   // do_bonds_to_hydrogens = 1;  // Fix off, is this sensible?  Valgrind
                               // // complained (rightly) about jump
                               // 80-90 lines below.
                               // (do_bonds_to_hydrogens had be unitialized).
   // 
                               // 20060812: Actually, I now doubt that
                               // this is uninitialized.  The value is
                               // set from/in the constructor from the
                               // molecule_class_info_t value
                               // draw_hydrogens_flag, which is set in
                               // the molecule_class_info_t
                               // constructors.  So I don't know what
                               // is going on... I will comment out
                               // this line now and see if valgrind
                               // still complains.
                               // 
                               // OK, I just checked. Valgrind no
                               // longer complains if I comment this line.
   // 20070407 Back here again! Valgrind *does* complain about
   // uninitialized do_bonds_to_hydrogens when I enable NCS ghosts.
   
   
   // and disulfides bond flag:
   do_disulfide_bonds_flag = 1;


   if (SelAtom.n_selected_atoms <= 0) return;

   // count the number of hydrogen atoms and non-hydrogen atoms:
   //
//    std::string element;
//    int n_H = 0;
//    int n_non_H = 0;
//    for (int i=0; i<SelAtom.n_selected_atoms; i++) {
//       element = SelAtom.atom_selection[i]->element;
//       if (element == " H" || element == " D") {
// 	 n_H++;
//       } else {
// 	 n_non_H++;
//       }
//    }

   // Now, let's not forget that some atoms don't have contacts, so
   // whenever we find a contact for an atom, we mark it with
   // UserDefinedData "found bond".
   // 

   int uddHnd = SelAtom.mol->RegisterUDInteger(UDR_ATOM, "found bond");
   short int have_udd_atoms = 1;
   if (uddHnd<0)  {
      std::cout << " atom bonding registration failed.\n";
      have_udd_atoms = 0;
   } else {
      for (int i=0; i<SelAtom.n_selected_atoms; i++) { 
	 SelAtom.atom_selection[i]->PutUDData(uddHnd,0);
      } 
   }

   // ------------------
   //   FIXME
   // ------------------
   // 
   // We need to loop over each model.  Currently the SelAtom is not
   // selected on model.
   //
   // consider a vector (over all models) of these:
   // 
   // class model_bond_atom_info_t {
   //       std::vector<PCAtom> hydrogen_atoms_;
   //       std::vector<PCAtom> non_hydrogen_atoms_;
   //    public:
   //    PPCAtom     Hydrogen_atoms();
   //    PPCAtom non_Hydrogen_atoms();
   //    int n_H();
   //    int n_non_H();
   // }
   
   int imodel; 
   int n_models = SelAtom.mol->GetNumberOfModels(); // models start at number 1
   std::vector<coot::model_bond_atom_info_t> atom_stuff_vec(n_models+1);
   
   for (int i=0; i<SelAtom.n_selected_atoms; i++) {
      imodel = SelAtom.atom_selection[i]->GetModelNum();
      if ((imodel <= n_models) && (imodel > 0)) { 
	 atom_stuff_vec[imodel].add_atom(SelAtom.atom_selection[i]);
      }
   }



   for (int imodel=1; imodel<=n_models; imodel++) { 

      PPCAtom     Hydrogen_atoms = atom_stuff_vec[imodel].Hydrogen_atoms();
      PPCAtom non_Hydrogen_atoms = atom_stuff_vec[imodel].non_Hydrogen_atoms();
      int n_non_H = atom_stuff_vec[imodel].n_non_H();
      int n_H = atom_stuff_vec[imodel].n_H();

      construct_from_atom_selection(SelAtom,
				    non_Hydrogen_atoms, n_non_H,
				    non_Hydrogen_atoms, n_non_H,
				    min_dist, max_dist, atom_colour_type,
				    0, have_udd_atoms, uddHnd);

      if (do_bonds_to_hydrogens && (n_H > 0)) {

	 float H_min_dist = 0.7;
	 float H_max_dist = 1.42;

	 // H-H
	 construct_from_atom_selection(SelAtom,
				       Hydrogen_atoms, n_H,
				       Hydrogen_atoms, n_H,
				       H_min_dist, H_max_dist, atom_colour_type,
				       0, have_udd_atoms, uddHnd);

	 // H-X
	 construct_from_atom_selection(SelAtom,
				       non_Hydrogen_atoms, n_non_H,
				       Hydrogen_atoms, n_H,
				       H_min_dist, H_max_dist, atom_colour_type,
				       1, have_udd_atoms, uddHnd);
      }

      construct_from_model_links(SelAtom.mol->GetModel(imodel), atom_colour_type);
      
      // std::cout << "DEBUG:: (post) SelAtom: mol, n_selected_atoms "
      // << SelAtom.mol << " " << SelAtom.n_selected_atoms << std::endl;

      if (do_disulfide_bonds_flag)
	 do_disulphide_bonds(SelAtom, imodel);

      // now we have dealt with the bonded atoms, lets find the non-bonded
      // atoms...

      if (have_udd_atoms) {
   
	 // for atoms with no neighbour (contacts):
	 coot::Cartesian small_vec_x(star_size, 0.0, 0.0);
	 coot::Cartesian small_vec_y(0.0, star_size, 0.0);
	 coot::Cartesian small_vec_z(0.0, 0.0, star_size);

	 int ic; // changed by reference (UDData)
	 int col;

	 for (int i=0; i<n_non_H; i++) {
      
	    if (non_Hydrogen_atoms[i]->GetUDData(uddHnd, ic) == UDDATA_Ok) {
	       if ((ic == 0) ||
		   (!strcmp(non_Hydrogen_atoms[i]->element, " S")) ||
		   (!strcmp(non_Hydrogen_atoms[i]->element, "SE")) ||
		   (!strcmp(non_Hydrogen_atoms[i]->element, " P"))) {

// 		  std::cout << " No contact for " << non_Hydrogen_atoms[i]
// 			    << std::endl;
	       
		  // no contact found or was Sulphur, or Phosphor

		  // So, was this a seleno-methione?
		  //
		  CResidue *atom_residue_p = non_Hydrogen_atoms[i]->residue;
		  if (atom_residue_p) {
		     std::string resname = non_Hydrogen_atoms[i]->GetResName();
		     if ((is_from_symmetry_flag == 0) &&
			 (resname == "MSE" || resname == "MET" || resname == "MSO"
			  || resname == "CYS" )) {
			handle_MET_or_MSE_case(non_Hydrogen_atoms[i], uddHnd, atom_colour_type);
		     } else {
			std::string ele = non_Hydrogen_atoms[i]->element;
			if (ele == "CL" || ele == "BR" || ele == " S" ||  ele == " I"
			    || ele == "Cl" || ele == "Br" 
			    || ele == "PT" 
			    || ele == "AS" || ele == " P" || ele == "AU" || ele == "HG") {
			   handle_long_bonded_atom(non_Hydrogen_atoms[i], uddHnd, atom_colour_type);
			}
		     }
		  } else {
		     std::cout << "INFO:: trapped atom without residue in non-bonded atom check: "
			       << non_Hydrogen_atoms[i] << std::endl;
		  }
	       }
	    } else {
	       std::cout << "missing UDData for atom "
			 << non_Hydrogen_atoms[i] << "\n";
	    }
	 }


	 // Make the stars...
	 // 
	 for (int i=0; i<n_non_H; i++) {
	    if (non_Hydrogen_atoms[i]->GetUDData(uddHnd, ic) == UDDATA_Ok) {
	       if (ic == 0) {
		  // no contact found
		  col = atom_colour(non_Hydrogen_atoms[i], atom_colour_type);
		  coot::Cartesian atom(non_Hydrogen_atoms[i]->x,
				       non_Hydrogen_atoms[i]->y,
				       non_Hydrogen_atoms[i]->z);
		  
		  addBond(col, atom+small_vec_x, atom-small_vec_x);
		  addBond(col, atom+small_vec_y, atom-small_vec_y);
		  addBond(col, atom+small_vec_z, atom-small_vec_z);
	       }
	    }
	 }


	 if (do_bonds_to_hydrogens && (n_H > 0)) {
	    for (int i=0; i<n_H; i++) {
	       if (Hydrogen_atoms[i]->GetUDData(uddHnd, ic) == UDDATA_Ok) {
		  if (ic == 0) {
	       
		     // no contact found
		     col = atom_colour(Hydrogen_atoms[i], atom_colour_type);
		     coot::Cartesian atom(Hydrogen_atoms[i]->x,
					  Hydrogen_atoms[i]->y,
					  Hydrogen_atoms[i]->z);
	       
		     addBond(col, atom+small_vec_x, atom-small_vec_x);
		     addBond(col, atom+small_vec_y, atom-small_vec_y);
		     addBond(col, atom+small_vec_z, atom-small_vec_z);
		  }
	       }
	    }
	 }
      }

      delete [] Hydrogen_atoms;
      delete [] non_Hydrogen_atoms;
   }

   add_zero_occ_spots(SelAtom);
   add_atom_centres(SelAtom, atom_colour_type);
}

void
Bond_lines_container::handle_MET_or_MSE_case(PCAtom mse_atom,
					     int udd_handle, 
					     int atom_colour_type,
					     coot::my_atom_colour_map_t *atom_colour_map_p) {

   // std::cout << "Handling MET/MSE case for atom " << mse_atom << std::endl;
   
   std::string atom_name(mse_atom->name);
   std::string residue_name(mse_atom->GetResName());
   if (residue_name == "MET" || residue_name == "MSE" || residue_name == "MSO") { 
      if (atom_name == "SE  " || atom_name == " SD ") {
	 int col = atom_colour(mse_atom, atom_colour_type, atom_colour_map_p);

	 // We need to add special bonds SE -> CE and SE -> CG.
	 PPCAtom residue_atoms;
	 int nResidueAtoms;
	 mse_atom->residue->GetAtomTable(residue_atoms, nResidueAtoms);
	 for (int i=0; i<nResidueAtoms; i++) {
	    std::string table_atom_name(residue_atoms[i]->name);
	    if (table_atom_name == " CG " ||
		table_atom_name == " CE " ) {
	       // mse_atom or met_atom now of course.
	       coot::Cartesian cart_at1(mse_atom->x, mse_atom->y, mse_atom->z);
	       coot::Cartesian cart_at2(residue_atoms[i]->x,
					residue_atoms[i]->y,
					residue_atoms[i]->z);
	    
	       std::string altconf1 = mse_atom->altLoc;
	       std::string altconf2 = residue_atoms[i]->altLoc;
	       if ( (altconf1=="") || (altconf2=="") || (altconf1==altconf2) ) {
		  coot::Cartesian bond_mid_point = cart_at1.mid_point(cart_at2);
		  int colc = atom_colour(residue_atoms[i], atom_colour_type, atom_colour_map_p);
		  // int colc = atom_colour_type; // just to check

		  float bond_length = (cart_at1 - cart_at2).amplitude();
		  if (bond_length < 3.0) { // surely not longer than this?
		     addBond(col,  cart_at1, bond_mid_point); // The SE->mid_pt bond
		     addBond(colc, bond_mid_point, cart_at2); // The Cx->mid_pt bond
		     // mark atom as bonded.
		     residue_atoms[i]->PutUDData(udd_handle, 1);
		     mse_atom->PutUDData(udd_handle, 1);
		  }
	       }
	    }
	 }
      }
   }
   if (residue_name == "CYS") {
      int col = atom_colour(mse_atom, atom_colour_type, atom_colour_map_p);
      // We need to add special bonds SE -> CE and SE -> CG.
      PPCAtom residue_atoms;
      int nResidueAtoms;
      mse_atom->residue->GetAtomTable(residue_atoms, nResidueAtoms);
      for (int i=0; i<nResidueAtoms; i++) {
	 std::string table_atom_name(residue_atoms[i]->name);
	 if (table_atom_name == " CB ") {
	    coot::Cartesian cart_at1(mse_atom->x, mse_atom->y, mse_atom->z);
	    coot::Cartesian cart_at2(residue_atoms[i]->x,
				     residue_atoms[i]->y,
				     residue_atoms[i]->z);

	    std::string altconf1 = mse_atom->altLoc;
	    std::string altconf2 = residue_atoms[i]->altLoc;
	    if ( (altconf1=="") || (altconf2=="") || (altconf1==altconf2) ) {
	       coot::Cartesian bond_mid_point = cart_at1.mid_point(cart_at2);
	       int colc = atom_colour(residue_atoms[i], atom_colour_type, atom_colour_map_p);
	       addBond(col,  cart_at1, bond_mid_point);
	       addBond(colc, bond_mid_point, cart_at2);
	       // mark atom as bonded.
	       residue_atoms[i]->PutUDData(udd_handle, 1);
	       mse_atom->PutUDData(udd_handle, 1);
	    }
	 }
      }
   }
}

void
Bond_lines_container::handle_long_bonded_atom(PCAtom atom,
					      int udd_handle,
					      int atom_colour_type,
					      coot::my_atom_colour_map_t *atom_colour_map_p) {

   float bond_limit = 2.16; // A S-S bonds are 2.05A.  So we've added
			    // some wiggle room (2.1 was too short for
			    // some dictionary S-S).
   
   std::string atom_name(atom->name);
   std::string residue_name(atom->GetResName());
   std::string element(atom->element);
   CResidue *res = atom->residue;

   // std::cout << "handling long bonds for " << atom << std::endl;
   // std::cout << " ele name :" << element << ":" << std::endl;
   
   if (atom_name == "AS  ")
      bond_limit = 2.4;
   if (element == "AU")
      bond_limit = 2.4;
   if (element == "AS")
      bond_limit = 2.4;
   if (element == "HG")
      bond_limit = 2.4;
   if (element == " I")
      bond_limit = 2.3; // C-I is 2.13 according to wikipedia
   float  bl2 = bond_limit * bond_limit;
   float h_bl2 = 1.8 * 1.8; // 20100208 all bonds to hydrogens are less than 1.8A ? (guess)
   short int bond_added_flag = 0;

   if (res) {
      // do the bonding by hand:
      coot::Cartesian atom_pos(atom->x, atom->y, atom->z);
      int col = atom_colour(atom, atom_colour_type, atom_colour_map_p);
      PPCAtom residue_atoms = 0;
      int nResidueAtoms;
      res->GetAtomTable(residue_atoms, nResidueAtoms);
      for (int i=0; i<nResidueAtoms; i++) {
	 if (residue_atoms[i] != atom) {
	    coot::Cartesian res_atom_pos(residue_atoms[i]->x,
					 residue_atoms[i]->y,
					 residue_atoms[i]->z);

	    // We compared squard bond distances (so that we don't
	    // have to take the square root for everything of course).
	    //
	    std::string res_atom_ele = residue_atoms[i]->element;
	    float len2 = (atom_pos - res_atom_pos).amplitude_squared();
	    if (((len2 <   bl2) && (res_atom_ele != " H")) ||
		((len2 < h_bl2) && (res_atom_ele == " H"))) {
	       std::string altconf1 = atom->altLoc;
	       std::string altconf2 = residue_atoms[i]->altLoc;
	       if ( (altconf1=="") || (altconf2=="") || (altconf1==altconf2) ) {
		  coot::Cartesian bond_mid_point = atom_pos.mid_point(res_atom_pos);
		  int colc = atom_colour(residue_atoms[i], atom_colour_type, atom_colour_map_p);
		  addBond(col,  atom_pos, bond_mid_point);
		  addBond(colc, bond_mid_point, res_atom_pos);
		  bond_added_flag = 1;
		  // mark the atom as bonded.
		  atom->PutUDData(udd_handle, 1);
		  residue_atoms[i]->PutUDData(udd_handle, 1);
	       }
	    }
	 }
      }
   }
   
   if (!bond_added_flag) {
      // bond it like a single atom then:
      
      float star_size = 0.4;
      coot::Cartesian small_vec_x(star_size, 0.0, 0.0);
      coot::Cartesian small_vec_y(0.0, star_size, 0.0);
      coot::Cartesian small_vec_z(0.0, 0.0, star_size);

      int col = atom_colour(atom, atom_colour_type, atom_colour_map_p);
      coot::Cartesian atom_pos(atom->x, atom->y, atom->z);
      addBond(col, atom_pos+small_vec_x, atom_pos-small_vec_x);
      addBond(col, atom_pos+small_vec_y, atom_pos-small_vec_y);
      addBond(col, atom_pos+small_vec_z, atom_pos-small_vec_z);
   }
}

 



// This finds bonds between a residue and the protein (in SelAtom).
// It is used for the environment bonds box.
// 
// What do we do about drawing hydrogens?
// 
Bond_lines_container::Bond_lines_container(const atom_selection_container_t &SelAtom,
					   PPCAtom residue_atoms,
					   int n_residue_atoms,
					   coot::protein_geometry *protein_geom_p, // modifiable
					   short int residue_is_water_flag,
					   bool draw_env_distances_to_hydrogens_flag,
					   float min_dist,
					   float max_dist) {

   if (0) 
      std::cout << "Environment distances NO symm" << std::endl;
   do_bonds_to_hydrogens = 1;  // added 20070629
   
   b_factor_scale = 1.0;
   int ncontacts;
   PSContact contact = NULL;
   // initialize each colour in the Bond_lines_container
   //
   if (bonds.size() == 0) { 
      for (int i=0; i<10; i++) { 
	 Bond_lines a(i);
	 bonds.push_back(a);
      }
   }
   

   // seqDist has to be 0 here or else failure to find contacts to new
   // waters added as pointer atoms. Why?  Don't know - a bug in mmdb,
   // I suspect.
   //
   SelAtom.mol->SeekContacts(residue_atoms, n_residue_atoms,
			     SelAtom.atom_selection, SelAtom.n_selected_atoms,
			     min_dist,
			     max_dist,
			     0,  // seqDist (in same residue allowed)
			     contact, ncontacts);

   if (0) {  // debugging seqDist
      std::cout << " DEBUG:: there are " << n_residue_atoms << " residue atoms "
		<< " and " << SelAtom.n_selected_atoms << " mol atoms\n";
      for (int iat=0; iat<n_residue_atoms; iat++)
	 std::cout << "residue atoms: " << iat << "/" << n_residue_atoms
		   << " " << residue_atoms[iat] << std::endl;
      for (int iat=0; iat<SelAtom.n_selected_atoms; iat++)
	 std::cout << "    mol atoms: " << iat << "/" << SelAtom.n_selected_atoms
		   << " " << SelAtom.atom_selection[iat] << std::endl;
   }

   
   if (ncontacts > 0) {
      for (int i=0; i<ncontacts; i++) {
	 if ( draw_these_residue_contacts(residue_atoms[contact[i].id1]->GetResidue(),
					  SelAtom.atom_selection[contact[i].id2]->GetResidue(),
					  protein_geom_p)
	      || residue_is_water_flag) {
	    
	    coot::Cartesian atom_1(residue_atoms[ contact[i].id1 ]->x,
				   residue_atoms[ contact[i].id1 ]->y,
				   residue_atoms[ contact[i].id1 ]->z);
	    coot::Cartesian atom_2(SelAtom.atom_selection[ contact[i].id2 ]->x,
				   SelAtom.atom_selection[ contact[i].id2 ]->y,
				   SelAtom.atom_selection[ contact[i].id2 ]->z);
	    std::string ele1 = residue_atoms[ contact[i].id1 ]->element;
	    std::string ele2 = SelAtom.atom_selection[ contact[i].id2 ]->element;

	    if (0) { // debug
	       std::cout << " DEBUG:: add env dist "
			 << residue_atoms[ contact[i].id1 ] << " to "
			 << SelAtom.atom_selection[ contact[i].id2 ]
			 << std::endl;
	    }
	    if (draw_env_distances_to_hydrogens_flag ||
		((ele1 != " H") && (ele2 != " H"))) { 
	       if (ele1 == " C")
		  addBond(0, atom_1, atom_2);
	       else {
		  if (ele2 == " C") { 
		     addBond(0, atom_1, atom_2);
		  } else { 
		     if (ele1 == " H" && ele2 == " H") { 
			addBond(0, atom_1, atom_2);
		     } else { 
			addBond(1, atom_1, atom_2);
		     }
		  }
	       }
	    }
	 } 
      }
      delete [] contact;
   }
}

// Shall we draw environment bonds between this_residue and env_residue?
//
// No, if the residues are next to each other in sequence in the same
// chain and are a polymer.
//     
// Otherwise, yes.
// 
bool
Bond_lines_container::draw_these_residue_contacts(CResidue *this_residue,
						  CResidue *env_residue,
						  coot::protein_geometry *protein_geom_p // modifiable
						  ) {

   if (this_residue != env_residue) { 
      std::string ch1(this_residue->GetChainID());
      std::string ch2(env_residue->GetChainID());
      if ((abs(this_residue->GetSeqNum() - env_residue->GetSeqNum()) > 1)
	  || (ch1 != ch2)) {
	 return 1;
      } else {
	 // are we in a polymer? if so, no draw.
	 //
	 std::string this_res_type = this_residue->GetResName();
	 std::string env_residue_res_type = env_residue->GetResName();
	 if (protein_geom_p->linkable_residue_types_p(this_res_type, env_residue_res_type)) {
	    return 0;
	 } else {
	    return 1;
	 }
      }
   } else {
      return 0;
   } 
}



// This finds bonds between a residue and the protein (in SelAtom).
// It is used for the environment bonds box.
//
// It seems that this will draw bond to hydrogens, even if they are
// not being displayed.  I guess that we need to pass
// bonds_to_hydrogens flag.  And above.
// 
Bond_lines_container::Bond_lines_container(const atom_selection_container_t &SelAtom,
					   PPCAtom residue_atoms,
					   int nResidueAtoms,
					   float min_dist,
					   float max_dist,
					   bool draw_env_distances_to_hydrogens_flag,
					   short int do_symm) {

   // std::cout << "Environment distances with symm" << std::endl;

   do_bonds_to_hydrogens = 1;  // added 20070629

   b_factor_scale = 1.0;
   if (bonds.size() == 0) { 
      for (int i=0; i<10; i++) { 
	 Bond_lines a(i);
	 bonds.push_back(a);
      }
   }
   
   if (min_dist> max_dist) {
      float tmp = max_dist;
      max_dist = min_dist;
      min_dist = tmp;
   }
   
   for (int iresatom=0; iresatom< nResidueAtoms; iresatom++) { 
      coot::Cartesian res_atom_pos(residue_atoms[iresatom]->x,
				   residue_atoms[iresatom]->y,
				   residue_atoms[iresatom]->z);
      
      molecule_extents_t mol_extents(SelAtom, max_dist);
      std::vector<std::pair<symm_trans_t, Cell_Translation> > boxes =
	 mol_extents.which_boxes(res_atom_pos, SelAtom);
      
      for (unsigned int ibox=0; ibox<boxes.size(); ibox++) { 
	 PPCAtom translated = trans_sel(SelAtom, boxes[ibox]);
	 
	 for (int it=0; it<SelAtom.n_selected_atoms; it++) { 
	    
	    coot::Cartesian symm_atom(translated[it]->x,
				      translated[it]->y,
				      translated[it]->z);
	    
	    float d = coot::Cartesian::LineLength(symm_atom, res_atom_pos);
 	    // std::cout << translated[it] << " d = " << d << std::endl;
	    if (d < max_dist && d >= min_dist) {
	       std::string ele1 = residue_atoms[iresatom]->element;
	       std::string ele2 = translated[it]->element;

	       std::cout << "adding bond: " << residue_atoms[iresatom] << " -> "
			 << translated[it] << std::endl;

	       if (draw_env_distances_to_hydrogens_flag ||
		   ((ele1 != " H") && (ele2 != " H"))) { 
		  if (ele1 == " C")
		     addBond(0, symm_atom, res_atom_pos);
		  else {
		     if (ele2 == " C") { 
			addBond(0, symm_atom, res_atom_pos);
		     } else { 
			if (ele1 == " H" && ele2 == " H") { 
			   addBond(0, symm_atom, res_atom_pos);
			} else { 
			   addBond(1, symm_atom, res_atom_pos);
			}
		     }
		  }
	       }
	    }
	 }
	 // valgrind suggestion: 050803 - deleting translated[i] too.
	 // c.f. pick usage/deletion of translation atoms
	 for (int i=0; i<SelAtom.n_selected_atoms; i++)
	    delete translated[i];
	 delete [] translated;
      }
   }
   // std::cout << "returning from Bond_lines_container (symm env) constructor \n" ;
}

// This *looks* like a kludge - symm_trans is forced into a vector.
// But it is not a kludge. addSymmetry needs to take a vector, because
// it's needed for CA symm addition too.
// 
//
graphical_bonds_container
Bond_lines_container::addSymmetry(const atom_selection_container_t &SelAtom,
				  coot::Cartesian point,
				  float symm_distance,
				  const std::vector<std::pair<symm_trans_t, Cell_Translation> > &symm_trans,
				  short int symmetry_as_ca_flag,
				  short int symmetry_whole_chain_flag) {

   graphical_bonds_container gbc; 

   if (symmetry_as_ca_flag == 1) {
      gbc = addSymmetry_calphas(SelAtom, point, symm_distance, symm_trans);

   } else {

      if (symmetry_whole_chain_flag) {
	 gbc = addSymmetry_whole_chain(SelAtom, point,
				       symm_distance, symm_trans);
      } else {
	 
      
	 PSContact contact = NULL;
	 int ncontacts;

	 if (SelAtom.n_selected_atoms > 0) { 
	    PCAtom point_atom_p = new CAtom;
	    point_atom_p->SetCoordinates(point.get_x(), point.get_y(),
					 point.get_z(), 1.0, 99.9);

	    for(unsigned int ii=0; ii<symm_trans.size(); ii++) { // i.e. boxes

	       // 	 cout << "---------------------------- "
	       // 	      << symm_trans[ii]  
	       // 	      << " ---------------------------- "
	       // 	      << endl;
	 
	       // now get an atom selection where all the atoms are moved
	       // by symm_trans[ii]
	       //
	       // Create a selection trans_selection that is a copy of
	       // SelAtom.atom_selection that has had symmetry
	       // transformation symm_trans[ii] applied to it.
	       // 
	       PPCAtom trans_selection = trans_sel(SelAtom, symm_trans[ii]); // heavyweight

	       contact = NULL; // assigned next, deleted below.
	       SelAtom.mol->SeekContacts(point_atom_p, trans_selection,
					 SelAtom.n_selected_atoms,
					 0.0, symm_distance,
					 0,  // seqDist
					 contact, ncontacts);
	 
	       if (ncontacts > 0 ) {

		  // the atom_selection of Contact_Sel is allocated.
		  // We give it below.
		  // 
		  atom_selection_container_t Contact_Sel = 
		     ContactSel(trans_selection, contact, ncontacts); 
		  Contact_Sel.mol = SelAtom.mol;

		  construct_from_asc(Contact_Sel, 0.01, 1.95,
				     coot::COLOUR_BY_ATOM_TYPE, 1);
		  gbc = make_graphical_symmetry_bonds();

		  // Now give back the atom_selection, (but not the atoms
		  // because they are created by trans_sel and are given
		  // back later.     giveback_1
		  //
		  delete [] Contact_Sel.atom_selection;
	    
	       } else {

		  // do a constructor for no symmetry bonds
		  // and reset graphical_symmetry_bonds.

		  no_symmetry_bonds();

	       }

	       // Tidy up.
	       // 
	       if (trans_selection) { 
		  for (int j=0; j<SelAtom.n_selected_atoms; j++)
		     // delete each of the atoms
		     delete trans_selection[j];
		  delete [] trans_selection;
	       }

	       delete [] contact;
	    }
	    delete point_atom_p;
	 }
      }
   }
   return gbc; 
}

// FYI: there is only one element to symm_trans, the is called from
// them addSymmetry_vector_symms wrapper
graphical_bonds_container
Bond_lines_container::addSymmetry_whole_chain(const atom_selection_container_t &SelAtom,
					      const coot::Cartesian &point,
					      float symm_distance, 
					      const std::vector <std::pair<symm_trans_t, Cell_Translation> > &symm_trans) {

   graphical_bonds_container gbc;
   mat44 my_matt;
   mat44 identity_matt;
   //   PSContact contact = NULL;
   // int ncontacts;
   // long i_contact_group = 1;

   for (int im=0; im<4; im++) { 
      for (int jm=0; jm<4; jm++) {
	 identity_matt[im][jm] = 0.0;
      }
   }
   for (int im=0; im<4; im++)
      identity_matt[im][im] = 1.0;  

   // There should be only one of these in the new (may 2004) system:
   for (unsigned int ii=0; ii<symm_trans.size(); ii++) {

      mat44 mol_to_origin_matt;
      SelAtom.mol->GetTMatrix(mol_to_origin_matt, 0,
			      -symm_trans[ii].second.us,
			      -symm_trans[ii].second.vs,
			      -symm_trans[ii].second.ws);
      
      int ierr = SelAtom.mol->GetTMatrix(my_matt,
					 symm_trans[ii].first.isym(),
					 symm_trans[ii].first.x(),
					 symm_trans[ii].first.y(),
					 symm_trans[ii].first.z());

      if (ierr != 0) {
	 std::cout << "!!!!!!!! something BAD with CMMDBCryst.GetTMatrix"
	      << std::endl;
      }
      
      int nmodels = SelAtom.mol->GetNumberOfModels();
      for (int imodel=1; imodel<=nmodels; imodel++) {
	 // Here we want to get a selection of all atoms that have
	 // been translated to this symm_trans, then calculate bonds
	 // on them.
	 CAtom **transsel = new PCAtom[SelAtom.n_selected_atoms];
	 for (int i=0; i<SelAtom.n_selected_atoms; i++) {
	    transsel[i] = new CAtom;
	    transsel[i]->Copy(SelAtom.atom_selection[i]);
	    transsel[i]->residue = SelAtom.atom_selection[i]->residue;
	    transsel[i]->Transform(mol_to_origin_matt);
	    transsel[i]->Transform(my_matt);
	 }

	 // So now we have transsel.  We need to calculate bonds from
	 // these:

	 short int atom_colour_type = coot::COLOUR_BY_ATOM_TYPE;
	 construct_from_atom_selection(SelAtom,
				       transsel, SelAtom.n_selected_atoms,
				       transsel, SelAtom.n_selected_atoms,
				       0.1, 1.8,
				       atom_colour_type,
				       0, 1, SelAtom.UDDAtomIndexHandle);

	 for (int i=0; i<SelAtom.n_selected_atoms; i++)
	    delete transsel[i]; // delete the atoms
	 delete [] transsel;  // delete the atom list
      }
   }
   if (symm_trans.size() > 0)
      gbc = make_graphical_symmetry_bonds();
   return gbc;
}

 
// 
graphical_bonds_container
Bond_lines_container::addSymmetry_with_mmdb(const atom_selection_container_t &SelAtom,
					    coot::Cartesian point,
					    float symm_distance,
					    const std::vector<std::pair<symm_trans_t, Cell_Translation> > &symm_trans,
					    short int symmetry_as_ca_flag) {

   graphical_bonds_container gbc; 

   if (symmetry_as_ca_flag == 1) {
      gbc = addSymmetry_calphas(SelAtom, point, symm_distance, symm_trans);

   } else {
      
      PSContact contact;
      int ncontacts;

      if (SelAtom.n_selected_atoms > 0) { 
	 PCAtom point_atom_p = new CAtom;
	 point_atom_p->SetCoordinates(point.get_x(), point.get_y(),
				      point.get_z(), 1.0, 99.9);

	 for(unsigned int ii=0; ii<symm_trans.size(); ii++) { // i.e. boxes

	    // 	 cout << "---------------------------- "
	    // 	      << symm_trans[ii]  
	    // 	      << " ---------------------------- "
	    // 	      << endl;
	 
	    // now get an atom selection where all the atoms are moved
	    // by symm_trans[ii]
	    //
	    // Create a selection trans_selection that is a copy of
	    // SelAtom.atom_selection that has had symmetry
	    // transformation symm_trans[ii] applied to it.
	    // 
	    PPCAtom trans_selection = trans_sel(SelAtom, symm_trans[ii]);

	    contact = NULL;
	    SelAtom.mol->SeekContacts(point_atom_p, trans_selection,
				      SelAtom.n_selected_atoms,
				      0.0, symm_distance,
				      0,  // seqDist
				      contact, ncontacts);
	 
	    // cout << symm_trans[ii] << " had " << ncontacts
	    // << " symmetry atom within "
	    // << symm_distance << "A radius." << endl;

	    if (ncontacts > 0 ) {

	       // the atom_selection of Contact_Sel is allocated,
	       // where do we give it up?
	       // 
	       atom_selection_container_t Contact_Sel = 
		  ContactSel(trans_selection, contact, ncontacts); 
	       Contact_Sel.mol = SelAtom.mol;

	       construct_from_asc(Contact_Sel, 0.01, 1.95, coot::COLOUR_BY_ATOM_TYPE, 1);
	       gbc = make_graphical_symmetry_bonds();

	       // Now give back the atom_selection, (but not the atoms
	       // because they are created by trans_sel and are given
	       // back later.     giveback_1
	       //
	       delete [] Contact_Sel.atom_selection;
	    
	    } else {

	       // do a constructor for no symmetry bonds
	       // and reset graphical_symmetry_bonds.

	       no_symmetry_bonds();

	    }

	    // Tidy up.
	    // 
	    if (trans_selection) { 
	       for (int j=0; j<SelAtom.n_selected_atoms; j++)
		  // delete each of the atoms
		  delete trans_selection[j];
	       delete [] trans_selection;
	    }

	    delete [] contact;
	 }
	 delete point_atom_p;
      }
   }
   return gbc; 
}

// symmetry with colour-by-symmetry-operator:
std::vector<std::pair<graphical_bonds_container, std::pair<symm_trans_t, Cell_Translation> > >
Bond_lines_container::addSymmetry_vector_symms(const atom_selection_container_t &SelAtom,
	    coot::Cartesian point,
	    float symm_distance,
	    const std::vector<std::pair<symm_trans_t, Cell_Translation> > &symm_trans,
	    short int symmetry_as_ca_flag,
            short int symmetry_whole_chain_flag,
	    short int draw_hydrogens_flag) {

   std::vector<std::pair<graphical_bonds_container, std::pair<symm_trans_t, Cell_Translation> > > r;
   do_bonds_to_hydrogens = draw_hydrogens_flag;

   for (unsigned int i=0; i<symm_trans.size(); i++) {
      std::pair<symm_trans_t, Cell_Translation>  this_symm_trans = symm_trans[i];
      std::vector<std::pair<symm_trans_t, Cell_Translation> > this_symm_trans_vec;
      this_symm_trans_vec.push_back(this_symm_trans);
      r.push_back(std::pair<graphical_bonds_container,
		  std::pair<symm_trans_t, Cell_Translation > > (addSymmetry(SelAtom, 
									    point, symm_distance,
									    this_symm_trans_vec,
									    symmetry_as_ca_flag,
									    symmetry_whole_chain_flag), 
								symm_trans[i]));
   }
   return r;
}



graphical_bonds_container
Bond_lines_container::addSymmetry_calphas(const atom_selection_container_t &SelAtom,
					  const coot::Cartesian &point,
					  float symm_distance,
					  const std::vector<std::pair<symm_trans_t, Cell_Translation> > &symm_trans) {

   graphical_bonds_container gbc;
   // std::cout << "DEBUG:: There are " << symm_trans.size() << "
   // symm_trans units" << std::endl;
   mat44 my_matt;
   CAtom t_atom1;
   CAtom t_atom2;

   for (unsigned int ii=0; ii<symm_trans.size(); ii++) {
      CMMDBCryst *cryst_p =  (CMMDBCryst *) &SelAtom.mol->get_cell();

      int err = cryst_p->GetTMatrix(my_matt, symm_trans[ii].first.isym(), symm_trans[ii].first.x(),
				    symm_trans[ii].first.y(), symm_trans[ii].first.z());

      mat44 mol_to_origin_matt;
      cryst_p->GetTMatrix(mol_to_origin_matt, 0,
			  -symm_trans[ii].second.us,
			  -symm_trans[ii].second.vs,
			  -symm_trans[ii].second.ws);
      
      if (err != 0) {
	 cout << "!!!!!!!!!!!!!! something BAD with CMMDBCryst.GetTMatrix"
	      << endl;
      }
      
      int nmodels = SelAtom.mol->GetNumberOfModels();
      for (int imodel=1; imodel<=nmodels; imodel++) {
	 int nchains = SelAtom.mol->GetNumberOfChains(imodel);
	 for (int ichain=0; ichain<nchains; ichain++) {
	    PCChain chain = SelAtom.mol->GetChain(imodel, ichain);
	    int nres = chain->GetNumberOfResidues();
	    for (int ires=0; ires<(nres-1); ires++) {
	       // std::cout << "residue " << ires << " to  " << ires+1 << std::endl;
	       // one residue to the next...
	       PCResidue res1 = chain->GetResidue(ires);
	       PCResidue res2 = chain->GetResidue(ires+1);

	       // Search through the atoms of each residue looking for CAs:
	       std::vector<CAtom *> ca_this;
	       std::vector<CAtom *> ca_next;

	       PPCAtom residue_atoms;
	       int n_residue_atoms;
	       res1->GetAtomTable(residue_atoms, n_residue_atoms);
	       for (int iat=0; iat<n_residue_atoms; iat++) {
		  if (residue_atoms[iat]) { 
		     if (std::string(residue_atoms[iat]->name) == " CA " ||
			 std::string(residue_atoms[iat]->name) == " P  ") {
			ca_this.push_back(residue_atoms[iat]);
		     }
		  }
	       }
	       res2->GetAtomTable(residue_atoms, n_residue_atoms);
	       for (int iat=0; iat<n_residue_atoms; iat++) {
		  if (residue_atoms[iat]) { 
		     if (std::string(residue_atoms[iat]->name) == " CA " ||
			 std::string(residue_atoms[iat]->name) == " P  ") {
			ca_next.push_back(residue_atoms[iat]);
		     }
		  }
	       }
	       // so now we have the 2 vectors filled with CA atoms of
	       // this and the next residue.  What are the connections?
	       // 
	       if (ca_this.size() > 0) {
		  if (ca_next.size() > 0) {
		     for (int iat=0; iat<ca_this.size(); iat++) { 
			std::string altconf1 = ca_this[iat]->altLoc;
			for (int jat=0; jat<ca_next.size(); jat++) { 
			   std::string altconf2 = ca_next[jat]->altLoc;

			   if ((altconf1 == altconf2) ||
			       (altconf1 == "") ||
			       (altconf2 == "")) {
			      coot::Cartesian ca_1(ca_this[iat]->x,
						   ca_this[iat]->y,
						   ca_this[iat]->z);
			      coot::Cartesian ca_2(ca_next[jat]->x,
						   ca_next[jat]->y,
						   ca_next[jat]->z);

			      double len = (ca_1 - ca_2).amplitude();
			      // CA-CA or P-P
			      if (((len < 4.7) && (len > 2.4)) ||
				  ((len<8) && (len>5))) {
				 
				 t_atom1.Copy(ca_this[iat]);
				 t_atom2.Copy(ca_next[jat]);
				 t_atom1.Transform(mol_to_origin_matt);
				 t_atom2.Transform(mol_to_origin_matt);
				 t_atom1.Transform(my_matt);
				 t_atom2.Transform(my_matt);

				 coot::Cartesian atom1(t_atom1.x, t_atom1.y, t_atom1.z);
				 coot::Cartesian atom2(t_atom2.x, t_atom2.y, t_atom2.z);
				 
				 addBond(0, atom1, atom2); // these are C atoms.
			      }
			   }
			}
		     }
		  } 
	       }
	    }
	 }
      }
   }

   if (symm_trans.size() > 0)
      gbc = make_graphical_symmetry_bonds();

   return gbc;
}

std::vector<std::pair<graphical_bonds_container, symm_trans_t> >
Bond_lines_container::add_NCS(const atom_selection_container_t &SelAtom,
	    coot::Cartesian point,
	    float symm_distance,
            std::vector<std::pair<coot::coot_mat44, symm_trans_t> > &strict_ncs_mats,
            short int symmetry_as_ca_flag,
            short int symmetry_whole_chain_flag) {

   std::vector<std::pair<graphical_bonds_container, symm_trans_t> > r;
   for (unsigned int i=0; i<strict_ncs_mats.size(); i++) {
      r.push_back(std::pair<graphical_bonds_container, symm_trans_t> (add_NCS_molecule(SelAtom, point, symm_distance, strict_ncs_mats[i], symmetry_as_ca_flag, symmetry_whole_chain_flag), strict_ncs_mats[i].second));
   }
   return r;
}

graphical_bonds_container
Bond_lines_container::add_NCS_molecule(const atom_selection_container_t &SelAtom,
				       const coot::Cartesian &point,
				       float symm_distance,
				       const std::pair<coot::coot_mat44, symm_trans_t> &strict_ncs_mat,
				       short int symmetry_as_ca_flag,
				       short int symmetry_whole_chain_flag) {

   graphical_bonds_container gbc;
   mat44 m;

   if (symmetry_as_ca_flag) {
      gbc = add_NCS_molecule_calphas(SelAtom, point, symm_distance, strict_ncs_mat);
   } else { 
      if (symmetry_whole_chain_flag) {
	 gbc = add_NCS_molecule_whole_chain(SelAtom, point, symm_distance, strict_ncs_mat);
      } else { 
   
	 m[0][0] = strict_ncs_mat.first.m[0].v4[0];
	 m[0][1] = strict_ncs_mat.first.m[0].v4[1];
	 m[0][2] = strict_ncs_mat.first.m[0].v4[2];
	 m[1][0] = strict_ncs_mat.first.m[1].v4[0];
	 m[1][1] = strict_ncs_mat.first.m[1].v4[1];
	 m[1][2] = strict_ncs_mat.first.m[1].v4[2];
	 m[2][0] = strict_ncs_mat.first.m[2].v4[0];
	 m[2][1] = strict_ncs_mat.first.m[2].v4[1];
	 m[2][2] = strict_ncs_mat.first.m[2].v4[2];
	 m[0][3] = strict_ncs_mat.first.m[0].v4[3];  // t
	 m[1][3] = strict_ncs_mat.first.m[1].v4[3];  // t
	 m[2][3] = strict_ncs_mat.first.m[2].v4[3];  // t
	 m[3][0] = strict_ncs_mat.first.m[3].v4[0];
	 m[3][1] = strict_ncs_mat.first.m[3].v4[1];
	 m[3][2] = strict_ncs_mat.first.m[3].v4[2];
	 m[3][3] = strict_ncs_mat.first.m[3].v4[3];

	 // gbc = make_graphical_symmetry_bonds();
      }
   }
   return gbc;
}

graphical_bonds_container
Bond_lines_container::add_NCS_molecule_whole_chain(const atom_selection_container_t &SelAtom,
						   const coot::Cartesian &point,
						   float symm_distance,
						   const std::pair<coot::coot_mat44, symm_trans_t> &strict_ncs_mat) {

   graphical_bonds_container gbc;
   mat44 my_matt;
   
   my_matt[0][0] = strict_ncs_mat.first.m[0].v4[0];
   my_matt[0][1] = strict_ncs_mat.first.m[0].v4[1];
   my_matt[0][2] = strict_ncs_mat.first.m[0].v4[2];
   my_matt[1][0] = strict_ncs_mat.first.m[1].v4[0];
   my_matt[1][1] = strict_ncs_mat.first.m[1].v4[1];
   my_matt[1][2] = strict_ncs_mat.first.m[1].v4[2];
   my_matt[2][0] = strict_ncs_mat.first.m[2].v4[0];
   my_matt[2][1] = strict_ncs_mat.first.m[2].v4[1];
   my_matt[2][2] = strict_ncs_mat.first.m[2].v4[2];
   my_matt[0][3] = strict_ncs_mat.first.m[0].v4[3];  // t
   my_matt[1][3] = strict_ncs_mat.first.m[1].v4[3];  // t
   my_matt[2][3] = strict_ncs_mat.first.m[2].v4[3];  // t
   my_matt[3][0] = strict_ncs_mat.first.m[3].v4[0];
   my_matt[3][1] = strict_ncs_mat.first.m[3].v4[1];
   my_matt[3][2] = strict_ncs_mat.first.m[3].v4[2];
   my_matt[3][3] = strict_ncs_mat.first.m[3].v4[3];

   int nmodels = SelAtom.mol->GetNumberOfModels();
   for (int imodel=1; imodel<=nmodels; imodel++) {
      // Here we want to get a selection of all atoms that have
      // been translated to this symm_trans, then calculate bonds
      // on them.
      CAtom **transsel = new PCAtom[SelAtom.n_selected_atoms];
      for (int i=0; i<SelAtom.n_selected_atoms; i++) {
	 transsel[i] = new CAtom;
	 transsel[i]->Copy(SelAtom.atom_selection[i]);
	 transsel[i]->residue = SelAtom.atom_selection[i]->residue;
	 transsel[i]->Transform(my_matt);
      }
      // So now we have transsel.  We need to calculate bonds from
      // these:
      
      short int atom_colour_type = coot::COLOUR_BY_ATOM_TYPE;
      construct_from_atom_selection(SelAtom,
				    transsel, SelAtom.n_selected_atoms,
				    transsel, SelAtom.n_selected_atoms,
				    0.1, 1.8,
				    atom_colour_type,
				    0, 1, SelAtom.UDDAtomIndexHandle);
//       for (int i=0; i<SelAtom.n_selected_atoms; i++) {
// 	 std::cout << "atom sel[" << i << "] " << SelAtom.atom_selection[i] << "\n";
// 	 std::cout << "transsel[" << i << "] " << transsel[i] << "\n";
//       }
      gbc = make_graphical_symmetry_bonds(); 
      for (int i=0; i<SelAtom.n_selected_atoms; i++)
	 delete transsel[i];
      delete [] transsel;
      
   }
   return gbc;
}

graphical_bonds_container
Bond_lines_container::add_NCS_molecule_calphas(const atom_selection_container_t &SelAtom,
					       const coot::Cartesian &point,
					       float symm_distance,
					       const std::pair<coot::coot_mat44, symm_trans_t> &strict_ncs_mat) {

   graphical_bonds_container gbc;

   return gbc;

}

// We have a set of contacts in contact for the atom selection trans_sel.
//
// We want to return an atom selection that is just the contact atoms.
// We also need the number of atom (which is the the same as the
// number of atom contacts), so we fill in the 2 parts of an
// atom_selection_container_t.
//
atom_selection_container_t
Bond_lines_container::ContactSel(PPCAtom trans_sel,
				 PSContact contact, int ncontacts) const {

   // We need to sort the contacts.
   // 
   // If we were an Object Oriented Programmer, we would expect to
   // sort contacts with contact->Sort().  Does that work with mmdb?
   // 
   // Hah, you're kidding me, right?  Grumble...
   //
   SortContacts(contact, ncontacts, CNSORT_2INC);
   //
   int id2, last_id2 = -1;
   int n_contact_atoms = 0;
   //
   // keep in mind that there is a contact id2->id1 as well as id1->id2.

   atom_selection_container_t TransSel;

   // You might ask: where does this get deleted?  It does get deleted: giveback_1
   // 
   TransSel.atom_selection = new PCAtom[ncontacts]; // 

   for (int icon=0; icon<ncontacts; icon++) {
      // use sorted contacts
      id2 = contact[icon].id2;
      TransSel.atom_selection[n_contact_atoms] = trans_sel[id2];
      n_contact_atoms++;
   }
   
   TransSel.n_selected_atoms = n_contact_atoms;
   if (n_contact_atoms > ncontacts){ 
      cout << "disaster n_contact_atoms (" << n_contact_atoms
	   << ") > ncontacts (" << ncontacts << ")" << endl;
   }
   return TransSel;

}


// See above comments re const CMMDBCryst &.
// 
PPCAtom
Bond_lines_container::trans_sel(atom_selection_container_t AtomSel,
				const std::pair<symm_trans_t, Cell_Translation> &symm_trans) const {
   mat44 my_matt;
   //
   // modify my_matt;
   // 
   // GetTMatrix is not const:  Arrrrrrrggggh!  Ghod this is *so* frustrating.
   //
   // I want to say "simply":
   // AtomSel.mol->get_cell().GetTMatrix(my_matt, symm_trans.isym(), symm_trans.x(),
   //                                    symm_trans.y(), symm_trans.z());
   // so this is the ugly hack work around...
   //
   // Note that we have to move to a pointer because CMMDBCryst
   // destructor will trash an CMMDBCryst returned from get_cell()
   // (which is actually a const ref to the real one).
   //
   // And we do (CMMDBCryst *) casting because &AtomSel.mol->get_cell();
   // returns (as I just said) a const ref and the compiler cannot
   // convert const `CMMDBCryst *' to `CMMDBCryst *' (sensibly enough).
   //
   // And we can't type cryst_p as const because (see above) GetTMatrix is
   // not a const function and can't use a const object.
   //
   // Fun eh?
   // 
   CMMDBCryst *cryst_p =  (CMMDBCryst *) &AtomSel.mol->get_cell();

   int err = cryst_p->GetTMatrix(my_matt, symm_trans.first.isym(), symm_trans.first.x(),
				 symm_trans.first.y(), symm_trans.first.z());

   if (err != 0) {
      cout << "!!!!!!!!!!!!!! something BAD with CMMDBCryst.GetTMatrix"
	   << endl;
   }
   
   // cout << "using symm_trans: " << symm_trans << endl;

   //cout << "applying my_matt in trans_sel: " << endl;
//    cout << my_matt[0][0] << " "  << my_matt[0][1] << " "
// 	<< my_matt[0][2] << " "  << my_matt[0][3] << " "  << endl
// 	<< my_matt[1][0] << " "  << my_matt[1][1] << " "
// 	<< my_matt[1][2] << " "  << my_matt[1][3] << " "  << endl
// 	<< my_matt[2][0] << " "  << my_matt[2][1] << " "
// 	<< my_matt[2][2] << " "  << my_matt[2][3] << " "  << endl
// 	<< my_matt[3][0] << " "  << my_matt[3][1] << " "
// 	<< my_matt[3][2] << " "  << my_matt[3][3] << " "  << endl;


   mat44 mol_to_origin_matt;
   cryst_p->GetTMatrix(mol_to_origin_matt, 0,
		       -symm_trans.second.us,
		       -symm_trans.second.vs,
		       -symm_trans.second.ws);
		       
   // Whoah!  Big alloc!  Given back in addSymmetry
   // 
   //cout << "allocating translation space for " << AtomSel.n_selected_atoms
   //	<< " atoms" << endl;
   // 
   // 
   PPCAtom trans_selection = new PCAtom[AtomSel.n_selected_atoms];
   for (int ii=0; ii<AtomSel.n_selected_atoms; ii++) {

      trans_selection[ii] = new CAtom;

//       trans_selection[ii]->SetCoordinates(AtomSel.atom_selection[ii]->x,
// 					  AtomSel.atom_selection[ii]->y,
// 					  AtomSel.atom_selection[ii]->z,
// 					  1.0, 99.9);
//       trans_selection[ii]->SetAtomName(   AtomSel.atom_selection[ii]->name);
//       trans_selection[ii]->SetElementName(AtomSel.atom_selection[ii]->element);
//       trans_selection[ii]->SetResidue(    AtomSel.atom_selection[ii]->GetResidue());
//       trans_selection[ii]->tempFactor =   AtomSel.atom_selection[ii]->tempFactor;


      
      // can't set seqNum:      
      //trans_selection[ii]->seqNum     =   AtomSel.atom_selection[ii]->GetSeqNum();
      // I don't know why, but this doesn't work.
      // trans_selection[ii]->SetChain(      AtomSel.atom_selection[ii]->GetChain());

      trans_selection[ii]->Copy(AtomSel.atom_selection[ii]);
      trans_selection[ii]->residue = AtomSel.atom_selection[ii]->residue;
      trans_selection[ii]->Transform(mol_to_origin_matt);
      trans_selection[ii]->Transform(my_matt);
   }
   return trans_selection;
}

void
Bond_lines_container::do_disulphide_bonds(atom_selection_container_t SelAtom,
					  int imodel) {

   // Lets make a new selection using mol.  We step back and sidestep
   // to a different atom selection.
   // 
   PPCAtom Sulfur_selection;
   int n_sulfurs;
   mat44 my_matt;
   PSContact contact = NULL;
   int ncontacts = 0; // initially no S contacts.
   long i_contact_group = 1;
   int col; 

   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;
      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

   
   int selHnd2 = SelAtom.mol->NewSelection();

   // model 1
   // Note that now we force the resname to be CYS and atom name to be SG
   SelAtom.mol->SelectAtoms(selHnd2, imodel,"*", ANY_RES, "*", ANY_RES, "*",
			    "CYS"," SG ","S","*" );

   SelAtom.mol->GetSelIndex(selHnd2, Sulfur_selection, n_sulfurs);

   if (n_sulfurs > 0) { 
      SelAtom.mol->SeekContacts(Sulfur_selection, n_sulfurs,
				Sulfur_selection, n_sulfurs,
				0.01, 3.0, // min, max dist.
				0,         // in same res also.
				contact, ncontacts, 
				0, &my_matt, i_contact_group);
   }

//    if (verbose_reporting)
//       cout << "Found " << ncontacts/2 << " disulfides" << endl;

   if (ncontacts > 0) {
      for (int i=0; i< ncontacts; i++) {
	 
	 if ( contact[i].id2 > contact[i].id1 ) {

	    std::string aloc_1(Sulfur_selection[ contact[i].id1 ]->altLoc);
	    std::string aloc_2(Sulfur_selection[ contact[i].id2 ]->altLoc);

	    if ( (aloc_1=="") || (aloc_2=="") || (aloc_1==aloc_2) ) {
	       coot::Cartesian atom_1(Sulfur_selection[ contact[i].id1 ]->x,
				Sulfur_selection[ contact[i].id1 ]->y,
				Sulfur_selection[ contact[i].id1 ]->z);

	       coot::Cartesian atom_2(Sulfur_selection[ contact[i].id2 ]->x,
				Sulfur_selection[ contact[i].id2 ]->y,
				Sulfur_selection[ contact[i].id2 ]->z);

	       col = atom_colour(Sulfur_selection[ contact[i].id1 ],
				 coot::DISULFIDE_COLOUR);
	       
	       if (! ((Sulfur_selection[ contact[i].id1 ]->GetSeqNum() ==
		       Sulfur_selection[ contact[i].id2 ]->GetSeqNum()) &&
		      (Sulfur_selection[ contact[i].id1 ]->GetChainID() ==
		       Sulfur_selection[ contact[i].id2 ]->GetChainID()))) { 
		  addBond(col, atom_1, atom_2);
	       }
	    }
	 }
      }
      delete [] contact;
   }
   SelAtom.mol->DeleteSelection(selHnd2);
}

// reset bonds (delete if necessary) any graphics bonds.
// 
void Bond_lines_container::no_symmetry_bonds() {

}


Bond_lines_container::Bond_lines_container(symm_keys key) {

   do_bonds_to_hydrogens = 1;  // added 20070629
   b_factor_scale = 1.0;
   if (key == NO_SYMMETRY_BONDS) {

      no_symmetry_bonds();
      
   } else { 
      std::cout << "Bond_lines_container::Bond_lines_container(symm_keys key)"
		<< " no such key as " << key << std::endl;
   }
}
      
   

void Bond_lines_container::check_static() const {

	graphical_bonds_container pot; 
	
        cout << "check: num_colours:"     << pot.num_colours << endl;
	cout << "check: bonds:"           << pot.bonds_ << endl;
        cout << "check: bonds::numlines " << pot.bonds_[12].num_lines << endl;

}

graphical_bonds_container 
Bond_lines_container::make_graphical_bonds() const {
   return make_graphical_bonds(1); // allow simple-minded thinning
}


graphical_bonds_container 
Bond_lines_container::make_graphical_bonds_no_thinning() const {
   return make_graphical_bonds(0); // no thinning
} 

graphical_bonds_container 
Bond_lines_container::make_graphical_bonds(bool thinning_flag) const {

   graphical_bonds_container box;

   box.num_colours = bonds.size();
   box.bonds_ = new Lines_list[bonds.size()];

   // i is the colour index
   int ibs = bonds.size();
   for (int i=0; i<ibs; i++) {

      box.bonds_[i].num_lines = bonds[i].size();
      box.bonds_[i].pair_list = new coot::CartesianPair[bonds[i].size()];
      for (int j=0; j<bonds[i].size(); j++) { 
	 box.bonds_[i].pair_list[j] = bonds[i][j];
	 if (thinning_flag)
	    if (j == HYDROGEN_GREY_BOND)
	       box.bonds_[j].thin_lines_flag = 1;
      }
   }
   box.add_zero_occ_spots(zero_occ_spot);
   box.add_atom_centres(atom_centres, atom_centres_colour);
   return box;
}

graphical_bonds_container
Bond_lines_container::make_graphical_symmetry_bonds() const {
 
   graphical_bonds_container box;
   box.num_colours = bonds.size();
   box.bonds_ = NULL;
   // box.bonds_ = new Lines_list[bonds.size()]; // surely not!
   box.symmetry_has_been_created = 1;

   // This block can never happen!  box is created here and the
   // constructor sets symmetry_bonds_ to NULL.
   // 
   // delete the old bonds
   // 
   if (box.symmetry_bonds_ != NULL) {
      //
      //cout << "deleting symmetry bonds (bonds.size() is " << bonds.size()
      //	   << ")" << endl;
      for (int i = 0; i < box.num_colours; i++) {
	 if (box.symmetry_bonds_[i].num_lines > 0) {
	    //cout << "box.symmetry_bonds_[" << i << "].num_lines is "
	    //<< box.symmetry_bonds_[i].num_lines << endl;
	    delete [] box.symmetry_bonds_[i].pair_list;
	 }
      }
   }
      
   box.symmetry_bonds_ = new Lines_list[bonds.size()];
//    std::cout << "allocating symmetry bonds " << box.symmetry_bonds_
// 	     << " for " << bonds.size() << " symmetry bonds" << std::endl;

   int num_lines;
   for (int i = 0; i < box.num_colours; i++) {

      // i is the colour
      //
      num_lines =  bonds[i].size();
      box.symmetry_bonds_[i].num_lines = num_lines;

      if (num_lines > 0 ){ 
      
	 box.symmetry_bonds_[i].pair_list = new coot::CartesianPair[bonds[i].size()];

	 int bis = bonds[i].size();
	 for (int j=0; j<bis; j++) {
	    box.symmetry_bonds_[i].pair_list[j] = bonds[i][j];
	 }
      }
   }
   return box; 
}

void
Bond_lines_container::check() const {
   //
   cout << "Bond_lines_container::check() bonds.size() " << bonds.size() << endl;
   if (bonds.size() > 0) {
      cout <<  "Bond_lines_container::check() bonds[0].size(): "
	   << bonds[0].size() << endl; 
   }
   if (bonds.size() > 1) {
      cout <<  "Bond_lines_container::check() bonds[1].size(): "
	   << bonds[1].size() << endl; 
   }
}


// initialise
// graphical_bonds_container::bonds

//
void
Bond_lines_container::write(std::string filename) const {


   cout << "Write bonds to file: " << filename.c_str() << endl;
   
   std::ofstream bondsout(filename.c_str()); 
   if (! bondsout) {
      // error
      cout << "Could not open " << filename.c_str() << " for some reason\n";
   } else { 
   
      for (unsigned int i = 0; i < bonds.size(); i++) {
	 
	 bondsout<< bonds[i].size() << " bonds of colour " << i << endl;

	 int bis=bonds[i].size();
	 for (int j = 0; j<bis ; j++) {

	    // This gets it pass CC, strangely
            bondsout << bonds[i].GetStart(j); 
            bondsout << " to ";
            bondsout << bonds[i].GetFinish(j) << endl;
	 }
      }
   }
}

coot::CartesianPair
Bond_lines::operator[](int i) const {

   return points[i];
}


coot::Cartesian
Bond_lines::GetStart(int i) const {

   return points[i].getStart();
}

coot::Cartesian
Bond_lines::GetFinish(int i) const {

   return points[i].getFinish();
} 

Bond_lines_container::Bond_lines_container(int col) {

   do_bonds_to_hydrogens = 1;  // added 20070629
   b_factor_scale = 1.0;
   std::cout << "Strange Bond_lines_container(int col)" << std::endl;
   Bond_lines a;
   bonds.push_back(a);
}

//
void
Bond_lines_container::addBond(int col,
			      const coot::Cartesian &start,
			      const coot::Cartesian &end) {

   coot::CartesianPair pair(start,end);
   bonds[col].add_bond(pair);
}

//
void
Bond_lines_container::add_dashed_bond(int col,
				      const coot::Cartesian &start,
				      const coot::Cartesian &end,
				      bool is_half_bond_flag) {

   int n_dash = 16;
   if (is_half_bond_flag)
      n_dash = 8;

   for (int idash=0; idash<n_dash; idash+=2) {
      float frac_1 = float(idash  )/float(n_dash);
      float frac_2 = float(idash+1)/float(n_dash);
      coot::Cartesian this_start(start + (end-start).by_scalar(frac_1));
      coot::Cartesian this_end(  start + (end-start).by_scalar(frac_2));
      coot::CartesianPair pair(this_start,this_end);
      bonds[col].add_bond(pair);
   } 
}

//
void
Bond_lines::add_bond(coot::CartesianPair pair) {
   points.push_back(pair);
} 

//
Bond_lines::Bond_lines() {

   // This gets called when we resize a Bond_lines_container's bonds array.
   // 
   // std::cout << "nothing much" << std::endl;
}

//
Bond_lines::Bond_lines(int col) {
   colour = col;
}

int 
Bond_lines::size() const {
   return points.size();
}



// ------------------------------------------------------------------
//            Alpha Carbon Trace, C-alpha calpha
// ------------------------------------------------------------------


// The distances are the minimum and maximum distances to look check
// between for ca-ca (pseudo) bond.  Typically 3.6-3.8 (tight) or
// 0.0-5.0 (loose).
// 
void
Bond_lines_container::do_Ca_bonds(atom_selection_container_t SelAtom,
				  float min_dist, float max_dist) {

   udd_has_ca_handle = SelAtom.mol->RegisterUDInteger (UDR_RESIDUE, "has CA");
   if (!udd_has_ca_handle) {
      std::cout << "ERROR getting udd_has_ca_handle\n";
   }
   coot::my_atom_colour_map_t atom_colour_map_in;
   coot::my_atom_colour_map_t atom_colour_map = 
      do_Ca_or_P_bonds_internal(SelAtom, " CA ", atom_colour_map_in,
				min_dist, max_dist, coot::COLOUR_BY_CHAIN);
   do_Ca_or_P_bonds_internal(SelAtom, " P  ",  atom_colour_map,
			     0.1,      7.5, coot::COLOUR_BY_CHAIN);
}


// where bond_colour_type is one of 
// enum bond_colour_t { COLOUR_BY_CHAIN=0, COLOUR_BY_SEC_STRUCT=1};
// 
coot::my_atom_colour_map_t 
Bond_lines_container::do_Ca_or_P_bonds_internal(atom_selection_container_t SelAtom,
						const char *backbone_atom_id,
						coot::my_atom_colour_map_t atom_colour_map,
						float min_dist, float max_dist, int bond_colour_type) {

//    std::cout << "DEBUG:: do_Ca_or_P_bonds_internal with bond_colour_type "
// 	     << bond_colour_type << std::endl;

   PPCAtom Ca_selection = NULL; 
   int n_ca;
   mat44 my_matt;
   int ncontacts;
   long i_contact_group = 1;
   PSContact contact = NULL;
   int col;

   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;
  
   int selHnd_ca = SelAtom.mol->NewSelection();
   const char *ele = "C";
   if (! strcmp(backbone_atom_id, " P  "))
      ele = "P";
   
   SelAtom.mol->SelectAtoms(selHnd_ca, 0,"*",ANY_RES,"*",ANY_RES,"*",
			    "*", backbone_atom_id, ele,"*" );

   SelAtom.mol->GetSelIndex(selHnd_ca, Ca_selection, n_ca);
   int atom_colours_udd = -1; // unset/bad
   if (bond_colour_type == coot::COLOUR_BY_RAINBOW)
      atom_colours_udd = set_rainbow_colours(selHnd_ca, SelAtom.mol);

   SelAtom.mol->SeekContacts(Ca_selection, n_ca,
			     Ca_selection, n_ca,
			     min_dist, max_dist, // min, max dist.
			     0,         // in same res also.
			     contact, ncontacts, 
			     0, &my_matt, i_contact_group);

   // std::cout << "Found " << ncontacts/2 << " Ca-Ca links" << std::endl;

   // Note that contact can be null when ncontacts is 2.

   // 
   short int *found_contact = new short int[SelAtom.n_selected_atoms];
   for (int i=0; i<SelAtom.n_selected_atoms; i++)
      found_contact[i] = 0;

   if (ncontacts > 0) {

      if (contact != NULL) {
	 int istat;
	 // zero the array.
	 for (int i=0; i<SelAtom.n_selected_atoms; i++) found_contact[i] = 0;
      
	 for (int i=0; i< ncontacts; i++) {
	    std::string chainid1;
	    std::string chainid2;
	    if ( contact[i].id2 >  contact[i].id1 ) {

	       int res_1 = Ca_selection[ contact[i].id1 ]->GetSeqNum();
	       int res_2 = Ca_selection[ contact[i].id2 ]->GetSeqNum();

	       if (udd_has_ca_handle >= 1) {
		  istat = Ca_selection [ contact[i].id1 ]->residue->PutUDData(udd_has_ca_handle, 1);
 		  if (istat == UDDATA_WrongUDRType)
 		     std::cout << "ERROR::  UDDATA_WrongUDRType in do_Ca_bonds 1" << std::endl;
		  istat = Ca_selection [ contact[i].id2 ]->residue->PutUDData(udd_has_ca_handle, 1);
 		  if (istat == UDDATA_WrongUDRType)
 		     std::cout << "ERROR::  UDDATA_WrongUDRType in do_Ca_bonds 2" << std::endl;
	       }
	    
	       // this +/- 1 residue test
	       if (abs(res_1 - res_2) < 2) {

		  std::string altloc1 = Ca_selection[ contact[i].id1 ]->altLoc;
		  std::string altloc2 = Ca_selection[ contact[i].id2 ]->altLoc;
		  
		  if (Ca_selection[ contact[i].id1 ]->GetModel() ==
		      Ca_selection[ contact[i].id2 ]->GetModel()) {


		     if (altloc1 == "" || altloc2 == "" || altloc1 == altloc2) { 
 
			found_contact[contact[i].id1] = 1;  // true
			found_contact[contact[i].id2] = 1;
		     
			coot::Cartesian ca_1(Ca_selection[ contact[i].id1 ]->x,
					     Ca_selection[ contact[i].id1 ]->y,
					     Ca_selection[ contact[i].id1 ]->z);
		     
			coot::Cartesian ca_2(Ca_selection[ contact[i].id2 ]->x,
					     Ca_selection[ contact[i].id2 ]->y,
					     Ca_selection[ contact[i].id2 ]->z);
		     
			chainid1 = (Ca_selection[ contact[i].id1 ]->GetChainID());
			chainid2 = (Ca_selection[ contact[i].id2 ]->GetChainID());

			if (chainid1 == chainid2) {

// 			   std::cout << "bond_colour_type: " << bond_colour_type << " vs "
// 				     << Bond_lines_container::COLOUR_BY_B_FACTOR << " and "
// 				     << coot::COLOUR_BY_SEC_STRUCT
// 				     << std::endl;
			      
			   if (bond_colour_type == Bond_lines_container::COLOUR_BY_B_FACTOR) {

			      coot::Cartesian bond_mid_point = ca_1.mid_point(ca_2);
			      col = atom_colour(Ca_selection[ contact[i].id1 ],
						coot::COLOUR_BY_B_FACTOR);
			      addBond(col, ca_1, bond_mid_point);
			      col = atom_colour(Ca_selection[ contact[i].id2 ],
						coot::COLOUR_BY_B_FACTOR);
			      addBond(col, bond_mid_point, ca_2);
			      
			   } else {
			      
			      if (bond_colour_type == coot::COLOUR_BY_SEC_STRUCT)
				 col = atom_colour(Ca_selection[ contact[i].id1 ], bond_colour_type);
			      else
				 if (bond_colour_type == coot::COLOUR_BY_RAINBOW) {
				    if (atom_colours_udd > 0) {
				       realtype f;
				       if (Ca_selection[contact[i].id1]->GetUDData(atom_colours_udd, f) == UDDATA_Ok) {
					  col = atom_colour_map.index_for_rainbow(f);
				       } else {
					  col = 0;
				       }
				    } else {
				       col = 0;
				    }
				 } else {
				    col = atom_colour_map.index_for_chain(chainid1);
				 }
			      bonds_size_colour_check(col);
			      addBond(col, ca_1, ca_2);
			   }
			}
		     }
		  }
	       }
	    }
	 }
	 delete [] contact;
      }
   }


   // for atoms with no neighbour (contacts):
   coot::Cartesian small_vec_x(0.5, 0.0, 0.0);
   coot::Cartesian small_vec_y(0.0, 0.5, 0.0);
   coot::Cartesian small_vec_z(0.0, 0.0, 0.5);

   for (int i=0; i<n_ca; i++) {

      if ( found_contact[i] == 0 ) {
	 // no contact found
	 col = atom_colour(Ca_selection[i], coot::COLOUR_BY_ATOM_TYPE);
	 coot::Cartesian atom(Ca_selection[i]->x,
			Ca_selection[i]->y,
			Ca_selection[i]->z);

	 addBond(col, atom+small_vec_x, atom-small_vec_x);
	 addBond(col, atom+small_vec_y, atom-small_vec_y);
	 addBond(col, atom+small_vec_z, atom-small_vec_z);
      }
   }

   delete [] found_contact;
   return atom_colour_map;
}


int
Bond_lines_container::set_rainbow_colours(int selHnd_ca, CMMDBManager *mol) {

   PPCAtom Ca_selection = NULL; 
   int n_ca;
   mol->GetSelIndex(selHnd_ca, Ca_selection, n_ca);
   int udd_handle = mol->RegisterUDReal(UDR_ATOM, "rainbow circle point");
   if (udd_handle > 0) { 


      int n_models = mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) { 
	 CModel *model_p = mol->GetModel(imod);
	 int nchains = model_p->GetNumberOfChains();
	 for (int ich=0; ich<nchains; ich++) {

	    CChain *chain_p = model_p->GetChain(ich);
	    const char *chain_id = chain_p->GetChainID();
	    int selhnd_chain = mol->NewSelection();
	    mol->SelectAtoms(selhnd_chain, 0, chain_id , ANY_RES, "*", ANY_RES, "*",
			     "*", " CA ", " C", "*");

	    int n_selected_atoms;
	    PPCAtom atom_selection = NULL;
	    mol->GetSelIndex(selhnd_chain, atom_selection, n_selected_atoms);
// 	    std::cout << "DEBUG:: rainbow " << n_selected_atoms
// 		      << " selected in chain " << chain_id << std::endl;
	    std::pair<short int, int> min_pair = coot::util::min_resno_in_chain(chain_p);
	    std::pair<short int, int> max_pair = coot::util::max_resno_in_chain(chain_p);
	    if (min_pair.first && max_pair.first) {
	       float range = max_pair.second - min_pair.second;
	       for (int iat=0; iat<n_selected_atoms; iat++) {
		  float chain_pos =
		     float(atom_selection[iat]->GetSeqNum() - min_pair.second - 0.01)/range;
		  atom_selection[iat]->PutUDData(udd_handle, chain_pos);
// 		  std::cout << "DEBUG:: atom " << atom_selection[iat]
// 			    << " has rainbow colour " << chain_pos << "\n";
	       }
	    }
	    mol->DeleteSelection(selhnd_chain);
	 }
      }
   }
   return udd_handle;
}


int
Bond_lines_container::atom_colour(CAtom *at, int bond_colour_type,
				  coot::my_atom_colour_map_t *atom_colour_map_p) {

   int col = 0;
   // coot::my_atom_colour_map_t atom_colour_map;
   // std::cout << "In atom-colour with b_factor_scale: " << b_factor_scale << std::endl;

   if (bond_colour_type == coot::COLOUR_BY_CHAIN) {
      if (atom_colour_map_p) { 
	 col = atom_colour_map_p->index_for_chain(std::string(at->GetChainID()));
	 std::cout << " atom_colour_map->index_for_chain(\"" << at->GetChainID()
		   << "\") returns " << col << std::endl;
      }
   } else { 
      if (bond_colour_type == coot::COLOUR_BY_SEC_STRUCT) { 
	 int sse = at->residue->SSE;
	 switch (sse)  {
	 case SSE_None: 
	    col = 0;
	    break;
	 case SSE_Strand:
	    col = 1;
	   break;
	 case SSE_Bulge:  
	    col = 1;
	    break;
	 case SSE_3Turn:  
	    col = 2;
	    break;
	 case SSE_4Turn:  
	    col = 2;
	    break;
	 case SSE_5Turn:  
	    col = 2;
	    break;
	 case SSE_Helix:  
	    col = 2;
	    break;
	 default:
	    col = 3;
	 }
      } else { 
	 if (bond_colour_type == coot::COLOUR_BY_ATOM_TYPE) {
	    std::string element = at->element;

	    if (element == " C") {
	       return YELLOW_BOND;
	    } else {
	       if (element == " N") {
		  return BLUE_BOND;
	       } else {
		  if (element == " O") {
		     return RED_BOND;
		  } else {
		     if (element == " S") {
			return GREEN_BOND;
		     } else {
			if (element == " H") {
			   return HYDROGEN_GREY_BOND;
			}
		     }
		  }
	       }
	    }
	    return GREY_BOND;
	 } else {

	    if (bond_colour_type == coot::COLOUR_BY_CHAIN_C_ONLY) {
	       std::string element = at->element;

	       if (element == " C") {
		  if (atom_colour_map_p) {
		     int l_col = atom_colour_map_p->index_for_chain(std::string(at->GetChainID()));
		     return l_col;
		  } else {
		     std::cout << "ERROR:: Null atom_colour_map_p with COLOUR_BY_CHAIN_C_ONLY mode"
			       << std::endl;
		     return col;
		  } 
	       } else {
		  if (element == " N") {
		     return BLUE_BOND;
		  } else {
		     if (element == " O") {
			return RED_BOND;
		     } else {
			if (element == " S") {
			   return GREEN_BOND;
			} else {
			   if (element == " H") {
			      return HYDROGEN_GREY_BOND;
			   }
			}
		     }
		  }
	       }
	       return GREY_BOND;	       

	    } else { 

	       if (bond_colour_type == coot::DISULFIDE_COLOUR) {
		  return GREEN_BOND;
	       } else {
		  if (bond_colour_type == coot::COLOUR_BY_OCCUPANCY) {
		     if (at->occupancy > 0.95) {
			return BLUE_BOND;
		     } else {
			if (at->occupancy < 0.05) {
			   return RED_BOND;
			} else {
			   if (at->occupancy > 0.7) {
			      return CYAN_BOND;
			   } else {
			      if (at->occupancy > 0.45) {
				 return GREEN_BOND;
			      } else { 
				 if (at->occupancy > 0.25) {
				    return YELLOW_BOND;
				 } else {
				    return ORANGE_BOND;
				 }
			      }
			   }
			}
		     }
		  } else {
		     if (bond_colour_type == coot::COLOUR_BY_B_FACTOR) {
			float scaled_b = at->tempFactor*b_factor_scale;
			if (scaled_b < 10.0) {
			   return BLUE_BOND;
			} else {
			   if (scaled_b > 80.0) {
			      return RED_BOND;
			   } else {
			      if (scaled_b < 22.0) {
				 return CYAN_BOND;
			      } else {
				 if (scaled_b < 36.0) {
				    return GREEN_BOND;
				 } else {
				    if (scaled_b < 48.0) {
				       return YELLOW_BOND;
				    } else {
				       if (scaled_b < 62.0) {
					  return ORANGE_BOND;
				       }
				    }
				 }
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   return col;
}

void
Bond_lines_container::do_Ca_plus_ligands_bonds(atom_selection_container_t SelAtom,
					       float min_dist, float max_dist) { 
   //do_Ca_plus_ligands_bonds(SelAtom, min_dist, max_dist, coot::COLOUR_BY_ATOM_TYPE);
   do_Ca_plus_ligands_bonds(SelAtom, min_dist, max_dist, coot::COLOUR_BY_CHAIN);
}

void
Bond_lines_container::do_Ca_plus_ligands_bonds(atom_selection_container_t SelAtom,
					       float min_dist, float max_dist, 
					       int atom_colour_type) {

   std::cout << "do_Ca_plus_ligands_bonds with atom_colour_type "
	     << atom_colour_type << std::endl;
   
   CModel *model_p = SelAtom.mol->GetModel(1);
   if (model_p) {
      int istat;
      udd_has_ca_handle = SelAtom.mol->RegisterUDInteger (UDR_RESIDUE, "has CA");
      int nchains = model_p->GetNumberOfChains();
      CResidue *residue_p;
      CChain   *chain_p;
      for (int ichain=0; ichain<nchains; ichain++) { 
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    if (residue_p) { 
	       istat = residue_p->PutUDData(udd_has_ca_handle, 0);
	       if (istat == UDDATA_WrongUDRType) {
		  std::cout << "ERROR::  UDDATA_WrongUDRType in "
			    << "do_Ca_plus_ligands_bonds" << std::endl;
	       } 
	    }
	 }
      }

      coot::my_atom_colour_map_t acm;
      do_Ca_or_P_bonds_internal(SelAtom, " CA ", acm, min_dist, max_dist, atom_colour_type);

      // do_Ca_plus_ligands_bonds has set udd_has_ca_handle on
      // residues that have CAs.  Now let's run through the residues
      // again looking for residues that don't have this handle set.
      // We should do normal bonds for those residues (if they aren't
      // water, of course).

      std::vector<CAtom *> ligand_atoms;
      for (int ichain=0; ichain<nchains; ichain++) { 
	 chain_p = model_p->GetChain(ichain);
	 int nres = chain_p->GetNumberOfResidues();
	 int ic;
	 for (int ires=0; ires<nres; ires++) { 
	    residue_p = chain_p->GetResidue(ires);
	    if (residue_p) { 
	       if (residue_p->GetUDData(udd_has_ca_handle, ic) == UDDATA_Ok) { 
		  if (ic == 0) {
		     // Residue was not rendered as CA, needs normal bonds
		     std::string resname(residue_p->name);
		     if (resname != "WAT" && resname != "HOH") {
			int natoms = residue_p->GetNumberOfAtoms();
			for (int iat=0; iat<natoms; iat++) { 
			   ligand_atoms.push_back(residue_p->GetAtom(iat));
			}
		     }
		  }
	       }
	    }
	 }
      }
      if (ligand_atoms.size() > 0) { 
	 PCAtom *ligand_atoms_selection = new PCAtom[ligand_atoms.size()];
	 for(unsigned int iat=0; iat<ligand_atoms.size(); iat++) { 
	    ligand_atoms_selection[iat] = ligand_atoms[iat];
	 }
	 add_ligand_bonds(SelAtom, ligand_atoms_selection, ligand_atoms.size());
	 delete [] ligand_atoms_selection;
      } 
   }
}

void 
Bond_lines_container::do_normal_bonds_no_water(const atom_selection_container_t &asc_in,
					       float min_dist, 
					       float max_dist) {

   atom_selection_container_t asc = asc_in;

   // Now make a new atom selection that excludes WAT and HOH by using SKEY_XOR
   int newSelectionHandle = asc.mol->NewSelection();
   asc.SelectionHandle = newSelectionHandle;
   std::string solvent_res = "WAT,HOH";
   
   // We need to select all atoms here first, or crash when going back to all-atom view.
   asc.mol->SelectAtoms(asc.SelectionHandle, 0, "*",
			ANY_RES, "*",
			ANY_RES, "*",
			"*", "*", "*", "*");

   asc.mol->Select(asc.SelectionHandle, STYPE_ATOM, 0, "*",
		   ANY_RES, "*",
		   ANY_RES, "*",
		   (char *) solvent_res.c_str(), "*", "*", "*",
		   SKEY_XOR);

   asc.mol->GetSelIndex(asc.SelectionHandle, asc.atom_selection, asc.n_selected_atoms);
   // std::cout << "after water selection: n_selected_atoms: " << asc.n_selected_atoms << std::endl;
   construct_from_asc(asc, min_dist, max_dist, coot::COLOUR_BY_ATOM_TYPE, 0);
   asc.mol->DeleteSelection(asc.SelectionHandle);
}


int 
Bond_lines_container::add_ligand_bonds(const atom_selection_container_t &SelAtom, 
				       PPCAtom ligand_atoms_selection, 
				       int n_ligand_atoms) {

   int ibond = 0;
   atom_selection_container_t asc = SelAtom;
   asc.atom_selection = ligand_atoms_selection;
   asc.n_selected_atoms = n_ligand_atoms;
   construct_from_asc(asc, 0.01, 1.9, coot::COLOUR_BY_ATOM_TYPE, 0);

   return ibond;

}

void 
Bond_lines_container::do_colour_sec_struct_bonds(const atom_selection_container_t &asc,
						 float min_dist, float max_dist) { 

   if (asc.n_selected_atoms > 0) { 
      CModel *model_p = asc.mol->GetModel(1);
      model_p->CalcSecStructure(1);
      construct_from_asc(asc, 0.01, 1.9, coot::COLOUR_BY_SEC_STRUCT, 0);
   }
}
  

void 
Bond_lines_container::do_Ca_plus_ligands_colour_sec_struct_bonds(const atom_selection_container_t &asc,
								 float min_dist, float max_dist) { 
   if (asc.n_selected_atoms > 0) { 
      CModel *model_p = asc.mol->GetModel(1);
      int aminoSelHnd = -1;
#ifdef HAVE_MMDB_WITH_CISPEP
      model_p->CalcSecStructure(1, aminoSelHnd);
#else
      model_p->CalcSecStructure(1);
#endif // HAVE_MMDB_WITH_CISPEP      
      do_Ca_plus_ligands_bonds(asc, min_dist, max_dist, coot::COLOUR_BY_SEC_STRUCT);
   }
}



// Hah!  I'd written this so long ago that I'd forgotten that I'd done
// it.  I have reimplemented symmetry Calphas: addSymmetry_calphas().
//
// This function no longer works.
void
Bond_lines_container::do_symmetry_Ca_bonds(atom_selection_container_t SelAtom,
					   symm_trans_t symm_trans){

   Cell_Translation cell_trans;
   PPCAtom trans_ca_selection = NULL; // trans_sel(SelAtom, symm_trans);
   int n_ca;
   mat44 my_matt;
   int ncontacts;
   long i_contact_group = 1;
   PSContact contact = NULL;
   int col;

   // adjust my_matt to the symm_trans:
   //
   CMMDBCryst *cryst_p =  (CMMDBCryst *) &SelAtom.mol->get_cell();

   int err = cryst_p->GetTMatrix(my_matt, symm_trans.isym(), symm_trans.x(),
				 symm_trans.y(), symm_trans.z());

   if (err != 0) {
      cout << "!!!!!!!!!!!!!! something BAD with CMMDBCryst.GetTMatrix"
	   << endl;
   }
   
   int selHnd_ca = SelAtom.mol->NewSelection();

   SelAtom.mol->SelectAtoms(selHnd_ca, 0, "*",ANY_RES,"*",ANY_RES,"*",
			    "*"," CA ","C","*" );
   
   SelAtom.mol->GetSelIndex(selHnd_ca, trans_ca_selection, n_ca);

   SelAtom.mol->SeekContacts(trans_ca_selection, n_ca,
			     trans_ca_selection, n_ca,
			     0.01, 5.0, // min, max dist.
			     0,         // in same res also.
			     contact, ncontacts, 
			     0, &my_matt, i_contact_group);

   cout << "Found " << ncontacts/2 << " Ca-Ca links" << endl;

   if (ncontacts > 0) {
      for (int i=0; i< ncontacts; i++) {
	 if ( contact[i].id2 >  contact[i].id1 ) {

	    coot::Cartesian ca_1(trans_ca_selection[ contact[i].id1 ]->x,
			   trans_ca_selection[ contact[i].id1 ]->y,
			   trans_ca_selection[ contact[i].id1 ]->z);
 
	    coot::Cartesian ca_2(trans_ca_selection[ contact[i].id2 ]->x,
			   trans_ca_selection[ contact[i].id2 ]->y,
			   trans_ca_selection[ contact[i].id2 ]->z);

	    col = 0; // Carbon = yellow is the first color
	    addBond(col, ca_1, ca_2);
	 }
      }
   }
   delete [] contact;

}

// I add Colour by Segment at last (28 Oct 2003)
//
// 2001130 Now we pass the change_c_only_flag at the request of Phil.
// So if change_c_only_flag is 0, then we have a single colour for
// each chain.
//
// If change_c_only_flag is 1, then we want only Carbons of the chain
// to be coloured by chain.  The other atoms should be coloured by
// atom type.
// 
void
Bond_lines_container::do_colour_by_chain_bonds(const atom_selection_container_t &asc,
					       int draw_hydrogens_flag,
					       short int change_c_only_flag) {

   if (change_c_only_flag) {
      do_colour_by_chain_bonds_change_only(asc, draw_hydrogens_flag);
      return;
   }

   PSContact contact = NULL;
   int ncontacts;
   long i_contact_group = 1;

   float max_dist = 1.9;
   float min_dist = 0.01; // As in the constructor
			  // Bond_lines_container::Bond_lines_container(const
			  // atom_selection_container_t &SelAtom, int
			  // do_disulphide_bonds_in, int
			  // do_bonds_to_hydrogens_in)

   // matrix stuff
   mat44 my_matt;
   CSymOps symm;

   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;
      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;
   int col = 0; // atom (segment) colour

   int n_models = asc.mol->GetNumberOfModels();
   
   for (int imodel=1; imodel<=n_models; imodel++) { 

      PPCAtom atom_selection = 0;
      int n_selected_atoms = 0;
      contact = NULL;

      // make a new atom selection, based on the model.
      int SelectionHandle = asc.mol->NewSelection();
      asc.mol->SelectAtoms (SelectionHandle, imodel, "*",
			    ANY_RES, // starting resno, an int
			    "*", // any insertion code
			    ANY_RES, // ending resno
			    "*", // ending insertion code
			    "*", // any residue name
			    "*", // atom name
			    "*", // elements
			    "*"  // alt loc.
			    );
      
      asc.mol->GetSelIndex(SelectionHandle, atom_selection, n_selected_atoms);
      
      asc.mol->SeekContacts(atom_selection, n_selected_atoms,
			    atom_selection, n_selected_atoms,
			    min_dist, max_dist, // min, max distances
			    0,        // seqDist 0 -> in same res also
			    contact, ncontacts,
			    0, &my_matt, i_contact_group);
      
      coot::my_atom_colour_map_t atom_colour_map;
      int res1, res2;

      // Now, let's not forget that some atoms don't have contacts, so
      // whenever we find a contact for an atom, we mark it with
      // UserDefinedData "found bond".
      // 
      int uddHnd = asc.mol->RegisterUDInteger ( UDR_ATOM,"found bond" );
      if (uddHnd<0)  {
	 std::cout << " atom bonding registration failed.\n";
      } else {
	 for (int i=0; i<n_selected_atoms; i++)
	    atom_selection[i]->PutUDData(uddHnd,0);
      }

      if (ncontacts > 0) { 

	 CAtom *at1 = 0;
	 CAtom *at2 = 0;
	 std::string element1;
	 std::string element2;
      
	 for (int i=0; i< ncontacts; i++) {
	    if (contact[i].id2 > contact[i].id1) {

	       at1 = atom_selection[ contact[i].id1 ];
	       at2 = atom_selection[ contact[i].id2 ];
	    
	       res1 = at1->GetSeqNum();
	       res2 = at2->GetSeqNum();

	       if (abs(res1 - res2) < 2) { 

		  std::string segid1(at1->GetChainID());
		  std::string segid2(at2->GetChainID());
		  col = atom_colour_map.index_for_chain(segid1); 

		  if (segid1 == segid2) {

		     element1 = at1->element;
		     element2 = at2->element;
		     if ( (draw_hydrogens_flag == 1) ||
			  (element1 != " H" && element1 != " D" &&
			   element2 != " H" && element2 != " D") ) { 

			coot::Cartesian atom_1(at1->x, at1->y, at1->z);
			coot::Cartesian atom_2(at2->x, at2->y, at2->z);
		  
			// alternate location test
			// 	    
			std::string aloc_1(at1->altLoc);
			std::string aloc_2(at2->altLoc);
			// 
			if (aloc_1 == "" || aloc_2 == "" || aloc_1 == aloc_2) {
			   bonds_size_colour_check(col);
			   addBond(col, atom_1, atom_2);
		     
			   if (uddHnd>=0) {
			      at1->PutUDData(uddHnd,1);
			      at2->PutUDData(uddHnd,1);
			   }
			}
		     } else {
			// It was a hydrogen (or bonded to Hydrogen).
			// Mark it as bonded (we don't want to see single
			// unbonded (stared) hydorgens.
			if (uddHnd>=0) {
			   at1->PutUDData(uddHnd,1);
			   at2->PutUDData(uddHnd,1);
			}
		     } 
		  }
	       }
	    }
	 }
	 delete [] contact;
      }

      if (uddHnd>=0) {
    
	 // for atoms with no neighbour (contacts):
	 coot::Cartesian small_vec_x(0.5, 0.0, 0.0);
	 coot::Cartesian small_vec_y(0.0, 0.5, 0.0);
	 coot::Cartesian small_vec_z(0.0, 0.0, 0.5);

	 int ic; // changed by reference;
	 int col;
	 for (int i=0; i<n_selected_atoms; i++) { 
	    if ( atom_selection[i]->GetUDData(uddHnd, ic) == UDDATA_Ok ) {
	       if ((ic == 0) ||
		   (!strcmp(atom_selection[i]->element, " S")) ||
		   (!strcmp(atom_selection[i]->element, "SE")) ||
		   (!strcmp(atom_selection[i]->element, " P"))) {

		  std::string segid(atom_selection[i]->GetChainID());
		  col = atom_colour_map.index_for_chain(segid);

		  // 		  std::cout << " No contact for " << non_Hydrogen_atoms[i]
		  // 			    << std::endl;
	       
		  // no contact found or was Sulphur, or Phosphor

		  // So, was this a seleno-methione?
		  //
		  CResidue *atom_residue_p = atom_selection[i]->residue;
		  if (atom_residue_p) {
		     std::string resname = atom_selection[i]->GetResName();
		     if (1 &&
			 (resname == "MSE" || resname == "MET" || resname == "MSO"
			  || resname == "CYS" )) {
			handle_MET_or_MSE_case(atom_selection[i], uddHnd, col);
		     } else {
			std::string ele = atom_selection[i]->element;
			if (ele == "CL" || ele == "BR" || ele == " S" ||  ele == " I"
			    || ele == "Cl" || ele == "Br" 
			    || ele == "PT" 
			    || ele == "AS" || ele == " P" || ele == "AU" || ele == "HG") {
			   handle_long_bonded_atom(atom_selection[i], uddHnd, col);
			}
		     }
		  } else {
		     std::cout << "INFO:: trapped atom without residue in non-bonded atom check: "
			       << atom_selection[i] << std::endl;
		  }
	       }
	    }
	 }
      }
      asc.mol->DeleteSelection(SelectionHandle);
   }
   add_zero_occ_spots(asc);
   int atom_colour_type = coot::COLOUR_BY_CHAIN;
   add_atom_centres(asc, atom_colour_type);
}

void
Bond_lines_container::do_colour_by_chain_bonds_change_only(const atom_selection_container_t &asc,
							   int draw_hydrogens_flag) {

   // std::cout << "debug:: colour by chain, carbons only" << std::endl;
   float max_dist = 1.9;
   float min_dist = 0.01; // As in the constructor
			  // Bond_lines_container::Bond_lines_container(const
			  // atom_selection_container_t &SelAtom, int
			  // do_disulphide_bonds_in, int
			  // do_bonds_to_hydrogens_in)
   PSContact contact = NULL;
   int ncontacts;
   long i_contact_group = 1;

   // matrix stuff
   mat44 my_matt;
   CSymOps symm;
   int col = 0; // atom (segment) colour
   coot::my_atom_colour_map_t atom_colour_map;

   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;
      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

   int n_models = asc.mol->GetNumberOfModels();
   for (int imodel=1; imodel<=n_models; imodel++) { 

      PPCAtom atom_selection = 0;
      int n_selected_atoms = 0;
      contact = NULL;

      // make a new atom selection, based on the model.
      int SelectionHandle = asc.mol->NewSelection();
      asc.mol->SelectAtoms (SelectionHandle, imodel, "*",
			    ANY_RES, // starting resno, an int
			    "*", // any insertion code
			    ANY_RES, // ending resno
			    "*", // ending insertion code
			    "*", // any residue name
			    "*", // atom name
			    "*", // elements
			    "*"  // alt loc.
			    );
      
      asc.mol->GetSelIndex(SelectionHandle, atom_selection, n_selected_atoms);
      

      asc.mol->SeekContacts(atom_selection, n_selected_atoms,
			    atom_selection, n_selected_atoms,
			    min_dist, max_dist, // min, max distances
			    0,        // seqDist 0 -> in same res also
			    contact, ncontacts,
			    0, &my_matt, i_contact_group);

      int uddHnd = asc.mol->RegisterUDInteger (UDR_ATOM, "found bond");
      if (uddHnd<0)  {
	 std::cout << " atom bonding registration failed.\n";
      } else {
	 for (int i=0; i<n_selected_atoms; i++)
	    atom_selection[i]->PutUDData(uddHnd, 0);
      }

      if (ncontacts > 0) { 

	 CAtom *at1 = 0;
	 CAtom *at2 = 0;
	 std::string element1;
	 std::string element2;
	 int res1, res2;
	 int atom_colour_type = coot::COLOUR_BY_ATOM_TYPE;
      
	 for (int i=0; i< ncontacts; i++) {
	    if (contact[i].id2 > contact[i].id1) {

	       at1 = atom_selection[ contact[i].id1 ];
	       at2 = atom_selection[ contact[i].id2 ];
	    
	       res1 = at1->GetSeqNum();
	       res2 = at2->GetSeqNum();

	       if (abs(res1 - res2) < 2) { 

		  std::string segid1(at1->GetChainID());
		  std::string segid2(at2->GetChainID());
		  col = atom_colour_map.index_for_chain(segid1); 

		  if (segid1 == segid2) {

		     element1 = at1->element;
		     element2 = at2->element;
		     if ( (draw_hydrogens_flag == 1) ||
			  (element1 != " H" && element1 != " D" &&
			   element2 != " H" && element2 != " D") ) { 

			coot::Cartesian atom_1(at1->x, at1->y, at1->z);
			coot::Cartesian atom_2(at2->x, at2->y, at2->z);
		  
			// alternate location test
			// 	    
			std::string aloc_1(at1->altLoc);
			std::string aloc_2(at2->altLoc);
			// 
			if (aloc_1 == "" || aloc_2 == "" || aloc_1 == aloc_2) {

			   if (element1 != element2) {
			   
			      // Bonded to different atom elements.
			      //
			      coot::Cartesian bond_mid_point = atom_1.mid_point(atom_2);
			      if (element1 != " C") {
				 col = atom_colour(atom_selection[ contact[i].id1 ], atom_colour_type);
				 bonds_size_colour_check(col);
				 addBond(col, atom_1, bond_mid_point);
			      } else {
				 bonds_size_colour_check(col);
				 addBond(col, atom_1, bond_mid_point);
			      }
			      if (element2 != " C") {
				 col = atom_colour(atom_selection[ contact[i].id2 ], atom_colour_type);
				 bonds_size_colour_check(col);
				 addBond(col, atom_2, bond_mid_point);
			      } else {
				 bonds_size_colour_check(col);
				 addBond(col, atom_2, bond_mid_point);
			      }
			   
			   } else { 

			      if (element1 == " C") { 
				 bonds_size_colour_check(col);
				 addBond(col, atom_1, atom_2);
			      } else {
				 col = atom_colour(atom_selection[ contact[i].id1 ], atom_colour_type);
				 bonds_size_colour_check(col);
				 addBond(col, atom_1, atom_2);
			      } 
			   }

			   // we drew a bond.  Mark it up.
			   if (uddHnd>=0) {
			      at1->PutUDData(uddHnd,1);
			      at2->PutUDData(uddHnd,1);
			   }
			}
		     } else {
			// It was a hydrogen (or bonded to Hydrogen).
			// Mark it as bonded (we don't want to see single
			// unbonded (stared) hydorgens.
			if (uddHnd>=0) {
			   at1->PutUDData(uddHnd,1);
			   at2->PutUDData(uddHnd,1);
			}
		     } 
		  }
	       }
	    }
	 }
	 delete [] contact;
      }


      if (uddHnd>=0) {
    
	 // for atoms with no neighbour (contacts):
	 coot::Cartesian small_vec_x(0.5, 0.0, 0.0);
	 coot::Cartesian small_vec_y(0.0, 0.5, 0.0);
	 coot::Cartesian small_vec_z(0.0, 0.0, 0.5);

	 int ic; // changed by reference;
	 int col;
	 int atom_colour_type = coot::COLOUR_BY_CHAIN_C_ONLY; // the atoms connected to the SE (say) are Cs.
	 for (int i=0; i<n_selected_atoms; i++) { 
	    if ( atom_selection[i]->GetUDData(uddHnd, ic) == UDDATA_Ok ) {
	       if ((ic == 0) ||
		   (!strcmp(atom_selection[i]->element, " S")) ||
		   (!strcmp(atom_selection[i]->element, "SE")) ||
		   (!strcmp(atom_selection[i]->element, " P"))) {

		  std::string segid(atom_selection[i]->GetChainID());
		  col = atom_colour_map.index_for_chain(segid);

		  // 		  std::cout << " No contact for " << non_Hydrogen_atoms[i]
		  // 			    << std::endl;
	       
		  // no contact found or was Sulphur, or Phosphor

		  // So, was this a seleno-methione?
		  //
		  CResidue *atom_residue_p = atom_selection[i]->residue;
		  if (atom_residue_p) {
		     std::string resname = atom_selection[i]->GetResName();
		     if (1 &&
			 (resname == "MSE" || resname == "MET" || resname == "MSO"
			  || resname == "CYS" )) {
			handle_MET_or_MSE_case(atom_selection[i], uddHnd, atom_colour_type,
					       &atom_colour_map);
		     } else {
			std::string ele = atom_selection[i]->element;
			if (ele == "CL" || ele == "BR" || ele == " S" ||  ele == " I"
			    || ele == "Cl" || ele == "Br" 
			    || ele == "PT" 
			    || ele == "AS" || ele == " P" || ele == "AU" || ele == "HG") {
			   handle_long_bonded_atom(atom_selection[i], uddHnd, atom_colour_type);
			}
		     }
		  } else {
		     std::cout << "INFO:: trapped atom without residue in non-bonded atom check: "
			       << atom_selection[i] << std::endl;
		  }
	       }
	    }
	 }
      }

      asc.mol->DeleteSelection(SelectionHandle);
   }
   add_zero_occ_spots(asc);
   int atom_colour_type = coot::COLOUR_BY_CHAIN;
   add_atom_centres(asc, atom_colour_type);
}

void
Bond_lines_container::do_colour_by_molecule_bonds(const atom_selection_container_t &asc,
						  int draw_hydrogens_flag) { 

   PSContact contact = NULL;
   int ncontacts;
   long i_contact_group = 1;

   float max_dist = 1.9;
   float min_dist = 0.01; // As in the constructor
			  // Bond_lines_container::Bond_lines_container(const
			  // atom_selection_container_t &SelAtom, int
			  // do_disulphide_bonds_in, int
			  // do_bonds_to_hydrogens_in)

   // matrix stuff
   mat44 my_matt;
   CSymOps symm;

   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;
      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;
   int col = 0; // atom (segment) colour


   int n_models = asc.mol->GetNumberOfModels();
   
   for (int imodel=1; imodel<=n_models; imodel++) { 

      PPCAtom atom_selection = 0;
      int n_selected_atoms = 0;
      contact = NULL;

      // make a new atom selection, based on the model.
      int SelectionHandle = asc.mol->NewSelection();
      asc.mol->SelectAtoms (SelectionHandle, imodel, "*",
			    ANY_RES, // starting resno, an int
			    "*", // any insertion code
			    ANY_RES, // ending resno
			    "*", // ending insertion code
			    "*", // any residue name
			    "*", // atom name
			    "*", // elements
			    "*"  // alt loc.
			    );
      
      asc.mol->GetSelIndex(SelectionHandle, atom_selection, n_selected_atoms);
      
      asc.mol->SeekContacts(atom_selection, n_selected_atoms,
			    atom_selection, n_selected_atoms,
			    min_dist, max_dist, // min, max distances
			    0,        // seqDist 0 -> in same res also
			    contact, ncontacts,
			    0, &my_matt, i_contact_group);

      coot::my_atom_colour_map_t atom_colour_map;
      int res1, res2;
      // Now, let's not forget that some atoms don't have contacts, so
      // whenever we find a contact for an atom, we mark it with
      // UserDefinedData "found bond".
      // 
      int uddHnd = asc.mol->RegisterUDInteger (UDR_ATOM, "found bond");
      if (uddHnd<0)  {
	 std::cout << " atom bonding registration failed.\n";
      } else {
	 for (int i=0; i<n_selected_atoms; i++)
	    atom_selection[i]->PutUDData(uddHnd,0);
      }

      if (ncontacts > 0) { 
      
	 std::string element1;
	 std::string element2;
	 for (int i=0; i< ncontacts; i++) {
	    if (contact[i].id2 > contact[i].id1) {

	       CAtom *at1 = atom_selection[ contact[i].id1 ];
	       CAtom *at2 = atom_selection[ contact[i].id2 ];
      
	       res1 = at1->GetSeqNum();
	       res2 = at2->GetSeqNum();

	       if (abs(res1 - res2) < 2) { 
		  coot::Cartesian atom_1(at1->x, at1->y, at1->z);
		  coot::Cartesian atom_2(at2->x, at2->y, at2->z);

		  element1 = at1->element;
		  element2 = at2->element;
		  if ( (draw_hydrogens_flag == 1) ||
		       (element1 != " H" && element1 != " D" &&
			element2 != " H" && element2 != " D") ) { 
		  
		     // alternate location test
		     // 	    
		     std::string aloc_1(at1->altLoc);
		     std::string aloc_2(at2->altLoc);
		     // 
		     if (aloc_1 == "" || aloc_2 == "" || aloc_1 == aloc_2) {
			bonds_size_colour_check(col);
			addBond(col, atom_1, atom_2);
		     
			if (uddHnd>=0) {
			   at1->PutUDData(uddHnd,1);
			   at2->PutUDData(uddHnd,1);
			}
		     }
		  } else {
		     // It was a hydrogen (or bonded to Hydrogen).
		     // Mark it as bonded (we don't want to see single
		     // unbonded (stared) hydorgens.
		     if (uddHnd>=0) {
			at1->PutUDData(uddHnd,1);
			at2->PutUDData(uddHnd,1);
		     }
		  }
	       }
	    }
	 }
	 delete [] contact;
	 contact = NULL;
      }

      if (uddHnd>=0) {
    
	 // for atoms with no neighbour (contacts):
	 coot::Cartesian small_vec_x(0.5, 0.0, 0.0);
	 coot::Cartesian small_vec_y(0.0, 0.5, 0.0);
	 coot::Cartesian small_vec_z(0.0, 0.0, 0.5);

	 int ic; // changed by reference;
	 int col;
	 for (int i=0; i<n_selected_atoms; i++) { 
	    if ( atom_selection[i]->GetUDData(uddHnd,ic) == UDDATA_Ok ) {
	       if (ic == 0) {
		  std::string segid(atom_selection[i]->GetChainID());
		  col = atom_colour_map.index_for_chain(segid);
		  bonds_size_colour_check(col);
		  coot::Cartesian atom(atom_selection[i]->x,
				       atom_selection[i]->y,
				       atom_selection[i]->z);
	       
		  addBond(col, atom+small_vec_x, atom-small_vec_x); 
		  addBond(col, atom+small_vec_y, atom-small_vec_y);
		  addBond(col, atom+small_vec_z, atom-small_vec_z);
	       }
	    }
	 }
      }
      asc.mol->DeleteSelection(SelectionHandle);
   }
   add_zero_occ_spots(asc);
   int atom_colour_type = coot::COLOUR_BY_CHAIN;
   add_atom_centres(asc, atom_colour_type);
}


void
Bond_lines_container::add_zero_occ_spots(const atom_selection_container_t &SelAtom) {

   zero_occ_spot.resize(0);

   for (int i=0; i<SelAtom.n_selected_atoms; i++) { 
      if (SelAtom.atom_selection[i]->occupancy < 0.01 &&
	  SelAtom.atom_selection[i]->occupancy > -1) { // shelx occ test
	 // we don't want to see atoms with occupancy -61 from a shelx ins
	 // file with zero occupancy spots.
	 zero_occ_spot.push_back(coot::Cartesian(SelAtom.atom_selection[i]->x,
						 SelAtom.atom_selection[i]->y,
						 SelAtom.atom_selection[i]->z));
      }
   }
}

void
Bond_lines_container::add_atom_centres(const atom_selection_container_t &SelAtom,
				       int atom_colour_type) {

   atom_centres.clear();
   atom_centres_colour.clear();

   for (int i=0; i<SelAtom.n_selected_atoms; i++) {
      bool is_H_flag = 0;
      if (std::string(SelAtom.atom_selection[i]->element) == " H")
	 is_H_flag = 1;
      if (do_bonds_to_hydrogens || (do_bonds_to_hydrogens == 0 && (!is_H_flag))) {
	 std::pair<bool, coot::Cartesian> p(is_H_flag,
					    coot::Cartesian(SelAtom.atom_selection[i]->x,
							    SelAtom.atom_selection[i]->y,
							    SelAtom.atom_selection[i]->z));
	 atom_centres.push_back(p);
	 atom_centres_colour.push_back(atom_colour(SelAtom.atom_selection[i], atom_colour_type));
      }
   }
}




void
graphical_bonds_container::add_zero_occ_spots(const std::vector<coot::Cartesian> &spots) { 

   n_zero_occ_spot = spots.size();

   if (n_zero_occ_spot > 0) {
      zero_occ_spot = new coot::Cartesian[n_zero_occ_spot];
      for (int j=0; j<n_zero_occ_spot; j++) { 
	 zero_occ_spot[j] = spots[j];
      }
   }
}

void
graphical_bonds_container::add_atom_centres(const std::vector<std::pair<bool,coot::Cartesian> > &centres,
					    const std::vector<int> &colours) {

   if (colours.size() != centres.size()) {
      std::cout << "ERROR!! colours.size() != centres.size() in add_atom_centres\n";
   }
   n_atom_centres_ = centres.size();
   atom_centres_ = new std::pair<bool, coot::Cartesian>[n_atom_centres_];
   atom_centres_colour_ = new int[n_atom_centres_];
   for (int i=0; i<n_atom_centres_; i++) {
      atom_centres_[i] = centres[i];
      atom_centres_colour_[i] = colours[i];
   }
}


