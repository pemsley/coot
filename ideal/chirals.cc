/* ideal/chirals.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 by The University of York
 * Copyright 2008, 2009, 2010  by The University of Oxford
 * Copyright 2013 by Medical Research Council
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
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

#include <string.h> // for strcmp


#include <algorithm> // for sort
#include <stdexcept>

#include "compat/coot-sysdep.h"

//
#include "coot-utils/coot-coord-extras.hh"  // is_nucleotide_by_dict

#include "simple-restraint.hh"
#include "utils/logging.hh"
extern logging logger;

// #include "mmdb.h" // for printing of mmdb::Atom pointers as info not raw
                     // pointers.  Removed. Too much (linking issues in)
                     // Makefile pain.



// Before we ask this question, we need to match residues which chiral
// restraints.
//
// i.e we need a function that takes a mmdb::Manager * and returns
// pairs of simple_restraint, mmdb::Residue, that is done in
//  molecule_class_info-other.cc
// 
std::vector<std::pair<short int, coot::atom_spec_t> >
coot::is_inverted_chiral_atom_p(const coot::dict_chiral_restraint_t &chiral_restraint,
				mmdb::Residue *res) {

   short int ibad=0;
   coot::atom_spec_t chiral_atom;
   std::vector<std::pair<short int, coot::atom_spec_t> > v;

   int i_no_res_atoms = res->GetNumberOfAtoms();
   
   for (int iat1=0; iat1<i_no_res_atoms; iat1++) {
      std::string pdb_atom_name1(res->GetAtom(iat1)->name);
      if (pdb_atom_name1 == chiral_restraint.atom_id_1_4c()) {
	 
	 for (int iat2=0; iat2<i_no_res_atoms; iat2++) {
	    std::string pdb_atom_name2(res->GetAtom(iat2)->name);
	    if (pdb_atom_name2 == chiral_restraint.atom_id_2_4c()) {
	       
	       for (int iat3=0; iat3<i_no_res_atoms; iat3++) {
		  std::string pdb_atom_name3(res->GetAtom(iat3)->name);
		  if (pdb_atom_name3 == chiral_restraint.atom_id_3_4c()) {
		     
		     for (int iatc=0; iatc<i_no_res_atoms; iatc++) {
			std::string pdb_atom_namec(res->GetAtom(iatc)->name);
			if (pdb_atom_namec == chiral_restraint.atom_id_c_4c()) {

			   // Now, do they have corresponding
			   // altconfs?  The bonded atoms must have
			   // the same altconf as each other (or blank).
			   
			   // 
			   // e.g. for Chiral atom CB in THR: CA   CB,A OG1,A CG2,A (1)
			   //                                 CA   CB   OG1,A CG2,A (2)
			   //                 but not:        CA   CB   OG1,A CG2,B (3)

			   std::string chiral_alt_conf(res->GetAtom(iatc)->altLoc);
			   std::string altLoc1(res->GetAtom(iat1)->altLoc);
			   std::string altLoc2(res->GetAtom(iat2)->altLoc);
			   std::string altLoc3(res->GetAtom(iat3)->altLoc);
			   short int matching_altlocs = 0;

			   // These should catch most cases, 
			   // but there may be some that they don't.

			   if ( chiral_alt_conf != "" &&  
			       (altLoc1 == "" || altLoc1 == chiral_alt_conf) &&
			       (altLoc2 == "" || altLoc2 == chiral_alt_conf) &&
			       (altLoc3 == "" || altLoc3 == chiral_alt_conf) ) { 
			      matching_altlocs = 1;
			   } 
			   if ( chiral_alt_conf == "" &&  
			       (altLoc1 == "" ) &&
			       (altLoc2 == "" ) &&
			       (altLoc3 == "" ) ) { 
			      matching_altlocs = 1;
			   }

			   // We need to check like this to get rid of case (3):
			   if ( chiral_alt_conf == "" &&  
			       (altLoc1 == altLoc2 || altLoc1 == "" || altLoc2 == "") &&
			       (altLoc1 == altLoc3 || altLoc1 == "" || altLoc3 == "") &&
			       (altLoc2 == altLoc3 || altLoc2 == "" || altLoc3 == "") ) { 
			      matching_altlocs = 1;
			   }

			   if (matching_altlocs) { 

			      clipper::Coord_orth centre(res->GetAtom(iatc)->x,
							 res->GetAtom(iatc)->y,
							 res->GetAtom(iatc)->z);
			      clipper::Coord_orth a1(res->GetAtom(iat1)->x,
						     res->GetAtom(iat1)->y,
						     res->GetAtom(iat1)->z);
			      clipper::Coord_orth a2(res->GetAtom(iat2)->x,
						     res->GetAtom(iat2)->y,
						     res->GetAtom(iat2)->z);
			      clipper::Coord_orth a3(res->GetAtom(iat3)->x,
						     res->GetAtom(iat3)->y,
						     res->GetAtom(iat3)->z);
			   
			      clipper::Coord_orth a = a1 - centre;
			      clipper::Coord_orth b = a2 - centre;
			      clipper::Coord_orth c = a3 - centre;
			      double cv = clipper::Coord_orth::dot(a, clipper::Coord_orth::cross(b,c));

			      chiral_atom = atom_spec_t(res->GetChainID(),
                                                        res->GetSeqNum(),
                                                        res->GetInsCode(),
                                                        res->GetAtom(iatc)->name,
                                                        res->GetAtom(iatc)->altLoc);

			      if (cv*chiral_restraint.volume_sign < 0) {
// 				 std::cout << "DEBUG:: " << res->name << " "
// 					   << res->GetSeqNum() << " " << pdb_atom_namec
// 					   << " has bad chiral volume\n";
				 ibad = 1;
			      } else { 
// 				 std::cout << "DEBUG:: " << res->name << " "
// 					   << res->GetSeqNum() << " " << pdb_atom_namec
// 					   << " has good chiral volume\n";
				 ibad = 0;
			      }
			      v.push_back(std::pair<short int, coot::atom_spec_t> (ibad, chiral_atom));
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   
   // return std::pair<short int, coot::atom_spec_t> (ibad, chiral_atom);
   return v;
}


// Return a list of bad chiral volumes for this molecule:
// 
// Return also a flag for the status of this test, were there any
// residues for which we didn't find restraints?  The flag is the
// number of residue names in first part of the returned pair.
// 0 is good.
// 
std::pair<std::vector<std::string> , std::vector <coot::atom_spec_t> >
coot::inverted_chiral_volumes(int imol,
			      mmdb::Manager *mol, protein_geometry *geom_p,
			      int cif_dictionary_read_number) {

   std::vector <coot::atom_spec_t> v;
   int restraints_status = 1;
   std::vector<std::string> unknown_types_vec; 

   int n_models = mol->GetNumberOfModels();
   for (int imod=1; imod<=n_models; imod++) { 
      
      mmdb::Model *model_p = mol->GetModel(imod);
      mmdb::Chain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      if (nchains <= 0) { 
	 std::cout << "bad nchains in molecule " << nchains
		   << std::endl;
      } else {
	 for (int ichain=0; ichain<nchains; ichain++) {
	    chain_p = model_p->GetChain(ichain);
	    if (chain_p == NULL) {  
	       // This should not be necessary. It seem to be a
	       // result of mmdb corruption elsewhere - possibly
	       // DeleteChain in update_molecule_to().
	       std::cout << "NULL chain in ... " << std::endl;
	    } else { 
	       int nres = chain_p->GetNumberOfResidues();
	       mmdb::PResidue residue_p;
	       // mmdb::Atom *at;
	       for (int ires=0; ires<nres; ires++) { 
		  residue_p = chain_p->GetResidue(ires);
		  int n_atoms = residue_p->GetNumberOfAtoms();
		  if (n_atoms > 3) {
		     std::string residue_type(residue_p->name);
		     if (residue_type == "UNK")
			residue_type = "ALA";
		     if (! geom_p->have_dictionary_for_residue_type(residue_type,
								    imol,
								    cif_dictionary_read_number)) {
			std::cout << "WARNING::! Failed to find restraint for residue type "
				  << residue_type << std::endl;
			restraints_status = 0;
			// only add the residue_type if it is not already in the vector
			short int irfound = 0;
			for (unsigned int ir=0; ir<unknown_types_vec.size(); ir++)
			   if (unknown_types_vec[ir] == residue_type) {
			      irfound = 1;
			      break;
			   }
			if (irfound == 0) 
			   unknown_types_vec.push_back(residue_type);
		     } else { 
			std::vector<coot::dict_chiral_restraint_t> chiral_restraints = 
			   geom_p->get_monomer_chiral_volumes(std::string(residue_p->name), imol);
			coot::dict_chiral_restraint_t chiral_restraint;
			for (unsigned int irestr=0; irestr<chiral_restraints.size(); irestr++) { 
			   chiral_restraint = chiral_restraints[irestr];
			   if (! chiral_restraint.is_a_both_restraint()) {
			      std::vector<std::pair<short int, coot::atom_spec_t> > c = 
				 coot::is_inverted_chiral_atom_p(chiral_restraint, residue_p);
			      for (unsigned int ibad=0; ibad<c.size(); ibad++) { 
				 if (c[ibad].first) {
                                    if (false) // I don't need to see this.
                                       std::cout << "INFO:: found bad chiral atom: "
                                                 << chain_p->GetChainID() << " "
                                                 << residue_p->GetSeqNum() << " "
                                                 << residue_p->GetInsCode() << " "
                                                 << chiral_restraint.atom_id_c_4c() << " "
                                                 << c[ibad].second.alt_conf << std::endl;

				    v.push_back(coot::atom_spec_t(chain_p->GetChainID(),
								  residue_p->GetSeqNum(),
								  residue_p->GetInsCode(),
								  chiral_restraint.atom_id_c_4c(),
								  c[ibad].second.alt_conf));
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
   std::pair<std::vector<std::string>, std::vector<coot::atom_spec_t> > pair(unknown_types_vec, v);
   return pair;
}


			     
void
coot::fix_chiral_atom_maybe (const simple_restraint &chiral_restraint,
			     gsl_vector *v) {
   int idx = 3*(chiral_restraint.atom_index_centre);
   clipper::Coord_orth centre(gsl_vector_get(v, idx),
			      gsl_vector_get(v, idx+1),
			      gsl_vector_get(v, idx+2));

   idx = 3*(chiral_restraint.atom_index_1);
   clipper::Coord_orth a1(gsl_vector_get(v, idx),
			      gsl_vector_get(v, idx+1),
			      gsl_vector_get(v, idx+2));
   idx = 3*(chiral_restraint.atom_index_2);
   clipper::Coord_orth a2(gsl_vector_get(v, idx),
			      gsl_vector_get(v, idx+1),
			      gsl_vector_get(v, idx+2));
   idx = 3*(chiral_restraint.atom_index_3);
   clipper::Coord_orth a3(gsl_vector_get(v, idx),
			      gsl_vector_get(v, idx+1),
			      gsl_vector_get(v, idx+2));

   clipper::Coord_orth a = a1 - centre;
   clipper::Coord_orth b = a2 - centre;
   clipper::Coord_orth c = a3 - centre;

   double cv = clipper::Coord_orth::dot(a, clipper::Coord_orth::cross(b,c));
   double volume_sign = chiral_restraint.chiral_volume_sign;

   std::cout << "DEBUG:::::::: Fix chiral maybe :::::: "
	     << volume_sign*cv << std::endl;

   if (volume_sign*cv < 0.0) {
      // we have a mismatch between restraint and actual, needs fixing
      std::cout << "Atom index " << chiral_restraint.atom_index_centre
		<< " is undergoing chiral centre inversion\n";
      fix_chiral_atom_internal(chiral_restraint, v);
   }
}

void
coot::fix_chiral_atom_internal(const simple_restraint &chiral_restraint,
			       gsl_vector *v) {

   // we need to find the chiral plane (i.e. the plane of the
   // non-centre atoms of the chiral restraint).
   //
   // We need to find the distance of the chiral atom from the chiral
   // plane (the vector is normal to the plane - and the plane
   // equation defines the normal).
   //
   // We will then move the atom to the other side of the plane.
   //

   //Recall d = Ax+By+Cz-D, where d is the distance from the
   //plane of the chiral atom.
   //

   // So we have 3 points, a b and c.  The cross-product b-a x c-a
   // gives us the normal the the plane.

   int idx = 3*(chiral_restraint.atom_index_centre);
   clipper::Coord_orth centre(gsl_vector_get(v, idx),
			      gsl_vector_get(v, idx+1),
			      gsl_vector_get(v, idx+2));

   idx = 3*(chiral_restraint.atom_index_1);
   clipper::Coord_orth a1(gsl_vector_get(v, idx),
			  gsl_vector_get(v, idx+1),
			  gsl_vector_get(v, idx+2));
   idx = 3*(chiral_restraint.atom_index_2);
   clipper::Coord_orth a2(gsl_vector_get(v, idx),
			  gsl_vector_get(v, idx+1),
			  gsl_vector_get(v, idx+2));
   idx = 3*(chiral_restraint.atom_index_3);
   clipper::Coord_orth a3(gsl_vector_get(v, idx),
			  gsl_vector_get(v, idx+1),
			  gsl_vector_get(v, idx+2));

   clipper::Coord_orth a1_2 = a1 - a2;
   clipper::Coord_orth a2_3 = a2 - a3;
   clipper::Coord_orth a3_1 = a3 - a1;

//     clipper::Coord_orth normal =
//        clipper::Coord_orth(b_a.y()*c_a.z() - b_a.z()*c_a.y(),
//  			  b_a.z()*c_a.x() - b_a.x()*c_a.z(),
//  			  b_a.x()*c_a.y() - b_a.y()*c_a.x());

   clipper::Coord_orth scaled_normal =
      clipper::Coord_orth(a1.y()*a2_3.z() + a2.y()*a3_1.z() + a3.y()*a1_2.z(),
			  a1.z()*a2_3.x() + a2.z()*a3_1.x() + a3.z()*a1_2.x(),
			  a1.x()*a2_3.y() + a2.x()*a3_1.y() + a3.x()*a1_2.y());
   double scaled_D =   (
			a1.x()*(a2.y()*a3.z()-a3.y()*a2.z()) + 
		        a2.x()*(a3.y()*a1.z()-a1.y()*a3.z()) +
			a3.x()*(a1.y()*a2.z()-a2.y()*a1.z()));

   // normalize:
   double len_s_n = sqrt(scaled_normal.lengthsq());
   double one_over_len_sn = 1.0/len_s_n;

   clipper::Coord_orth normal = one_over_len_sn * scaled_normal;
   double D = one_over_len_sn * scaled_D;

   std::cout << "normal now: " << normal.format() << "D: " << D << "\n";
   
   //d = Ax+By+Cz-D
   double d = normal.x()*centre.x() + normal.y()*centre.y() + normal.z()*centre.z() - D;

   std::cout << "d is " << d << " for atom index "
	     << chiral_restraint.atom_index_centre << "\n";

   // Cool, now lets move the centre atom to the other side of the plane
   // by 0.5A

   double factor;
   if (d<0)
      factor = d - 0.5;
   else
      factor = d + 0.5;
   clipper::Coord_orth shift = - factor*normal;

   std::cout  << "DEBUG::  CHIRAL: shifting atom index "
	      << chiral_restraint.atom_index_centre << " by "
	      << shift.format() << "\n";

   idx = 3*(chiral_restraint.atom_index_centre);
   gsl_vector_set(v, idx,   gsl_vector_get(v, idx  ) + shift.x());
   gsl_vector_set(v, idx+1, gsl_vector_get(v, idx+1) + shift.y());
   gsl_vector_set(v, idx+2, gsl_vector_get(v, idx+2) + shift.z());
}

void
coot::restraints_container_t::fix_chiral_atoms_maybe(gsl_vector *s) {

   if (restraints_usage_flag & coot::CHIRAL_VOLUME_MASK) {
      for(int i=0; i<size(); i++) {
	 {
	    const simple_restraint &rest = restraints_vec[i];
	    if ( restraints_vec[i].restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) {
	       fix_chiral_atom_maybe(restraints_vec[i], s);
	    }
	 }
      }
   }
}

bool
coot::restraints_container_t::check_pushable_chiral_hydrogens(gsl_vector *v) {

   bool state = 0; // none
   if (restraints_usage_flag & coot::CHIRAL_VOLUME_MASK) {
      for (int i=0; i<size(); i++) {
	       {
	          const simple_restraint &rest = restraints_vec[i];
	          if ( rest.restraint_type == coot::CHIRAL_VOLUME_RESTRAINT) {
	             if (rest.chiral_hydrogen_index != coot::UNSET_INDEX) {
		              // so we have a single H attached to this chiral centre.

		              // is the hydrogen on the wrong side of the chiral centre?
		              //
		              bool val = chiral_hydrogen_needs_pushing(rest, v);
		              // std::cout << "::::  chiral_hydrogen_needs_pushing() returned " << val << std::endl;
		              if (val) {
		                 push_chiral_hydrogen(rest, v);
		                 state = 1;
		                 break; // only do one at a time.
		              }
	             }
	          }
	       }
      }
   }
   return state;
} 

// a single H on this chiral centre is on the wrong side of the
// chiral centre?
// (not sure that this is the best name for this function)
//
// Note: chiral_restraint.chiral_hydrogen_index must be valid!
// 
bool
coot::restraints_container_t::chiral_hydrogen_needs_pushing(const coot::simple_restraint &chiral_restraint,
							    const gsl_vector *v) const {

   int n_bad_angles = 0; 

   int n_angles = 0;
   double angle_distortion = 0;
   double angle_distortion_max = 35; // if more than this, then the hydrogen is pushable.

   // we do nothing if the chiral hydrogen is part of the definition
   // of the chiral centre (because the pushing routine determines the
   // position of the hydrogen from the position of the non-hydrogen
   // atoms connected to the chiral centre).
   // 
   if (chiral_restraint.chiral_hydrogen_index != coot::UNSET_INDEX) {
      if (chiral_restraint.chiral_hydrogen_index == chiral_restraint.atom_index_1 ||
	  chiral_restraint.chiral_hydrogen_index == chiral_restraint.atom_index_2 ||
	  chiral_restraint.chiral_hydrogen_index == chiral_restraint.atom_index_2) {
	 return 0;
      }
   } 
   
      
   // now check, does this have an inverted chiral centre?
   //
   // bool icc = has_inverted_chiral_centre(chiral_restraint, v);
   
   bool tiny_cv = has_tiny_chiral_centre_volume(chiral_restraint, v);

   if (! tiny_cv)
      return 0; // don't push the chiral hydrogen.


   // OK, so perhaps the chiral volume was just the right side of 0,
   // so let's check the angles

   // second check, are there highly distorted angle terms?
   // 
   for (int i=0; i<size(); i++) {
      if (restraints_usage_flag & coot::ANGLES_MASK) {
	 if ( (*this)[i].restraint_type == coot::ANGLE_RESTRAINT) {
	    const coot::simple_restraint &restraint = (*this)[i];

	    // now select the angle terms involving the hydrogen
	    bool flag = 0;
	    if (chiral_restraint.atom_index_centre == restraint.atom_index_2) {
	       if (chiral_restraint.chiral_hydrogen_index == restraint.atom_index_1 ||
		   chiral_restraint.chiral_hydrogen_index == restraint.atom_index_3) {
		  flag = 1;
	       }
	    }

	    if (flag) { // i.e. angle term for hydrogen
	       angle_distortion = coot::distortion_score_angle(restraint, v);
	       if (angle_distortion > angle_distortion_max) { 
		  std::cout << "::angle distortion for restraint " << i << ":  "
			    << angle_distortion << std::endl;
		  n_bad_angles++;
	       }
	    }
	 }
      }
   }

   // consider the 2 positions:
   // 
   // 1) at the end of refinement, the chiral volume signs match (barely) and we have
   //    distorted angles
   //
   // 2) if we have a starting structure that was correct for the initial (other)
   //   chirality (that is, before we flipped the chiral restraint on the atom), then we
   //   start off with clearly incorrect chiral volume sign (but the angles are not bad).
   //   We don't want to flip the H in that case.

   bool needs_pushing = 0;

   if ((n_bad_angles > 1) && tiny_cv)
      needs_pushing = 1;

   if (0)
      std::cout << ".... chiral_hydrogen_needs_pushing() restraint " << chiral_restraint 
		<< " returns " << needs_pushing
		<< " with " << n_bad_angles << " bad angles and tiny_cv " << tiny_cv << std::endl;
   return needs_pushing;
} 

bool
coot::restraints_container_t::has_inverted_chiral_centre(const coot::simple_restraint &chiral_restraint,
							 const gsl_vector *v) const {

   bool r = 0;
   int idx = 3*(chiral_restraint.atom_index_centre);
   clipper::Coord_orth centre(gsl_vector_get(v, idx),
			      gsl_vector_get(v, idx+1),
			      gsl_vector_get(v, idx+2));
   
   idx = 3*(chiral_restraint.atom_index_1);
   clipper::Coord_orth a1(gsl_vector_get(v, idx),
			  gsl_vector_get(v, idx+1),
			  gsl_vector_get(v, idx+2));
   idx = 3*(chiral_restraint.atom_index_2);
   clipper::Coord_orth a2(gsl_vector_get(v, idx),
			  gsl_vector_get(v, idx+1),
			  gsl_vector_get(v, idx+2));
   idx = 3*(chiral_restraint.atom_index_3);
   clipper::Coord_orth a3(gsl_vector_get(v, idx),
			  gsl_vector_get(v, idx+1),
			  gsl_vector_get(v, idx+2));
   
   clipper::Coord_orth a = a1 - centre;
   clipper::Coord_orth b = a2 - centre;
   clipper::Coord_orth c = a3 - centre;
   double cv = clipper::Coord_orth::dot(a, clipper::Coord_orth::cross(b,c));
   
   double distortion = cv - chiral_restraint.target_chiral_volume;
   double distort_2  = distortion * distortion / (chiral_restraint.sigma * chiral_restraint.sigma);

   // either both positve or both negative is fine...
   if (cv*chiral_restraint.chiral_volume_sign < 0)
      r = true; // yes, a problem

   if (false) {
      std::cout << "debug in has_inverted_chiral_centre() cv = " << cv << " target: " 
		<< chiral_restraint.target_chiral_volume << std::endl;
      std::cout << "debug in has_inverted_chiral_centre() distortion = " << distortion << " " 
		<< ", sigma = " << chiral_restraint.sigma << " distort_2 = " << distort_2 << std::endl;
      std::cout << "debug in has_inverted_chiral_centre() result: " << r << std::endl;
   }
   return r;

}

bool
coot::restraints_container_t::has_tiny_chiral_centre_volume(const coot::simple_restraint &chiral_restraint,
							    const gsl_vector *v) const {

   bool r = false;
   int idx = 3*(chiral_restraint.atom_index_centre);
   clipper::Coord_orth centre(gsl_vector_get(v, idx),
			      gsl_vector_get(v, idx+1),
			      gsl_vector_get(v, idx+2));
   
   idx = 3*(chiral_restraint.atom_index_1);
   clipper::Coord_orth a1(gsl_vector_get(v, idx),
			  gsl_vector_get(v, idx+1),
			  gsl_vector_get(v, idx+2));
   idx = 3*(chiral_restraint.atom_index_2);
   clipper::Coord_orth a2(gsl_vector_get(v, idx),
			  gsl_vector_get(v, idx+1),
			  gsl_vector_get(v, idx+2));
   idx = 3*(chiral_restraint.atom_index_3);
   clipper::Coord_orth a3(gsl_vector_get(v, idx),
			  gsl_vector_get(v, idx+1),
			  gsl_vector_get(v, idx+2));
   
   clipper::Coord_orth a = a1 - centre;
   clipper::Coord_orth b = a2 - centre;
   clipper::Coord_orth c = a3 - centre;
   double cv = clipper::Coord_orth::dot(a, clipper::Coord_orth::cross(b,c));
   
   // double distortion = cv - chiral_restraint.target_chiral_volume;
   // double distort_2  = distortion * distortion / (chiral_restraint.sigma * chiral_restraint.sigma);

   double chiral_fraction = fabs(cv/chiral_restraint.target_chiral_volume);

   // std::cout << "     chiral_fraction " << chiral_fraction << std::endl; // typically 1

   if (chiral_fraction < 0.4)
      r = true; 
   
   return r;

}

void
coot::restraints_container_t::push_chiral_hydrogen(const simple_restraint &chiral_restraint, gsl_vector *v) {
   
   bool r = 0;
   if (chiral_restraint.chiral_hydrogen_index != coot::UNSET_INDEX) {

      int idx = 3*(chiral_restraint.atom_index_centre);
      clipper::Coord_orth centre(gsl_vector_get(v, idx),
				 gsl_vector_get(v, idx+1),
				 gsl_vector_get(v, idx+2));
   
      idx = 3*(chiral_restraint.atom_index_1);
      clipper::Coord_orth a1(gsl_vector_get(v, idx),
			     gsl_vector_get(v, idx+1),
			     gsl_vector_get(v, idx+2));
      idx = 3*(chiral_restraint.atom_index_2);
      clipper::Coord_orth a2(gsl_vector_get(v, idx),
			     gsl_vector_get(v, idx+1),
			     gsl_vector_get(v, idx+2));
      idx = 3*(chiral_restraint.atom_index_3);
      clipper::Coord_orth a3(gsl_vector_get(v, idx),
			     gsl_vector_get(v, idx+1),
			     gsl_vector_get(v, idx+2));

      idx = 3*(chiral_restraint.chiral_hydrogen_index);
      clipper::Coord_orth h_current_pos(gsl_vector_get(v, idx),
					gsl_vector_get(v, idx+1),
					gsl_vector_get(v, idx+2));
      
      clipper::Coord_orth a = a1 - centre;
      clipper::Coord_orth b = a2 - centre;
      clipper::Coord_orth c = a3 - centre;

      clipper::Coord_orth mid_point_sum = a1 + a2 + a3;
      clipper::Coord_orth mid_point(mid_point_sum.x()/3.0, mid_point_sum.y()/3.0, mid_point_sum.z()/3.0);

      clipper::Coord_orth dv_u(clipper::Coord_orth(centre - mid_point).unit());

      clipper::Coord_orth new_h_pos = centre + 1.09 * dv_u;

      std::cout << "::INFO pushing H "
		<< coot::atom_spec_t(atom[chiral_restraint.chiral_hydrogen_index]) 
		<< " on " << coot::atom_spec_t(atom[chiral_restraint.atom_index_centre])
		<< " from " << h_current_pos.format()
		<< " to " << new_h_pos.format() << std::endl;
      
      idx = 3*chiral_restraint.chiral_hydrogen_index;
      gsl_vector_set(v, idx,     new_h_pos.x());
      gsl_vector_set(v, idx + 1, new_h_pos.y());
      gsl_vector_set(v, idx + 2, new_h_pos.z());
   }
   
}

bool
coot::restraints_container_t::check_through_ring_bonds(gsl_vector *v) {

   // actually this only checks for a very long bond (distortion)

   bool status = false;
   for (int i=0; i<size(); i++) {
      if (restraints_usage_flag & coot::BONDS_MASK) { 
	 if ( (*this)[i].restraint_type == coot::BOND_RESTRAINT) {
	    bool p = bond_is_very_long((*this)[i], v);
	    if (p) {
// 	       std::cout << "    restraint " << i << " " << (*this)[i]
// 			 << " has a very long bond " << std::endl;
	    }
	    
	 }
      }
   }
   return status;
}


bool
coot::restraints_container_t::bond_is_very_long(const coot::simple_restraint &bond_restraint,
						const gsl_vector *v) const {
   bool status = false;
   int idx = 3*(bond_restraint.atom_index_1);
   clipper::Coord_orth a1(gsl_vector_get(v, idx),
			  gsl_vector_get(v, idx+1),
			  gsl_vector_get(v, idx+2));
   idx = 3*(bond_restraint.atom_index_2);
   clipper::Coord_orth a2(gsl_vector_get(v, idx),
			  gsl_vector_get(v, idx+1),
			  gsl_vector_get(v, idx+2));
   double l = clipper::Coord_orth::length(a1,a2);
   double bit = l - bond_restraint.target_value;

   if (bit > 1.0) {
      // std::cout << "long bond  bit: " << bit << std::endl;
      status = true;
   }
   

   return status;
}



// make restraints and get distortions. chiral_volume_limit_for_outlier
// should/might be about 2.0.
// the second is the chiral atom and the distortion score
std::pair<std::vector<std::string> , std::vector <std::pair<coot::atom_spec_t, double> > >
coot::distorted_chiral_volumes(int imol, mmdb::Manager *mol, protein_geometry *geom_p,
                               int cif_dictionary_read_number,
                               double chiral_volume_limit_for_outlier) {

   auto make_local_residues = [] (mmdb::Manager *mol) {

      std::vector<std::pair<bool,mmdb::Residue *> > lr;
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  lr.push_back(std::make_pair(false, residue_p));
               }
            }
         }
      }
      return lr;
   };

   std::vector<std::string> missing_types;
   std::vector<std::pair<atom_spec_t, double> > d_specs;

   std::vector<std::pair<bool,mmdb::Residue *> > local_residues = make_local_residues(mol);
   std::vector<mmdb::Link> links;
   std::vector<coot::atom_spec_t> fixed_atom_specs;
   clipper::Xmap<float> dummy_map;
   restraints_container_t restraints(local_residues, links, *geom_p, mol, fixed_atom_specs, &dummy_map);
   restraint_usage_Flags flags = coot::CHIRAL_VOLUMES;
   bool do_trans_peptide_restraints = false;
   bool do_link_restraints = false;
   bool do_flank_restraints = false;
   pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
   int n_threads = 1;
   ctpl::thread_pool tp(n_threads);
   restraints.thread_pool(&tp, n_threads);

   // this is making nbc restraints I think - stop it.
   int n_restraints = restraints.make_restraints(imol, *geom_p, flags, 1, do_trans_peptide_restraints,
                                                0.0, false, true, true, false, pseudos, do_link_restraints, do_flank_restraints);

   std::cout << "--------------------- distorted_chiral_volumes() made " << n_restraints << " restraints" << std::endl;

   if (n_restraints > 0) {
      geometry_distortion_info_container_t gdic = restraints.geometric_distortions();
      for (std::size_t id=0; id<gdic.geometry_distortion.size(); id++) {
         const auto &rest = gdic.geometry_distortion[id];
         if (rest.distortion_score > chiral_volume_limit_for_outlier) {
            if (rest.atom_indices.size() == 4) {
               mmdb:: Atom *at = restraints.get_atom(rest.atom_indices[0]); // the chiral centre of course
               if (at)
                  d_specs.push_back(std::make_pair(atom_spec_t(at), rest.distortion_score));
            }
         }
      }
   }
   return std::make_pair(missing_types, d_specs);
}
