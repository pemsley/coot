/* analysis/cablam.cc
 * 
 * Copyright 2013 by the Medical Research Council
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

#include <map>
#include <clipper/core/coords.h>

#include "coot-utils/coot-coord-utils.hh"
#include "cablam.hh"

clipper::Coord_orth
coot::cablam::get_closest_CA_CA_approach(const coot::torsion_atom_quad &quad) const {

   // this is not the function I need.  I need
   //
   // get_closest_CA_CA_approach(CA_pos_p, CA_pos_t, O_pos_p);

   clipper::Coord_orth CA_p(quad.atom_1->x, quad.atom_1->y, quad.atom_1->z);
   clipper::Coord_orth CA_t(quad.atom_2->x, quad.atom_2->y, quad.atom_2->z);
   clipper::Coord_orth  O_t(quad.atom_4->x, quad.atom_4->y, quad.atom_4->z);

   clipper::Coord_orth PT = CA_t - CA_p;
   clipper::Coord_orth PT_unit(PT.unit());
   clipper::Coord_orth PO =  O_t - CA_p;

   double PO_length = sqrt(PO.lengthsq());
   
   double cos_alpha_p = (PT * PO)/(PO_length * sqrt(PT.lengthsq()));

   double PC_length = cos_alpha_p * PO_length;
   clipper::Coord_orth Closest_approach = CA_p + PC_length * PT_unit;
   return Closest_approach;
} 

clipper::Coord_orth
coot::cablam::get_closest_CA_CA_approach(const clipper::Coord_orth &CA_pos_p,
					 const clipper::Coord_orth &CA_pos_t,
					 const clipper::Coord_orth &O_pos_p) const {

   clipper::Coord_orth PT = CA_pos_t - CA_pos_p;
   clipper::Coord_orth PO = O_pos_p - CA_pos_p;
   clipper::Coord_orth PT_unit(PT.unit());

   double PO_length = sqrt(PO.lengthsq());
   double PT_length = sqrt(PT.lengthsq());
   double cos_alpha = clipper::Coord_orth::dot(PT, PO)/(PO_length * PT_length);

   double PC_length = cos_alpha * PO_length;
   
   clipper::Coord_orth Closest_approach = CA_pos_p + PC_length * PT_unit;
   return Closest_approach;

}


coot::cablam::cablam(mmdb::PResidue *residues, int n_sel_residues) {

   // there are 4 angles/torsions of interest:
   //
   // nu:    torsion O(i-1) cl_ap(i-1) cl_ap(i) O(i)
   // mu_in: torsion CA(i-2) CA(i-1) CA(i) CA(i+1)
   // mu_out torsion CA(i-1) CA(i) CA(i+1) CA(i+2)
   // angle CA(i-2) CA(i-1) CA(i)
   // where
   // cl_ap(i-1) (aka X(i-i)) is the closest_approach of the O to the CA(i-1)-CA(i) vector
   // 

   std::map<mmdb::Residue *, torsion_atom_quad> residue_CA_CA_quads;
   std::map<mmdb::Residue *, torsion_atom_quad> residue_OX_XO_quads;

   for (int ires=2; ires<n_sel_residues; ires++) {

      if ((ires+2) < n_sel_residues) {

	 // we dont calculate nu if we cant calculate mu_out - heyho
	 
	 mmdb::Residue *res_minus_2 = residues[ires-2];
	 mmdb::Residue *res_minus_1 = residues[ires-1];
	 mmdb::Residue *res_0       = residues[ires];
	 mmdb::Residue *res_plus_1  = residues[ires+1];
	 mmdb::Residue *res_plus_2  = residues[ires+2];

	 // mu_in:  CA(i-2) CA(i-1) CA(i) CA(i+1) 
	 // mu_out: CA(i) CA(i-1) CA(i) CA(i+1) CA(i+2)

	 // test that we are in the same chain
	 // test that N_res_t is close to C_res_p
	 // test that residue numbers are consecutive.

	 bool proceed = false;
	 if ((res_minus_2->GetSeqNum()+1) == res_minus_1->GetSeqNum()) {
	    if ((res_minus_1->GetSeqNum()+1) == res_0->GetSeqNum()) {
	       if ((res_0->GetSeqNum()+1) == res_plus_1->GetSeqNum()) {
		  if ((res_plus_1->GetSeqNum()+1) == res_plus_2->GetSeqNum()) {

		     if (res_minus_2->GetChain() == res_plus_2->GetChain()) {
			if (res_minus_2->GetChain() == res_plus_2->GetChain()) {
			   if (res_0->GetChain() == res_plus_2->GetChain()) {
			      if (res_plus_1->GetChain() == res_plus_2->GetChain()) {
				 proceed = true;
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }

	 if (proceed) {

	    bool proceed_2 = false;
	    mmdb::Atom *CA_m_2 = res_minus_2->GetAtom(" CA ");
	    mmdb::Atom *CA_m_1 = res_minus_1->GetAtom(" CA ");
	    mmdb::Atom *CA_0   = res_0->GetAtom(" CA ");
	    mmdb::Atom *CA_p_1 = res_plus_1->GetAtom(" CA ");
	    mmdb::Atom *CA_p_2 = res_plus_2->GetAtom(" CA ");

	    if (CA_m_2 && CA_m_1 && CA_0 && CA_p_1 && CA_p_2) { 

	       std::string ac_m_2 = CA_m_2->altLoc;
	       std::string ac_m_1 = CA_m_1->altLoc;
	       std::string ac_0   = CA_0->altLoc;
	       std::string ac_p_1 = CA_p_1->altLoc;
	       std::string ac_p_2 = CA_p_2->altLoc;

	       if (ac_m_2 == "") { 
		  if (ac_m_1 == "") { 
		     if (ac_0 == "") { 
			if (ac_p_1 == "") { 
			   if (ac_p_2 == "") {
			      proceed_2 = true;
			   }
			}
		     }
		  }
	       }
	    }

	    if (proceed_2) { 

	       // check N and C are close
	       mmdb::Atom *C_m_2 = res_minus_2->GetAtom(" C  ");
	       mmdb::Atom *C_m_1 = res_minus_1->GetAtom(" C  ");
	       mmdb::Atom *C_0   = res_0->GetAtom(" C  ");
	       mmdb::Atom *C_p_1 = res_plus_1->GetAtom(" C  ");
	       mmdb::Atom *C_p_2 = res_plus_2->GetAtom(" C  ");
	       mmdb::Atom *N_m_2 = res_minus_2->GetAtom(" N  ");
	       mmdb::Atom *N_m_1 = res_minus_1->GetAtom(" N  ");
	       mmdb::Atom *N_0   = res_0->GetAtom(" N  ");
	       mmdb::Atom *N_p_1 = res_plus_1->GetAtom(" N  ");
	       mmdb::Atom *N_p_2 = res_plus_2->GetAtom(" N  ");

	       clipper::Coord_orth co_C_m_2 = coot::co(C_m_2);
	       clipper::Coord_orth co_C_m_1 = coot::co(C_m_1);
	       clipper::Coord_orth co_C_0   = coot::co(C_0);
	       clipper::Coord_orth co_C_p_1 = coot::co(C_p_1);
	       clipper::Coord_orth co_C_p_2 = coot::co(C_p_2);

	       clipper::Coord_orth co_N_m_2 = coot::co(N_m_2);
	       clipper::Coord_orth co_N_m_1 = coot::co(N_m_1);
	       clipper::Coord_orth co_N_0   = coot::co(N_0);
	       clipper::Coord_orth co_N_p_1 = coot::co(N_p_1);
	       clipper::Coord_orth co_N_p_2 = coot::co(N_p_2);

	       double l_1 = clipper::Coord_orth::length(co_C_m_2, co_N_m_1);
	       double l_2 = clipper::Coord_orth::length(co_C_m_1, co_N_0);
	       double l_3 = clipper::Coord_orth::length(co_C_0,   co_N_p_1);
	       double l_4 = clipper::Coord_orth::length(co_C_p_1, co_N_p_2);

	       bool proceed_3 = false;
	       if ((l_1 > 1.3) && (l_1 < 1.4)) {
		  if ((l_2 > 1.3) && (l_2 < 1.4)) {
		     if ((l_3 > 1.3) && (l_3 < 1.4)) {
			if ((l_4 > 1.3) && (l_4 < 1.4)) {
			   proceed_3 = true;
			}
		     }
		  }
	       }

	       if (proceed_3) {

		  mmdb::Atom *O_m_1 = res_minus_1->GetAtom(" O  ");
		  mmdb::Atom *O_0   = res_0->GetAtom(" O  ");

		  if (O_m_1 && O_0) {
	       
		     coot::atom_quad mu_in (CA_m_2, CA_m_1, CA_0,   CA_p_1);
		     coot::atom_quad mu_out(CA_m_1, CA_0,   CA_p_1, CA_p_2);

		     clipper::Coord_orth co_O_m_1 = coot::co(O_m_1);
		     clipper::Coord_orth co_O_0   = coot::co(O_0);

		     clipper::Coord_orth co_CA_m_1 = coot::co(CA_m_1);
		     clipper::Coord_orth co_CA_0   = coot::co(CA_0);
		     clipper::Coord_orth co_CA_p_1 = coot::co(CA_p_1);

		     clipper::Coord_orth cl_ap_1 = get_closest_CA_CA_approach(co_CA_m_1, co_CA_0,   co_O_m_1);
		     clipper::Coord_orth cl_ap_2 = get_closest_CA_CA_approach(co_CA_0,   co_CA_p_1, co_O_0);

		     double nu = clipper::Coord_orth::torsion(co_O_m_1, cl_ap_1, cl_ap_2, co_O_0);

		     std::cout << atom_spec_t(CA_0) << " " << mu_in.torsion() << " " << mu_out.torsion()
			       << " " << nu * 180/M_PI << std::endl;

		  }
	       }
	    }
	 }
      }
   }
}
