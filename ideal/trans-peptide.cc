/* ideal/trans-peptide.cc
 * 
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

#include "simple-restraint.hh"
int
coot::restraints_container_t::add_link_trans_peptide(mmdb::Residue *first,
						     mmdb::Residue *second,
						     short int is_fixed_first,
						     short int is_fixed_second,
						     const coot::protein_geometry &geom) {
   int n_torsion = 0;

   mmdb::PPAtom first_sel;
   mmdb::PPAtom second_sel;
   int n_first_res_atoms, n_second_res_atoms;
   int n_atom_1, n_atom_2, n_atom_3, n_atom_4;
   
   first->GetAtomTable(first_sel,   n_first_res_atoms); 
   second->GetAtomTable(second_sel, n_second_res_atoms);

   // assigned to either first_sel or second_sel when atom_1_comp_id
   // (etc.) are known.
   mmdb::PPAtom atom_1_sel, atom_2_sel, atom_3_sel, atom_4_sel;

   if (n_first_res_atoms <= 0) {
      std::cout << "no atoms in first residue!? " << std::endl;
   }
   if (n_second_res_atoms <= 0) {
      std::cout << "no atoms in second residue!? " << std::endl;
   }

   std::vector<bool> fixed_flags(4);
   fixed_flags[0] = 0;  // not fixed
   fixed_flags[1] = 0;
   fixed_flags[2] = 0;
   fixed_flags[3] = 0;

   atom_1_sel = first_sel;
   atom_2_sel = first_sel;
   atom_3_sel = second_sel;
   atom_4_sel = second_sel;

   n_atom_1 = n_first_res_atoms;
   n_atom_2 = n_first_res_atoms;
   n_atom_3 = n_second_res_atoms;
   n_atom_4 = n_second_res_atoms;
   fixed_flags[0] = is_fixed_first;
   fixed_flags[1] = is_fixed_first;
   fixed_flags[2] = is_fixed_second;
   fixed_flags[3] = is_fixed_second;

   for (int ifat=0; ifat<n_atom_1; ifat++) { 
      std::string pdb_atom_name_1(atom_1_sel[ifat]->name);
	       
      if (pdb_atom_name_1 == " CA ") {
	 for (int isat=0; isat<n_atom_2; isat++) { 
	    std::string pdb_atom_name_2(atom_2_sel[isat]->name);
		     
	    if (pdb_atom_name_2 == " C  ") {
	       for (int itat=0; itat<n_atom_3; itat++) { 
		  std::string pdb_atom_name_3(atom_3_sel[itat]->name);
			   
		  if (pdb_atom_name_3 == " N  ") {
		     for (int iffat=0; iffat<n_atom_4; iffat++) {
			std::string pdb_atom_name_4(atom_4_sel[iffat]->name);
				 			   
			if (pdb_atom_name_4 == " CA ") {
				    
			   int index1 = get_asc_index(atom_1_sel[ifat]->name,
						      atom_1_sel[ifat]->altLoc,
						      atom_1_sel[ifat]->residue->seqNum,
						      atom_1_sel[ifat]->GetInsCode(),
						      atom_1_sel[ifat]->GetChainID());
			
			   int index2 = get_asc_index(atom_2_sel[isat]->name,
						      atom_2_sel[isat]->altLoc,
						      atom_2_sel[isat]->residue->seqNum,
						      atom_2_sel[isat]->GetInsCode(),
						      atom_2_sel[isat]->GetChainID());
				    
			   int index3 = get_asc_index(atom_3_sel[itat]->name,
						      atom_3_sel[itat]->altLoc,
						      atom_3_sel[itat]->residue->seqNum,
						      atom_3_sel[itat]->GetInsCode(),
						      atom_3_sel[itat]->GetChainID());

			   int index4 = get_asc_index(atom_4_sel[iffat]->name,
						      atom_4_sel[iffat]->altLoc,
						      atom_4_sel[iffat]->residue->seqNum,
						      atom_4_sel[iffat]->GetInsCode(),
						      atom_4_sel[iffat]->residue->GetChainID());

			   if (false) 
			      std::cout << "trans-peptide restraint.... "
					<< " from atoms \n    "
					<< atom_1_sel[ifat]->name << " " 
					<< atom_1_sel[ifat]->GetSeqNum() << "\n    " 
					<< atom_2_sel[isat]->name << " " 
					<< atom_2_sel[isat]->GetSeqNum() << "\n    " 
					<< atom_3_sel[itat]->name << " " 
					<< atom_3_sel[itat]->GetSeqNum() << "\n    " 
					<< atom_4_sel[iffat]->name << " " 
					<< atom_4_sel[iffat]->GetSeqNum() << "\n";

			   // if the angle is currently trans, the we should add a
			   // trans peptide restraint.
				    
			   clipper::Coord_orth a1 = co(atom_1_sel[ifat]);
			   clipper::Coord_orth a2 = co(atom_2_sel[isat]);
			   clipper::Coord_orth a3 = co(atom_3_sel[itat]);
			   clipper::Coord_orth a4 = co(atom_4_sel[iffat]);

			   double omega = clipper::Coord_orth::torsion(a1, a2, a3, a4);
			   double dist = clipper::Coord_orth::length(a2,a3);

			   if (false)
			      std::cout << "... current omega " << omega * 180/M_PI
					<< " and N-C dist " << dist << std::endl;
				 
			   if (dist < 2.0) {
			      if ((omega > 0.5 * M_PI) || (omega < -0.5 * M_PI)) {

				 std::vector<bool> other_fixed_flags = make_fixed_flags(index1, index2, index3, index4);
				 for (unsigned int ii=0; ii<other_fixed_flags.size(); ii++)
				    if (other_fixed_flags[ii])
				       fixed_flags[ii] = true;
				 double target_omega = 180.0;
				 double esd = 2.0; // 5.0 lets slip 72A in 2bmd to trans
				 add(TRANS_PEPTIDE_RESTRAINT, index1, index2, index3, index4,
				     fixed_flags,
				     target_omega,
				     esd,
				     1.2, // dummy value
				     1);
				 n_torsion++;
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
   return n_torsion; 
}

//
// Return the distortion score from a single torsion restraint.
// 
// can throw a std::runtime_error if there is a problem calculating the torsion.
// 
double
coot::distortion_score_trans_peptide(const coot::simple_restraint &restraint,
				     const gsl_vector *v) {

   // First calculate the torsion:
   // theta = arctan(E/G); 
   // where E = a.(bxc) and G = -a.c + (a.b)(b.c)

   int idx; 

   idx = 3*(restraint.atom_index_1);
   clipper::Coord_orth P1(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(restraint.atom_index_2); 
   clipper::Coord_orth P2(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(restraint.atom_index_3); 
   clipper::Coord_orth P3(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));
   idx = 3*(restraint.atom_index_4); 
   clipper::Coord_orth P4(gsl_vector_get(v,idx), 
			  gsl_vector_get(v,idx+1), 
			  gsl_vector_get(v,idx+2));

   clipper::Coord_orth a = P2 - P1; 
   clipper::Coord_orth b = P3 - P2; 
   clipper::Coord_orth c = P4 - P3;

   double al = sqrt(clipper::Coord_orth::dot(a,a));
   double bl = sqrt(clipper::Coord_orth::dot(b,b));
   double cl = sqrt(clipper::Coord_orth::dot(c,c));
   double cos_a1 = clipper::Coord_orth::dot(a,b)/(al*bl);
   double cos_a2 = clipper::Coord_orth::dot(b,c)/(bl*cl);

   // instabilty when the P2-P3-P4 or P1-P2-p3 line is linear. Give up with the derivatives
   // similar escape in the derivatives
   //
   if (cos_a1 > 0.9 || cos_a2> 0.9) {
      return 0;
   } else {

      // b*b * [ a.(bxc)/b ]
      double E = clipper::Coord_orth::dot(a,clipper::Coord_orth::cross(b,c)) *
	 sqrt( b.lengthsq() );

      // b*b * [ -a.c+(a.b)(b.c)/(b*b) ] = -a.c*b*b + (a.b)(b.c)
      double G = - clipper::Coord_orth::dot(a,c)*b.lengthsq()
	 + clipper::Coord_orth::dot(a,b)*clipper::Coord_orth::dot(b,c);

      double theta = clipper::Util::rad2d(atan2(E,G));
      if (false) {
	 if (clipper::Util::isnan(theta)) {
	    std::string mess = "WARNING: distortion_score_torsion() observed torsion theta is a NAN!";
	    throw std::runtime_error(mess);
	 }
      }

      double diff = theta - restraint.target_value;

      // the target is 180, so if theta is -179, then we want
      // to add 360
      // if theta is -110, we want to add 360, -> 250 - diff is 70

      if (diff > 180)
	 diff -= 360;
      if (diff < -180)
	 diff += 360;

      if (false) { // debug
	 double pen = diff*diff/(restraint.sigma * restraint.sigma);
	 std::cout << "in distortion_trans_peptide theta (calc): " << theta
		   << " periodicity " << restraint.periodicity
		   << " target "      << restraint.target_value
		   << " diff: " << diff << " ";
	 std::cout << "sigma= " << restraint.sigma
		   << " weight= " << pow(restraint.sigma,-2.0) << " ";
	 std::cout << "score " << pen;
	 if (false) {
	    std::cout << " " << P1.format();
	    std::cout << " " << P2.format();
	    std::cout << " " << P3.format();
	    std::cout << " " << P4.format();
	 }
	 std::cout << "\n";
      }
      return diff*diff/(restraint.sigma * restraint.sigma);
   }
}

// Add in the torsion gradients
//
void coot::my_df_trans_peptides(const gsl_vector *v, 
				void *params, 
				gsl_vector *df) {
   
   int n_trans_peptide_restr = 0; 
   int idx; 

   // first extract the object from params 
   //
   restraints_container_t *restraints = (restraints_container_t *)params;

   if (restraints->restraints_usage_flag & TRANS_PEPTIDE_MASK) { 

      int restraints_size = restraints->size();
      for (int i=0; i< restraints_size; i++) {

	 const simple_restraint &restraint = restraints->at(i);
      
	 if (restraint.restraint_type == TRANS_PEPTIDE_RESTRAINT) {

	    n_trans_peptide_restr++;

	    idx = 3*(restraint.atom_index_1); 
	    clipper::Coord_orth P1(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));
	    idx = 3*(restraint.atom_index_2); 
	    clipper::Coord_orth P2(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));
	    idx = 3*(restraint.atom_index_3); 
	    clipper::Coord_orth P3(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));
	    idx = 3*(restraint.atom_index_4); 
	    clipper::Coord_orth P4(gsl_vector_get(v,idx), 
				   gsl_vector_get(v,idx+1), 
				   gsl_vector_get(v,idx+2));

	    try {

	       // if the bond angles are (near) linear, then the distortion gradients
	       // are zero.
	       //
	       distortion_torsion_gradients_t dtg =
		  fill_distortion_torsion_gradients(P1, P2, P3, P4);

	       if (! dtg.zero_gradients) {

		  double diff = dtg.theta - restraint.target_value;
		  // because trans restraints - 180, (see distortion score notes)
		  if (diff > 180)
		     diff -= 360;
		  if (diff < -180)
		     diff += 360;

		  if (false)
		     std::cout << "in df_trans_peptide: dtg.theta is " << dtg.theta 
			       <<  " and target is " << restraint.target_value 
			       << " and diff is " << diff 
			       << " and periodicity: " << restraint.periodicity << std::endl;

		  double tt = tan(clipper::Util::d2rad(dtg.theta));
		  double trans_peptide_scale = (1.0/(1+tt*tt)) *
		     clipper::Util::rad2d(1.0);

		  double weight = 1/(restraint.sigma * restraint.sigma);

		  // 	       std::cout << "trans_peptide weight: " << weight << std::endl;
		  // 	       std::cout << "trans_peptide_scale : " << trans_peptide_scale << std::endl; 
		  // 	       std::cout << "diff          : " << trans_peptide_scale << std::endl; 	       

		  double xP1_contrib = 2.0*diff*dtg.dD_dxP1*trans_peptide_scale * weight;
		  double xP2_contrib = 2.0*diff*dtg.dD_dxP2*trans_peptide_scale * weight;
		  double xP3_contrib = 2.0*diff*dtg.dD_dxP3*trans_peptide_scale * weight;
		  double xP4_contrib = 2.0*diff*dtg.dD_dxP4*trans_peptide_scale * weight;

		  double yP1_contrib = 2.0*diff*dtg.dD_dyP1*trans_peptide_scale * weight;
		  double yP2_contrib = 2.0*diff*dtg.dD_dyP2*trans_peptide_scale * weight;
		  double yP3_contrib = 2.0*diff*dtg.dD_dyP3*trans_peptide_scale * weight;
		  double yP4_contrib = 2.0*diff*dtg.dD_dyP4*trans_peptide_scale * weight;

		  double zP1_contrib = 2.0*diff*dtg.dD_dzP1*trans_peptide_scale * weight;
		  double zP2_contrib = 2.0*diff*dtg.dD_dzP2*trans_peptide_scale * weight;
		  double zP3_contrib = 2.0*diff*dtg.dD_dzP3*trans_peptide_scale * weight;
		  double zP4_contrib = 2.0*diff*dtg.dD_dzP4*trans_peptide_scale * weight;

		  if (false) {
		     std::cout << "trans-peptide " << dtg.theta << " derivs: "
			       << dtg.dD_dxP1 << " " << dtg.dD_dyP1 << " " << dtg.dD_dzP1 << " "
			       << dtg.dD_dxP2 << " " << dtg.dD_dyP2 << " " << dtg.dD_dzP2 << " "
			       << dtg.dD_dxP3 << " " << dtg.dD_dyP3 << " " << dtg.dD_dzP3 << " "
			       << dtg.dD_dxP4 << " " << dtg.dD_dyP4 << " " << dtg.dD_dzP4 << " "
			       << dtg.dD_dxP1 << " " << dtg.dD_dyP1 << " " << dtg.dD_dzP1 << " "
			       << std::endl;
		  }
		  if (false) {
		     std::cout << "trans-peptide " << dtg.theta << " delta: "
			       << xP1_contrib << " " << yP1_contrib << " " << zP1_contrib << " "
			       << xP2_contrib << " " << yP2_contrib << " " << zP2_contrib << " "
			       << xP3_contrib << " " << yP3_contrib << " " << zP3_contrib << " "
			       << xP4_contrib << " " << yP4_contrib << " " << zP4_contrib << " "
			       << " prior-dfs "
			       << gsl_vector_get(df, 3*restraint.atom_index_1) << " "
			       << gsl_vector_get(df, 3*restraint.atom_index_1+1) << " "
			       << gsl_vector_get(df, 3*restraint.atom_index_1+2) << " "
			       << gsl_vector_get(df, 3*restraint.atom_index_2) << " "
			       << gsl_vector_get(df, 3*restraint.atom_index_2+1) << " "
			       << gsl_vector_get(df, 3*restraint.atom_index_2+2) << " "
			       << gsl_vector_get(df, 3*restraint.atom_index_3) << " "
			       << gsl_vector_get(df, 3*restraint.atom_index_3+1) << " "
			       << gsl_vector_get(df, 3*restraint.atom_index_3+2) << " "
			       << gsl_vector_get(df, 3*restraint.atom_index_4) << " "
			       << gsl_vector_get(df, 3*restraint.atom_index_4+1) << " "
			       << gsl_vector_get(df, 3*restraint.atom_index_4+2) << " "
			       << std::endl;
		  }

		  if (! restraint.fixed_atom_flags[0]) { 
		     idx = 3*(restraint.atom_index_1);
		     *gsl_vector_ptr(df, idx  ) += xP1_contrib;
		     *gsl_vector_ptr(df, idx+1) += yP1_contrib;
		     *gsl_vector_ptr(df, idx+2) += zP1_contrib;
		  }

		  if (! restraint.fixed_atom_flags[1]) { 
		     idx = 3*(restraint.atom_index_2);
		     *gsl_vector_ptr(df, idx  ) += xP2_contrib;
		     *gsl_vector_ptr(df, idx+1) += yP2_contrib;
		     *gsl_vector_ptr(df, idx+2) += zP2_contrib;
		  }

		  if (! restraint.fixed_atom_flags[2]) { 
		     idx = 3*(restraint.atom_index_3);
		     *gsl_vector_ptr(df, idx  ) += xP3_contrib;
		     *gsl_vector_ptr(df, idx+1) += yP3_contrib;
		     *gsl_vector_ptr(df, idx+2) += zP3_contrib;
		  }

		  if (! restraint.fixed_atom_flags[3]) { 
		     idx = 3*(restraint.atom_index_4);
		     *gsl_vector_ptr(df, idx  ) += xP4_contrib;
		     *gsl_vector_ptr(df, idx+1) += yP4_contrib;
		     *gsl_vector_ptr(df, idx+2) += zP4_contrib;
		  }
	       }
	    }
	    catch (const std::runtime_error &rte) {
	       std::cout << "Caught runtime_error" << rte.what() << std::endl;
	    } 
	 } 
      }
   }
}
