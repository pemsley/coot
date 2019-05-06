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

#include <iomanip>

#include "simple-restraint.hh"
int
coot::restraints_container_t::add_link_trans_peptide(mmdb::Residue *first,
						     mmdb::Residue *second,
						     short int is_fixed_first,
						     short int is_fixed_second,
						     const coot::protein_geometry &geom) {
   int n_trans_peptide_torsion = 0;

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

				 if (none_are_fixed_p(fixed_flags)) {

				    if (false) {
				       std::cout << "debug:: making trans-peptide restraint with fixed flags: ";
				       for (std::size_t jj=0; jj<fixed_flags.size(); jj++)
					  std::cout << " " << fixed_flags[jj];
				       std::cout << std::endl;
				    }
				    double target_omega = 180.0;
				    double esd = 2.0; // 5.0 lets slip 72A in 2bmd to trans
				    // esd = 0.2; // 20171228 trial - does this beat the plane restraints?
				    add(TRANS_PEPTIDE_RESTRAINT, index1, index2, index3, index4,
					fixed_flags,
					target_omega,
					esd,
					1.2, // dummy value
					1);
				    n_trans_peptide_torsion++;
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
   return n_trans_peptide_torsion;
}

//
// Return the distortion score from a single torsion restraint.
// 
// can throw a std::runtime_error if there is a problem calculating the torsion.
// 
double
coot::distortion_score_trans_peptide(const int &restraint_index,
				     const coot::simple_restraint &restraint,
				     const gsl_vector *v) {

   int idx;

   // checked: P1 is CA_1
   //          P2 is C_1
   //          P3 is N_2
   //          P4 is CA_2

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

   // mid-point is a misnomer here
   // I mean "point at the "closest-approach" fraction -
   // which is close to the mid-point, but not quite.

   // i.e. if closest_approach_fraction_CA_CA is 0.9, that means
   // it is a lot closer to CA_2 than CA_1
   //
   double closest_approach_fraction_CA_CA = 0.5;
   double closest_approach_fraction_C_N   = 0.5;
   double best_closest_approach = 0.055; // or whatever it was

   // correct values: 0.473, 0.407, and 0.055
   // 20181108 tried them - they made things worse.
   // Using 0.055 alone was an improvement - I need to revisit
   // these numbers - maybe they are round the "wrong" way.

   double q_CA_CA = 1.0 - closest_approach_fraction_CA_CA;
   double q_C_N   = 1.0 - closest_approach_fraction_C_N;

   clipper::Coord_orth mid_pt_1 = q_CA_CA * P1 + closest_approach_fraction_CA_CA * P4;
   clipper::Coord_orth mid_pt_2 = q_C_N   * P2 + closest_approach_fraction_C_N * P3;

   double dist_sqrd = (mid_pt_2-mid_pt_1).lengthsq();
   // dd is the distance from the "mid-points" to the expected distance
   // between "midpoints" for an ideal trans-peptide

   double delta = sqrt(dist_sqrd) - best_closest_approach;

   double trans_pep_dist_scale_factor = 4000.0; // needs tweaking

   double d = trans_pep_dist_scale_factor * delta * delta;

   return d;
}

// Add in the torsion gradients
//
void coot::my_df_trans_peptides(const gsl_vector *v, 
				void *params, 
				gsl_vector *df) {

   restraints_container_t *restraints = static_cast<restraints_container_t *>(params);

   if (restraints->restraints_usage_flag & coot::TRANS_PEPTIDE_MASK) {

      for (unsigned int i=restraints->restraints_limits_trans_peptide.first;
	   i<=restraints->restraints_limits_trans_peptide.second; i++) {

         const simple_restraint &restraint = restraints->at(i);
	 if (restraint.restraint_type == coot::TRANS_PEPTIDE_RESTRAINT) {

	    int idx;

	    // checked:    P1 is CA_1
	    //             P2 is C_1
	    //             P3 is N_2
	    //             P4 is CA_2

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

	    // mid-point is a misnomer here
	    // I mean "point at the "closest-approach" fraction -
	    // which is close to the real mid-point, but not quite.

	    // i.e. if closest_approach_fraction_CA_CA is 0.9, that means
	    // it is a lot closer to CA_2 than CA_1
	    //
            double closest_approach_fraction_CA_CA = 0.5;
            double closest_approach_fraction_C_N   = 0.5;
            double best_closest_approach = 0.055; // or whatever it was

	    const double &p_CA_CA = closest_approach_fraction_CA_CA; // shorthand aliases
	    const double &p_C_N   = closest_approach_fraction_C_N;

	    double q_CA_CA = 1.0 - closest_approach_fraction_CA_CA;
	    double q_C_N   = 1.0 - closest_approach_fraction_C_N;

	    clipper::Coord_orth mid_pt_1 = q_CA_CA * P1 + closest_approach_fraction_CA_CA * P4;
	    clipper::Coord_orth mid_pt_2 = q_C_N   * P2 + closest_approach_fraction_C_N * P3;

	    double dist_sqrd = (mid_pt_2-mid_pt_1).lengthsq();

	    double trans_pep_dist_scale_factor = 4000.0; // needs tweaking
	    double weight = trans_pep_dist_scale_factor;

	    // d is the distance from the "mid-points" to the expected distance
	    // between "midpoints" for an ideal trans-peptide
	    double b = sqrt(dist_sqrd);
	    double delta = b - best_closest_approach;

	    double dS_ddelta = weight * 2.0 * delta;
	    double db_da = 0.5 / b;

	    // std::cout << "b: " << b << " delta: " << delta << std::endl;

	    double constant_part = dS_ddelta * db_da;

	    double xP1_contrib = constant_part * 2.0 * q_CA_CA * ( mid_pt_1.x() - mid_pt_2.x());
	    double yP1_contrib = constant_part * 2.0 * q_CA_CA * ( mid_pt_1.y() - mid_pt_2.y());
	    double zP1_contrib = constant_part * 2.0 * q_CA_CA * ( mid_pt_1.z() - mid_pt_2.z());

	    double xP2_contrib = constant_part * 2.0 * q_C_N * ( mid_pt_2.x() - mid_pt_1.x());
	    double yP2_contrib = constant_part * 2.0 * q_C_N * ( mid_pt_2.y() - mid_pt_1.y());
	    double zP2_contrib = constant_part * 2.0 * q_C_N * ( mid_pt_2.z() - mid_pt_1.z());

	    double xP3_contrib = constant_part * 2.0 * p_C_N * ( mid_pt_2.x() - mid_pt_1.x());
	    double yP3_contrib = constant_part * 2.0 * p_C_N * ( mid_pt_2.y() - mid_pt_1.y());
	    double zP3_contrib = constant_part * 2.0 * p_C_N * ( mid_pt_2.z() - mid_pt_1.z());

	    double xP4_contrib = constant_part * 2.0 * p_CA_CA * ( mid_pt_1.x() - mid_pt_2.x());
	    double yP4_contrib = constant_part * 2.0 * p_CA_CA * ( mid_pt_1.y() - mid_pt_2.y());
	    double zP4_contrib = constant_part * 2.0 * p_CA_CA * ( mid_pt_1.z() - mid_pt_2.z());

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
   }
}
