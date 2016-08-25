
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

   std::vector<bool> fixed_flag(4);
   fixed_flag[0] = 0;  // not fixed
   fixed_flag[1] = 0;
   fixed_flag[2] = 0;
   fixed_flag[3] = 0;

   atom_1_sel = first_sel;
   atom_2_sel = first_sel;
   atom_3_sel = second_sel;
   atom_4_sel = second_sel;

   n_atom_1 = n_first_res_atoms;
   n_atom_2 = n_first_res_atoms;
   n_atom_3 = n_second_res_atoms;
   n_atom_4 = n_second_res_atoms;
   fixed_flag[0] = is_fixed_first;
   fixed_flag[1] = is_fixed_first;
   fixed_flag[2] = is_fixed_second;
   fixed_flag[3] = is_fixed_second;

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

				 double target_omega = 180.0;
				 double esd = 2.0; // 5.0 lets slip 72A in 2bmd to trans
				 add(TRANS_PEPTIDE_RESTRAINT, index1, index2, index3, index4,
				     fixed_flag,
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

   // b*b * [ a.(bxc)/b ]
   double E = clipper::Coord_orth::dot(a,clipper::Coord_orth::cross(b,c)) *
      sqrt( b.lengthsq() );

   // b*b * [ -a.c+(a.b)(b.c)/(b*b) ] = -a.c*b*b + (a.b)(b.c)
   double G = - clipper::Coord_orth::dot(a,c)*b.lengthsq()
      + clipper::Coord_orth::dot(a,b)*clipper::Coord_orth::dot(b,c);

   double theta = clipper::Util::rad2d(atan2(E,G));
   if (theta < 0)
      theta += 360;

   if (false) { 
      if (clipper::Util::isnan(theta)) {
	 std::string mess = "WARNING: distortion_score_torsion() observed torsion theta is a NAN!";
	 throw std::runtime_error(mess);
      }
   }

   double diff = theta - restraint.target_value;

   if (false) { // debug 
      double pen = diff*diff/(restraint.sigma * restraint.sigma);
	 std::cout << "distortion_trans_peptide theta (calc): " << theta 
		   << " periodicity " << restraint.periodicity
		   << " target "      << restraint.target_value
		   << " diff: " << diff << endl ;
      std::cout << "in distortion_trans_peptide: sigma = " << restraint.sigma
		<< ", weight=" << pow(restraint.sigma,-2.0)
		<< " and diff is " << diff << std::endl;
      std::cout << "distortion score trans-peptide "
		<< diff*diff/(restraint.sigma * restraint.sigma) << " "
		<< std::endl;
   }
      
   return diff*diff/(restraint.sigma * restraint.sigma);
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

	 const simple_restraint &restraint = (*restraints)[i];
      
	 if ( restraint.restraint_type == TRANS_PEPTIDE_RESTRAINT) {

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
	       distortion_torsion_gradients_t dtg =
		  fill_distortion_torsion_gradients(P1, P2, P3, P4);

	       if (true) {

		  if (dtg.theta < 0)
		     dtg.theta += 360;
		  double diff = dtg.theta - restraint.target_value;
		  
		  if (false) 
		     std::cout << "in df_trans_peptide: dtg.theta is " << dtg.theta 
			       <<  " and target is " << restraint.target_value 
			       << " and diff is " << diff 
			       << " and periodicity: " << restraint.periodicity <<  endl;

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
