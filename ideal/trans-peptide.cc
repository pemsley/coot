
#include "simple-restraint.hh"
int
coot::restraints_container_t::add_link_trans_peptide(std::string link_type,
						     mmdb::Residue *first,
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

   for (int i=0; i<geom.link_size(); i++) {
      if (geom.link(i).link_id == link_type) {

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
				       double esd = 5.0;
				       add(TORSION_RESTRAINT, index1, index2, index3, index4,
					   fixed_flag,
					   target_omega,
					   esd,
					   1.2, // dummy value
					   1);
				       n_torsion++;
				       std::cout << "!!!!!!!!!!!!!!!! added peptide trans restraint"
						 << std::endl;
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
   return n_torsion; 
}
