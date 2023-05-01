/* ligand/test-ligand.cc
 * 
 * Copyright 2007, 2008 by The University of Oxford
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms ofn the GNU General Public License as published by
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
 * 02110-1301, USA.
*/

#include "clipper/core/rotation.h"
#include "coords/mmdb-extras.h"
#include "coords/mmdb.h"
#include "coot-utils/coot-coord-utils.hh"
#include "torsion-general.hh"
#include <dirent.h>

#include "rotamer.hh"



#include <algorithm> // for reverse()

class compare_stats_t {
 public:
   float angle;
   float sum_deviance;
   int n_atoms;
   compare_stats_t() {
      angle = -999.9;
      sum_deviance = 999.9;
      n_atoms = 0;
   }
   compare_stats_t(float ang_in, float dev_in, int n_atoms_in) {
      angle = ang_in;
      sum_deviance = dev_in;
      n_atoms = n_atoms_in;
   }
};

// return 1 on rtop being identity matrix
bool identity_matrix_test(const clipper::RTop_orth &rtop) {

   clipper::Coord_orth t(rtop.trn());
   double lensq = t.lengthsq();
   //std::cout << "lensq " << lensq << std::endl;
   if (lensq > 0.01) return 0;
   clipper::Rotation rot(rtop.rot());
   double abs_angle = 180.0/M_PI * rot.abs_angle();
   std::cout << "abs_angle " << abs_angle << std::endl;
   if (abs_angle > 0.01) return 0;
   return 1;
}

// return the matrix rotation angle
compare_stats_t
compare_by_overlap(std::string chain_id, int resno, mmdb::Manager *mol1, mmdb::Manager *mol2) {

   compare_stats_t stats;
   std::vector<coot::lsq_range_match_info_t> matches;
   matches.push_back(coot::lsq_range_match_info_t(resno, resno, chain_id,
						  resno, resno, chain_id,
						  COOT_LSQ_ALL));

   std::pair<short int, clipper::RTop_orth> p = coot::util::get_lsq_matrix(mol1, mol2, matches, 1);

   if (!p.first) {
      std::cout << "   Couldn't get matrix" << std::endl;
      return stats;
   } else {
      std::cout << "Got operator: \n" << p.second.format() << std::endl;
      if (identity_matrix_test(p.second)) {
	 std::cout << "   Identity matrix detected.  Overlap not tested." << std::endl
		   << "   Failed" << std::endl;
	 return stats;
      }
	 
      coot::util::transform_mol(mol2, p.second);
      float atom_dist_diff = 0.0;
      int n_atoms_matched = 0;

      int imod = 1;
      
      mmdb::Model *model_p_1 = mol1->GetModel(imod);
      mmdb::Chain *chain_p_1;
      int nchains_1 = model_p_1->GetNumberOfChains();
      for (int ichain_1=0; ichain_1<nchains_1; ichain_1++) {
	 chain_p_1 = model_p_1->GetChain(ichain_1);
	 int nres_1 = chain_p_1->GetNumberOfResidues();
	 mmdb::PResidue residue_p_1;
	 mmdb::Atom *at_1;
	 for (int ires_1=0; ires_1<nres_1; ires_1++) { 
	    residue_p_1 = chain_p_1->GetResidue(ires_1);
	    int n_atoms_1 = residue_p_1->GetNumberOfAtoms();
	 
	    for (int iat_1=0; iat_1<n_atoms_1; iat_1++) {
	       at_1 = residue_p_1->GetAtom(iat_1);
	       coot::atom_spec_t s(at_1);

	       // other molecule:
	       mmdb::Model *model_p_2 = mol2->GetModel(imod);
	       mmdb::Chain *chain_p_2;
	       int nchains_2 = model_p_2->GetNumberOfChains();
	       for (int ichain_2=0; ichain_2<nchains_2; ichain_2++) {
		  chain_p_2 = model_p_2->GetChain(ichain_2);
		  if (std::string(chain_p_1->GetChainID()) == std::string(chain_p_2->GetChainID())) {
		     int nres_2 = chain_p_2->GetNumberOfResidues();
		     mmdb::PResidue residue_p_2;
		     mmdb::Atom *at_2;
		     for (int ires_2=0; ires_2<nres_2; ires_2++) { 
			residue_p_2 = chain_p_2->GetResidue(ires_2);
			if (residue_p_1->GetSeqNum() == residue_p_2->GetSeqNum()) { 
			   int n_atoms_2 = residue_p_2->GetNumberOfAtoms();
			   
			   for (int iat_2=0; iat_2<n_atoms_2; iat_2++) {
			      at_2 = residue_p_2->GetAtom(iat_2);
			      if (s.matches_spec(at_2)) {
				 n_atoms_matched++;
				 ::Cartesian p1(at_1->x, at_1->y, at_1->z);
				 ::Cartesian p2(at_2->x, at_2->y, at_2->z);
				 float d = LineLength(p1, p2);
				 atom_dist_diff += d;
			      }
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
      std::cout << "  Number of atoms matched: " << n_atoms_matched << " deviation sum: "
		<< atom_dist_diff << std::endl;
      clipper::Rotation rot(p.second.rot());
      float abs_angle = 180.0/M_PI * rot.abs_angle();
      stats.angle = abs_angle;
      stats.n_atoms = n_atoms_matched;
      stats.sum_deviance = atom_dist_diff;
   }
   return stats;
}

int test_torsion_general(atom_selection_container_t asc, std::string pdb_filename) { 

   if (!asc.read_success) {
	 std::cout << "   Problem reading " << pdb_filename << std::endl;
	 return 2; 
   }
   bool tested = 0;
   std::string target_chain_id = "B"; // but note that atom specs
   // are always blank, because
   // chain id gets lost on
   // create_mmdbmanager_from_residue()
   int target_resno = 81;
   bool reverse_torsion = 0;
   mmdb::Model *model_p = asc.mol->GetModel(1);
   if (! model_p) {
      std::cout << "   Model not found " << std::endl;
   } else { 
      mmdb::Chain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 if (target_chain_id == chain_p->GetChainID()) { 
	    int nres = chain_p->GetNumberOfResidues();
	    mmdb::PResidue residue_p;
	    mmdb::Atom *at;
	    for (int ires=0; ires<nres; ires++) { 
	       residue_p = chain_p->GetResidue(ires);
	       if (residue_p->GetSeqNum() == target_resno) { 
		  std::vector<std::string> atom_names;
		  atom_names.push_back(" C  ");
		  atom_names.push_back(" CA ");
		  atom_names.push_back(" CB ");
		  atom_names.push_back(" CG ");
		  std::vector<coot::atom_spec_t> torsion_general_atom_specs;
		  std::vector<coot::atom_spec_t> reverse_atom_specs;
		  std::vector<std::string>::iterator it;
		  for (it = atom_names.begin(); it != atom_names.end(); it++) {
		     coot::atom_spec_t as("", target_resno, "", *it, "");
		     torsion_general_atom_specs.push_back(as);
		     reverse_atom_specs.push_back(as);
		  }
		  std::reverse(reverse_atom_specs.begin(), reverse_atom_specs.end());
		  mmdb::Manager *res_mol_1 =
		     coot::util::create_mmdbmanager_from_residue(residue_p);
		  mmdb::Manager *res_mol_2 =
		     coot::util::create_mmdbmanager_from_residue(residue_p);

		  mmdb::Residue *res_copy_1 = coot::util::get_residue("", target_resno, "", res_mol_1);
		  mmdb::Residue *res_copy_2 = coot::util::get_residue("", target_resno, "", res_mol_2);

		  if (! res_copy_1) {
		     std::cout << "   Error can't find residue in new molecule" << std::endl;
		     return 2;
		  }
		  if (! res_copy_2) {
		     std::cout << "   Error can't find residue in new molecule" << std::endl;
		     return 2;
		  }
		     
		  coot::torsion_general tg_1(res_copy_1, res_mol_1, torsion_general_atom_specs);
		  coot::torsion_general tg_2(res_copy_2, res_mol_2,         reverse_atom_specs);
		  Tree tree1 = tg_1.GetTree();
		  Tree tree2 = tg_2.GetTree();
		  double diff = 20;
		  int istat_1 = tg_1.change_by(diff, &tree1); // fiddle with tree
		  int istat_2 = tg_2.change_by(diff, &tree2);

		  mmdb::PPAtom residue_atoms;
		  int n_residue_atoms;
		  res_copy_1->GetAtomTable(residue_atoms, n_residue_atoms);
		  std::vector<::Cartesian > coords;

		  res_mol_1->WritePDBASCII("rotated-1.pdb");
		  res_mol_2->WritePDBASCII("rotated-2.pdb");

		  if (istat_1 != 0) {
		     std::cout << "  Error rotating coordinates 1" << std::endl;
		     return 2;
		  } else {
		     if (istat_2 != 0) {
			std::cout << "  Error rotating coordinates 2" << std::endl;
			return 2;
		     } else {
			tested=1;
			// So in 1 we have turned the top fragment by 20 degrees
			// in 2 we have turned the bottom by 20 degrees.
			//
			// after overlapping, the structures should be very similar
			compare_stats_t stats =
			   compare_by_overlap("", target_resno, res_mol_1, res_mol_2);
			if (fabs(diff - stats.angle) > 0.001) {
			   std::cout << "  Error in angle of rotation matrix" << std::endl;
			   return 2;
			} else {
			   if (stats.n_atoms < 6) {
			      std::cout << "  Error:: not enough atoms matched" << std::endl;
			      return 2;
			   } else {
			      if (stats.sum_deviance > 0.01) {
				 std::cout << "  Error:: atoms too different" << std::endl;
				 return 2;
			      } else { 
				 std::cout << "  Similar rotation, low deviance: pass"
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
      if (!tested) {
	 std::cout << "   Failed to find residue \"" << target_chain_id << "\" "
		   << target_resno << " in input coordinates" << std::endl;
	 return 2;
      }
   }
   return 0;
}

namespace coot { 
   bool matches_data_name(const std::string &file_str) {
      bool m = 1; 
      return m; 
   }
}

void get_rotamer_probabilities(const std::string &rotamer_probability_dir) {

   std::vector<std::string> prob_data_list; 
   DIR *ref_struct_dir = opendir(rotamer_probability_dir.c_str()); 

   if (ref_struct_dir == NULL) { 
      std::cout << "An error occured on opendir" << std::endl;
   } else { 

      std::cout << rotamer_probability_dir
		<< " successfully opened" << std::endl; 

      // loop until the end of the filelist (readdir returns NULL)
      // 
      struct dirent *dir_ent; 

      while(1) { 
	 dir_ent = readdir(ref_struct_dir); 

	 if (dir_ent != NULL) { 

	    std::string file_str(dir_ent->d_name); 
	    if (coot::matches_data_name(file_str)) {

	       // Construct a file in a unix directory - non-portable?
	       // 
	       std::string s(rotamer_probability_dir);
	       s += "/"; 
	       s += file_str; 

	       prob_data_list.push_back(s); 
	    }
	 } else { 
	    break;
	 }
      }
      closedir(ref_struct_dir); 
   } 
   // return prob_data_list; 
} 

void rotamer_tables() {

   std::string dir = "../../coot/rama-data";
   coot::rotamer_probability_tables tables;
   tables.fill_tables(dir);

}

int test_torsion(const std::string &pdb_filename) {

   atom_selection_container_t asc = get_atom_selection(pdb_filename, true, false, false);
   int retval = test_torsion_general(asc, pdb_filename);
   return retval;
}

#include "ligand.hh"

int test_spherical(const std::string &pdb_filename,
		   const std::string &mtz_filename) {

   int r = 0;

   atom_selection_container_t asc = get_atom_selection(pdb_filename, true, false, false);
   clipper::Xmap<float> xmap;
   std::string f_col   = "FWT";
   std::string phi_col = "PHWT";
   bool use_weights = 0;
   bool is_diff_map = 0;
   coot::ligand lig;
   bool stat = false;
   try {
      stat = lig.map_fill_from_mtz(mtz_filename, f_col, phi_col, "",
				   use_weights, is_diff_map);
   }
   catch (const clipper::Message_fatal &m) {
      std::cout << "ERROR:: reading " << mtz_filename << std::endl;
   }
   if (! stat) {
      std::cout << "ERROR: in filling map from mtz file: " << mtz_filename
		<< std::endl;
   } else {
      int imod = 1;
      mmdb::Model *model_p = asc.mol->GetModel(imod);
      if (model_p) {
	 std::pair<float, float> mdoa = lig.mean_and_variance_where_the_atoms_are(asc.mol);
	 std::cout << "here with mdoa " << mdoa.first << " " << mdoa.second << std::endl;
	 if (mdoa.first > 0.0) {
	    int n_chains = model_p->GetNumberOfChains();
	    for (int ichain=0; ichain<n_chains; ichain++) {
	       mmdb::Chain *chain_p = model_p->GetChain(ichain);
	       int nres = chain_p->GetNumberOfResidues();
	       for (int ires=0; ires<nres; ires++) {
		  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
		  std::string rn(residue_p->GetResName());
		  if (rn == "HOH") {
		     int n_atoms = residue_p->GetNumberOfAtoms();
		     if (n_atoms == 1) {
			for (int iat=0; iat<n_atoms; iat++) {
			   mmdb::Atom *at = residue_p->GetAtom(iat);
			   if (! at->isTer()) {
			      clipper::Coord_orth pt = coot::co(at);
			      coot::ligand::spherical_density_score_t non_spherical_score =
				 lig.spherical_density_score(pt, mdoa.first);
			      std::cout << "water " << pdb_filename <<  " " << coot::atom_spec_t(at)
					<< " density " << non_spherical_score.density_at_position
					<< " non-spherical-score " << non_spherical_score.non_spherical_score
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
   return r;
}


int main(int argc, char **argv) {

   int r = 1;

   if (false) {
      if (argc < 2) {
	 std::cout << "Usage: Need to pass a pdb file name" << std::endl;
      } else {
	 std::string pdb_filename = argv[1];
	 r = test_torsion(pdb_filename);
      }
   }

   if (false)
      rotamer_tables();

   if (true) {
      if (argc < 3) {
	 std::cout << "Usage: Need to pass a pdb file name and map" << std::endl;
      } else {
	 std::string pdb_filename = argv[1];
	 std::string mtz_filename = argv[2];
	 r = test_spherical(pdb_filename, mtz_filename);
      }
   }

   return r;
}
