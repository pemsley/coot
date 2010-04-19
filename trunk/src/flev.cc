/* src/flev.cc
 * 
 * Copyright 2010 The University of Oxford
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

/*  ----------------------------------------------------------------------- */
/*               Flattened Ligand Environment View  Interface               */
/*  ----------------------------------------------------------------------- */

#include "c-interface-ligands.hh"
#include "mmdb-extras.h"
#include "mmdb.h"

#include "graphics-info.h"
#include "c-interface.h"


void fle_view_internal(int imol, const char *chain_id, int res_no, const char *ins_code, 
		       int imol_ligand_fragment, 
		       const char *prodrg_output_flat_mol_file_name,
		       const char *prodrg_output_flat_pdb_file_name) {

   float residues_near_radius = 4.5;
   
   // Plan:
   //
   // read in the prodrg output flat pdb file.
   // 
   // lsq match from the residue at imol/chain_id/resno/ins_code to that the first
   // residue in the prodrg flat molecule.  Get the RTop.
   //
   // Get the residues in imol that are within 4A (say) of the residue
   // at imol/chain_id/resno/ins_code.  Find the centres of these
   // residue (construct a helper class)
   // 
   // class fle_residues_helper_t { 
   //    centre
   //    residue spec
   //    residue name
   // }
   //
   //
   // Rotate centres by RTop
   //
   // Write a file of lines like:
   //
   // RES flat_pos_x flat_pos_y TYPE spec
   // 
   // (flat_pos_x, flat_pos_y are in the coordinate system of the
   // flattened ligand)
   //
   // Additionally, a RES card can have the immediately following record:
   //
   // BOND to_atom_x to_atom_y colour_name
   //

   atom_selection_container_t flat = get_atom_selection(prodrg_output_flat_pdb_file_name);
   if (flat.read_success) {
      if (is_valid_model_molecule(imol_ligand_fragment)) { 
	 CMMDBManager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
	 CResidue  *res_ref = coot::util::get_residue(chain_id, res_no, ins_code, mol);
	 CResidue *flat_res = coot::util::get_residue("", 1, "", flat.mol); // prodrg attribs
	 if (res_ref) { 

	    CMMDBManager *ligand_mol =
	       graphics_info_t::molecules[imol_ligand_fragment].atom_sel.mol;
	    int every_nth = 1;
	    std::vector<coot::lsq_range_match_info_t> matches;
	    coot::lsq_range_match_info_t match(1, 1, "", res_no, res_no, chain_id,
					       COOT_LSQ_ALL);
	    matches.push_back(match);
	    std::pair<short int, clipper::RTop_orth> lsq_mat = 
	       coot::util::get_lsq_matrix(flat.mol, ligand_mol, matches, every_nth);

	    if (lsq_mat.first) { 
	       std::vector<CResidue *> residues =
		  coot::residues_near_residue(res_ref, mol, residues_near_radius);

	       // residues needs to be filtered to remove waters that
	       // are not connected to a protein atom.

	       std::vector<CResidue *> filtered_residues =
		  coot::filter_residues_by_solvent_contact(res_ref, mol, residues, 3.6);

	       // for the atoms in the ligand only, for the moment.
	       // More would require a class containing with and
	       // without solvent accessibilites for each residue.
	       // Maybe that can be another function.  And that would
	       // need to consider neighbours of neighbours, perhaps
	       // done with a larger radius.
	       //
	       coot::dots_representation_info_t dots;
	       std::vector<std::pair<coot::atom_spec_t, float> > s_a_v = 
		  dots.solvent_accessibilities(res_ref, filtered_residues);

	       // for the ligand environment residues:
	       std::vector<coot::solvent_exposure_difference_helper_t> sed = 
		  dots.solvent_exposure_differences(res_ref, filtered_residues);

	       write_solvent_accessibilities(s_a_v, flat_res);

	       if (0)
		  for (unsigned int i=0; i<s_a_v.size(); i++)
		     std::cout << "   " << i << " " << s_a_v[i].first << " " << s_a_v[i].second
			       << std::endl;

	       //
	       std::map<std::string, std::string> name_map =
		  coot::make_flat_ligand_name_map(flat_res);

	       // Now what are the bonds between the ligand and residues?
	       // 
	       // a vector of fle-ligand-bond which contain
	       // ligand-atom-name residue-spec and bond type
	       // (acceptor/donor).
	       //
	       std::vector<coot::fle_ligand_bond_t> bonds_to_ligand =
		  coot::get_fle_ligand_bonds(res_ref, filtered_residues, name_map);

	       if (0) 
		  std::cout << "Found ================== " << bonds_to_ligand.size()
			    << " ==================== bonds to ligand " << std::endl;

	       if (filtered_residues.size()) {
		  std::vector<coot::fle_residues_helper_t> centres(filtered_residues.size());
		  for (unsigned int ires=0; ires<filtered_residues.size(); ires++) { 
		     CResidue *res_copy = coot::util::deep_copy_this_residue(filtered_residues[ires]);
		     std::string res_name = filtered_residues[ires]->GetResName();
		     coot::util::transform_atoms(res_copy, lsq_mat.second);
		     std::pair<bool, clipper::Coord_orth> c =
			coot::util::get_residue_centre(res_copy);
		     if (c.first) {
			if (0) 
			   std::cout << "DEBUG:: creating fle_centre with centre "
				     << c.second.format() << std::endl;
			coot::fle_residues_helper_t fle_centre(c.second,
							       coot::residue_spec_t(filtered_residues[ires]),
							       res_name);
			centres[ires] = fle_centre;
		     } else {
			std::cout << "WARNING:: failed to get residue centre for "
				  << coot::residue_spec_t(res_copy) << std::endl;
		     }
		     delete res_copy;
		  }


		  write_fle_centres(centres, bonds_to_ligand, sed, flat_res);
		  
	       }
	    }
	 }
      }
   } 
}

std::map<std::string, std::string>
coot::make_flat_ligand_name_map(CResidue *flat_res) {

   double bond_to_H_dist = 1.1;

   double b2Hd2 = bond_to_H_dist * bond_to_H_dist;
   std::map<std::string, std::string> map;
   PPCAtom residue_atoms = 0;
   int n_residue_atoms;
   flat_res->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      CAtom *at_i = residue_atoms[iat];
      std::string ele_i = at_i->element;
      clipper::Coord_orth pt_i(at_i->x, at_i->y, at_i->z);
      if (ele_i == " H") { 
	 for (int jat=0; jat<n_residue_atoms; jat++) {
	    if (iat != jat) { 
	       CAtom *at_j = residue_atoms[jat];
	       std::string ele_j = at_j->element;
	       if (ele_j != " H") { 
		  clipper::Coord_orth pt_j(at_j->x, at_j->y, at_j->z);
		  if ((pt_i - pt_j).lengthsq() < b2Hd2) {
		     map[at_j->name] = at_i->name;
		     break;
		  } 
	       }
	    }
	 }
      }
   }

//    std::map<std::string, std::string>::const_iterator it;
//    std::cout << "=== name map: === " << std::endl;
//    for (it=map.begin(); it!=map.end(); it++) {
//       std::cout << "  :" << it->first << ": :" << it->second << ":" << std::endl;
//    }
//    std::cout << "=== === " << std::endl;

   return map;
} 


void
coot::write_fle_centres(const std::vector<fle_residues_helper_t> &v,
			const std::vector<coot::fle_ligand_bond_t> &bonds_to_ligand,
			const std::vector<coot::solvent_exposure_difference_helper_t> &sed, 
			CResidue *res_flat) {

   std::string file_name = "coot-tmp-fle-view-residue-info.txt";
   std::ofstream of(file_name.c_str());
   if (!of) {
      std::cout << "failed to open output file " << file_name << std::endl;
   } else {

      for (unsigned int i=0; i<v.size(); i++) {
	 if (v[i].is_set) { 
	    of << "RES " << v[i].centre.x() << " " << v[i].centre.y() << " "
	       << v[i].centre.z() << " "
	       << v[i].residue_name << " " << v[i].spec.chain << v[i].spec.resno << "\n";

	    // now, was there a bond with this res spec?
	    // 
	    // if so find the x y coordinates of the atom with this
	    // name in the res_flat and write them out.
	    // 
	    for (unsigned int ib=0; ib<bonds_to_ligand.size(); ib++) { 
	       if (bonds_to_ligand[ib].res_spec == v[i].spec) {
		  PPCAtom residue_atoms = 0;
		  int n_residue_atoms;
		  res_flat->GetAtomTable(residue_atoms, n_residue_atoms);
		  for (unsigned int iat=0; iat<n_residue_atoms; iat++) {
		     CAtom *at = residue_atoms[iat];
		     std::string atom_name(at->name);
		     if (atom_name == bonds_to_ligand[ib].ligand_atom_name) {
			// the coordinates of the matching flat ligand
			// atom (done like solvent accessibilites).
			of << "BOND " << at->name
			   << " length: " << bonds_to_ligand[ib].bond_length
			   << " type: " << bonds_to_ligand[ib].bond_type
			   << "\n";
			break;
		     }
		  }
	       }
	    }


	    // Is there a solvent accessibility associated with this residue?
	    // 
	    for (unsigned int ised=0; ised<sed.size(); ised++) { 
	       if (sed[ised].res_spec == v[i].spec) {
		  of << "SOLVENT_ACCESSIBILITY_DIFFERENCE "
		     << sed[ised].exposure_fraction_holo << " " 
		     << sed[ised].exposure_fraction_apo << "\n";
	       }
	    }
	    
	 } else {
	    std::cout << "WARNING centre " << i << " was not set " << std::endl;
	 } 
      }
   }

}

// the reference_residue is the flat residue
void
coot::write_solvent_accessibilities(const std::vector<std::pair<coot::atom_spec_t, float> > &sav,
				    CResidue *reference_residue) {

   // for each of the atoms in reference_residue, find the spec in sav
   // and write out the coordinates and the atom name and the
   // solvent_accessibility.

   std::string file_name = "coot-tmp-fle-view-solvent-accessibilites.txt";
   std::ofstream of(file_name.c_str());
   if (!of) {
      std::cout << "failed to open output file " << file_name << std::endl;
   } else {
      PPCAtom residue_atoms = 0;
      int n_residue_atoms;
      reference_residue->GetAtomTable(residue_atoms, n_residue_atoms);
      for (unsigned int i=0; i<n_residue_atoms; i++) {
	 CAtom *at = residue_atoms[i];
	 std::string atom_name(at->GetAtomName());
	 for (unsigned int j=0; j<sav.size(); j++) {
	    if (atom_name == sav[j].first.atom_name) {
	       of << "ATOM:" << atom_name << " " << at->x << " " <<  at->y << " " << at->z
		  << " " << sav[j].second << "\n";
	       break;
	    }
	 }
      }
   }
}

// if there is a name map, use it, otherwise just bond to the found
// atom.  Using the name map allows bond to hydrogens hanging off of
// ligand atoms (e.g. the H on an O or an H on an N). The Hs are not
// in res_ref (often/typically), but are in the flat residue.
// 
std::vector<coot::fle_ligand_bond_t>
coot::get_fle_ligand_bonds(CResidue *res_ref,
			   const std::vector<CResidue *> &residues,
			   const std::map<std::string, std::string> &name_map) {
   std::vector<coot::fle_ligand_bond_t> v;

   // do it by hand, rather than create a molecule and use SeekContacts();
   double bond_length = 3.3;

   double bl_2 = bond_length * bond_length;
   PPCAtom ref_residue_atoms = 0;
   int n_ref_residue_atoms;
   res_ref->GetAtomTable(ref_residue_atoms, n_ref_residue_atoms);
   for (unsigned int iat=0; iat<n_ref_residue_atoms; iat++) {
      std::string ele_ref = ref_residue_atoms[iat]->element;
      if (ele_ref != " C") { 
	 clipper::Coord_orth ref_pt(ref_residue_atoms[iat]->x,
				    ref_residue_atoms[iat]->y,
				    ref_residue_atoms[iat]->z);
	 for (unsigned int ires=0; ires<residues.size(); ires++) {
	    PPCAtom residue_atoms = 0;
	    int n_residue_atoms;
	    residues[ires]->GetAtomTable(residue_atoms, n_residue_atoms);
	    for (unsigned int jat=0; jat<n_residue_atoms; jat++) {
	       std::string ele = residue_atoms[jat]->element;
	       if (ele != " C") { 
		  clipper::Coord_orth pt(residue_atoms[jat]->x,
					 residue_atoms[jat]->y,
					 residue_atoms[jat]->z);
		  double diff_2 = (ref_pt - pt).lengthsq();
		  if (diff_2 < bl_2) {
		     std::string atom_name = ref_residue_atoms[iat]->name;
		     double bl = sqrt(diff_2);
		     coot::residue_spec_t res_spec(residues[ires]);

		     // donor/acceptor from the residue to the ligand
		     // 
		     int bond_type = fle_ligand_bond_t::H_BOND_DONOR_SIDECHAIN;

		     // it's not a donor from a mainchain oxygen
		     //
		     std::string residue_atom_name(residue_atoms[jat]->name);
		     if (ele == " O")
			bond_type = fle_ligand_bond_t::H_BOND_ACCEPTOR_SIDECHAIN;
		     if (ele == " N")
			bond_type = fle_ligand_bond_t::H_BOND_DONOR_SIDECHAIN;
		     if (residue_atom_name == " O  ")
			bond_type = fle_ligand_bond_t::H_BOND_ACCEPTOR_MAINCHAIN;
		     if (residue_atom_name == " N  ")
			bond_type = fle_ligand_bond_t::H_BOND_DONOR_MAINCHAIN;
		     if (residue_atom_name == " H  ")
			bond_type = fle_ligand_bond_t::H_BOND_DONOR_MAINCHAIN;

		     // note, can't use [] because name_map is const
		     std::map<std::string, std::string>::const_iterator it =
			name_map.find(atom_name);
		     if (it != name_map.end()) { 
			atom_name = it->second;
			bond_type = fle_ligand_bond_t::H_BOND_ACCEPTOR_SIDECHAIN;
			if (coot::is_main_chain_p(residue_atoms[jat]))
			   bond_type = fle_ligand_bond_t::H_BOND_ACCEPTOR_MAINCHAIN;
		     }

		     coot::fle_ligand_bond_t bond(atom_name, bond_type, bl, res_spec);
		     v.push_back(bond);
		     break;
		  }
	       }
	    }
	 }
      }
   }
   return v;
}
