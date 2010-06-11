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
#include "lbg-graph.hh"

#include "coot-h-bonds.hh"

std::ostream&
coot::operator<< (std::ostream& s, const pi_stacking_instance_t &stack) {

   s << "[" << stack.type << " " << coot::residue_spec_t(stack.res) << " "
     << stack.overlap_score << " ligand-atom-name :" <<  stack.ligand_cationic_atom_name
     << ": ";
   for (unsigned int i=0; i<stack.ligand_ring_atom_names.size(); i++)
      s << "  :" << stack.ligand_ring_atom_names[i] << ":   " ;
   s << "]";
   return s;
}

void fle_view_internal(int imol, const char *chain_id, int res_no, const char *ins_code, 
		       int imol_ligand_fragment, 
		       const char *prodrg_output_flat_mol_file_name,
		       const char *prodrg_output_flat_pdb_file_name,
		       const char *prodrg_output_3d_pdb_file_name,
		       const char *prodrg_output_dict_cif_file_name) {
   
   float residues_near_radius = 4.5;
   graphics_info_t g;
   coot::protein_geometry *geom_p = g.Geom_p();
   
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

	 std::string res_name_res_ref; // just in case res_ref is NULL
	 if (res_ref)
	    res_name_res_ref = res_ref->GetResName();

	 // We only want to read a cif dictionary if is is not one of
	 // the standard ones (this is because "standard" residues
	 // have corrrect group for linking (typically "L-peptide")
	 // but the file from PRODRG is always "monomer" - and that
	 // means bad news if we subsequently try to look for covalent
	 // bonds.
	 //
	 if (! (coot::standard_residue_name_p(res_name_res_ref)))
	    handle_cif_dictionary(prodrg_output_dict_cif_file_name);
	 
	 
	 if (res_ref) { 
	    
	    std::string ligand_res_name(res_ref->GetResName());
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
	       
	       // for just the ligand atoms.
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
	       // std::vector<coot::fle_ligand_bond_t> bonds_to_ligand =
	       // coot::get_fle_ligand_bonds(res_ref, filtered_residues, name_map);

	       // using new (20100522) h_bond class (based on geometry)
	       //
 	       graphics_info_t g;
 	       std::vector<coot::fle_ligand_bond_t> bonds_to_ligand = 
 		  coot::get_fle_ligand_bonds(res_ref, filtered_residues, mol,
 					     name_map, *geom_p);

	       if (1) 
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
		  
		  std::pair<bool, coot::dictionary_residue_restraints_t> p = 
		     geom_p->get_monomer_restraints(ligand_res_name);
		  
		  if (! p.first) {
		     std::cout << "Failed to get monomer_restraints for PRODRG residue"
			       << std::endl;
		  } else {
		     
		     coot::pi_stacking_container_t pi_stack_info(p.second, filtered_residues, res_ref);
		     
		     write_fle_centres(centres, bonds_to_ligand, sed, pi_stack_info, flat_res);
		     
		  }
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


coot::pi_stacking_container_t::pi_stacking_container_t(const coot::dictionary_residue_restraints_t &monomer_restraints,
						       
						       const std::vector<CResidue *> &residues,
						       CResidue *res_ref) {

   bool debug = 0;

   // --------------------------------------------------------------------
   //          ligand ring systems
   // --------------------------------------------------------------------
   
   // get a list of aromatic bonds, so that they can be used to find
   // aromatic rings.
   std::vector<std::vector<std::string> > ring_list = get_ligand_aromatic_ring_list(monomer_restraints);

   // float pi_overlap_thresh = 0.2;
   float pi_overlap_thresh = 0.0015; // play value
   // float pi_overlap_thresh = 0.2; 
   for (unsigned int iring=0; iring<ring_list.size(); iring++) {
      try {
	 std::pair<clipper::Coord_orth, clipper::Coord_orth> ligand_ring_pi_pts = 
	    get_ring_pi_centre_points(ring_list[iring], res_ref);

	 if (debug) {
	    std::cout << " ligand ring ";
	    for (unsigned int iat=0; iat<ring_list[iring].size(); iat++)
	       std::cout << ring_list[iring][iat] << "  ";
	    
	    std::cout << " points " << ligand_ring_pi_pts.first.format() << " "
		      << ligand_ring_pi_pts.second.format() << std::endl;
	 }
	 
	 for (unsigned int ires=0; ires<residues.size(); ires++) {

	    if (debug) {
	       std::string res_name(residues[ires]->GetResName());
	       std::cout << "==== Environment residue " << coot::residue_spec_t(residues[ires])
			 << " " << res_name << std::endl;
	    }

	    // return a pair that is the score and the stacking type
	    std::pair<float, int> pi_overlap_1 =
	       get_pi_overlap_to_ligand_ring(residues[ires], ligand_ring_pi_pts.first);
	    std::pair<float, int> pi_overlap_2 =
	       get_pi_overlap_to_ligand_ring(residues[ires], ligand_ring_pi_pts.second);

	    if (debug) 
	       std::cout << "   Overlaps:  score "
			 << pi_overlap_1.first << " type: " << pi_overlap_1.second << "  score: "
			 << pi_overlap_2.first << " type: " << pi_overlap_2.second << std::endl;
	       
	    if (pi_overlap_1.first > pi_overlap_thresh) {
	       coot::pi_stacking_instance_t st(residues[ires],
					       pi_overlap_1.second,
					       ring_list[iring]);
	       st.overlap_score = pi_overlap_1.first;
	       stackings.push_back(st);
	    }
	    if (pi_overlap_2.first > pi_overlap_thresh) {
	       coot::pi_stacking_instance_t st(residues[ires],
					       pi_overlap_2.second,
					       ring_list[iring]);
	       st.overlap_score = pi_overlap_2.first;
	       stackings.push_back(st);
	    }
	 }
      }
      catch (std::runtime_error rte) {
	 std::cout << "WARNING:: " << rte.what() << std::endl;
      }
   }

   // --------------------------------------------------------------------
   //          ligand cations (-> protein ring systems)
   // --------------------------------------------------------------------
   
   // Now, we need to test the ligand cations against the environment
   // residues that have pi system (not ligand ring systems).
   //
   // OK, so what are the ligand cations?

   std::vector<std::pair<std::string, clipper::Coord_orth> > cation_points =
      get_ligand_cations(res_ref, monomer_restraints);

   for (unsigned int icat=0; icat<cation_points.size(); icat++) {
      for (unsigned int ires=0; ires<residues.size(); ires++) {
	 // get_ligand_cation_residue_pi_overlap
	 float score = get_pi_overlap_to_ligand_cation(residues[ires], cation_points[icat].second);

	 if (score > pi_overlap_thresh) { 
	    // add a stacking to stackings.
	    coot::pi_stacking_instance_t stacking(residues[ires], cation_points[icat].first);
	    stacking.overlap_score = score;
	    stackings.push_back(stacking);
	 }
      }
   }

}

std::vector<std::vector<std::string> >
coot::pi_stacking_container_t::get_ligand_aromatic_ring_list(const coot::dictionary_residue_restraints_t &monomer_restraints) const {

   // get a list of aromatic bonds, so that they can be used to find
   // aromatic rings.
   // 
   std::vector<std::pair<std::string, std::string> > bonds;
   for (unsigned int irest=0; irest<monomer_restraints.bond_restraint.size(); irest++) {
      if (monomer_restraints.bond_restraint[irest].type() == "aromatic") {
	 std::pair<std::string, std::string> p(monomer_restraints.bond_restraint[irest].atom_id_1_4c(),
					       monomer_restraints.bond_restraint[irest].atom_id_2_4c());
	 bonds.push_back(p);
      }
   }
   
   coot::aromatic_graph_t arom(bonds);
   std::vector<std::vector<std::string> > ring_list = arom.ring_list();

   if (0) {
      std::cout << "----------- " << ring_list.size() << " rings ---------- " << std::endl;
      for (unsigned int i=0; i<ring_list.size(); i++) {
	 std::cout << "ring " << i << "\n   ";
	 for (unsigned int j=0; j<ring_list[i].size(); j++) { 
	    std::cout << ring_list[i][j] << "  ";
	 }
	 std::cout << std::endl;
      }
   }
   return ring_list;
}


// by search through res_ref
std::vector<std::pair<std::string, clipper::Coord_orth> >
coot::pi_stacking_container_t::get_ligand_cations(CResidue *res_ref,
						 const coot::dictionary_residue_restraints_t &monomer_restraints) const {

   std::vector<std::pair<std::string, clipper::Coord_orth> > v;
   int n_residue_atoms;
   PPCAtom residue_atoms = NULL;
   res_ref->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) { 
      std::string ele(residue_atoms[iat]->element);
      if (ele == " N") {
	 // how many bonds does this N have?
	 int n_bonds = 0;
	 std::string atom_name(residue_atoms[iat]->name);
	 for (unsigned int ibond=0; ibond<monomer_restraints.bond_restraint.size(); ibond++) { 
	    if ((monomer_restraints.bond_restraint[ibond].atom_id_1_4c() == atom_name) || (monomer_restraints.bond_restraint[ibond].atom_id_2_4c() == atom_name)) {
	       if (monomer_restraints.bond_restraint[ibond].type() == "single")
		  n_bonds++;
	       if (monomer_restraints.bond_restraint[ibond].type() == "double")
		  n_bonds += 2;
	       if (monomer_restraints.bond_restraint[ibond].type() == "triple")
		  n_bonds += 3;
	    } 
	 }

	 if (n_bonds > 3) { // i.e. 4
	    clipper::Coord_orth pt(residue_atoms[iat]->x,
				   residue_atoms[iat]->y,
				   residue_atoms[iat]->z);
	    std::pair<std::string, clipper::Coord_orth> p(atom_name, pt);
	    v.push_back(p);
	 }
      }
   }
   return v;
}


// can throw an exception because it calls get_ring_pi_centre_points()
//
// should return the stacking type, e.g. PI_CATION_STACKING.
// 
std::pair<float, int>
coot::pi_stacking_container_t::get_pi_overlap_to_ligand_ring(CResidue *res,
							     const clipper::Coord_orth &ligand_pi_point) const {

   float pi_pi_score = 0;
   float pi_cation_score = 0;
   
   std::string res_name(res->GetResName());
   int stacking_type = coot::pi_stacking_instance_t::PI_PI_STACKING;

   // First test all the ring systems in the residue
   // 
   std::vector<std::vector<std::string> > atom_names = ring_atom_names(res_name);
   for (unsigned int iring=0; iring<atom_names.size(); iring++) {
      std::pair<clipper::Coord_orth, clipper::Coord_orth> residue_pi_points =
	 get_ring_pi_centre_points(atom_names[iring], res);
      float score_1 = overlap_of_pi_spheres(ligand_pi_point, residue_pi_points.first,
					    0.78, -1, 0.78, -1);
      float score_2 = overlap_of_pi_spheres(ligand_pi_point, residue_pi_points.second,
					    0.78, -1, 0.78, -1);
      if (score_1 > pi_pi_score)
	 pi_pi_score = score_1;
      if (score_2 > pi_pi_score)
	 pi_pi_score = score_2;
   }

   // Now test all the cation-pi interactions
   // 
   std::vector<clipper::Coord_orth> cation_atom_point = get_cation_atom_positions(res);

   // std::cout << "DEBUG:: there are " << cation_atom_point.size()
   // << " cation atom points" << std::endl;
   
   for (unsigned int icat=0; icat<cation_atom_point.size(); icat++) {
      pi_cation_score += overlap_of_cation_pi(ligand_pi_point, cation_atom_point[icat]);
   }

   float score = pi_pi_score;
   if (pi_cation_score > pi_pi_score) {
      score = pi_cation_score;
      stacking_type = coot::pi_stacking_instance_t::PI_CATION_STACKING;
   }
   
   return std::pair<float, int> (score, stacking_type);
}

// Return the best score of the ligand cation overlap to any of the
// ring systems in res (if any of course, typically this function just
// falls through returning 0.0).
// 
float
coot::pi_stacking_container_t::get_pi_overlap_to_ligand_cation(CResidue *res,
							       const clipper::Coord_orth &pt) const {

   float score = 0.0;
   std::string res_name(res->GetResName());
   std::vector<std::vector<std::string> > atom_names = ring_atom_names(res_name);
   for (unsigned int iring=0; iring<atom_names.size(); iring++) {
      std::pair<clipper::Coord_orth, clipper::Coord_orth> residue_pi_points =
	 get_ring_pi_centre_points(atom_names[iring], res);
      float pi_cation_score_1 = overlap_of_cation_pi(pt, residue_pi_points.first);
      float pi_cation_score_2 = overlap_of_cation_pi(pt, residue_pi_points.second );
      if (pi_cation_score_1 > score)
	 score = pi_cation_score_1;
      if (pi_cation_score_2 > score)
	 score = pi_cation_score_2;
   }
   return score;
}


// the sum of the product of the function: m1 exp [ m2 d^2 ]. Where d
// is the is distance of this particular point from either the control
// point (the case of a ring system) or the atom holding the charge (in
// the case of cation in cation-pi interactions).
// 
// m1, m2 pairs from Clarke & Labute:
// 0.78 -1.0
// 7.8  -0.05
// 
float
coot::pi_stacking_container_t::overlap_of_pi_spheres(const clipper::Coord_orth &pt_1,
						     const clipper::Coord_orth &pt_2,
						     const double &m1_pt_1, const double &m2_pt_1,
						     const double &m1_pt_2, const double &m2_pt_2) const {

   // 0.78 exp(-|r-V|^2) where r is a point in space and V is the
   // normally extended ring centre point.  Return the "integral of
   // the function over space".
   
   double score = 0;
   double grid_size = 3.0;
   double grid_step = 0.2;
   double fudge = 5; // so that ideal overlap gives score of 1.0;
   
   for (double delta_x = -grid_size; delta_x<=grid_size; delta_x += grid_step) {
      for (double delta_y = -grid_size; delta_y<=grid_size; delta_y += grid_step) {
	 for (double delta_z = -grid_size; delta_z<=grid_size; delta_z += grid_step) {
	    clipper::Coord_orth t(pt_1.x() + delta_x,
				  pt_1.y() + delta_y,
				  pt_1.z() + delta_z);
	    double d_1_sqd = (pt_1 - t).lengthsq();
	    double d_2_sqd = (pt_2 - t).lengthsq();
	    double score_1 = 0.0;
	    double score_2 = 0.0;
	    if (d_1_sqd < 12.0)
	       score_1 = m1_pt_1 * exp(m2_pt_1 * d_1_sqd);
	    if (d_2_sqd < 12.0)
	       score_2 = m1_pt_2 * exp(m2_pt_2 * d_2_sqd);
	    // std::cout << "scores: " << d_1_sqd << " -> " << score_1 << "     "
	    // << d_2_sqd << " -> " << score_2 << "\n";
	    score += score_1 * score_2;
	 }
      }
   }
   score *= fudge * (grid_step * grid_step * grid_step);
   return score; // double to float
}


// Can throw an exception if not all the ring_atom_names were found in
// res_ref, or lsq plane normal goes bad.
// 
std::pair<clipper::Coord_orth, clipper::Coord_orth>
coot::pi_stacking_container_t::get_ring_pi_centre_points(const std::vector<std::string> &ring_atom_names,
							 CResidue *res_ref) const {
   // dummy points, overwritten.
   clipper::Coord_orth pt_1(0,0,0);
   clipper::Coord_orth pt_2(0,0,0);

   int n_residue_atoms;
   PPCAtom residue_atoms = NULL;
   res_ref->GetAtomTable(residue_atoms, n_residue_atoms);

   // fill this vector
   std::vector<clipper::Coord_orth> aromatic_plane_points;
   for (unsigned int iring_at=0; iring_at<ring_atom_names.size(); iring_at++) {
      for (int iat=0; iat<n_residue_atoms; iat++) {
	 std::string atom_name(residue_atoms[iat]->name);
	 if (atom_name == ring_atom_names[iring_at]) {
	    clipper::Coord_orth at_pt(residue_atoms[iat]->x,
				      residue_atoms[iat]->y,
				      residue_atoms[iat]->z);
	    aromatic_plane_points.push_back(at_pt);
	    break;
	 }
      }
   }
   if (aromatic_plane_points.size() != ring_atom_names.size()) {
      coot::residue_spec_t spec(res_ref);
      std::string mess = "Not all aromatic atoms were found in residue ";
      mess += spec.chain; 
      mess += " "; 
      mess += spec.resno; 
      throw std::runtime_error(mess);
   } else {
      // can throw an exception
      std::pair<clipper::Coord_orth, clipper::Coord_orth> p =
	 ring_centre_and_normal(aromatic_plane_points);
      double plane_offset = 1; // Angstroms plane normal offset, Clark & Labute 2007
      pt_1 = p.first + plane_offset * p.second;
      pt_2 = p.first - plane_offset * p.second;
      if (0) { 
	 std::cout << "   Centre: " << p.first.x() << " " << p.first.y() << " " << p.first.z() << std::endl;
	 std::cout << "   Normal: " << p.second.x() << " " << p.second.y() << " " << p.second.z() << std::endl;
	 std::cout << " Pi point: " << pt_1.x() << " " << pt_1.y() << " " << pt_1.z() << std::endl;
      }
   }
   return std::pair<clipper::Coord_orth, clipper::Coord_orth> (pt_1, pt_2);
}

// can throw an exception if not enough points found in pts.
//
std::pair<clipper::Coord_orth, clipper::Coord_orth>
coot::pi_stacking_container_t::ring_centre_and_normal(const std::vector<clipper::Coord_orth> &pts) const {

   clipper::Coord_orth centre(0,0,0);
   clipper::Coord_orth normal(0,0,0);

   if (pts.size() < 3) {
      std::string mess = "in ring_centre_and_normal() not enough point to calculate lsq plane";
      throw runtime_error(mess);
   }

   double sum_x = 0, sum_y = 0, sum_z = 0;
   for (unsigned int ipl=0; ipl<pts.size(); ipl++) {
      sum_x += pts[ipl].x();
      sum_y += pts[ipl].y();
      sum_z += pts[ipl].z();
   }

   double divisor = 1.0/double(pts.size());
   centre = clipper::Coord_orth(sum_x*divisor, sum_y*divisor, sum_z*divisor);
   clipper::Matrix<double> mat(3,3);
   for (unsigned int ipl=0; ipl<pts.size(); ipl++) {
      mat(0,0) += (pts[ipl].x()-centre.x()) * (pts[ipl].x()-centre.x());
      mat(1,1) += (pts[ipl].y()-centre.y()) * (pts[ipl].y()-centre.y());
      mat(2,2) += (pts[ipl].z()-centre.z()) * (pts[ipl].z()-centre.z());
      mat(0,1) += (pts[ipl].x()-centre.x()) * (pts[ipl].y()-centre.y());
      mat(0,2) += (pts[ipl].x()-centre.x()) * (pts[ipl].z()-centre.z());
      mat(1,2) += (pts[ipl].y()-centre.y()) * (pts[ipl].z()-centre.z());
   }
   mat(1,0) = mat(0,1);
   mat(2,0) = mat(0,2);
   mat(2,1) = mat(1,2);
   std::vector<double> eigens = mat.eigen(false);
   if (0) 
      std::cout << "   we get eigen values: "
		<< eigens[0] << "  "
		<< eigens[1] << "  "
		<< eigens[2] << std::endl;

   int eigen_index = 0;
   if (eigens[1] < eigens[eigen_index])
      eigen_index = 1;
   if (eigens[2] < eigens[eigen_index])
      eigen_index = 2;

   clipper::Coord_orth eigen_vec(mat(0, eigen_index),
				 mat(1, eigen_index),
				 mat(2, eigen_index));

   double sum_sq = 1e-20;
   sum_sq += eigen_vec.x() * eigen_vec.x();
   sum_sq += eigen_vec.y() * eigen_vec.y();
   sum_sq += eigen_vec.z() * eigen_vec.z();

   normal = clipper::Coord_orth(eigen_vec.x()/sum_sq,
				eigen_vec.y()/sum_sq,
				eigen_vec.z()/sum_sq);

   return std::pair<clipper::Coord_orth, clipper::Coord_orth> (centre, normal);
}

// TRP (for example) has 2 rings, so we have to return a vector
//
// Don't forget RNA and DNA! FIXME
//
std::vector<std::vector<std::string> >
coot::pi_stacking_container_t::ring_atom_names(const std::string &residue_name) const {

   std::vector<std::vector<std::string> > v_outer;
   if ((residue_name == "PHE") || (residue_name == "TYR") || (residue_name == "PTY")) {
      std::vector<std::string> v;
      v.push_back(" CG ");
      v.push_back(" CZ ");
      v.push_back(" CD1");
      v.push_back(" CD2");
      v.push_back(" CE1");
      v.push_back(" CE2");
      v_outer.push_back(v);
   }
   if (residue_name == "TRP") {
      std::vector<std::string> v;
      v.push_back(" CG ");
      v.push_back(" CD1");
      v.push_back(" CD2");
      v.push_back(" NE1");
      v.push_back(" CE2");
      v_outer.push_back(v);
      v.clear();
      v.push_back(" CD2");
      v.push_back(" CE2");
      v.push_back(" CE3");
      v.push_back(" CZ3");
      v.push_back(" CZ2");
      v.push_back(" CH2");
      v_outer.push_back(v);
   }

   if ((residue_name == "DA") ||   // or RNA equivalent
       (residue_name == "Ad") ||
       (residue_name == "Ar")) { 
      std::vector<std::string> v;
      v.push_back(" C5 ");
      v.push_back(" C4 ");
      v.push_back(" N9 ");
      v.push_back(" C8 ");
      v.push_back(" N7 ");
      v_outer.push_back(v);
      v.clear();
      v.push_back(" C6 ");
      v.push_back(" C5 ");
      v.push_back(" C4 ");
      v.push_back(" N3 ");
      v.push_back(" C2 ");
      v.push_back(" N1 ");
   }
       
   if ((residue_name == "DG") ||   // or RNA equivalent
       (residue_name == "Gd") ||
       (residue_name == "Gr")) { 
      std::vector<std::string> v;
      v.push_back(" C5 ");
      v.push_back(" C4 ");
      v.push_back(" N9 ");
      v.push_back(" C8 ");
      v.push_back(" N7 ");
      v_outer.push_back(v);
      v.clear();
      v.push_back(" C6 ");
      v.push_back(" C5 ");
      v.push_back(" C4 ");
      v.push_back(" N3 ");
      v.push_back(" C2 ");
      v.push_back(" N1 ");
   }

   if ((residue_name == "DC") ||   // or RNA equivalent
       (residue_name == "Cd") ||
       (residue_name == "Cr")) { 
      std::vector<std::string> v;
      v.push_back(" N1 ");
      v.push_back(" C2 ");
      v.push_back(" N3 ");
      v.push_back(" C4 ");
      v.push_back(" C5 ");
      v.push_back(" C6 ");
      v_outer.push_back(v);
   }
   
   if ((residue_name == "DT") ||   // or RNA equivalent
       (residue_name == "Td") ||
       (residue_name == "Tr")) { 
      std::vector<std::string> v;
      v.push_back(" N1 ");
      v.push_back(" C2 ");
      v.push_back(" N3 ");
      v.push_back(" C4 ");
      v.push_back(" C5 ");
      v.push_back(" C6 ");
      v_outer.push_back(v);
   }
       
   return v_outer;
}


std::vector<clipper::Coord_orth>
coot::pi_stacking_container_t::get_cation_atom_positions(CResidue *res) const {
   
   std::vector<clipper::Coord_orth> v;

   std::string res_name(res->GetResName());

   if (res_name == "LYS") {
      PPCAtom residue_atoms = 0;
      int n_residue_atoms;
      res->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int i=0; i<n_residue_atoms; i++) {
	 std::string atom_name(residue_atoms[i]->name);
	 if (atom_name == " NZ ") {
	    clipper::Coord_orth pt(residue_atoms[i]->x,
				   residue_atoms[i]->y,
				   residue_atoms[i]->z);
	    v.push_back(pt);
	 }
      }
   }

   if (res_name == "ARG") {
      PPCAtom residue_atoms = 0;
      int n_residue_atoms;
      res->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int i=0; i<n_residue_atoms; i++) {
	 std::string atom_name(residue_atoms[i]->name);
	 if ((atom_name == " NH1") ||
	     (atom_name == " NH2")) {
	    clipper::Coord_orth pt(residue_atoms[i]->x,
				   residue_atoms[i]->y,
				   residue_atoms[i]->z);
	    v.push_back(pt);
	 }
      }
   }
   

   return v;
}

float
coot::pi_stacking_container_t::overlap_of_cation_pi(const clipper::Coord_orth &ligand_pi_point,
						    const clipper::Coord_orth &cation_atom_point) const {

   // multipliers follow the order of the atom/control points
   // 
   float score = overlap_of_pi_spheres(ligand_pi_point, cation_atom_point, 0.78, -1, 7.8, -0.05);

   // std::cout << "overlap_of_cation_pi() calling overlap_of_pi_spheres() returns "
   // << score << std::endl;
   
   return score;
}
 

void
coot::write_fle_centres(const std::vector<fle_residues_helper_t> &v,
			const std::vector<coot::fle_ligand_bond_t> &bonds_to_ligand,
			const std::vector<coot::solvent_exposure_difference_helper_t> &sed,
			const coot::pi_stacking_container_t &stack_info,
			CResidue *res_flat) {

   if (0) {
      std::cout << "DEBUG:: in write_fle_centres() " << stack_info.stackings.size()
		<< " stackings:\n";
      for (unsigned int i=0; i<stack_info.stackings.size(); i++)
      std::cout << "    " << i << ":  " << stack_info.stackings[i] << std::endl;
      std::cout << std::endl;
   }
   
   std::string file_name = "coot-tmp-fle-view-residue-info.txt";
   std::ofstream of(file_name.c_str());
   if (!of) {
      std::cout << "failed to open output file " << file_name << std::endl;
   } else {
      
      for (unsigned int i=0; i<v.size(); i++) {
	 if (v[i].is_set) { 
	    of << "RES " << v[i].centre.x() << " " << v[i].centre.y() << " "
	       << v[i].centre.z() << " "
	       << v[i].residue_name << " " << v[i].spec.chain << v[i].spec.resno;

	    // Add the qualifying distance to protein if this is a water atom.
	    // 
	    if (v[i].residue_name == "HOH") { 
	       for (unsigned int ib=0; ib<bonds_to_ligand.size(); ib++) { 
		  if (bonds_to_ligand[ib].res_spec == v[i].spec) {
		     of << " water->protein-length " << bonds_to_ligand[ib].water_protein_length;
		     break;
		  }
	       }
	       of << "\n";
	    } else { 
	       of << "\n";
	    }


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

	    // Is there a stacking interaction associated with this residue?
	    //
	    // There could, in fact, be several stackings associated
	    // with this residue.  We want to pick the best one.  So
	    // we need to run through the list of stacking first,
	    // looking for the best.
	    //
	    // 
	    double best_stacking_score = -1.0;
	    coot::pi_stacking_instance_t best_stacking(NULL, "");
	    
	    for (unsigned int istack=0; istack<stack_info.stackings.size(); istack++) {
	       coot::residue_spec_t spec(stack_info.stackings[istack].res);
	       if (spec == v[i].spec) {
		  if (stack_info.stackings[istack].overlap_score > best_stacking_score) {
		     best_stacking_score = stack_info.stackings[istack].overlap_score;
		     best_stacking = stack_info.stackings[istack];
		  }
	       }
	    }

	    if (best_stacking_score > 0.0) {
	       std::string type = "    pi-pi"; // leading spaces
					       // important in format
					       // of atom names

	       // Recall that CATION_PI_STACKING is for ligand ring
	       // systems, PI_CATION_STACKING is for protein residue
	       // ring systems.
	       // 
	       if (best_stacking.type == coot::pi_stacking_instance_t::CATION_PI_STACKING)
		  type = "cation-pi";
	       if (best_stacking.type == coot::pi_stacking_instance_t::PI_CATION_STACKING)
		  type = "pi-cation";
	       if (0) // debug
		  of << "# Best stacking for RES " 
		     << v[i].residue_name << " " << v[i].spec.chain
		     << v[i].spec.resno << " is " << best_stacking.type << " "
		     << best_stacking_score << std::endl;
	       of << "STACKING " << type << " ";
	       switch (best_stacking.type) {
	       case coot::pi_stacking_instance_t::PI_PI_STACKING: 
		  for (unsigned int jat=0;
		       jat<best_stacking.ligand_ring_atom_names.size();
		       jat++) {
		     of << best_stacking.ligand_ring_atom_names[jat] << "  ";
		  }
		  break;
		  
	       case coot::pi_stacking_instance_t::CATION_PI_STACKING:
		  of << best_stacking.ligand_cationic_atom_name;
		  break;
		  
	       case coot::pi_stacking_instance_t::PI_CATION_STACKING:
		  for (unsigned int jat=0;
		       jat<best_stacking.ligand_ring_atom_names.size();
		       jat++) {
		     of << best_stacking.ligand_ring_atom_names[jat] << "  ";
		  }
		  break;
	       }
	       of << "\n";
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

// This function is currently not used.
// 
// if there is a name map, use it, otherwise just bond to the found
// atom.  Using the name map allows bond to hydrogens hanging off of
// ligand atoms (e.g. the H on an O or an H on an N). The Hs are not
// in res_ref (often/typically), but are in the flat residue.
// 
std::vector<coot::fle_ligand_bond_t>
coot::get_fle_ligand_bonds(CResidue *ligand_res,
			   const std::vector<CResidue *> &residues,
			   const std::map<std::string, std::string> &name_map) {
   std::vector<coot::fle_ligand_bond_t> v;

   // do it by hand, rather than create a molecule and use SeekContacts();
   double bond_length = 3.3;

   double bl_2 = bond_length * bond_length;
   PPCAtom ligand_residue_atoms = 0;
   int n_ligand_residue_atoms;
   ligand_res->GetAtomTable(ligand_residue_atoms, n_ligand_residue_atoms);
   for (unsigned int iat=0; iat<n_ligand_residue_atoms; iat++) {
      std::string ele_ligand = ligand_residue_atoms[iat]->element;
      if (ele_ligand != " C") { 
	 clipper::Coord_orth ref_pt(ligand_residue_atoms[iat]->x,
				    ligand_residue_atoms[iat]->y,
				    ligand_residue_atoms[iat]->z);
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
		     std::string ligand_atom_name = ligand_residue_atoms[iat]->name;
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

		     // We don't want to map between hydroxyl oxygens -> it attached H
		     // when the residue to which it has a bond is a metal!  (and in
		     // that case, we should probably remove the Hs from the ligand
		     // anyway - but not today).
		     //
		     if (is_a_metal(residues[ires])) {
			bond_type = fle_ligand_bond_t::METAL_CONTACT_BOND;
		     } else {
			// note, can't use [] because name_map is const
			std::map<std::string, std::string>::const_iterator it =
			   name_map.find(ligand_atom_name);
			if (it != name_map.end()) {
			   // If the map happens, that's presumably because we found a H
			   // attached to an N (or an H attached to an O), either way, we
			   // are sitting now on an H.
			   ligand_atom_name = it->second;

			   // is this OK?
			   // 
			   ele = " H";
			   
			   bond_type = fle_ligand_bond_t::H_BOND_ACCEPTOR_SIDECHAIN;
			   if (coot::is_main_chain_p(residue_atoms[jat]))
			      bond_type = fle_ligand_bond_t::H_BOND_ACCEPTOR_MAINCHAIN;
			}
		     }

		     // that was not enough.  We still need to do some
		     // more rejections.  Here's a start:
		     if (! ((ele_ligand == " O") && (ele == " O"))) {
			
			coot::fle_ligand_bond_t bond(ligand_atom_name, bond_type, bl, res_spec);
			v.push_back(bond);
			break;
		     }
		  }
	       }
	    }
	 }
      }
   }
   return v;
}

// Use coot::h_bonds class to generate ligands.  We do that by creating a synthetic
// temporary  molecule and atom selections.
// 
std::vector<coot::fle_ligand_bond_t>
coot::get_fle_ligand_bonds(CResidue *ligand_res,
			   const std::vector<CResidue *> &residues,
			   CMMDBManager *mol, 
			   const std::map<std::string, std::string> &name_map,
			   const coot::protein_geometry &geom) {
   std::vector<coot::fle_ligand_bond_t> v; // returned value

   bool debug = 0;
   std::vector<CResidue *> rv = residues;
   rv.push_back(ligand_res);

   std::pair<bool, CMMDBManager *> m = coot::util::create_mmdbmanager_from_residue_vector(rv);
   coot::residue_spec_t ligand_spec(ligand_res);

   if (m.first) { 
      int SelHnd_all = m.second->NewSelection();
      int SelHnd_lig = m.second->NewSelection();
      m.second->SelectAtoms(SelHnd_all, 0, "*", ANY_RES, "*", ANY_RES, "*", "*", "*", "*", "*");
      m.second->SelectAtoms(SelHnd_lig, 0, ligand_spec.chain.c_str(),
			    ligand_spec.resno, ligand_spec.insertion_code.c_str(),
			    ligand_spec.resno, ligand_spec.insertion_code.c_str(),
			    "*", "*", "*", "*");

      // -----------------------
      //   hydrogen bonds 
      // -----------------------
      
      coot::h_bonds hb;
      std::vector<coot::h_bond> hbonds = hb.get(SelHnd_lig, SelHnd_all, m.second, geom);

      std::cout << "DEBUG:: get_fle_ligand_bonds from h_bonds class found "
		<< hbonds.size() << " bonds." << std::endl;

      for (unsigned int i=0; i<hbonds.size(); i++) { 
	 std::cout << coot::atom_spec_t(hbonds[i].donor) << "..."
		   << coot::atom_spec_t(hbonds[i].acceptor) << " with ligand donor flag "
		   << hbonds[i].ligand_atom_is_donor << std::endl;
	 std::string ligand_atom_name = hbonds[i].acceptor->name;
	 coot::residue_spec_t res_spec(hbonds[i].donor);
	 CResidue *ligand_residue = hbonds[i].donor->residue;
	 int bond_type = fle_ligand_bond_t::get_bond_type(hbonds[i].donor,
							  hbonds[i].acceptor,
							  hbonds[i].ligand_atom_is_donor);

	 if (hbonds[i].ligand_atom_is_donor) {
	    ligand_atom_name = hbonds[i].donor->name;
	    res_spec = coot::residue_spec_t(hbonds[i].acceptor);
	    ligand_residue = hbonds[i].acceptor->residue;
	 }

	 // Now, in 3D (pre-prodrgification) we don't have (polar) Hs on the ligand
	 // (but we do in 2D), the map allows transfer from the ligand O or N to the
	 // polar H in FLEV.
	 // 
	 std::map<std::string, std::string>::const_iterator it = name_map.find(ligand_atom_name);
	 if (it != name_map.end()) {
	    // If the map happens, that's presumably because we found a H
	    // attached to an N (or an H attached to an O), either way, we
	    // are sitting now on an H.
	    ligand_atom_name = it->second;
	 }

	 if (debug) 
	    std::cout << "constructing fle ligand bond " << ligand_atom_name
		      << " " << bond_type << " " << hbonds[i].dist << " " 
		      << res_spec << std::endl;
	 
	 coot::fle_ligand_bond_t bond(ligand_atom_name, bond_type, hbonds[i].dist, res_spec);
	 std::string residue_name = ligand_residue->GetResName();
	 if (residue_name == "HOH")
	    bond.water_protein_length = find_water_protein_length(ligand_residue, mol);
	 v.push_back(bond);
      }

      // -----------------------
      //   covalent bonds 
      // -----------------------

      std::vector<coot::fle_ligand_bond_t> covalent_bonds =
	 get_covalent_bonds(m.second, SelHnd_lig, SelHnd_all, ligand_spec, geom);
      for (unsigned int i=0; i<covalent_bonds.size(); i++)
	 v.push_back(covalent_bonds[i]);


      // -----------------------
      //   metal bonds 
      // -----------------------

      std::vector<coot::fle_ligand_bond_t> metal_bonds = get_metal_bonds(ligand_res, residues);
      for (unsigned int i=0; i<metal_bonds.size(); i++)
	 v.push_back(metal_bonds[i]);

      
      // -----------------------
      //   clean up 
      // -----------------------
      
      m.second->DeleteSelection(SelHnd_lig);
      m.second->DeleteSelection(SelHnd_all);
      delete m.second;

   }

   return v;
}

std::vector<coot::fle_ligand_bond_t>
coot::get_covalent_bonds(CMMDBManager *mol,
			 int SelHnd_lig,
			 int SelHnd_all,
			 const residue_spec_t &ligand_spec,
			 const protein_geometry &geom) {
   
   std::vector<coot::fle_ligand_bond_t> v;
   int SelHnd_local = mol->NewSelection();
   mol->SelectAtoms(SelHnd_local, 0, "*", ANY_RES, "*", ANY_RES, "*", "*", "*", "*", "*");
   mol->Select(SelHnd_local, STYPE_ATOM, 0, ligand_spec.chain.c_str(),
	       ligand_spec.resno, ligand_spec.insertion_code.c_str(),
	       ligand_spec.resno, ligand_spec.insertion_code.c_str(),
	       "*", "*", "*", "*", SKEY_XOR);

   // now find contacts:
   // 
   PSContact pscontact = NULL;
   int n_contacts;
   long i_contact_group = 1;
   mat44 my_matt;
   CSymOps symm;
   for (int i=0; i<4; i++) 
      for (int j=0; j<4; j++) 
	 my_matt[i][j] = 0.0;      
   for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

   PPCAtom lig_atom_selection = 0;
   int n_lig_atoms;
   mol->GetSelIndex(SelHnd_lig, lig_atom_selection, n_lig_atoms);

   PPCAtom other_atom_selection = 0;
   int n_other_atoms;
   mol->GetSelIndex(SelHnd_local, other_atom_selection, n_other_atoms);

   realtype min_dist = 0.1;
   realtype max_dist = 2.3; //  even S-S is shorter than this, I think
   realtype cno_max_dist = 1.8;

   mol->SeekContacts(  lig_atom_selection,   n_lig_atoms,
		     other_atom_selection, n_other_atoms,
		     min_dist, max_dist, // min, max distances
		     0,        // seqDist 0 -> in same res also
		     pscontact, n_contacts,
		       0, &my_matt, i_contact_group);

   std::vector<std::pair<CResidue *, CResidue *> > contacting_pairs_vec;

   if (n_contacts > 0) {
      if (pscontact) {
	 for (int i=0; i<n_contacts; i++) {
	    CAtom *at_1 =   lig_atom_selection[pscontact[i].id1];
	    CAtom *at_2 = other_atom_selection[pscontact[i].id2];
// 	    std::cout << "DEUBG:: Covalent test "
// 		      << coot::atom_spec_t(at_1) << "..."
// 		      << coot::atom_spec_t(at_2) << std::endl;
	    std::pair<CResidue *, CResidue *> pair(at_1->GetResidue(),
						   at_2->GetResidue());

	    std::string ele_1 = at_1->element;
	    std::string ele_2 = at_2->element;
	    clipper::Coord_orth pt_1(at_1->x, at_1->y, at_1->z);
	    clipper::Coord_orth pt_2(at_2->x, at_2->y, at_2->z);
	    double d = (pt_1-pt_2).lengthsq();

	    double dist_for_bond = max_dist;
	    if (((ele_1 == " C") || (ele_1 == " N") || (ele_1 == "O")) &&
		((ele_2 == " C") || (ele_2 == " N") || (ele_2 == "O")))
	       dist_for_bond = cno_max_dist;

	    if (d < dist_for_bond) { 

	       // only add this pair if it is not already in the list:
	       if (std::find(contacting_pairs_vec.begin(), contacting_pairs_vec.end(), pair) ==
		   contacting_pairs_vec.end()) {
		  contacting_pairs_vec.push_back(pair);
		  coot::residue_spec_t res_spec(at_2->GetResidue());
		  int bond_type = coot::fle_ligand_bond_t::BOND_COVALENT;
		  std::string ligand_atom_name(at_1->name);
		  coot::fle_ligand_bond_t bond(ligand_atom_name, bond_type, d, res_spec);
		  v.push_back(bond);
	       }
	    }
	 }
      }
   }

   mol->DeleteSelection(SelHnd_local);
   return v;
}

std::vector<coot::fle_ligand_bond_t>
coot::get_metal_bonds(CResidue *ligand_residue, const std::vector<CResidue *> &residues) {

   // a non-Hydrogen, non-Carbon ligand atom, that is.
   double max_dist_metal_to_ligand_atom = 3.5; // pass this parameter?
   
   std::vector<coot::fle_ligand_bond_t> v;

   double best_dist_sqrd = max_dist_metal_to_ligand_atom * max_dist_metal_to_ligand_atom;
   std::string ligand_atom_name; // goes with best_dist_sqrd
   
   PPCAtom ligand_residue_atoms = 0;
   int n_ligand_residue_atoms;
   ligand_residue->GetAtomTable(ligand_residue_atoms, n_ligand_residue_atoms);
   for (unsigned int i=0; i<residues.size(); i++) { 
      if (coot::is_a_metal(residues[i])) {
	 PPCAtom residue_atoms = 0;
	 int n_residue_atoms;
	 residues[i]->GetAtomTable(residue_atoms, n_residue_atoms);
	 for (unsigned int irat=0; irat<n_residue_atoms; irat++) {
	    for (unsigned int ilat=0; ilat<n_residue_atoms; ilat++) {
	       std::string ele(residue_atoms[irat]->element);
	       if ((ele == " H") || (ele == " C")) { 
		  clipper::Coord_orth pt_1(ligand_residue_atoms[ilat]->x,
					   ligand_residue_atoms[ilat]->y,
					   ligand_residue_atoms[ilat]->z);
		  clipper::Coord_orth pt_2(residue_atoms[irat]->x,
					   residue_atoms[irat]->y,
					   residue_atoms[irat]->z);
		  double d2 = (pt_1-pt_2).clipper::Coord_orth::lengthsq();
		  if (d2 < best_dist_sqrd) {
		     best_dist_sqrd = d2;
		     ligand_atom_name = ligand_residue_atoms[ilat]->name;
		  }
	       }
	    }
	 }
	 if (best_dist_sqrd < max_dist_metal_to_ligand_atom * max_dist_metal_to_ligand_atom) {
	    coot::fle_ligand_bond_t bond(ligand_atom_name,
					 coot::fle_ligand_bond_t::METAL_CONTACT_BOND,
					 sqrt(best_dist_sqrd),
					 coot::residue_spec_t(residues[i]));
	    v.push_back(bond);
	 }
      }
   }

   return v;
}


// should be in coot-utils perhaps?
//
bool
coot::is_a_metal(CResidue *res) {

   bool r = 0;
   std::string res_name = res->GetResName();
   if (res_name == "MG")
      return 1;
   if (res_name == "CA")
      return 1;
   if (res_name == "MN")
      return 1;
   if (res_name == "FE")
      return 1;
   if (res_name == "FE")
      return 1;
   if (res_name == "K")
      return 1;
   if (res_name == "NA")
      return 1;
   if (res_name == "CO")
      return 1;
   if (res_name == "NI")
      return 1;
   if (res_name == "CU")
      return 1;
   if (res_name == "ZN")
      return 1;
   if (res_name == "RU")
      return 1;
   if (res_name == "PT")
      return 1;
   if (res_name == "AU")
      return 1;
   if (res_name == "AG")
      return 1;

   return 0;

}


// return 100 if no other contact found (strange!)
// 
double
coot::find_water_protein_length(CResidue *ligand_residue, CMMDBManager *mol) {

   double dist = 100;

   double dist_sqrd = dist * dist;
   double dist_sqrd_init = dist_sqrd;
   
   PPCAtom ligand_residue_atoms = 0;
   int n_ligand_residue_atoms;
   ligand_residue->GetAtomTable(ligand_residue_atoms, n_ligand_residue_atoms);


   int imod = 1;
   CModel *model_p = mol->GetModel(imod);
   CChain *chain_p;
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      CResidue *residue_p;
      for (int ires=0; ires<nres; ires++) { 
	 residue_p = chain_p->GetResidue(ires);
	 if (ligand_residue != residue_p) {
	    std::string residue_name(residue_p->GetResName());
	    if (residue_name != "HOH") { 
	       PPCAtom residue_atoms = 0;
	       int n_residue_atoms;
	       residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
	       for (unsigned int il=0; il<n_ligand_residue_atoms; il++) {
		  for (unsigned int irat=0; irat<n_residue_atoms; irat++) {
		     std::string ele(residue_atoms[irat]->element);
		     if ((ele == " O") || (ele == " N")) { 
			clipper::Coord_orth pt_1(ligand_residue_atoms[il]->x,
						 ligand_residue_atoms[il]->y,
						 ligand_residue_atoms[il]->z);
			clipper::Coord_orth pt_2(residue_atoms[irat]->x,
						 residue_atoms[irat]->y,
						 residue_atoms[irat]->z);
			double d2 = (pt_1-pt_2).clipper::Coord_orth::lengthsq();
			if (d2 < dist_sqrd) {
			   dist_sqrd = d2;
			}
		     }
		  }
	       }
	    }
	 }
      }
   }

   if (dist_sqrd < dist_sqrd_init)  // usually is.
      dist = sqrt(dist_sqrd);

   return dist;
}

   

bool
coot::standard_residue_name_p(const std::string &rn) {

   std::vector<std::string> s;

   s.push_back("ALA");
   s.push_back("ASP");
   s.push_back("ASN");
   s.push_back("CYS");
   s.push_back("GLN");
   s.push_back("GLY");
   s.push_back("GLU");
   s.push_back("PHE");
   s.push_back("HIS");
   s.push_back("ILE");
   s.push_back("LYS");
   s.push_back("LEU");
   s.push_back("MET");
   s.push_back("MSE");
   s.push_back("PRO");
   s.push_back("ARG");
   s.push_back("SER");
   s.push_back("THR");
   s.push_back("VAL");
   s.push_back("TRP");
   s.push_back("TYR");

   s.push_back("AR");
   s.push_back("AD");
   s.push_back("CR");
   s.push_back("CD");
   s.push_back("GR");
   s.push_back("GD");
   s.push_back("TD");

   if (std::find(s.begin(), s.end(), rn) == s.end())
      return 0;
   else
      return 1;
}

// The substitution contour:
// 
// The substitution contour depends on the hydridisation of the atoms
// of the ligands.  We need to "fire cannonballs" along the direction
// of the hydrogens and see when they bash into the atoms of the
// protein.
//
// Most (I guess) hydrogens of ligands are "riding", by that, I mean
// that the position of the hydrogen is determined by the geometry of
// the atom to which it is attached (e.g. the H on an N that has 2
// single bonds to 2 carbons).
//
// So, PRODRG will give us a ligand that has hydrogens on it - good
// show (we will need to run it in "regular" mode ("MINI PREP", say)
// to generate such hydrogens on the model).
//
// We need to determine which are the riding hydrogens and which are
// the rotatable ones. The riding hydrogens give us the cannonball
// vector trivially. For the rotable ones, we should draw a number
// from a probability distribution defined by the parameters in the
// cif dictionary.
// 
// How do we know which are riding and which are not?
//
// Here is my guess: riding hydrogens do not have torsions for them in
// the dictionary.  Rotable hydrogens do.  Note though, that PRODRG
// hydrogen names (in the cif file) are problematic.  You can't rely
// on name matching of hydrogens between the PDB and the dictionary.
//
// So here is what we have to do:
//
// i) From the dictionary, look through the list of non-Hydrogen atoms
// to see if there are hydrogens bonded to the atom.
// 
//    ii) If so, does this non-Hydrogen atom appear as atom 2 in a 
//    torsion restraint with a hydrogen as atom 1 or appear as atom 3
//    in a torsion restraint with a hydrogen as atom 4?
//
//        iii) If so, then this is a torsionable hydrogen, we need to
//        generate vectors by random sampling from the probability
//        distribution of this torsion.
//
//        However, if it is not, then this is a riding hydrogen.  Add
//        the vector from the ligand atom to the hydrogen as a
//        cannonball vector associated with that ligand atom.
//
// OK, at the end of the hybridisation analysis, we want a vector of
// ligand "heavy" (in this case, non-Hydrogen) atoms and associated
// with each atom is a vector (i.e. list) of cannonball vectors
// (i.e. directions).  It is not yet clear to me if we need anything
// else, so maybe a vector of pairs of ligand atoms and a vector of
// Coord_orths will do.
//
// This vector of pairs will in turn be processed into a "bash"
// distance for each atom (how far can we go from the atom before
// bashing in the atoms of the protein?)  Note that, in doing so, we
// need to be able to specify "infinite bash distance". Let's make a
// trivial class for bash distance.
//
// Then these distances get overlayed onto the ligand grid.
//
// The the ligand grid gets contoured.
//
// And plotted to the canvas.

std::vector<std::pair<std::string, int> >
coot::get_prodrg_hybridizations(const std::string &prodrg_log_file_name) {

   std::cout << "in get_prodrg_hybridizations() "  << std::endl;
   std::vector<std::pair<std::string, int> > v;
   std::ifstream f(prodrg_log_file_name.c_str());
   if (! f) {
      std::cout << "WARNING:: failed to open " << prodrg_log_file_name << std::endl;
   } else {
      std::string line;
      bool hybridizations_flag = 0;
      while (std::getline(f, line)) {
	 std::vector<std::string> words =
	    coot::util::split_string_no_blanks(line, " ");
	 if (words.size() > 1) {
	    if (words[0] == "DETERMINING") {
	       if (words[1] == "HYBRIDISATION") {
		  hybridizations_flag = 1;
	       }
	    }
	    if (words[0] == "PRODRG>") {
	       hybridizations_flag = 0;
	    }
	 }
	 if (hybridizations_flag) {
	    if (line.find("qual = ") != std::string::npos ) {
	       std::string atom_name=line.substr(1,4);
	       std::string hybridization_str = line.substr(13,3);
	       std::cout << ":" << atom_name <<": :" << hybridization_str
			 << ":" << std::endl;
	       int hybridization = coot::SP_HYBRIDIZATION;
	       if (hybridization_str == "sp3")
		  hybridization = coot::SP3_HYBRIDIZATION;
	       if (hybridization_str == "sp2")
		  hybridization = coot::SP2_HYBRIDIZATION;
	       std::pair<std::string, int> p(atom_name, hybridization);
	       v.push_back(p);
	    }
	 }
      }
   }
   return v;
}

std::vector<std::pair<CAtom *, std::vector<clipper::Coord_orth> > >
coot::get_cannonball_vectors(CResidue *ligand_res_3d,
			     const coot::dictionary_residue_restraints_t &monomer_restraints) {

   std::vector<std::pair<CAtom *, std::vector<clipper::Coord_orth> > > v;

   return v;
}
