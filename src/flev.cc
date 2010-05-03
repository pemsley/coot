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


void fle_view_internal(int imol, const char *chain_id, int res_no, const char *ins_code, 
		       int imol_ligand_fragment, 
		       const char *prodrg_output_flat_mol_file_name,
		       const char *prodrg_output_flat_pdb_file_name,
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


		  handle_cif_dictionary(prodrg_output_dict_cif_file_name);
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

   float pi_overlap_thresh = 0.2;
   for (unsigned int iring=0; iring<ring_list.size(); iring++) {
      try {
	 std::pair<clipper::Coord_orth, clipper::Coord_orth> ligand_ring_pi_pts = 
	    get_ring_pi_centre_points(ring_list[iring], res_ref);
	 
	 for (unsigned int ires=0; ires<residues.size(); ires++) {
	    std::string res_name(residues[ires]->GetResName());

	    // Don't forget DNA and RNA!
	    // 
	    if (1) {
	       std::cout << "==== Environment residue " << coot::residue_spec_t(residues[ires])
			 << " " << res_name << std::endl;
	       double pi_overlap_1 = get_pi_overlap(residues[ires], ligand_ring_pi_pts.first);
	       double pi_overlap_2 = get_pi_overlap(residues[ires], ligand_ring_pi_pts.second);
	       std::cout << "   Overlaps: " << pi_overlap_1 << " " << pi_overlap_2
			 << std::endl;
	       if (pi_overlap_1 > pi_overlap_thresh) {
		  coot::pi_stacking_instance_t st(residues[ires],
						  coot::pi_stacking_instance_t::PI_PI_STACKING,
						  ring_list[iring]);
		  stackings.push_back(st);
	       }
	       if (pi_overlap_2 > pi_overlap_thresh) {
		  coot::pi_stacking_instance_t st(residues[ires],
						  coot::pi_stacking_instance_t::PI_PI_STACKING,
						  ring_list[iring]);
		  stackings.push_back(st);
	       }
	    }
	 }
      }
      catch (std::runtime_error rte) {
	 std::cout << "WARNING:: " << rte.what() << std::endl;
      }
   }
}

// can throw an exception because it calls get_ring_pi_centre_points()
//
// should return the stacking type, e.g. PI_CATION_STACKING.
// 
float
coot::pi_stacking_container_t::get_pi_overlap(CResidue *res,
					      clipper::Coord_orth &ligand_pi_point) const {

   float score = 0;
   std::string res_name(res->GetResName());
   std::vector<std::vector<std::string> > atom_names = ring_atom_names(res_name);
   for (unsigned int iring=0; iring<atom_names.size(); iring++) {
      std::pair<clipper::Coord_orth, clipper::Coord_orth> residue_pi_points =
	 get_ring_pi_centre_points(atom_names[iring], res);
      float score_1 = overlap_of_pi_spheres(ligand_pi_point, residue_pi_points.first);
      float score_2 = overlap_of_pi_spheres(ligand_pi_point, residue_pi_points.second);
      if (score_1 > score)
	 score = score_1;
      if (score_2 > score)
	 score = score_2;
   }
   return score;
}

float
coot::pi_stacking_container_t::overlap_of_pi_spheres(const clipper::Coord_orth &pt_1,
						     const clipper::Coord_orth &pt_2) const {

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
	       score_1 = 0.78 * exp(-d_1_sqd);
	    if (d_2_sqd < 12.0)
	       score_2 = 0.78 * exp(-d_2_sqd);
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
   return v_outer;
}


void
coot::write_fle_centres(const std::vector<fle_residues_helper_t> &v,
			const std::vector<coot::fle_ligand_bond_t> &bonds_to_ligand,
			const std::vector<coot::solvent_exposure_difference_helper_t> &sed,
			const coot::pi_stacking_container_t &stack_info,
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

	    // Is there a stacking interaction associated with this residue?
	    for (unsigned int istack=0; istack<stack_info.stackings.size(); istack++) {
	       coot::residue_spec_t spec(stack_info.stackings[istack].res);
	       std::string type = "    pi-pi"; // leading spaces important in format of atom names
	       if (stack_info.stackings[istack].type ==
		   coot::pi_stacking_instance_t::PI_CATION_STACKING)
  		           type = "pi-cation";
	       if (spec == v[i].spec) {
		  of << "STACKING " << type << " ";
		  for (unsigned int jat=0;
		       jat<stack_info.stackings[istack].ligand_ring_atom_names.size();
		       jat++) {
		     of << stack_info.stackings[istack].ligand_ring_atom_names[jat]
			<< "  ";
		  }
		  of << "\n";
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
