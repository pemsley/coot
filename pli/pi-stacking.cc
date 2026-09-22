/* lbg/pi-stacking.cc
 * 
 * Copyright 2010 by the University of Oxford
 * Copyright 2014 by Medical Research Council
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


#include <cmath>
#include <cstring> // otherwise strchr() problems when using clang/Mac,
                   // when including mmdb2/mmdb_manager.h at the top (or
                   // after Python.h?)

#include <mmdb2/mmdb_manager.h>
#include <clipper/core/coords.h>

#include "geometry/protein-geometry.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "geometry/residue-and-atom-specs.hh"
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include "lidia-core/rdkit-interface.hh"
#endif

#include "pi-stacking.hh"

pli::pi_stacking_container_t::pi_stacking_container_t(const coot::dictionary_residue_restraints_t &monomer_restraints,
						       const std::vector<mmdb::Residue *> &residues,
						       mmdb::Residue *res_ref) {

   // get a list of aromatic bonds, so that they can be used to find
   // aromatic rings.
   std::vector<std::vector<std::string> > aromatic_ring_list_atom_names =
      get_aromatic_ring_list(monomer_restraints);

   init(monomer_restraints, residues, res_ref, aromatic_ring_list_atom_names);

}


#ifdef MAKE_ENHANCED_LIGAND_TOOLS
pli::pi_stacking_container_t::pi_stacking_container_t(const coot::dictionary_residue_restraints_t &monomer_restraints,
						       const std::vector<mmdb::Residue *> &residues,
						       mmdb::Residue *res_ref, const RDKit::ROMol &mol) {

   // get a list of aromatic bonds, so that they can be used to find
   // aromatic rings.
   std::vector<std::vector<std::string> > aromatic_ring_list_atom_names =
      get_aromatic_ring_list(monomer_restraints, mol);

   init(monomer_restraints, residues, res_ref, aromatic_ring_list_atom_names);
}
#endif // MAKE_ENHANCED_LIGAND_TOOLS

void
pli::pi_stacking_container_t::init(const coot::dictionary_residue_restraints_t &monomer_restraints,
				    const std::vector<mmdb::Residue *> &residues,
				    mmdb::Residue *res_ref,
				    const std::vector<std::vector<std::string> > &aromatic_ring_list) {

   bool debug = false;

   // float pi_overlap_thresh = 0.0015; // play value
   float pi_pi_overlap_thresh = 0.1;    // HLZ in 3LTW
   float pi_cation_overlap_thresh = 30;  // ZZG in 2wot is 27, close but spurious interaction,
                                         // a bit more than that then.

   for (unsigned int iring=0; iring<aromatic_ring_list.size(); iring++) {
      try {
	 std::pair<clipper::Coord_orth, clipper::Coord_orth> ligand_ring_pi_pts = 
	    get_ring_pi_centre_points(aromatic_ring_list[iring], res_ref);

	 if (debug) {
	    std::cout << "========= ligand ring ";
	    for (unsigned int iat=0; iat<aromatic_ring_list[iring].size(); iat++)
	       std::cout << aromatic_ring_list[iring][iat] << "  ";
	    
	    std::cout << " ====== points " << ligand_ring_pi_pts.first.format() << " "
		      << ligand_ring_pi_pts.second.format() << std::endl;
	 }
	 
	 for (unsigned int ires=0; ires<residues.size(); ires++) {

	    if (debug) {
	       std::string res_name(residues[ires]->GetResName());
	       std::cout << "   ==== Environment residue " << coot::residue_spec_t(residues[ires])
			 << " " << res_name << std::endl;
	    }

	    // score both pi points (above and below the ligand ring) and keep the
	    // better one - one stacking per (ligand ring, residue).
	    overlap_result_t pi_overlap_1 =
	       get_pi_overlap_to_ligand_ring(residues[ires], ligand_ring_pi_pts.first);
	    overlap_result_t pi_overlap_2 =
	       get_pi_overlap_to_ligand_ring(residues[ires], ligand_ring_pi_pts.second);

	    if (debug)
	       std::cout << "    protein cation:ligand ring: Overlaps:  score "
			 << pi_overlap_1.score << " type: " << pi_overlap_1.type << "  score: "
			 << pi_overlap_2.score << " type: " << pi_overlap_2.type << std::endl;

	    const overlap_result_t &best = (pi_overlap_2.score > pi_overlap_1.score) ? pi_overlap_2 : pi_overlap_1;

	    float thresh = -1;
	    if (best.type == pi_stacking_instance_t::PI_PI_STACKING)
	       thresh = pi_pi_overlap_thresh;
	    if (best.type == pi_stacking_instance_t::PI_CATION_STACKING)
	       thresh = pi_cation_overlap_thresh;

	    if (thresh > 0.0 && best.score > thresh) {
	       pi_stacking_instance_t st(residues[ires], best.type, aromatic_ring_list[iring]);
	       st.overlap_score = best.score;
	       st.residue_ring_atom_names = best.residue_ring_atom_names;
	       st.residue_cation_atom_name = best.residue_cation_atom_name;
	       fill_geometry(st, res_ref);
	       if (debug)
		  std::cout << "adding a stacking " << st << std::endl;
	       stackings.push_back(st);
	    } else {
	       // no overlap: is it a T-shaped (edge-to-face) stack?
	       std::vector<std::string> e2f_ring =
		  edge_to_face_ring_by_geometry(aromatic_ring_list[iring], res_ref, residues[ires]);
	       if (! e2f_ring.empty()) {
		  pi_stacking_instance_t st(residues[ires], pi_stacking_instance_t::PI_PI_STACKING,
					    aromatic_ring_list[iring]);
		  st.overlap_score = (best.type == pi_stacking_instance_t::PI_PI_STACKING) ? best.score : 0.0f;
		  st.residue_ring_atom_names = e2f_ring;
		  fill_geometry(st, res_ref);
		  if (debug)
		     std::cout << "adding an edge-to-face stacking " << st << std::endl;
		  stackings.push_back(st);
	       }
	    }
	 }
      }
      catch (const std::runtime_error &rte) {
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
	 overlap_result_t r = get_pi_overlap_to_ligand_cation(residues[ires], cation_points[icat].second);
	 if (r.score > pi_cation_overlap_thresh) {
	    // add a stacking to stackings.
	    pi_stacking_instance_t stacking(residues[ires], cation_points[icat].first);
	    stacking.overlap_score = r.score;
	    stacking.residue_ring_atom_names = r.residue_ring_atom_names;
	    fill_geometry(stacking, res_ref);
	    stackings.push_back(stacking);
	 }
      }
   }
}


#ifdef MAKE_ENHANCED_LIGAND_TOOLS
std::vector<std::vector<std::string> >
pli::pi_stacking_container_t::get_aromatic_ring_list(const RDKit::ROMol &mol_in) const {

   // get a list of aromatic bonds, so that they can be used to find
   // aromatic rings.  To do that, we need to set the aromaticity of
   // the molecule, so make a copy.
   //
   RDKit::RWMol mol(mol_in);
   RDKit::MolOps::cleanUp(mol);
   RDKit::MolOps::setAromaticity(mol);

   // coot::debug_rdkit_molecule(&mol);
   
   std::vector<std::vector<std::string> > ring_list;
   RDKit::RingInfo *ri = mol.getRingInfo(); // don't delete 
   if (! ri->isInitialized()) {
   } else {
      const std::vector<std::vector<int> > &bond_rings = ri->bondRings();
      for (unsigned int ibr=0; ibr<bond_rings.size(); ibr++) {
	 bool all_aromatic = true;
	 for (unsigned int ib=0; ib<bond_rings[ibr].size(); ib++) {
	    const RDKit::Bond *bond_p = mol.getBondWithIdx(bond_rings[ibr][ib]);
	    if (! bond_p->getIsAromatic()) {
	       all_aromatic = false;
	    }
	    
	    int idx_1 = bond_p->getBeginAtomIdx();
	    int idx_2 = bond_p->getEndAtomIdx();
	    std::string at_name_1;
	    std::string at_name_2;

	    try { 
	       mol[idx_1]->getProp("name", at_name_1);
	       mol[idx_2]->getProp("name", at_name_2);
	    }
	    catch (const KeyErrorException &kee) {
	       std::cout << "kee " << kee.what() << std::endl;
	       
	    } 
	 }
	 // std::cout << std::endl;

	 
	 if (all_aromatic) {
	    std::vector<std::string> ring_atom_names;
	    for (unsigned int ib=0; ib<bond_rings[ibr].size(); ib++) {
	       const RDKit::Bond *bond_p = mol.getBondWithIdx(bond_rings[ibr][ib]);
	       int idx_1 = bond_p->getBeginAtomIdx();
	       int idx_2 = bond_p->getEndAtomIdx();
	       std::string at_name_1;
	       std::string at_name_2;
	       try { 
		  mol[idx_1]->getProp("name", at_name_1);
		  mol[idx_2]->getProp("name", at_name_2);
	       }
	       catch (const KeyErrorException &kee) {
		  // std::cout << "kee " << kee.what() << std::endl;
	       }
	       if (!at_name_1.empty() && !at_name_2.empty()) {
		  if (std::find(ring_atom_names.begin(),
				ring_atom_names.end(),
				at_name_1) == ring_atom_names.end())
		     ring_atom_names.push_back(at_name_1);
		  if (std::find(ring_atom_names.begin(),
				ring_atom_names.end(),
				at_name_2) == ring_atom_names.end())
		     ring_atom_names.push_back(at_name_2);
	       }
	    }
	    if (ring_atom_names.size() > 0)
	       ring_list.push_back(ring_atom_names);
	 } 
      }
   }
   return ring_list;
}
#endif // MAKE_ENHANCED_LIGAND_TOOLS


#ifdef MAKE_ENHANCED_LIGAND_TOOLS
std::vector<std::vector<std::string> >
pli::pi_stacking_container_t::get_aromatic_ring_list(const coot::dictionary_residue_restraints_t &monomer_restraints,
						      const RDKit::ROMol &rdkit_mol) const { 

   // Get a list of aromatic bonds, so that they can be used to find
   // aromatic rings.
   // 
   std::vector<std::vector<std::string> > ring_list;
   if (monomer_restraints.ligand_has_aromatic_bonds_p())
      ring_list = monomer_restraints.get_ligand_aromatic_ring_list();
   else
      // it was (or may have been) a kekulized dictionary.  Use rdkit
      // to get the aromatic ring list
      ring_list = get_aromatic_ring_list(rdkit_mol);
   return ring_list;
}
#endif


std::vector<std::vector<std::string> >
pli::pi_stacking_container_t::get_aromatic_ring_list(const coot::dictionary_residue_restraints_t &monomer_restraints) const {

   // get a list of aromatic bonds, so that they can be used to find
   // aromatic rings.
   // 
   std::vector<std::vector<std::string> > ring_list;
   if (monomer_restraints.ligand_has_aromatic_bonds_p())
      ring_list = monomer_restraints.get_ligand_aromatic_ring_list();
   return ring_list;
}

// by search through res_ref
std::vector<std::pair<std::string, clipper::Coord_orth> >
pli::pi_stacking_container_t::get_ligand_cations(mmdb::Residue *res_ref,
						 const coot::dictionary_residue_restraints_t &monomer_restraints) const {

   std::vector<std::pair<std::string, clipper::Coord_orth> > v;
   int n_residue_atoms;
   mmdb::PPAtom residue_atoms = NULL;
   res_ref->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) { 
      std::string ele(residue_atoms[iat]->element);
      if (ele == " N") {
	 // how many bonds does this N have?
	 int n_bonds = 0;
	 std::string atom_name(residue_atoms[iat]->name);
	 for (unsigned int ibond=0; ibond<monomer_restraints.bond_restraint.size(); ibond++) {
	    const coot::dict_bond_restraint_t &br = monomer_restraints.bond_restraint[ibond];
	    if (br.atom_id_1_4c() == atom_name) { 
	       std::string other_atom_name = br.atom_id_2_4c();
	       mmdb::Atom *at = res_ref->GetAtom(other_atom_name.c_str());
	       if (at) { 
		  if (br.type() == "single")
		     n_bonds++;
		  if (br.type() == "double")
		     n_bonds += 2;
		  if (br.type() == "triple")
		     n_bonds += 3;
	       }
	    }

	    if (br.atom_id_2_4c() == atom_name) {
	       std::string other_atom_name = br.atom_id_1_4c();
	       mmdb::Atom *at = res_ref->GetAtom(other_atom_name.c_str());
	       if (at) { 
		  if (br.type() == "single")
		     n_bonds++;
		  if (br.type() == "double")
		     n_bonds += 2;
		  if (br.type() == "triple")
		     n_bonds += 3;
	       }
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
   if (0) { // debug
      std::cout << "debug:: in coot::pi_stacking_container_t::get_ligand_cations() found " << v.size()
		<< " ligand cation points " << std::endl;
      for (unsigned int i=0; i<v.size(); i++)
	 std::cout << "   :" << v[i].first << ":  " << v[i].second.format() << std::endl;
   }
   return v;
}


// can throw an exception because it calls get_ring_pi_centre_points()
//
// should return the stacking type, e.g. PI_CATION_STACKING.
// 
pli::pi_stacking_container_t::overlap_result_t
pli::pi_stacking_container_t::get_pi_overlap_to_ligand_ring(mmdb::Residue *res,
							     const clipper::Coord_orth &ligand_pi_point) const {

   float pi_pi_score = 0;
   float pi_cation_score = 0;
   std::string res_name(res->GetResName());
   std::vector<std::string> best_ring_atom_names;
   std::string best_cation_atom_name;

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
      float score = (score_2 > score_1) ? score_2 : score_1;
      if (score > pi_pi_score) {
	 pi_pi_score = score;
	 best_ring_atom_names = atom_names[iring];
      }
   }

   // Now test all the cation-pi interactions (the score is the sum over the
   // residue's cation atoms; we report the atom that contributed most)
   // 
   std::vector<std::pair<std::string, clipper::Coord_orth> > cation_atom_point = get_cation_atom_positions(res);
   float best_cation_contribution = 0.0;
   for (unsigned int icat=0; icat<cation_atom_point.size(); icat++) {
      float c = overlap_of_cation_pi(ligand_pi_point, cation_atom_point[icat].second);
      pi_cation_score += c;
      if (c > best_cation_contribution) {
	 best_cation_contribution = c;
	 best_cation_atom_name = cation_atom_point[icat].first;
      }
   }

   overlap_result_t r;
   r.score = pi_pi_score;
   r.type = pi_stacking_instance_t::PI_PI_STACKING;
   r.residue_ring_atom_names = best_ring_atom_names;
   if (pi_cation_score > pi_pi_score) {
      r.score = pi_cation_score;
      r.type = pi_stacking_instance_t::PI_CATION_STACKING;
      r.residue_ring_atom_names.clear();
      r.residue_cation_atom_name = best_cation_atom_name;
   }
   return r;
}

// Return the best score of the ligand cation overlap to any of the
// ring systems in res (if any of course, typically this function just
// falls through returning 0.0).
// 
pli::pi_stacking_container_t::overlap_result_t
pli::pi_stacking_container_t::get_pi_overlap_to_ligand_cation(mmdb::Residue *res,
                                                              const clipper::Coord_orth &pt) const {

   overlap_result_t r;
   r.type = pi_stacking_instance_t::CATION_PI_STACKING;
   std::string res_name(res->GetResName());
   std::vector<std::vector<std::string> > atom_names = ring_atom_names(res_name);
   for (unsigned int iring=0; iring<atom_names.size(); iring++) {
      std::pair<clipper::Coord_orth, clipper::Coord_orth> residue_pi_points =
	 get_ring_pi_centre_points(atom_names[iring], res);
      float pi_cation_score_1 = overlap_of_cation_pi(pt, residue_pi_points.first);
      float pi_cation_score_2 = overlap_of_cation_pi(pt, residue_pi_points.second );
      float score = (pi_cation_score_2 > pi_cation_score_1) ? pi_cation_score_2 : pi_cation_score_1;
      if (score > r.score) {
	 r.score = score;
	 r.residue_ring_atom_names = atom_names[iring];
      }
   }
   return r;
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
pli::pi_stacking_container_t::overlap_of_pi_spheres(const clipper::Coord_orth &pt_1,
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
pli::pi_stacking_container_t::get_ring_pi_centre_points(const std::vector<std::string> &ring_atom_names,
							 mmdb::Residue *res_ref) const {
   // dummy points, overwritten.
   clipper::Coord_orth pt_1(0,0,0);
   clipper::Coord_orth pt_2(0,0,0);

   int n_residue_atoms;
   mmdb::PPAtom residue_atoms = NULL;
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
      mess += spec.chain_id; 
      mess += " "; 
      mess += spec.res_no; 
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
pli::pi_stacking_container_t::ring_centre_and_normal(const std::vector<clipper::Coord_orth> &pts) const {

   clipper::Coord_orth centre(0,0,0);
   clipper::Coord_orth normal(0,0,0);

   if (pts.size() < 3) {
      std::string mess = "in ring_centre_and_normal() not enough point to calculate lsq plane";
      throw std::runtime_error(mess);
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

   double norm = sqrt(sum_sq);
   normal = clipper::Coord_orth(eigen_vec.x()/norm,
				eigen_vec.y()/norm,
				eigen_vec.z()/norm);

   return std::pair<clipper::Coord_orth, clipper::Coord_orth> (centre, normal);
}

// TRP (for example) has 2 rings, so we have to return a vector
//
// Don't forget RNA and DNA! FIXME
//
std::vector<std::vector<std::string> >
pli::pi_stacking_container_t::ring_atom_names(const std::string &residue_name) const {

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
   if (residue_name == "HIS") {
      std::vector<std::string> v;
      v.push_back(" CG ");
      v.push_back(" ND1");
      v.push_back(" CD2");
      v.push_back(" CE1");
      v.push_back(" NE2");
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
       (residue_name == "A") ||
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
       (residue_name == "G") ||
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
       (residue_name == "C") ||
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
       (residue_name == "T") ||
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

   if ((residue_name == "DU") ||   // or RNA equivalent
       (residue_name == "U") ||
       (residue_name == "Ud") ||
       (residue_name == "Ur")) { 
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
       (residue_name == "T") ||
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


std::vector<std::pair<std::string, clipper::Coord_orth> >
pli::pi_stacking_container_t::get_cation_atom_positions(mmdb::Residue *res) const {

   std::vector<std::pair<std::string, clipper::Coord_orth> > v;

   std::string res_name(res->GetResName());

   if (res_name == "LYS") {
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      res->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int i=0; i<n_residue_atoms; i++) {
	 std::string atom_name(residue_atoms[i]->name);
	 if (atom_name == " NZ ") {
	    clipper::Coord_orth pt(residue_atoms[i]->x,
				   residue_atoms[i]->y,
				   residue_atoms[i]->z);
	    v.push_back(std::make_pair(atom_name, pt));
	 }
      }
   }

   if (res_name == "ARG") {
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      res->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int i=0; i<n_residue_atoms; i++) {
	 std::string atom_name(residue_atoms[i]->name);
	 if ((atom_name == " NH1") ||
	     (atom_name == " NH2")) {
	    clipper::Coord_orth pt(residue_atoms[i]->x,
				   residue_atoms[i]->y,
				   residue_atoms[i]->z);
	    v.push_back(std::make_pair(atom_name, pt));
	 }
      }
   }
   

   return v;
}

float
pli::pi_stacking_container_t::overlap_of_cation_pi(const clipper::Coord_orth &ligand_pi_point,
						    const clipper::Coord_orth &cation_atom_point) const {

   // multipliers follow the order of the atom/control points
   //
   float score = overlap_of_pi_spheres(ligand_pi_point, cation_atom_point, 0.78, -1, 7.8, -0.05);

   // std::cout << "overlap_of_cation_pi() calling overlap_of_pi_spheres() returns "
   // << score << std::endl;

   return score;
}


std::vector<clipper::Coord_orth>
pli::pi_stacking_container_t::ring_atom_positions(const std::vector<std::string> &names,
                                                  mmdb::Residue *res) const {
   std::vector<clipper::Coord_orth> pts;
   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms = 0;
   res->GetAtomTable(residue_atoms, n_residue_atoms);
   for (unsigned int i=0; i<names.size(); i++)
      for (int iat=0; iat<n_residue_atoms; iat++)
         if (names[i] == std::string(residue_atoms[iat]->name)) {
            pts.push_back(clipper::Coord_orth(residue_atoms[iat]->x, residue_atoms[iat]->y, residue_atoms[iat]->z));
            break;
         }
   return pts;
}

std::vector<std::string>
pli::pi_stacking_container_t::edge_to_face_ring_by_geometry(const std::vector<std::string> &ligand_ring_atom_names,
                                                            mmdb::Residue *res_ref, mmdb::Residue *res) const {

   const double centroid_dist_max = 6.5;  // ProLIF EdgeToFace defaults
   const double plane_angle_min   = 50.0;
   const double n2c_angle_max     = 30.0;
   const double intersect_radius  = 1.5;

   std::vector<std::string> result;
   std::vector<clipper::Coord_orth> lig_pts = ring_atom_positions(ligand_ring_atom_names, res_ref);
   if (lig_pts.size() < 3) return result;

   auto folded_angle = [] (const clipper::Coord_orth &a, const clipper::Coord_orth &b) {
      double c = clipper::Coord_orth::dot(a.unit(), b.unit());
      if (c < 0.0) c = -c;
      if (c > 1.0) c = 1.0;
      return clipper::Util::rad2d(acos(c));
   };

   try {
      std::pair<clipper::Coord_orth, clipper::Coord_orth> lig = ring_centre_and_normal(lig_pts);
      std::vector<std::vector<std::string> > res_rings = ring_atom_names(std::string(res->GetResName()));
      double best_d = centroid_dist_max + 1.0;
      for (unsigned int iring=0; iring<res_rings.size(); iring++) {
	 std::vector<clipper::Coord_orth> res_pts = ring_atom_positions(res_rings[iring], res);
	 if (res_pts.size() < 3) continue;
	 std::pair<clipper::Coord_orth, clipper::Coord_orth> rr = ring_centre_and_normal(res_pts);
	 clipper::Coord_orth cc = rr.first - lig.first;
	 double d = sqrt(cc.lengthsq());
	 if (d > centroid_dist_max) continue;
	 if (folded_angle(lig.second, rr.second) < plane_angle_min) continue;
	 double a1 = folded_angle(lig.second, cc);
	 double a2 = folded_angle(rr.second, cc);
	 if (a1 > n2c_angle_max && a2 > n2c_angle_max) continue;

	 // the line where the two ring planes meet: direction n1 x n2, through
	 // a point satisfying n1.p = n1.c1, n2.p = n2.c2, (n1 x n2).p = 0
	 const clipper::Coord_orth &n1 = lig.second;
	 const clipper::Coord_orth &n2 = rr.second;
	 clipper::Coord_orth dir(clipper::Coord_orth::cross(n1, n2));
	 clipper::Mat33<double> A(n1.x(),  n1.y(),  n1.z(),
				  n2.x(),  n2.y(),  n2.z(),
				  dir.x(), dir.y(), dir.z());
	 if (fabs(A.det()) < 1e-6) continue; // parallel planes
	 clipper::Vec3<double> rhs(clipper::Coord_orth::dot(n1, lig.first),
				   clipper::Coord_orth::dot(n2, rr.first),
				   0.0);
	 clipper::Vec3<double> pv = A.inverse() * rhs;
	 clipper::Coord_orth p(pv[0], pv[1], pv[2]);
	 // projection of the ligand centroid onto that line
	 clipper::Coord_orth u(dir.unit());
	 clipper::Coord_orth q = p + clipper::Coord_orth(clipper::Coord_orth::dot(lig.first - p, u) * u);
	 double dq_lig = sqrt((q - lig.first).lengthsq());
	 double dq_res = sqrt((q - rr.first).lengthsq());
	 if (dq_lig > intersect_radius && dq_res > intersect_radius) continue;

	 if (d < best_d) {
	    best_d = d;
	    result = res_rings[iring];
	 }
      }
   }
   catch (const std::runtime_error &rte) {
      result.clear();
   }
   return result;
}

// Geometry for comparison with other tools (ProLIF, PLIP): centroid distance,
// plane angle and normal-to-centroid angle. Angles are folded into 0..90.
//
void
pli::pi_stacking_container_t::fill_geometry(pi_stacking_instance_t &st, mmdb::Residue *res_ref) const {

   auto ring_points = [this] (const std::vector<std::string> &names, mmdb::Residue *res) {
      return ring_atom_positions(names, res);
   };
   auto atom_position = [] (const std::string &name, mmdb::Residue *res, clipper::Coord_orth &pos) {
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms = 0;
      res->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++)
         if (name == std::string(residue_atoms[iat]->name)) {
            pos = clipper::Coord_orth(residue_atoms[iat]->x, residue_atoms[iat]->y, residue_atoms[iat]->z);
            return true;
         }
      return false;
   };
   auto folded_angle = [] (const clipper::Coord_orth &a, const clipper::Coord_orth &b) {
      double c = clipper::Coord_orth::dot(a.unit(), b.unit());
      if (c < 0.0) c = -c;
      if (c > 1.0) c = 1.0;
      return clipper::Util::rad2d(acos(c));
   };

   try {
      if (st.type == pi_stacking_instance_t::PI_PI_STACKING) {
         std::vector<clipper::Coord_orth> lig_pts = ring_points(st.ligand_ring_atom_names, res_ref);
         std::vector<clipper::Coord_orth> res_pts = ring_points(st.residue_ring_atom_names, st.res);
         if (lig_pts.size() < 3 || res_pts.size() < 3) return;
         std::pair<clipper::Coord_orth, clipper::Coord_orth> lig = ring_centre_and_normal(lig_pts);
         std::pair<clipper::Coord_orth, clipper::Coord_orth> res = ring_centre_and_normal(res_pts);
         clipper::Coord_orth cc = res.first - lig.first;
         st.centroid_distance = sqrt(cc.lengthsq());
         st.plane_angle = folded_angle(lig.second, res.second);
         double a1 = folded_angle(lig.second, cc);
         double a2 = folded_angle(res.second, cc);
         st.normal_to_centroid_angle = (a2 < a1) ? a2 : a1;
      }
      if (st.type == pi_stacking_instance_t::PI_CATION_STACKING) {
         // ligand ring, residue cation atom
         std::vector<clipper::Coord_orth> lig_pts = ring_points(st.ligand_ring_atom_names, res_ref);
         clipper::Coord_orth cat;
         if (lig_pts.size() < 3) return;
         if (! atom_position(st.residue_cation_atom_name, st.res, cat)) return;
         std::pair<clipper::Coord_orth, clipper::Coord_orth> lig = ring_centre_and_normal(lig_pts);
         clipper::Coord_orth cc = cat - lig.first;
         st.centroid_distance = sqrt(cc.lengthsq());
         st.normal_to_centroid_angle = folded_angle(lig.second, cc);
      }
      if (st.type == pi_stacking_instance_t::CATION_PI_STACKING) {
         // ligand cation atom, residue ring
         std::vector<clipper::Coord_orth> res_pts = ring_points(st.residue_ring_atom_names, st.res);
         clipper::Coord_orth cat;
         if (res_pts.size() < 3) return;
         if (! atom_position(st.ligand_cationic_atom_name, res_ref, cat)) return;
         std::pair<clipper::Coord_orth, clipper::Coord_orth> res = ring_centre_and_normal(res_pts);
         clipper::Coord_orth cc = cat - res.first;
         st.centroid_distance = sqrt(cc.lengthsq());
         st.normal_to_centroid_angle = folded_angle(res.second, cc);
      }
   }
   catch (const std::runtime_error &rte) {
      // leave the geometry unset
   }
}

std::vector<pli::pi_stacking_instance_t>
pli::get_pi_stackings(mmdb::Residue *ligand_residue, mmdb::Manager *mol,
                      const coot::protein_geometry &geom, int imol,
                      float residues_near_radius) {

   std::vector<pi_stacking_instance_t> v;
   if (! ligand_residue) return v;
   std::string res_name(ligand_residue->GetResName());
   std::pair<bool, coot::dictionary_residue_restraints_t> rp = geom.get_monomer_restraints(res_name, imol);
   if (! rp.first) {
      std::cout << "WARNING:: get_pi_stackings(): no dictionary for " << res_name << std::endl;
      return v;
   }
   std::vector<mmdb::Residue *> residues = coot::residues_near_residue(ligand_residue, mol, residues_near_radius);

   if (rp.second.ligand_has_aromatic_bonds_p()) {
      pi_stacking_container_t psc(rp.second, residues, ligand_residue);
      v = psc.stackings;
   } else {
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
      // a kekulized dictionary: let RDKit find the aromatic rings
      try {
         RDKit::RWMol rdkm = coot::rdkit_mol_sanitized(ligand_residue, imol, geom);
         pi_stacking_container_t psc(rp.second, residues, ligand_residue, rdkm);
         v = psc.stackings;
      }
      catch (const std::exception &e) {
         std::cout << "WARNING:: get_pi_stackings(): RDKit molecule failed for " << res_name
                   << ": " << e.what() << std::endl;
      }
#else
      pi_stacking_container_t psc(rp.second, residues, ligand_residue); // no rings found
      v = psc.stackings;
#endif
   }
   return v;
}

std::ostream&
pli::operator<<(std::ostream& s, const pi_stacking_instance_t &stack) {

   std::string st_type = "UNKNOWN";
   if (stack.type == pi_stacking_instance_t::NO_STACKING)        st_type = "NO_STACKING";
   if (stack.type == pi_stacking_instance_t::PI_PI_STACKING)     st_type = "PI_PI_STACKING";
   if (stack.type == pi_stacking_instance_t::PI_CATION_STACKING) st_type = "PI_CATION_STACKING";
   if (stack.type == pi_stacking_instance_t::CATION_PI_STACKING) st_type = "CATION_PI_STACKING";
   
   s << "[" << st_type << " " << coot::residue_spec_t(stack.res) << " "
     << stack.overlap_score << " ligand-atom-name :"
     <<  stack.ligand_cationic_atom_name
     << ": ";
   for (unsigned int i=0; i<stack.ligand_ring_atom_names.size(); i++)
      s << "  :" << stack.ligand_ring_atom_names[i] << ":   " ;
   s << "]";
   return s;
}
