
#include <map>
#include <algorithm>

// these are needed to compile molecule-compile-info.h:
// 
#include "Cartesian.h"
#include "mmdb-extras.h"
#include "mmdb-crystal.h"

// morphing
#include "coot-coord-utils.hh"
#include "rigid-body.hh"

#include "molecule-class-info.h"

// return an index of the new restraint
int
molecule_class_info_t::add_extra_bond_restraint(coot::atom_spec_t atom_1,
						coot::atom_spec_t atom_2,
						double bond_dist, double esd) {
   int r = -1; // unset
   CAtom *at_1 = get_atom(atom_1);
   CAtom *at_2 = get_atom(atom_2);
   if (at_1) {
      int atom_index = -1;
      at_1->GetUDData(atom_sel.UDDAtomIndexHandle, atom_index); // set atom_index
      atom_1.int_user_data = atom_index;
   }
   if (at_2) {
      int atom_index = -1;
      at_2->GetUDData(atom_sel.UDDAtomIndexHandle, atom_index); // set atom_index
      atom_2.int_user_data = atom_index;
   }
   if (at_1 && at_2) { 
      coot::extra_restraints_t::extra_bond_restraint_t bond(atom_1, atom_2, bond_dist, esd);
      extra_restraints.bond_restraints.push_back(bond);
      update_extra_restraints_representation();
      r = extra_restraints.bond_restraints.size() -1;
   }
   return r;
}


void molecule_class_info_t::remove_extra_bond_restraint(coot::atom_spec_t atom_1, coot::atom_spec_t atom_2) {

   int n_bonds_pre = extra_restraints.bond_restraints.size();
   std::vector<coot::extra_restraints_t::extra_bond_restraint_t>::iterator it;
   for (it=extra_restraints.bond_restraints.begin(); it != extra_restraints.bond_restraints.end(); it++) { 
      if (((it->atom_1 == atom_1) &&
	   (it->atom_2 == atom_2)) ||
	  ((it->atom_2 == atom_1) &&
	   (it->atom_1 == atom_2))) {
	 extra_restraints.bond_restraints.erase(it);
	 std::cout << "deleted extra bond restraint " << atom_1 << " to " << atom_2 << std::endl;
	 // break;
      }
   }
   int n_bonds_post = extra_restraints.bond_restraints.size();
   std::cout << "DEBUG:: pre: " << n_bonds_pre << " post " << n_bonds_post << std::endl;
   update_extra_restraints_representation();
}

void
molecule_class_info_t::add_refmac_extra_restraints(const std::string &file_name) {

   coot::extra_restraints_t r;
   r.read_refmac_distance_restraints(file_name);
   extra_restraints.add_restraints(r);
   std::cout << "in add_refmac_extra_restraints we have " << extra_restraints.bond_restraints.size()
	     << " bond restraints " << std::endl;
   update_extra_restraints_representation();
}

void
molecule_class_info_t::delete_extra_restraints_for_residue(const coot::residue_spec_t &rs) {

   unsigned int pre_n = extra_restraints.bond_restraints.size(); 
   extra_restraints.delete_restraints_for_residue(rs);
   unsigned int post_n = extra_restraints.bond_restraints.size();
   if (post_n != pre_n)
      update_extra_restraints_representation();
}

void
molecule_class_info_t::delete_extra_restraints_worse_than(const double &n_sigma) {

   unsigned int pre_n = extra_restraints.bond_restraints.size();

   // the real dist, with atom specs for keys
   // (to be used by the erasor)
   // 
   std::map<std::pair<coot::atom_spec_t, coot::atom_spec_t>, double> dist_map;
   std::map<coot::atom_spec_t, CAtom *> atom_map;
   std::map<coot::atom_spec_t, CAtom *>::const_iterator it_1;
   std::map<coot::atom_spec_t, CAtom *>::const_iterator it_2;
   // first fill the dist_map and fill the atom_map as you do so.
   for (unsigned int i=0; i<extra_restraints.bond_restraints.size(); i++) {
      coot::extra_restraints_t::extra_bond_restraint_t &br = extra_restraints.bond_restraints[i];
      CAtom *at_1 = NULL;
      CAtom *at_2 = NULL;
      it_1 = atom_map.find(br.atom_1);
      it_2 = atom_map.find(br.atom_2);
      if (it_1 == atom_map.end()) {
	 at_1 = get_atom(br.atom_1);
	 atom_map[br.atom_1] = at_1;
      } else {
	 at_1 = it_1->second; // most of the hits, I hope
      }
      // and the same for 2:
      if (it_2 == atom_map.end()) {
	 at_2 = get_atom(br.atom_2);
	 atom_map[br.atom_2] = at_2;
      } else {
	 at_2 = it_2->second; 
      }
      if (at_1 && at_2) {
	 double dx = at_1->x - at_2->x;
	 double dy = at_1->y - at_2->y;
	 double dz = at_1->z - at_2->z;
	 double dist_sq = dx*dx + dy*dy + dz*dz;
	 if (dist_sq < 0) dist_sq = 0;
	 dist_map[std::pair<coot::atom_spec_t, coot::atom_spec_t>(br.atom_1, br.atom_2)] = sqrt(dist_sq);
      }
   }

   extra_restraints.bond_restraints.erase(std::remove_if(extra_restraints.bond_restraints.begin(), extra_restraints.bond_restraints.end(), coot::extra_restraints_t::bond_erasor(dist_map, n_sigma)), extra_restraints.bond_restraints.end());

   unsigned int post_n = extra_restraints.bond_restraints.size();
   if (post_n != pre_n)
      update_extra_restraints_representation();
   std::cout << "INFO deleted : " << pre_n - post_n << " extra bond restraints" << std::endl;
} 


void
molecule_class_info_t::clear_extra_restraints() {
   extra_restraints.clear();
} 

// return an index of the new restraint
int
molecule_class_info_t::add_extra_torsion_restraint(coot::atom_spec_t atom_1,
						   coot::atom_spec_t atom_2,
						   coot::atom_spec_t atom_3,
						   coot::atom_spec_t atom_4,
						   double torsion_angle, double esd, int period) {

   CAtom *at_1 = get_atom(atom_1);
   CAtom *at_2 = get_atom(atom_2);
   CAtom *at_3 = get_atom(atom_3);
   CAtom *at_4 = get_atom(atom_4);
   if (at_1) {
      int atom_index = -1;
      at_1->GetUDData(atom_sel.UDDAtomIndexHandle, atom_index); // set atom_index
      atom_1.int_user_data = atom_index;
   }
   if (at_2) {
      int atom_index = -1;
      at_2->GetUDData(atom_sel.UDDAtomIndexHandle, atom_index); // set atom_index
      atom_2.int_user_data = atom_index;
   }
   if (at_3) {
      int atom_index = -1;
      at_3->GetUDData(atom_sel.UDDAtomIndexHandle, atom_index); // set atom_index
      atom_3.int_user_data = atom_index;
   }
   if (at_4) {
      int atom_index = -1;
      at_4->GetUDData(atom_sel.UDDAtomIndexHandle, atom_index); // set atom_index
      atom_4.int_user_data = atom_index;
   }
   coot::extra_restraints_t::extra_torsion_restraint_t tors(atom_1, atom_2,
							    atom_3,atom_4,
							    torsion_angle, esd, period);
   extra_restraints.torsion_restraints.push_back(tors);
   update_extra_restraints_representation();
   return extra_restraints.torsion_restraints.size() -1;
}


int
molecule_class_info_t::morph_fit_all(const clipper::Xmap<float> &xmap_in, float transformation_average_radius) {

   int success = 0;

   // construct minimol fragments 

   int n_neighb=2; // either side of central residue
   int imod = 1;
   CModel *model_p = atom_sel.mol->GetModel(imod);
   CChain *chain_p;
   int n_chains = model_p->GetNumberOfChains();
   // store the local origin too.
   std::map<CResidue *, morph_rtop_triple> rtop_map;
   for (int ichain=0; ichain<n_chains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();

      for (int ires=n_neighb; ires<(nres-n_neighb); ires++) { // residue-in-chain loop

	 std::vector<CResidue *> fragment_residues;
	 CResidue *residue_p = chain_p->GetResidue(ires);
	 std::pair<bool, clipper::Coord_orth> local_centre = residue_centre(residue_p);
	 if (local_centre.first) { 
	    std::cout << "\rINFO:: Getting RTops for " << coot::residue_spec_t(residue_p);
	    std::cout.flush();
	    for (int ifragres=-n_neighb; ifragres<=n_neighb; ifragres++) {
	       CResidue *r = chain_p->GetResidue(ires+ifragres);
	       if (!r)
		  break;
	       fragment_residues.push_back(r);
	    }

	    coot::minimol::fragment f;
	    for (unsigned int ifr=0; ifr<fragment_residues.size(); ifr++) {
	       f.addresidue(coot::minimol::residue(fragment_residues[ifr]), false);
	    } 
	    coot::minimol::molecule m(f);
	    // m.write_file("working-fragment.pdb",10);
	    // returns the rtop to move m into map.
	    std::pair<bool, clipper::RTop_orth> rtop = coot::get_rigid_body_fit_rtop(&m, local_centre.second, xmap_in);
	    morph_rtop_triple rt(local_centre.second, rtop);
	    rtop_map[residue_p] = rt; 
	 }
      }
      std::cout << std::endl;
   }

   std::map<CResidue *, morph_rtop_triple>::const_iterator it;
   if (rtop_map.size()) {
      success = 1;
      make_backup();

      std::map<CResidue *, clipper::RTop_orth> simple_shifts;
      std::map<CResidue *, clipper::RTop_orth> smooth_shifts;
      
      for (it=rtop_map.begin(); it!=rtop_map.end(); it++) {
	 CResidue *this_residue = it->first;
	 if (it->second.valid) {

	    // Morphing step is super-fast
	    // std::cout << "\rINFO:: Morphing " << coot::residue_spec_t(this_residue);
	    // std::cout.flush();
	    
	    std::vector<std::pair<clipper::RTop_orth, float> > rtops;
	    // std::cout << "this residue:\n" << it->second.second.format() << std::endl;
	    rtops.push_back(std::pair<clipper::RTop_orth,float>(it->second.rtop, 1));
 	    std::vector<CResidue *> neighb_residues =
 	       coot::residues_near_residue(this_residue, atom_sel.mol, transformation_average_radius);

  	    for (unsigned int i_n_res=0; i_n_res<neighb_residues.size(); i_n_res++) { 
	       std::map<CResidue *, morph_rtop_triple>::const_iterator it_for_neighb =
		  rtop_map.find(neighb_residues[i_n_res]);
	       if (it_for_neighb != rtop_map.end()) { 
		  if (it_for_neighb->second.valid) {
		     std::pair<clipper::RTop_orth, float> p(it_for_neighb->second.rtop, 1.0);
		     rtops.push_back(p);
		     if (0)
			std::cout << "adding rtop for " << coot::residue_spec_t(neighb_residues[i_n_res]) << "\n"
				  << it_for_neighb->second.rtop.format() << std::endl;
		  }
	       }
	    }

	    // pre-local shifts and quaternion-based rtop averaging
	    // morph_residue_atoms_by_average_rtops(this_residue, rtops);

	    coot::util::quaternion q(0,0,0,0);
	    clipper::RTop_orth rtop = q.centroid_rtop(rtops);
	    // if (0)
	    // std::cout << "yields averaged RTop:\n" << rtop.format() << std::endl;
	    
	    transform_by_internal(rtop, it->first);
	    simple_shifts[this_residue] = it->second.rtop; // save just to view them
	    smooth_shifts[this_residue] = rtop;
	    
	    
	 } else {
	    std::cout << "no RTop for " << coot::residue_spec_t(it->first) << std::endl;
	 } 
      }
      std::cout << std::endl;
      atom_sel.mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked();

      // morph_show_shifts(simple_shifts, smooth_shifts);
      
   }
   return success;
}

// I fail to make a function that does a good "average" of RTops,
// so do it long-hand by generating sets of coordinates by applying
// each rtop to each atom - weights are transfered in the second part of the pair.
//
// This doesn't do backups or unsaved changes marking of course.
void
molecule_class_info_t::morph_residue_atoms_by_average_rtops(CResidue *residue_p,
							    const std::vector<std::pair<clipper::RTop_orth, float> > &rtops) {

   PPCAtom residue_atoms = 0;
   int n_residue_atoms;
   if (rtops.size()) { 
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (unsigned int iat=0; iat<n_residue_atoms; iat++) {
	 CAtom *at = residue_atoms[iat];
	 clipper::Coord_orth pt(at->x, at->y, at->z);
	 clipper::Coord_orth sum_transformed_pts(0,0,0);
	 double sum_weights = 0.0;
	 for (unsigned int i_rtop=0; i_rtop<rtops.size(); i_rtop++) { 
	    clipper::Coord_orth t_pt = pt.transform(rtops[i_rtop].first);
	    double weight = rtops[i_rtop].second;
	    sum_weights += weight;
	    sum_transformed_pts += t_pt*weight;
	 }
	 if (sum_weights > 0.0) {
	    double inv_weight = 1.0/sum_weights;
	    clipper::Coord_orth new_pt(sum_transformed_pts.x() * inv_weight,
				       sum_transformed_pts.y() * inv_weight,
				       sum_transformed_pts.z() * inv_weight);
	    at->x = new_pt.x();
	    at->y = new_pt.y();
	    at->z = new_pt.z();
	 }
      }
   }
}



void
molecule_class_info_t::morph_show_shifts(const std::map<CResidue *, clipper::RTop_orth> &simple_shifts,
					 const std::map<CResidue *, clipper::RTop_orth> &smooth_shifts) const {

   // write a file
   
   std::map<CResidue *, clipper::RTop_orth>::const_iterator it;
   std::ofstream f("morph-shifts.scm");

   std::string ss;
   ss = "(define simple-shifts (new-generic-object-number \"simple-shifts\"))";
   f << ss << "\n";
   ss = "(define smooth-shifts (new-generic-object-number \"smooth-shifts\"))";
   f << ss << "\n";
   ss = "(set-display-generic-object simple-shifts 1)";
   f << ss << "\n";
   ss = "(set-display-generic-object smooth-shifts 1)";
   f << ss << "\n";

   for (it=simple_shifts.begin(); it!=simple_shifts.end(); it++) {
      CResidue *r = it->first;
      std::pair<bool, clipper::Coord_orth> rc = residue_centre(r);
      CAtom *C_alpha = r->GetAtom(" CA ");
      if (! C_alpha)
	 C_alpha = r->GetAtom(" P  ");
      if (C_alpha) {
	 clipper::Coord_orth ca_pos(C_alpha->x, C_alpha->y, C_alpha->z);
	 if (rc.first) {
	    std::string s;
	    std::string line_colour = "yellow";
	    std::string ball_colour = "yellow";

	    clipper::Coord_orth to_pos = ca_pos.transform(it->second);
	 
	    s += "(to-generic-object-add-line  simple-shifts ";
	    s += "\"";
	    s += line_colour;
	    s += "\"";
	    s += "  2 ";
	    s += coot::util::float_to_string(ca_pos.x());
	    s += " ";
	    s += coot::util::float_to_string(ca_pos.y());
	    s += " ";
	    s += coot::util::float_to_string(ca_pos.z());
	    s += " ";
	    s += coot::util::float_to_string(to_pos.x());
	    s += " ";
	    s += coot::util::float_to_string(to_pos.y());
	    s += " ";
	    s += coot::util::float_to_string(to_pos.z());
	    s += " ";
	    s += ")";
	    f << s << "\n";
	    s = "";
	    s += "(to-generic-object-add-point simple-shifts ";
	    s += "\"";
	    s += ball_colour;
	    s += "\"";
	    s += " 12                   ";
	    s += coot::util::float_to_string(to_pos.x());
	    s += " ";
	    s += coot::util::float_to_string(to_pos.y());
	    s += " ";
	    s += coot::util::float_to_string(to_pos.z());
	    s += " ";
	    s += ")";
	    f << s << "\n";
	 }
      }
   }
   for (it=smooth_shifts.begin(); it!=smooth_shifts.end(); it++) {
      CResidue *r = it->first;
      std::pair<bool, clipper::Coord_orth> rc = residue_centre(r);
      CAtom *C_alpha = r->GetAtom(" CA ");
      if (! C_alpha)
	 C_alpha = r->GetAtom(" P  ");
      if (C_alpha) {
	 clipper::Coord_orth ca_pos(C_alpha->x, C_alpha->y, C_alpha->z);
	 if (rc.first) {
	    std::string s;
	    std::string line_colour = "red";
	    std::string ball_colour = "red";

	    clipper::Coord_orth to_pos = ca_pos.transform(it->second);
	 
	    s += "(to-generic-object-add-line  smooth-shifts ";
	    s += "\"";
	    s += line_colour;
	    s += "\"";
	    s += " 2 ";
	    s += coot::util::float_to_string(ca_pos.x());
	    s += " ";
	    s += coot::util::float_to_string(ca_pos.y());
	    s += " ";
	    s += coot::util::float_to_string(ca_pos.z());
	    s += " ";
	    s += coot::util::float_to_string(to_pos.x());
	    s += " ";
	    s += coot::util::float_to_string(to_pos.y());
	    s += " ";
	    s += coot::util::float_to_string(to_pos.z());
	    s += " ";
	    s += ")";
	    f << s << "\n";
	    s = "";
	    s += "(to-generic-object-add-point smooth-shifts ";
	    s += "\"";
	    s += ball_colour;
	    s += "\"";
	    s += " 12                 ";
	    s += coot::util::float_to_string(to_pos.x());
	    s += " ";
	    s += coot::util::float_to_string(to_pos.y());
	    s += " ";
	    s += coot::util::float_to_string(to_pos.z());
	    s += " ";
	    s += ")";
	    f << s << "\n";
	 }
      }
   }
   f.close();
}

