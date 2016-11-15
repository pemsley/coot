/* coot-utils/atom-overlaps.cc
 * 
 * Copyright 2015 by Medical Research Council
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


#include "atom-overlaps.hh"
#include "coot-coord-utils.hh"
// #include "coot-h-bonds.hh"

coot::atom_overlaps_container_t::atom_overlaps_container_t(mmdb::Residue *res_central_in,
							   const std::vector<mmdb::Residue *> &neighbours_in,
							   mmdb::Manager *mol_in,
							   const protein_geometry *geom_p_in) {
   geom_p = geom_p_in;
   res_central = res_central_in;
   neighbours = neighbours_in;
   mol = mol_in;
   clash_spike_length = 0.5;
   init();

}

// is this used?
coot::atom_overlaps_container_t::atom_overlaps_container_t(mmdb::Residue *res_central_in,
							   mmdb::Residue *neighbour,
							   mmdb::Manager *mol_in,
							   const protein_geometry *geom_p_in) {
   geom_p = geom_p_in;
   res_central = res_central_in;
   neighbours.push_back(neighbour);
   mol = mol_in;
   clash_spike_length = 0.5;
   init();

}

// clash spike_length should be 0.5;
coot::atom_overlaps_container_t::atom_overlaps_container_t(mmdb::Residue *res_central_in,
							   const std::vector<mmdb::Residue *> &neighbours_in,
							   mmdb::Manager *mol_in,
							   const protein_geometry *geom_p_in,
							   double clash_spike_length_in,
							   double probe_radius_in) {
   probe_radius = probe_radius_in;
   geom_p = geom_p_in;
   res_central = res_central_in;
   neighbours = neighbours_in;
   mol = mol_in;
   clash_spike_length = clash_spike_length_in;
   init();

}


void
coot::atom_overlaps_container_t::init() {

   have_dictionary = false; // initially.

   if (res_central) {

      std::string cres_name = res_central->GetResName();
      std::pair<bool, dictionary_residue_restraints_t> d =
	 geom_p->get_monomer_restraints(cres_name, protein_geometry::IMOL_ENC_ANY);
      if (! d.first) {
	 std::cout << "Failed to get dictionary for " << cres_name << std::endl;
      } else {
	 // Happy path
	 central_residue_dictionary = d.second;

	 if (false)
	    std::cout << "central_residue_dictionary has " << central_residue_dictionary.atom_info.size()
		      << " atoms and " << central_residue_dictionary.bond_restraint.size()
		      << " bond restraints " << std::endl;

	 neighb_dictionaries.resize(neighbours.size());
	 have_dictionary = true;
	 for (unsigned int i=0; i<neighbours.size(); i++) {
	    std::string residue_name = neighbours[i]->GetResName();
	    d = geom_p->get_monomer_restraints(residue_name, protein_geometry::IMOL_ENC_ANY);
	    if (! d.first) {
	       std::cout << "WARNING:: Overlap fail. Failed to get dictionary for name "
			 << residue_name << std::endl;
	       have_dictionary = false;
	       break;
	    } else {
	       // happy path
	       neighb_dictionaries[i] = d.second;
	    }
	 }
      }

      if (have_dictionary) {
	 fill_ligand_atom_neighbour_map(); // and add radius
	 mark_donors_and_acceptors();
      }
   }
}

void
coot::atom_overlaps_container_t::mark_donors_and_acceptors() {

   // now mark donors and acceptors.
   //
   udd_h_bond_type_handle = mol->RegisterUDInteger(mmdb::UDR_ATOM, "hb_type");

   // 
   mmdb::PAtom *central_residue_atoms = 0;
   int n_central_residue_atoms;
   res_central->GetAtomTable(central_residue_atoms, n_central_residue_atoms);
   for (int iat=0; iat<n_central_residue_atoms; iat++) { 
      mmdb::Atom *at = central_residue_atoms[iat];
      std::string atom_name(at->name);
      std::string ele = at->element;
      if (ele == " H") {
	 // Hydrogens have energy type "H" from Refmac and acedrg, that doesn't
	 // tell us if this atom is a donor hydrogen.
	 // So, find the atom to which the H is attached and if that is a donor then this
	 // is a hydrogen bond hydrogen.
	 std::string heavy_neighb_of_H_atom =
	    central_residue_dictionary.get_bonded_atom(atom_name);
	 if (! heavy_neighb_of_H_atom.empty()) {
	    std::string neigh_energy_type = central_residue_dictionary.type_energy(heavy_neighb_of_H_atom);
	    energy_lib_atom neighb_ela = geom_p->get_energy_lib_atom(neigh_energy_type);
	    hb_t neighb_hb_type = neighb_ela.hb_type;
	    if (neighb_hb_type == coot::HB_DONOR) {
	       // std::cout << "----- adding ligand HB_HYDROGEN udd " << atom_spec_t(at) << std::endl;
	       at->PutUDData(udd_h_bond_type_handle, coot::HB_HYDROGEN); // hb_t -> int
	    }
	    if (neighb_hb_type == coot::HB_BOTH) {
	       // std::cout << "----- adding ligand HB_HYDROGEN udd " << atom_spec_t(at) << std::endl;
	       at->PutUDData(udd_h_bond_type_handle, coot::HB_HYDROGEN); // hb_t -> int
	    }
	 }
      } else {
	 std::string energy_type = central_residue_dictionary.type_energy(atom_name);
	 energy_lib_atom ela = geom_p->get_energy_lib_atom(energy_type);
	 hb_t hb_type = ela.hb_type;
	 at->PutUDData(udd_h_bond_type_handle, hb_type); // hb_t -> int
      }
   }

   for (unsigned int i=0; i<neighbours.size(); i++) {
      const dictionary_residue_restraints_t &dict = neighb_dictionaries[i];
      mmdb::PAtom *residue_atoms = 0;
      int n_residue_atoms;
      neighbours[i]->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) { 
	 mmdb::Atom *n_at = residue_atoms[iat];
	 std::string atom_name(n_at->name);
	 std::string ele = n_at->element;
	 if (ele == " H") {
	    // as above
	    std::string heavy_neighb_of_H_atom = dict.get_bonded_atom(atom_name);
	    if (! heavy_neighb_of_H_atom.empty()) {
	       std::string neigh_energy_type = dict.type_energy(heavy_neighb_of_H_atom);
	       energy_lib_atom neighb_ela = geom_p->get_energy_lib_atom(neigh_energy_type);
	       hb_t neighb_hb_type = neighb_ela.hb_type;
	       if (neighb_hb_type == coot::HB_DONOR) {
		  // std::cout << "----- adding env HB_HYDROGEN udd " << atom_spec_t(n_at) << std::endl;
		  n_at->PutUDData(udd_h_bond_type_handle, coot::HB_HYDROGEN); // hb_t -> int
	       }
	       if (neighb_hb_type == coot::HB_BOTH) {
		  // std::cout << "----- adding env HB_HYDROGEN udd " << atom_spec_t(n_at) << std::endl;
		  n_at->PutUDData(udd_h_bond_type_handle, coot::HB_HYDROGEN); // hb_t -> int
	       }
	    }
	 } else {
	    std::string atom_name(n_at->name);
	    std::string energy_type = neighb_dictionaries[i].type_energy(atom_name);
	    energy_lib_atom ela = geom_p->get_energy_lib_atom(energy_type);
	    hb_t hb_type = ela.hb_type;
	    n_at->PutUDData(udd_h_bond_type_handle, hb_type); // hb_t -> int
	 }
      }
   }
}

// and radius
void
coot::atom_overlaps_container_t::fill_ligand_atom_neighbour_map() {

   mmdb::realtype max_dist = 2.3;
   if (mol) {
      mmdb::Contact *pscontact = NULL;
      int n_contacts;
      float min_dist = 0.01;
      long i_contact_group = 1;
      mmdb::mat44 my_matt;
      mmdb::SymOps symm;
      for (int i=0; i<4; i++)
	 for (int j=0; j<4; j++)
	    my_matt[i][j] = 0.0;
      for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms;
      res_central->GetAtomTable(residue_atoms, n_residue_atoms);

      mol->SeekContacts(residue_atoms, n_residue_atoms,
			residue_atoms, n_residue_atoms,
			0, max_dist,
			0, // in same residue
			pscontact, n_contacts,
			0, &my_matt, i_contact_group); // makes reverses also
      if (n_contacts > 0) {
	 if (pscontact) {
	    for (int i=0; i<n_contacts; i++) {
	       mmdb::Atom *neighb_at = residue_atoms[pscontact[i].id2];
	       double radius = get_vdw_radius_ligand_atom(neighb_at);
	       std::pair<mmdb::Atom *, double> p(neighb_at, radius);
	       ligand_atom_neighbour_map[pscontact[i].id1].push_back(p);
	    }
	 }
      }
   }
}




void
coot::atom_overlaps_container_t::make_overlaps() {

   if (! have_dictionary) {
      std::cout << "No dictionary" << std::endl;
      return;
   }

   // This is completely non-clever about finding which atoms are close
   // to each other - it checks all neighbour atoms against all
   // central-residue atoms.

   double dist_crit = 2.3; // Everything that is bonded is less than this
   double dist_crit_sqrd = (2*dist_crit) * (2*dist_crit);

   mmdb::PAtom *central_residue_atoms = 0;
   int n_central_residue_atoms;
   res_central->GetAtomTable(central_residue_atoms, n_central_residue_atoms);

   for (int j=0; j<n_central_residue_atoms; j++) {

      mmdb::Atom *cr_at = central_residue_atoms[j];
      // std::cout << "Surface points for ligand atom " << atom_spec_t(cr_at) << std::endl;
      clipper::Coord_orth co_cr_at = co(cr_at);
      double r_1 = get_vdw_radius_ligand_atom(cr_at);

      for (unsigned int i=0; i<neighbours.size(); i++) { 
	 mmdb::PAtom *residue_atoms = 0;
	 int n_residue_atoms;
	 neighbours[i]->GetAtomTable(residue_atoms, n_residue_atoms);
	 for (int iat=0; iat<n_residue_atoms; iat++) { 
	    mmdb::Atom *n_at = residue_atoms[iat];
	    clipper::Coord_orth co_n_at = co(n_at);

	    double ds = (co_cr_at - co_n_at).lengthsq();
	    if (ds < dist_crit_sqrd) {

	       double r_2 = get_vdw_radius_neighb_atom(n_at, i);
	       double d = sqrt(ds);

	       // first is yes/no, second is if the H is on the ligand
	       // 
	       std::pair<bool, bool> might_be_h_bond_flag = is_h_bond_H_and_acceptor(cr_at, n_at);

	       if (d < (r_1 + r_2 + probe_radius)) { 
		  double o = get_overlap_volume(d, r_2, r_1);
		  bool h_bond_flag = false;
		  if (might_be_h_bond_flag.first) {
		     h_bond_flag = true;
		     if (might_be_h_bond_flag.second) { 
			std::cout << atom_spec_t(cr_at) << "   " << " and " << atom_spec_t(n_at)
				  << " r_1 " << r_1 << " and r_2 " << r_2  <<  " and d " << d 
				  << " overlap " << o << " IS H-Bond (ligand donor)" << std::endl;
		     } else { 
			   std::cout << atom_spec_t(cr_at) << "   " << " and " << atom_spec_t(n_at)
				     << " r_1 " << r_1 << " and r_2 " << r_2  <<  " and d " << d 
				     << " overlap " << o << " IS H-Bond (ligand acceptor)" << std::endl;
			   
		     }
		  } else { 
		     std::cout << atom_spec_t(cr_at) << "   " << " and " << atom_spec_t(n_at)
			       << " r_1 " << r_1 << " and r_2 " << r_2  <<  " and d " << d 
			       << " overlap " << o << std::endl;
		  }
		  atom_overlap_t ao(j, cr_at, n_at, r_1, r_2, o);
		  // atom_overlap_t ao2(-1, n_at, cr_at, r_2, r_1, o);
		  ao.is_h_bond = h_bond_flag;
		  overlaps.push_back(ao);
		  // overlaps.push_back(ao2);
	       } else {
		  if (might_be_h_bond_flag.first) { 
		     if (d < (dist_crit + 0.5)) {
			std::cout << atom_spec_t(cr_at) << "   " << " and " << atom_spec_t(n_at)
				  << " r_1 " << r_1 << " and r_2 " << r_2  <<  " and d " << d
				  << " but might be h-bond anyway (is this strange?)"
				  << std::endl;
			double o = 0;
			// atom_overlap_t ao(cr_at, n_at, r_1, r_2, o);
			// ao.is_h_bond = true;
			// overlaps.push_back(ao);
		     }
		  }
	       }
	    }
	 }
      }
   }
}

// first is yes/no, second is if the H is on the ligand
// 
std::pair<bool, bool>
coot::atom_overlaps_container_t::is_h_bond_H_and_acceptor(mmdb::Atom *ligand_atom,
							  mmdb::Atom *env_atom
							  // const double &d
							  ) const {

   bool status = false;
   bool H_on_ligand = false;
   unsigned int n_h = 0;
   bool h_on_ligand_atom = false;
   bool h_on_env_atom    = false;

   int hb_1 = -1;
   int hb_2 = -1;
   if (ligand_atom->GetUDData(udd_h_bond_type_handle, hb_1) == mmdb::UDDATA_Ok) { 
      if (env_atom->GetUDData(udd_h_bond_type_handle, hb_2) == mmdb::UDDATA_Ok) {

	 if (false) // testing
	    std::cout << "    hb_1 " << hb_1 << " for " << coot::atom_spec_t(ligand_atom)
		      << " and hb_2 " << hb_2 << " for neighb-atom " << coot::atom_spec_t(env_atom)
		      << std::endl;

	 if (hb_1 == HB_HYDROGEN) {
	    if (hb_2 == HB_ACCEPTOR || hb_2 == HB_BOTH) {
	       status = true;
	       H_on_ligand = true;
	    }
	 } 

	 if (hb_1 == HB_ACCEPTOR || hb_1 == HB_BOTH) {
	    if (hb_2 == HB_HYDROGEN) {
	       status = true;
	       H_on_ligand = false;
	    }
	 }
      } else {
	 // it's not bad.  Some H atoms dont have hydrogen UDD.
	 if (false)
	    std::cout << "   bad get of uud_h_bond info for env atom " << coot::atom_spec_t(env_atom)
		      << std::endl;
      }
   } else {
      if (false)
	 std::cout << "   bad get of uud_h_bond info for ligand atom " << coot::atom_spec_t(ligand_atom)
		   << std::endl;
   }
   return std::pair<bool, bool> (status, H_on_ligand);
}

// in A^3
double
coot::atom_overlaps_container_t::get_overlap_volume(const double &d, const double &r_1, const double &r_2) const {

   // V = π /(12d) * (r1+r2−d) ^2 * (d^2+2d(r1+r2)−3*(r1−r2)^2)

   double V = (M_PI/(12*d)) * (r_1+r_2-d) * (r_1+r_2-d) * (d*d + 2.0*d*(r_1+r_2) - 3*(r_1-r_2)*(r_1-r_2));
   return V;
} 

double
coot::atom_overlaps_container_t::get_vdw_radius_ligand_atom(mmdb::Atom *at) {

   double r = 2.5;

   std::map<mmdb::Atom *, double>::const_iterator it = central_residue_atoms_vdw_radius_map.find(at);
   if (it == central_residue_atoms_vdw_radius_map.end()) {
      // we need to add it then

      // What's the energy type of Atom at?
      //
      // PDBv3 FIXME - change from 4-char
      std::string te = central_residue_dictionary.type_energy(at->GetAtomName());
      if (! te.empty()) {
	 std::map<std::string, double>::const_iterator it_type =
	    type_to_vdw_radius_map.find(te);
	 if (it_type == type_to_vdw_radius_map.end()) {
	    // didn't find it. so look it up and add it.
	    r = geom_p->get_energy_lib_atom(te).vdw_radius;
	    hb_t t = geom_p->get_energy_lib_atom(te).hb_type;

	    if (false)
	       std::cout << "setting map: type_to_vdw_radius_h_bond_type_map[" << te << "] to ("
			 << r << "," << t << ")" << std::endl;

	    type_to_vdw_radius_map[te] = r;
	 } else {
	    r = it_type->second;
	 }
	 central_residue_atoms_vdw_radius_map[at] = r;
      } else {
	 std::cout << "failed to find type-energy for atom " << atom_spec_t(at) << std::endl;
      }
   } else {
      if (false)
	 std::cout << "radius for atom " << atom_spec_t(at) << " was found in map: value: "
		   << it->second << std::endl;
      r = it->second;
   }

   return r;
}


double
coot::atom_overlaps_container_t::get_vdw_radius_neighb_atom(int idx_neigh_atom) const {

   // no index checking (a bit cowboy?)
   //
   double r = neighb_atom_radius[idx_neigh_atom];
   return r;

}

double
coot::atom_overlaps_container_t::get_vdw_radius_neighb_atom(mmdb::Atom *at, unsigned int idx_res) {

   double r = 1.5;

   std::map<mmdb::Atom *, double>::const_iterator it = neighbour_atoms_vdw_radius_map.find(at);
   if (it == neighbour_atoms_vdw_radius_map.end()) {
      // we need to add it then

      // What's the energy type of Atom at?
      //
      std::string te = neighb_dictionaries[idx_res].type_energy(at->GetAtomName());
      std::map<std::string, double>::const_iterator it_type = type_to_vdw_radius_map.find(te);
      if (it_type == type_to_vdw_radius_map.end()) {
	 // didn't find te in types map. so look it up from the dictionary and add to the types map
	 r = geom_p->get_energy_lib_atom(te).vdw_radius;
	 hb_t t = geom_p->get_energy_lib_atom(te).hb_type;
	 type_to_vdw_radius_map[te] = r;
      } else {
	 r = it_type->second;
      }
      neighbour_atoms_vdw_radius_map[at] = r;
   } else {
      r = it->second;
   }

   return r;
} 



coot::hb_t 
coot::atom_overlaps_container_t::get_h_bond_type(mmdb::Atom *at) {

   hb_t type = HB_UNASSIGNED;
   std::string atom_name = at->name;
   std::string res_name = at->GetResName();
   type = geom_p->get_h_bond_type(atom_name, res_name, protein_geometry::IMOL_ENC_ANY); // heavyweight

   return type;
} 


void
coot::atom_overlaps_container_t::contact_dots_for_overlaps() const {

   double spike_length = 0.5;
   double clash_dist = 0.4;

   double dot_density = 0.35;
   dot_density = 1.0;

   double   phi_step = 5.0 * (M_PI/180.0);
   double theta_step = 5.0 * (M_PI/180.0);
   if (dot_density > 0.0) {
      phi_step   /= dot_density;
      theta_step /= dot_density;
   }

   for (unsigned int i=0; i<overlaps.size(); i++) {

      // std::cout << "considering overlap idx: " << i << std::endl;

      clipper::Coord_orth pt_at_1(overlaps[i].atom_1->x,
				  overlaps[i].atom_1->y,
				  overlaps[i].atom_1->z);
      clipper::Coord_orth pt_at_2(overlaps[i].atom_2->x,
				  overlaps[i].atom_2->y,
				  overlaps[i].atom_2->z);
      const double &r_1 = overlaps[i].r_1;
      const double &r_2 = overlaps[i].r_2;
      const int &idx =    overlaps[i].ligand_atom_index;
      double r_2_sqrd = r_2 * r_2;
      double r_2_plus_prb_squard  = r_2_sqrd + 2 * r_2 * probe_radius + probe_radius * probe_radius;
      double r_2_minux_prb_squard = r_2_sqrd - 2 * r_2 * probe_radius + probe_radius * probe_radius;
      //
      bool done_atom_name = false;

      bool even = true;
      for (double theta=0; theta<M_PI; theta+=theta_step) {
	 double phi_step_inner = phi_step + 0.1 * pow(theta-0.5*M_PI, 2);
	 for (double phi=0; phi<2*M_PI; phi+=phi_step_inner) {
	    if (even) {
	       clipper::Coord_orth pt(r_1*cos(phi)*sin(theta),
				      r_1*sin(phi)*sin(theta),
				      r_1*cos(theta));
	       clipper::Coord_orth pt_at_surface = pt + pt_at_1;
	       double d_sqrd = (pt_at_2 - pt_at_surface).lengthsq();

	       // std::cout << "comparing " << sqrt(d_sqrd) << " vs " << sqrt(r_2_plus_prb_squard)
	       // << " with r_2 " << r_2 << std::endl;

	       if (d_sqrd > r_2_plus_prb_squard) {

		  bool draw_it = ! is_inside_another_ligand_atom(idx, pt_at_surface);

		  if (false) // debugging
		     if (std::string(overlaps[i].atom_1->name) != " HO3")
			draw_it = false;

		  if (draw_it) {

		     std::cout << "considering overlap idx: " << i << " " << atom_spec_t(overlaps[i].atom_1)
			       << " to " << atom_spec_t(overlaps[i].atom_2) << std::endl;

		     clipper::Coord_orth vect_to_pt_1 = pt_at_1 - pt_at_surface;
		     clipper::Coord_orth vect_to_pt_1_unit(vect_to_pt_1.unit());
		     clipper::Coord_orth pt_spike_inner = pt_at_surface + spike_length * vect_to_pt_1_unit;

		     std::string type = "wide-contact";
		     bool only_pt = true; // not a spike

		     if (d_sqrd < r_2_sqrd)
			type = "close-contact";

		     if (d_sqrd < (r_2_sqrd - 2 * r_2 * clash_dist + clash_dist * clash_dist)) {
			type = "clash";
			only_pt = false;
		     }

		     if (overlaps[i].is_h_bond) {
			type = "H-bond";
		     }

		     if (! done_atom_name) {
			std::cout << "   spike for atom " << coot::atom_spec_t(overlaps[i].atom_1)
				  << std::endl;
			done_atom_name = true;
		     }

		     if (only_pt)
			pt_spike_inner = pt_at_surface;

		     // on the surface of atom_1 inside the sphere of atom_2
		     std::cout << "spike "
			       << type << " "
			       << pt_at_surface.x() << " "
			       << pt_at_surface.y() << " "
			       << pt_at_surface.z() << " to "
			       << pt_spike_inner.x() << " "
			       << pt_spike_inner.y() << " "
			       << pt_spike_inner.z()
			       << " theta " << theta << " phi " << phi
			       << std::endl;
		  }
	       }
	    }
	 }
      }
   }
}

bool
coot::atom_overlaps_container_t::is_inside_another_ligand_atom(int idx,
							       const clipper::Coord_orth &dot_pt) const {
   bool r = false;

   if (idx >= 0) {

      const std::vector<std::pair<mmdb::Atom *, double> > &v = ligand_atom_neighbour_map.find(idx)->second;
      for (unsigned int i=0; i<v.size(); i++) {
	 clipper::Coord_orth pt = co(v[i].first);
	 double dist_sqrd = (dot_pt - pt).lengthsq();

	 const double &radius_other = v[i].second;
	 if (dist_sqrd < radius_other * radius_other) {
	    r = true;
	    break;
	 }
      }
   }
   return r;
}
bool
coot::atom_overlaps_container_t::is_inside_another_ligand_atom(int idx,
							       const clipper::Coord_orth &probe_pos,
							       const clipper::Coord_orth &dot_pt) const {

   // for better speed, change the atom -> radius map to atom index -> radius vector
   // (and directly look up the radius).

   bool consider_probe_also = false; // if this is true, then don't make surface points for inner cusps
   
   bool r = false;

   if (idx >= 0) {

      const std::vector<std::pair<mmdb::Atom *, double> > &v = ligand_atom_neighbour_map.find(idx)->second;
      for (unsigned int i=0; i<v.size(); i++) {
	 clipper::Coord_orth pt = co(v[i].first);
	 double dist_sqrd = (dot_pt - pt).lengthsq(); // should be r_1 squared, right?
	 double radius_other = v[i].second;
	 radius_other += probe_radius;
	 if (dist_sqrd < radius_other * radius_other) {
	    r = true;
	    break;
	 }
      }
   }
   return r;
}


coot::atom_overlaps_dots_container_t
coot::atom_overlaps_container_t::contact_dots() {

   atom_overlaps_dots_container_t ao;
   mmdb::realtype max_dist = 4.0; // max distance for an interaction
   if (mol) {
      mmdb::Contact *pscontact = NULL;
      int n_contacts;
      float min_dist = 0.01;
      long i_contact_group = 1;
      mmdb::mat44 my_matt;
      mmdb::SymOps symm;
      for (int i=0; i<4; i++)
	 for (int j=0; j<4; j++)
	    my_matt[i][j] = 0.0;
      for (int i=0; i<4; i++) my_matt[i][i] = 1.0;

      mmdb::Atom **ligand_residue_atoms = 0;
      int n_ligand_residue_atoms;
      res_central->GetAtomTable(ligand_residue_atoms, n_ligand_residue_atoms);

      std::vector<residue_spec_t> env_residue_specs(neighbours.size());
      for (unsigned int i=0; i<neighbours.size(); i++)
	 env_residue_specs[i] = residue_spec_t(neighbours[i]);
      add_residue_neighbour_index_to_neighbour_atoms();

      int mask_mode = 0; // all atoms
      int i_sel_hnd_env_atoms = specs_to_atom_selection(env_residue_specs, mol, mask_mode);

      mmdb::Atom **env_residue_atoms = 0;
      int n_env_residue_atoms;
      mol->GetSelIndex(i_sel_hnd_env_atoms, env_residue_atoms, n_env_residue_atoms);
      setup_env_residue_atoms_radii(i_sel_hnd_env_atoms);

      mol->SeekContacts(ligand_residue_atoms, n_ligand_residue_atoms,
			env_residue_atoms, n_env_residue_atoms,
			0, max_dist,
			1, // 0: in same residue also?
			pscontact, n_contacts,
			0, &my_matt, i_contact_group);
      if (n_contacts > 0) {
	 if (pscontact) {
	    // std::cout << "n_contacts: " << n_contacts << std::endl;
	    for (int i=0; i<n_contacts; i++)
	       ligand_to_env_atom_neighbour_map[pscontact[i].id1].push_back(pscontact[i].id2);

	    // std::map<int, std::vector<mmdb::Atom *> >::const_iterator it;
	    std::map<int, std::vector<int> >::const_iterator it;
	    for (it= ligand_to_env_atom_neighbour_map.begin();
		 it!=ligand_to_env_atom_neighbour_map.end();
		 it++) {

	       mmdb::Atom *cr_at = ligand_residue_atoms[it->first];

	       // std::cout << "Surface points for " << atom_spec_t(cr_at) << std::endl;

	       clipper::Coord_orth pt_at_1 = co(cr_at);

	       double dot_density = 0.35;
	       dot_density = 1.05;
	       if (std::string(cr_at->element) == " H")
		  dot_density *=0.66; // so that surface dots on H atoms don't appear (weirdly) more fine

	       double   phi_step = 5.0 * (M_PI/180.0);
	       double theta_step = 5.0 * (M_PI/180.0);
	       if (dot_density > 0.0) {
		  phi_step   /= dot_density;
		  theta_step /= dot_density;
	       }

	       double r_1 = get_vdw_radius_ligand_atom(cr_at);

	       bool even = true;
	       for (double theta=0; theta<M_PI; theta+=theta_step) {
		  double phi_step_inner = phi_step + 0.1 * pow(theta-0.5*M_PI, 2);
		  for (double phi=0; phi<2*M_PI; phi+=phi_step_inner) {
		     if (even) {
			clipper::Coord_orth pt(r_1*cos(phi)*sin(theta),
					       r_1*sin(phi)*sin(theta),
					       r_1*cos(theta));
			clipper::Coord_orth pt_at_surface = pt + pt_at_1;
			bool draw_it = ! is_inside_another_ligand_atom(it->first, pt_at_surface);

			// if (std::string(cr_at->name) != " O1 ") draw_it = false; // debugging

			if (draw_it) {

			   // is the point on the surface of cr_at inside the sphere
			   // of an environment atom?
			   // If so, we want to know which one to which it was closest

			   double biggest_overlap = -1; // should be positive if we get a hit
			   mmdb::Atom *atom_with_biggest_overlap = 0;
			   double r_2_for_biggest_overlap = 0;

			   if (false) { // debug
			      std::cout << ":::: atom " << atom_spec_t(cr_at)
					<< " has radius " << r_1 << " and " << it->second.size()
					<< " env neighbours" << std::endl;
			      for (unsigned int ii=0; ii<it->second.size(); ii++)
				 std::cout << "   " << atom_spec_t(env_residue_atoms[it->second[ii]])
					   << std::endl;
			   }

			   for (unsigned int jj=0; jj<it->second.size(); jj++) {

			      mmdb::Atom *neighb_atom = env_residue_atoms[it->second[jj]];
			      double r_2 = get_vdw_radius_neighb_atom(it->second[jj]);
			      double r_2_sqrd = r_2 * r_2;
			      double r_2_plus_prb_squard =
				 r_2_sqrd + 2 * r_2 * probe_radius + probe_radius * probe_radius;

			      clipper::Coord_orth pt_na = co(neighb_atom);
			      double d_sqrd = (pt_na - pt_at_surface).lengthsq();

			      if (false)
				 std::cout << " for atom "
					   << atom_spec_t(env_residue_atoms[it->second[jj]])
					   << " comparing " << d_sqrd << " vs "
					   << r_2_plus_prb_squard << " with r_2 " << r_2
					   << " probe_radius " << probe_radius << std::endl;

			      if (d_sqrd < r_2_plus_prb_squard) {

				 // OK it was close to something.
				 double delta_d_sqrd = r_2_plus_prb_squard - d_sqrd;
				 if (delta_d_sqrd > biggest_overlap) {
				    biggest_overlap = delta_d_sqrd;
				    atom_with_biggest_overlap = neighb_atom;
				    r_2_for_biggest_overlap = r_2;
				 }
			      }
			   }

			   if (atom_with_biggest_overlap) {
			      double d_surface_pt_to_atom_sqrd =
				 (co(atom_with_biggest_overlap)-pt_at_surface).lengthsq();
			      double d_surface_pt_to_atom = sqrt(d_surface_pt_to_atom_sqrd);
			      double overlap_delta = r_2_for_biggest_overlap - d_surface_pt_to_atom;

			      // first is yes/no, second is H-is-on-ligand?
			      std::pair<bool, bool> might_be_h_bond_flag =
				 is_h_bond_H_and_acceptor(cr_at, atom_with_biggest_overlap);
			      bool is_h_bond = false;
			      if (might_be_h_bond_flag.first)
				 is_h_bond = true;
			      std::string c_type = overlap_delta_to_contact_type(overlap_delta, is_h_bond);

			      if (false)
				 std::cout << "spike "
					   << c_type << " "
					   << pt_at_surface.x() << " "
					   << pt_at_surface.y() << " "
					   << pt_at_surface.z() << " to "
					   << pt_at_surface.x() << " "
					   << pt_at_surface.y() << " "
					   << pt_at_surface.z()
					   << " theta " << theta << " phi " << phi
					   << std::endl;

			      if (c_type != "clash") {
				 ao.dots[c_type].push_back(pt_at_surface);
			      } else {
				 clipper::Coord_orth vect_to_pt_1 = pt_at_1 - pt_at_surface;
				 clipper::Coord_orth vect_to_pt_1_unit(vect_to_pt_1.unit());
				 clipper::Coord_orth pt_spike_inner =
				    pt_at_surface + clash_spike_length * vect_to_pt_1_unit;
				 std::pair<clipper::Coord_orth, clipper::Coord_orth> p(pt_at_surface,
										       pt_spike_inner);
				 ao.spikes.positions.push_back(p);
			      }

			   } else {

			      // no environment atom was close to this ligand atom, so just add
			      // a surface point

			      if (false)
				 std::cout << "spike-surface "
					   << pt_at_surface.x() << " "
					   << pt_at_surface.y() << " "
					   << pt_at_surface.z() << std::endl;

			      ao.dots["vdw-surface"].push_back(pt_at_surface);
			   }
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   return ao;
}

// return H-bond, or wide-contact or close-contact or small-overlap or big-overlap or clash
//
std::string
coot::atom_overlaps_container_t::overlap_delta_to_contact_type(double delta, bool is_h_bond) const {

// from the Word 1999 paper:
// pale green for H-bonds
// green (narrow gaps) or yellow (slight overlaps, < 0.2) for good contacts
// blue for wider gaps > 0.25
// orange and red for unfavourable (0.25 to 0.4)
// hot pink for >= 0.4

   std::string r = "wide-contact";

   // std::cout << "overlap-delta " << delta << " " << is_h_bond << std::endl;

   if (is_h_bond) {
      delta -= 0.7;
      if (delta > 0.4)
	 r = "clash";
      else
	 r = "H-bond";
   } else {
      if (delta > -0.1)         // Word: -0.25, // -0.15 allows too much green
	 r = "close-contact";
      if (delta > 0.07)          // Word: 0
	 r = "small-overlap";
      if (delta > 0.33)           // Word: 0.2  // 0.15, 0.18, 0.2, 0.25, 0.3 allows too much red
	 r = "big-overlap";
      if (delta > 0.4)          // Word: 0.4 // 0.3, 0.35 too much clash
	 r = "clash";
   }
   return r;
}

void
coot::atom_overlaps_container_t::add_residue_neighbour_index_to_neighbour_atoms() {

   udd_residue_index_handle = mol->RegisterUDInteger(mmdb::UDR_ATOM, "neighb-residue-index");
   for (unsigned int i=0; i<neighbours.size(); i++) {
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms;
      neighbours[i]->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
	 mmdb::Atom *at = residue_atoms[iat];
	 at->PutUDData(udd_residue_index_handle, int(i));
      }
   }
}



// fill std::vector<double> env_residue_radii;
void
coot::atom_overlaps_container_t::setup_env_residue_atoms_radii(int i_sel_hnd_env_atoms) {

   double r = 1.5;
   mmdb::Atom **env_residue_atoms = 0;
   int n_env_residue_atoms;
   mol->GetSelIndex(i_sel_hnd_env_atoms, env_residue_atoms, n_env_residue_atoms);
   neighb_atom_radius.resize(n_env_residue_atoms);

   for (int i=0; i<n_env_residue_atoms; i++) {
      mmdb::Atom *at = env_residue_atoms[i];
      mmdb::Residue *res = at->residue;
      int residue_index;
      at->GetUDData(udd_residue_index_handle, residue_index);
      const dictionary_residue_restraints_t &rest = neighb_dictionaries[residue_index];
      if (false) // debugging residue indexing
	 std::cout << "residue name " << res->GetResName() << " with comp_id "
		   << rest.residue_info.comp_id << std::endl;
      std::string te = rest.type_energy(at->GetAtomName());
      std::map<std::string, double>::const_iterator it_type = type_to_vdw_radius_map.find(te);
      if (it_type == type_to_vdw_radius_map.end()) {
	 // didn't find te in types map. so look it up from the dictionary and add to the types map
	 r = geom_p->get_energy_lib_atom(te).vdw_radius;
	 hb_t t = geom_p->get_energy_lib_atom(te).hb_type;
	 type_to_vdw_radius_map[te] = r;
      } else {
	 r = it_type->second;
      }
      neighb_atom_radius[i] = r;
   }
}
