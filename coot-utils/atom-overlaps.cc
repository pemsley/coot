

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
   init();

}

coot::atom_overlaps_container_t::atom_overlaps_container_t(mmdb::Residue *res_central_in,
							   mmdb::Residue *neighbour,
							   mmdb::Manager *mol,
							   const protein_geometry *geom_p_in) {
   geom_p = geom_p_in;
   res_central = res_central_in;
   neighbours.push_back(neighbour);
   init();
   
}


void
coot::atom_overlaps_container_t::init() {

   have_dictionary = false; // initially.
   
   if (res_central) { 

      std::string cres_name = res_central->GetResName();
      std::pair<bool, dictionary_residue_restraints_t> d =
	 geom_p->get_monomer_restraints(cres_name);
      if (! d.first) {
	 std::cout << "Failed to get dictionary for " << cres_name << std::endl;
      } else {
	 // Happy path
	 central_residue_dictionary = d.second;

	 neighb_dictionaries.resize(neighbours.size());
	 have_dictionary = true;
	 for (unsigned int i=0; i<neighbours.size(); i++) {
	    std::string residue_name = neighbours[i]->GetResName();
	    d = geom_p->get_monomer_restraints(residue_name);
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
		  std::cout << "----- adding HB_HYDROGEN udd " << atom_name << std::endl;
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
	    mmdb::PAtom *residue_atoms = 0;
	    int n_residue_atoms;
	    neighbours[i]->GetAtomTable(residue_atoms, n_residue_atoms);
	    for (int iat=0; iat<n_residue_atoms; iat++) { 
	       mmdb::Atom *n_at = residue_atoms[iat];
	       std::string atom_name(n_at->name);
	       std::string energy_type = neighb_dictionaries[i].type_energy(atom_name);
	       energy_lib_atom ela = geom_p->get_energy_lib_atom(energy_type);
	       hb_t hb_type = ela.hb_type;
	       n_at->PutUDData(udd_h_bond_type_handle, hb_type); // hb_t -> int
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
	       std::pair<bool, bool> might_be_h_bond_flag = is_h_bond_H_and_acceptor(cr_at, n_at, d);

	       if (d < (r_1 + r_2)) { 
		  double o = get_overlap_volume(d, r_2, r_1);
		  if (might_be_h_bond_flag.first) { 
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
		  atom_overlap_t ao(cr_at, n_at, o);
		  // ao.is_h_bond = h_bond_flag;
		  overlaps.push_back(ao);
	       } else {
		  if (might_be_h_bond_flag.first) { 
		     if (d < (dist_crit + 0.5)) {
			std::cout << atom_spec_t(cr_at) << "   " << " and " << atom_spec_t(n_at)
				  << " r_1 " << r_1 << " and r_2 " << r_2  <<  " and d " << d
				  << " but might be h-bond anyway " << std::endl;
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
							  mmdb::Atom *env_atom,
							  const double &d) const {

   bool status = false;
   bool H_on_ligand = false;
   unsigned int n_h = 0;
   bool h_on_ligand_atom = false;
   bool h_on_env_atom    = false;

   int hb_1 = -1;
   int hb_2 = -1;
   if (ligand_atom->GetUDData(udd_h_bond_type_handle, hb_1) == mmdb::UDDATA_Ok) { 
      if (env_atom->GetUDData(udd_h_bond_type_handle, hb_2) == mmdb::UDDATA_Ok) {
	 // std::cout << "hb_1 " << hb_1 << " and hb_2 " << hb_2 << std::endl;
	 if (hb_1 == HB_HYDROGEN) {
	    if (hb_2 == HB_ACCEPTOR || hb_2 == HB_BOTH) {
	       status = true;
	       H_on_ligand = true;
	    }
	 } 
      }

      if (hb_1 == HB_ACCEPTOR || hb_1 == HB_BOTH) {
	 if (hb_2 == HB_HYDROGEN) {
	    status = true;
	    H_on_ligand = false;
	 }
      }
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

   double r = 1.5;

   std::map<mmdb::Atom *, double>::const_iterator it = central_residue_atoms_vdw_radius_map.find(at);
   if (it == central_residue_atoms_vdw_radius_map.end()) {
      // we need to add it then

      // What's the energy type of Atom at?
      //
      std::string te = central_residue_dictionary.type_energy(at->GetAtomName());
      std::map<std::string, double>::const_iterator it_type =
	 type_to_vdw_radius_map.find(te);
      if (it_type == type_to_vdw_radius_map.end()) {
	 // didn't find it. so look it up and add it.
	 r = geom_p->get_energy_lib_atom(te).vdw_radius;
	 hb_t t = geom_p->get_energy_lib_atom(te).hb_type;
	 // std::cout << "setting map: type_to_vdw_radius_h_bond_type_map[" << te << "] to ("
	 // << r << "," << t << ")" << std::endl;
	 type_to_vdw_radius_map[te] = r;
      } else {
	 r = it_type->second;
      }
      central_residue_atoms_vdw_radius_map[at] = r;
   } else {
      r = it->second;
   } 

   return r;
}



double
coot::atom_overlaps_container_t::get_vdw_radius_neighb_atom(mmdb::Atom *at, unsigned int ires) {

   double r = 1.5;

   std::map<mmdb::Atom *, double>::const_iterator it = neighbour_atoms_vdw_radius_map.find(at);
   if (it == neighbour_atoms_vdw_radius_map.end()) {
      // we need to add it then

      // What's the energy type of Atom at?
      //
      std::string te = neighb_dictionaries[ires].type_energy(at->GetAtomName());
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
   type = geom_p->get_h_bond_type(atom_name, res_name); // heavyweight

   return type;
} 

