
#include <stdexcept>

#include <sys/types.h> // for stating
#include <sys/stat.h>

#include "protein-geometry.hh"

void
coot::protein_geometry::add_energy_lib_atom(    const coot::energy_lib_atom    &atom) {
   energy_lib.atom_map[atom.type] = atom;
} 

void
coot::protein_geometry::add_energy_lib_bond(const coot::energy_lib_bond &bond) {
   energy_lib.bonds.push_back(bond);
}


void
coot::protein_geometry::add_energy_lib_angle(   const coot::energy_lib_angle   &angle) {
   energy_lib.angles.push_back(angle);
}

void
coot::protein_geometry::add_energy_lib_torosion(const coot::energy_lib_torsion &torsion) {
   energy_lib.torsions.push_back(torsion);
}

void
coot::protein_geometry::read_energy_lib(const std::string &file_name) {

   struct stat buf;
   int istat = stat(file_name.c_str(), &buf);
   if (istat != 0) {
      std::cout << "WARNING:: energy lib " << file_name << " not found.\n";
      return;
   }
   
   CMMCIFFile ciffile;
   int ierr = ciffile.ReadMMCIFFile(file_name.c_str());
   if (ierr!=CIFRC_Ok) {
      std::cout << "dirty mmCIF file? " << file_name.c_str() << std::endl;
      std::cout << "    Bad CIFRC_Ok on ReadMMCIFFile" << std::endl;
      std::cout << "    " << GetErrorDescription(ierr) << std::endl;
      char        err_buff[1000];
      std::cout <<  "CIF error rc=" << ierr << " reason:" << 
	 GetCIFMessage (err_buff,ierr) << std::endl;
   } else {
      std::cout << "There are " << ciffile.GetNofData() << " data in "
		<< file_name << std::endl;
      for(int idata=0; idata<ciffile.GetNofData(); idata++) { 
	 PCMMCIFData data = ciffile.GetCIFData(idata);
	 // if (std::string(data->GetDataName()).substr(0,5) == "_lib_atom") {
	 // energy_lib_atoms(mmCIFLoop);

	 // std::cout << "read_energy_lib: " << data->GetDataName() << std::endl;

	 if (std::string(data->GetDataName()) == "energy") {
	    for (int icat=0; icat<data->GetNumberOfCategories(); icat++) { 
      
	       PCMMCIFCategory cat = data->GetCategory(icat);
	       std::string cat_name(cat->GetCategoryName());
	       
	       // std::cout << "DEBUG:: init_link is handling " << cat_name << std::endl;
	       
	       PCMMCIFLoop mmCIFLoop = data->GetLoop(cat_name.c_str());
	       
	       if (mmCIFLoop == NULL) { 
		  std::cout << "null loop" << std::endl; 
	       } else {
		  int n_chiral = 0;
		  if (cat_name == "_lib_atom")
		     add_energy_lib_atoms(mmCIFLoop);
		  if (cat_name == "_lib_bond")
		     add_energy_lib_bonds(mmCIFLoop);
	       }
	    }
	 } 
      }
   }
}


void
coot::protein_geometry::add_energy_lib_atoms(PCMMCIFLoop mmCIFLoop) {

   // note that that:
   // if (ierr) {
   //    xxx = -1
   // }
   // 
   // construction is necessary because if the GetXXX() fails then the
   // int/float is set to 0 or some such by the function.

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {
      std::string type;
      realtype weight = -1;
      int hb_type = coot::energy_lib_atom::HB_UNASSIGNED;
      realtype vdw_radius = -1;
      realtype vdwh_radius = -1; // with implicit hydrogen, I presume
      realtype ion_radius = -1;
      std::string element;
      int valency = -1;
      int sp_hybridisation = -1;
      int ierr;
      int ierr_tot = 0;

      char *s;

      s = mmCIFLoop->GetString("type", j, ierr);
      ierr_tot += ierr;
      if (s) type = s;
         
      // This can fail (to set weight - we still have a useful atom description).
      //
      ierr = mmCIFLoop->GetReal(weight, "weight", j);
      if (ierr) {
	 weight = -1;
      }

      s = mmCIFLoop->GetString("hb_type", j, ierr);
      ierr_tot += ierr;
      if (s) {
	 std::string ss(s);
	 if (ss == "D")
	    hb_type = coot::energy_lib_atom::HB_DONOR;
	 if (ss == "A")
	    hb_type = coot::energy_lib_atom::HB_ACCEPTOR;
	 if (ss == "B")
	    hb_type = coot::energy_lib_atom::HB_BOTH;
	 if (ss == "N")
	    hb_type = coot::energy_lib_atom::HB_NEITHER;
	 if (ss == "H")
	    hb_type = coot::energy_lib_atom::HB_HYDROGEN;
      }

      // This can fail (to set vdw_radius - we still have a useful atom description).
      //
      ierr = mmCIFLoop->GetReal(vdw_radius, "vdw_radius",j);
      if (ierr) {
	 vdw_radius = -1;
      }

      // This can fail (to set vdwh_radius - we still have a useful atom description).
      //
      ierr = mmCIFLoop->GetReal(vdwh_radius, "vdwh_radius",j);
      if (ierr) {
	 vdwh_radius = -1;
      }

      // This can fail (to set ion_radius - we still have a useful atom description).
      //
      ierr = mmCIFLoop->GetReal(ion_radius, "ion_radius", j);
      if (ierr) {
	 ion_radius = -1;
      }

      s = mmCIFLoop->GetString("element", j, ierr);
      ierr_tot += ierr;
      if (s)
	 element = s;
      
      ierr = mmCIFLoop->GetInteger(valency, "valency", j);
      ierr_tot += ierr;
      
      // This can fail (to set the hybridisation - we still have a useful atom description).
      //
      ierr = mmCIFLoop->GetInteger(sp_hybridisation, "sp", j);
      if (ierr) {
	 sp_hybridisation = -1;
      }

      if (ierr_tot == 0) {
	 coot::energy_lib_atom at(type, hb_type, weight, vdw_radius, vdwh_radius,
				  ion_radius, element, valency, sp_hybridisation);
	 // std::cout << "DEBUG:: adding energy atom: " << at << std::endl;
	 add_energy_lib_atom(at);
      } 
      
   }
}

std::ostream&
coot::operator<<(std::ostream &s, const energy_lib_atom &at) {

   s << "[type: " << at.type << " weight: " << at.weight << " vdw_radius: " << at.vdw_radius
     << " vdwh_radius: " << at.vdwh_radius << " ion_radius: " << at.ion_radius
     << " element: " << at.element
     << " valency: " << at.valency << " sp_hybridisation: " << at.sp_hybridisation
     << "]";
   return s;
}

std::ostream&
coot::operator<<(std::ostream &s, const energy_lib_bond &bond) {

   s << "[type: " << bond.atom_type_1 << " " << bond.atom_type_2 << " "
     << bond.length << " " << bond.esd << "]";
   return s;
}


void
coot::protein_geometry::add_energy_lib_bonds(PCMMCIFLoop mmCIFLoop) {

   for (int j=0; j<mmCIFLoop->GetLoopLength(); j++) {
      
      std::string atom_type_1;
      std::string atom_type_2;
      std::string type;
      realtype spring_const;
      realtype length;
      realtype value_esd; 
      int ierr;
      int ierr_tot = 0;

      char *s;

      s = mmCIFLoop->GetString("atom_type_1", j, ierr);
      ierr_tot += ierr;
      if (s) atom_type_1 = s;

      s = mmCIFLoop->GetString("atom_type_2", j, ierr);
      ierr_tot += ierr;
      if (s) atom_type_1 = s;

      s = mmCIFLoop->GetString("type", j, ierr);
      ierr_tot += ierr;
      if (s) type = s;
      
      ierr = mmCIFLoop->GetReal(spring_const, "const", j);
      if (ierr) {
	 spring_const = 420;
      }
      
      ierr = mmCIFLoop->GetReal(length, "length", j);
      ierr_tot += ierr;

      ierr = mmCIFLoop->GetReal(value_esd, "value_esd", j);
      if (ierr) {
	 value_esd = 0.02;
      }

      if (ierr_tot == 0) {
	 coot::energy_lib_bond bond(atom_type_1, atom_type_1, type, spring_const, length, value_esd);
	 add_energy_lib_bond(bond);
      }
   } 
} 




int
coot::protein_geometry::get_h_bond_type(const std::string &atom_name, const std::string &monomer_name) const {

   bool debug = 0;  // before debugging this, is ener_lib.cif being
		    // read correctly?
   
   // this is heavy!
   // 
   std::pair<bool, coot::dictionary_residue_restraints_t> r =
      get_monomer_restraints(monomer_name);

   int hb_type = coot::energy_lib_atom::HB_UNASSIGNED;

   if (! r.first) {
      std::string m = "No dictionary for monomer_type: ";
      m += monomer_name;
      std::cout << m << std::endl;
   } else {
      for (unsigned int i=0; i<r.second.atom_info.size(); i++) {
	 if (r.second.atom_info[i].atom_id_4c == atom_name) { 
	    std::string type = r.second.atom_info[i].type_energy;
	    if (type.length()) {
	       std::map<std::string, coot::energy_lib_atom>::const_iterator it = 
		  energy_lib.atom_map.find(type);
	       if (it != energy_lib.atom_map.end()) { 
		  hb_type = it->second.hb_type;
		  if (debug)
		     std::cout << "DEBUG:: found hb_type " << hb_type << " for " << atom_name
			       << " given energy type " << type << std::endl;
	       }
	    }
	    break;
	 }
      }
   } 

   if (debug)
      if (hb_type == coot::energy_lib_atom::HB_UNASSIGNED)
	 std::cout << " failed to get_h_bond_type for " << atom_name << " in " << monomer_name
		   << std::endl;
	 
   return hb_type;

} 


// Find the non_bonded_contact distance, get_nbc_dist()
// 
// Return a pair, if not found the first is 0.  Look up in the energy_lib.
// 
std::pair<bool, double>
coot::protein_geometry::get_nbc_dist(const std::string &energy_type_1,
				     const std::string &energy_type_2) const {

   std::pair<bool, double> r(false, 0);
   std::map<std::string, coot::energy_lib_atom>::const_iterator it_1 = energy_lib.atom_map.find(energy_type_1);
   std::map<std::string, coot::energy_lib_atom>::const_iterator it_2 = energy_lib.atom_map.find(energy_type_2);

   if (it_1 != energy_lib.atom_map.end()) { 
      if (it_2 != energy_lib.atom_map.end()) {
	 r.first = true;
	 r.second = it_1->second.vdw_radius + it_2->second.vdw_radius;

	 // hack in a correction factor - that makes atoms fit into density  :-/
	 // 
	 // 0.9 is too much.

	 // There is an algorithm problem.
	 // 
	 // I need to reject atoms that have bond (done), angle and
	 // torsion interactions from NBC interactions, I think.  Then
	 // I can set r.second to 1.0 maybe.
	 // 
	 r.second *= 0.84;

	 // ring atoms should not be NBCed to each other.  Not sure
	 // that 5 atom rings need to be excluded in this manner.
	 // 
	 if (it_1->first == "CR15" || it_1->first == "CR16" || it_1->first == "CR1"  ||
	     it_1->first == "CR6"  || it_1->first == "CR5"  || it_1->first == "CR5"  ||
	     it_1->first == "CR56" || it_1->first == "CR5"  || it_1->first == "CR66" ||
	     it_1->first == "NPA"  || it_1->first == "NPB"  || it_1->first == "NRD5" ||
	     it_1->first == "NRD6" || it_1->first == "NR15" || it_1->first == "NR16" ||
	     it_1->first == "NR6"  || it_1->first == "NR5") {

	    if (it_2->first == "CR15" || it_2->first == "CR16" || it_2->first == "CR1"  ||
		it_2->first == "CR6"  || it_2->first == "CR5"  || it_2->first == "CR5"  ||
		it_2->first == "CR56" || it_2->first == "CR5"  || it_2->first == "CR66" ||
		it_2->first == "NPA"  || it_2->first == "NPB"  || it_2->first == "NRD5" ||
		it_2->first == "NRD6" || it_2->first == "NR15" || it_2->first == "NR16" ||
		it_2->first == "NR6"  || it_2->first == "NR5") {
	       r.second = 2.2;
	    }
	 }


	 // hydrogen bonds can be closer
	 // 
	 if ((it_1->second.hb_type == coot::energy_lib_atom::HB_DONOR ||
	      it_1->second.hb_type == coot::energy_lib_atom::HB_BOTH  ||
	      it_1->second.hb_type == coot::energy_lib_atom::HB_HYDROGEN) &&
	     (it_2->second.hb_type == coot::energy_lib_atom::HB_ACCEPTOR ||
	      it_2->second.hb_type == coot::energy_lib_atom::HB_BOTH)) { 
	    r.second -= 0.7;
	    // actual hydrogens to aceptors can be shorter still 
	    if (it_1->second.hb_type == coot::energy_lib_atom::HB_HYDROGEN)
	       r.second -=0.3;
	 } 

	 if ((it_2->second.hb_type == coot::energy_lib_atom::HB_DONOR ||
	      it_2->second.hb_type == coot::energy_lib_atom::HB_BOTH  ||
	      it_2->second.hb_type == coot::energy_lib_atom::HB_HYDROGEN) &&
	     (it_1->second.hb_type == coot::energy_lib_atom::HB_ACCEPTOR ||
	      it_1->second.hb_type == coot::energy_lib_atom::HB_BOTH)) { 
	    r.second -= 0.7;
	    // as above
	    if (it_1->second.hb_type == coot::energy_lib_atom::HB_HYDROGEN)
	       r.second -=0.3;
	 }
      }
   }
   return r;
} 


coot::energy_lib_bond
coot::energy_lib_t::get_bond(const std::string &energy_type_1,
			     const std::string &energy_type_2) const {

   energy_lib_bond bond;
   std::map<std::string, energy_lib_atom>::const_iterator it_1 = atom_map.find(energy_type_1);
   std::map<std::string, energy_lib_atom>::const_iterator it_2 = atom_map.find(energy_type_2);

   if (it_1 != atom_map.end() && it_2 != atom_map.end()) {
      for (unsigned int ibond=0; ibond<bonds.size(); ibond++) { 
	 std::cout << "energy bond " << ibond << " " << bonds[ibond] << std::endl;
      }

   } else {
      std::cout << "WARNING:: in get_bond() failed to find energy types given " << energy_type_1
		<< " and " << energy_type_2 << std::endl;
   } 

   return bond;
} 
