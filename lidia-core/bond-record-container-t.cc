
#include <fstream>
#include <iomanip>
#include <list>

#include "rdkit-interface.hh"
#include "GraphMol/Rings.h"
#include "GraphMol/RingInfo.h"
#include "boost/dynamic_bitset.hpp"

#include "utils/coot-utils.hh"

#include "cod-atom-types.hh"
#include "bond-record-container-t.hh"

#include "coot-utils/coot-coord-utils.hh"

// can throw std::runtime_error.
// 
// The input atom names are directly from the pdbx cif file, so the
// don't contain spaces, so the input at_name_in should not contain
// spaces.
unsigned int
cod::bond_record_container_t::get_atom_index(const std::string &at_name_in,
					     const RDKit::RWMol &mol) const {

   // check "name" property

   unsigned int idx = 0;
   bool found = false;
   
   unsigned int n_mol_atoms = mol.getNumAtoms();
   for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
      RDKit::ATOM_SPTR at_p = mol[iat];
      try {
	 std::string name = "";
	 at_p->getProp("name", name);
	 // std::cout << "   debug " << iat << " \"" << name << "\"" << std::endl;
	 if (name == at_name_in) {
	    found = true;
	    idx = iat;
	 }
      }
      catch (const KeyErrorException &kee) {
      }
   }

   if (! found) {
      std::cout << "get_atom_index() throwing rte for atom \"" << at_name_in
		<< "\""<< std::endl;
      std::string m(std::string("atom name \"") + at_name_in +
		    std::string("\" not found in dictionary atom name list"));
      throw std::runtime_error(m);
   } 
   return idx;
}


// can throw std::runtime_error
unsigned int
cod::bond_record_container_t::get_atom_index(const std::string &at_name_1,
					     const coot::dictionary_residue_restraints_t &rest) const {


   unsigned idx = 0;
   bool found = false;
   for (unsigned int iat=0; iat<rest.atom_info.size(); iat++) { 
      if (rest.atom_info[iat].atom_id_4c == at_name_1) {
	 found = true;
	 break;
      }
   }

   if (! found) {
      std::string m(std::string("atom name ") + at_name_1 +
		    std::string(" not found in dictionary atom name list"));
      throw std::runtime_error(m);
   }
   return idx;
}



bool
cod::bond_record_container_t::write_atom_type_indices(const std::string &file_name) const {

   bool status = false;

   std::ofstream f(file_name.c_str());
   if (f) {
      std::map<std::string, unsigned int>::const_iterator it_map = atom_types_map.begin();
      unsigned int n = atom_types_map.size();
      for (unsigned int i=0; i<n; i++) {
	 f << std::setw(6) << it_map->second << " " << it_map->first << "\n";
	 it_map++;
      }
   }
   f.close();

   return status;
}


bool
cod::bond_record_container_t::write(const std::string file_name) const {

   write_atom_type_indices("atom-indices.tab");

   bool status = false;
   std::ofstream f(file_name.c_str());

   if (f) {
      for (unsigned int i=0; i<bonds.size(); i++) {
	 const bond_table_record_t &btr = bonds[i];
	 // f << bonds[i] << "\n";
	 // btr.write(f, max_atom_type_width);
	 // dangerous, but should not fail
	 btr.write(f,
		   atom_types_map.find(btr.cod_type_1)->second,
		   atom_types_map.find(btr.cod_type_2)->second);
      }
      f.close();
      status = true;
   }
   return status;
}

std::string::size_type
cod::bond_record_container_t::get_max_atom_type_width() const {

   std::size_t m = 0;
   for (unsigned int i=0; i<bonds.size(); i++) { 
      std::string::size_type s1 = bonds[i].cod_type_1.length();
      std::string::size_type s2 = bonds[i].cod_type_2.length();

      if (s1 > m) m = s1;
      if (s2 > m) m = s2;
   }

   return m;
}


bool
cod::bond_record_container_t::read_acedrg_table(const std::string &file_name) {


   bool status = false;
   std::ifstream f(file_name.c_str());

   if (f) {

      std::cout << "opened " << file_name << std::endl;
      std::string line;

      while (std::getline(f, line)) {
	 std::vector<std::string> bits = coot::util::split_string_no_blanks(line);

	 if (bits.size() == 18) {
	    const std::string &cod_type_1  = bits[10];
	    const std::string &cod_type_2  = bits[11];
	 
	    // const std::string &mean_str    = bits[12];
	    // const std::string &std_dev_str = bits[13];
	    // const std::string &count_str   = bits[14];

	    const std::string &mean_str    = bits[15];
	    const std::string &std_dev_str = bits[16];
	    const std::string &count_str   = bits[17];

	    try {
	       double mean   = coot::util::string_to_float(mean_str);
	       double stddev = coot::util::string_to_float(std_dev_str);
	       int count     = coot::util::string_to_int(count_str);

	       bond_table_record_t rec(cod_type_1, cod_type_2, mean, stddev, count);
	       add(rec);
	       status = true;
	    }
	    catch (const std::runtime_error &rte) {
	       std::cout << "error converting " << rte.what() << std::endl;
	    }
	 }
      }
   } else {
      std::cout << "failed to open " << file_name << std::endl;
   } 

   return status;
}



// return the consolidated table
//
void
cod::bond_record_container_t::read_acedrg_table_dir(const std::string &dir_name) {

   std::string glob_pattern = "*.table";
   std::vector<std::string> tables = coot::util::glob_files(dir_name, glob_pattern);

   for (unsigned int i=0; i<tables.size(); i++) {
      const std::string &table = tables[i];
      bond_record_container_t single;
      bool success = single.read_acedrg_table(table);
      if (success)
	 add_table(single);
   }

   // this should be a brc fuction
   
   std::cout << "stored " << size() << " bond records" << std::endl;
   std::cout << "-- pre-sort " << std::endl;
   sort();
   std::cout << "-- post-sort " << std::endl;

   fill_atom_map();
}


void
cod::bond_record_container_t::fill_atom_map() {

   // now fill the atom map
   // v1
   //
   if (false) { 
      std::list<std::string> l;
      std::cout << "pre  fill atom map" << std::endl;
      std::list<std::string>::const_iterator it;
      for (unsigned int i=0; i<bonds.size(); i++) {
	 // it = l.find(bonds[i].cod_type_1);
	 it = std::find(l.begin(), l.end(), bonds[i].cod_type_1);
	 if (it == l.end())
	    l.push_back(bonds[i].cod_type_1);
	 it = std::find(l.begin(), l.end(), bonds[i].cod_type_2);
	 if (it == l.end())
	    l.push_back(bonds[i].cod_type_2);
      }
      std::cout << "post fill atom map" << std::endl;
   }

   std::set<std::string> set;
   std::cout << "pre  fill atom set" << std::endl;
   for (unsigned int i=0; i<bonds.size(); i++) {
      set.insert(bonds[i].cod_type_1);
      set.insert(bonds[i].cod_type_2);
   }
   std::cout << "post fill atom set " << set.size() << std::endl;

   unsigned int n = set.size();
   std::set<std::string>::const_iterator it;
   unsigned int idx = 0;
   std::cout << "pre  fill atom map" << std::endl;
   for (it=set.begin(); it!=set.end(); it++) {
      // std::cout << "    " << *it << std::endl;
      atom_types_map[*it] = idx;
      idx++;
   }
   std::cout << "post fill atom map" << std::endl;

   std::map<std::string, unsigned int>::const_iterator it_map;
   it_map = atom_types_map.begin();
}

bool
cod::bond_record_container_t::read(const std::string &atom_type_indices_file_name,
				   const std::string &bonds_file_name) {

   bool success = true;
   std::vector<std::string> types = read_atom_type_indices(atom_type_indices_file_name);
   if (types.size())
      success = read_bonds(bonds_file_name, types);
   return success;
} 


std::vector<std::string>
cod::bond_record_container_t::read_atom_type_indices(const std::string &atom_type_indices_file_name) const {

   std::vector<std::string> types(300000);
   bool status = true;
   std::ifstream f(atom_type_indices_file_name.c_str());
   if (f) {
      std::string line;
      try { 
	 while (std::getline(f, line)) {
	    int idx = coot::util::string_to_int(line.substr(0, 6));
	    types[idx] = line.substr(7);
	 }
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "Error while reading " << atom_type_indices_file_name
		   << " - failed to parse: " << line << " " << rte.what() << std::endl;
      }
   } else {
      status = false;
   }
   return types;
}

bool
cod::bond_record_container_t::read_bonds(const std::string &bonds_file_name,
					 const std::vector<std::string> &atom_types) {

   // this is the slow function.
   //
   // It takes 3 seconds to run,
   // of which 1s is the bonds_map[][] line,
   // running getline() on the file takes 0.2s,
   // 1.8 s is spend on string -> number
   // the atom types look up takes ~0.1s each (really?)
   
   bool status = true;
   std::ifstream f(bonds_file_name.c_str());
   if (f) {
      // bonds.reserve(563711); // nice speedup.
      std::string line;
      try { 
	 while (std::getline(f, line)) {

	    if (true) {
	       double bl      = coot::util::string_to_double(line.substr( 0,10));
	       double std_dev = coot::util::string_to_double(line.substr(11,20));
	       int count      = coot::util::string_to_int(line.substr(21, 6));
	       int idx_1      = coot::util::string_to_int(line.substr(28, 6));
	       int idx_2      = coot::util::string_to_int(line.substr(36, 6));

	       std::string type_1 = atom_types[idx_1];
	       std::string type_2 = atom_types[idx_2];

	       if (false)
		  std::cout << "parsed "
			    << std::setw(10) << bl << " "
			    << std::setw(10) << std_dev << " "
			    << std::setw( 8) << count
			    << " idx_1: " << idx_1 << "    idx_2: " << idx_2 << "   "
			    << type_1 << "   " << type_2 << "\n";

	       if (type_2 < type_1)
		  std::swap(type_1, type_2);
	       bond_table_record_t bond(type_1, type_2, bl, std_dev, count);

	       bonds_map[type_1][type_2] = bond;
	    }
	 }
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "Error while reading " << bonds_file_name
		   << " - failed to parse: " << line << " " << rte.what() << std::endl;
      }

      if (false) { 
	 int n_forward = 0;
	 int n_backward = 0;
	 int n_equal = 0;
	 int n_many  = 0;

	 for (unsigned int i=0; i<bonds.size(); i++) { 
	    if (bonds[i].cod_type_1 > bonds[i].cod_type_2)
	       n_forward++;
	    if (bonds[i].cod_type_1 < bonds[i].cod_type_2)
	       n_backward++;
	    if (bonds[i].cod_type_1 == bonds[i].cod_type_2)
	       n_equal++;
	    if (bonds[i].count > 7)
	       n_many++;
	 }

	 std::cout << "           n_forward  " << n_forward << std::endl;
	 std::cout << "           n_backward " << n_backward << std::endl;
	 std::cout << "           n_equal    " << n_equal    << std::endl;
	 std::cout << "           n_many     " << n_many     << std::endl;
      }
      
      if (false)
	 for (unsigned int i=0; i<20; i++)
	    std::cout << "   " << bonds[i] << std::endl;
	 

   } else {
      status = false;
   }

   return status;
}

void
cod::bond_record_container_t::check() const {

   std::cout << "start check " << std::endl;
   std::map<std::string, std::vector<bond_table_record_t> > bonds_map;
   for (unsigned int i=0; i<bonds.size(); i++)
      bonds_map[bonds[i].cod_type_1].push_back(bonds[i]);
   std::cout << " done check - 0" << std::endl;

   std::map<std::string, std::vector<bond_table_record_t> >::const_iterator it_bonds_map;

   std::map<unsigned int, unsigned int> counts;
   std::map<unsigned int, unsigned int>::const_iterator it_counts;
   
   for (it_bonds_map=bonds_map.begin(); it_bonds_map!=bonds_map.end(); it_bonds_map++)
      counts[it_bonds_map->second.size()]++;
   std::cout << " done check - 1" << std::endl;

   if (false)
      for (it_counts=counts.begin(); it_counts!=counts.end(); it_counts++)
	 // second is the frequency
	 std::cout << "   " << it_counts->first << " " << it_counts->second
		   << std::endl;

   
   std::cout << " done check all" << std::endl;
}

std::vector<bool>
cod::bond_record_container_t::get_is_hydrogen_flags(const RDKit::RWMol &rdkm) const {

   unsigned int n_mol_atoms = rdkm.getNumAtoms();
   std::vector<bool> flags(n_mol_atoms);
   for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
      RDKit::ATOM_SPTR at_p = rdkm[iat];
      try {
	 if (at_p->getAtomicNum() == 1) {
	    flags[iat] = true;
	 } else {
	    flags[iat] = false;
	 }
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "   " << rte.what() << std::endl;
      }
   }
   return flags;
}

// atom types are generated from the atoms of rest (which contains hydrogens)
// 
void
cod::bond_record_container_t::validate(mmdb::Residue *res,
				       const coot::dictionary_residue_restraints_t &rest) const {

   if (res) { 
      std::string res_name = res->GetResName();
      std::cout << "validate: " << res_name << " " << rest.residue_info.comp_id
		<< std::endl;

      if (res_name == rest.residue_info.comp_id) {

	 try {
	    
	    RDKit::RWMol rdkm = coot::rdkit_mol(rest);

	    // coot::debug_rdkit_molecule(&rdkm);

	    atom_types_t t;
	    std::vector<std::string> v = t.get_cod_atom_types(rdkm, true);

	    if (rdkm.getNumAtoms() == v.size()) {
	       unsigned int n_mol_atoms = rdkm.getNumAtoms();
	       for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
		  RDKit::ATOM_SPTR at_p = rdkm[iat];
		  try {
		     std::string name;
		     at_p->getProp("name", name);
		     std::cout << "    " << iat << " " << name
			       << "  " << v[iat] << std::endl;
		  }
		  catch (const KeyErrorException &kee) {
		     std::cout << "   " << iat << " " << kee.what()
			       << std::endl;
		  }
	       }
	    }

	    if (false)
	       for (unsigned int ii=0; ii<v.size(); ii++)
		  std::cout << "get_cod_atom_types() " << ii << " " << v[ii]
			    << std::endl;

	    std::vector<bool> v_H = get_is_hydrogen_flags(rdkm);

	    // check that v is the same length as residue has atoms
	    //
	    if ((rdkm.getNumAtoms() == v.size()) && (v.size() == v_H.size())) {

	       for (unsigned int ibond=0; ibond<rest.bond_restraint.size(); ibond++) {
		  const coot::dict_bond_restraint_t &bond = rest.bond_restraint[ibond];

		  std::string at_name_1_4c = bond.atom_id_1_4c();
		  std::string at_name_2_4c = bond.atom_id_2_4c();

		  std::string at_name_1 = coot::util::remove_whitespace(bond.atom_id_1());
		  std::string at_name_2 = coot::util::remove_whitespace(bond.atom_id_2());

		  if (false) { 
		     std::cout << "in validate(): at_name_1 is \"" << at_name_1
			       << "\"" << std::endl;
		     std::cout << "               at_name_2 is \"" << at_name_2
			       << "\"" << std::endl;
		  }

		  try {
		     
		     // use atom indices and look in v. Use non-whitespace names
		     // 
		     unsigned int atom_idx_1 = get_atom_index(at_name_1, rdkm);
		     unsigned int atom_idx_2 = get_atom_index(at_name_2, rdkm);

		     if (false)
			std::cout << "atom indices "
				  << atom_idx_1 << " " << atom_idx_2 << std::endl;

		     if (! v_H[atom_idx_1] && ! v_H[atom_idx_2]) {

			std::string cod_type_1 = v[atom_idx_1];
			std::string cod_type_2 = v[atom_idx_2];

			if (cod_type_2 < cod_type_1)
			   std::swap(cod_type_1, cod_type_2);
		     
			bond_table_record_t cod_bond = 
			   get_cod_bond_from_table(cod_type_1, cod_type_2);
		     
			double dist =
			   get_bond_distance_from_model(at_name_1_4c, at_name_2_4c, res);

			std::cout << "compare: model: " << dist
				  << " vs tables: " << cod_bond.mean << " +/- "
				  << cod_bond.std_dev << " "
				  << cod_type_1 << "  " << cod_type_2
				  << std::endl;
		     }
		  }
		  catch (const std::runtime_error &rte) {
		     std::cout << "No bond for atom names: \""
			       << at_name_1 << "\" \""
			       << at_name_2 << "\": " << rte.what() << std::endl;
		  }
	       }
	    } else {
	       std::cout << "mismatch between number of atoms in molecule and "
			 << "COD types list" << std::endl;
	    } 
	 }
	 catch(const std::runtime_error &rte) {
	    std::cout << "error::" << rte.what() << std::endl;
	 }
      } else {
	 std::cout << "Mismatch in residue name vs dictionary comp_id "
		   << res_name << " " << rest.residue_info.comp_id
		   << std::endl;
      } 
   } else {
      std::cout << "Null res" << std::endl;
   } 

}


cod::bond_table_record_t
cod::bond_record_container_t::get_cod_bond_from_table(const std::string &cod_type_1,
						      const std::string &cod_type_2) const {

   std::map<std::string, std::map<std::string, bond_table_record_t> >::const_iterator it;

   bond_table_record_t bond;
   
   it = bonds_map.find(cod_type_1);
   if (it != bonds_map.end()) {
      std::map<std::string, bond_table_record_t>::const_iterator it_i;
      it_i = it->second.find(cod_type_2);
      if (it_i != it->second.end()) {
	 bond = it_i->second;
      } else {
	 std::string m("missing cod_type_2 " + cod_type_2);
	 throw(std::runtime_error(m));
      }
   } else {
      std::string m("missing cod_type_1 " + cod_type_1);
      throw(std::runtime_error(m));
   }
   return bond;

}

// This function may be called for bonds to Hydrogen atoms - they are in the dictionary
// but often not in the model.
//
double
cod::bond_record_container_t::get_bond_distance_from_model(const std::string &at_name_1,
							   const std::string &at_name_2,
							   mmdb::Residue *residue_p) const {

   double d = 0;
   mmdb::PPAtom residue_atoms = 0;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);

   mmdb::Atom *at_1 = NULL;
   mmdb::Atom *at_2 = NULL;

   if (n_residue_atoms) {
      for (int iat=0; iat<n_residue_atoms; iat++) { 
	 mmdb::Atom *at_p = residue_atoms[iat];
	 std::string atom_name = at_p->GetAtomName();
	 if (atom_name == at_name_1)
	    at_1 = at_p;
	 if (atom_name == at_name_2)
	    at_2 = at_p;
      }

      if (at_1 && at_2) {

	 clipper::Coord_orth pt_1 = coot::co(at_1);
	 clipper::Coord_orth pt_2 = coot::co(at_2);

	 d = clipper::Coord_orth::length(pt_1, pt_2);
      }
   }

   return d;

}
