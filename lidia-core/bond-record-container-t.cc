/* lidia-core/bond-record-container-t.cc
 * 
 * Copyright 2016 by Medical Research Council
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

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

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

// #include "coot-utils/coot-coord-utils.hh" out of order now

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
      const RDKit::Atom *at_p = mol[iat];
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
      std::map<std::string, unsigned int>::const_iterator it_map = level_4_atom_types_map.begin();
      unsigned int n = level_4_atom_types_map.size();
      for (unsigned int i=0; i<n; i++) {
	 f << std::setw(6) << it_map->second << " " << it_map->first << "\n";
	 it_map++;
      }
   }
   f.close();

   return status;
}


bool
cod::bond_record_container_t::write(const std::string &atom_indices_file_name,
				    const std::string &file_name) const {

   write_atom_type_indices(atom_indices_file_name);

   bool status = false;
   std::ofstream f(file_name.c_str());

   if (f) {
      for (unsigned int i=0; i<bonds.size(); i++) {
	 const bond_table_record_t &btr = bonds[i];
	 // f << bonds[i] << "\n";
	 // btr.write(f, max_atom_type_width);
	 // dangerous, but should not fail
	 btr.write(f,
		   level_4_atom_types_map.find(btr.cod_type_1.level_4)->second,
		   level_4_atom_types_map.find(btr.cod_type_2.level_4)->second);
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
      std::string::size_type s1 = bonds[i].cod_type_1.level_4.length();
      std::string::size_type s2 = bonds[i].cod_type_2.level_4.length();

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
      int line_number = 0;

      while (std::getline(f, line)) {
	 line_number++;
	 std::vector<std::string> bits = coot::util::split_string_no_blanks(line);

	 if (bits.size() == 18) {
	    const std::string &cod_type_1_hash_str = bits[0];
	    const std::string &cod_type_2_hash_str = bits[1];

	    const std::string &hybridization_token = bits[2];
	    const std::string &ring_token = bits[3];

	    // 20170613 capture the degree:colon types too
	    const std::string &atom_1_colon_type  = bits[4];
	    const std::string &atom_2_colon_type  = bits[5];

	    // 20170612 extra-electron types now
	    const std::string &cod_type_1_level_2_str  = bits[6];
	    const std::string &cod_type_2_level_2_str  = bits[7];

	    atom_level_2_type cod_type_atom_1_l2(cod_type_1_level_2_str);
	    atom_level_2_type cod_type_atom_2_l2(cod_type_2_level_2_str);
	    
	    const std::string &cod_type_1_level_4  = bits[10];
	    const std::string &cod_type_2_level_4  = bits[11];

	    std::string cod_type_1_level_3 =
	       atom_type_t::level_4_type_to_level_3_type(cod_type_1_level_4);
	    std::string cod_type_2_level_3 =
	       atom_type_t::level_4_type_to_level_3_type(cod_type_2_level_4);

	    if (false) // debug
	       std::cout << "debug:: making an cod-atom from "
			 << cod_type_1_hash_str << " "
			 << atom_1_colon_type << "   "
			 << cod_type_atom_1_l2.string() << "    "
			 << cod_type_1_level_3 << "    "
			 << cod_type_1_level_4 << "\n";
	    atom_type_t cod_type_1(cod_type_1_hash_str, atom_1_colon_type,
				   cod_type_atom_1_l2,
				   cod_type_1_level_3,  cod_type_1_level_4);
	    atom_type_t cod_type_2(cod_type_2_hash_str, atom_2_colon_type,
				   cod_type_atom_2_l2,
				   cod_type_2_level_3,  cod_type_2_level_4);

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
	       rec.file_name = file_name;
	       rec.line_number = line_number;
	       rec.hybridization_token = hybridization_token;
	       rec.ring_token = ring_token;
	       add(rec);
	       status = true;
	    }
	    catch (const std::runtime_error &rte) {
	       std::cout << "error converting " << rte.what() << std::endl;
	    }
	 } else {
	    std::cout << "from " << file_name << " reject " << line << std::endl;
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
      std::string fn = coot::util::file_name_non_directory(table);
      if (fn != "bond_idx.table") { 
	 bond_record_container_t single;
	 bool success = single.read_acedrg_table(table);
	 if (success)
	    add_table(single);
      }
   }

   std::cout << "stored " << size() << " bond records" << std::endl;

   // why do I sort?
   std::cout << "-- pre-sort " << std::endl;
   sort();
   std::cout << "-- post-sort " << std::endl;

   std::cout << "--  pre-fill bonds map " << std::endl;
   fill_bonds_map();
   std::cout << "-- post-fill bonds map " << std::endl;
   std::cout << "--  pre-fill atoms map " << std::endl;
   // fill_level_4_atom_map(); // what does this do, I mean, why is the atom map needed?
   fill_cod_atom_type_map();
   std::cout << "-- post-fill atoms map " << std::endl;
}

void
cod::bond_record_container_t::fill_cod_atom_type_map() {

   cod_atom_type_set.clear();
   for (unsigned int i=0; i<bonds.size(); i++) {
      if (false)
	 std::cout << i << " inserting " << bonds[i].cod_type_1.level_4 << " " << bonds[i].cod_type_1.level_3
		   << " " << bonds[i].cod_type_1.level_2.string() << " " << bonds[i].cod_type_1.hash_value
		   << std::endl;
      cod_atom_type_set.insert(bonds[i].cod_type_1);
      cod_atom_type_set.insert(bonds[i].cod_type_2);
   }
   std::set<atom_type_t>::const_iterator it;
   int idx = 0;
   for (it=cod_atom_type_set.begin(); it!=cod_atom_type_set.end(); it++) {
      cod_atom_type_map.insert(std::pair<atom_type_t, int> (*it, idx));
      idx++;
   }
}

// old, just level_4 strings
void
cod::bond_record_container_t::fill_level_4_atom_map() {

   // now fill the atom map (with level_4 types)

   std::set<std::string> set;
   std::cout << "pre  fill atom set" << std::endl;
   for (unsigned int i=0; i<bonds.size(); i++) {
      set.insert(bonds[i].cod_type_1.level_4);
      set.insert(bonds[i].cod_type_2.level_4);
   }
   std::cout << "post fill atom set " << set.size() << std::endl;

   unsigned int n = set.size();
   std::set<std::string>::const_iterator it;
   unsigned int idx = 0;
   std::cout << "pre  fill atom map" << std::endl;
   for (it=set.begin(); it!=set.end(); it++) {
      // std::cout << "    " << *it << std::endl;
      level_4_atom_types_map[*it] = idx;
      idx++;
   }
   std::cout << "post fill atom map" << std::endl;

   std::map<std::string, unsigned int>::const_iterator it_map;
   it_map = level_4_atom_types_map.begin();
}

void
cod::bond_record_container_t::fill_bonds_map() {

   for (unsigned int i=0; i<bonds.size(); i++) {
      const atom_type_t &c1 = bonds[i].cod_type_1;
      const atom_type_t &c2 = bonds[i].cod_type_2;

      std::string c1l2 = c1.level_2.string();
      std::string c2l2 = c2.level_2.string();

      if (std::find(bonds_map[c1l2][c2l2][c1.level_3][c2.level_3].begin(),
		    bonds_map[c1l2][c2l2][c1.level_3][c2.level_3].end(),
		    bonds[i]) == bonds_map[c1l2][c2l2][c1.level_3][c2.level_3].end())
	 bonds_map[c1l2][c2l2][c1.level_3][c2.level_3].push_back(bonds[i]);
      else
	 std::cout << "WARNING " << bonds[i] << " was already in the bonds_map" << std::endl;

      // we don't need reverse lookup for bonds between atoms of the same type.
      if (! (c1 == c2)) {
	 bond_table_record_t reverse_bond = bonds[i];
	 std::swap(reverse_bond.cod_type_1, reverse_bond.cod_type_2);
	 if (std::find(bonds_map[c2l2][c1l2][c2.level_3][c1.level_3].begin(),
		       bonds_map[c2l2][c1l2][c2.level_3][c1.level_3].end(),
		       reverse_bond) == bonds_map[c2l2][c1l2][c2.level_3][c1.level_3].end())
	    bonds_map[c2l2][c1l2][c2.level_3][c1.level_3].push_back(reverse_bond);
	 else
	    std::cout << "WARNING " << reverse_bond << " was already in the bonds_map" << std::endl;
      }
   }
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
	       double std_dev = coot::util::string_to_double(line.substr(12,21));
	       int count      = coot::util::string_to_int(line.substr(22, 7));
	       int idx_1      = coot::util::string_to_int(line.substr(29, 7));
	       int idx_2      = coot::util::string_to_int(line.substr(37, 7));

	       std::string type_1_level_4 = atom_types[idx_1];
	       std::string type_2_level_4 = atom_types[idx_2];

	       atom_type_t type_1(type_1_level_4); // hash and level-2 are empty
	       atom_type_t type_2(type_2_level_4);

	       if (false)
		  std::cout << "parsed "
			    << std::setw(10) << bl << " "
			    << std::setw(10) << std_dev << " "
			    << std::setw( 8) << count
			    << " idx_1: " << idx_1 << "    idx_2: " << idx_2 << "   "
			    << type_1.level_4 << "   " << type_2.level_4 << "\n";

	       if (type_2 < type_1)
		  std::swap(type_1, type_2);
	       
	       bond_table_record_t bond(type_1, type_2, bl, std_dev, count);

	       std::string l3_type_1 = type_1.level_3;
	       std::string l3_type_2 = type_2.level_3;

	       // this function (read_bonds()) doesn't read the level-2 types
	       // bonds_map[l3_type_1][l3_type_2].push_back(bond);
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
	    if (bonds[i].cod_type_2 < bonds[i].cod_type_1)
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
   
   std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, std::vector<bond_table_record_t> > > > >::const_iterator it_l2_a;
   std::map<std::string, std::map<std::string, std::map<std::string, std::vector<bond_table_record_t> > > >::const_iterator it_l2_b;
   
   for (it_l2_a=bonds_map.begin(); it_l2_a!=bonds_map.end(); it_l2_a++) {
      std::string key_l2_a = it_l2_a->first;

      for (it_l2_b=it_l2_a->second.begin(); it_l2_b!=it_l2_a->second.end(); it_l2_b++) {
	 std::string key_l2_b = it_l2_b->first;
	 std::cout << "   l2 types " << key_l2_a << " " << key_l2_b << std::endl;
      }
   }
   std::cout << " done check()" << std::endl;
}

std::vector<bool>
cod::bond_record_container_t::get_is_hydrogen_flags(const RDKit::RWMol &rdkm) const {

   unsigned int n_mol_atoms = rdkm.getNumAtoms();
   std::vector<bool> flags(n_mol_atoms);
   for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
      const RDKit::Atom *at_p = rdkm[iat];
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
	    coot::rdkit_mol_sanitize(rdkm);

	    // coot::debug_rdkit_molecule(&rdkm);

	    atom_types_t t;
	    std::vector<cod::atom_type_t> v = t.get_cod_atom_types(rdkm, true);

	    if (rdkm.getNumAtoms() == v.size()) {
	       unsigned int n_mol_atoms = rdkm.getNumAtoms();

	       std::cout << "---- validate() types table ----- " << std::endl;
	       for (unsigned int iat=0; iat<n_mol_atoms; iat++) {
		  RDKit::Atom *at_p = rdkm[iat];
		  try {
		     std::string name;
		     at_p->getProp("name", name);
		     std::cout << "    " << iat << " " << name
			       << "  " << v[iat].hash_value
			       << "  \"" << v[iat].level_2.string() << "\""
			       << "  \"" << v[iat].level_3 << "\""
			       << "  \"" << v[iat].level_4 << "\"" << std::endl;
		  }
		  catch (const KeyErrorException &kee) {
		     std::cout << "   " << iat << " " << kee.what()
			       << std::endl;
		  }
	       }
	       std::cout << "------------------- " << std::endl;
	    }


	    // check(); // debug

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

			if (true)
			   std::cout << "Bond for atoms: " << at_name_1
				     << " " << at_name_2 << std::endl;

			atom_type_t cod_type_1 = v[atom_idx_1];
			atom_type_t cod_type_2 = v[atom_idx_2];

			if (cod_type_2.level_4 < cod_type_1.level_4)
			   std::swap(cod_type_1, cod_type_2);
		     
			bond_table_record_t cod_bond = 
			   get_cod_bond_from_table(cod_type_1, cod_type_2);
		     
			double dist =
			   get_bond_distance_from_model(at_name_1_4c, at_name_2_4c, res);
			
			double z = 9999;
			if (cod_bond.std_dev > 0)
			   z = std::abs((cod_bond.mean-dist)/cod_bond.std_dev);

			std::cout << "  compare: " << std::setw(4) << at_name_1 << " "
				  << std::setw(4) << at_name_2
				  << std::setw(8) << " model: "
				  << std::setw(8) << dist
				  << " vs tables: "
				  << std::setw(8) << cod_bond.mean << " +/- "
				  << std::setw(6) << cod_bond.std_dev
				  << " counts: " << cod_bond.count
				  << "    z = "
				  << std::setw(8) << z
				  << std::endl;

			   // << cod_type_1.level_3 << "  " << cod_type_1.level_4 << " "
			   // << cod_type_2.level_3 << "  " << cod_type_2.level_4
		     }
		  }
		  catch (const std::runtime_error &rte) {
		     std::cout << "   No bond for atom names: \""
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
cod::bond_record_container_t::get_cod_bond_from_table(const cod::atom_type_t &cod_type_1,
						      const cod::atom_type_t &cod_type_2) const {

   cod::bond_table_record_t bond;
   bool found_bond = false;

   if (true) {
      std::cout << "  get_cod_bond_from_table() using "
		<< cod_type_1.level_2.string() << "   "
		<< cod_type_1.level_3 << std::endl;
      std::cout << "                                  "
		<< cod_type_2.level_2.string() << "   "
		<< cod_type_2.level_3 << std::endl;
   }
      
   std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, std::vector<bond_table_record_t> > > > >::const_iterator it_l2_a;
   std::map<std::string, std::map<std::string, std::map<std::string, std::vector<bond_table_record_t> > > >::const_iterator it_l2_b;
   
   std::map<std::string, std::map<std::string, std::vector<bond_table_record_t> > >::const_iterator it_l3_a;

   std::map<std::string, std::vector<bond_table_record_t> >::const_iterator it_l3_b;
   
   it_l2_a = bonds_map.find(cod_type_1.level_2.string());
   if (it_l2_a != bonds_map.end()) {
      it_l2_b = it_l2_a->second.find(cod_type_2.level_2.string());
      if (it_l2_b != it_l2_a->second.end()) {
	 it_l3_a = it_l2_b->second.find(cod_type_1.level_3);
	 if (it_l3_a != it_l2_b->second.end()) {
	    it_l3_b = it_l3_a->second.find(cod_type_2.level_3);
	    if (it_l3_b != it_l3_a->second.end()) {

	       // OK, normal case, 2 level-3 hits
	       // unsigned int approx_level = 0;
	       bond_table_record_t::approximation_level_t approx_level = bond_table_record_t::FULL_ATOM;
	       bond = make_bond_from_level_3_vector(cod_type_1, cod_type_2, it_l3_b->second,
						    approx_level);
	       found_bond = true;

	    } else {

	       // we need to accumulate a level_3 vector from level_2s.
	       //
	       // unsigned int approx_level = 1;
	       bond_table_record_t::approximation_level_t approx_level = bond_table_record_t::EXTRA_ELECTRON;
	       bond = make_bond_from_level_2_map(cod_type_1, cod_type_2, it_l2_b->second, approx_level);
	       found_bond = true;
	       
	       // std::string m("missing cod_type_2 level_3 " + cod_type_1.level_2.string());
	       // throw(std::runtime_error(m));
	    }
	 } else {

	    // unsigned int approx_level = 1;
	    bond_table_record_t::approximation_level_t approx_level = bond_table_record_t::MAIN_SECTION;
	    bond = make_bond_from_level_2_map(cod_type_1, cod_type_2, it_l2_b->second, approx_level);
	    found_bond = true;

	    // t3_miss_diagnose(cod_type_1, cod_type_2); // debugging
	    // std::string m("missing cod_type_1 level_3 " + cod_type_1.level_3);
	    // throw(std::runtime_error(m));
	 }
      } else {
	 std::string m("missing cod_type_2 level_2 " + cod_type_2.level_2.string());
	 throw(std::runtime_error(m));
      }
   } else {
      std::string m("missing cod_type_1 level_2 " + cod_type_1.level_2.string());
      throw(std::runtime_error(m));
   }
   return bond;
}

void
cod::bond_record_container_t::add_table(const bond_record_container_t &brc) {
   for (unsigned int i=0; i<brc.bonds.size(); i++) {
      std::vector<bond_table_record_t>::const_iterator it;

      if (false) { // debugging.
	 // they are unique
	 it = std::find(bonds.begin(), bonds.end(), brc.bonds[i]);
	 if (it == bonds.end()) {
	    bonds.push_back(brc.bonds[i]);
	 } else {
	    std::cout << "Not adding " << brc.bonds[i] << " --already-exists-- "
		      << *it << std::endl;
	 }
      }

      // real
      bonds.push_back(brc.bonds[i]);

   }
}


// select or average (consolidate)
//
cod::bond_table_record_t
cod::bond_record_container_t::make_bond_from_level_3_vector(const cod::atom_type_t &cod_type_1,
							    const cod::atom_type_t &cod_type_2,
							    const std::vector<cod::bond_table_record_t> &v,
							    bond_table_record_t::approximation_level_t approx_level) const {

   cod::bond_table_record_t b = v[0];
   unsigned int min_count_for_l4_match = 4; // counts should be at least this value
                                            // Fei's value is 4, I used to have 6 (guess).

   // when we are consolidating bonds, what should be the sum of counts so that we
   // think that we've properly constructed a bond.
   // 
   unsigned int min_count_sum_for_merge = 4;

   if (v.size() > 1) {

      bool found_bond = false;
      // first, is there an exact level_4 match?
      //
      for (unsigned int i=0; i<v.size(); i++) { 
	 if (cod_type_1 == v[i].cod_type_1) {
	    if (cod_type_2 == v[i].cod_type_2) {
	       if (v[1].count >= min_count_for_l4_match) {
		  b = v[i];
		  found_bond = true;
		  break;
	       }
	    }
	 }
      }

      if (! found_bond) {

	 // local bonds: collect a set of bond records that have exact
	 // match for cod_type_1
	 // 
	 std::vector<bond_table_record_t> local_bonds;
	 for (unsigned int i=0; i<v.size(); i++) { 
	    if (cod_type_1 == v[i].cod_type_1) {
	       // do we need a "distance" metric here before we add this?
	       local_bonds.push_back(v[i]);
	    }
	 }
	 unsigned int n_sum = 0;
	 for (unsigned int j=0; j<local_bonds.size(); j++)
	    n_sum += local_bonds[j].count;
	 if (n_sum > min_count_sum_for_merge) {
	    b = consolidate_bonds(cod_type_1, cod_type_2, local_bonds, approx_level);
	    found_bond = true;
	 }
      }

      if (! found_bond) {

	 // likewise but the other testing the other atom
	 std::vector<bond_table_record_t> local_bonds;
	 for (unsigned int i=0; i<v.size(); i++) { 
	    if (cod_type_2 == v[i].cod_type_2) {
	       local_bonds.push_back(v[i]);
	    }
	 }
	 unsigned int n_sum = 0;
	 for (unsigned int j=0; j<local_bonds.size(); j++)
	    n_sum += local_bonds[j].count;
	 if (n_sum > min_count_sum_for_merge) {
	    b = consolidate_bonds(cod_type_1, cod_type_2, local_bonds, approx_level);
	    found_bond = true;
	 }
      }

      if (! found_bond)
	 b = consolidate_bonds(cod_type_1, cod_type_2, v, approx_level);
   }
   return b;
}

// generalization
//
cod::bond_table_record_t
cod::bond_record_container_t::make_bond_from_level_2_map(const atom_type_t &cod_type_1,
							 const atom_type_t &cod_type_2,
							 const std::map<std::string, std::map<std::string, std::vector<bond_table_record_t> > > &l3_map,
							 cod::bond_table_record_t::approximation_level_t approx_level) const {

   bond_table_record_t bond;
   std::vector<bond_table_record_t> v;
   std::map<std::string, std::map<std::string, std::vector<bond_table_record_t> > >::const_iterator it_l3_a_local;
   std::map<std::string, std::vector<bond_table_record_t> >::const_iterator it_l3_b_local;

   for (it_l3_a_local=l3_map.begin();
	it_l3_a_local != l3_map.end();
	it_l3_a_local++) {
      for (it_l3_b_local=it_l3_a_local->second.begin();
	   it_l3_b_local != it_l3_a_local->second.end();
	   it_l3_b_local++) {
	 for (unsigned int i=0; i<it_l3_b_local->second.size(); i++) {
	    const bond_table_record_t &bond = it_l3_b_local->second[i];
	    v.push_back(bond);
	 }
      }
   }

   bond = make_bond_from_level_3_vector(cod_type_1, cod_type_2, v, approx_level);

   return bond;

}


cod::bond_table_record_t
cod::bond_record_container_t::consolidate_bonds(const cod::atom_type_t &cod_type_1,
						const cod::atom_type_t &cod_type_2,
						const std::vector<cod::bond_table_record_t> &lb,
						bond_table_record_t::approximation_level_t approx_level) const {

   double mean_sum = 0;
   double var_sum = 0;
   unsigned int n_sum = 0;
   for (unsigned int j=0; j<lb.size(); j++) {
      const double &sd = lb[j].std_dev;
      mean_sum += lb[j].mean * lb[j].count;
      var_sum  += sd * sd * lb[j].count;
      n_sum += lb[j].count;
   }

   double mean = mean_sum/double(n_sum);
   double var  = var_sum/double(n_sum);
   double sd = sqrt(var);
   return bond_table_record_t(cod_type_1, cod_type_2, mean, sd, n_sum, approx_level);
}
	       

void
cod::bond_record_container_t::t3_miss_diagnose(const atom_type_t &cod_type_1,
					       const atom_type_t &cod_type_2) const {

   std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, std::vector<bond_table_record_t> > > > >::const_iterator it_l2_a;
   std::map<std::string, std::map<std::string, std::map<std::string, std::vector<bond_table_record_t> > > >::const_iterator it_l2_b;
   
   std::map<std::string, std::map<std::string, std::vector<bond_table_record_t> > >::const_iterator it_l3_a;

   std::map<std::string, std::vector<bond_table_record_t> >::const_iterator it_l3_b;

   bond_table_record_t bond_found;
   bool found_bond = false;

   bool vector_found = false;
   
   for (it_l2_a=bonds_map.begin(); it_l2_a!=bonds_map.end(); it_l2_a++) {
      for (it_l2_b=it_l2_a->second.begin(); it_l2_b!=it_l2_a->second.end(); it_l2_b++) {
	 for (it_l3_a=it_l2_b->second.begin(); it_l3_a != it_l2_b->second.end(); it_l3_a++) {
	    for (it_l3_b=it_l3_a->second.begin(); it_l3_b != it_l3_a->second.end(); it_l3_b++) {
	       for (unsigned int i=0; i<it_l3_b->second.size(); i++) {
		  const bond_table_record_t &bond = it_l3_b->second[i];
		  if (bond.cod_type_1.l3_match(cod_type_1)) {

		     std::cout << "found-1 " << bond.cod_type_1.level_3 << std::endl;
		     if (bond.cod_type_2.l3_match(cod_type_2)) {
			std::cout << "found-2 " << cod_type_2.level_3 << std::endl;
			bond_found = bond;
			found_bond = true;
			break;
		     } else {
			std::cout << "  test for cod_type_2 " << bond.cod_type_2.level_3
				  << " is not \n                      " << cod_type_2.level_3
				  << std::endl;
		     }
		  }
		  if (bond.cod_type_1.l3_match(cod_type_1)) {
		     if (bond.cod_type_2.l3_match(cod_type_2)) {
			vector_found = true;
		     }
		  }

		  // reverse indexing

		  if (bond.cod_type_1.l3_match(cod_type_2)) {
		     std::cout << "Reverse Found-1 " << bond.cod_type_1.level_3 << std::endl;
		     if (bond.cod_type_2.l3_match(cod_type_1)) {
			std::cout << "Reverse Found-2 " << bond.cod_type_2.level_3
				  << "  !!! ::::::::::::::::::: " << std::endl;
		     } else {
			std::cout << "  reverse test for cod_type_2 " << bond.cod_type_2.level_3
				  << " is not \n                              "
				  << cod_type_1.level_3 << std::endl;
		     }
		  }
		  
	       }
	       if (found_bond)
		  break;
	    }
	    if (found_bond)
	       break;
	 }
	 if (found_bond)
	    break;
      }
      if (found_bond)
	 break;
   }

   if (found_bond) {
      std::cout << "::::::: Hmmmmm t3_miss_diagnose() found bond "
		<< std::endl;
   } else {
      std::cout << "::::::: Hmmmmm t3_miss_diagnose() bond not found "
		<< std::endl;
   }
   
   if (vector_found) {
      std::cout << "::::::: Hmmmmm t3_miss_diagnose() found vector "
		<< std::endl;
   } else {
      std::cout << "::::::: Hmmmmm t3_miss_diagnose() vector not found "
		<< std::endl;
   }
   
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

	 clipper::Coord_orth pt_1 = co(at_1);
	 clipper::Coord_orth pt_2 = co(at_2);

	 d = clipper::Coord_orth::length(pt_1, pt_2);
      }
   }

   return d;

}

#ifdef USE_SQLITE3

// create the database with tables (but not data)
//
sqlite3 *
cod::bond_record_container_t::make_sqlite_db(const std::string &db_file_name) {

   sqlite3 *db = NULL;
   if (! coot::file_exists(db_file_name)) {
      int rc = sqlite3_open(db_file_name.c_str(), &db);
      char *zErrMsg = 0;
      std::string command;

      if (rc == 0) {
	 // happy path
	 std::vector<std::string> commands;
	 std::string create_cmd =
	    std::string("CREATE TABLE COD_TYPE_4_INDICES (") +
	    std::string("level_4_atom_type NVCHAR(200) PRIMARY KEY, level_3_atom_type NVCHAR(200), level_2_atom_type NVCHAR(100), colon_degree_atom_type NVCHAR(12), hash_code INT, atom_index INT )");

	 rc = sqlite3_exec(db, create_cmd.c_str(), db_callback, 0, &zErrMsg);
	 if (rc) {
	    if (zErrMsg) {
	       std::cout << "rc for " << create_cmd << " " << rc << " " << zErrMsg << std::endl;
	    } else {
	       std::cout << "rc for " << create_cmd << " " << rc << " " << std::endl;
	    }
	    db = NULL; // something bad happend
	 } else {
	    // I think the ring info is a boolean (not nvchar(10)).
	    std::cout << "create table COD_TYPE_4_INDICES success" << std::endl;
	    create_cmd = "CREATE TABLE COD_TYPE_4_BONDS ";
	    create_cmd += "(atom_index_1 INT, atom_index_2 INT, hybridization NVCHAR(20), ring NVCHAR(10), mean REAL, sd REAL, count INT, PRIMARY KEY (atom_index_1, atom_index_2, hybridization, ring) )";
	    rc = sqlite3_exec(db, create_cmd.c_str(), db_callback, 0, &zErrMsg);
	    if (rc) {
	       if (zErrMsg) {
		  std::cout << "rc for " << create_cmd << " " << rc << " " << zErrMsg << std::endl;
	       } else {
		  std::cout << "rc for " << create_cmd << " " << rc << " " << std::endl;
	       }
	       db = NULL; // something bad happend
	    } else {
	       std::cout << "create table COD_TYPE_4_BONDS success" << std::endl;
	    }
	 }
      }
   } else {
      std::cout << "WARNING:: database file " << db_file_name << " already exists" << std::endl;
   }
   return db;
}
#endif // USE_SQLITE3

void
cod::bond_record_container_t::make_db(const std::string &db_file_name) {

#ifdef USE_SQLITE3
   sqlite3 *db = NULL;
   if (! coot::file_exists(db_file_name))
      db = make_sqlite_db(db_file_name); // database empty but with ligands table

   if (db) {
      char *zErrMsg = 0;
      sqlite3_exec(db, "BEGIN", db_callback, 0, &zErrMsg);
      // fill level 4 types
      db_add_level_4_types(db);
      
      sqlite3_exec(db, "END", db_callback, 0, &zErrMsg);
   } else {
      std::cout << "WARNING:: empty db in make_db() " << std::endl;
   }
#endif // USE_SQLITE3
}

#ifdef USE_SQLITE3
// static
int
cod::bond_record_container_t::db_callback(void *NotUsed, int argc, char **argv, char **azColName) {
   int status = 0;
   for (int i=0; i<argc; i++) {
      std::cout << " db_callback(): " << i << " " << argv[1];
   }
   std::cout << std::endl;
   return status;
}
#endif // USE_SQLITE3

#ifdef USE_SQLITE3
bool
cod::bond_record_container_t::db_add_level_4_types(sqlite3 *db) {

   bool status = true;
   std::string table_name = "COD_TYPE_4_INDICES";

   // std::map<std::string, unsigned int>::const_iterator it_map = level_4_atom_types_map.begin();
   // unsigned int n = level_4_atom_types_map.size();
   // std::cout << "in db_add_level_4_types() adding " << n << " atom types" << std::endl;

   std::map<atom_type_t, int>::const_iterator it_map = cod_atom_type_map.begin();
   unsigned int n = cod_atom_type_map.size();

   for (unsigned int i=0; i<n; i++) {
      std::string level_4 = it_map->first.level_4;
      std::string level_3 = it_map->first.level_3;
      std::string level_2 = it_map->first.level_2.string();
      std::string degree_colon_type = it_map->first.neighb_degrees_str();
      int hash_code = it_map->first.hash_value;
      std::string cmd = "INSERT INTO " + table_name + " ";
      cmd += "(level_4_atom_type, level_3_atom_type, level_2_atom_type, colon_degree_atom_type, hash_code, atom_index) ";
      cmd += "VALUES (";
      cmd += coot::util::single_quote(level_4, "'");
      cmd += ", ";
      cmd += coot::util::single_quote(level_3, "'");
      cmd += ", ";
      cmd += coot::util::single_quote(level_2, "'");
      cmd += ", ";
      cmd += coot::util::single_quote(degree_colon_type, "'");
      cmd += ", ";
      cmd += coot::util::int_to_string(hash_code);
      cmd += ", ";
      cmd += coot::util::int_to_string(it_map->second);
      cmd += ");";
      it_map++;
      char *zErrMsg = 0;
      int rc = sqlite3_exec(db, cmd.c_str(), db_callback, 0, &zErrMsg);
      if (rc !=  SQLITE_OK) {
	 if (zErrMsg) {
	    std::cout << "ERROR: processing command: " << cmd << " " << zErrMsg << std::endl;
	 } else {
	    std::cout << "ERROR when processing command: " << cmd << std::endl;
	    sqlite3_free(zErrMsg);
	 }
      } else {
	 // was OK!
	 status = true;
	 // std::cout << "cmd " << cmd << " ran OK " << std::endl;
      }
   }

   // ----------------------- bonds (indexed by atom type index) ----------------------

   table_name = "COD_TYPE_4_BONDS";

   std::map<std::string, std::map<std::string, std::map<std::string, std::map<std::string, std::vector<bond_table_record_t> > > > >::const_iterator it_1;
   std::map<std::string, std::map<std::string, std::map<std::string, std::vector<bond_table_record_t> > > >::const_iterator it_2;
   std::map<std::string, std::map<std::string, std::vector<bond_table_record_t> > >::const_iterator it_3;
   std::map<std::string, std::vector<bond_table_record_t> >::const_iterator it_4;
   std::vector<bond_table_record_t>::const_iterator it_5;

   unsigned int bond_index = 0;
   if (false) {
      for (it_1 = bonds_map.begin(); it_1 != bonds_map.end(); it_1++) {
	 for (it_2 = it_1->second.begin(); it_2 != it_1->second.end(); it_2++) {
	    for (it_3 = it_2->second.begin(); it_3 != it_2->second.end(); it_3++) {
	       for (it_4 = it_3->second.begin(); it_4 != it_3->second.end(); it_4++) {
		  for (it_5 = it_4->second.begin(); it_5 != it_4->second.end(); it_5++) {
		     const bond_table_record_t &btr = *it_5;
		     unsigned int index_1 = cod_atom_type_map[btr.cod_type_1];
		     unsigned int index_2 = cod_atom_type_map[btr.cod_type_2];
		     std::cout << ":::bond "
			       << btr.cod_type_1.level_4 << " "
			       << btr.cod_type_1.level_3 << " "
			       << btr.cod_type_1.level_2.string() << " "
			       << btr.cod_type_1.hash_value
			       << " index " << index_1 << " "
			       << btr.cod_type_2.level_4 << " "
			       << btr.cod_type_2.level_3 << " "
			       << btr.cod_type_2.level_2.string() << " "
			       << btr.cod_type_2.hash_value
			       << " index " << index_2
			       << "\n";
		  }
	       }
	    }
	 }
      }
   }
   // std::cout << "restarting it_1 from begin " << bonds_map.size() << std::endl;
   unsigned int it_1_count = 0;
   for (it_1 = bonds_map.begin(); it_1 != bonds_map.end(); it_1++) {
      it_1_count++;
      // std::cout << "restarting it_2 from begin " << it_1->second.size() << std::endl;
      unsigned int it_2_count = 0;
      for (it_2 = it_1->second.begin(); it_2 != it_1->second.end(); it_2++) {
	 it_2_count++;
	 // std::cout << "restarting it_3 from begin " << it_2->second.size() << std::endl;
	 unsigned int it_3_count = 0;
	 for (it_3 = it_2->second.begin(); it_3 != it_2->second.end(); it_3++) {
	    it_3_count++;
	    // std::cout << "restarting it_4 from begin " << it_3->second.size() << std::endl;
	    unsigned int it_4_count = 0;
	    for (it_4 = it_3->second.begin(); it_4 != it_3->second.end(); it_4++) {
	       it_4_count++;
	       // std::cout << "restarting it_5 from begin " << it_4->second.size() << std::endl;
	       unsigned int it_5_count = 0;
	       for (it_5 = it_4->second.begin(); it_5 != it_4->second.end(); it_5++) {
		  it_5_count++;

		  const bond_table_record_t &btr = *it_5;
		  unsigned int index_1 = cod_atom_type_map[btr.cod_type_1];
		  unsigned int index_2 = cod_atom_type_map[btr.cod_type_2];
		  if (false)
		     std::cout << "insert from bond map " << index_1 << " " << index_2 << " "
			       << btr << " counts " << it_1_count << " " << it_2_count << " "
			       << it_3_count << " " << it_4_count << " " << it_5_count << " "
			       << std::endl;
		  std::string cmd = "INSERT INTO " + table_name + " ";
		  cmd += ("(atom_index_1, atom_index_2, hybridization, ring, mean, sd, count)");
		  cmd += " VALUES (";
		  cmd += coot::util::int_to_string(index_1);
		  cmd += ", ";
		  cmd += coot::util::int_to_string(index_2);
		  cmd += ", ";
		  cmd += coot::util::single_quote(btr.hybridization_token, "'");
		  cmd += ", ";
		  cmd += coot::util::single_quote(btr.ring_token, "'");
		  cmd += ", ";
		  cmd += coot::util::float_to_string_using_dec_pl(btr.mean, 6);
		  cmd += ", ";
		  cmd += coot::util::float_to_string_using_dec_pl(btr.std_dev, 6);
		  cmd += ", ";
		  cmd += coot::util::int_to_string(btr.count);
		  cmd += ");";
      
		  char *zErrMsg = 0;
		  int rc = sqlite3_exec(db, cmd.c_str(), db_callback, 0, &zErrMsg);
		  if (rc !=  SQLITE_OK) {
		     if (zErrMsg) {
			std::cout << "ERROR: processing command: " << cmd << " " << zErrMsg
				  << " for "
				  << btr.cod_type_1.level_4 << " "
				  << btr.cod_type_2.level_4
				  << btr.hybridization_token << " "
				  << btr.ring_token << " from "
				  << btr.file_name << " line " << btr.line_number
				  << " counts "
				  << it_1_count << " "
				  << it_2_count << " "
				  << it_3_count << " "
				  << it_4_count << " "
				  << it_5_count << std::endl;
		     } else {
			std::cout << "ERROR when processing command: " << cmd << std::endl;
			sqlite3_free(zErrMsg);
		     }
		  } else {
		     // was OK!
		     status = true;
		     if (false)
			std::cout << "insert of            " << index_1 << " " << index_2 << " "
				  << btr << " was OK" << std::endl;
		  }

		  bond_index++;
	       }
	    }
	 }
      }
   }

   return status;
}
#endif // USE_SQLITE3


#endif // MAKE_ENHANCED_LIGAND_TOOLS
