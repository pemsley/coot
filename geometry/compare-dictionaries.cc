/* geometry/compare-dictionaries.cc
 * 
 * Copyright 2013 by Medical Research Council
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

#include "compat/coot-getopt.h"

#include "utils/coot-utils.hh"
#include "protein-geometry.hh"
// remember that coot-utils headers are not allowed here


mmdb::Manager *
create_mmdbmanager_from_residue(mmdb::Residue *res) {

   mmdb::Manager *mol = NULL;

   if (res) {
      mol = new mmdb::Manager;
      mmdb::Model *model_p = new mmdb::Model;
      mmdb::Chain *chain_p = new mmdb::Chain;
      chain_p->AddResidue(res);
      model_p->AddChain(chain_p);
      mol->AddModel(model_p);
      if (mol) {
	 chain_p->SetChainID(res->GetChainID());
      } else {
	 chain_p->SetChainID("");
      }
   }
   return mol;
}


int
compare_dictionaries(const std::string &type,
		     const std::string &file_name_1,
		     const std::string &file_name_2,
		     double bond_length_tolerance,
		     double bond_esd_tolerance,
		     double angle_tolerance,
		     double angle_esd_tolerance,
		     bool compare_hydrogens,
		     bool output_energy_types,
		     bool quiet) {

   int status = 0;
   // bool found_type_in_dictionary_1 = true; // not used
   bool found_type_in_dictionary_2 = true;

   coot::protein_geometry pg_1;
   coot::protein_geometry pg_2;

   pg_1.set_verbose(false);
   pg_2.set_verbose(false);

   pg_1.init_refmac_mon_lib(file_name_1, 0);
   pg_2.init_refmac_mon_lib(file_name_2, 0);

   std::pair<bool, coot::dictionary_residue_restraints_t> r1 = 
      pg_1.get_monomer_restraints(type, coot::protein_geometry::IMOL_ENC_ANY);
   std::pair<bool, coot::dictionary_residue_restraints_t> r2 = 
      pg_2.get_monomer_restraints(type, coot::protein_geometry::IMOL_ENC_ANY);


   if (!r1.first) { 
      std::cout << "Failed to find restraints for type " << type << " in "
		<< file_name_1 << std::endl;
      status = 1;
   } else {
      if (!r2.first) { 
	 std::cout << "Failed to find restraints for type " << type << " in "
		   << file_name_2 << std::endl;
	 found_type_in_dictionary_2 = false;
	 status = 1;
      } else {
	 // Happy path

	 if (output_energy_types) { 
	    std::cout << "Comparing dictionaries:: " << file_name_1 << " vs. " << file_name_2
		      << std::endl;
	    
	    // Needs to be done in a directory that knows about RDKit
	    // (which currently geometry is not).
	    // 
	    // RDKit::RWMol m = rdkit_mol(r1.second);
	    // std::vector<std::string> cod_name = get_cod_atom_types(m);
	 } 


	 //
	 // if compare_status is true, they matched.
	 bool compare_status = r1.second.compare(r2.second,
						 bond_length_tolerance,
						 bond_esd_tolerance,
						 angle_tolerance,
						 angle_esd_tolerance,
						 compare_hydrogens,
						 output_energy_types,
						 quiet);
	 status = !compare_status; // invert for unix return value (0 happy)
      }
   }

   if (! found_type_in_dictionary_2) {
      std::cout << "There are " << pg_1.size() << " types in pg_1" << std::endl;

      std::cout << "Trying Graph matching..." << std::endl;

      if (pg_2.size() == 1) {
	 std::string ref_type = pg_2[0].second.comp_id();
	 std::cout << "getting dictionary for " << ref_type << " from " << file_name_2 << std::endl;
	 r2 = pg_2.get_monomer_restraints(ref_type, coot::protein_geometry::IMOL_ENC_ANY);
	 if (r2.first) { // should be!
	    mmdb::Residue *residue_p = NULL; // for the moment.
	    std::string new_comp_id = type;
	    std::string new_compound_name = "some-new-name-here";

	    std::cout << "calling match_to_reference " << r2.second.comp_id() << " vs "
		      << r1.second.comp_id() << std::endl;
	    coot::dictionary_match_info_t dmi =
	       r2.second.match_to_reference(r1.second, residue_p,
					    new_comp_id, new_compound_name);

	    std::cout << "got " << dmi.n_matches << " matches and had "
		      << r1.second.number_of_non_hydrogen_atoms() << " non-H atoms" << std::endl;

	    if (dmi.n_matches == r1.second.number_of_non_hydrogen_atoms()) {
	       std::string cif_dict_out = "matching-dictionary.cif";
	       dmi.dict.write_cif(cif_dict_out);
	       bool get_idealized = true;
	       float b_factor = 30;
	       mmdb::Residue *residue_p = dmi.dict.GetResidue(get_idealized, b_factor);
	       mmdb::Manager *m = create_mmdbmanager_from_residue(residue_p);
	       m->WritePDBASCII("new-dict.pdb");
	    }
	 }
      }
   }
   
   return status;
}

void print_help(std::string cmd) {
   std::cout << "Usage: " << cmd << " "
      "        --help     this help\n" << 
      "        --quiet    do not report satisfatory matches\n" << 
      "        --type     residue type to match\n" << 
      "        --dict-1   file with reference dictionary\n"
      "        --dict-2   file with comparison dictionary\n" << 
      "        --bond-length-tol\n" << 
      "        --bond-length-esd-tol\n" << 
      "        --angle-tol\n" << 
      "        --angle-esd-tol\n" << 
      "        --full     output comp_ids, filename and energy types also\n" << 
      "        --include-hydrogens (default: geometry excludes hydrogens)\n" << 
      std::endl;
   
}



int main(int argc, char **argv) {

   int status = 0;
   bool quiet = false;
   bool compare_hydrogens = false;
   bool output_energy_types = false; // and comp_id

   if (argc < 4) {
      print_help(argv[0]);
   } else {

      std::string type;
      std::string file_name_1;
      std::string file_name_2;

      double bond_length_tolerance = 0.01;
      double bond_esd_tolerance    = 0.01;
      double angle_tolerance       = 1;
      double angle_esd_tolerance   = 1;
      
      const char *optstr = "q";
      struct option long_options[] = {
	 {"include-hydrogens",    0, 0, 0},
	 {"Hs",    0, 0, 0},
	 {"help",    0, 0, 0},
	 {"quiet",   0, 0, 0},
	 {"full",    0, 0, 0},
	 {"type",    1, 0, 0},
	 {"dict-1",  1, 0, 0},
	 {"dict-2",  1, 0, 0},
	 {"bond-length-tol",      1, 0, 0},
	 {"bond-length-esd-tol",  1, 0, 0},
	 {"angle-tol",            1, 0, 0},
	 {"angle-esd-tol",        1, 0, 0},
	 {0, 0, 0, 0}
      };

      int ch;
      int option_index = 0;
      while ( -1 != 
	      (ch = coot_getopt_long(argc, argv, optstr, long_options, &option_index))) {

	 switch(ch) {
	    
	 case 0:
	    if (coot_optarg) {
	       std::string arg_str = long_options[option_index].name;
	       if (arg_str == "type")
		  type = coot_optarg;
	       if (arg_str == "dict-1")
		  file_name_1 = coot_optarg;
	       if (arg_str == "dict-2")
		  file_name_2 = coot_optarg;
	       if (arg_str == "bond-length-tol") { 
		  try {
		     bond_length_tolerance = coot::util::string_to_float(coot_optarg);
		  }
		  catch (const std::runtime_error &rte) {
		     std::cout << "bad number " << coot_optarg << std::endl;
		  }
	       }
	       if (arg_str == "bond-length-esd-tol") { 
		  try {
		     bond_esd_tolerance = coot::util::string_to_float(coot_optarg);
		  }
		  catch (const std::runtime_error &rte) {
		     std::cout << "bad number " << coot_optarg << std::endl;
		  }
	       }
	       if (arg_str == "angle-tol") { 
		  try {
		     angle_tolerance = coot::util::string_to_float(coot_optarg);
		  }
		  catch (const std::runtime_error &rte) {
		     std::cout << "bad number " << coot_optarg << std::endl;
		  }
	       }
	       if (arg_str == "angle-esd-tol") { 
		  try {
		     angle_esd_tolerance = coot::util::string_to_float(coot_optarg);
		  }
		  catch (const std::runtime_error &rte) {
		     std::cout << "bad number " << coot_optarg << std::endl;
		  }
	       }

	       
	    } else {
	       std::string arg_str = long_options[option_index].name;
	       if (arg_str == "help") { 
		  print_help(argv[0]);
	       }
	       if (arg_str == "quiet") { 
		  quiet = true;
	       }
	       if (arg_str == "include-hydrogens") { 
		  compare_hydrogens = true;
	       }
	       if (arg_str == "Hs") { 
		  compare_hydrogens = true;
	       }
	       if (arg_str == "full") { 
		  output_energy_types = true;
	       }
	    }
	    break;

	 case 'q':
	    quiet = true;
	 }
      }

      if (type.empty()) {
	 std::cout << "missing type" << std::endl;
      } else { 
	 if (file_name_1.empty()) {
	    std::cout << "missing dict-1" << std::endl;
	 } else { 
	    if (file_name_2.empty()) {
	       std::cout << "missing dict-2" << std::endl;
	    } else {
	       // Happy path
	       status = compare_dictionaries(type, file_name_1, file_name_2,
					     bond_length_tolerance,
					     bond_esd_tolerance,
					     angle_tolerance,
					     angle_esd_tolerance,
					     compare_hydrogens,
					     output_energy_types,
					     quiet);
	    }
	 }
      }
   }
   return status;
}
