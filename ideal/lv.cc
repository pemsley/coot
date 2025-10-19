/*
 * ideal/lv.cc
 * 
 * Copyright 2021 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 */


#include "compat/coot-getopt.h"
#include "utils/coot-utils.hh"
#include "simple-restraint.hh"

void show_usage() {

   std::string prog_name = "ligand-validation";

   std::cout << "Usage: " << prog_name << "\n"
             << "       --pdbin pdb-in-filename\n"
             << "       --chain-id ligand-chain-id\n"
             << "       --res-no ligand-residue-number\n"
             << "       --dictionary cif-dictionary-file-name\n"
             << "       --pocket       # include interactions to the binding pocket\n"
             << std::endl;
}

void print_version() {

   std::cout << "coot-ligand-validation version " << VERSION << std::endl;

}

class input_data_t {
public:
   input_data_t() {
      is_good = false;
      use_pocket = false;
   }
   bool is_good;
   bool use_pocket; // i.e. use NBC from ligand neighbours
   std::string pdb_file_name;
   std::string dictionary;
   std::string chain_id;
   std::string res_no_ligand_str;
   std::string res_no_protein_str;
};

input_data_t get_input_details(int argc, char **argv) {

   input_data_t d;

   int ch;
   int option_index = 0;

   const char *optstr = "i:h:f:p:o:m:1:2:c:w";

   struct option long_options[] = {
      {"version",  0, 0, 0},
      {"help",     0, 0, 0},
      {"pocket",   0, 0, 0}, // use protein NBCs in the analysis
      {"pdbin",    1, 0, 0},
      {"res-no",   1, 0, 0},
      {"chain-id", 1, 0, 0},
      {"dictin",   1, 0, 0},
      {"dictionary", 1, 0, 0},
   };

   while (-1 !=
          (ch = coot_getopt_long(argc, argv, optstr, long_options, &option_index))) {

      // std::cout << "debug::: ch: " << ch << std::endl;

      switch (ch) {
      case 0:
         // if (! coot_optarg) std::cout << "No coot_optarg" << std::endl;
         // std::cout << "coot_optarg " << coot_optarg << std::endl;
         if (coot_optarg) {
            std::string arg_str = long_options[option_index].name;
            if (arg_str == "pdbin") {
               d.pdb_file_name = coot_optarg;
            }
            if (arg_str == "chain-id") {
               d.chain_id = coot_optarg;
            }
            if (arg_str == "res-no") {
               d.res_no_ligand_str = coot_optarg;
            }
            if (arg_str == "res-no-protein") {
               d.res_no_protein_str = coot_optarg;
            }
            if (arg_str == "dictionary") {
               d.dictionary = coot_optarg;
            }
            if (arg_str == "dictin") {
               d.dictionary = coot_optarg;
            }
         } else {
            // long argument without parameter:
            std::string arg_str(long_options[option_index].name);
            if (arg_str == "help") {
               show_usage();
               exit(0);
            }
            if (arg_str == "version") {
               print_version();
               exit(0);
            }
            if (arg_str == "pocket") {
               d.use_pocket = true;
            }
         }
         break;

      case 'h':
         { show_usage(); exit(0); }
         break;

      case 'v':
         { print_version(); exit(0); }
         break;
      }
   } // end while

   if (! d.dictionary.empty())
      if (! d.pdb_file_name.empty())
         if (! d.res_no_ligand_str.empty())
            d.is_good = true;

   return d;
}

#include "coot-utils/atom-overlaps.hh"

void
validate_ligand(const std::string &pdb_file_name,
                mmdb::Residue *residue_p, mmdb::Manager *mol, bool include_environment_contacts_flag,
                coot::protein_geometry  &geom) {

   coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS;
   bool do_residue_internal_torsions = true;
   bool do_trans_peptide_restraints = false;
   coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
   int imol = 0; // dummy
   std::vector<std::pair<bool,mmdb::Residue *> > residues;
   residues.push_back(std::pair<bool, mmdb::Residue *>(false, residue_p));

   std::cout << "INFO:: include_environment_contacts_flag " << include_environment_contacts_flag << std::endl;
   if (include_environment_contacts_flag) {
   } else {
      // We need to create a molecule (and replace mol) that is just the ligand then.
      std::vector<mmdb::Residue *> res_vec;
      res_vec.push_back(residue_p);
      std::pair<bool, mmdb::Manager *> mol_new = coot::util::create_mmdbmanager_from_residue_vector(res_vec, mol);
      if (mol_new.first) {
         mol = mol_new.second;
      } else {
         std::cout << "ERROR:: Failure to isolate ligand " << std::endl;
         return;
      }
   }

   std::string monomer_type = residue_p->GetResName();
   int res_no = residue_p->GetSeqNum();
   const std::string chain_id = residue_p->GetChainID();
   std::vector<coot::atom_spec_t> fixed_atom_specs;
   std::vector<mmdb::Link> links;
   coot::restraints_container_t restraints(residues, links, geom, mol, fixed_atom_specs, 0);

   int n_threads = 1;
   ctpl::thread_pool tp(n_threads);
   restraints.thread_pool(&tp, n_threads);
   restraints.set_quiet_reporting();
   int nrestraints = restraints.make_restraints(imol,
                                                geom, flags,
                                                do_residue_internal_torsions,
                                                do_trans_peptide_restraints,
                                                0.0, 0, true, true, false,
                                                pseudos);

   if (nrestraints > 0) {
      restraints.analyze_for_bad_restraints(coot::restraints_container_t::BAD_RESTRAINT_ANALYSIS_INCLUDE_ANGLES);
      coot::geometry_distortion_info_container_t gdic = restraints.geometric_distortions();
      double d = gdic.distortion();
      const double k = 350;  // ratio derived from mon-lib analysis
      double ceu = d/(k);
      std::cout << "INFO:: " << pdb_file_name << " " << chain_id << " " << res_no << " " << monomer_type
                << " ligand distortion score:   " << ceu << " ceus" << std::endl;
   }

   if (include_environment_contacts_flag) {
      std::vector<mmdb::Residue *> neighbours = coot::residues_near_residue(residue_p, mol, 5);
      coot::atom_overlaps_container_t ao(residue_p, neighbours, mol, &geom);
      ao.make_overlaps();
      float os = ao.score();
      std::cout << "INFO:: " << pdb_file_name << " " << chain_id << " " << res_no << " " << monomer_type
                << " ligand atom overlap score: " << os << " A^3 / 1000 atoms" << std::endl;
   }

}

#include "coot-utils/atom-selection-container.hh"

void
read_data_validate_ligand(const input_data_t &input_data) {

   if (coot::file_exists(input_data.pdb_file_name)) {

      // ReadCoorFile doesn't add atom indexes UDData. So I get warning messages in
      // create_mmdbmanager_from_residue_vector(). So let's use get_atom_selection
      //mmdb::Manager *mol = new mmdb::Manager;
      // mol->ReadCoorFile(input_data.pdb_file_name.c_str());
      bool allow_dups = true;
      bool verbose = false;
      bool convert_flag = false;
      std::string pdb_file_name = input_data.pdb_file_name;
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, false, allow_dups, verbose);
      mmdb::Manager *mol = asc.mol;

      try {
         int res_no = coot::util::string_to_int(input_data.res_no_ligand_str);
         coot::residue_spec_t spec(input_data.chain_id, res_no, "");
         mmdb::Residue *residue_p = coot::util::get_residue(spec, mol);
         if (residue_p) {
            if (coot::file_exists(input_data.dictionary)) {

               coot::protein_geometry geom;
               geom.set_verbose(false);
               if (input_data.use_pocket)
                  geom.init_standard();
               int read_number = 44;
               std::vector<std::string> nsr = coot::util::non_standard_residue_types_in_molecule(mol);
               for (const auto &res_type : nsr)
                  geom.try_dynamic_add(res_type, read_number++);

               int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
               geom.init_refmac_mon_lib(input_data.dictionary, read_number, imol_enc);

               bool include_environment_contacts_flag = false;
               if (input_data.use_pocket)
                  include_environment_contacts_flag = true;
               validate_ligand(pdb_file_name, residue_p, mol, include_environment_contacts_flag, geom);

            } else {
                  std::cout << "Failed to find the restraints file " << input_data.dictionary << std::endl;
            }
         } else {
            std::cout << "residue not found " << input_data.pdb_file_name << " " << input_data.chain_id
                      << " " << input_data.res_no_ligand_str << std::endl;
         }
      }
      catch (const std::runtime_error &e) {
         std::cout << "ERROR::" << e.what() << std::endl;
      }
   }
}

int main(int argc, char **argv) {

   int status = 0;

   input_data_t input_data = get_input_details(argc, argv);

   if (input_data.is_good) {
         read_data_validate_ligand(input_data);
   } else {
      std::cout << "Incomprehensible input " << std::endl;
      show_usage();
   }

   return status;

}
