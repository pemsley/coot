/* ideal/extra-restraints.cc
 * 
 * Copyright 2010  by The University of Oxford
 * Copyright 2013, 2015 by Medical Research Council
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */


#include <stdexcept>
#include <algorithm>
#include <fstream>
#include <optional>

#include "simple-restraint.hh"
#include "utils/logging.hh"
extern logging logger;


void
coot::extra_restraints_t::add_start_pos_restraint(const atom_spec_t &atom_1_in, double esd_in) {
   bool already_exists = false;
   for (unsigned int i=0; i<start_pos_restraints.size(); i++) {
      if (start_pos_restraints[i].atom_1 == atom_1_in) {
         start_pos_restraints[i].esd = esd_in;
         already_exists = true;
         break;
      }
   }
   if (! already_exists) {
      extra_start_pos_restraint_t e(atom_1_in, esd_in);
      start_pos_restraints.push_back(e);
   }
}

int
coot::extra_restraints_t::read_refmac_extra_restraints(const std::string &file_name) {

   class torsion_data_t {
   public:
      atom_spec_t spec_1;
      atom_spec_t spec_2;
      atom_spec_t spec_3;
      atom_spec_t spec_4;
      float value;
      float sigma;
      torsion_data_t() : value(0), sigma(0) {}
      torsion_data_t(const atom_spec_t &a1, const atom_spec_t &a2, const atom_spec_t &a3, const atom_spec_t &a4) :
         spec_1(a1), spec_2(a2), spec_3(a3), spec_4(a4), value(0), sigma(0) {}
   };

   // parse_with_new_style_bonds_parser(words);
   auto get_bond_data = [] (const std::vector<std::string> &words_in) -> std::optional<extra_bond_restraint_t> {

      extra_bond_restraint_t br;
      unsigned int n_words = words_in.size();
      std::vector<std::string> words = words_in;
      for (std::string &word : words)
         word = util::upcase(word);
      unsigned int i_word = 0;

      // ---------------- atom 1 ---------------------

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "EXTE") {
               i_word++;
            }  else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "DIST") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "FIRS") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "CHAI") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         br.atom_1.chain_id = words[i_word];
         i_word++;
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "RESI") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         br.atom_1.res_no = util::string_to_int(words[i_word]);
         i_word++;
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 2) {
            if (words[i_word].substr(0,3) == "INS") {
               i_word++;
               if (i_word < n_words) {
                  br.atom_1.ins_code = words[i_word];
                  i_word++;
               }
            }
         }
      }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "ATOM") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         br.atom_1.atom_name = words[i_word];
         i_word++;
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 2) {
            if (words[i_word].substr(0,3) == "INS") {
               i_word++;
               if (i_word < n_words) {
                  br.atom_1.ins_code = words[i_word];
                  i_word++;
               }
            }
         }
      } else { return {}; }

      // ---------------- atom 2 ---------------------

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "SECO") {
               i_word++;
            }  else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "CHAI") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         br.atom_2.chain_id = words[i_word];
         i_word++;
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "RESI") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         br.atom_2.res_no = util::string_to_int(words[i_word]);
         i_word++;
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 2) {
            if (words[i_word].substr(0,3) == "INS") {
               i_word++;
               if (i_word < n_words) {
                  br.atom_2.ins_code = words[i_word];
                  i_word++;
               }
            }
         }
      }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "ATOM") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         br.atom_2.atom_name = words[i_word];
         i_word++;
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 2) {
            if (words[i_word].substr(0,3) == "INS") {
               i_word++;
               if (i_word < n_words) {
                  br.atom_2.ins_code = words[i_word];
                  i_word++;
               }
            }
         }
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "VALU") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         br.bond_dist = util::string_to_float(words[i_word]);
         i_word++;
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "SIGM") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         br.esd = util::string_to_float(words[i_word]);
         i_word++;
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "TYPE") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         br.type = util::string_to_int(words[i_word]);
         i_word++;
      } else { return {}; }

      return br;
   };

   // TIL to specify the return type of a lambda function
   auto get_torsion_data = [] (const std::vector<std::string> &words_in) -> std::optional<torsion_data_t> {

      torsion_data_t td;
      unsigned int n_words = words_in.size();
      std::vector<std::string> words = words_in;
      for (std::string &word : words)
         word = util::upcase(word);
      unsigned int i_word = 0;

      // ---------------- atom 1 ---------------------

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "EXTE") {
               i_word++;
            }  else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "TORS") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "FIRS") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "CHAI") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         td.spec_1.chain_id = words[i_word];
         i_word++;
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "RESI") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         td.spec_1.res_no = util::string_to_int(words[i_word]);
         i_word++;
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 2) {
            if (words[i_word].substr(0,3) == "INS") {
               i_word++;
               if (i_word < n_words) {
                  td.spec_1.ins_code = words[i_word];
                  i_word++;
               }
            }
         }
      }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "ATOM") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         td.spec_1.atom_name = words[i_word];
         i_word++;
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 2) {
            if (words[i_word].substr(0,3) == "INS") {
               i_word++;
               if (i_word < n_words) {
                  td.spec_1.ins_code = words[i_word];
                  i_word++;
               }
            }
         }
      }


      // ---------------- atom 2 ---------------------

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "NEXT") {
               i_word++;
            }  else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "CHAI") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         td.spec_2.chain_id = words[i_word];
         i_word++;
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "RESI") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         td.spec_2.res_no = util::string_to_int(words[i_word]);
         i_word++;
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 2) {
            if (words[i_word].substr(0,3) == "INS") {
               i_word++;
               if (i_word < n_words) {
                  td.spec_2.ins_code = words[i_word];
                  i_word++;
               }
            }
         }
      }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "ATOM") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         td.spec_2.atom_name = words[i_word];
         i_word++;
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 2) {
            if (words[i_word].substr(0,3) == "INS") {
               i_word++;
               if (i_word < n_words) {
                  td.spec_2.ins_code = words[i_word];
                  i_word++;
               }
            }
         }
      }

      // ---------------- atom 3 ---------------------

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "NEXT") {
               i_word++;
            }  else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "CHAI") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         td.spec_3.chain_id = words[i_word];
         i_word++;
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "RESI") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         td.spec_3.res_no = util::string_to_int(words[i_word]);
         i_word++;
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 2) {
            if (words[i_word].substr(0,3) == "INS") {
               i_word++;
               if (i_word < n_words) {
                  td.spec_3.ins_code = words[i_word];
                  i_word++;
               }
            }
         }
      }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "ATOM") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         td.spec_3.atom_name = words[i_word];
         i_word++;
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 2) {
            if (words[i_word].substr(0,3) == "INS") {
               i_word++;
               if (i_word < n_words) {
                  td.spec_3.ins_code = words[i_word];
                  i_word++;
               }
            }
         }
      }

      // ---------------- atom 4 ---------------------


      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "NEXT") {
               i_word++;
            }  else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "CHAI") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         td.spec_4.chain_id = words[i_word];
         i_word++;
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "RESI") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         td.spec_4.res_no = util::string_to_int(words[i_word]);
         i_word++;
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 2) {
            if (words[i_word].substr(0,3) == "INS") {
               i_word++;
               if (i_word < n_words) {
                  td.spec_4.ins_code = words[i_word];
                  i_word++;
               }
            }
         }
      }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "ATOM") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         td.spec_4.atom_name = words[i_word];
         i_word++;
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 2) {
            if (words[i_word].substr(0,3) == "INS") {
               i_word++;
               if (i_word < n_words) {
                  td.spec_4.ins_code = words[i_word];
                  i_word++;
               }
            }
         }
      }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "VALU") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         td.value = util::string_to_float(words[i_word]);
         i_word++;
      } else { return {}; }

      if (i_word < n_words) {
         if (words[i_word].length() > 3) {
            if (words[i_word].substr(0,4) == "SIGM") {
               i_word++;
            } else { return {}; }
         } else { return {}; }
      } else { return {}; }

      if (i_word < n_words) {
         td.sigma = util::string_to_float(words[i_word]);
         return td;
      } else { return {}; }

      return {};
   };

      // std::vector<extra_bond_restraint_t> bond_restraints;
      // std::vector<extra_angle_restraint_t> angle_restraints;
      // std::vector<extra_torsion_restraint_t> torsion_restraints;
      // std::vector<extra_geman_mcclure_restraint_t> geman_mcclure_restraints;
      // std::vector<extra_start_pos_restraint_t> start_pos_restraints;
      // std::vector<extra_target_position_restraint_t> target_position_restraints;
      // std::vector<parallel_planes_t> parallel_plane_restraints;

   std::map<std::string, int> restraints_counts_map_pre;
   restraints_counts_map_pre["bond"]    =    bond_restraints.size();
   restraints_counts_map_pre["angle"]   =   angle_restraints.size();
   restraints_counts_map_pre["torsion"] = torsion_restraints.size();
   restraints_counts_map_pre["geman_mcclure"]   =   geman_mcclure_restraints.size();
   restraints_counts_map_pre["start_pos"]       =       start_pos_restraints.size();
   restraints_counts_map_pre["target_position"] = target_position_restraints.size();
   restraints_counts_map_pre["parallel_plane"]  =  parallel_plane_restraints.size();

   if (file_exists(file_name)) {
      std::string line;
      std::vector<std::string> lines;
      std::ifstream f(file_name.c_str());
      while (std::getline(f, line)) {
         lines.push_back(line);
      }

      for (unsigned int i=0; i<lines.size(); i++) {
         const std::string line = lines[i];
         std::vector<std::string> words = util::split_string_no_blanks(lines[i], " ");
         if (matches_bond_template_p(words)) {
            try {
               std::string chain_id_1 = words[4];
               std::string ins_code_1 = words[8];
               std::string atm_name_1 = words[10];
               int res_no_1 = util::string_to_int(words[6]);

               std::string chain_id_2 = words[13];
               std::string ins_code_2 = words[17];
               std::string atm_name_2 = words[19];
               int res_no_2 = util::string_to_int(words[15]);

               if (ins_code_1 == ".") ins_code_1 = "";
               if (ins_code_2 == ".") ins_code_2 = "";

               double d = util::string_to_float(words[21]);
               double e = util::string_to_float(words[23]);

               std::string n1 = atom_id_mmdb_expand(atm_name_1);
               std::string n2 = atom_id_mmdb_expand(atm_name_2);

               atom_spec_t spec_1(chain_id_1, res_no_1, ins_code_1, n1, "");
               atom_spec_t spec_2(chain_id_2, res_no_2, ins_code_2, n2, "");
               extra_bond_restraint_t br(spec_1, spec_2, d, e);
               bond_restraints.push_back(br);
            }
            catch (const std::runtime_error &rte) {
               // silently ignore
               std::cout << "WARNING:: rte on : " << lines[i] << std::endl;
            }

         } else {

            // parse_with_new_style_bonds_parser(words);
            std::optional<extra_bond_restraint_t> new_style_bonds_parser = get_bond_data(words);

            if (new_style_bonds_parser) {

               bond_restraints.push_back(new_style_bonds_parser.value());

            } else {

               if (matches_angle_template_p(words)) {

                  try {
                     std::string chain_id_1 = words[4];
                     std::string ins_code_1 = words[8];
                     std::string atm_name_1 = words[10];
                     int res_no_1 = util::string_to_int(words[6]);

                     std::string chain_id_2 = words[13];
                     std::string ins_code_2 = words[17];
                     std::string atm_name_2 = words[19];
                     int res_no_2 = util::string_to_int(words[15]);

                     std::string chain_id_3 = words[22];
                     std::string ins_code_3 = words[26];
                     std::string atm_name_3 = words[28];
                     int res_no_3 = util::string_to_int(words[24]);

                     if (ins_code_1 == ".") ins_code_1 = "";
                     if (ins_code_2 == ".") ins_code_2 = "";
                     if (ins_code_3 == ".") ins_code_3 = "";

                     double d = util::string_to_float(words[30]);
                     double e = util::string_to_float(words[32]);

                     std::string n1 = atom_id_mmdb_expand(atm_name_1);
                     std::string n2 = atom_id_mmdb_expand(atm_name_2);
                     std::string n3 = atom_id_mmdb_expand(atm_name_3);

                     atom_spec_t spec_1(chain_id_1, res_no_1, ins_code_1, n1, "");
                     atom_spec_t spec_2(chain_id_2, res_no_2, ins_code_2, n2, "");
                     atom_spec_t spec_3(chain_id_3, res_no_3, ins_code_3, n3, "");
                     extra_angle_restraint_t ar(spec_1, spec_2, spec_3, d, e);
                     angle_restraints.push_back(ar);

                  }
                  catch (const std::runtime_error &rte) {
                     std::cout << "WARNING:: rte on : " << lines[i] << std::endl;
                  }

               } else {

                  if (matches_torsion_template_p(words)) {

                     try {
                        std::string chain_id_1 = words[4];
                        std::string ins_code_1 = words[8];
                        std::string atm_name_1 = words[10];
                        int res_no_1 = util::string_to_int(words[6]);

                        std::string chain_id_2 = words[13];
                        std::string ins_code_2 = words[17];
                        std::string atm_name_2 = words[19];
                        int res_no_2 = util::string_to_int(words[15]);

                        std::string chain_id_3 = words[22];
                        std::string ins_code_3 = words[26];
                        std::string atm_name_3 = words[28];
                        int res_no_3 = util::string_to_int(words[24]);

                        std::string chain_id_4 = words[31];
                        std::string ins_code_4 = words[35];
                        std::string atm_name_4 = words[37];
                        int res_no_4 = util::string_to_int(words[33]);

                        if (ins_code_1 == ".") ins_code_1 = "";
                        if (ins_code_2 == ".") ins_code_2 = "";
                        if (ins_code_3 == ".") ins_code_3 = "";
                        if (ins_code_4 == ".") ins_code_4 = "";

                        double d = util::string_to_float(words[39]);
                        double e = util::string_to_float(words[41]);

                        std::string n1 = atom_id_mmdb_expand(atm_name_1);
                        std::string n2 = atom_id_mmdb_expand(atm_name_2);
                        std::string n3 = atom_id_mmdb_expand(atm_name_3);
                        std::string n4 = atom_id_mmdb_expand(atm_name_4);

                        atom_spec_t spec_1(chain_id_1, res_no_1, ins_code_1, n1, "");
                        atom_spec_t spec_2(chain_id_2, res_no_2, ins_code_2, n2, "");
                        atom_spec_t spec_3(chain_id_3, res_no_3, ins_code_3, n3, "");
                        atom_spec_t spec_4(chain_id_4, res_no_4, ins_code_4, n4, "");
                        extra_torsion_restraint_t tr(spec_1, spec_2, spec_3, spec_4, d, e, 1);
                        torsion_restraints.push_back(tr);

                     }
                     catch (const std::runtime_error &rte) {
                        std::cout << "WARNING:: rte on : " << lines[i] << std::endl;
                     }

                  } else {

                     std::optional<torsion_data_t> torsion_data = get_torsion_data(words);

                     if (torsion_data) {
                        extra_torsion_restraint_t tr(torsion_data.value().spec_1, torsion_data.value().spec_2,
                                                     torsion_data.value().spec_3, torsion_data.value().spec_4,
                                                     torsion_data.value().value, torsion_data.value().sigma, 1);
                        if (false) // debugging
                           std::cout << "read torsion "
                                     << torsion_data.value().spec_1 << " " << torsion_data.value().spec_2 << " "
                                     << torsion_data.value().spec_3 << " " << torsion_data.value().spec_4 << " "
                                     << torsion_data.value().value  << " " << torsion_data.value().sigma
                                     << std::endl;

                        torsion_restraints.push_back(tr);

                     } else {

                        parallel_planes_t ppr(lines[i]); // try to parse the line and make a restraint
                        if (ppr.matches) {
                           // add parallel plane (aka "stacking") restraint
                           parallel_plane_restraints.push_back(ppr);
                        } else {
                           logger.log(log_t::INFO, logging::ltw("Failed to match this restraint to templates: "), logging::ltw(line));
                           // std::cout << "INFO:: Failed to match this restraint to templates:\n   "
                           //           << line << std::endl;
                        }
                     }
                  }
               }
            }
         }
      }
   }
   std::map<std::string, int> restraints_counts_map_post;
   restraints_counts_map_post["bond"]    =    bond_restraints.size();
   restraints_counts_map_post["angle"]   =   angle_restraints.size();
   restraints_counts_map_post["torsion"] = torsion_restraints.size();
   restraints_counts_map_post["geman_mcclure"]   =   geman_mcclure_restraints.size();
   restraints_counts_map_post["start_pos"]       =       start_pos_restraints.size();
   restraints_counts_map_post["target_position"] = target_position_restraints.size();
   restraints_counts_map_post["parallel_plane"]  =  parallel_plane_restraints.size();

   int n_delta = 0;
   for (auto t : std::vector{"bond", "angle", "torsion", "geman_mcclure", "start_pos", "target_position", "parallel_plane"}) {
      n_delta += restraints_counts_map_post[t] - restraints_counts_map_pre[t];
   }
   return n_delta;
}

void
coot::extra_restraints_t::delete_restraints_for_residue(const residue_spec_t &rs) {

   bond_restraints.erase(std::remove(bond_restraints.begin(), bond_restraints.end(), rs), bond_restraints.end());
}

bool
coot::extra_restraints_t::matches_torsion_template_p(const std::vector<std::string> &words) const {

   bool status = false;

   if (words.size() == 42 || words.size() == 44) { // allow type field at the end
      std::vector<std::string> u(words.size());
      for (unsigned int i=0; i<words.size(); i++)
         u[i] = coot::util::upcase(words[i]);
      if (u[0].length() > 3) {
         if (u[0].substr(0,4) == "EXTE") {
            if (u[1].length() > 3) {
               if (u[1].substr(0,4) == "TORS") {
                  if (u[2].length() > 3) {
                     if (u[2].substr(0,4) == "FIRS") {
                        if (u[3].length() > 3) {
                           if (u[3].substr(0,4) == "CHAI") {
                              if (u[5].length() > 3) {
                                 if (u[5].substr(0,4) == "RESI") {
                                    if (u[7].length() > 2) {
                                       if (u[7].substr(0,3) == "INS") {
                                          if (u[9].length() > 3) {
                                             if (u[9].substr(0,4) == "ATOM") {
                                                if (u[11].length() > 3) {
                                                   if (u[11].substr(0,4) == "NEXT") {
                                                      if (u[12].length() > 3) {
                                                         if (u[12].substr(0,4) == "CHAI") {
                                                            if (u[14].length() > 3) {
                                                               if (u[14].substr(0,4) == "RESI") {
                                                                  if (u[16].length() > 2) {
                                                                     if (u[16].substr(0,3) == "INS") {
                                                                        if (u[18].length() > 3) {
                                                                           if (u[18].substr(0,4) == "ATOM") {
                                                                              if (u[20].length() > 3) {
                                                                                 if (u[20].substr(0,4) == "NEXT") {
                                                                                    if (u[21].length() > 3) {
                                                                                       if (u[21].substr(0,4) == "CHAI") {
                                                                                          if (u[23].length() > 3) {
                                                                                             if (u[23].substr(0,4) == "RESI") {
                                                                                                if (u[25].length() > 2) {
                                                                                                   if (u[25].substr(0,3) == "INS") {
                                                                                                      if (u[27].length() > 3) {
                                                                                                         if (u[27].substr(0,4) == "ATOM") {
                                                                                                            if (u[29].length() > 3) {
                                                                                                               if (u[29].substr(0,4) == "NEXT") {
                                                                                                                  if (u[30].length() > 3) {
                                                                                                                     if (u[30].substr(0,4) == "CHAI") {
                                                                                                                        if (u[32].length() > 3) {
                                                                                                                           if (u[32].substr(0,4) == "RESI") {
                                                                                                                              if (u[34].length() > 2) {
                                                                                                                                 if (u[34].substr(0,3) == "INS") {
                                                                                                                                    if (u[36].length() > 3) {
                                                                                                                                       if (u[36].substr(0,4) == "ATOM") {
                                                                                                                                          if (u[38].length() > 3) {
                                                                                                                                             if (u[38].substr(0,4) == "VALU") {
                                                                                                                                                if (u[40].length() > 3) {
                                                                                                                                                   if (u[40].substr(0,4) == "SIGM") {
                                                                                                                                                      status = true;
                                                                                                                                                   }
                                                                                                                                                }
                                                                                                                                             }
                                                                                                                                          }
                                                                                                                                       }
                                                                                                                                    }
                                                                                                                                 }
                                                                                                                              }
                                                                                                                           }
                                                                                                                        }
                                                                                                                     }
                                                                                                                  }
                                                                                                               }
                                                                                                            }
                                                                                                         }
                                                                                                      }
                                                                                                   }
                                                                                                }
                                                                                             }
                                                                                          }
                                                                                       }
                                                                                    }
                                                                                 }
                                                                              }
                                                                           }
                                                                        }
                                                                     }
                                                                  }
                                                               }
                                                            }
                                                         }
                                                      }
                                                   }
                                                }
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

   if (false) { // debugging
      for (unsigned int i=0; i<words.size(); i++) {
         std::cout << " " << words[i] <<  " ";
      }
      std::cout << "status: " << status << std::endl;
   }
   return status;
}

bool
coot::extra_restraints_t::matches_angle_template_p(const std::vector<std::string> &words) const {

   bool status = false;

   if (words.size() == 33 || words.size() == 35) { // allow type field at the end
      std::vector<std::string> u(words.size());
      for (unsigned int i=0; i<words.size(); i++)
         u[i] = coot::util::upcase(words[i]);
      if (u[0].length() > 3) {
         if (u[0].substr(0,4) == "EXTE") {
            if (u[1].length() > 3) {
               if (u[1].substr(0,4) == "ANGL") {
                  if (u[2].length() > 3) {
                     if (u[2].substr(0,4) == "FIRS") {
                        if (u[3].length() > 3) {
                           if (u[3].substr(0,4) == "CHAI") {
                              if (u[5].length() > 3) {
                                 if (u[5].substr(0,4) == "RESI") {
                                    if (u[7].length() > 2) {
                                       if (u[7].substr(0,3) == "INS") {
                                          if (u[9].length() > 3) {
                                             if (u[9].substr(0,4) == "ATOM") {
                                                if (u[11].length() > 3) {
                                                   if (u[11].substr(0,4) == "NEXT") {
                                                      if (u[12].length() > 3) {
                                                         if (u[12].substr(0,4) == "CHAI") {
                                                            if (u[14].length() > 3) {
                                                               if (u[14].substr(0,4) == "RESI") {
                                                                  if (u[16].length() > 2) {
                                                                     if (u[16].substr(0,3) == "INS") {
                                                                        if (u[18].length() > 3) {
                                                                           if (u[18].substr(0,4) == "ATOM") {
                                                                              if (u[20].length() > 3) {
                                                                                 if (u[20].substr(0,4) == "NEXT") {
                                                                                    if (u[21].length() > 3) {
                                                                                       if (u[21].substr(0,4) == "CHAI") {
                                                                                          if (u[23].length() > 3) {
                                                                                             if (u[23].substr(0,4) == "RESI") {
                                                                                                if (u[25].length() > 2) {
                                                                                                   if (u[25].substr(0,3) == "INS") {
                                                                                                      if (u[27].length() > 3) {
                                                                                                         if (u[27].substr(0,4) == "ATOM") {
                                                                                                            if (u[29].length() > 3) {
                                                                                                               if (u[29].substr(0,4) == "VALU") {
                                                                                                                  if (u[31].length() > 3) {
                                                                                                                     if (u[31].substr(0,4) == "SIGM") {
                                                                                                                        status = true;
                                                                                                                     }
                                                                                                                  }
                                                                                                               }
                                                                                                            }
                                                                                                         }
                                                                                                      }
                                                                                                   }
                                                                                                }
                                                                                             }
                                                                                          }
                                                                                       }
                                                                                    }
                                                                                 }
                                                                              }
                                                                           }
                                                                        }
                                                                     }
                                                                  }
                                                               }
                                                            }
                                                         }
                                                      }
                                                   }
                                                }
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

   if (false) { // debugging
      for (unsigned int i=0; i<words.size(); i++) {
         std::cout << " " << words[i] <<  " ";
      }
      std::cout << "status: " << status << std::endl;
   }
   return status;
}


bool
coot::extra_restraints_t::matches_bond_template_p(const std::vector<std::string> &words) const {

   bool status = false;
   if (words.size() >= 24 || words.size() == 26) { // allow type field at the end
      std::vector<std::string> u(words.size());
      for (unsigned int i=0; i<words.size(); i++)
         u[i] = coot::util::upcase(words[i]);
      if (u[0].length() > 3) {
         if (u[0].substr(0,4) == "EXTE") {
            if (u[1].length() > 3) {
               if (u[1].substr(0,4) == "DIST") {
                  if (u[2].length() > 3) {
                     if (u[2].substr(0,4) == "FIRS") {
                        if (u[3].length() > 3) {
                           if (u[3].substr(0,4) == "CHAI") {
                              if (u[5].length() > 3) {
                                 if (u[5].substr(0,4) == "RESI") {
                                    if (u[7].length() > 2) {
                                       if (u[7].substr(0,3) == "INS") {
                                          if (u[9].length() > 3) {
                                             if (u[9].substr(0,4) == "ATOM") {
                                                if (u[11].length() > 3) {
                                                   if (u[11].substr(0,4) == "SECO") {
                                                      if (u[12].length() > 3) {
                                                         if (u[12].substr(0,4) == "CHAI") {
                                                            if (u[14].length() > 3) {
                                                               if (u[14].substr(0,4) == "RESI") {
                                                                  if (u[16].length() > 2) {
                                                                     if (u[16].substr(0,3) == "INS") {
                                                                        if (u[18].length() > 3) {
                                                                           if (u[18].substr(0,4) == "ATOM") {
                                                                              if (u[20].length() > 3) {
                                                                                 if (u[20].substr(0,4) == "VALU") {
                                                                                    if (u[22].length() > 3) {
                                                                                       if (u[22].substr(0,4) == "SIGM") {
                                                                                          status = true;
                                                                                       }
                                                                                    }
                                                                                 }
                                                                              }
                                                                           }
                                                                        }
                                                                     }
                                                                  }
                                                               }
                                                            }
                                                         }
                                                      }
                                                   }
                                                }
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   } else {
      if (false)
         std::cout << "WARNING:: matches_bond_template_p() test: not 24 words" << std::endl;
   }
   return status;
}



void
coot::restraints_container_t::add_extra_restraints(int imol,
                                                   const std::string &description,
                                                   const extra_restraints_t &extra_restraints,
                                                   const protein_geometry &geom) {

   // std::cout << "HHHHHHHHHere in add_extra_restraints() - from " << description << std::endl;

   if (false) {
      logger.log(log_t::INFO, logging::ltw("extra restraints origin: "), logging::ltw(description));
      // std::cout << "INFO:: extra restraints origin: " << description << std::endl;
      std::cout << "--------------------- in add_extra_restraints() we have "
                << restraints_vec.size() << " standard restraints "
                << std::endl;
      std::cout << "--------------------- in add_extra_restraints() adding "
                << extra_restraints.bond_restraints.size() << " extra bond restraints "
                << std::endl;
      std::cout << "--------------------- in add_extra_restraints() adding "
                << extra_restraints.geman_mcclure_restraints.size() << " extra GM restraints "
                << std::endl;
      std::cout << "--------------------- in add_extra_restraints() adding "
                << extra_restraints.angle_restraints.size() << " extra angle restraints "
                << std::endl;
      std::cout << "--------------------- in add_extra_restraints() adding "
                << extra_restraints.torsion_restraints.size() << " extra torsion restraints "
                << std::endl;
      std::cout << "--------------------- in add_extra_restraints() parallel-plane adding "
                << extra_restraints.parallel_plane_restraints.size() << " pp restraints "
                << std::endl;
      std::cout << "--------------------- in add_extra_restraints() target-position adding "
                << extra_restraints.target_position_restraints.size() << " position restraints "
                << std::endl;
   }

   add_extra_bond_restraints(extra_restraints);
   add_extra_geman_mcclure_restraints(extra_restraints);
   add_extra_angle_restraints(extra_restraints);
   add_extra_torsion_restraints(extra_restraints);
   add_extra_start_pos_restraints(extra_restraints);
   add_extra_target_position_restraints(extra_restraints);
   add_extra_parallel_plane_restraints(imol, extra_restraints, geom);
   make_restraint_types_index_limits();

   post_add_new_restraints();
}

void
coot::restraints_container_t::fill_old_to_new_index_vector() {

   bool debug = false;

   int hnd = mol->GetUDDHandle(mmdb::UDR_ATOM, "old atom index");
   old_atom_index_to_new_atom_index.resize(ATOM_INDEX_MAX, -1);
   for (int i=0; i<n_atoms; i++) {
      int idx;
      int ierr = atom[i]->GetUDData(hnd, idx);
      if (ierr == mmdb::UDDATA_Ok) {
         if ((idx >= 0) && (idx < ATOM_INDEX_MAX))
            if (debug)
               std::cout << "debug:: fill_old_to_new_index_vector(): converting index " << idx << " " << i << std::endl;
         old_atom_index_to_new_atom_index[idx] = i;
      } else {
         if (debug)
            std::cout << "debug:: fill_old_to_new_index_vector(): mmdb error atom " << i << std::endl;
      }
   }
}

void
coot::restraints_container_t::add_extra_target_position_restraints(const extra_restraints_t &extra_restraints) {

   // it doesn't make sense to add a target position restraint for a fixed atom

   for (unsigned int i=0; i<extra_restraints.target_position_restraints.size(); i++) {
      const extra_restraints_t::extra_target_position_restraint_t &pr =
         extra_restraints.target_position_restraints[i];

      mmdb::Residue *residue_p = NULL;
      mmdb::Atom *at = 0;
      bool fixed = false;
      if (from_residue_vector) {
         residue_spec_t res_spec(pr.atom_spec);
         for (unsigned int ir=0; ir<residues_vec.size(); ir++) {
            if (residue_spec_t(residues_vec[ir].second) == res_spec) {
               residue_p = residues_vec[ir].second;
               fixed = residues_vec[ir].first;
               break;
            }
         }
      } else {
         // Fill me
      }

      if (! fixed) {

         if (residue_p) {
            // set "at" from atoms in residue_p
            mmdb::Atom **residue_atoms = 0;
            int n_residue_atoms;
            residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
            for (int iat=0; iat<n_residue_atoms; iat++) {
               std::string atom_name(residue_atoms[iat]->name);
               if (atom_name == extra_restraints.target_position_restraints[i].atom_spec.atom_name) {
                  std::string alt_loc(residue_atoms[iat]->altLoc);
                  if (alt_loc == extra_restraints.target_position_restraints[i].atom_spec.alt_conf) {
                     at = residue_atoms[iat];
                     break;
                  }
               }
            }
         }

         if (at) {
            int atom_index = -1;
            at->GetUDData(udd_atom_index_handle, atom_index);
            // std::cout << "add target position for atom index " << atom_index << std::endl;
            add_user_defined_target_position_restraint(TARGET_POS_RESTRAINT, atom_index, pr.atom_spec, pr.pos, pr.weight);
         }
      }
   }
}

// can I find the atoms using the atom indices from the original molecule?
bool
coot::restraints_container_t::try_add_using_old_atom_indices(const extra_restraints_t::extra_bond_restraint_t &ebr) {
   bool success = false;
   bool debug = false;

   int idx_1_old = ebr.atom_1.int_user_data; // where do these get set? (they seem to work though!)
   int idx_2_old = ebr.atom_2.int_user_data;

   if (debug)
      std::cout << "debug:: idx_1_old: " << idx_1_old << " idx_2_old: " << idx_2_old << "\n";

   if (idx_1_old >= 0 && idx_1_old < ATOM_INDEX_MAX) {
      if (idx_2_old >= 0 && idx_2_old < ATOM_INDEX_MAX) {
         int index_1 = old_atom_index_to_new_atom_index[idx_1_old];
         int index_2 = old_atom_index_to_new_atom_index[idx_2_old];
         if (debug)
            std::cout << "debug:: index_1: " << index_1 << " index_2: " << index_2 << "\n";
         if ((index_1 != -1) && (index_2 != -1)) {
            std::vector<bool> fixed_flags = make_fixed_flags(index_1, index_2);

            add(BOND_RESTRAINT, index_1, index_2, fixed_flags, ebr.bond_dist, ebr.esd, 1.2 /* dummy value */);

            if (debug)
               std::cout << "DEBUG:: extra bond restraint done fast\n";
            success = true;
         }
      }
   }
   return success;
}

bool
coot::restraints_container_t::try_add_using_old_atom_indices(const extra_restraints_t::extra_geman_mcclure_restraint_t &ebr) {

   bool success = false;
   bool debug = false;

   int idx_1_old = ebr.atom_1.int_user_data; // where do these get set? (they seem to work though!)
   int idx_2_old = ebr.atom_2.int_user_data;

   if (debug)
      std::cout << "debug:: idx_1_old: " << idx_1_old << " idx_2_old: " << idx_2_old << "\n";

   if (idx_1_old >= 0 && idx_1_old < ATOM_INDEX_MAX) {
      if (idx_2_old >= 0 && idx_2_old < ATOM_INDEX_MAX) {
         int index_1 = old_atom_index_to_new_atom_index[idx_1_old];
         int index_2 = old_atom_index_to_new_atom_index[idx_2_old];
         if (debug)
            std::cout << "debug:: index_1: " << index_1 << " index_2: " << index_2 << "\n";
         if ((index_1 != -1) && (index_2 != -1)) {
            std::vector<bool> fixed_flags = make_fixed_flags(index_1, index_2);
            add_geman_mcclure_distance(GEMAN_MCCLURE_DISTANCE_RESTRAINT, index_1, index_2, fixed_flags,
                                       ebr.bond_dist, ebr.esd);
            if (debug)
               std::cout << "DEBUG:: GM done fast\n";
            success = true;
         }
      }
   }

   if (false)
      std::cout << "try_add_using_old_atom_indices() added status for GM restraint " << success << std::endl;
   return success;

}



void
coot::restraints_container_t::add_extra_bond_restraints(const extra_restraints_t &extra_restraints) {

   int n_extra_bond_restraints = 0;
   // don't add the restraint if both the residues are fixed.
   //

   fill_old_to_new_index_vector();
   // pre-calculate the residue specs for speed
   std::vector<residue_spec_t> residues_vec_residue_specs(residues_vec.size());
   for (unsigned int ir=0; ir<residues_vec.size(); ir++)
      residues_vec_residue_specs[ir] = residue_spec_t(residues_vec[ir].second);

   for (unsigned int i=0; i<extra_restraints.bond_restraints.size(); i++) {

      mmdb::Residue *r_1 = NULL;
      mmdb::Residue *r_2 = NULL;
      mmdb::Atom *at_1 = 0;
      mmdb::Atom *at_2 = 0;
      bool fixed_1 = false;
      bool fixed_2 = false;

      const extra_restraints_t::extra_bond_restraint_t &ebr = extra_restraints.bond_restraints[i];

      bool done = try_add_using_old_atom_indices(ebr);
      if (done) continue;
      if (from_residue_vector) {
         residue_spec_t br_res_atom_1(ebr.atom_1);
         residue_spec_t br_res_atom_2(ebr.atom_2);
         for (unsigned int ir=0; ir<residues_vec.size(); ir++) {
            if (!r_1) {
               if (br_res_atom_1 == residues_vec_residue_specs[ir]) {
                  r_1 = residues_vec[ir].second;
                  fixed_1 = residues_vec[ir].first;
               }
            }
            if (! r_2) {
               if (br_res_atom_2 == residues_vec_residue_specs[ir]) {
                  r_2 = residues_vec[ir].second;
                  fixed_2 = residues_vec[ir].first;
               }
            }
            if (r_1 && r_2) break;
         }
      } else {

         // bleugh.
         int selHnd = mol->NewSelection();  // d
         mol->Select (selHnd, mmdb::STYPE_RESIDUE, 1,       // .. TYPE, iModel
                      chain_id_save.c_str(), // Chain(s)
                      istart_res, "*", // starting res
                      iend_res,   "*", // ending   res
                      "*",  // residue name
                      "*",  // Residue must contain this atom name?
                      "*",  // Residue must contain this Element?
                      "*",  // altLocs
                      mmdb::SKEY_NEW // selection key
                      );
         int nSelResidues_local = 0;
         mmdb::PPResidue SelResidue_local= 0;
         mol->GetSelIndex (selHnd, SelResidue_local, nSelResidues_local);
         for (int ir=0; ir<nSelResidues_local; ir++) {
            if (coot::residue_spec_t(ebr.atom_1) == coot::residue_spec_t(SelResidue_local[ir])) {
               r_1 = SelResidue_local[ir];
               fixed_1 = fixed_check(ir);
            }
            if (coot::residue_spec_t(ebr.atom_2) == coot::residue_spec_t(SelResidue_local[ir])) {
               r_2 = SelResidue_local[ir];
               fixed_1 = fixed_check(ir);
            }
         } 
         mol->DeleteSelection(selHnd);

      }
      
      if (r_1 && r_2) {
         if (! (fixed_1 && fixed_2)) {
            mmdb::PPAtom residue_atoms_1 = 0;
            mmdb::PPAtom residue_atoms_2 = 0;
            int n_residue_atoms_1;
            int n_residue_atoms_2;
            r_1->GetAtomTable(residue_atoms_1, n_residue_atoms_1);
            r_2->GetAtomTable(residue_atoms_2, n_residue_atoms_2);

            for (int iat=0; iat<n_residue_atoms_1; iat++) { 
               std::string atom_name_1(residue_atoms_1[iat]->name);
               if (atom_name_1 == ebr.atom_1.atom_name) {
                  std::string alt_loc_1(residue_atoms_1[iat]->altLoc);
                  if (alt_loc_1 == ebr.atom_1.alt_conf) {
                     at_1 = residue_atoms_1[iat];
                     break;
                  }
               }
            }
            for (int iat=0; iat<n_residue_atoms_2; iat++) { 
               std::string atom_name_2(residue_atoms_2[iat]->name);
               if (atom_name_2 == ebr.atom_2.atom_name) {
                  std::string alt_loc_2(residue_atoms_2[iat]->altLoc);
                  if (alt_loc_2 == ebr.atom_2.alt_conf) {
                     at_2 = residue_atoms_2[iat];
                     break;
                  }
               }
            }

            if (at_1 && at_2) {
               int index_1 = -1;
               int index_2 = -1;
               at_1->GetUDData(udd_atom_index_handle, index_1);
               at_2->GetUDData(udd_atom_index_handle, index_2);
               if ((index_1 != -1) && (index_2 != -1)) {
                  std::vector<bool> fixed_flags = make_fixed_flags(index_1, index_2);

                  add(BOND_RESTRAINT, index_1, index_2, fixed_flags, ebr.bond_dist, ebr.esd, 1.2 /* dummy value */);

                  n_extra_bond_restraints++;

                  //mark these atoms as bonded so that we don't add a non-bonded restraint between them.
                  // 20170423 - But we *do* want NBC between these atoms...
                  //          otherwise we get horrid crunching.
                  // 20191119 However, if there is a *bond* restraint (not a GM restraint)
                  // then we don't also want a NBC.

                  bonded_atom_indices[index_1].insert(index_2);
                  bonded_atom_indices[index_2].insert(index_1);
               }
            }
         }
      }
   }
   if (false)
      logger.log(log_t::INFO, logging::ltw("--------------------------  made "), logging::ltw(n_extra_bond_restraints), logging::ltw(" extra bond restraints"));
      // std::cout << "INFO:: --------------------------  made " << n_extra_bond_restraints
      //           << " extra bond restraints" << std::endl;
}

void
coot::restraints_container_t::add_extra_geman_mcclure_restraints(const extra_restraints_t &extra_restraints) {

   int n_extra_bond_restraints = 0;
   // don't add the restraint if both the residues are fixed.
   //

   fill_old_to_new_index_vector();
   // pre-calculate the residue specs for speed
   std::vector<residue_spec_t> residues_vec_residue_specs(residues_vec.size());
   for (unsigned int ir=0; ir<residues_vec.size(); ir++)
      residues_vec_residue_specs[ir] = residue_spec_t(residues_vec[ir].second);

   for (unsigned int i=0; i<extra_restraints.geman_mcclure_restraints.size(); i++) {
      mmdb::Residue *r_1 = NULL;
      mmdb::Residue *r_2 = NULL;
      bool fixed_1 = false;
      bool fixed_2 = false;

      const extra_restraints_t::extra_geman_mcclure_restraint_t &ebr = extra_restraints.geman_mcclure_restraints[i];

      // this is the fast/right way to do it if it works
      bool done = try_add_using_old_atom_indices(ebr);
      if (done) {
         n_extra_bond_restraints++;
         continue;
      }

      if (from_residue_vector) {
         residue_spec_t br_res_atom_1(extra_restraints.geman_mcclure_restraints[i].atom_1);
         residue_spec_t br_res_atom_2(extra_restraints.geman_mcclure_restraints[i].atom_2);
         for (unsigned int ir=0; ir<residues_vec.size(); ir++) {
            if (!r_1) {
               if (br_res_atom_1 == residues_vec_residue_specs[ir]) {
                  r_1 = residues_vec[ir].second;
                  fixed_1 = residues_vec[ir].first;
               }
            }
            if (! r_2) {
               if (br_res_atom_2 == residues_vec_residue_specs[ir]) {
                  r_2 = residues_vec[ir].second;
                  fixed_2 = residues_vec[ir].first;
               }
            }
            if (r_1 && r_2) break;
         }
      } else {

         // bleugh.
         int selHnd = mol->NewSelection();  // d
         mol->Select (selHnd, mmdb::STYPE_RESIDUE, 1,       // .. TYPE, iModel
                      chain_id_save.c_str(), // Chain(s)
                      istart_res, "*", // starting res
                      iend_res,   "*", // ending   res
                      "*",  // residue name
                      "*",  // Residue must contain this atom name?
                      "*",  // Residue must contain this Element?
                      "*",  // altLocs
                      mmdb::SKEY_NEW // selection key
                      );
         int nSelResidues_local = 0;
         mmdb::PPResidue SelResidue_local= 0;
         mol->GetSelIndex (selHnd, SelResidue_local, nSelResidues_local);
         for (int ir=0; ir<nSelResidues_local; ir++) {
            if (coot::residue_spec_t(ebr.atom_1) ==
                coot::residue_spec_t(SelResidue_local[ir])) {
               r_1 = SelResidue_local[ir];
               fixed_1 = fixed_check(ir);
            }
            if (coot::residue_spec_t(ebr.atom_2) ==
                coot::residue_spec_t(SelResidue_local[ir])) {
               r_2 = SelResidue_local[ir];
               fixed_1 = fixed_check(ir);
            }
         } 
         mol->DeleteSelection(selHnd);
      }

      if (r_1 && r_2) {
         if (! (fixed_1 && fixed_2)) {
            mmdb::Atom *at_1 = 0;
            mmdb::Atom *at_2 = 0;

            mmdb::PPAtom residue_atoms_1 = 0;
            mmdb::PPAtom residue_atoms_2 = 0;
            int n_residue_atoms_1;
            int n_residue_atoms_2;
            r_1->GetAtomTable(residue_atoms_1, n_residue_atoms_1);
            r_2->GetAtomTable(residue_atoms_2, n_residue_atoms_2);

            for (int iat=0; iat<n_residue_atoms_1; iat++) { 
               std::string atom_name_1(residue_atoms_1[iat]->name);
               if (atom_name_1 == ebr.atom_1.atom_name) {
                  std::string alt_loc_1(residue_atoms_1[iat]->altLoc);
                  if (alt_loc_1 == ebr.atom_1.alt_conf) {
                     at_1 = residue_atoms_1[iat];
                     break;
                  }
               }
            }
            for (int iat=0; iat<n_residue_atoms_2; iat++) { 
               std::string atom_name_2(residue_atoms_2[iat]->name);
               if (atom_name_2 == ebr.atom_2.atom_name) {
                  std::string alt_loc_2(residue_atoms_2[iat]->altLoc);
                  if (alt_loc_2 == ebr.atom_2.alt_conf) {
                     at_2 = residue_atoms_2[iat];
                     break;
                  }
               }
            }

            if (at_1 && at_2) {
               int index_1 = -1; 
               int index_2 = -1;
               at_1->GetUDData(udd_atom_index_handle, index_1);
               at_2->GetUDData(udd_atom_index_handle, index_2);
               if ((index_1 != -1) && (index_2 != -1)) { 
                  std::vector<bool> fixed_flags = make_fixed_flags(index_1, index_2);

                   add_geman_mcclure_distance(GEMAN_MCCLURE_DISTANCE_RESTRAINT, index_1, index_2, fixed_flags,
                                              ebr.bond_dist, ebr.esd);
                  n_extra_bond_restraints++;

               }
            } 
         } 
      } 
   }

   // hopfully they have been added using the fast method,
   if (false)
      logger.log(log_t::INFO, logging::ltw("--------------------------  add_extra_geman_mcclure_restraints() made "), logging::ltw(n_extra_bond_restraints), logging::ltw(" extra GM restraints"));
      // std::cout << "INFO:: --------------------------  add_extra_geman_mcclure_restraints() made "
      //           << n_extra_bond_restraints << " extra GM restraints" << std::endl;

}


void
coot::restraints_container_t::add_extra_torsion_restraints(const extra_restraints_t &extra_restraints) {

   // std::cout << "in add_extra_torsion_restraints() << extra_restraints.torsion_restraints "
   // << extra_restraints.torsion_restraints.size() << std::endl;

   for (unsigned int i=0; i<extra_restraints.torsion_restraints.size(); i++) {

      if (false) {
         std::cout << "in add_extra_torsion_restraints() loop " << i << std::endl;
         std::cout << " "
                   << extra_restraints.torsion_restraints[i].atom_1 << "\n "
                   << extra_restraints.torsion_restraints[i].atom_2 << "\n "
                   << extra_restraints.torsion_restraints[i].atom_3 << "\n "
                   << extra_restraints.torsion_restraints[i].atom_4 << std::endl;
      }

      mmdb::Residue *r_1 = NULL;
      mmdb::Residue *r_2 = NULL;
      mmdb::Residue *r_3 = NULL;
      mmdb::Residue *r_4 = NULL;
      mmdb::Atom *at_1 = 0;
      mmdb::Atom *at_2 = 0;
      mmdb::Atom *at_3 = 0;
      mmdb::Atom *at_4 = 0;
      bool fixed_1 = 0;
      bool fixed_2 = 0;
      bool fixed_3 = 0;
      bool fixed_4 = 0;
      if (from_residue_vector) {
         for (unsigned int ir=0; ir<residues_vec.size(); ir++) {
            if (coot::residue_spec_t(extra_restraints.torsion_restraints[i].atom_1) ==
                coot::residue_spec_t(residues_vec[ir].second)) {
               r_1 = residues_vec[ir].second;
               fixed_1 = residues_vec[ir].first;
            }
            if (coot::residue_spec_t(extra_restraints.torsion_restraints[i].atom_2) ==
                coot::residue_spec_t(residues_vec[ir].second)) {
               r_2 = residues_vec[ir].second;
               fixed_2 = residues_vec[ir].first;
            }
            if (coot::residue_spec_t(extra_restraints.torsion_restraints[i].atom_3) ==
                coot::residue_spec_t(residues_vec[ir].second)) {
               r_3 = residues_vec[ir].second;
               fixed_3 = residues_vec[ir].first;
            }
            if (coot::residue_spec_t(extra_restraints.torsion_restraints[i].atom_4) ==
                coot::residue_spec_t(residues_vec[ir].second)) {
               r_4 = residues_vec[ir].second;
               fixed_4 = residues_vec[ir].first;
            }
         }
      } else {

         // bleugh.
         int selHnd = mol->NewSelection(); // d
         mol->Select (selHnd, mmdb::STYPE_RESIDUE, 1,       // .. TYPE, iModel
                      chain_id_save.c_str(), // Chain(s)
                      istart_res, "*", // starting res
                      iend_res,   "*", // ending   res
                      "*",  // residue name
                      "*",  // Residue must contain this atom name?
                      "*",  // Residue must contain this Element?
                      "*",  // altLocs
                      mmdb::SKEY_NEW // selection key
                      );
         int nSelResidues_local = 0;
         mmdb::PPResidue SelResidue_local= 0;
         mol->GetSelIndex (selHnd, SelResidue_local, nSelResidues_local);
         for (int ir=0; ir<nSelResidues_local; ir++) {
            if (coot::residue_spec_t(extra_restraints.torsion_restraints[i].atom_1) ==
                coot::residue_spec_t(SelResidue_local[ir]))
               r_1 = SelResidue_local[ir];
            if (coot::residue_spec_t(extra_restraints.torsion_restraints[i].atom_2) ==
                coot::residue_spec_t(SelResidue_local[ir]))
               r_2 = SelResidue_local[ir];
            if (coot::residue_spec_t(extra_restraints.torsion_restraints[i].atom_3) ==
                coot::residue_spec_t(SelResidue_local[ir]))
               r_3 = SelResidue_local[ir];
            if (coot::residue_spec_t(extra_restraints.torsion_restraints[i].atom_4) ==
                coot::residue_spec_t(SelResidue_local[ir]))
               r_4 = SelResidue_local[ir];
         }
         mol->DeleteSelection(selHnd);
      }

      if (r_1 && r_2 && r_3 && r_4) {
         mmdb::PPAtom residue_atoms = 0;
         int n_residue_atoms;
         r_1->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            if (coot::atom_spec_t(residue_atoms[iat]) == extra_restraints.torsion_restraints[i].atom_1) {
               at_1 = residue_atoms[iat];
               break;
            } 
         }
         residue_atoms = 0; // just to be safe
         r_2->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            if (coot::atom_spec_t(residue_atoms[iat]) == extra_restraints.torsion_restraints[i].atom_2) {
               at_2 = residue_atoms[iat];
               break;
            } 
         }
         residue_atoms = 0; // just to be safe
         r_3->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            if (coot::atom_spec_t(residue_atoms[iat]) == extra_restraints.torsion_restraints[i].atom_3) {
               at_3 = residue_atoms[iat];
               break;
            } 
         }
         residue_atoms = 0; // just to be safe
         r_4->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            if (coot::atom_spec_t(residue_atoms[iat]) == extra_restraints.torsion_restraints[i].atom_4) {
               at_4 = residue_atoms[iat];
               break;
            } 
         }

         if (at_1 && at_2 && at_3 && at_4) {
            int index_1 = -1; 
            int index_2 = -1;
            int index_3 = -1; 
            int index_4 = -1;
            at_1->GetUDData(udd_atom_index_handle, index_1);
            at_2->GetUDData(udd_atom_index_handle, index_2);
            at_3->GetUDData(udd_atom_index_handle, index_3);
            at_4->GetUDData(udd_atom_index_handle, index_4);
            if ((index_1 != -1) && (index_2 != -1) && (index_3 != -1) && (index_4 != -1)) { 
               std::vector<bool> fixed_flags = make_fixed_flags(index_1, index_2, index_3, index_4);
               if (fixed_1) fixed_flags[0] = 1;
               if (fixed_2) fixed_flags[1] = 1;
               if (fixed_3) fixed_flags[2] = 1;
               if (fixed_4) fixed_flags[3] = 1;

               if (true)
                  std::cout << "DEBUG:: adding user-defined torsion restraint with fixed flags: "
                            << "[" << index_1 << " " << coot::atom_spec_t(atom[index_1]) << " " << fixed_flags[0] << "]  " 
                            << "[" << index_2 << " " << coot::atom_spec_t(atom[index_2]) << " " << fixed_flags[1] << "]  " 
                            << "[" << index_3 << " " << coot::atom_spec_t(atom[index_3]) << " " << fixed_flags[2] << "]  " 
                            << "[" << index_4 << " " << coot::atom_spec_t(atom[index_4]) << " " << fixed_flags[3] << "]  " 
                            << std::endl;

               add_user_defined_torsion_restraint(TORSION_RESTRAINT,
                                                  index_1, index_2, index_3, index_4,
                                                  fixed_flags,
                                                  extra_restraints.torsion_restraints[i].torsion_angle,
                                                  extra_restraints.torsion_restraints[i].esd,
                                                  1.2, // dummy value
                                                  extra_restraints.torsion_restraints[i].period);
            }
         } 
      }
   }
}

void
coot::restraints_container_t::add_extra_parallel_plane_restraints(int imol,
                                                                  const extra_restraints_t &extra_restraints,
                                                                  const protein_geometry &geom) {

//    std::cout << "------ in add_extra_parallel_plane_restraints() " << extra_restraints.parallel_plane_restraints.size()
//              << " pp restraints " << std::endl;

   for (unsigned int i=0; i<extra_restraints.parallel_plane_restraints.size(); i++) {
      std::vector<int> plane_1_atom_indices;
      std::vector<int> plane_2_atom_indices;
      const parallel_planes_t &r = extra_restraints.parallel_plane_restraints[i];
      mmdb::Residue *r_1 = util::get_residue(r.plane_1_atoms.res_spec, mol);
      mmdb::Residue *r_2 = util::get_residue(r.plane_2_atoms.res_spec, mol);

      if (false) {
         std::cout << "------ in add_extra_parallel_plane_restraints() extracting 1 " << r.plane_1_atoms.res_spec
                   << std::endl;
         std::cout << "------ in add_extra_parallel_plane_restraints() extracting 2 " << r.plane_2_atoms.res_spec
                   << std::endl;
         std::cout << "------ in add_extra_parallel_plane_restraints() extracting from mol " << mol << std::endl;
         std::cout << "------ in add_extra_parallel_plane_restraints() " << r_1 << " " << r_2 << std::endl;
      }

      if (r_1 && r_2) {

         bool fixed_1 = 0;
         bool fixed_2 = 0;

         // 20131112 OK, so the extra restraints have non-spaced
         // names, we need (at the moment at least) to look up the
         // 4-char names so that we can match the atom names from the
         // PDB.  For that, we need the dictionary for the residue
         // types of the selected residues.
         //
         //
         std::string res_type_1 = r_1->GetResName();
         std::string res_type_2 = r_2->GetResName();
         std::pair<bool, dictionary_residue_restraints_t> dri_1 = geom.get_monomer_restraints(res_type_1, imol);
         std::pair<bool, dictionary_residue_restraints_t> dri_2 = geom.get_monomer_restraints(res_type_2, imol);

         if (dri_1.first && dri_2.first) {


            if (from_residue_vector) {
               for (unsigned int ir=0; ir<residues_vec.size(); ir++) {
                  if (residues_vec[ir].second == r_1) 
                     fixed_1 = residues_vec[ir].first;
                  if (residues_vec[ir].second == r_2) 
                     fixed_2 = residues_vec[ir].first;
               }
            } 
         
            mmdb::PPAtom residue_atoms = 0;
            int n_residue_atoms;

            // add to plane_1_atom_indices
            // 
            r_1->GetAtomTable(residue_atoms, n_residue_atoms);
            for (unsigned int i_rest_at=0; i_rest_at<r.plane_1_atoms.atom_names.size(); i_rest_at++) {
               std::string plane_atom_expanded_name =
                  dri_1.second.atom_name_for_tree_4c(r.plane_1_atoms.atom_names[i_rest_at]);
               for (int iat=0; iat<n_residue_atoms; iat++) {
                  mmdb::Atom *at = residue_atoms[iat];
                  std::string atom_name(at->name);
                  std::string alt_conf(at->altLoc);
                  if (plane_atom_expanded_name == atom_name) {
                     if (r.plane_1_atoms.alt_conf == alt_conf) {
                        int idx = -1;
                        if (at->GetUDData(udd_atom_index_handle, idx) == mmdb::UDDATA_Ok) { 
                           if (idx != -1) {
                              plane_1_atom_indices.push_back(idx);
                              if (0)
                                 std::cout << "adding plane-1 parallel plane atom " << atom_spec_t(at)
                                           << " which has idx " << idx << std::endl;
                           }
                        } else {
                           std::cout << "no udd_atom_index_handle for " <<  atom_spec_t(at) << std::endl;
                        } 
                     }
                  }
               }
            }
            // same for plane_2_atom_indices
            //
            residue_atoms = 0;
            r_2->GetAtomTable(residue_atoms, n_residue_atoms);
            for (unsigned int i_rest_at=0; i_rest_at<r.plane_2_atoms.atom_names.size(); i_rest_at++) {
               std::string plane_atom_expanded_name =
                  dri_2.second.atom_name_for_tree_4c(r.plane_2_atoms.atom_names[i_rest_at]);
               for (int iat=0; iat<n_residue_atoms; iat++) {
                  mmdb::Atom *at = residue_atoms[iat];
                  std::string atom_name(at->name);
                  std::string alt_conf(at->altLoc);
                  // std::cout << "testing :" << plane_atom_expanded_name << ": vs :" << atom_name << ":" << std::endl;
                  if (plane_atom_expanded_name == atom_name) {
                     if (r.plane_2_atoms.alt_conf == alt_conf) {
                        int idx = -1;
                        if (at->GetUDData(udd_atom_index_handle, idx) == mmdb::UDDATA_Ok) { 
                           if (idx != -1) {
                              plane_2_atom_indices.push_back(idx);
                              if (0)
                                 std::cout << "adding plane-2 parallel plane atom " << atom_spec_t(at)
                                           << " which has idx " << idx << std::endl;
                           }
                        } else {
                           std::cout << "no udd_atom_index_handle for " <<  atom_spec_t(at) << std::endl;
                        }
                     }
                  }
               }
            }

            if (plane_1_atom_indices.size() > 3) {
               if (plane_2_atom_indices.size() > 3) {

                  std::vector<bool> fixed_atoms_plane_1 = make_fixed_flags(plane_1_atom_indices);
                  std::vector<bool> fixed_atoms_plane_2 = make_fixed_flags(plane_2_atom_indices);

                  simple_restraint sr(PARALLEL_PLANES_RESTRAINT,
                                      plane_1_atom_indices,
                                      plane_2_atom_indices,
                                      fixed_atoms_plane_1,
                                      fixed_atoms_plane_2,
                                      r.target_angle, r.sigma_combined_planes);

                  restraints_vec.push_back(sr);
                  if (0)
                     std::cout << "after pp restraints with sigma " << sr.sigma << " from " << r.sigma_combined_planes
                               << "   after push restraints_vec is now of size "
                               << restraints_vec.size() << std::endl;
               }
            }
         }
      }
   }
}



// We want to interpolate proSMART restraints from start to final model.
// We have proSMART restraints for both models.
// 
// Let's return a list of extra bond restraint indices that are
// between the bond restraints of this (presumably start) and
// final.
// 
std::vector<std::pair<unsigned int, unsigned int> >
coot::extra_restraints_t::find_pair_indices(const extra_restraints_t &final) const {

   std::vector<std::pair<unsigned int, unsigned int> > v;

   for (unsigned int ibond=0; ibond<bond_restraints.size(); ibond++) { 
      const extra_bond_restraint_t &br_i = bond_restraints[ibond];
      bool found = false;
      for (unsigned int jbond=0; jbond<final.bond_restraints.size(); jbond++) {
         const extra_bond_restraint_t &br_j = final.bond_restraints[jbond];


         // std::cout << "comparing:  " << br_i.atom_1 << "  vs    "
         // << br_j.atom_1 << std::endl;

         if (br_i.atom_1 == br_j.atom_1) {
            if (br_i.atom_2 == br_j.atom_2) {
               std::pair<unsigned int, unsigned int> p(ibond, jbond);
               v.push_back(p);
               found = true;
               break;
            }
         }
         if (br_i.atom_1 == br_j.atom_2) {
            if (br_i.atom_2 == br_j.atom_1) {
               std::pair<unsigned int, unsigned int> p(ibond, jbond);
               v.push_back(p);
               found = true;
               break;
            }
         }
      }
   }
   return v;
}

void
coot::extra_restraints_t::write_interpolated_restraints(const extra_restraints_t &final,
                                                        unsigned int n_path_points,
                                                        std::string file_name_stub) const {

   if (n_path_points <= 2)
      return;

   unsigned int i_end = n_path_points - 1;

   // when path i_path = 0 is the start (this)
   // 
   // when path i_path = i_end = n_path_points - 1 is final
   // 
   std::vector<std::pair<unsigned int, unsigned int> > indices = find_pair_indices(final);

   for (unsigned int i_path=0; i_path<n_path_points; i_path++) {
      std::string file_name = file_name_stub + util::int_to_string(i_path) + ".rst";
      std::ofstream f(file_name.c_str());
      if (f) {
         double frac = double(i_path)/double(i_end);
         for (unsigned int ii=0; ii<indices.size(); ii++) { 
            const std::pair<unsigned int, unsigned int> &pi = indices[ii];
            write_interpolated_restraints(f, final.bond_restraints, frac,
                                          pi.first, pi.second);
         }
      }
      f.close();
   }
}

void
coot::extra_restraints_t::write_interpolated_restraints(std::ofstream &f,
                                                        const std::vector<coot::extra_restraints_t::extra_bond_restraint_t> &final_bond_restraints,
                                                        double frac,
                                                        unsigned int idx_1,
                                                        unsigned int idx_2) const {

   // extra restraints are also written in user-defined-restraints.scm

   const extra_bond_restraint_t &br_1 = bond_restraints[idx_1];
   const extra_bond_restraint_t &br_2 = final_bond_restraints[idx_2];

   double delta = (br_2.bond_dist - br_1.bond_dist);
   double d = br_1.bond_dist + frac * delta;
   double esd = br_1.esd; // simple
                                
   f << "EXTE DIST FIRST ";
   f << "CHAIN ";
   if (br_1.atom_1.chain_id == " " || br_1.atom_1.chain_id == "")
      f << ".";
   else
      f << br_1.atom_1.chain_id;
   f << " RESI ";
   f << util::int_to_string(br_1.atom_1.res_no);
   f << " INS ";
   if (br_1.atom_1.ins_code == " " || br_1.atom_1.ins_code == "")
      f << ".";
   else
      f << br_1.atom_1.ins_code;
   f << " ATOM ";
   f << br_1.atom_1.atom_name;
   f << " ";

   f << " SECOND ";
   
   f << "CHAIN ";
   if (br_1.atom_2.chain_id == " " || br_1.atom_2.chain_id == "")
      f << ".";
   else
      f << br_1.atom_2.chain_id;
   f << " RESI ";
   f << util::int_to_string(br_1.atom_2.res_no);
   f << " INS ";
   if (br_1.atom_2.ins_code == " " || br_1.atom_2.ins_code == "")
      f << ".";
   else
      f << br_1.atom_2.ins_code;
   f << " ATOM ";
   f << br_1.atom_2.atom_name;
   f << " ";

   f << " VALUE ";
   f << d;
   f << " SIGMA ";
   f << esd;
   f << "\n";

} 


void
coot::extra_restraints_t::write_interpolated_models_and_restraints(const extra_restraints_t &final,
                                                                   mmdb::Manager *mol_1, // corresponds to this
                                                                   mmdb::Manager *mol_2, // corresponds to final
                                                                   unsigned int n_path_points,
                                                                   std::string file_name_stub) const {

   if (n_path_points <= 2)
      return;

   // when path i_path = 0 is the start (this)
   // 
   // when path i_path = i_end = n_path_points - 1 is final
   //
   if (mol_1) { 
      if (mol_2) {

         mmdb::Manager *mol_running = new mmdb::Manager;
         mol_running->Copy(mol_1, mmdb::MMDBFCM_All);
         
         std::map<mmdb::Atom *, clipper::Coord_orth> matching_atoms_1 = position_point_map(mol_running, mol_1);
         std::map<mmdb::Atom *, clipper::Coord_orth> matching_atoms_2 = position_point_map(mol_running, mol_2);
         
         logger.log(log_t::INFO, logging::ltw("found "), logging::ltw(static_cast<unsigned int>(matching_atoms_1.size())), logging::ltw(" (1) matching atoms"));
         // std::cout << "INFO:: found " << matching_atoms_1.size() << " (1) matching atoms " << std::endl;
         logger.log(log_t::INFO, logging::ltw("found "), logging::ltw(static_cast<unsigned int>(matching_atoms_2.size())), logging::ltw(" (2) matching atoms"));
         // std::cout << "INFO:: found " << matching_atoms_2.size() << " (2) matching atoms " << std::endl;

         if (matching_atoms_1.size() && matching_atoms_2.size()) {
            write_interpolated_restraints(final, n_path_points, file_name_stub);
            write_interpolated_models(mol_running, matching_atoms_1, matching_atoms_2, n_path_points, file_name_stub);
         }
      }
   }
}

std::map<mmdb::Atom *, clipper::Coord_orth>
coot::extra_restraints_t::position_point_map(mmdb::Manager *mol_running,
                                             mmdb::Manager *mol_ref) const {

   std::map<mmdb::Atom *, clipper::Coord_orth> matching_atoms;
   
   if (mol_running) { 
      if (mol_ref) { 
         
         for(int imod = 1; imod<=mol_running->GetNumberOfModels(); imod++) {
            mmdb::Model *model_1_p = mol_running->GetModel(imod);
            if (model_1_p) {

               for(int jmod = 1; jmod<=mol_ref->GetNumberOfModels(); jmod++) {
                  mmdb::Model *model_2_p = mol_ref->GetModel(jmod);
                  if (model_2_p) {

                     mmdb::Chain *chain_1_p;
                     int n_chains_1 = model_1_p->GetNumberOfChains();
                     for (int ichain=0; ichain<n_chains_1; ichain++) {
                        chain_1_p = model_1_p->GetChain(ichain);
                        if (chain_1_p) {
                           std::string chain_1_id = chain_1_p->GetChainID();

                           mmdb::Chain *chain_2_p;
                           int n_chains_2 = model_2_p->GetNumberOfChains();
                           for (int jchain=0; jchain<n_chains_2; jchain++) {
                              chain_2_p = model_2_p->GetChain(jchain);
                              if (chain_2_p) {
                                 std::string chain_2_id = chain_2_p->GetChainID();
                                 if (chain_1_id == chain_2_id) { 

                                    int nres_1 = chain_1_p->GetNumberOfResidues();
                                    mmdb::Residue *residue_1_p;
                                    mmdb::Atom *at_1_p;
                                    for (int ires=0; ires<nres_1; ires++) {
                                       residue_1_p = chain_1_p->GetResidue(ires);
                                       int res_no_1 = residue_1_p->GetSeqNum();
                                       std::string ins_code_1 = residue_1_p->GetInsCode();

                                       int nres_2 = chain_2_p->GetNumberOfResidues();
                                       mmdb::Residue *residue_2_p;
                                       mmdb::Atom *at_2_p;
                                       for (int jres=0; jres<nres_2; jres++) {
                                          residue_2_p = chain_2_p->GetResidue(jres);
                                          int res_no_2 = residue_2_p->GetSeqNum();
                                          std::string ins_code_2 = residue_2_p->GetInsCode();

                                          if (res_no_2 == res_no_1) {
                                             if (ins_code_2 == ins_code_1) { 

                                                int n_atoms_1 = residue_1_p->GetNumberOfAtoms();
                                                for (int iat=0; iat<n_atoms_1; iat++) {
                                                   at_1_p = residue_1_p->GetAtom(iat);
                                                   std::string atom_name_1(at_1_p->name);
                                                   std::string alt_conf_1(at_1_p->altLoc);

                                                   int n_atoms_2 = residue_2_p->GetNumberOfAtoms();
                                                   for (int jat=0; jat<n_atoms_2; jat++) {
                                                      at_2_p = residue_2_p->GetAtom(jat);
                                                      std::string atom_name_2(at_2_p->name);
                                                      std::string alt_conf_2(at_2_p->altLoc);
                                                      
                                                      if (atom_name_2 == atom_name_1) {
                                                         if (alt_conf_2 == alt_conf_1) {
                                                            
                                                            clipper::Coord_orth co = coot::co(at_2_p);
                                                            matching_atoms[at_1_p] = co;
                                                            break;
                                                         }
                                                      }
                                                   }
                                                }
                                                break; // there won't be another residue that matches
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return matching_atoms;
} 

void
coot::extra_restraints_t::write_interpolated_models(mmdb::Manager *mol_running,
                                                    const std::map<mmdb::Atom *, clipper::Coord_orth> &matching_atoms_1,
                                                    const std::map<mmdb::Atom *, clipper::Coord_orth> &matching_atoms_2,
                                                    unsigned int n_path_points,
                                                    std::string file_name_stub) const {

   unsigned int i_end = n_path_points - 1;

   logger.log(log_t::INFO, logging::ltw("number of interpolation points: "), logging::ltw(n_path_points));
   // std::cout << "INFO:: number of interpolation points: " << n_path_points << std::endl;
   for (unsigned int i_path=0; i_path<n_path_points; i_path++) {
      std::string file_name = file_name_stub + util::int_to_string(i_path) + ".pdb";
      double frac = double(i_path)/double(i_end);

      for(int imod = 1; imod<=mol_running->GetNumberOfModels(); imod++) {
         mmdb::Model *model_p = mol_running->GetModel(imod);
         if (model_p) { 
            mmdb::Chain *chain_p;
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               chain_p = model_p->GetChain(ichain);
               int nres = chain_p->GetNumberOfResidues();
               mmdb::Residue *residue_p;
               mmdb::Atom *at;
               for (int ires=0; ires<nres; ires++) { 
                  residue_p = chain_p->GetResidue(ires);
                  int n_atoms = residue_p->GetNumberOfAtoms();
                  for (int iat=0; iat<n_atoms; iat++) {
                     at = residue_p->GetAtom(iat);
                     atom_spec_t spec(at);
                     std::map<mmdb::Atom *, clipper::Coord_orth>::const_iterator it_1;
                     std::map<mmdb::Atom *, clipper::Coord_orth>::const_iterator it_2;
                     it_1 = matching_atoms_1.find(at);
                     it_2 = matching_atoms_2.find(at);
                     if (it_1 != matching_atoms_1.end()) {
                        if (it_2 != matching_atoms_2.end()) {
                           const clipper::Coord_orth &pt_1 = it_1->second;
                           const clipper::Coord_orth &pt_2 = it_2->second;
                           clipper::Coord_orth pt(pt_1 + (pt_2 - pt_1) * frac);
                           at->x = pt.x();
                           at->y = pt.y();
                           at->z = pt.z();
                        } else {
                           std::cout << "failed to find spec for it_2 " << spec << std::endl;
                        } 
                     } else {
                        std::cout << "failed to find spec for it_1 " << spec << std::endl;
                     }
                  }
               }
            }
         }
      }
      mol_running->FinishStructEdit();
      mol_running->WritePDBASCII(file_name.c_str());
   }
}
