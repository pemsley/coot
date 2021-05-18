/* coot-utils/residue-and-atom-specs.hh
 * 
 * Copyright 2011, 2012 by The University of Oxford
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

#include "glyco-torsions.hh"

int main(int argc, char **argv) {

   if (argc <= 2) {
      std::cout << "Usage: " << argv[0] << " template-pdb link-type new-comp-id (generate mode)\n";
      std::cout << "   or: " << argv[0] << " test link-type new-comp-id \n";
   } else {

      mmdb::InitMatType();

      if (std::string(argv[1]) == "test") {

         // use the link-by-torsion reference files to make a new residue (i.e. to test
         // the function)

         if (argc > 3) {
            std::string link_type = argv[2];
            std::string new_residue_type = argv[3];

            std::string file_name = "link-by-torsion-to-pyranose-core-" + link_type + ".tab";
            if (! coot::file_exists(file_name)) {
               std::cout << "ERROR: file not found " << file_name << std::endl;
            } else {
               coot::link_by_torsion_t l(file_name);
               l.new_residue_type = new_residue_type;

               // Get a base residue
               mmdb::Manager *mol = new mmdb::Manager;
               std::string pdb_file_name =
                  "pdb-templates/pyranose-pyranose-via-" + link_type + ".pdb";

               if (link_type == "NAG-ASN")
                  pdb_file_name = "pdb-templates/ASN-NAG-via-NAG-ASN.pdb";

               if (! coot::file_exists(pdb_file_name)) {
                  std::cout << "ERROR:: pdb file " << pdb_file_name << " does not exist"
                            << std::endl;
               } else {
                  int err = mol->ReadPDBASCII(pdb_file_name.c_str());
                  int cr = mol->CrystReady();
                  std::cout << "INFO:: read " << pdb_file_name << " gives code " << err
                            << " and cryst status: " << cr << std::endl;
                  mmdb::Residue *base_residue_p = coot::util::get_first_residue(mol);
                  if (! base_residue_p) {
                     std::cout << "ERROR:: no base residue " << link_type << " in " << file_name
                               << std::endl;
                  } else {
                     std::string decor_file_name = new_residue_type + "-decorations.tab";
                     if (! coot::file_exists(decor_file_name)) {
                        std::cout << "No file " << decor_file_name << std::endl;
                     } else {
                        coot::link_by_torsion_t decor(decor_file_name);
                        if (! decor.filled()) {
                           std::cout << "Decorations not filled from " << decor_file_name
                                     << std::endl;
                        } else {
                           l.add(decor);
                           l.set_new_residue_number(1);
                           mmdb::Residue *r = l.make_residue(base_residue_p);
                           if (r) {
                              mmdb::Manager *mol = coot::util::create_mmdbmanager_from_residue(r);
                              std::string output_pdb_file_name = "output-" + new_residue_type + ".pdb";
                              mol->WritePDBASCII(output_pdb_file_name.c_str());
                           }
                        }
                     }
                  }
               }
            }
         }

      } else {

         // make the link-by-torsion reference files

         std::cout << "--- generate mode --- " << std::endl;

         std::string file_name = argv[1];
         std::string link_type = argv[2]; // e.g. "ALPHA1-6";
         std::string new_residue_type = "MAN";
         if (argc == 4)
            new_residue_type = argv[3];

         mmdb::Manager *mol = new mmdb::Manager;
         int status = mol->ReadPDBASCII(file_name.c_str());
         if (status != mmdb::Error_NoError) {
            std::cout << "ERROR:: on reading " << file_name << std::endl;
         } else {
            int cr = mol->CrystReady();
            std::cout << "INFO:: read " << file_name << " gives code " << status
                      << " and cryst status: " << cr << "\n";
            std::cout << "INFO:: read " << file_name << " OK " << std::endl;
            std::pair<mmdb::Residue *, mmdb::Residue *> p = coot::link_by_torsion_t::get_residue_pair(mol);
            if (! p.first || !p.second) {
               std::cout << "Failed to get residue pair from " << file_name << std::endl;
            } else {
               // generate link torsions, write link torsions
               std::cout << "INFO:: Getting names for link_type " << link_type << std::endl;
               coot::link_by_torsion_base_t to_core = coot::get_names_for_link_type(link_type);
               if (! to_core.filled()) {
                  std::cout << "ERROR:: failed to get names for " << link_type << std::endl;
               } else {
                  coot::link_by_torsion_t l_to_core(to_core, p.first, p.second);
                  if (! l_to_core.filled()) {
                     std::cout << "ERROR:: failed to fill link to core params for "
                               << link_type << std::endl;
                  } else {
                     std::string file_name = "link-by-torsion-to-pyranose-core-" + link_type + ".tab";
                     l_to_core.write(file_name);

                     // and the decorations on the core:
                     coot::link_by_torsion_base_t decor = coot::get_decorations(new_residue_type);
                     coot::link_by_torsion_t l_decor(decor, p.first, p.second);
                     if (! l_decor.filled()) {
                        std::cout << "ERROR:: No decorations for "
                                  << new_residue_type << ", i.e. needs coot::" <<  new_residue_type
                                  << "_decorations() function to be written" << std::endl;
                     } else {
                        std::string file_name =  new_residue_type + "-decorations.tab";
                        l_decor.write(file_name);
                     }
                  }
               }
            }
         }
      }
   }
   return 0;
}
//
