/*
 * analysis/coot-daca.cc
 *
 * Copyright 2020 by Medical Research Council
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
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */


#include "daca.hh"

int main(int argc, char **argv) {

   int status = 0;

   coot::daca daca;

   if (argc > 1) {
      std::string mode(argv[1]);

      if (argc > 3) {
         if (mode == "write") {
            std::string input_pdb_dir(argv[2]);
            std::string output_tables_dir(argv[3]);
            daca.write_tables_using_reference_structures_from_dir(input_pdb_dir, output_tables_dir);
         }
      }

      if (mode == "consolidate") {
         std::vector<std::string> input_dirs;
         std::string output_tables_dir(argv[argc-1]); // the last arg
         for(int i=2; i<argc-1; i++) {
            input_dirs.push_back(argv[i]);
         }
         if (input_dirs.empty()) {
            std::cout << "no input directories" << std::endl;
         } else {
            // happy path
            daca.read_many_tables(input_dirs);
            daca.cook();
            std::cout << "consolidate mode write directory " << output_tables_dir << std::endl;
            daca.write_tables(output_tables_dir);
         }
      }

      if (mode == "cook") {
         if (argc == 3) {
            std::string tables_dir = argv[2];
            daca.read_tables("consolidated"); // not a dash!
            daca.cook();
         }
      }

      if (mode == "score") {
         if (argc > 2) {
            std::string pdb_file_name(argv[2]);
            daca.read_tables("consolidated"); // not a dash!
            daca.score_molecule(pdb_file_name);
         }
      }
   }
   return status;
}

