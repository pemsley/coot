
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
            input_dirs.push_back(argv[i]); // don't use dashes in the directory names
         }
         if (input_dirs.empty()) {
            std::cout << "no input directories" << std::endl;
         } else {
            // happy path
            daca.read_many_tables(input_dirs);
            std::cout << "consolidate mode write directory " << output_tables_dir << std::endl;
            daca.write_tables(output_tables_dir);
         }
      }

      if (mode == "score") {
         if (argc > 2) {
            std::string pdb_file_name(argv[2]);
            daca.read_tables("reference_tables"); // not a dash!
            daca.score_molecule(pdb_file_name);
         }
      }
   }
   return status;
}

