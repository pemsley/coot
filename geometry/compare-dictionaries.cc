#ifdef __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#else
#define __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#undef __GNU_LIBRARY__
#endif

#include "protein-geometry.hh"

int
compare_dictionaries(const std::string &type,
		     const std::string &file_name_1,
		     const std::string &file_name_2,
		     bool quiet) {

   int status = 0;
   std::string dir = ".";

   coot::protein_geometry pg_1;
   coot::protein_geometry pg_2;

   pg_1.set_verbose(false);
   pg_2.set_verbose(false);

   pg_1.init_refmac_mon_lib(file_name_1, 0);
   pg_2.init_refmac_mon_lib(file_name_2, 0);

   std::pair<bool, coot::dictionary_residue_restraints_t> r1 = 
      pg_1.get_monomer_restraints(type);
   std::pair<bool, coot::dictionary_residue_restraints_t> r2 = 
      pg_2.get_monomer_restraints(type);

   if (!r1.first) { 
      std::cout << "Failed to find restraints for type " << type << " in "
		<< file_name_1 << std::endl;
      status = 1;
   } else {
      if (!r2.first) { 
	 std::cout << "Failed to find restraints for type " << type << " in "
		   << file_name_2 << std::endl;
	 status = 1;
      } else {
	 // Happy path
	 //
	 // if compare_status is true, they matched.
	 bool compare_status = r1.second.compare(r2.second, quiet);
	 status = !compare_status; // invert for unix return value (0 happy)
      }
   }
   return status;
} 

int main(int argc, char **argv) {

   int status = 0;
   bool quiet = false;

   if (argc < 4) {
      std::cout << "Usage: " << argv[0] << " type dict-file-name-1 dict-file-name2"
		<< std::endl;
   } else {

      std::string type;
      std::string file_name_1;
      std::string file_name_2;
      
      const char *optstr = "q";
      struct option long_options[] = {
	 {"quiet",   0, 0, 0},
	 {"type",    1, 0, 0},
	 {"dict-1",  1, 0, 0},
	 {"dict-2",  1, 0, 0},
	 {0, 0, 0, 0}
      };

      int ch;
      int option_index = 0;
      while ( -1 != 
	      (ch = getopt_long(argc, argv, optstr, long_options, &option_index))) { 

	 switch(ch) {
	    
	 case 0:
	    if (optarg) {
	       std::string arg_str = long_options[option_index].name;
	       if (arg_str == "type")
		  type = optarg;
	       if (arg_str == "dict-1")
		  file_name_1 = optarg;
	       if (arg_str == "dict-2")
		  file_name_2 = optarg;
	    } else {
	       std::string arg_str = long_options[option_index].name;
	       if (arg_str == "quiet") { 
		  quiet = true;
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
	       status = compare_dictionaries(type, file_name_1, file_name_2, quiet);
	    }
	 }
      }
   }
   return status;
}
