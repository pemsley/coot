#ifdef __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#else
#define __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#undef __GNU_LIBRARY__
#endif

#include "utils/coot-utils.hh"
#include "protein-geometry.hh"

int
compare_dictionaries(const std::string &type,
		     const std::string &file_name_1,
		     const std::string &file_name_2,
		     double bond_length_tolerance,
		     double bond_esd_tolerance,
		     double angle_tolerance,
		     double angle_esd_tolerance,
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
	 bool compare_status = r1.second.compare(r2.second,
						 bond_length_tolerance,
						 bond_esd_tolerance,
						 angle_tolerance,
						 angle_esd_tolerance,
						 quiet);
	 status = !compare_status; // invert for unix return value (0 happy)
      }
   }
   return status;
}

void print_help(std::string cmd) {
   std::cout << "Usage: " << cmd << " "
      "--help     this help\n" << 
      "--quiet    do not report satisfatory matches\n" << 
      "--type     residue tupe to match\n" << 
      "--dict-1   file with reference dictionary\n"
      "--dict-2   file with comparison dictionary\n" << 
      "--bond-length-tol\n" << 
      "--bond-length-esd-tol" << 
      "--angle-tol" << 
      "--angle-esd-tol" << 
      std::endl;
   
} 

int main(int argc, char **argv) {

   int status = 0;
   bool quiet = false;

   if (argc < 4) {
      print_help(argv[0]);
   } else {

      std::string type;
      std::string file_name_1;
      std::string file_name_2;

      double bond_length_tolerance = 0.01;
      double bond_esd_tolerance    = 0.005;
      double angle_tolerance       = 0.3;
      double angle_esd_tolerance   = 0.15;
      
      const char *optstr = "q";
      struct option long_options[] = {
	 {"help",   0, 0, 0},
	 {"quiet",   0, 0, 0},
	 {"type",    1, 0, 0},
	 {"dict-1",  1, 0, 0},
	 {"dict-2",  1, 0, 0},
	 {"bond-length-tol",  1, 0, 0},
	 {"bond-length-esd-tol",  1, 0, 0},
	 {"angle-tol",  1, 0, 0},
	 {"angle-esd-tol",  1, 0, 0},
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
	       if (arg_str == "bond-length-tol") { 
		  try {
		     bond_length_tolerance = coot::util::string_to_float(optarg);
		  }
		  catch (const std::runtime_error &rte) {
		     std::cout << "bad number " << optarg << std::endl;
		  }
	       }
	       if (arg_str == "bond-length-esd-tol") { 
		  try {
		     bond_esd_tolerance = coot::util::string_to_float(optarg);
		  }
		  catch (const std::runtime_error &rte) {
		     std::cout << "bad number " << optarg << std::endl;
		  }
	       }
	       if (arg_str == "angle-tol") { 
		  try {
		     angle_tolerance = coot::util::string_to_float(optarg);
		  }
		  catch (const std::runtime_error &rte) {
		     std::cout << "bad number " << optarg << std::endl;
		  }
	       }
	       if (arg_str == "angle-esd-tol") { 
		  try {
		     angle_esd_tolerance = coot::util::string_to_float(optarg);
		  }
		  catch (const std::runtime_error &rte) {
		     std::cout << "bad number " << optarg << std::endl;
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
					     quiet);
	    }
	 }
      }
   }
   return status;
}
