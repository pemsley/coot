
#include "protein-geometry.hh"

int main(int argc, char **argv) {

   int status = 0;

   if (argc < 4) {
      std::cout << "Usage: " << argv[0] << " type dict-file-name-1 dict-file-name2"
		<< std::endl;
   } else {
      std::string type = argv[1];
      std::string file_name_1 = argv[2];
      std::string file_name_2 = argv[3];
      std::string dir = ".";

      coot::protein_geometry pg_1;
      coot::protein_geometry pg_2;

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
	    status = r1.second.compare(r2.second);
	 }
      }
   }
   return status;
}
