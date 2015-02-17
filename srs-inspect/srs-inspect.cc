
#include <sys/types.h>
#include <sys/stat.h> 
#include <unistd.h>

#include <iostream>
#include <string>

#include "ccp4srs/ccp4srs_manager.h"


    
bool
file_exists(const std::string &filename) {

   struct stat s;
   int fstat = stat(filename.c_str(), &s);
   if ( fstat == -1 ) { // file not exist
      return 0;
   } else {
      return 1;
   }
}


bool
exists_in_refmac(const std::string &id,
		 const std::string &monomer_dir) {

   bool exists = false;
   if (id.length()) {
      char c_sub_dir = tolower(id[0]);
      std::string full_path = monomer_dir;
      full_path += "/";
      full_path += c_sub_dir;
      full_path += "/";
      full_path +=  id;
      full_path +=  ".cif";

      if (file_exists(full_path)) { 
	 // std::cout << "checking " << full_path << " exists "  << std::endl;
	 exists = true;
      } else { 
	 // std::cout << "         " << full_path << " does not exist"  << std::endl;
      } 
      
   } 
   return exists;

} 

int main(int argc, char **argv) {

   int status = 0;
   std::string srs_dir;
   std::string refmac_monomer_dir = "monomers";

   if (argc > 1)
      refmac_monomer_dir = argv[1];

   if (argc > 1)
      srs_dir = argv[1];

   ccp4srs::Manager *srs = new ccp4srs::Manager();
   int rc = srs->loadIndex(srs_dir.c_str());
   if (rc != ccp4srs::CCP4SRS_Ok)  {
      std::cout << "SRS files not found: " << srs_dir << std::endl;
      status = 1;
      
   } else {

      int l = srs->n_entries();
      for (int i=0; i<l ;i++)  {
	 ccp4srs::Monomer  *Monomer = srs->getMonomer(i, NULL);
	 if (Monomer)  {
	    // std::cout << "   " << Monomer->ID() << std::endl;
	    std::string id = Monomer->ID();

	    if (id.length()) { 

	       if (! exists_in_refmac(id, refmac_monomer_dir)) {
		  // std::cout << "   Missing in Refmac: " << id << std::endl;
		  std::cout << id << "\n";
	       }
	    }
	 }
      }
   }

   delete srs;
   return status;
}
