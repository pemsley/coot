
#include <clipper/ccp4/ccp4_map_io.h>

#include "../utils/coot-utils.hh"
#include "trace.hh"

int main(int argc, char **argv) {

   bool debug = false;

   std::string map_file_name = "trace-test.map";

   if (argc > 1)
      map_file_name = argv[1];
   
   if (coot::file_exists(map_file_name)) {

      try { 
	 clipper::CCP4MAPfile file;
	 clipper::Xmap<float> xmap;
	 file.open_read(map_file_name);
	 file.import_xmap(xmap);
	 file.close_read();

	 coot::trace t(xmap);


	 // test from a pdb file
	 std::string test_pdb_file_name = "test-trace-template.pdbxx";
	 // test with a null moll or flood mol
	 if (coot::file_exists(test_pdb_file_name)) {

	    mmdb::Manager *mol = new mmdb::Manager;
	    mmdb::ERROR_CODE err = mol->ReadCoorFile(test_pdb_file_name.c_str());
	    if (! err) {
	       std::cout << "running with test mol: " << std::endl;
	       t.test_model(mol);
	    }
	 } else {

	    t.action();

	 } 
      }
      
      // problem reading the map, perhaps?
      // 
      catch (const clipper::Message_fatal &mess) {
	 std::cout << "ERROR:: " << mess.text() << std::endl;
      } 
   }
   return 0;
}
