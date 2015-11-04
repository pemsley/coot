
#include <clipper/ccp4/ccp4_map_io.h>

#include "../utils/coot-utils.hh"
#include "trace.hh"

int main(int argc, char **argv) {

   bool debug = false;

   std::string map_file_name = "trace-test.map";

   if (argc > 1)
      map_file_name = argv[1];
   
   if (coot::file_exists(map_file_name)) { 
      clipper::CCP4MAPfile file;
      clipper::Xmap<float> xmap;
      file.open_read(map_file_name);
      file.import_xmap(xmap);
      file.close_read();

      if (debug) {
	 clipper::CCP4MAPfile mapout;
	 mapout.open_write("duplicate.map");
	 mapout.export_xmap(xmap);
	 mapout.close_write();
      }
      
      coot::trace t(xmap);
      t.action();
   }

   return 0;

}
