#include <string>
#include <iostream>
#include "segmap.hh"

int main(int argc, char **argv) {

   bool done = false;

   if (argc == 3) {
      // generate the stats from the sampled maps.
      std::string a1(argv[1]);
      if (a1 == "proc") {
	 std::string map_file_name(argv[2]);
	 clipper::CCP4MAPfile file;
	 try {
	    file.open_read(map_file_name);
	    clipper::Xmap<float> xmap;
	    file.import_xmap(xmap);
	    coot::segmap s(xmap);
	    s.proc();
	    done = true;
	 }
	 catch (const clipper::Message_base &exc) {
	    std::cout << "WARNING:: failed to open " << map_file_name << std::endl;
	 }
      }
   }

   return 0;
}
