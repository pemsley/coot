
#include <vector>
#include <string>
#include <clipper/core/xmap.h>
#include <clipper/ccp4/ccp4_map_io.h>
#include "tube-finder.hh"

int main(int argc, char **argv) {

   int status = 0;

   if (argc> 1) {
      std::string map_file_name = argv[1];
      clipper::CCP4MAPfile file;
      clipper::Xmap<float> xmap;
      file.open_read(map_file_name);
      file.import_xmap(xmap);
      coot::tube_finder_t tf(xmap);
      std::vector<clipper::Coord_orth> positions = tf.get_positions();
      for (unsigned int i=0; i<positions.size(); i++) {
         std::cout << i << " " << positions[i].format() << std::endl;
      }
   }

   return status;
}
