
#include <clipper/ccp4/ccp4_map_io.h>
#include "occlusion.hh"
#include "coot-map-utils.hh"

#include "density-contour/CIsoSurface.h"



int main(int argc, char **argv) {

   int status = 0;

   std::string file_name = "test.map";
   if (argc == 2)
      file_name = argv[1];

   try {
      clipper::CCP4MAPfile file;
      clipper::Xmap<float> xmap;
      file.open_read(file_name);
      file.import_xmap(xmap);

      // now contour the xmap

      CIsoSurface<float> my_isosurface;
      coot::CartesianPairInfo v;
      float contour_level = 0.8;
      float dy_radius = 6.0;
      coot::Cartesian centre(30, 30, 30);
      int isample_step = 1;
      bool is_em_map = false;

      v = my_isosurface.GenerateSurface_from_Xmap(xmap,
                                                  contour_level,
                                                  dy_radius, centre,
                                                  isample_step,
                                                  0,1,1,
                                                  is_em_map);

   }
   catch (const clipper::Message_base &exc) {
      std::cout << "WARNING:: failed to open " << file_name << std::endl;
   }





   return status;

}
