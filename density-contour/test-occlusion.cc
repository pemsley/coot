
#include <chrono>
#include <clipper/ccp4/ccp4_map_io.h>
#include "coot-utils/coot-map-utils.hh"
#include "occlusion.hh"

#include "density-contour/CIsoSurface.h"



int main(int argc, char **argv) {

   int status = 0;

   std::string file_name = "test.map";
   if (argc == 2)
      file_name = argv[1];

   try {
      std::cout << "Reading " << file_name << std::endl;

      auto tp_0 = std::chrono::high_resolution_clock::now();
      clipper::CCP4MAPfile file;
      clipper::Xmap<float> xmap;
      file.open_read(file_name);
      file.import_xmap(xmap);

      // now contour the xmap

      auto tp_1 = std::chrono::high_resolution_clock::now();
      CIsoSurface<float> my_isosurface;
      coot::CartesianPairInfo v;
      float contour_level = 0.8;
      float dy_radius = 18.0;
      coot::Cartesian centre(58, 0, 14);
      int isample_step = 1;
      bool is_em_map = false;

      coot::density_contour_triangles_container_t tri_con;
      tri_con = my_isosurface.GenerateTriangles_from_Xmap(xmap,
                                                          contour_level,
                                                          dy_radius, centre,
                                                          isample_step);
      auto tp_2 = std::chrono::high_resolution_clock::now();

      // tri_con contains points, normals and triangle indices
      std::cout << "Got " << tri_con.points.size() << " points" << std::endl;
      std::cout << "Got " << tri_con.normals.size() << " normals" << std::endl;
      std::cout << "Got " << tri_con.point_indices.size() << " triangles" << std::endl;

      // convert (for now)
      auto tp_3 = std::chrono::high_resolution_clock::now();
      std::vector<coot::augmented_position> positions;
      positions.resize(tri_con.points.size());
      unsigned int n = tri_con.points.size();
      for (unsigned int i=0; i<n; i++) {
         const clipper::Coord_orth pos  = tri_con.points[i];
         const clipper::Coord_orth norm = tri_con.normals[i];
         positions[i] = coot::augmented_position(pos, norm);
      }
      auto tp_4 = std::chrono::high_resolution_clock::now();

      // calculate occlusions
      coot::set_occlusions(positions);
      auto tp_5 = std::chrono::high_resolution_clock::now();

      auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
      auto d21 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_2 - tp_1).count();
      auto d32 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_2).count();
      auto d43 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_4 - tp_3).count();
      auto d54 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_5 - tp_4).count();
      std::cout << "import " << d10 <<  " tri_cons " << d21 << " screen: " << d32 << " convert " << d43 << " set_occlusions " << d54
                << " ms " << std::endl;

   }
   catch (const clipper::Message_base &exc) {
      std::cout << "WARNING:: failed to open " << file_name << std::endl;
   }



   return status;

}
