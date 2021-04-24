//
#include <cstring>
#include <string.h>
#include "coot-hole.hh"

int main(int argc, char **argv) {

   if (argc > 1) {
      std::string file_name = argv[1];
      int  error_count;
      char error_buf[500];
      mmdb::InitMatType();
      mmdb::Manager *m = new mmdb::Manager;
      mmdb::ERROR_CODE err = m->ReadCoorFile(file_name.c_str());
      if (err) {
         std::cout << "There was an error reading " << file_name.c_str() << ".\n";
         std::cout << "ERROR " << err << " READ: "
                   << mmdb::GetErrorDescription(err) << std::endl;
         m->GetInputBuffer(error_buf, error_count);
         if (error_count >= 0) {
            std::cout << " LINE #" << error_count << "\n " << error_buf << std::endl;
         }
      } else {

         // clipper::Coord_orth from_pt(4.88335275650024, -0.676153182983398, 17.5295238494873);
         // clipper::Coord_orth   to_pt(2.77147531509399, -1.79189586639404, -11.4192752838135);

         // 1alz
         clipper::Coord_orth from_pt(21.5, 14.5, 11.5);
         clipper::Coord_orth   to_pt(21.5, 52.5, 11.5);

         // 1c4d
         //from_pt = clipper::Coord_orth(-0.36,  21.2, 7.74);
         //to_pt   = clipper::Coord_orth(-0.36, -4.48, 6.80);

         // 30gc
         // from_pt = clipper::Coord_orth(0.0, 0,  38);
         // to_pt   = clipper::Coord_orth(0.0, 0, -23);

         // energy lib:
         coot::protein_geometry geom;
         geom.init_standard();
         geom.try_dynamic_add("FOR", 1);
         geom.try_dynamic_add("DLE", 2);
         geom.try_dynamic_add("DVA", 2);
         geom.try_dynamic_add("ETA", 2);
         geom.try_dynamic_add("EOH", 2);

         coot::hole hole(m, from_pt, to_pt, geom);
         hole.generate();
      }
   }
   return 0;
}
