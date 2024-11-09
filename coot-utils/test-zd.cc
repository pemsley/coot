
#include <clipper/ccp4/ccp4_map_io.h>
#include "coot-map-utils.hh"

int main(int argc, char **argv) {

   auto make_difference_maps = [] (const clipper::Xmap<float> &xmap_zde,
                           std::vector<clipper::Xmap<float> > &xmaps) {

      const clipper::Xmap<float> &xmap_0 = xmaps[0];
      const clipper::Xmap<float> &xmap_1 = xmaps[1];
      for(unsigned int ii=0; ii<xmaps.size(); ii++) {
         const clipper::Xmap<float> &xmap_frame = xmaps[ii];
         std::pair<clipper::Xmap<float>, float> dm =
            coot::util::difference_map(xmap_frame, xmap_0, 1.0);
            int idx = ii + 1; // ii is zero-indexed.
            std::string file_name = std::string("frame-") + std::to_string(idx) + std::string(".map");

            std::cout << "############### map " << file_name << " has rmsd " << dm.second << std::endl;
            clipper::CCP4MAPfile file;
            file.open_write(file_name);
            file.export_xmap(dm.first);
      }

   };

   int status = 0;

   if (argc > 3) {
      std::vector<clipper::Xmap<float> > xmaps;
      std::vector<clipper::Xmap<float> *> xmap_ps;
      clipper::Xmap<float> xmap_mask;

      try {
         clipper::CCP4MAPfile file;
         file.open_read(argv[1]);
         file.import_xmap(xmap_mask);

         for (int i=2; i<argc; i++) {
            clipper::CCP4MAPfile file;
            clipper::Xmap<float> xmap;
            file.open_read(argv[i]);
            file.import_xmap(xmap);
            xmaps.push_back(xmap);
         }
         for(unsigned int ii=0; ii<xmaps.size(); ii++)
            xmap_ps.push_back(&xmaps[ii]);
         clipper::Xmap<float> xmap = coot::util::real_space_zero_dose_extrapolation(xmap_ps, xmap_mask);
         clipper::CCP4MAPfile outmapfile;
         outmapfile.open_write("zde.map");
         outmapfile.export_xmap(xmap);
         outmapfile.close_write();

         make_difference_maps(xmap, xmaps);
      }
      catch (const std::runtime_error &rte) {
         std::cout << "runtime_error " << rte.what() << std::endl;
      }


   }

   return status;
}
