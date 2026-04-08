
#include <string>
#include <filesystem>
#include <clipper/ccp4/ccp4_mtz_io.h>
#include <clipper/ccp4/ccp4_map_io.h>
#include "born.hh"

void
coot::born(const std::string &mtz_file_name, const std::string &f_col, const std::string &phi_col, const std::string &map_name) {

  auto make_import_datanames = [] (const std::string &f_col, const std::string &phi_col) {
     return "/*/*/[" + f_col + " " + phi_col + "]";
  };

   bool debug = true;
   std::filesystem::path mtz_file_path(mtz_file_name);
   if (std::filesystem::exists(mtz_file_path)) {
      try {
         clipper::CCP4MTZfile mtzin;
         mtzin.open_read(mtz_file_name);
         clipper::Spacegroup spgr = mtzin.spacegroup();
         clipper::Cell       cell = mtzin.cell();
         clipper::Resolution reso = mtzin.resolution();
         clipper::HKL_info hkls, hklp;
         typedef clipper::HKL_data_base::HKL_reference_index HRI;

         clipper::HKL_data<clipper::data32::F_phi> fphi(hkls);
         std::string datanames = make_import_datanames(f_col, phi_col);
         std::cout << "datanames: " << datanames << std::endl;
         mtzin.import_hkl_data(fphi, datanames);
         mtzin.close_read();

         hklp.init(spgr, cell, reso, true);

         // make new coeffs
         int count = 0;
         clipper::HKL_data<clipper::data32::F_phi> fphi_new( hklp );
         for (HRI ih = fphi.first(); !ih.last(); ih.next()) {
             const clipper::data32::F_phi &f = fphi[ih.hkl()];
             fphi_new[ih].phi() = 0.0 ;
             if ( !f.missing() ) {
                fphi_new[ih].f() = f.f()*f.f();
                count++;
             }
         }
         if (debug)
            std::cout << "DEBUG:: processed " << count << " reflections" << std::endl;

         // xmap.spacegroup(), xmap.cell(), gs
         // clipper::Resolution reso(2.0);
         clipper::Grid_sampling gs(spgr, cell, reso);
         clipper::Xmap<float> xmap(spgr, cell, gs);
         xmap.fft_from(fphi_new);

         clipper::CCP4MAPfile mapout;
         mapout.open_write("born.map");
         mapout.export_xmap(xmap);
         mapout.close_write();

      }

      catch (const std::runtime_error &rte) {
         std::cout << "ERROR" << rte.what() << std::endl;
      }
   }

}
