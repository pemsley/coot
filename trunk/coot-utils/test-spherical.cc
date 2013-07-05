#include <iostream>
#include <string>
#include <vector>

#include <clipper/ccp4/ccp4_map_io.h>

#include <gsl/gsl_sf_legendre.h>

void calc_almn(const clipper::Xmap<float> &xmap, clipper::Coord_orth &pos) {

   float r = 10;
   std::cout << "almn here" << std::endl;
   
   clipper::Grid_sampling grid = xmap.grid_sampling();
   clipper::Cell cell = xmap.cell();
   clipper::Grid_range gr(cell, grid, r);
   clipper::Coord_grid g, g0, g1;
   clipper::Coord_orth c0, c1;
   typedef clipper::Xmap<float>::Map_reference_coord MRC;
   MRC i0, iu, iv, iw;
   g0 = g + gr.min();
   g1 = g + gr.max();
   i0 = MRC(xmap, g0);
   std::cout << "grid_range " << gr.format() << std::endl;
   for (iu = i0; iu.coord().u() <= g1.u(); iu.next_u()) {
      for ( iv = iu; iv.coord().v() <= g1.v(); iv.next_v() ) { 
	 for ( iw = iv; iw.coord().w() <= g1.w(); iw.next_w() ) {
	    c0 = iw.coord_orth() - pos;
	    float rho = xmap[iw];
	    std::cout << "c0: " << c0.format() << " rho: " << rho << std::endl;
	 }
      }
   }
} 

void test_spherical() {
   
   int l = 20;
   int m = 8;
   int n_bins = 100;
   
   for (unsigned int iv=0; iv<n_bins; iv++) { 
      double x = -1 + 2* double(iv)/double(n_bins);
      gsl_sf_result result;
      int success = gsl_sf_legendre_sphPlm_e(l, m, x, &result);
      // std::cout << "success: " << success << std::endl;
      std::cout << "result:  " << result.val << std::endl;
   }
}

int main(int argc, char **argv) {


   if (argc > 1) {
      std::string file_name(argv[1]);
      clipper::CCP4MAPfile file;
      try {
	 file.open_read(file_name);
	 clipper::Grid_sampling fgs = file.grid_sampling();
	 clipper::Xmap<float> xmap;
	 file.import_xmap(xmap);
	 std::cout << "map grid sampling " << xmap.grid_sampling().format() << std::endl;
	 clipper::Coord_orth pos(30,10,20);
	 calc_almn(xmap, pos);
      }
      catch (const clipper::Message_base &exc) {
	 std::cout << "WARNING:: failed to open " << file_name << std::endl;
      } 
   } 
   return 0;
} 
