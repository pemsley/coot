
#include <chrono>
#include "coot-utils/coot-map-utils.hh"
#include "mini-texture.hh"

mini_texture_t::mini_texture_t(const clipper::Xmap<float> &xmap, int section_index,
                               float data_value_for_bottom, float data_value_for_top) {

   auto tp_start = std::chrono::high_resolution_clock::now();

   data_value_for_top_of_range    = data_value_for_top;
   data_value_for_bottom_of_range = data_value_for_bottom;

   float f_range = data_value_for_top - data_value_for_bottom;

   clipper::Grid_sampling gs = xmap.grid_sampling();
   // std::cout << "mini_texture_t  constructor: "  << gs.format() << std::endl;
   int image_data_size = gs.nu() * gs.nv() * 4;
   image_data = new unsigned char[image_data_size];
   // auto tp_t_id2 = std::chrono::high_resolution_clock::now();
   // auto deltaid = std::chrono::duration_cast<std::chrono::milliseconds>(tp_t_id2 - tp_t_id1);


   clipper::Cell cell = xmap.cell();
   x_size = cell.a();
   y_size = cell.b();

   float z_frac = static_cast<float>(section_index) / static_cast<float>(gs.nw());
   z_position = z_frac * cell.c();

   clipper::Coord_grid cg_0(0,0,section_index);
   clipper::Coord_grid cg_1(gs.nu()-1, gs.nv()-1, section_index);
   clipper::Grid_map grid(cg_0, cg_1);
   clipper::Xmap_base::Map_reference_coord ix( xmap, grid.min()), iu, iv, iw;
   int nv = gs.nv();
   int nu = gs.nu();

   auto tp_pre_grid = std::chrono::high_resolution_clock::now();
   int c_u = 0;
   for ( iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() ) {
      int c_v = 0;
      for ( iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
         int c_w = 0;
         for ( iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
            const float &f = xmap[iw];
            float f_in_range = (f-data_value_for_bottom)/f_range;
            if (f_in_range < 0.0) f_in_range = 0.0;
            if (f_in_range > 1.0) f_in_range = 1.0;
            int y_coord = c_v;
            int x_coord = c_u;
            int idx = 4 * (y_coord + nv * x_coord);
            // idx = 4 * (x_coord + nv * y_coord);
            if (idx >= image_data_size) {
               std::cout << "out of index " << idx << " " << image_data_size << std::endl;
            } else {
               image_data[idx]   = static_cast<unsigned char>(255 * f_in_range);
               image_data[idx+1] = static_cast<unsigned char>(255 * f_in_range);
               image_data[idx+2] = static_cast<unsigned char>(255 * f_in_range);
               image_data[idx+3] = 255;
               // if (image_data[idx] > 200) { std::cout << idx << std::endl;}
            }
            c_w++;
         }
         c_v++;
      }
      c_u++;
   }
   width  = gs.nv();
   height = gs.nu();

   if (false) {
      auto tp_post_grid = std::chrono::high_resolution_clock::now();
      auto delta2  = std::chrono::duration_cast<std::chrono::milliseconds>(tp_post_grid - tp_pre_grid);
      auto deltaa  = std::chrono::duration_cast<std::chrono::milliseconds>(tp_post_grid - tp_start);
      // auto deltatr = std::chrono::duration_cast<std::chrono::milliseconds>(tp_t_1 - tp_start);
      std::cout << "constructor all:             " << deltaa.count()   << " ms" << std::endl;
      std::cout << "constructor grid time:       " << delta2.count()   << " ms" << std::endl;
      // std::cout << "constructor image-data time: " << deltaid.count()  << " ms" << std::endl;
      // std::cout << "constructor test_range:      " << deltatr.count()  << " ms" << std::endl;
   }

}

mini_texture_t::~mini_texture_t() {
   // delete [] image_data;
}

void
mini_texture_t::clear() {
   delete [] image_data;
}
