
#include "coot-utils/coot-map-utils.hh"
#include "mini-texture.hh"

mini_texture_t::mini_texture_t(const clipper::Xmap<float> &xmap, int section_id) {

   std::pair<float, float> mv = coot::util::mean_and_variance(xmap);
   float mean = mv.first;
   float std_dev = std::sqrt(mv.second);

   // colour map info needs to be passed
   // for now we will do grey-scale
   float data_value_for_top    =  4.0;  // full white
   float data_value_for_bottom = -1.0;  // full black

   data_value_for_top    = mean + 2.5f * std_dev;
   data_value_for_bottom = mean - std_dev;
   float f_range = data_value_for_top - data_value_for_bottom;

   clipper::Grid_sampling gs = xmap.grid_sampling();
   std::cout << "mini_texture_t  constructor: "  << gs.format() << std::endl;
   image_data = new unsigned char[gs.nu() * gs.nv() * 4];

   clipper::Coord_grid cg_0(0,0,section_id);
   clipper::Coord_grid cg_1(gs.nu()-1, gs.nv()-1, section_id);
   clipper::Grid_map grid(cg_0, cg_1);
   clipper::Xmap_base::Map_reference_coord ix( xmap, grid.min()), iu, iv, iw;
   int nv = gs.nv();

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
            int idx = 4 * (c_v + nv * c_u);
            // std::cout << "idx " << idx << " idx_old " << idx_old << " c_u " << c_u << " nu " << nu << " c_v " << c_v << " nv " << nv << "\n";
            image_data[idx]   = static_cast<unsigned char>(255 * f_in_range);
            image_data[idx+1] = static_cast<unsigned char>(255 * f_in_range);
            image_data[idx+2] = static_cast<unsigned char>(255 * f_in_range);
            image_data[idx+3] = 255;
            // if (image_data[idx] > 200) { std::cout << idx << std::endl;}
            c_w++;
         }
         c_v++;
      }
      c_u++;
   }

   width  = gs.nv();
   height = gs.nu();
}

mini_texture_t::~mini_texture_t() {
   // delete [] image_data;
}

void
mini_texture_t::clear() {
   delete [] image_data;
}
