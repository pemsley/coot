
#include "coot-map-heavy.hh"

// the inner method is the same as is used for gaussian surface
clipper::Xmap<float> coot::util::make_gaussian_atom_map_for_mask(const clipper::Xmap<float> &map_ref,
                                                                 mmdb::Manager *mol,
                                                                 int atom_selection_handle,
                                                                 float sigma, float box_radius) {

   auto place_atom_in_grid = [] (const clipper::Coord_orth &pt,
                                 clipper::Xmap<float> &xmap,
                                 float sigma,
                                 float box_radius) {

      // std::cout << "place_atom_in_grid " << pt.format() << std::endl;
      clipper::Coord_frac centre_f = pt.coord_frac(xmap.cell());
      float box_radius_sqrd = box_radius * box_radius;

      clipper::Coord_frac box0(
                               centre_f.u() - box_radius/xmap.cell().descr().a(),
                               centre_f.v() - box_radius/xmap.cell().descr().b(),
                               centre_f.w() - box_radius/xmap.cell().descr().c() );
      clipper::Coord_frac box1(
                               centre_f.u() + box_radius/xmap.cell().descr().a(),
                               centre_f.v() + box_radius/xmap.cell().descr().b(),
                               centre_f.w() + box_radius/xmap.cell().descr().c() );

      clipper::Grid_map grid( box0.coord_grid(xmap.grid_sampling()),
                              box1.coord_grid(xmap.grid_sampling()));

      clipper::Xmap_base::Map_reference_coord ix( xmap, grid.min() ), iu, iv, iw;
      for (iu = ix; iu.coord().u() <= grid.max().u(); iu.next_u() )  {
         for (iv = iu; iv.coord().v() <= grid.max().v(); iv.next_v() ) {
            for (iw = iv; iw.coord().w() <= grid.max().w(); iw.next_w() ) {
               clipper::Coord_grid c_g = iw.coord();
               clipper::Coord_frac c_f = c_g.coord_frac(xmap.grid_sampling());
               clipper::Coord_orth c_o = c_f.coord_orth(xmap.cell());
               // float plength = clipper::Coord_orth::length(c_o, pt);
               float plength_sqrd = (c_o - pt).lengthsq();
               float z_sqrd = plength_sqrd/(sigma*sigma);
               float x_prime = - z_sqrd;
               if (plength_sqrd < box_radius_sqrd) {
                  // int idx = float_to_expo_index(x_prime);
                  // float v = expo[idx];
                  float atomic_number_scaling = 1.0;
                  float v = atomic_number_scaling * expf(x_prime);
                  xmap[iw] += v;
                  if (false)
                     std::cout << "adding " << v << " to " << c_g.format() << " from "
                               << x_prime << " using z " << std::sqrt(z_sqrd) << " x_prime " << x_prime
                               << std::endl;
               }
            }
         }
      }
   };

   clipper::Xmap<float> xmap(map_ref.spacegroup(), map_ref.cell(), map_ref.grid_sampling());

   mmdb::PPAtom atom_selection = 0;
   int n_atoms = 0;
   mol->GetSelIndex(atom_selection_handle, atom_selection, n_atoms);

   if (n_atoms > 0) {
      if (atom_selection) {
         for (int i=0; i<n_atoms; i++) {
            mmdb::Atom *at = atom_selection[i];
            if (! at->isTer()) {
               clipper::Coord_orth pos(at->x, at->y, at->z);
               place_atom_in_grid(pos, xmap, sigma, box_radius); // change xmap
            }
         }
      }
   }

   return xmap;
}

