
#include <cmath>
#include <iostream>
#include <algorithm>

#include <clipper/core/map_interp.h>
#include "coot-coord-utils.hh"
#include "rebox-map.hh"

coot::util::reboxed_map_t
coot::util::rebox_map(const clipper::Xmap<float> &xmap_in,
                      mmdb::Manager *mol,
                      int SelectionHandle,
                      float border,
                      int n_pixels_per_edge) {

   reboxed_map_t result;

   if (!mol) {
      result.message = "Null molecule pointer";
      return result;
   }

   // 1. Get atom extents
   std::pair<clipper::Coord_orth, clipper::Coord_orth> ext = util::extents(mol, SelectionHandle);
   clipper::Coord_orth p_min = ext.first;
   clipper::Coord_orth p_max = ext.second;

   // check we got something sensible
   if (p_min.x() > p_max.x()) {
      result.message = "No atoms in selection";
      return result;
   }

   // 2. Add border
   p_min -= clipper::Coord_orth(border, border, border);
   p_max += clipper::Coord_orth(border, border, border);

   // 3. Determine the cube edge length (max of the three extents)
   double dx = p_max.x() - p_min.x();
   double dy = p_max.y() - p_min.y();
   double dz = p_max.z() - p_min.z();
   double edge_length = std::max({dx, dy, dz});

   // 4. Work out the grid spacing from the input map
   clipper::Cell orig_cell = xmap_in.cell();
   clipper::Grid_sampling orig_gs = xmap_in.grid_sampling();
   double gs_a = orig_cell.a() / static_cast<double>(orig_gs.nu());
   double gs_b = orig_cell.b() / static_cast<double>(orig_gs.nv());
   double gs_c = orig_cell.c() / static_cast<double>(orig_gs.nw());
   double grid_spacing = (gs_a + gs_b + gs_c) / 3.0;

   std::cout << "INFO:: rebox_map: input molecule extents "
             << ext.first.format() << " to " << ext.second.format() << std::endl;
   std::cout << "INFO:: rebox_map: padded extents dx=" << dx << " dy=" << dy << " dz=" << dz << std::endl;
   std::cout << "INFO:: rebox_map: input map grid spacing " << grid_spacing << " A" << std::endl;

   // 5. Determine grid size
   //    The edge length is always determined by the atom extents.
   //    If n_pixels_per_edge is specified, the grid spacing is adjusted to fit.
   //    Otherwise, the grid count is derived from the edge length and input map spacing.
   int n_grid;
   if (n_pixels_per_edge > 0) {
      n_grid = n_pixels_per_edge;
   } else {
      n_grid = static_cast<int>(std::ceil(edge_length / grid_spacing));
      if (n_grid < 1) n_grid = 1;
   }

   // 6. Centre the cube on the atom selection centre
   clipper::Coord_orth centre(0.5 * (ext.first.x() + ext.second.x()),
                              0.5 * (ext.first.y() + ext.second.y()),
                              0.5 * (ext.first.z() + ext.second.z()));

   double half_edge = 0.5 * edge_length;
   clipper::Coord_orth new_origin(centre.x() - half_edge,
                                  centre.y() - half_edge,
                                  centre.z() - half_edge);

   // The offset is the translation from the new map origin to the old coordinate frame.
   // To place the model into the new map frame, subtract this offset from the atom coords.
   clipper::Coord_orth offset = new_origin;

   // 7. Create new P1 cubic cell and Xmap
   clipper::Spacegroup sg_p1(clipper::Spacegroup::P1);
   clipper::Cell_descr cell_descr(edge_length, edge_length, edge_length, 90.0, 90.0, 90.0);
   clipper::Cell new_cell(cell_descr);
   clipper::Grid_sampling new_gs(n_grid, n_grid, n_grid);

   clipper::Xmap<float> new_xmap;
   new_xmap.init(sg_p1, new_cell, new_gs);

   // 8. Fill the new map by interpolating from the original
   clipper::Xmap_base::Map_reference_index ix;
   for (ix = new_xmap.first(); !ix.last(); ix.next()) {
      // grid coord in new map
      clipper::Coord_grid cg = ix.coord();
      // orthogonal position in new map frame
      clipper::Coord_frac cf = cg.coord_frac(new_gs);
      clipper::Coord_orth co = cf.coord_orth(new_cell);
      // translate back to original coordinate frame
      clipper::Coord_orth co_orig(co.x() + offset.x(),
                                  co.y() + offset.y(),
                                  co.z() + offset.z());
      // convert to fractional in original map
      clipper::Coord_frac cf_orig = co_orig.coord_frac(orig_cell);
      // interpolate (Interp_cubic wants Coord_map, not Coord_frac)
      clipper::Coord_map cm_orig = cf_orig.coord_map(orig_gs);
      float density = 0.0f;
      clipper::Interp_cubic::interp(xmap_in, cm_orig, density);
      new_xmap[ix] = density;
   }

   result.xmap = new_xmap;
   result.offset = offset;
   result.success = true;
   result.message = "Reboxed map: " + std::to_string(n_grid) + " grid points per edge, "
                  + "edge length " + std::to_string(edge_length) + " A";

   std::cout << "INFO:: rebox_map: " << result.message << std::endl;
   std::cout << "INFO:: rebox_map: offset " << offset.format() << std::endl;

   return result;
}
