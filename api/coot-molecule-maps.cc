/*
 * api/coot-molecule-maps.cc
 * 
 * Copyright 2020 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */


#include <thread>
#include <iostream>
#include <iomanip>

#include <clipper/ccp4/ccp4_map_io.h>

#include "utils/coot-utils.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "density-contour/occlusion.hh"
#include "density-contour/transfer-occlusions.hh"
#include "coot-molecule.hh"

// 20221126-PE set this for now. When it is restored, jiggle_fit_multi_thread_func_1 and jiggle_fit_multi_thread_func_2
// will need to be transfered.
// It is set (not configured) in ideal/simple-restraint.hh
#undef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

std::atomic<bool> coot::molecule_t::draw_vector_sets_lock(false);

void
coot::molecule_t::set_map_colour(coot::colour_holder ch) {
   map_colour = ch;
}

int
coot::molecule_t::write_map(const std::string &file_name) const {

   int status = 0;

   if (! xmap.is_null()) {
      clipper::CCP4MAPfile mapout;
      mapout.open_write(file_name);
      mapout.export_xmap(xmap);
      mapout.close_write();
      status = 1;
   }

   return status;

}

bool
coot::molecule_t::is_EM_map() const {

   bool ret_is_em = false;

   if (false)
      std::cout << "in coot::molecule::is_EM_map() A " << std::endl;

   if (has_xmap()) {
      if (is_em_map_cached_flag == 1) { // -1 means unset
         ret_is_em = true;
      }
   }
   return ret_is_em;
}


short int
coot::molecule_t::is_em_map_cached_state() {

   if (is_em_map_cached_flag == -1) {

      if (has_xmap()) { // FIXME - need to test for NXmap too.
         bool is_em = is_EM_map();
         is_em_map_cached_flag = is_em;
      }
   }
   return is_em_map_cached_flag;
}

void
coot::molecule_t::scale_map(float scale_factor) {

   // can be multi-threaded?

   if (has_xmap()) {
      clipper::Xmap_base::Map_reference_index ix;
      for (ix = xmap.first(); !ix.last(); ix.next() )
	 xmap[ix] *= scale_factor;
   }
}


void
coot::molecule_t::set_map_is_difference_map(bool state) {
   xmap_is_diff_map = state;
}

void
coot::molecule_t::associate_data_mtz_file_with_map(const std::string &data_mtz_file_name,
                                                   const std::string &f_col, const std::string &sigf_col,
                                                   const std::string &r_free_col) {

   // where should they be stored?
   refmac_mtz_filename = data_mtz_file_name;
   refmac_fobs_col     = f_col;
   refmac_sigfobs_col  = sigf_col;
   refmac_r_free_col = r_free_col;
   refmac_r_free_flag_sensible = true; // don't call this function unless this is true
                                       // ideally, this should be tested. sfcalc_genmaps_using_bulk_solvent()
                                       // crashes if this flag is not set true.
                                       // (that's not ideal, just to be clear) - it is however the current
                                       // state of things.

}



void
coot::molecule_t::update_map_triangles_using_thread_pool(float radius, coot::Cartesian centre, float contour_level, ctpl::thread_pool *thread_pool_p) {

   auto local_gensurf_and_add_vecs_threaded_workpackage = [] (unsigned int thread_index,
                                                              const clipper::Xmap<float> *xmap_p,
                                                              float contour_level, float dy_radius,
                                                              coot::Cartesian centre,
                                                              int isample_step,
                                                              int iream_start, int n_reams,
                                                              bool is_em_map,
                                                              std::vector<coot::density_contour_triangles_container_t> *draw_vector_sets_p,
                                                              std::atomic<unsigned int> &done_count_for_threads) {

      gensurf_and_add_vecs_threaded_workpackage(xmap_p, contour_level, dy_radius, centre, isample_step, iream_start, n_reams,
                                                is_em_map, draw_vector_sets_p);
      done_count_for_threads++;
   };

   if (!xmap.is_null()) {

      bool is_em_map = false;
      if (is_em_map_cached_state() == 1) {
         is_em_map = true;
      }

      int isample_step = 1;
      clear_draw_vecs();
      int n_reams = coot::get_max_number_of_threads() - 1;
      if (n_reams < 1) n_reams = 1;
      std::atomic<unsigned int> done_count_for_threads(0);
      for (int ii=0; ii<n_reams; ii++) {
         thread_pool_p->push(local_gensurf_and_add_vecs_threaded_workpackage,
                             &xmap, contour_level, radius, centre,
                             isample_step, ii, n_reams, is_em_map,
                             &draw_vector_sets, std::ref(done_count_for_threads));
      }
      while (done_count_for_threads < static_cast<unsigned int>(n_reams))
         std::this_thread::sleep_for(std::chrono::microseconds(3));

      if (xmap_is_diff_map) {
         clear_diff_map_draw_vecs();
         done_count_for_threads = 0;
         for (int ii=0; ii<n_reams; ii++) {
            thread_pool_p->push(local_gensurf_and_add_vecs_threaded_workpackage,
                                &xmap, -contour_level, radius, centre,
                                isample_step, ii, n_reams, is_em_map,
                                &draw_diff_map_vector_sets, std::ref(done_count_for_threads));
         }
         while (done_count_for_threads < static_cast<unsigned int>(n_reams))
            std::this_thread::sleep_for(std::chrono::microseconds(3));
      }
   }
}

void
coot::molecule_t::update_map_triangles(float radius, coot::Cartesian centre, float contour_level) {

      // std::cout   << "DEBUG:: update_map_triangles() at center: " << centre << std::endl;
   // std::cout   << "DEBUG:: update_map_triangles() g.zoom: " << g.zoom << std::endl;

   // duck out of doing map OpenGL map things if we are not in gui mode
   // (for figure making, from jupyter (say) in the future, this is probably not the right
   // thing to do.

   // if (! graphics_info_t::use_graphics_interface_flag) return;

   CIsoSurface<float> my_isosurface;
   coot::CartesianPairInfo v;
   int isample_step = 1;

   bool is_em_map = false;
   if (is_em_map_cached_state() == 1) {
      is_em_map = true;
   }


   // 20221008-PE this block were statics in the graphics_info_t (the molecules container)
   bool dynamic_map_resampling = false;
   bool dynamic_map_size_display = false;
   bool is_dynamically_transformed_map_flag = false;
   float zoom = 100;
   float dynamic_map_zoom_offset = 0.0;

   if (dynamic_map_resampling == 1)
      // isample_step = 1 + int (0.009*g.zoom);
      isample_step = 1 + int (0.009*(zoom + dynamic_map_zoom_offset));

   if (isample_step > 15)
      isample_step = 15;

   // for critical points of size display and resampling being different:
   //
   float dy_radius = radius;
   if (dynamic_map_size_display == 1) {
      if (isample_step <= 15 )
         dy_radius *= float(isample_step);
      else
         dy_radius *= 15.0;
   }

   //
   if (isample_step <= 0) {
      std::cout << "WARNING:: Bad zoom   ("<< zoom << "):  setting isample_step to 1" << std::endl;
      isample_step = 1;
   }
   if (dy_radius <= 0.0) {
      std::cout << "WARNING:: Bad radius (" << dy_radius << ") setting to 10" << std::endl;
      dy_radius = 10.0;
   }

   // dynamically transformed maps get their vectors from molecule B
   // (we are looking at molecule at atoms in molecule A) which then
   // have a the inverse of that transformation applied to them.
   //
   // But note that to get from the centre in the A molecule to the
   // corresponding centre in the B molecule we need to apply the
   // *inverse* of the transformation in the map_ghost_info.

   if (is_dynamically_transformed_map_flag) {
      clipper::Coord_orth c(centre.x(), centre.y(), centre.z());
      clipper::Coord_orth ct = c.transform(map_ghost_info.rtop.inverse());
      centre = coot::Cartesian(ct.x(), ct.y(), ct.z());
   }

   if (!xmap.is_null()) {

      clear_draw_vecs();
      std::vector<std::thread> threads;
      int n_reams = coot::get_max_number_of_threads() - 1;
      if (n_reams < 1) n_reams = 1;

      for (int ii=0; ii<n_reams; ii++) {
         threads.push_back(std::thread(gensurf_and_add_vecs_threaded_workpackage,
                                       &xmap, contour_level, dy_radius, centre,
                                       isample_step, ii, n_reams, is_em_map,
                                       &draw_vector_sets));
      }
      for (int ii=0; ii<n_reams; ii++)
         threads[ii].join();

      threads.clear();
      if (xmap_is_diff_map) {
         clear_diff_map_draw_vecs();
         for (int ii=0; ii<n_reams; ii++) {
            threads.push_back(std::thread(gensurf_and_add_vecs_threaded_workpackage,
                                          &xmap, -contour_level, dy_radius, centre,
                                          isample_step, ii, n_reams, is_em_map,
                                          &draw_diff_map_vector_sets));
         }
         for (int ii=0; ii<n_reams; ii++)
            threads[ii].join();

      }

      if (is_dynamically_transformed_map_flag) {
         for (unsigned int ii=0; ii<draw_vector_sets.size(); ii++) {
            // needs the type changing?       FIXME
            // dynamically_transform(draw_vector_sets[ii]);
         }
      }

      // post_process_map_triangles();

      if (false) {
         for (std::size_t i=0; i<draw_vector_sets.size(); i++) {
            coot::density_contour_triangles_container_t &tri_con = draw_vector_sets[i];
            std::vector<coot::augmented_position> positions(tri_con.points.size());
            unsigned int n = draw_vector_sets[i].points.size();
            for (unsigned int j=0; j<n; j++) {
               const clipper::Coord_orth &pos  = tri_con.points[j];
               const clipper::Coord_orth &norm = tri_con.normals[j];
               positions[j] = coot::augmented_position(pos, norm);
            }
            coot::set_occlusions(positions); // crash, related to range
            coot::transfer_occlusions(positions, &draw_vector_sets[i]);
         }
      }

   }
}


// this function is outside the coot_molecule class

void gensurf_and_add_vecs_threaded_workpackage(const clipper::Xmap<float> *xmap_p,
                                               float contour_level, float dy_radius,
                                               coot::Cartesian centre,
                                               int isample_step,
                                               int iream_start, int n_reams,
                                               bool is_em_map,
                                               std::vector<coot::density_contour_triangles_container_t> *draw_vector_sets_p) {

   try {
      CIsoSurface<float> my_isosurface;

      bool use_vertex_gradients_for_map_normals_flag = false; // give user control
      coot::density_contour_triangles_container_t tri_con =
        my_isosurface.GenerateTriangles_from_Xmap(std::cref(*xmap_p),
                                                  contour_level, dy_radius, centre, isample_step,
                                                  iream_start, n_reams, is_em_map,
                                                  use_vertex_gradients_for_map_normals_flag);

      // we are about to put the triangles into draw_vectors, so get the lock to
      // do that, so that the threads don't try to change draw_vectors at the same time.
      //
      bool unlocked = false;
      while (! coot::molecule_t::draw_vector_sets_lock.compare_exchange_weak(unlocked, true)) {
         std::this_thread::sleep_for(std::chrono::microseconds(10));
         unlocked = false;
      }

      // no longer dynamically change the size of draw_vector_sets
      // clear_draw_vecs will set the size to zero. If we find a element with size 0,
      // replace that one, rather than adding to draw_vector_sets
      //
      // draw_vector_sets_p->push_back(v);
      //
      bool done = false;
      for (unsigned int i=0; i<draw_vector_sets_p->size(); i++) {
         // std::cout << "gensurf_and_add_vecs_threaded_workpackage() checking i " << i << " data "
         // << draw_vector_sets_p->at(i).data << " size " << draw_vector_sets_p->at(i).size
         // << std::endl;
         if (draw_vector_sets_p->at(i).empty()) {
            // std::cout << "   replacing set at " << i << " data" << v.data << " size " << v.size
            // << std::endl;
            // perhaps I can std::move this? I don't need tri_con after this.

            // draw_vector_sets_p->at(i) = tri_con;
            std::move(tri_con.points.begin(), tri_con.points.end(), std::back_inserter(draw_vector_sets_p->at(i).points));
            std::move(tri_con.normals.begin(), tri_con.normals.end(), std::back_inserter(draw_vector_sets_p->at(i).normals));
            std::move(tri_con.point_indices.begin(), tri_con.point_indices.end(), std::back_inserter(draw_vector_sets_p->at(i).point_indices));
            done = true;
            break;
         }
      }
      if (! done) {
         // OK, let's push this one back then
         // std::cout << "gensurf_and_draw_vecs_threaded_workpackage() adding another draw vector set, "
         // << "current size " << draw_vector_sets_p->size() << " with " << v.data << " " << v.size
         // << std::endl;
         draw_vector_sets_p->push_back(tri_con);
      }

      coot::molecule_t::draw_vector_sets_lock = false; // unlock
   }
   catch (const std::out_of_range &oor) {
      std::cout << "ERROR:: contouring threaded workpackage " << oor.what() << std::endl;
   }
}

void
coot::molecule_t::clear_draw_vecs() {

   bool unlocked = false;
   while (!draw_vector_sets_lock.compare_exchange_weak(unlocked, true)) {
      std::this_thread::sleep_for(std::chrono::microseconds(10));
      unlocked = false;
   }
   for (std::size_t i=0; i<draw_vector_sets.size(); i++) {
      draw_vector_sets[i].clear();
   }
   draw_vector_sets_lock = false; // unlock

}

void
coot::molecule_t::clear_diff_map_draw_vecs() {

   bool unlocked = false;
   while (!draw_vector_sets_lock.compare_exchange_weak(unlocked, true)) {
      std::this_thread::sleep_for(std::chrono::microseconds(10));
      unlocked = false;
   }
   for (std::size_t i=0; i<draw_diff_map_vector_sets.size(); i++) {
      draw_diff_map_vector_sets[i].clear();
   }
   draw_vector_sets_lock = false; // unlock

}

coot::simple_mesh_t
coot::molecule_t::get_map_contours_mesh(clipper::Coord_orth position, float radius, float contour_level,
                                        bool use_thread_pool, ctpl::thread_pool *thread_pool_p) {

   // std::cout << "!!! ##### get_map_contours_mesh() for imol " << imol_no << std::endl;

   auto coord_orth_to_glm = [] (const clipper::Coord_orth &co) {
      return glm::vec3(co.x(), co.y(), co.z());
   };

   coot::simple_mesh_t m; // initially status is good (1).

   coot::Cartesian p(position.x(), position.y(), position.z());

   bool show_timings = true;
   auto tp_0 = std::chrono::high_resolution_clock::now();

   if (use_thread_pool)
      update_map_triangles_using_thread_pool(radius, p, contour_level, thread_pool_p);
   else
      update_map_triangles(radius, p, contour_level);
   auto tp_1 = std::chrono::high_resolution_clock::now();
   auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
   if (show_timings)
      std::cout << "Timings: map contouring " << d10 << " milliseconds" << std::endl;

   // now convert the contents of the draw-vector sets to a simple_mesh_t.

   // coot::colour_holder map_colour(0.3, 0.3, 0.8); // now is a class member
   if (xmap_is_diff_map)
      map_colour = coot::colour_holder(0.4, 0.8, 0.4);

   try {

      auto &vertices  = m.vertices;
      auto &triangles = m.triangles;

      std::vector<coot::density_contour_triangles_container_t>::const_iterator it;
      glm::vec4 col(map_colour.red, map_colour.green, map_colour.blue, 1.0f);
      for (it=draw_vector_sets.begin(); it!=draw_vector_sets.end(); ++it) {
         const coot::density_contour_triangles_container_t &tri_con(*it);
         unsigned int idx_base = vertices.size();
         for (unsigned int i=0; i<tri_con.points.size(); i++) {
            glm::vec3 pos    = coord_orth_to_glm(tri_con.points[i]);
            glm::vec3 normal = coord_orth_to_glm(- tri_con.normals[i]); // reverse normal.
            coot::api::vnc_vertex vert(pos, normal, col);
            vertices.push_back(vert);
         }
         for (unsigned int i=0; i<tri_con.point_indices.size(); i++) {
            g_triangle tri(tri_con.point_indices[i].pointID[0],
                           tri_con.point_indices[i].pointID[1],
                           tri_con.point_indices[i].pointID[2]);
            tri.rebase(idx_base);
            triangles.push_back(tri);

            // removed map triangle centres block here.

         }
      }

      if (xmap_is_diff_map) {
         glm::vec4 diff_map_col = glm::vec4(0.8, 0.4, 0.4, 1.0f);
         for (it=draw_diff_map_vector_sets.begin(); it!=draw_diff_map_vector_sets.end(); ++it) {
            const coot::density_contour_triangles_container_t &tri_con(*it);
            unsigned int idx_base = vertices.size();
            for (unsigned int i=0; i<tri_con.points.size(); i++) {
               glm::vec3 pos    = coord_orth_to_glm(tri_con.points[i]);
               glm::vec3 normal = coord_orth_to_glm(tri_con.normals[i]); // non-reverse normal!
               coot::api::vnc_vertex vert(pos, normal, diff_map_col);
               vertices.push_back(vert);
            }
            for (unsigned int i=0; i<tri_con.point_indices.size(); i++) {
               g_triangle tri(tri_con.point_indices[i].pointID[0],
                              tri_con.point_indices[i].pointID[1],
                              tri_con.point_indices[i].pointID[2]);
               tri.rebase(idx_base);
               triangles.push_back(tri);

               // removed map triangle centres block here.

            }
         }
      }
   }

   catch (const std::bad_alloc &ba) {
      std::cout << "WARNING:: in get_map_contours_mesh() bad alloc. " << std::endl;
      std::cout << "WARNING:: " << ba.what() << std::endl;
      m.clear();
      m.status = 0;
   }

   catch (...) {
      std::cout << "WARNING:: in get_map_contours_mesh() caught something else!" << std::endl;
      m.clear();
      m.status = 0;
   }


   const auto &vertices = m.vertices;

   if (false) {
      // I want to see if the vertices are greater than the cell (219x219x219)
      for (unsigned int i=0; i<vertices.size(); i++) {
         if (vertices[i].pos.x > 219.0 || vertices[i].pos.y > 219.0 || vertices[i].pos.z > 219.0) {
            std::cout << "ERROR:: vertex " << i << " " << vertices[i].pos.x << " "
                      << vertices[i].pos.y << " " << vertices[i].pos.z << std::endl;
         }
         // test for negative values also:
         if (vertices[i].pos.x < -1.0 || vertices[i].pos.y < -1.0 || vertices[i].pos.z < -1.0) {
            std::cout << "ERROR:: vertex " << i << " " << vertices[i].pos.x << " "
                      << vertices[i].pos.y << " " << vertices[i].pos.z << std::endl;
         }

         // print out everything:
         std::cout << "INFO::  vertex " << i << " " << vertices[i].pos.x << " "
                   << vertices[i].pos.y << " " << vertices[i].pos.z << std::endl;

      }
   }


   return m;
}

coot::simple_mesh_t
coot::molecule_t::get_map_contours_mesh_using_other_map_for_colours(const clipper::Coord_orth &position, float radius, float contour_level,
                                                                    const clipper::Xmap<float> &other_map) {

   auto coord_orth_to_glm = [] (const clipper::Coord_orth &co) {
      return glm::vec3(co.x(), co.y(), co.z());
   };

   auto clipper_to_cartesian = [] (const clipper::Coord_orth &c) {
      return Cartesian(c.x(), c.y(), c.z()); };

   auto cpos = clipper_to_cartesian(position);
   update_map_triangles(radius, cpos, contour_level);

   coot::simple_mesh_t m; // initially status is good (1).
   auto &vertices  = m.vertices;
   auto &triangles = m.triangles;
   std::vector<coot::density_contour_triangles_container_t>::const_iterator it;
   for (it=draw_vector_sets.begin(); it!=draw_vector_sets.end(); ++it) {
      const coot::density_contour_triangles_container_t &tri_con(*it);
      unsigned int idx_base = vertices.size();
      for (unsigned int i=0; i<tri_con.points.size(); i++) {
         glm::vec3 pos    = coord_orth_to_glm(tri_con.points[i]);
         glm::vec3 normal = coord_orth_to_glm(-tri_con.normals[i]); // reverse normals
         clipper::Coord_orth clipper_pos(pos.x, pos.y, pos.z);
         glm::vec4 col = position_to_colour_using_other_map(clipper_pos, other_map);
         api::vnc_vertex vert(pos, normal, col);
         vertices.push_back(vert);
      }
      for (unsigned int i=0; i<tri_con.point_indices.size(); i++) {
         g_triangle tri(tri_con.point_indices[i].pointID[0],
                        tri_con.point_indices[i].pointID[1],
                        tri_con.point_indices[i].pointID[2]);
         tri.rebase(idx_base);
         triangles.push_back(tri);
      }
   }
   return m;
}

coot::simple_mesh_t
coot::molecule_t::get_map_contours_mesh_using_other_map_for_colours(const clipper::Coord_orth &position,
								    float radius, float contour_level,
								    const user_defined_colour_table_t &udct,
								    const clipper::Xmap<float> &other_map) {

   auto coord_orth_to_glm = [] (const clipper::Coord_orth &co) {
      return glm::vec3(co.x(), co.y(), co.z());
   };

   auto clipper_to_cartesian = [] (const clipper::Coord_orth &c) {
      return Cartesian(c.x(), c.y(), c.z()); };

   auto frac_to_col = [] (float f, const user_defined_colour_table_t &udct) {
      glm::vec4 col(0.1f, 0.1f, 0.1f, 1.0f);
      for (unsigned int i=0; i<(udct.colour_table.size() -1); i++) {
	 if (udct.colour_table[i].colour_frac <= f) {
	    if (f <= udct.colour_table[i+1].colour_frac) {
	       user_defined_colour_table_t::colour_pair_t lower_point = udct.colour_table[i];
	       user_defined_colour_table_t::colour_pair_t upper_point = udct.colour_table[i+1];
	       float fraction = (f - lower_point.colour_frac) / (upper_point.colour_frac - lower_point.colour_frac);
	       float r = lower_point.colour.r + (upper_point.colour.r - lower_point.colour.r) * fraction;
	       float g = lower_point.colour.g + (upper_point.colour.g - lower_point.colour.g) * fraction;
	       float b = lower_point.colour.b + (upper_point.colour.b - lower_point.colour.b) * fraction;
	       col = glm::vec4(r,g,b,1.0f);
	    }
	 }
      }
      return col;
   };

   auto position_to_colour = [frac_to_col] (const clipper::Coord_orth &c,
					    const clipper::Xmap<float> &other_map,
					    const user_defined_colour_table_t &udct,
					    float min_value, float max_value) {
      float dv = util::density_at_point(other_map, c);
      float f = 0.0;
      if (dv < min_value) {
	 f = 0.0;
      } else {
	 if (dv > max_value) {
	    f = 1.0;
	 } else {
	    // in the range
	    float range = max_value - min_value;
	    float m = dv - min_value;
	    f = m/range;
	 }
      }
      return frac_to_col(f, udct);
   };

   auto cpos = clipper_to_cartesian(position);
   update_map_triangles(radius, cpos, contour_level);

   coot::simple_mesh_t m; // initially status is good (1).
   auto &vertices  = m.vertices;
   auto &triangles = m.triangles;
   std::vector<coot::density_contour_triangles_container_t>::const_iterator it;
   const float &min_value = other_map_for_colouring_min_value;
   const float &max_value = other_map_for_colouring_max_value;
   for (it=draw_vector_sets.begin(); it!=draw_vector_sets.end(); ++it) {
      const coot::density_contour_triangles_container_t &tri_con(*it);
      unsigned int idx_base = vertices.size();
      for (unsigned int i=0; i<tri_con.points.size(); i++) {
         glm::vec3 pos    = coord_orth_to_glm(tri_con.points[i]);
         glm::vec3 normal = coord_orth_to_glm(-tri_con.normals[i]); // reverse normals
         clipper::Coord_orth clipper_pos(pos.x, pos.y, pos.z);
         glm::vec4 col = position_to_colour(clipper_pos, other_map, udct, min_value, max_value);
         api::vnc_vertex vert(pos, normal, col);
         vertices.push_back(vert);
      }
      for (unsigned int i=0; i<tri_con.point_indices.size(); i++) {
         g_triangle tri(tri_con.point_indices[i].pointID[0],
                        tri_con.point_indices[i].pointID[1],
                        tri_con.point_indices[i].pointID[2]);
         tri.rebase(idx_base);
         triangles.push_back(tri);
      }
   }
   return m;
}


void
coot::molecule_t::set_other_map_for_colouring_min_max(float min_v, float max_v) {
   other_map_for_colouring_min_value = min_v;
   other_map_for_colouring_max_value = max_v;
}


glm::vec4
coot::molecule_t::position_to_colour_using_other_map(const clipper::Coord_orth &position,
                                                     const clipper::Xmap<float> &other_map_for_colouring) const {

   float dv = coot::util::density_at_point(other_map_for_colouring, position);
   float f = 0.0;
   const float &min_value = other_map_for_colouring_min_value;
   const float &max_value = other_map_for_colouring_max_value;
   if (dv < min_value) {
      f = 0.0;
   } else {
      if (dv > max_value) {
         f = 1.0;
      } else {
         // in the range
         float range = max_value - min_value;
         float m = dv - min_value;
         f = m/range;
      }
   }

   glm::vec4 col = fraction_to_colour(f);
   return col;
}

#include "colour-functions.hh"

glm::vec4
coot::molecule_t::fraction_to_colour(float fraction) const {

   // 20231101-PE innards reworked - now based on get_bond_colour_by_colour_wheel_position()

   float sat = radial_map_colour_saturation;
   // coot::colour_t cc(0.6+0.4*sat, 0.6-0.6*sat, 0.6-0.6*sat);
   std::vector<float> rgb(3);
   rgb[0] = 0.1f; rgb[1] =  0.1f; rgb[2] =  0.8f; // blue
   float max_colour = 30.0;
   // float rotation_size = float(fraction) * 36.0/max_colour - 0.2; pretty goo
   float rotation_size = float(fraction) * 36.0/max_colour - 0.28;
   rgb = rotate_rgb(rgb, rotation_size);

   glm::vec4 col(rgb[0], rgb[1], rgb[2], 1.0);

   return col;
}


#include "coot-utils/peak-search.hh"



// the molecule is passed so that the peaks are placed around the protein
std::vector<coot::molecule_t::interesting_place_t>
coot::molecule_t::difference_map_peaks(mmdb::Manager *mol, float n_rmsd) const {

   auto make_button_label = [] (unsigned int i, const std::pair<clipper::Coord_orth, float> &peak) {
      std::string s = "Index ";
      s += util::int_to_string(i);
      s += " Position: (";
      s += util::float_to_string_using_dec_pl(peak.first.x(), 2);
      s += std::string(", ");
      s += util::float_to_string_using_dec_pl(peak.first.y(), 2);
      s += std::string(", ");
      s += util::float_to_string_using_dec_pl(peak.first.z(), 2);
      s += std::string(") Height ");
      s += util::float_to_string_using_dec_pl(peak.second, 2);
      return s;
   };

   unsigned int n_limit = 100;

   std::vector<interesting_place_t> v;
   if (mol) {
      coot::peak_search ps(xmap);
      float rmsd = get_map_rmsd_approx();
      // float level = rmsd * n_rmsd; // not needed for skip_symmetry_check is false
      // this returns sorted peaks
      bool skip_symmetry_check = false;
      std::vector<std::pair<clipper::Coord_orth, float> > peaks = ps.get_peaks(xmap, mol, n_rmsd, true, true, skip_symmetry_check);
      for (unsigned int i=0; i<peaks.size(); i++) {
         const auto &peak = peaks[i];
         // difference_map_peaks_info_t dmp(peak.first, peak.second); // 20221105-PE as was, before generic type
         float badness = 7.0f * std::abs(peak.second/rmsd);
         std::string button_label = make_button_label(i, peak);
         interesting_place_t dmp("difference-map-peak", peak.first, button_label);
         dmp.residue_spec = get_residue_closest_to(mol, peak.first);
         dmp.set_feature_value(peak.second);
         dmp.set_badness_value(badness);
         v.push_back(dmp);
      }

      if (v.size() <= n_limit) {
         // sort them in numberical order (not absolute) - for the waterfall plot
         auto sorter_1 = [] (const interesting_place_t &i1, const interesting_place_t &i2) {
            return i1.feature_value < i2.feature_value;
         };
         std::sort(v.begin(), v.end(), sorter_1);
      } else {
         // first sort by absolute - then resize, then sort by value
         auto sorter_1 = [] (const interesting_place_t &i1, const interesting_place_t &i2) {
            return i1.feature_value < i2.feature_value;
         };
         auto sorter_2 = [] (const interesting_place_t &i1, const interesting_place_t &i2) {
            return fabsf(i2.feature_value) < fabsf(i1.feature_value);
         };
         std::sort(v.begin(), v.end(), sorter_2);
         v.resize(n_limit);
         std::sort(v.begin(), v.end(), sorter_1);
      }

   } else {
      std::cout << "ERROR:: " << __FUNCTION__ << "() null mol" << std::endl;
   }
   return v;
}


#include "coot-utils/xmap-stats.hh"

// map functions, return -1.0 on not-a-map
float
coot::molecule_t::get_map_mean() const {

   bool ignore_pseudo_zeros_for_map_stats = false; // set this to true for an EM map
   bool ipz = ignore_pseudo_zeros_for_map_stats;
   mean_and_variance<float> mv = map_density_distribution(xmap, 20, false, ipz);
   float m = mv.mean;
   return m;

}

// map functions, return -1.1 on not-a-map
float
coot::molecule_t::get_map_rmsd_approx() const {

   bool ignore_pseudo_zeros_for_map_stats = false; // set this to true for an EM map
   bool ipz = ignore_pseudo_zeros_for_map_stats;
   mean_and_variance<float> mv = map_density_distribution(xmap, 20, false, ipz);
   float rmsd = std::sqrt(mv.variance);
   return rmsd;

}

//! @return the suggested initial contour level. Return -1 on not-a-map
float
coot::molecule_t::get_suggested_initial_contour_level() const {

   float l = -1.0;

   if (is_valid_map_molecule()) {
      float rmsd = get_map_rmsd_approx();
      if (is_difference_map_p())
         l = 3.6 * rmsd;
      else
         l = 1.6 * rmsd;
      if (is_EM_map())
         l = 4.0 * rmsd;
   }
   return l;
}


bool
coot::molecule_t::is_difference_map_p() const {

   bool istat = false;
   if (is_valid_map_molecule())
      if (xmap_is_diff_map)
         istat = true;
   return istat;
}

// Return a pair.
//
// If first string of length 0 on error to construct dataname(s).
std::pair<std::string, std::string>
coot::molecule_t::make_import_datanames(const std::string &f_col_in,
                                        const std::string &phi_col_in,
                                        const std::string &weight_col_in,
                                        int use_weights) const {

   // If use_weights return 2 strings, else set something useful only for pair.first

   std::string f_col = f_col_in;
   std::string phi_col = phi_col_in;
   std::string weight_col = weight_col_in;

#ifdef WINDOWS_MINGW
   std::string::size_type islash_f   = coot::util::intelligent_debackslash(  f_col).find_last_of("/");
   std::string::size_type islash_phi = coot::util::intelligent_debackslash(phi_col).find_last_of("/");
#else
   std::string::size_type islash_f   =      f_col.find_last_of("/");
   std::string::size_type islash_phi =    phi_col.find_last_of("/");
#endif // MINGW

   short int label_error = 0;

   if (islash_f != std::string::npos) {
      // f_col is of form e.g. xxx/yyy/FWT
      if (f_col.length() > islash_f)
         f_col = f_col.substr(islash_f+1);
      else
         label_error = 1;
   }

   if (islash_phi != std::string::npos) {
      // phi_col is of form e.g. xxx/yyy/PHWT
      if (phi_col.length() > islash_phi)
         phi_col = phi_col.substr(islash_phi+1);
      else
         label_error = 1;
   }

   if (use_weights) {
      std::string::size_type islash_fom = weight_col.find_last_of("/");
      if (islash_fom != std::string::npos) {
         // weight_col is of form e.g. xxx/yyy/WT
         if (weight_col.length() > islash_fom)
            weight_col = weight_col.substr(islash_fom+1);
         else
            label_error = 1;
      }
   }


   std::pair<std::string, std::string> p("", "");

   if (!label_error) {
      std::string no_xtal_dataset_prefix= "/*/*/";
      if (use_weights) {
         p.first  = no_xtal_dataset_prefix + "[" +   f_col + " " +      f_col + "]";
         p.second = no_xtal_dataset_prefix + "[" + phi_col + " " + weight_col + "]";
      } else {
         p.first  = no_xtal_dataset_prefix + "[" +   f_col + " " + phi_col + "]";
      }
   }
   return p;
}



#include "clipper/ccp4/ccp4_mtz_io.h"
#include "clipper/cns/cns_hkl_io.h"
#include "clipper/cns/cns_map_io.h"
#include "clipper/core/hkl_compute.h"
#include "clipper/core/map_utils.h" // Map_stats
#include "clipper/core/resol_basisfn.h"
#include "clipper/core/resol_targetfn.h"
#include "clipper/mmdb/clipper_mmdb.h"
#include "clipper/clipper-phs.h"
#include "clipper/contrib/sfcalc_obs.h"
#include "clipper/contrib/sfscale.h"
#include "clipper/contrib/sfweight.h"

void
coot::molecule_t::fill_fobs_sigfobs() {

   bool show_timings = false;

   // set original_fobs_sigfobs_filled when done

   bool have_sensible_refmac_params = true; // 20221016-PE need to be set properly!

   if (have_sensible_refmac_params) {

      if (false)
         std::cout << "debug:: in fill_fobs_sigfobs() with original_fobs_sigfobs_filled " << original_fobs_sigfobs_filled
                   << " original_fobs_sigfobs_fill_tried_and_failed " << original_fobs_sigfobs_fill_tried_and_failed
                   << std::endl;

      // only try this once. If you try to import_hkl_data() when the original_fobs_sigfobs
      // already contains data, then crashiness.
      //

      if (! original_fobs_sigfobs_filled && ! original_fobs_sigfobs_fill_tried_and_failed) {

         auto tp_0 = std::chrono::high_resolution_clock::now();

         try {

            std::pair<std::string, std::string> p = make_import_datanames(Refmac_fobs_col(), Refmac_sigfobs_col(), "", 0);
            clipper::CCP4MTZfile *mtzin_p = new clipper::CCP4MTZfile; // original_fobs_sigfobs contains a pointer to
                                                                      // a cell in the crystals vector of a CCP4MTZfile.
                                                                      // The CCP4MTZfile goes out of score and takes
                                                                      // the crystal vector with it.
                                                                      // crystals is a vector of crystalinfo, which
                                                                      // is a structure that contains a MTZcrystal
                                                                      // which inherits from a Cell
                                                                      // Or something like that. ccp4_mtz_types.h,
                                                                      // ccp4_mtz_io.h and ccp4_mtz_io.cpp ::import_hkldata().
                                                                      // Anyway, something seems to go out of scope when
                                                                      // the molecule vector is resized. So
                                                                      // regenerate original_fobs_sigfobs from
                                                                      // the mtz file every time we need them.
                                                                      // This leak memory.  Meh... but better than
                                                                      // crashing. Likewise mtzin_p for R-free.
                                                                      // (20 each ms for RNAse dataset). 20210816-PE

                                                                      // Later note: now that original_fobs_sigfobs is a pointer
                                                                      // I probably don't need to mtzin object to be pointers.

            original_fobs_sigfobs_p = new clipper::HKL_data< clipper::datatypes::F_sigF<float> >;
            original_r_free_flags_p = new clipper::HKL_data< clipper::data32::Flag>;

            mtzin_p->open_read(Refmac_mtz_filename());
            mtzin_p->import_hkl_data(*original_fobs_sigfobs_p, p.first);
            mtzin_p->close_read();
            if (false)
               std::cout << "INFO:: fill_fobs_sigfobs(): reading " << Refmac_mtz_filename() << " provided "
                         << original_fobs_sigfobs_p->num_obs() << " data using data name: "
                         << p.first << std::endl;
            if (original_fobs_sigfobs_p->num_obs() > 10)
               original_fobs_sigfobs_filled = 1;
            else
               original_fobs_sigfobs_fill_tried_and_failed = true;

            // flags

            if (refmac_r_free_flag_sensible) {
               std::string dataname = "/*/*/[" + refmac_r_free_col + "]";
               // if refmac_r_free_col already has /x/y/Rfree - use that instead
               if (refmac_r_free_col.length() > 0) {
                  if (refmac_r_free_col[0] == '/') {
                     dataname = refmac_r_free_col;
                     dataname = "/*/*/[" + coot::util::file_name_non_directory(refmac_r_free_col) + "]";
                  }
               }
               if (false)
                  std::cout << "INFO:: About to read " << Refmac_mtz_filename() << " with dataname " << dataname << std::endl;
               clipper::CCP4MTZfile *mtzin_rfree_p = new clipper::CCP4MTZfile;
               mtzin_rfree_p->open_read(Refmac_mtz_filename());
               mtzin_rfree_p->import_hkl_data(*original_r_free_flags_p, dataname);
               mtzin_rfree_p->close_read();

               if (false)
                  std::cout << "INFO:: reading " << Refmac_mtz_filename() << " using dataname: " << dataname << " provided "
                            << original_r_free_flags_p->num_obs() << " R-free flags\n";
            } else {
               std::cout << "INFO:: no sensible R-free flag column label\n";
            }
         }
         catch (const clipper::Message_fatal &m) {
            std::cout << "ERROR:: bad columns " << m.text() << std::endl;
            have_sensible_refmac_params = false;
            original_fobs_sigfobs_filled = false;
            original_fobs_sigfobs_fill_tried_and_failed = true;
         }

         auto tp_1 = std::chrono::high_resolution_clock::now();
         auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
         if (show_timings)
            std::cout << "Timings: read mtz file and store data " << d10 << " milliseconds" << std::endl;
      }
   } else {
      std::cout << "DEBUG:: fill_fobs_sigfobs() no Fobs parameters\n";
   }
}




//! My ligands don't jiggle-jiggle...
//! Hey, what do you know, they actually do.
float
coot::molecule_t::fit_to_map_by_random_jiggle(const residue_spec_t &res_spec, const clipper::Xmap<float> &xmap, float map_rmsd,
                                              int n_trials, float jiggle_scale_factor) {

   float v = -1001.0;
   mmdb::Residue *residue_p = get_residue(res_spec);
   if (residue_p) {
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      bool use_biased_density_scoring = true;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      std::vector<mmdb::Chain *> chains; // empty - apply RTop to atoms of selection
      v = fit_to_map_by_random_jiggle(residue_atoms, n_residue_atoms, xmap, map_rmsd,
                                      n_trials, jiggle_scale_factor, use_biased_density_scoring, chains);
   } else {
      std::cout << "WARNING:: residue " << res_spec << " not found" << std::endl;
   }
   return v;
}

float
coot::molecule_t:: fit_to_map_by_random_jiggle_using_atom_selection(const std::string &cid, const clipper::Xmap<float> &xmap, float map_rmsd,
                                                                 int n_trials, float translation_scale_factor) {

   float v = -1001.0;
   if (is_valid_model_molecule()) {
      mmdb::PPAtom atoms = 0;
      int n_atoms;
      int selHnd = atom_sel.mol->NewSelection(); // d
      atom_sel.mol->Select(selHnd, mmdb::STYPE_ATOM, cid.c_str(), mmdb::SKEY_NEW);
      atom_sel.mol->GetSelIndex(selHnd, atoms, n_atoms);
      bool use_biased_density_scoring = true;
      std::vector<mmdb::Chain *> chains; // empty - apply RTop to atoms of selection
      v = fit_to_map_by_random_jiggle(atoms, n_atoms, xmap, map_rmsd, n_trials, translation_scale_factor, use_biased_density_scoring, chains);
      atom_sel.mol->DeleteSelection(selHnd);
   }
   return v;
}




#include "coot-utils/coot-map-utils.hh"
#include "coot-utils/coot-map-heavy.hh"

// called by above and split_water.
//
// chain_for_moving is default arg, with value 0.
//
// if chain_for_moving is not empty, apply the transformation
// the the atoms of chain_for_moving rather than to the atom of atom_selection
//
float
coot::molecule_t::fit_to_map_by_random_jiggle(mmdb::PPAtom atom_selection,
                                                   int n_atoms,
                                                   const clipper::Xmap<float> &xmap,
                                                   float map_sigma,
                                                   int n_trials,
                                                   float jiggle_scale_factor,
                                                   bool use_biased_density_scoring,
                                                   std::vector<mmdb::Chain *> chains_for_moving) {

   // 20221119-PE in api mode, compiled with CMake files, HAVE_BOOST_BASED_THREAD_POOL_LIBRARY is not set

   float v = 0;

   if (! atom_sel.mol) return v;

   std::vector<std::pair<std::string, int> > atom_numbers = coot::util::atomic_number_atom_list();
   if (n_trials <= 0)
      n_trials = 1;

   // set atoms so that we can get an initial score.
   std::vector<mmdb::Atom *> initial_atoms(n_atoms);
   std::vector<mmdb::Atom> direct_initial_atoms(n_atoms);
   for (int iat=0; iat<n_atoms; iat++)
      initial_atoms[iat] = atom_selection[iat];

   // We have to make a copy because direct_initial_atoms goes out of
   // scope and destroys the mmdb::Atoms (we don't want to take the
   // contents of the atom_selection out when we do that).
   //
   for (int iat=0; iat<n_atoms; iat++)
      direct_initial_atoms[iat].Copy(atom_selection[iat]);

   coot::minimol::molecule direct_mol(atom_selection, n_atoms, direct_initial_atoms);

   float (*density_scoring_function)(const coot::minimol::molecule &mol,
                                     const std::vector<std::pair<std::string, int> > &atom_number_list,
                                     const clipper::Xmap<float> &map) = coot::util::z_weighted_density_score_linear_interp;

   if (true)
      density_scoring_function = coot::util::z_weighted_density_score_nearest;

   // if (use_biased_density_scoring)
   //   density_scoring_function = coot::util::biased_z_weighted_density_score;


   // what residues are near to but not in atom_selection?
   //
   std::vector<mmdb::Residue *> neighbs;  // fill this
   //
   // we want to use residues_near_position(), so we want a list of residue that will be each of
   // the target residues for residues_near_residue().
   //
   std::vector<mmdb::Residue *> central_residues;
   for (int iat=0; iat<n_atoms; iat++) {
      mmdb::Atom *at = atom_selection[iat];
      mmdb::Residue *r = at->GetResidue();
      if (std::find(central_residues.begin(), central_residues.end(), r) == central_residues.end()) {
         central_residues.push_back(r);
      }
   }

   if (false)
      for (unsigned int ii=0; ii<central_residues.size(); ii++)
         std::cout << "            central residue: " << coot::residue_spec_t(central_residues[ii]) << std::endl;

   float radius = 4.0;
   for (unsigned int ires=0; ires<central_residues.size(); ires++) {
      mmdb::Residue *res_ref = central_residues[ires];
      std::pair<bool, clipper::Coord_orth> pt = util::get_residue_centre(res_ref);
      if (pt.first) {
         std::vector<mmdb::Residue *> r_residues =
         coot::residues_near_position(pt.second, atom_sel.mol, radius);
         for (unsigned int ii=0; ii<r_residues.size(); ii++) {
            if (std::find(neighbs.begin(), neighbs.end(), r_residues[ii]) == neighbs.end())
            if (std::find(central_residues.begin(), central_residues.end(), r_residues[ii]) == central_residues.end())
            neighbs.push_back(r_residues[ii]);
         }
      }
   }

   clipper::Xmap<float> xmap_masked = coot::util::mask_map(xmap, neighbs);

   // best score is the inital score (without the atoms being jiggled) (could be a good score!)
   //
   float initial_score = density_scoring_function(direct_mol, atom_numbers, xmap_masked);
   // float initial_score = coot::util::z_weighted_density_score(direct_mol, atom_numbers, xmap);
   // initial_score = coot::util::biased_z_weighted_density_score(direct_mol, atom_numbers, xmap);

   v = initial_score;
   float best_score = initial_score;

   std::cout << "---------------- initial_score " << initial_score << " ---------------" << std::endl;
   bool bested = false;
   coot::minimol::molecule best_molecule;
   clipper::RTop_orth best_rtop;

   // first, find the centre point.  We do that because otherwise we
   // do it lots of times in jiggle_atoms.  Inefficient.
   std::vector<double> p(3, 0.0);
   for (int iat=0; iat<n_atoms; iat++) {
      p[0] += atom_selection[iat]->x;
      p[1] += atom_selection[iat]->y;
      p[2] += atom_selection[iat]->z;
   }
   double fact = 1.0;
   if (n_atoms)
      fact = 1.0/float(n_atoms);
   clipper::Coord_orth centre_pt(p[0]*fact, p[1]*fact, p[2]*fact);

   std::vector<std::pair<clipper::RTop_orth, float> > trial_results(n_trials);
   bool do_multi_thread = false;

// 20240112-PE This needs an update -  see top of file
#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
   int n_threads = coot::get_max_number_of_threads();
   if (n_threads > 0)
      do_multi_thread = true;
#endif

   if (do_multi_thread) {

#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
      ctpl::thread_pool thread_pool(n_threads);
      try {
         unsigned int n_threads = coot::get_max_number_of_threads();

         for (int itrial=0; itrial<n_trials; itrial++) {

            auto tp_1 = std::chrono::high_resolution_clock::now();

            // throw into the mix a model that has been only small rotated/translate
            // or maybe nothing at all.

            thread_pool.push(jiggle_fit_multi_thread_func_1, itrial, n_trials, atom_selection, n_atoms,
                             initial_atoms, centre_pt, jiggle_scale_factor, atom_numbers,
                             &xmap_masked, // pointer arg, no std::ref()
                             density_scoring_function, &trial_results[itrial]);

            auto tp_2 = std::chrono::high_resolution_clock::now();
            auto d21 = std::chrono::duration_cast<std::chrono::microseconds>(tp_2 - tp_1).count();
            // not to self: it takes 40ms to copy a const xmap reference to the function.
            // question for self: was it actually a reference though? I suspect not, because
            // std::ref() was not in the code until I (just) added it.
            //

            // this is useful for debugging, but makes a mess
            if (false)
               std::cout << "pushing trial thread into pool: " << itrial << " " << d21
                         << " microseconds" << std::endl;
         }

         // wait for thread pool to finish jobs.
         bool wait_continue = true;
         while (wait_continue) {
            std::this_thread::sleep_for(std::chrono::milliseconds(200));
            if (thread_pool.n_idle() == thread_pool.size())
               wait_continue = false;
         }
      }

      catch (const std::bad_alloc &ba) {
         std::cout << "ERROR:: ------------------------ out of memory! ----------------- " << std::endl;
         std::cout << "ERROR:: " << ba.what() << std::endl;
      }
#endif
   } else {

      for (int itrial=0; itrial<n_trials; itrial++) {
         float annealing_factor = 1 - float(itrial)/float(n_trials);
         std::pair<clipper::RTop_orth, std::vector<mmdb::Atom> > jiggled_atoms =
         coot::util::jiggle_atoms(initial_atoms, centre_pt, jiggle_scale_factor);
         coot::minimol::molecule jiggled_mol(atom_selection, n_atoms, jiggled_atoms.second);
         if (false) { // debug solutions
            std::string jfn = "jiggled-" + std::to_string(itrial) + ".pdb";
            jiggled_mol.write_file(jfn, 20.0);
         }
         float this_score = density_scoring_function(jiggled_mol, atom_numbers, xmap_masked);
         std::pair<clipper::RTop_orth, float> p(jiggled_atoms.first, this_score);
         trial_results[itrial] = p;
      }
   }

   int n_for_rigid = int(float(n_trials) * 0.1);
   if (n_for_rigid > 10) n_for_rigid = 10;
   if (n_for_rigid == 0)  n_for_rigid = 1;

   if (false) {
      unsigned int n_top = 20;
      if (trial_results.size() < 20)
         n_top = trial_results.size();
      for (unsigned int i_trial=0; i_trial<n_top; i_trial++)
         std::cout << " debug pre-sort trial scores: " << i_trial << " " << trial_results[i_trial].second << std::endl;
   }

   auto trial_results_comparer = [] (const std::pair<clipper::RTop_orth, float> &a,
                                     const std::pair<clipper::RTop_orth, float> &b) {
      return (b.second < a.second);
   };

   std::sort(trial_results.begin(),
             trial_results.end(),
             trial_results_comparer);

   // sorted results (debugging)
   if (false) {
      unsigned int n_top = 20;
      if (trial_results.size() < 20)
         n_top = trial_results.size();
      for (unsigned int i_trial=0; i_trial<n_top; i_trial++)
         std::cout << " debug sorted trials: " << i_trial << " " << trial_results[i_trial].second << std::endl;
   }

   // Here grid-search each of top n_for_rigid solution, replacing
   // each by best of grid-search results.  {5,10} degrees x 3 angles?

   clipper::RTop_orth rtop_orth_identity;
   std::pair<clipper::RTop_orth, float> start_pair(rtop_orth_identity, 0);
   // these get updated in the upcoming loop
   std::vector<std::pair<clipper::RTop_orth, float> > post_fit_trial_results = trial_results;

   float best_score_so_far = -999999;

#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

   // fit and score best random jiggled results

   try {
      ctpl::thread_pool thread_pool(n_threads);
      for (int i_trial=0; i_trial<n_for_rigid; i_trial++) {
         // does the fitting
         thread_pool.push(jiggle_fit_multi_thread_func_2, direct_mol,
                          std::cref(xmap_masked), map_sigma,
                          centre_pt, atom_numbers,
                          trial_results[i_trial].second,
                          density_scoring_function,
                          &post_fit_trial_results[i_trial]);
      }

      // wait
      std::cout << "waiting for rigid-body fits..." << std::endl;
      bool wait_continue = true;
      while (wait_continue) {
         std::this_thread::sleep_for(std::chrono::milliseconds(200));
         if (thread_pool.n_idle() == thread_pool.size())
            wait_continue = false;
      }
   }
   catch (const std::bad_alloc &ba) {
      std::cout << "ERROR:: ------------------------ out of memory! ----------------- " << std::endl;
      std::cout << "ERROR:: " << ba.what() << std::endl;
   }

#else

   // non-threaded, fit and score top jiggled results

   for (int i_trial=0; i_trial<n_for_rigid; i_trial++) {
      coot::minimol::molecule  trial_mol = direct_mol;

      trial_mol.transform(trial_results[i_trial].first, centre_pt);
      coot::minimol::molecule fitted_mol = rigid_body_fit(trial_mol, xmap_masked, map_sigma);
      float this_score = density_scoring_function(fitted_mol, atom_numbers, xmap_masked);
      std::cout << "INFO:: Jiggle-fit: optimizing trial "
                << std::setw(3) << i_trial << ": prelim-score was "
                << std::setw(7) << trial_results[i_trial].second << " post-fit "
                << std::setw(5) << this_score;
      if (this_score > best_score_so_far) {
         best_score_so_far = this_score;
         if (this_score > initial_score) {
            std::cout << " ***";
         }
      }
      std::cout << std::endl;
      post_fit_trial_results[i_trial].second = this_score;
   }
#endif // HAVE_CXX_THREAD

   std::sort(post_fit_trial_results.begin(),
             post_fit_trial_results.end(),
             trial_results_comparer);

   if (true) {
      unsigned int n_top = 10;
      if (trial_results.size() < 10)
         n_top = trial_results.size();
      for (unsigned int i_trial=0; i_trial<n_top; i_trial++)
         std::cout << " debug:: raw        sorted trials: " << i_trial << " " << trial_results[i_trial].second << std::endl;
      for (int i_trial=0; i_trial<n_for_rigid; i_trial++)
         std::cout << " debug:: post-rigid sorted trials: " << i_trial << " " << post_fit_trial_results[i_trial].second << std::endl;
   }

   if (post_fit_trial_results[0].second > initial_score) {
      bested = true;
      best_rtop = post_fit_trial_results[0].first; // the rtop from before the rigid-body fitting
      coot::minimol::molecule  post_fit_mol = direct_mol;
      post_fit_mol.transform(post_fit_trial_results[0].first, centre_pt);
      coot::minimol::molecule fitted_mol = rigid_body_fit(post_fit_mol, xmap_masked, map_sigma);
      best_molecule = fitted_mol;

      float this_score = density_scoring_function(fitted_mol, atom_numbers, xmap_masked);
      std::cout << "INFO:: chose new molecule with score " << this_score << std::endl;
      best_score = this_score;
   }

   std::cout << "debug:: ............ here with bested: " << bested << std::endl;

   //
   if (bested) {
      make_backup("fit to map by random jiggle");
      std::cout << "INFO:: Improved fit from " << initial_score << " to " << best_score << std::endl;

      v = best_score;
      if (! best_molecule.is_empty()) {
         mmdb::Manager *mol = best_molecule.pcmmdbmanager();
         if (mol) {

            if (!chains_for_moving.empty()) {

               // move the atoms of chain for moving, not the atoms of the atom selection
               //
               // now fitted_mol contains the atoms of the atom selection fitted to density
               // We need to find the transformation from the current/original coordintes
               // to that fitted mol coordinates and then apply them to all the atom
               // in the chain

               std::cout << "DEBUG: we have " << chains_for_moving.size() << " chains for moving"
               << std::endl;

               for (unsigned int ich=0; ich<chains_for_moving.size(); ich++) {
                  mmdb::Chain *chain_for_moving = chains_for_moving[ich];
                  std::string chain_id = chain_for_moving->GetChainID();
                  std::cout << "DEBUG:: chain_for_moving " << chain_for_moving << " " << chain_id
                  << std::endl;
                  std::pair<int, int> mmr = coot::util::min_and_max_residues(chain_for_moving);
                  if (mmr.second >= mmr.first) {
                     std::vector<coot::lsq_range_match_info_t> matches;
                     coot::lsq_range_match_info_t match(mmr.first,
                                                        mmr.second, chain_id,
                                                        mmr.first,
                                                        mmr.second, chain_id,
                                                        lsq_t::MAIN);
                     matches.push_back(match);
                     mmdb::Manager *mol_1 = mol;
                     mmdb::Manager *mol_2 = atom_sel.mol;
                     std::pair<short int, clipper::RTop_orth> lsq_mat =
                        coot::util::get_lsq_matrix(mol_1, mol_2, matches, 1, false);
                     if (lsq_mat.first) {
                        const clipper::RTop_orth &rtop_of_fitted_mol = lsq_mat.second;
                        coot::util::transform_chain(chain_for_moving, rtop_of_fitted_mol);
                        std::cout << "DEBUG:: transforming chain " << chain_id << ":\n";
                        std::cout << rtop_of_fitted_mol.format() << "\n";

                     } else {
                        std::cout << "WARNING:: failed to make a rtop matrix" << std::endl;
                     }
                  }
               }

            } else {

               if (false) { // debugging fit.
                  int imod = 1;
                  mmdb::Model *model_p = mol->GetModel(imod);
                  if (model_p) {
                     int n_chains = model_p->GetNumberOfChains();
                     for (int ichain=0; ichain<n_chains; ichain++) {
                        mmdb::Chain *chain_p = model_p->GetChain(ichain);
                        int n_res = chain_p->GetNumberOfResidues();
                        for (int ires=0; ires<n_res; ires++) {
                           mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                           if (residue_p) {
                              int n_atoms = residue_p->GetNumberOfAtoms();
                              for (int iat=0; iat<n_atoms; iat++) {
                                 mmdb::Atom *at = residue_p->GetAtom(iat);
                                 if (! at->isTer()) {
                                    clipper::Coord_orth pt(at->x, at->y, at->z);
                                    float d = coot::util::density_at_point(xmap, pt);
                                    std::cout << coot::atom_spec_t(at) << " " << d << std::endl;
                                 }
                              }
                           }
                        }
                     }
                  }
               }

               // there were no chains for moving, so we apply the rtop to everything
               // with replace_coords().

               atom_selection_container_t asc_ligand = make_asc(mol);
               replace_coords(asc_ligand, false, true);
               // asc_ligand.mol->WritePDBASCII("asc_ligand.pdb");
            }

            // have_unsaved_changes_flag = 1;
            // make_bonds_type_checked();
         } else {
            std::cout << "ERROR:: fit_to_map_by_random_jiggle(): mol is null! " << std::endl;
         }
         delete mol;
      } else {
         std::cout << "ERROR:: fit_to_map_by_random_jiggle(): best_molecule is empty!" << std::endl;
      }
      // save_info.new_modification("fit_to_map_by_random_jiggle");
   } else {
      std::cout << " nothing better found " << std::endl;
   }
   return v;
}

#include "ligand/ligand.hh"

// return a fitted molecule
coot::minimol::molecule
coot::molecule_t::rigid_body_fit(const coot::minimol::molecule &mol_in,
                                 const clipper::Xmap<float> &xmap,
                                 float map_sigma) const {

   coot::ligand lig;
   lig.import_map_from(xmap, map_sigma);
   lig.install_ligand(mol_in);
   lig.find_centre_by_ligand(0); // don't test ligand size
   lig.set_map_atom_mask_radius(0.5);
   lig.set_dont_write_solutions();
   lig.set_dont_test_rotations();
   lig.set_acceptable_fit_fraction(0.1);
   lig.fit_ligands_to_clusters(1);
   unsigned int iclust = 0;
   unsigned int isol   = 0;
   coot::minimol::molecule moved_mol = lig.get_solution(isol, iclust);
   return moved_mol;
}


coot::util::map_molecule_centre_info_t
coot::molecule_t::get_map_molecule_centre() const {

   util::map_molecule_centre_info_t mc = util::map_molecule_centre(xmap);
   return mc;
}


//! @return the map histogram
coot::molecule_t::histogram_info_t
coot::molecule_t::get_map_histogram(unsigned int n_bins_in, float zoom_factor) const {

   float n_bins_fl = static_cast<float>(n_bins_in) * zoom_factor;
   int n_bins = static_cast<int>(n_bins_fl);
   mean_and_variance<float> mv = map_density_distribution(xmap, n_bins, false, false);
   float mean = mv.mean;
   float prev_range = mean - mv.min_density;
   float new_range = prev_range/zoom_factor;
   float new_min_density = mean - new_range;

   // how many bins are there between mv.min_density and new_min_density?
   unsigned int count = 0;
   float level = mv.min_density;
   while (level < new_min_density) {
      level += mv.bin_width;
      count++;
      // sanity
      if (count > 9999) break; // 20231023-PE needed a bigger limit for large zoom
   }

   // std::cout << "n_bins_in " << n_bins_in << " zoom_factor " << zoom_factor << " n_bins " << n_bins << std::endl;
   // std::cout << "Now create new_bins by removing the first " << count << " entries from mv.bins "
   // << "and limiting number of bins" << std::endl;

   std::vector<int> new_bins(n_bins_in, 0);
   for (unsigned int ibin=0; ibin<mv.bins.size(); ibin++) {
      int new_index = ibin - count;
      if (new_index >= 0) {
         if (new_index < static_cast<int>(n_bins_in))
            new_bins[new_index] = mv.bins[ibin];
      }
   }

   coot::molecule_t::histogram_info_t hi(new_min_density, mv.bin_width, new_bins);
   hi.mean = mean;
   hi.variance = mv.variance;
   return hi;
}

// just look at the vertices of the map - not the whole thing
coot::molecule_t::histogram_info_t
coot::molecule_t::get_map_vertices_histogram(const clipper::Xmap<float> &other_xmap,
					     const clipper::Coord_orth &pt,
					     float radius, float contour_level,
					     bool use_thread_pool, ctpl::thread_pool *thread_pool_p,
					     unsigned int n_bins) {


   coot::molecule_t::histogram_info_t hi;
   coot::simple_mesh_t mesh = get_map_contours_mesh(pt, radius, contour_level,
						    use_thread_pool, thread_pool_p);
   unsigned int n_points = 0;
   float max_d = -10e6f;
   float min_d =  10e6f;
   std::vector<float> d_values;
   d_values.reserve(40000);
   std::vector<coot::density_contour_triangles_container_t>::const_iterator it;
   for (it=draw_vector_sets.begin(); it!=draw_vector_sets.end(); ++it) {
      const coot::density_contour_triangles_container_t &tri_con(*it);
      for (unsigned int i=0; i<tri_con.points.size(); i++) {
	 const clipper::Coord_orth &co = tri_con.points[i];
	 float d = coot::util::density_at_point(other_xmap, co);
	 if (d > max_d) max_d = d;
	 if (d < min_d) min_d = d;
	 n_points++;
	 d_values.push_back(d);
      }
   }
   if (n_points > 0) {
      if (n_bins > 0) {
	 float range = max_d - min_d;
	 if (range > 0.0) {
	    float bin_width = range / static_cast<float>(n_bins);
	    float inv_range = 1.0f/range;
	    std::map<int, int> counts;
	    double sum = 0.0;
	    double sum_sq = 0.0;
	    int max_bin = 0;
	    for (unsigned int i=0; i<d_values.size(); i++) {
	       const float &d = d_values[i];
	       float frac = (d - min_d) * inv_range;
	       int bin = frac * static_cast<float>(n_bins);
	       if (bin > max_bin) max_bin = bin;
	       sum += d;
	       sum_sq += d * d;
	       counts[bin] += 1;
	    }
	    std::vector<int> counts_vec(max_bin + 1);
	    std::map<int, int>::const_iterator it;
	    for (it=counts.begin(); it!=counts.end(); ++it)
	       counts_vec[it->first] = it->second;
	    double mean = sum / static_cast<double>(n_points);
	    double var = sum_sq / static_cast<double>(n_points) - mean * mean;
	    if (var < 0.0) var = 0.0;
	    hi = histogram_info_t(min_d, bin_width, counts_vec);
	    hi.mean = mean;
	    hi.variance = var;
	 }
      }
   }
   return hi;

}

#include "coot-utils/diff-diff-map-peaks.hh"

std::vector<std::pair<clipper::Coord_orth, float> >
coot::molecule_t::get_updating_maps_diff_diff_map_peaks(const clipper::Coord_orth &screen_centre) const {

   clipper::Spacegroup sg = xmap.spacegroup();
   clipper::Cell cell = xmap.cell();
   if (false)
      std::cout << "debug:: in get_updating_maps_diff_diff_map_peaks() updating_maps_diff_diff_map_peaks"
                << " has size " << updating_maps_diff_diff_map_peaks.size() << std::endl;
   std::vector<std::pair<clipper::Coord_orth, float> > v2 =
      coot::move_peaks_to_around_position(screen_centre, sg, cell, updating_maps_diff_diff_map_peaks);

   return v2;
}

float
coot::molecule_t::get_density_at_position(const clipper::Coord_orth &pos) const {

   float f = util::density_at_point(xmap, pos);
   return f;
}

texture_as_floats_t
coot::molecule_t:: get_map_section_texture(int section_index, int axis,
                                           float data_value_for_bottom,
                                           float data_value_for_top) const {

   texture_as_floats_t t(xmap, section_index, axis, data_value_for_bottom, data_value_for_top);
   return t;
}


//! @return the number of section in the map along the give axis.
//! (0 for X-axis, 1 for y-axis, 2 for Z-axis).
//! return -1 on failure.
int
coot::molecule_t::get_number_of_map_sections(int axis_id) const {

   int n = -1;
   if (! xmap.is_null()) {
      clipper::Grid_sampling gs = xmap.grid_sampling();
      if (axis_id == 0) n = gs.nu();
      if (axis_id == 1) n = gs.nv();
      if (axis_id == 2) n = gs.nw();
   }
   return n;
}


double
coot::molecule_t::sum_density_for_atoms_in_residue(const std::string &cid,
                                                   const std::vector<std::string> &atom_names,
                                                   const clipper::Xmap<float> &xmap) const {

   double v = 0.0;
   mmdb::Residue *residue_p = cid_to_residue(cid);
   if (residue_p) {
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            std::string atom_name = at->GetAtomName();
            if (std::find(atom_names.begin(), atom_names.end(), atom_name) != atom_names.end()) {
               clipper::Coord_orth pos = co(at);
               float d = util::density_at_point(xmap, pos);
               v += static_cast<double>(d);
            }
         }
      }
   }
   return v;
}

