/* src/molecule-class-info-maps.cc
 *
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2007 by Paul Emsley
 * Copyright 2007, 2008, 2009, 2010 by The University of Oxford
 * Copyright 2013, 2014, 2015, 2016 by Medical Research Council
 *
 * Author: Paul Emsley
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
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

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include <epoxy/gl.h>

#include "compat/coot-sysdep.h"

// Having to set up the include files like this so that
// molecule-class-info.h can be parsed, is silly.

// For stat, mkdir:
#include <iomanip> // for std::setw

// is this a C++11 thing?
#include <functional> // std::ref() for GCC C++11 (not clang)

#include <sys/types.h>
#include <sys/stat.h>

#include <string>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>  // to_string()

#include <mmdb2/mmdb_manager.h>
#include "coords/mmdb-extras.h"
#include "coords/Cartesian.h"
#include "coords/mmdb-crystal.h"
#include "molecule-class-info.h"
#include "coot-utils/coot-coord-utils.hh"
#include "density-contour/CIsoSurface.h"

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

#include "clipper/clipper-cif.h"
#include "clipper/contrib/sfcalc.h"

#include "xmap-utils.h"
#include "coot-utils/xmap-stats.hh"

#include "graphics-info.h"
// #include <GL/glut.h> // needed (only?) for wirecube
#include "globjects.h" // for set_bond_colour()
#include "skeleton/graphical_skel.h"


// #include "coords/mmdb.h"

// for jiggle_fit
#include "coot-utils/coot-map-utils.hh"
#include "coot-utils/coot-map-heavy.hh"
#include "ligand/ligand.hh"


void
molecule_class_info_t::gtk3_draw() {

}

//
void
molecule_class_info_t::sharpen(float b_factor, bool try_gompertz, float gompertz_factor) {

   int n_data = 0;
   int n_tweaked = 0;
   int n_count = 0;
   bool verbose = false;
   bool debugging = false;

   if (debugging) {
      std::cout << "gompertz: " << try_gompertz << " " << gompertz_factor << std::endl;
   }


   bool do_gompertz = false;
   if (try_gompertz) {
      if (original_fobs_sigfobs_filled) {
         do_gompertz = 1;
      } else {
         if (have_sensible_refmac_params) {
            fill_fobs_sigfobs(); // sets original_fobs_sigfobs_filled
            if (have_sensible_refmac_params) {
               if (original_fobs_sigfobs_filled) {
                  do_gompertz = 1;
               } else {
                  std::cout << "WARNING:: Failure to read in F, sigF data" << std::endl;
               }
            }
         }
      }
   }

   if (original_fphis_filled) {

      clipper::HKL_info::HKL_reference_index hri;

      if (debugging)
      std::cout << "DEBUG:: sharpen: using saved " << original_fphis.num_obs()
      << " original data " << std::endl;

      if (debugging) {
         if (do_gompertz) {
            std::cout << "DEBUG:: do_gompertz: " << do_gompertz << " with "
            << original_fobs_sigfobs.num_obs() << " F,sigF reflections"
            << std::endl;
         } else {
            std::cout << "DEBUG:: no gompertz F/sigF scaling " << std::endl;
         }
      }


      if (debugging) {
         for (hri = original_fphis.first(); !hri.last(); hri.next()) {

            if (debugging)
            std::cout << "original_fphis: " << original_fphis[hri].f() << " "
            << hri.invresolsq() << std::endl;
            n_count++;
            if (n_count == 50)
            break;
         }
      }

      clipper::HKL_data< clipper::datatypes::F_phi<float> > fphis(original_fphis.spacegroup(),
      original_fphis.cell(),
      original_fphis.hkl_sampling());
      fphis = original_fphis;

      if (debugging) {
         n_count = 0;
         for (hri = fphis.first(); !hri.last(); hri.next()) {
            if (debugging)
            std::cout << "new fphis: " << fphis[hri].f() << " "
            << hri.invresolsq() << std::endl;
            n_count++;
            if (n_count == 50)
            break;
         }
      }

      if (debugging)
      std::cout << "INFO:: sharpening " << original_fphis.num_obs() << " "
      << fphis.num_obs() << " data " << std::endl;

      n_count = 0;
      int n_gompertz_count = 0;
      double gompertz_sum = 0.0; // for checking values
      for (hri = fphis.first(); !hri.last(); hri.next()) {
         n_data++;

	 // std::cout << " " << hri.invresolsq() << std::endl;

	 float f = fphis[hri].f();
	 if (! clipper::Util::is_nan(f)) {
	    float irs =  hri.invresolsq();
	    if (n_count < 50) {
	       n_count++;
	       if (debugging)
		  std::cout << hri.hkl().format() << " scale factor: e(-" << b_factor
			    << "*" << irs << ") = " << exp(-b_factor * irs)
			    << std::endl;
	    }
	    float gompertz_scale = 1.0;
	    if (do_gompertz) {
	       try {
		  clipper::datatypes::F_sigF<float> fsigf;
		  bool ok = original_fobs_sigfobs.get_data(hri.hkl(), fsigf);
		  if (ok) {
		     if (! clipper::Util::is_nan(fsigf.sigf())) {
			float ratio = fsigf.f()/fsigf.sigf();
			// gompertz function
			// y = ae^{be^{ct}}
			// a is 1.0 obviously.
			float b = -5;
			float c = 0.6586; // such that y is 0.5 for t = 3
			// c = -2; // test value
			gompertz_scale = exp(b*exp(-c*ratio));
			n_gompertz_count++;
			gompertz_sum += gompertz_scale;
		     }
		  }
	       }
	       catch (const clipper::Message_base &exc) {
		  std::cout << "WARNING:: Caught something in sharpen()" << std::endl;
	       }
	    }
	    fphis[hri].f() *= exp(-b_factor * irs * 0.25) * gompertz_scale; // 0.25 factor noted
	                                                                    // by Chang Liu.
	                                                                    // 20130112
	    n_tweaked++;
	 }
      }
      if (do_gompertz) {
	 if (n_gompertz_count)
	    std::cout << "INFO:: Average gompertz scale factor "
		      << gompertz_sum/double(n_gompertz_count)
		      << " from " << n_gompertz_count << " scaled reflections"
		      << std::endl;
	 else
	    std::cout << "WARNING:: no gompertz F/sig correction to reflections!"
		      << std::endl;
      }

      xmap.fft_from(fphis);

      float old_sigma = map_sigma_;
      mean_and_variance<float> mv = map_density_distribution(xmap, 40, false);
      map_mean_  = mv.mean;
      map_sigma_ = sqrt(mv.variance);
      map_max_   = mv.max_density;
      map_min_   = mv.min_density;
      sharpen_b_factor_ = b_factor;

      if (verbose) {
	 std::cout << "      Map mean: ........ " << map_mean_ << std::endl;
	 std::cout << "      Map sigma: ....... " << map_sigma_ << std::endl;
	 std::cout << "      Map maximum: ..... " << map_max_ << std::endl;
	 std::cout << "      Map minimum: ..... " << map_min_ << std::endl;
      }

      // dynamic contour level setting, (not perfect but better than
      // not compensating for the absolute level decreasing).
      if (old_sigma > 0)
	 contour_level *= map_sigma_/old_sigma;

      // update_map_colour_menu_manual(g.n_molecules, name_.c_str());
      // update_map_scroll_wheel_menu_manual(g.n_molecules, name_.c_str());

      update_map();
   }
}

void
molecule_class_info_t::clear_draw_vecs() {

   // crash on double free of the draw vectors. Not sure why. Let's add a lock

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
molecule_class_info_t::clear_diff_map_draw_vecs() {
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

// for negative the other map.
//
void
molecule_class_info_t::set_diff_map_draw_vecs(const coot::CartesianPair* c, int n) {
   // delete [] diff_map_draw_vectors;
   // diff_map_draw_vectors = c; n_diff_map_draw_vectors = n;

   // delete this function
}


void
molecule_class_info_t::update_map() {

   update_map_internal();
}


void
molecule_class_info_t::update_map_internal() {

   float radius = graphics_info_t::box_radius_xray;

   if (has_xmap()) {
      if (is_EM_map())
         radius = graphics_info_t::box_radius_em;

      if (false)
          std::cout << "in update_map_internal() " << radius << " vs x "
                    << graphics_info_t::box_radius_xray << " em "
                    << graphics_info_t::box_radius_em << " is-em: "
                    << is_EM_map() << std::endl;

      coot::Cartesian rc(graphics_info_t::RotationCentre_x(),
                         graphics_info_t::RotationCentre_y(),
                         graphics_info_t::RotationCentre_z());

      update_map_triangles(radius, rc);

      // 20200312-PE do stereo stuff here - removed for now.

   }
}

void
molecule_class_info_t::set_draw_solid_density_surface(bool state) {

   // is this function needed now?

   if (state)
      draw_it_for_map_standard_lines = false;
   else
      draw_it_for_map_standard_lines = true;

   update_map(); // gets solid triangles too.
}


// Create a new combo box for this newly created map.
//
// bleugh.  Using graphics_info_t here!?
//
void
molecule_class_info_t::update_map_in_display_control_widget() const {

   graphics_info_t g;

   std::string dmn = name_for_display_manager();
   if (g.display_control_window())
      display_control_map_combo_box(g.display_control_window(),
				    dmn.c_str(),
				    imol_no);

}


void
molecule_class_info_t::fill_fobs_sigfobs() {

   // set original_fobs_sigfobs_filled when done

   std::cout << "DEBUG:: in fill_fobs_sigfobs() with have_sensible_refmac_params "
            << have_sensible_refmac_params << std::endl;

   if (have_sensible_refmac_params) {

      // only try this once. If you try to import_hkl_data() when the original_fobs_sigfobs
      // already contains data, then crashiness.
      //
      if (! original_fobs_sigfobs_filled && ! original_fobs_sigfobs_fill_tried_and_failed) {

         try {

            std::pair<std::string, std::string> p =
               make_import_datanames(Refmac_fobs_col(), Refmac_sigfobs_col(), "", 0);
            clipper::CCP4MTZfile mtzin;
            mtzin.open_read(Refmac_mtz_filename());
            mtzin.import_hkl_data(original_fobs_sigfobs, p.first);
            mtzin.close_read();
            std::cout << "INFO:: reading " << Refmac_mtz_filename() << " provided "
                      << original_fobs_sigfobs.num_obs() << " data using data name: "
                      << p.first << std::endl;
            if (original_fobs_sigfobs.num_obs() > 10)
               original_fobs_sigfobs_filled = 1;
            else
               original_fobs_sigfobs_fill_tried_and_failed = true;

            // flags

            if (refmac_r_free_flag_sensible) {
               std::string dataname = "/*/*/[" + refmac_r_free_col + "]";
               // if refmac_r_free_col already has /x/y/Rfree - use that instead
               if (refmac_r_free_col.length() > 0)
                  if (refmac_r_free_col[0] == '/') {
                     dataname = refmac_r_free_col;
                     dataname = "/*/*/[" + coot::util::file_name_non_directory(refmac_r_free_col) + "]";
                  }
               std::cout << "About to read " << Refmac_mtz_filename() << " with dataname " << dataname << std::endl;
               mtzin.open_read(Refmac_mtz_filename());
               mtzin.import_hkl_data(original_r_free_flags, dataname);
               mtzin.close_read();

               std::cout << "INFO:: reading " << Refmac_mtz_filename()
                         << " using dataname: " << dataname << " provided "
                         << original_r_free_flags.num_obs() << " R-free flags\n";
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
      }
   } else {
      std::cout << "DEBUG:: fill_fobs_sigfobs() no Fobs parameters\n";
   }
}



// not a member of the class because of the burden it puts on the header: CIsoSurface is not
// needed to compile main.cc (or should not be)

#include "gensurf.hh"

//
void
molecule_class_info_t::update_map_triangles(float radius, coot::Cartesian centre) {

   // cut glass mode means
   // do_solid_surface_for_density
   // do_flat_shading_for_solid_density_surface is true

   CIsoSurface<float> my_isosurface;
   coot::CartesianPairInfo v;
   int isample_step = 1;
   graphics_info_t g;

   bool is_em_map = false;
   if (is_em_map_cached_state() == 1) {
      is_em_map = true;
   }

   // std::cout   << "DEBUG:: g.zoom: " << g.zoom << std::endl;

   if (g.dynamic_map_resampling == 1)
      // isample_step = 1 + int (0.009*g.zoom);
      isample_step = 1 + int (0.009*(g.zoom + g.dynamic_map_zoom_offset));

   if (isample_step > 15)
      isample_step = 15;

   // for critical points of size display and resampling being different:
   //
   float dy_radius = radius;
   if (g.dynamic_map_size_display == 1) {
      if (isample_step <= 15 )
         dy_radius *= float(isample_step);
      else
         dy_radius *= 15.0;
   }

   //
   if (isample_step <= 0) {
      std::cout << "WARNING:: Bad zoom   ("<< g.zoom
                << "):  setting isample_step to 1" << std::endl;
      isample_step = 1;
   }
   if (dy_radius <= 0.0) {
      std::cout << "WARNING:: Bad radius (" << dy_radius
                << ") setting to 10" << std::endl;
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

      // post_process_map_triangles();

      setup_glsl_map_rendering(); // turn tri_con into buffers.

      std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > vp =
         make_map_mesh();

      // Mesh gm(vp);
      // meshes.push_back(gm);
      // meshes.back().setup() here?                          
      

/*
Threaded difference map lines:
if (xmap_is_diff_map) { // do the negative level
   clear_diff_map_draw_vecs();
   std::vector<std::thread> threads;
   int n_reams = coot::get_max_number_of_threads() - 1;
   if (n_reams < 1) n_reams = 1;
   for (int ii=0; ii<n_reams; ii++) {
      int iream_start = ii;
      threads.push_back(std::thread(gensurf_and_add_vecs_threaded_workpackage,
                                    &xmap, -contour_level, dy_radius, centre,
                                    isample_step, iream_start, n_reams, is_em_map,
                                    &draw_diff_map_vector_sets));
   }
   for (int ii=0; ii<n_reams; ii++)
      threads[ii].join();
}

if (draw_it_for_solid_density_surface) {
   tri_con = my_isosurface.GenerateTriangles_from_Xmap(xmap,
                           contour_level, dy_radius, centre,
                           isample_step);


*/


/* single threaded triangles
      if (true) {
         std::cout << "calling my_isosurface.GenerateTriangles_from_Xmap()" << std::endl;
         tri_con = my_isosurface.GenerateTriangles_from_Xmap(xmap,
               contour_level, dy_radius, centre, isample_step, iream_start, n_reams, is_em_map);
         std::cout << "done my_isosurface.GenerateTriangles_from_Xmap()" << std::endl;

         if (xmap_is_diff_map) {
            tri_con_diff_map_neg = my_isosurface.GenerateTriangles_from_Xmap(xmap,
                  -contour_level,
                  dy_radius, centre,
                  isample_step);
         }
         setup_glsl_map_rendering(); // turn tri_con into buffers.
      }

*/
   }
}

void gensurf_and_add_vecs_threaded_workpackage(const clipper::Xmap<float> *xmap_p,
					       float contour_level, float dy_radius,
					       coot::Cartesian centre,
					       int isample_step,
					       int iream_start, int n_reams,
					       bool is_em_map,
					       std::vector<coot::density_contour_triangles_container_t> *draw_vector_sets_p) {

   try {
      CIsoSurface<float> my_isosurface;

      coot::density_contour_triangles_container_t tri_con =
        my_isosurface.GenerateTriangles_from_Xmap(std::cref(*xmap_p),
                                                  contour_level, dy_radius, centre, isample_step,
                                                  iream_start, n_reams, is_em_map);

      // we are about to put the triangles into draw_vectors, so get the lock to
      // do that, so that the threads don't try to change draw_vectors at the same time.
      //
      bool unlocked = false;
      while (! molecule_class_info_t::draw_vector_sets_lock.compare_exchange_weak(unlocked, true)) {
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

      molecule_class_info_t::draw_vector_sets_lock = false; // unlock
   }
   catch (const std::out_of_range &oor) {
      std::cout << "ERROR:: contouring threaded workpackage " << oor.what() << std::endl;
   }
}


void
molecule_class_info_t::post_process_map_triangles() {

   // average the normals of the vertices that are close.
   // note std::vector<coot::density_contour_triangles_container_t> draw_vector_sets;

   double min_dist = 0.03;
   double min_dist_sqrd = min_dist * min_dist;
   unsigned int n_reset = 0;

   for (unsigned int i=0; i<draw_vector_sets.size(); i++) {
      coot::density_contour_triangles_container_t &tri_con_i = draw_vector_sets[i];
      for (unsigned int ii=0; ii<tri_con_i.points.size(); ii++) {
         const clipper::Coord_orth &pt_i = tri_con_i.points[ii];
         std::vector<clipper::Coord_orth> neighb_normals;
         for (unsigned int j=0; j<draw_vector_sets.size(); j++) {
            const coot::density_contour_triangles_container_t &tri_con_j = draw_vector_sets[j];
            for (unsigned int jj=0; jj<tri_con_j.points.size(); jj++) {
               if (i == j && ii == jj ) continue;
               const clipper::Coord_orth &pt_j = tri_con_j.points[jj];
               double dd = (pt_i-pt_j).lengthsq();
               if (dd < min_dist_sqrd) {
                  neighb_normals.push_back(tri_con_j.normals[jj]);
               }
            }
         }
         if (! neighb_normals.empty()) {
            clipper::Coord_orth sum = tri_con_i.normals[ii];
            for (unsigned int in=0; in<neighb_normals.size(); in++)
               sum += neighb_normals[in];
            tri_con_i.normals[ii] = clipper::Coord_orth(sum.unit());
            n_reset++;
         }
      }
   }
   
   std::cout << "DEBUG:: n_reset " << n_reset << std::endl;

}

void
molecule_class_info_t::setup_map_cap(Shader *shader_p,
                                     const clipper::Coord_orth &base_point, // Bring it into this class.
                                     const clipper::Coord_orth &x_axis_uv, // Of the cap plane, of course.
                                     const clipper::Coord_orth &y_axis_uv,
                                     double x_axis_step_size,
                                     double y_axis_step_size,
                                     unsigned int n_x_axis_points,
                                     unsigned int n_y_axis_points) {

   // this line is completely vital! - But do I want gtk_gl_area_attach_buffers(gl_area) ?
   //
   gtk_gl_area_make_current(GTK_GL_AREA(graphics_info_t::glareas[0]));

   GLenum err = glGetError(); if (err) std::cout << "error in setup_map_cap() -- start -- " << err << std::endl;
   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > map_cap =
      make_map_cap(base_point, x_axis_uv, y_axis_uv, x_axis_step_size, y_axis_step_size,
                   n_x_axis_points, n_y_axis_points);

   shader_p->Use();
   Material material;

   // What was I trying to do here?
   // graphical_molecule gm;
   // graphical_molecule gm_cap(map_cap.first, map_cap.second);
   // gm.setup_simple_triangles(shader_p, material);
   // graphical_molecules.push_back(gm);
   // graphical_molecules.push_back(gm_cap);

   Mesh gm_cap(map_cap);
   meshes.push_back(gm_cap);
   meshes.back().setup(shader_p, material);
   
}


// there are called molecular meshes now
void
molecule_class_info_t::mesh_draw_normals(const glm::mat4 &mvp) {

   bool do_all = false;
   // there are molecular meshes now
   const float s = 0.1;
   if (do_all) {
      for (unsigned int i=0; i<meshes.size(); i++) {
         meshes[i].draw_normals(mvp, s);
      }
   } else {
      meshes.back().draw_normals(mvp, s);
   }
}


std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> >
molecule_class_info_t::make_map_mesh() {

   std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > vp;
   std::vector<s_generic_vertex> &vertices = vp.first;
   std::vector<g_triangle> &triangles = vp.second;

   std::vector<coot::density_contour_triangles_container_t>::const_iterator it;
   for (it=draw_vector_sets.begin(); it!=draw_vector_sets.end(); it++) {
      const coot::density_contour_triangles_container_t &tri_con(*it);
      // vertices
      int idx_base = vertices.size();
      for (unsigned int i=0; i<tri_con.points.size(); i++) {
         glm::vec3 p( tri_con.points[i].x(),  tri_con.points[i].y(),  tri_con.points[i].z());
         glm::vec3 n(tri_con.normals[i].x(), tri_con.normals[i].y(), tri_con.normals[i].z());
         glm::vec4 c(0.5, 0.5, 0.5, 1.0);
         s_generic_vertex g(p,n,c);
         vertices.push_back(g);
      }
      // triangles
      for (unsigned int i=0; i<tri_con.point_indices.size(); i++) {
         g_triangle t(tri_con.point_indices[i].pointID[0] + idx_base,
                      tri_con.point_indices[i].pointID[1] + idx_base,
                      tri_con.point_indices[i].pointID[2] + idx_base);
         if (triangles.size() < 10000)
            triangles.push_back(t);
      }
   }
   return vp;
}



void
molecule_class_info_t::sort_map_triangles(const clipper::Coord_orth &eye_position) {

   bool do_the_sort = false;
   if ((eye_position - previous_eye_position).lengthsq() > 0.0001) do_the_sort = true;

   if (! do_the_sort) return;

   for (unsigned int i=0; i<map_triangle_centres.size(); i++) {
      clipper::Coord_orth delta(map_triangle_centres[i].second.mid_point - eye_position);
      double dd = delta.lengthsq();
      map_triangle_centres[i].second.back_front_projection_distance = dd;
   }

   // this sign needs checking.
   auto map_triangle_sorter = [](const std::pair<int, TRIANGLE> &t1,
                                 const std::pair<int, TRIANGLE> &t2) {
                                 return (t1.second.back_front_projection_distance > t2.second.back_front_projection_distance);
                              };

   std::sort(map_triangle_centres.begin(), map_triangle_centres.end(), map_triangle_sorter);

   unsigned int n_triangle_centres = map_triangle_centres.size();

   int *indices_for_triangles = new int[3 * n_triangle_centres];
   for (unsigned int i=0; i<map_triangle_centres.size(); i++) {
      indices_for_triangles[3*i  ] = map_triangle_centres[i].second.pointID[0];
      indices_for_triangles[3*i+1] = map_triangle_centres[i].second.pointID[1];
      indices_for_triangles[3*i+2] = map_triangle_centres[i].second.pointID[2];
   }

   GLenum err = glGetError();
   // glDeleteBuffers(1, &m_IndexBuffer_for_map_triangles_ID); // doing this causes weird effects.
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_IndexBuffer_for_map_triangles_ID);
   err = glGetError(); if (err) std::cout << "GL error: sorting triangles: " << err << std::endl;
   // glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * n_indices_for_triangles,
   // &indices_for_triangles[0], GL_STATIC_DRAW);
   // perhaps don't use GL_STATIC_DRAW - use GL_DYNAMIC_DRAW
   glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * n_indices_for_triangles,
                &indices_for_triangles[0], GL_DYNAMIC_DRAW);
   err = glGetError(); if (err) std::cout << "GL error: sorting triangles: " << err << std::endl;

   delete [] indices_for_triangles;

   // for next time
   previous_eye_position = eye_position;
}



void
molecule_class_info_t::setup_glsl_map_rendering() {

   // This is called from update_map_triangles().

   bool debug = false;

   if (true) { // real map

      unsigned int sum_tri_con_points = 0;
      unsigned int sum_tri_con_normals = 0;
      unsigned int sum_tri_con_triangles = 0;
      std::vector<coot::density_contour_triangles_container_t>::const_iterator it;

      auto tp_0 = std::chrono::high_resolution_clock::now();
      for (it=draw_vector_sets.begin(); it!=draw_vector_sets.end(); it++) {
         const coot::density_contour_triangles_container_t &tri_con(*it);

         if (false) {
            std::cout << "------------------ setup_glsl_map_rendering() here A ------------" << std::endl;
            std::cout << "------------------ setup_glsl_map_rendering() here B ------------"
                      << tri_con.point_indices.size() << " triangles" << std::endl;
            std::cout << "------------------ setup_glsl_map_rendering() here C ------------"
                      << tri_con.normals.size() << " normals" << std::endl;
         }

         sum_tri_con_points    += tri_con.points.size();
         sum_tri_con_normals   += tri_con.normals.size();
         sum_tri_con_triangles += tri_con.point_indices.size();
      }
      if (xmap_is_diff_map) {
         for (it=draw_diff_map_vector_sets.begin(); it!=draw_diff_map_vector_sets.end(); it++) {
            const coot::density_contour_triangles_container_t &tri_con(*it);
            sum_tri_con_points    += tri_con.points.size();
            sum_tri_con_normals   += tri_con.normals.size();
            sum_tri_con_triangles += tri_con.point_indices.size();
         }
      }

      if (sum_tri_con_triangles > 0) {

         std::vector<int> idx_base_for_points_vec(   draw_vector_sets.size()); // and normals
         std::vector<int> idx_base_for_triangles_vec(draw_vector_sets.size());
         if (xmap_is_diff_map) {
            idx_base_for_points_vec.resize(draw_vector_sets.size() + draw_diff_map_vector_sets.size());
            idx_base_for_triangles_vec.resize(draw_vector_sets.size() + draw_diff_map_vector_sets.size());
         }
         for (unsigned int i=0; i<draw_vector_sets.size(); i++) {
            if (i==0) {
               idx_base_for_points_vec[i] = 0;
               idx_base_for_triangles_vec[i] = 0;
            } else {
               idx_base_for_points_vec[i]    = idx_base_for_points_vec[i-1]    + draw_vector_sets[i-1].points.size();
               idx_base_for_triangles_vec[i] = idx_base_for_triangles_vec[i-1] + draw_vector_sets[i-1].points.size();
            }
         }
         if (xmap_is_diff_map) {
            for (unsigned int i0=0; i0<draw_diff_map_vector_sets.size(); i0++) {
               unsigned int i = i0 + draw_vector_sets.size();
               const coot::density_contour_triangles_container_t &tri_con(draw_diff_map_vector_sets[i0]);
               try {
                  if (i0==0) {
                     idx_base_for_points_vec.at(i)    = idx_base_for_points_vec.at(i-1)    + draw_vector_sets.at(i-1).points.size();
                     idx_base_for_triangles_vec.at(i) = idx_base_for_triangles_vec.at(i-1) + draw_vector_sets.at(i-1).points.size();
                  } else {
                     idx_base_for_points_vec.at(i)    = idx_base_for_points_vec.at(i-1)    + draw_diff_map_vector_sets.at(i0-1).points.size();
                     idx_base_for_triangles_vec.at(i) = idx_base_for_triangles_vec.at(i-1) + draw_diff_map_vector_sets.at(i0-1).points.size();
                  }
               }
               catch (std::out_of_range &oor) {
                  std::cout << "ERROR:: caught out of range " << oor.what() << " at i " << i << std::endl;
               }
            }
         }

         // put the triangle mid-points into a single vector. Not sure that this is needed now.
         // There's a different mechanism (in the works) for triangle sorting.
         int idx_for_mid_points = 0;
         // std::vector<std::pair<int, TRIANGLE> > map_triangle_centres(sum_tri_con_triangles);
         map_triangle_centres.resize(sum_tri_con_triangles); // class member
         for (unsigned int i=0; i<draw_vector_sets.size(); i++) {
            const coot::density_contour_triangles_container_t &tri_con(draw_vector_sets[i]);
            int idx_base_for_triangles = idx_base_for_triangles_vec[i];
            for (std::size_t j=0; j<tri_con.point_indices.size(); j++) {
               map_triangle_centres[idx_for_mid_points].first = idx_for_mid_points;
               map_triangle_centres[idx_for_mid_points].second.mid_point = tri_con.point_indices[j].mid_point;
               map_triangle_centres[idx_for_mid_points].second.pointID[0] = idx_base_for_triangles + tri_con.point_indices[j].pointID[0];
               map_triangle_centres[idx_for_mid_points].second.pointID[1] = idx_base_for_triangles + tri_con.point_indices[j].pointID[1];
               map_triangle_centres[idx_for_mid_points].second.pointID[2] = idx_base_for_triangles + tri_con.point_indices[j].pointID[2];
               idx_for_mid_points++;
            }
         }

         // transfer the points
         float *points = new float[3 * sum_tri_con_points];

         int idx_points = 0;
         for (unsigned int i=0; i<draw_vector_sets.size(); i++) {
            const coot::density_contour_triangles_container_t &tri_con(draw_vector_sets[i]);
            for (std::size_t j=0; j<tri_con.points.size(); j++) {
               points[3*idx_points  ] = tri_con.points[j].x();
               points[3*idx_points+1] = tri_con.points[j].y();
               points[3*idx_points+2] = tri_con.points[j].z();
               idx_points++;
            }
         }
         if (xmap_is_diff_map) {
            for (unsigned int i0=0; i0<draw_diff_map_vector_sets.size(); i0++) {
               unsigned int i = i0 + draw_vector_sets.size();
               const coot::density_contour_triangles_container_t &tri_con(draw_diff_map_vector_sets[i0]);
               for (std::size_t j=0; j<tri_con.points.size(); j++) {
                  points[3*idx_points  ] = tri_con.points[j].x();
                  points[3*idx_points+1] = tri_con.points[j].y();
                  points[3*idx_points+2] = tri_con.points[j].z();
                  idx_points++;
               }
            }
         }

         n_vertices_for_map_VertexArray = 6 * sum_tri_con_triangles;
         n_indices_for_triangles    = 3 * sum_tri_con_triangles;
         n_indices_for_lines        = 6 * sum_tri_con_triangles;

         // "in-line" class
         class hashed_line_t {
         public:
            long hash;
            hashed_line_t(const unsigned int &i1, const unsigned int &i2) {
               if (i1 < i2)
                  hash = i1 + 1000000 * i2;
               else
                  hash = i2 + 1000000 * i1;
            }
            bool operator<(const hashed_line_t &h) const {
               return hash < h.hash;
            }
         };

         std::set<hashed_line_t> line_hash;
         unsigned int n_lines_added = 0;
         unsigned int n_lines_total = 0;

         int *indices_for_lines = new int[n_indices_for_lines];
         int idx_for_indices = 0;
         bool use_line_hash = false; // make user-settable parameter

         if (xmap_is_diff_map) use_line_hash = false; // too tricky indexing (for now)

         if (use_line_hash) {

            for (unsigned int i=0; i<draw_vector_sets.size(); i++) {
               const coot::density_contour_triangles_container_t &tri_con(draw_vector_sets[i]);
               int idx_base_for_points = idx_base_for_points_vec[i];
               for (std::size_t j=0; j<tri_con.point_indices.size(); j++) {

                  bool add_01 = true;
                  bool add_12 = true;
                  bool add_20 = true;
                  hashed_line_t hash_key_1 = hashed_line_t(idx_base_for_points + tri_con.point_indices[j].pointID[0],
                                                           idx_base_for_points + tri_con.point_indices[j].pointID[1]);
                  hashed_line_t hash_key_2 = hashed_line_t(idx_base_for_points + tri_con.point_indices[j].pointID[1],
                                                           idx_base_for_points + tri_con.point_indices[j].pointID[2]);
                  hashed_line_t hash_key_3 = hashed_line_t(idx_base_for_points + tri_con.point_indices[j].pointID[2],
                                                           idx_base_for_points + tri_con.point_indices[j].pointID[0]);

                  if (line_hash.find(hash_key_1) != line_hash.end()) add_01 = false; else line_hash.insert(hash_key_1);
                  if (line_hash.find(hash_key_2) != line_hash.end()) add_12 = false; else line_hash.insert(hash_key_2);
                  if (line_hash.find(hash_key_3) != line_hash.end()) add_20 = false; else line_hash.insert(hash_key_3);

                  if (add_01) {
                     indices_for_lines[idx_for_indices  ] = idx_base_for_points + tri_con.point_indices[j].pointID[0];
                     indices_for_lines[idx_for_indices+1] = idx_base_for_points + tri_con.point_indices[j].pointID[1];
                     idx_for_indices += 2;
                     n_lines_added++;
                  }
                  if (add_12) {
                     indices_for_lines[idx_for_indices  ] = idx_base_for_points + tri_con.point_indices[j].pointID[1];
                     indices_for_lines[idx_for_indices+1] = idx_base_for_points + tri_con.point_indices[j].pointID[2];
                     idx_for_indices += 2;
                     n_lines_added++;
                  }
                  if (add_20) {
                     indices_for_lines[idx_for_indices  ] = idx_base_for_points + tri_con.point_indices[j].pointID[2];
                     indices_for_lines[idx_for_indices+1] = idx_base_for_points + tri_con.point_indices[j].pointID[0];
                     idx_for_indices += 2;
                     n_lines_added++;
                  }
                  n_lines_total += 3;
               }
            }
            std::cout << "debug:: n_lines_total " << n_lines_total << " of which " << n_lines_added
                      << " were added " << std::endl;
         } else {

            for (unsigned int i=0; i<draw_vector_sets.size(); i++) {
               const coot::density_contour_triangles_container_t &tri_con(draw_vector_sets[i]);
               int idx_base_for_points = idx_base_for_points_vec[i];
               for (std::size_t j=0; j<tri_con.point_indices.size(); j++) {
                  indices_for_lines[6*idx_for_indices  ] = idx_base_for_points + tri_con.point_indices[j].pointID[0];
                  indices_for_lines[6*idx_for_indices+1] = idx_base_for_points + tri_con.point_indices[j].pointID[1];
                  indices_for_lines[6*idx_for_indices+2] = idx_base_for_points + tri_con.point_indices[j].pointID[1];
                  indices_for_lines[6*idx_for_indices+3] = idx_base_for_points + tri_con.point_indices[j].pointID[2];
                  indices_for_lines[6*idx_for_indices+4] = idx_base_for_points + tri_con.point_indices[j].pointID[2];
                  indices_for_lines[6*idx_for_indices+5] = idx_base_for_points + tri_con.point_indices[j].pointID[0];
                  idx_for_indices++;
               }
            }
            if (xmap_is_diff_map) {
               for (unsigned int i0=0; i0<draw_diff_map_vector_sets.size(); i0++) {
                  unsigned int i = i0 + draw_vector_sets.size();
                  const coot::density_contour_triangles_container_t &tri_con(draw_diff_map_vector_sets[i0]);
                  int idx_base_for_points = idx_base_for_points_vec[i];
                  for (std::size_t j=0; j<tri_con.point_indices.size(); j++) {
                     indices_for_lines[6*idx_for_indices  ] = idx_base_for_points + tri_con.point_indices[j].pointID[0];
                     indices_for_lines[6*idx_for_indices+1] = idx_base_for_points + tri_con.point_indices[j].pointID[1];
                     indices_for_lines[6*idx_for_indices+2] = idx_base_for_points + tri_con.point_indices[j].pointID[1];
                     indices_for_lines[6*idx_for_indices+3] = idx_base_for_points + tri_con.point_indices[j].pointID[2];
                     indices_for_lines[6*idx_for_indices+4] = idx_base_for_points + tri_con.point_indices[j].pointID[2];
                     indices_for_lines[6*idx_for_indices+5] = idx_base_for_points + tri_con.point_indices[j].pointID[0];
                     idx_for_indices++;
                  }
               }
            }
         }

         int *indices_for_triangles = new int[n_indices_for_triangles];
         int idx_for_triangles = 0;

         if (true) {
            for (unsigned int i=0; i<draw_vector_sets.size(); i++) {
               const coot::density_contour_triangles_container_t &tri_con(draw_vector_sets[i]);
               int idx_base_for_triangles = idx_base_for_triangles_vec[i];
               for (std::size_t i=0; i<tri_con.point_indices.size(); i++) {
                  indices_for_triangles[3*idx_for_triangles  ] = idx_base_for_triangles + tri_con.point_indices[i].pointID[0];
                  indices_for_triangles[3*idx_for_triangles+1] = idx_base_for_triangles + tri_con.point_indices[i].pointID[1];
                  indices_for_triangles[3*idx_for_triangles+2] = idx_base_for_triangles + tri_con.point_indices[i].pointID[2];
                  idx_for_triangles++;

                  // as we install these we should find out where the triangles centres are (or maybe before now)
                  // so that we can use the centres for sorting.
               }
            }
            if (xmap_is_diff_map) {
               for (unsigned int i0=0; i0<draw_diff_map_vector_sets.size(); i0++) {
                  unsigned int i = i0 + draw_vector_sets.size();
                  const coot::density_contour_triangles_container_t &tri_con(draw_diff_map_vector_sets[i0]);
                  int idx_base_for_triangles = idx_base_for_triangles_vec[i];
                  for (std::size_t i=0; i<tri_con.point_indices.size(); i++) {
                     indices_for_triangles[3*idx_for_triangles  ] = idx_base_for_triangles + tri_con.point_indices[i].pointID[0];
                     indices_for_triangles[3*idx_for_triangles+1] = idx_base_for_triangles + tri_con.point_indices[i].pointID[1];
                     indices_for_triangles[3*idx_for_triangles+2] = idx_base_for_triangles + tri_con.point_indices[i].pointID[2];
                     idx_for_triangles++;
                  }
               }
            }
         }

         // each index has a normal

         int n_normals = sum_tri_con_normals;
         float *normals = new float[3 * n_normals];
         int idx_for_normals = 0;
         for (unsigned int i0=0; i0<draw_vector_sets.size(); i0++) {
            const coot::density_contour_triangles_container_t &tri_con(draw_vector_sets[i0]);
            for (std::size_t i=0; i<tri_con.normals.size(); i++) {
               // std::cout << "real normal " << idx_for_normals << " " << tri_con.normals[i].format() << std::endl;
               normals[3*idx_for_normals  ] = tri_con.normals[i].x();
               normals[3*idx_for_normals+1] = tri_con.normals[i].y();
               normals[3*idx_for_normals+2] = tri_con.normals[i].z();
               idx_for_normals++;
            }
         }
         if (xmap_is_diff_map) {
            for (unsigned int i0=0; i0<draw_diff_map_vector_sets.size(); i0++) {
               const coot::density_contour_triangles_container_t &tri_con(draw_diff_map_vector_sets[i0]);
               for (std::size_t i=0; i<tri_con.normals.size(); i++) {
                  normals[3*idx_for_normals  ] = tri_con.normals[i].x();
                  normals[3*idx_for_normals+1] = tri_con.normals[i].y();
                  normals[3*idx_for_normals+2] = tri_con.normals[i].z();
                  idx_for_normals++;
               }
            }
         }

         // each vertex/point has a colour

         int n_colours = sum_tri_con_points;
         float *colours = new float[4 * n_colours];
         int idx_for_colours = 0;
         for (unsigned int i=0; i<draw_vector_sets.size(); i++) {
            const coot::density_contour_triangles_container_t &tri_con(draw_vector_sets[i]);
            // check here for difference map colours
            for (std::size_t j=0; j<tri_con.points.size(); j++) {
               if (radial_map_colouring_do_radial_colouring) {
                  if (idx_for_colours < n_colours) {
                     // Oh dear - indexing!
                     int idx_base = 3 * idx_for_colours;
                     clipper::Coord_orth co(points[idx_base], points[idx_base+1], points[idx_base+2]);
                     double dd = (co-radial_map_colour_centre).lengthsq();
                     double r = sqrt(dd);
                     GdkRGBA map_col = radius_to_colour(r, radial_map_colour_radius_min, radial_map_colour_radius_max);
                     colours[4*idx_for_colours  ] = map_col.red;
                     colours[4*idx_for_colours+1] = map_col.green;
                     colours[4*idx_for_colours+2] = map_col.blue;
                     colours[4*idx_for_colours+3] = 1.0f;
                  }
               } else {
                  if (colour_map_using_other_map_flag) {
                     int idx_base = 3 * idx_for_colours;
                     clipper::Coord_orth co(points[idx_base], points[idx_base+1], points[idx_base+2]);
                     GdkRGBA map_col = position_to_colour_using_other_map(co);
                     colours[4*idx_for_colours  ] = map_col.red;
                     colours[4*idx_for_colours+1] = map_col.green;
                     colours[4*idx_for_colours+2] = map_col.blue;
                     colours[4*idx_for_colours+3] = 1.0f;
                  } else {
                     // basic/standard single colour map
                     if ((4*idx_for_colours) < (4 * n_colours)) {
                        colours[4*idx_for_colours  ] = map_colour.red;
                        colours[4*idx_for_colours+1] = map_colour.green;
                        colours[4*idx_for_colours+2] = map_colour.blue;
                        colours[4*idx_for_colours+3] = 1.0;
                     } else {
                        std::cout << "oops indexing error for colours"
                                  << idx_for_colours << " " << n_colours << std::endl;
                     }
                  }
               }
               idx_for_colours++;
            }
         }
         if (xmap_is_diff_map) {
            for (unsigned int i0=0; i0<draw_diff_map_vector_sets.size(); i0++) {
               unsigned int i = i0 + draw_vector_sets.size();
               const coot::density_contour_triangles_container_t &tri_con(draw_diff_map_vector_sets[i0]);
               for (std::size_t j=0; j<tri_con.points.size(); j++) {
                  if ((4*idx_for_colours) < (4 * n_colours)) {
                     colours[4*idx_for_colours  ] = map_colour_negative_level.red;
                     colours[4*idx_for_colours+1] = map_colour_negative_level.green;
                     colours[4*idx_for_colours+2] = map_colour_negative_level.blue;
                     colours[4*idx_for_colours+3] = 1.0;;
                     idx_for_colours++;
                  } else {
                     std::cout << "oops indexing error for colours"
                               << idx_for_colours << " " << n_colours << std::endl;
                  }
               }
            }
         }

         // why is this needed? Is it needed? Done by caller? (or should be?)
         gtk_gl_area_make_current(GTK_GL_AREA(graphics_info_t::glareas[0]));

         glGenVertexArrays(1, &m_VertexArrayID_for_map);
         GLenum err = glGetError();
         if (err) std::cout << "############## in setup_glsl_map_rendering() error trying to bind vertex array "
                            << m_VertexArrayID_for_map << std::endl;
         glBindVertexArray(m_VertexArrayID_for_map);
         err = glGetError();
         if (err) std::cout << "############## in setup_glsl_map_rendering() error glBindVertexArray() "
                            << m_VertexArrayID_for_map << std::endl;
         err = glGetError();

         // positions
         glGenBuffers(1, &m_VertexBufferID);
         err = glGetError();
         // std::cout << "setup_glsl_map_rendering() glGenBuffers() err " << err << std::endl;
         glBindBuffer(GL_ARRAY_BUFFER, m_VertexBufferID);
         err = glGetError();
         // std::cout << "setup_glsl_map_rendering() glBindBuffer() err " << err << std::endl;
         int n_bytes_for_tri_con_points = sizeof(float) * 3 * sum_tri_con_points;
         // std::cout << "debug:: n_bytes_for_tri_con_points " << n_bytes_for_tri_con_points << std::endl;
         glBufferData(GL_ARRAY_BUFFER, n_bytes_for_tri_con_points, &points[0], GL_STATIC_DRAW);
         err = glGetError();
         // std::cout << "setup_glsl_map_rendering() glBufferData() err " << err << std::endl;
         glEnableVertexAttribArray(0);
         err = glGetError();
         // std::cout << "setup_glsl_map_rendering() glEnableVertexAttribArray() err " << err << std::endl;
         glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
         err = glGetError();
         // std::cout << "setup_glsl_map_rendering() glVertexAttribPointer() err " << err << std::endl;

         // normals
         glGenBuffers(1, &m_NormalBufferID);
         err = glGetError();
         // std::cout << "setup_glsl_map_rendering() glGenBuffers() for normals err " << err << std::endl;
         glBindBuffer(GL_ARRAY_BUFFER, m_NormalBufferID);
         err = glGetError();
         // std::cout << "setup_glsl_map_rendering() glBindBuffer() err " << err << std::endl;
         glBufferData(GL_ARRAY_BUFFER, n_normals * 3 * sizeof(float), &normals[0], GL_STATIC_DRAW);
         err = glGetError();
         // std::cout << "setup_glsl_map_rendering() glBufferData() for normals err " << err << std::endl;
         glEnableVertexAttribArray(1);
         err = glGetError();
         // std::cout << "setup_glsl_map_rendering() glEnableVertexAttribArray() err normals " << err << std::endl;
         glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
         err = glGetError();
         // std::cout << "setup_glsl_map_rendering() glVertexAttribPointer() err " << err << std::endl;


         // colours
         glGenBuffers(1, &m_ColourBufferID);
         err = glGetError();
         // std::cout << "setup_glsl_map_rendering() glGenBuffers() for colours err " << err << std::endl;
         glBindBuffer(GL_ARRAY_BUFFER, m_ColourBufferID);
         err = glGetError();
         // std::cout << "setup_glsl_map_rendering() glBindBuffer() err " << err << std::endl;
         glBufferData(GL_ARRAY_BUFFER, n_colours * 4 * sizeof(float), &colours[0], GL_STATIC_DRAW);
         err = glGetError();
         // std::cout << "setup_glsl_map_rendering() glBufferData() for colours err " << err << std::endl;
         glEnableVertexAttribArray(2);
         err = glGetError();
         // std::cout << "setup_glsl_map_rendering() glEnableVertexAttribArray() err colours " << err << std::endl;
         glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, 0, 0);
         err = glGetError();
         // std::cout << "setup_glsl_map_rendering() glVertexAttribPointer() err " << err << std::endl;


         // indices for map lines
         glGenBuffers(1, &m_IndexBuffer_for_map_lines_ID);
         err = glGetError();
         // std::cout << "setup_glsl_map_rendering() glGenBuffers() " << err << std::endl;
         glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_IndexBuffer_for_map_lines_ID);
         err = glGetError();
         // std::cout << "setup_glsl_map_rendering() glBindBuffer() " << err << std::endl;
         // Note to self: is n_vertices_for_VertexArray correct here?

         if (use_line_hash)
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * n_lines_added * 2,
                         &indices_for_lines[0], GL_STATIC_DRAW);
         else
            glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * n_indices_for_lines,
                         &indices_for_lines[0], GL_STATIC_DRAW);
         err = glGetError();
         // std::cout << "setup_glsl_map_rendering() glBufferData() " << err << std::endl;


         // indices for map triangles

         glGenBuffers(1, &m_IndexBuffer_for_map_triangles_ID);
         err = glGetError();
         // std::cout << "setup_glsl_map_rendering() glGenBuffers() for m_IndexBuffer_for_triangles_ID " << err << std::endl;
         glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_IndexBuffer_for_map_triangles_ID);
         err = glGetError();
         // std::cout << "setup_glsl_map_rendering() glBindBuffer() for triangles " << err << std::endl;
         glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * n_indices_for_triangles,
                      &indices_for_triangles[0], GL_STATIC_DRAW);
         err = glGetError();
         // std::cout << "setup_glsl_map_rendering() glBufferData() " << err << std::endl;

         delete [] points;
         delete [] indices_for_lines;
         delete [] indices_for_triangles;
         delete [] colours;

      }
      auto tp_1 = std::chrono::high_resolution_clock::now();
      auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
      // std::cout << "INFO:: Map triangle generation time " << d10 << " milliseconds" << std::endl;
   }

   if (false)
      std::cout << "------------------ setup_glsl_map_rendering() here -end- ------------ "
                << n_vertices_for_map_VertexArray << std::endl;
}


GdkRGBA
molecule_class_info_t::position_to_colour_using_other_map(const clipper::Coord_orth &position) {

   GdkRGBA c;
   c.red   = 0.0;
   c.green = 0.1;
   c.blue  = 0.0;
   c.alpha = 1.0;

   if (other_map_for_colouring_p) {
      float min_value = other_map_for_colouring_min_value;
      float max_value = other_map_for_colouring_max_value;
      float dv = coot::util::density_at_point(*other_map_for_colouring_p, position);

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

      c = fraction_to_colour(f);
      // std::cout << "fraction " << f << " col " << c.red << " " << c.green << " " << c.blue << std::endl;
   } else {
      return c;
   }

#if 0  // testing function

   // float v = coot::util::random()/static_cast<float>(RAND_MAX);

   clipper::Coord_orth pt_central(0,0,-20);
   double dd = (position - pt_central).lengthsq();
   float v = sqrt(dd) * 0.02;
   GdkRGBA c = fraction_to_colour(v);
#endif
   return c;
}



// not const because we sort in place the triangles of tri_con
void
molecule_class_info_t::draw_solid_density_surface(bool do_flat_shading) {

   // 20200312-PE empty this old function

   if (do_flat_shading)
      return; //

   if (draw_it_for_map) {
   }
}

void
molecule_class_info_t::display_solid_surface_triangles(const coot::density_contour_triangles_container_t &tc,
						       bool do_flat_shading) const {


   glBegin(GL_TRIANGLES);

   if (do_flat_shading) {
      for (unsigned int i=0; i<tc.point_indices.size(); i++) {

	 glNormal3f(tc.point_indices[i].normal_for_flat_shading.x(),
		    tc.point_indices[i].normal_for_flat_shading.y(),
		    tc.point_indices[i].normal_for_flat_shading.z());
	 glVertex3f(tc.points[tc.point_indices[i].pointID[0]].x(),
		    tc.points[tc.point_indices[i].pointID[0]].y(),
		    tc.points[tc.point_indices[i].pointID[0]].z());

	 glNormal3f(tc.point_indices[i].normal_for_flat_shading.x(),
		    tc.point_indices[i].normal_for_flat_shading.y(),
		    tc.point_indices[i].normal_for_flat_shading.z());
	 glVertex3f(tc.points[tc.point_indices[i].pointID[1]].x(),
		    tc.points[tc.point_indices[i].pointID[1]].y(),
		    tc.points[tc.point_indices[i].pointID[1]].z());

	 glNormal3f(tc.point_indices[i].normal_for_flat_shading.x(),
		    tc.point_indices[i].normal_for_flat_shading.y(),
		    tc.point_indices[i].normal_for_flat_shading.z());
	 glVertex3f(tc.points[tc.point_indices[i].pointID[2]].x(),
		    tc.points[tc.point_indices[i].pointID[2]].y(),
		    tc.points[tc.point_indices[i].pointID[2]].z());
      }

   } else {

      glShadeModel(GL_SMOOTH);
      bool opacity_experiment = false;

      /*
      coot::Cartesian rc(graphics_info_t::RotationCentre_x(),
			 graphics_info_t::RotationCentre_y(),
			 graphics_info_t::RotationCentre_z());

      float dist = 0.5 * graphics_info_t::zoom;
      GL_matrix glm;
      clipper::Coord_orth eye_dir(0,0,1);
      glm.from_quaternion(graphics_info_t::quat);
      clipper::Mat33<double> m = glm.to_clipper_mat();
      clipper::Coord_orth rot_dir(m * eye_dir);
      coot::Cartesian rot_dir_c(rot_dir.x(), rot_dir.y(), rot_dir.z());
      coot::Cartesian rot_dir_uv = rot_dir_c.unit();
      */

      for (unsigned int i=0; i<tc.point_indices.size(); i++) {

	 /*
	   if (opacity_experiment) {

	   // in fresnel mode, perhaps we need to change the specular too
	   // to make the glass more shiny? Something for the shader.

	   //
	   coot::Cartesian n =
	   tc.normals[tc.point_indices[i].pointID[0]] +
	   tc.normals[tc.point_indices[i].pointID[1]] +
	   tc.normals[tc.point_indices[i].pointID[2]];
	   coot::Cartesian n_uv = n.unit();

	   float cos_theta = coot::dot_product(n_uv, rot_dir_uv);
	   double opacity = pow(1.0 - pow(cos_theta, 6.0), 3);

	   GLfloat  mat_diffuse[]   = {float(map_colour[0][0]),
	   float(map_colour[0][1]),
	   float(map_colour[0][2]),
	   static_cast<GLfloat>(0.5 * opacity)};
	   GLfloat  mat_ambient[]   = {float(0.3*map_colour[0][0]),
	   float(0.3*map_colour[0][1]),
	   float(0.3*map_colour[0][2]),
	   static_cast<GLfloat>(opacity)};
	   glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,   mat_ambient);
	   glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mat_diffuse);
	   }
	 */

	 glNormal3f(tc.normals[tc.point_indices[i].pointID[0]].x(),
		    tc.normals[tc.point_indices[i].pointID[0]].y(),
		    tc.normals[tc.point_indices[i].pointID[0]].z());
	 glVertex3f(tc.points[tc.point_indices[i].pointID[0]].x(),
		    tc.points[tc.point_indices[i].pointID[0]].y(),
		    tc.points[tc.point_indices[i].pointID[0]].z());

	 glNormal3f(tc.normals[tc.point_indices[i].pointID[1]].x(),
		    tc.normals[tc.point_indices[i].pointID[1]].y(),
		    tc.normals[tc.point_indices[i].pointID[1]].z());
	 glVertex3f(tc.points[tc.point_indices[i].pointID[1]].x(),
		    tc.points[tc.point_indices[i].pointID[1]].y(),
		    tc.points[tc.point_indices[i].pointID[1]].z());

	 glNormal3f(tc.normals[tc.point_indices[i].pointID[2]].x(),
		    tc.normals[tc.point_indices[i].pointID[2]].y(),
		    tc.normals[tc.point_indices[i].pointID[2]].z());
	 glVertex3f(tc.points[tc.point_indices[i].pointID[2]].x(),
		    tc.points[tc.point_indices[i].pointID[2]].y(),
		    tc.points[tc.point_indices[i].pointID[2]].z());
      }
   }

   glEnd();
}

// is_neg is an optional arg
void
molecule_class_info_t::setup_density_surface_material(bool solid_mode, float opacity, bool is_neg) {

   if (solid_mode) {

      // normal solid

      GLfloat  ambientLight[] = { 0.3f, 0.3f, 0.3f, 1.0f };
      GLfloat  diffuseLight[] = { 0.4f, 0.4f, 0.4f, 1.0f };
      GLfloat specularLight[] = { 0.2f, 0.2f, 0.2f, 1.0f };

      // Assign created components to GL_LIGHT2
      glLightfv(GL_LIGHT2, GL_AMBIENT, ambientLight);
      glLightfv(GL_LIGHT2, GL_DIFFUSE, diffuseLight);
      glLightfv(GL_LIGHT2, GL_SPECULAR, specularLight);

      glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 100); // makes a difference.
      glDisable(GL_COLOR_MATERIAL);

      // narrowing from doubles to floats (there is no glMaterialdv).

      GLfloat  mat_specular[]  = {0.6f,  0.6f,  0.6f,  opacity}; // makes a difference
      GLfloat  mat_ambient[]   = {float(0.3*map_colour.red),
				                      float(0.3*map_colour.green),
				                      float(0.3*map_colour.blue), opacity};
      GLfloat  mat_diffuse[]   = {float(map_colour.red),
				                      float(map_colour.green),
				                      float(map_colour.blue), opacity};
      GLfloat  mat_shininess[] = {100}; // makes a difference

      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  mat_specular);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,   mat_ambient);
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mat_diffuse);

      if (is_neg) {
	      // override
         GLfloat  mat_specular[]  = {0.4,  0.4,  0.4,  opacity};
         GLfloat  mat_ambient[]   = {float(0.3*map_colour_negative_level.red),
                                     float(0.3*map_colour_negative_level.green),
				                         float(0.3*map_colour_negative_level.blue), opacity};
	      GLfloat  mat_diffuse[]   = {float(map_colour.red),
				                         float(map_colour.green),
				                         float(map_colour.blue), opacity};
	      GLfloat  mat_shininess[] = {100};

// 	 std::cout << " is_neg with map_colour: "
// 		   << map_colour[1][0] << " "
// 		   << map_colour[1][1] << " "
// 		   << map_colour[1][2] << " "
// 		   << std::endl;

         glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,  mat_specular);
         glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
         glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,   mat_ambient);
         glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,   mat_diffuse);

      }
   }
}


// modify v
void
molecule_class_info_t::dynamically_transform(coot::CartesianPairInfo v) {

   int s = v.size;
   for (int i=0; i<s; i++) {
      clipper::Coord_orth c1(v.data[i].getStart().x(),
			     v.data[i].getStart().y(),
			     v.data[i].getStart().z());
      clipper::Coord_orth c2(v.data[i].getFinish().x(),
			     v.data[i].getFinish().y(),
			     v.data[i].getFinish().z());
      clipper::Coord_orth ct1 = c1.transform(map_ghost_info.rtop);
      clipper::Coord_orth ct2 = c2.transform(map_ghost_info.rtop);
      v.data[i] = coot::CartesianPair(coot::Cartesian(ct1.x(), ct1.y(), ct1.z()),
				coot::Cartesian(ct2.x(), ct2.y(), ct2.z()));
   }

}

void
molecule_class_info_t::map_fill_from_mtz(const coot::mtz_to_map_info_t &mmi, const std::string &cwd, float sampling_rate) {

   map_fill_from_mtz(mmi.mtz_file_name, cwd, mmi.f_col, mmi.phi_col, mmi.w_col, mmi.use_weights, mmi.is_difference_map, sampling_rate);

}


//
void
molecule_class_info_t::map_fill_from_mtz(std::string mtz_file_name,
					 std::string cwd,
					 std::string f_col,
					 std::string phi_col,
					 std::string weight_col,
					 int use_weights,
					 int is_diff_map,
					 float sampling_rate) {

   short int use_reso_flag = 0;
   short int is_anomalous_flag = 0;
   map_fill_from_mtz_with_reso_limits(mtz_file_name,
				      cwd,
				      f_col,
				      phi_col,
				      weight_col,
				      use_weights,
				      is_anomalous_flag,
				      is_diff_map,
				      use_reso_flag, 0.0, 0.0, sampling_rate); // don't use these reso limits.

}


//
void
molecule_class_info_t::map_fill_from_mtz_with_reso_limits(std::string mtz_file_name,
							  std::string cwd,
							  std::string f_col,
							  std::string phi_col,
							  std::string weight_col,
							  int use_weights,
							  short int is_anomalous_flag,
							  int is_diff_map,
							  short int use_reso_limits,
							  float low_reso_limit,
							  float high_reso_limit,
							  float map_sampling_rate) {

   graphics_info_t g;

   // save for potential phase recombination in refmac later
   if (use_weights) {
      fourier_f_label = f_col;
      fourier_phi_label = phi_col;
      fourier_weight_label = weight_col; // magic label, we can go
					 // combining if this is not
					 // "";
//       std::cout << "DEBUG:: saving fourier_weight_label: " <<
// 	 fourier_weight_label << std::endl;
   }


   // std::cout << "DEBUG:: reso tinkering " << use_reso_limits << std::endl;
   clipper::Resolution user_resolution(high_reso_limit);
   clipper::Resolution fft_reso; // filled later

   //clipper::HKL_info myhkl;
   //clipper::MTZdataset mtzset;
   //clipper::MTZcrystal mtzxtl;

   long T0 = 0; // timer
   T0 = 0; // glutGet(GLUT_ELAPSED_TIME);

   clipper::CCP4MTZfile mtzin;
   mtzin.open_read( mtz_file_name );       // open new file
   clipper::HKL_data< clipper::datatypes::F_sigF<float> >  f_sigf_data;
   clipper::HKL_data< clipper::datatypes::Phi_fom<float> > phi_fom_data;
   clipper::HKL_data< clipper::datatypes::F_phi<float> >   fphidata;

   std::string mol_name = mtz_file_name + " ";
   mol_name += f_col;
   mol_name += " ";
   mol_name += phi_col;

   if (use_weights) {
      mol_name += " ";
      mol_name += weight_col;
   }

   if (use_reso_limits) {
      mol_name += " ";
      mol_name += g.float_to_string(low_reso_limit);
      mol_name += " ";
      mol_name += g.float_to_string(high_reso_limit);
   }

   initialize_map_things_on_read_molecule(mol_name, is_diff_map, is_anomalous_flag,
                                          g.swap_difference_map_colours);

   // If use weights, use both strings, else just use the first
   std::pair<std::string, std::string> p = make_import_datanames(f_col, phi_col, weight_col, use_weights);

   if (p.first.length() == 0) { // mechanism to signal an error
      std::cout << "ERROR:: fill_map.. - There was a column label error.\n";
   } else {

      if (use_weights) {
         // 	 std::cout << "DEBUG:: Importing f_sigf_data: " << p.first << std::endl;
         mtzin.import_hkl_data( f_sigf_data, p.first );
         // std::cout << "DEBUG:: Importing phi_fom_data: " << p.second << std::endl;
         mtzin.import_hkl_data(phi_fom_data, p.second);
         mtzin.close_read();
         fphidata.init( f_sigf_data.spacegroup(), f_sigf_data.cell(), f_sigf_data.hkl_sampling() );
         fphidata.compute(f_sigf_data, phi_fom_data,
                          clipper::datatypes::Compute_fphi_from_fsigf_phifom<float>());
       } else {
          // std::cout << "DEBUG:: Importing f_phi_data: " << p.first << std::endl;
          mtzin.import_hkl_data(fphidata, p.first);
          mtzin.close_read();
       }

      long T1 = 0; // glutGet(GLUT_ELAPSED_TIME);

      int n_reflections = fphidata.num_obs();
      std::cout << "INFO:: Number of observed reflections: " << n_reflections << "\n";
      if (n_reflections <= 0) {
         std::cout << "WARNING:: No reflections in mtz file!?" << std::endl;
      } else {
	 if (use_reso_limits) {
	    fft_reso = user_resolution;
	    filter_by_resolution(&fphidata, low_reso_limit, high_reso_limit);
	    data_resolution_ = high_reso_limit;
	 } else {
	    // fft_reso = myhkl.resolution();
	    // Kevin says do this instead:
	    //fft_reso = clipper::Resolution(1.0/sqrt(fphidata.invresolsq_range().max()));
	    fft_reso = fphidata.resolution();
	    data_resolution_ = 1.0/sqrt(fft_reso.invresolsq_limit());
	 }

	 if (is_anomalous_flag) {
	    fix_anomalous_phases(&fphidata);
	 }


	 cout << "INFO:: finding ASU unique map points with sampling rate "
   	      << map_sampling_rate	<< endl;
         clipper::Grid_sampling gs(fphidata.spacegroup(),
				   fphidata.cell(),
				   fft_reso,
				   map_sampling_rate);
	 cout << "INFO:: grid sampling..." << gs.format() << endl;
	 xmap.init( fphidata.spacegroup(), fphidata.cell(), gs); // 1.5 default
	 // 	 cout << "Grid..." << xmap.grid_sampling().format() << "\n";

	 long T2 = 0; // glutGet(GLUT_ELAPSED_TIME);
// 	 std::cout << "MTZ:: debug:: " << myhkl.spacegroup().symbol_hm() << " "
// 		   << myhkl.cell().descr().a() << " "
// 		   << myhkl.cell().descr().b() << " "
// 		   << myhkl.cell().descr().c() << " "
// 		   << clipper::Util::rad2d(myhkl.cell().descr().alpha()) << " "
// 		   << clipper::Util::rad2d(myhkl.cell().descr().beta ()) << " "
// 		   << clipper::Util::rad2d(myhkl.cell().descr().gamma()) << std::endl;
// 	 std::cout << "MTZ:: debug:: n_reflections: " << myhkl.num_reflections()
// 		   << std::endl;
// 	 int ncount = 0;
// 	 clipper::HKL_info::HKL_reference_index hri;
// 	 for (hri=fphidata.first(); !hri.last(); hri.next()) {
// 	    if (ncount < 500)
// 	       std::cout << " MTZ fphi: " << hri.hkl().h() << " "
// 			 << hri.hkl().k() << " " << hri.hkl().l() << " "
// 			 << fphidata[hri].f() << " "
// 			 << clipper::Util::rad2d(fphidata[hri].phi()) << std::endl;
// 	    ncount++;
// 	 }

	 // cout << "doing fft..." << endl;
	 xmap.fft_from( fphidata );                  // generate map
	 // cout << "done fft..." << endl;

	 long T3 = 0; // glutGet(GLUT_ELAPSED_TIME);
   #if 0
	 std::cout << "INFO:: " << float(T1-T0)/1000.0 << " seconds to read MTZ file\n";
	 std::cout << "INFO:: " << float(T2-T1)/1000.0 << " seconds to initialize map\n";
	 std::cout << "INFO:: " << float(T3-T2)/1000.0 << " seconds for FFT\n";
   #endif
   std::cout << "Fix Timers\n";
	 update_map_in_display_control_widget();

	 // Fill the class variables:
	 //   clipper::Map_stats stats(xmap);
	 //   map_mean_ = stats.mean();
	 //   map_sigma_ = stats.std_dev();

	 mean_and_variance<float> mv = map_density_distribution(xmap, 40, false, false);

	 save_mtz_file_name = mtz_file_name;
	 save_f_col = f_col;
	 save_phi_col = phi_col;
	 save_weight_col = weight_col;
	 save_use_weights = use_weights;
	 save_is_anomalous_map_flag = is_anomalous_flag;
	 save_is_diff_map_flag = is_diff_map;
	 save_high_reso_limit = high_reso_limit;
	 save_low_reso_limit = low_reso_limit;
	 save_use_reso_limits = use_reso_limits;

	 //
	 map_mean_  = mv.mean;
	 map_sigma_ = sqrt(mv.variance);
	 map_max_   = mv.max_density;
	 map_min_   = mv.min_density;

	 original_fphis_filled = 1;
	 original_fphis.init(fphidata.spacegroup(),fphidata.cell(),fphidata.hkl_sampling());
	 original_fphis = fphidata;

	 long T4 = 0; // glutGet(GLUT_ELAPSED_TIME);
	 std::cout << "INFO:: " << float(T4-T3)/1000.0 << " seconds for statistics\n";

	 std::cout << "      Map extents: ..... "
		   << xmap.grid_sampling().nu() << " "
		   << xmap.grid_sampling().nv() << " "
		   << xmap.grid_sampling().nw() << " " << std::endl;
	 std::cout << "      Map mean: ........ " << map_mean_ << std::endl;
	 std::cout << "      Map sigma: ....... " << map_sigma_ << std::endl;
	 std::cout << "      Map maximum: ..... " << map_max_ << std::endl;
	 std::cout << "      Map minimum: ..... " << map_min_ << std::endl;

	 set_initial_contour_level();

	 // update_map_colour_menu_manual(g.n_molecules, name_.c_str());
	 // update_map_scroll_wheel_menu_manual(g.n_molecules, name_.c_str());

	 update_map();
	 long T5 = 0; // glutGet(GLUT_ELAPSED_TIME);
	 //std::cout << "INFO:: " << float(T5-T4)/1000.0 << " seconds for contour map\n";
	 //std::cout << "INFO:: " << float(T5-T0)/1000.0 << " seconds in total\n";

	 // save state strings

	 std::string cwd = coot::util::current_working_dir();
	 std::string f1  = coot::util::intelligent_debackslash(mtz_file_name);
	 std::string f2  = coot::util::relativise_file_name(f1, cwd);
	 if (have_sensible_refmac_params) {
	    save_state_command_strings_.push_back("make-and-draw-map-with-refmac-params");
	    save_state_command_strings_.push_back(single_quote(f2));
	    save_state_command_strings_.push_back(single_quote(f_col));
	    save_state_command_strings_.push_back(single_quote(phi_col));
	    save_state_command_strings_.push_back(single_quote(weight_col));
	    save_state_command_strings_.push_back(g.int_to_string(use_weights));
	    save_state_command_strings_.push_back(g.int_to_string(is_diff_map));
	    save_state_command_strings_.push_back(g.int_to_string(1)); // have refmac params
	    save_state_command_strings_.push_back(single_quote(refmac_fobs_col));
	    save_state_command_strings_.push_back(single_quote(refmac_sigfobs_col));
	    save_state_command_strings_.push_back(single_quote(refmac_r_free_col));
	    save_state_command_strings_.push_back(g.int_to_string(refmac_r_free_flag_sensible));
	 } else {
	    if (save_use_reso_limits) {
	       save_state_command_strings_.push_back("make-and-draw-map-with-reso-with-refmac-params");
	       save_state_command_strings_.push_back(single_quote(f2));
	       save_state_command_strings_.push_back(single_quote(f_col));
	       save_state_command_strings_.push_back(single_quote(phi_col));
	       save_state_command_strings_.push_back(single_quote(weight_col));
	       save_state_command_strings_.push_back(g.int_to_string(use_weights));
	       save_state_command_strings_.push_back(g.int_to_string(is_diff_map));
	       save_state_command_strings_.push_back(g.int_to_string(0)); // have refmac params
	       save_state_command_strings_.push_back(single_quote(""));
	       save_state_command_strings_.push_back(single_quote(""));
	       save_state_command_strings_.push_back(single_quote(""));
	       save_state_command_strings_.push_back(g.int_to_string(0)); // sensible r-free
	       save_state_command_strings_.push_back(g.int_to_string(is_anomalous_flag));
	       save_state_command_strings_.push_back(g.int_to_string(save_use_reso_limits));
	       save_state_command_strings_.push_back(g.float_to_string( low_reso_limit));
	       save_state_command_strings_.push_back(g.float_to_string(high_reso_limit));
	    } else {
	       if (is_anomalous_flag) {
		  save_state_command_strings_.push_back("make-and-draw-map-with-reso-with-refmac-params");
		  save_state_command_strings_.push_back(single_quote(f2));
		  save_state_command_strings_.push_back(single_quote(f_col));
		  save_state_command_strings_.push_back(single_quote(phi_col));
		  save_state_command_strings_.push_back(single_quote(weight_col));
		  save_state_command_strings_.push_back(g.int_to_string(use_weights));
		  save_state_command_strings_.push_back(g.int_to_string(is_diff_map));
		  save_state_command_strings_.push_back(g.int_to_string(0)); // have refmac params
		  save_state_command_strings_.push_back(single_quote(""));
		  save_state_command_strings_.push_back(single_quote(""));
		  save_state_command_strings_.push_back(single_quote(""));
		  save_state_command_strings_.push_back(g.int_to_string(0)); // sensible r-free
		  save_state_command_strings_.push_back(g.int_to_string(is_anomalous_flag));
		  save_state_command_strings_.push_back(g.int_to_string(0)); // use reso limits
		  save_state_command_strings_.push_back(g.float_to_string(999.9));
		  save_state_command_strings_.push_back(g.float_to_string(1.2));
	       } else {
		  // bog standard.
		  save_state_command_strings_.push_back("make-and-draw-map");
		  save_state_command_strings_.push_back(single_quote(f2));
		  save_state_command_strings_.push_back(single_quote(f_col));
		  save_state_command_strings_.push_back(single_quote(phi_col));
		  save_state_command_strings_.push_back(single_quote(weight_col));
		  save_state_command_strings_.push_back(g.int_to_string(use_weights));
		  save_state_command_strings_.push_back(g.int_to_string(is_diff_map));
	       }
	    }
	 }
      }
   }
}



// return succes status, if mtz file is broken or empty, or
// non-existant, return 0.
//
bool
molecule_class_info_t::map_fill_from_cns_hkl(std::string cns_file_name,
					     std::string f_col,
					     int is_diff_map,
					     float map_sampling_rate)
{
   graphics_info_t g;

   try {
      long T0 = 0; // timer
      T0 = 0; // glutGet(GLUT_ELAPSED_TIME);

      clipper::CNS_HKLfile cnsin;
      cnsin.open_read( cns_file_name );       // open new file
      if (cnsin.cell().is_null() || cnsin.spacegroup().is_null()) {
	 std::cout << "WARNING:: Not an extended CNS file" << std::endl;
	 return 0;
      }
      clipper::HKL_sampling hklsam( cnsin.cell(), cnsin.resolution() );
      clipper::HKL_data< clipper::datatypes::F_phi<float> >
	 fphidata( cnsin.spacegroup(), cnsin.cell(), hklsam );
      cnsin.import_hkl_data( fphidata, f_col );
      cnsin.close_read();

      std::string mol_name = cns_file_name + " ";
      mol_name += f_col;

      original_fphis_filled = 1;
      original_fphis.init(fphidata.spacegroup(),fphidata.cell(),fphidata.hkl_sampling());
      original_fphis = fphidata;


      initialize_map_things_on_read_molecule(mol_name,
					     is_diff_map, false,
					     g.swap_difference_map_colours);
      long T1 = 0; // glutGet(GLUT_ELAPSED_TIME);

      int n_reflections = fphidata.num_obs();
      std::cout << "Number of OBSERVED reflections: " << n_reflections << "\n";
      if (n_reflections <= 0) {
	 std::cout << "WARNING:: No reflections in cns file!?" << std::endl;
	 return 0;
      }
      cout << "INFO:: finding ASU unique map points with sampling rate "
	   << map_sampling_rate	<< endl;
      clipper::Grid_sampling gs(fphidata.spacegroup(),
				fphidata.cell(),
				fphidata.resolution(),
				map_sampling_rate);
      cout << "INFO grid sampling..." << gs.format() << endl;
      xmap.init( fphidata.spacegroup(), fphidata.cell(), gs ); // 1.5 default
      // 	 cout << "Grid..." << xmap.grid_sampling().format() << "\n";

      long T2 = 0; // glutGet(GLUT_ELAPSED_TIME);

      // cout << "doing fft..." << endl;
      xmap.fft_from( fphidata );                  // generate map
      // cout << "done fft..." << endl;

      long T3 = 0; // glutGet(GLUT_ELAPSED_TIME);
      std::cout << "INFO:: " << float(T1-T0)/1000.0 << " seconds to read CNS file\n";
      std::cout << "INFO:: " << float(T2-T1)/1000.0 << " seconds to initialize map\n";
      std::cout << "INFO:: " << float(T3-T2)/1000.0 << " seconds for FFT\n";
      update_map_in_display_control_widget();

      // Fill the class variables:
      //   clipper::Map_stats stats(xmap);
      //   map_mean_ = stats.mean();
      //   map_sigma_ = stats.std_dev();

      mean_and_variance<float> mv = map_density_distribution(xmap, 40, false);

      save_mtz_file_name = cns_file_name;
      save_f_col = f_col;
      save_phi_col = "";
      save_weight_col = "";
      save_use_weights = 0;
      save_is_anomalous_map_flag = 0;
      save_is_diff_map_flag = is_diff_map;

      map_mean_  = mv.mean;
      map_sigma_ = sqrt(mv.variance);
      map_max_   = mv.max_density;
      map_min_   = mv.min_density;

      original_fphis.init(fphidata.spacegroup(),fphidata.cell(),fphidata.hkl_sampling());
      original_fphis = fphidata;

      long T4 = 0; // glutGet(GLUT_ELAPSED_TIME);
      std::cout << "INFO:: " << float(T4-T3)/1000.0 << " seconds for statistics\n";

      std::cout << "      Map mean: ........ " << map_mean_ << std::endl;
      std::cout << "      Map sigma: ....... " << map_sigma_ << std::endl;
      std::cout << "      Map maximum: ..... " << map_max_ << std::endl;
      std::cout << "      Map minimum: ..... " << map_min_ << std::endl;

      set_initial_contour_level();

      update_map();
      long T5 = 0; // glutGet(GLUT_ELAPSED_TIME);
      std::cout << "INFO:: " << float(T5-T4)/1000.0 << " seconds for contour map\n";
      std::cout << "INFO:: " << float(T5-T0)/1000.0 << " seconds in total\n";
      return 1;
   }
   catch (const clipper::Message_base &rte) {
      std::cout << "WARNING:: bad read of CNS HKL file " << cns_file_name << std::endl;
      return 0;
   }
}


void
molecule_class_info_t::set_refmac_save_state_commands(std::string mtz_file_name,
						      std::string f_col,
						      std::string phi_col,
						      std::string weight_col,
						      bool use_weights,
						      bool is_diff_map,
						      std::string refmac_fobs_col,
						      std::string refmac_sigfobs_col,
						      std::string refmac_r_free_col,
						      bool refmac_r_free_flag_sensible) {

   have_sensible_refmac_params = 1;
   save_state_command_strings_.clear();
   save_state_command_strings_.push_back("make-and-draw-map-with-refmac-params");
   save_state_command_strings_.push_back(single_quote(coot::util::intelligent_debackslash(mtz_file_name)));
   save_state_command_strings_.push_back(single_quote(f_col));
   save_state_command_strings_.push_back(single_quote(phi_col));
   save_state_command_strings_.push_back(single_quote(weight_col));
   save_state_command_strings_.push_back(coot::util::int_to_string(use_weights));
   save_state_command_strings_.push_back(coot::util::int_to_string(is_diff_map));
   save_state_command_strings_.push_back(coot::util::int_to_string(1)); // have refmac params
   save_state_command_strings_.push_back(single_quote(refmac_fobs_col));
   save_state_command_strings_.push_back(single_quote(refmac_sigfobs_col));
   save_state_command_strings_.push_back(single_quote(refmac_r_free_col));
   save_state_command_strings_.push_back(coot::util::int_to_string(refmac_r_free_flag_sensible));
}


std::vector<coot::atom_attribute_setting_help_t>
molecule_class_info_t::get_refmac_params() const {

   std::vector<coot::atom_attribute_setting_help_t> r;

   if (Have_sensible_refmac_params()) {
      r.push_back(coot::util::intelligent_debackslash(save_mtz_file_name));
      r.push_back(save_f_col);
      r.push_back(save_phi_col);
      r.push_back(save_weight_col);
      r.push_back(save_use_weights);
      r.push_back(save_is_diff_map_flag);
      r.push_back(1); // have refmac_params
      // r.push_back(refmac_mtz_filename); not sure if this should be given twice...
      r.push_back(refmac_fobs_col);
      r.push_back(refmac_sigfobs_col);
      r.push_back(refmac_r_free_col);
      r.push_back(refmac_r_free_flag_sensible);
   }
   return r;
}


void
molecule_class_info_t::fix_anomalous_phases(clipper::HKL_data< clipper::datatypes::F_phi<float> > *fphidata) const {

   for (clipper::HKL_info::HKL_reference_index hri = fphidata->first(); !hri.last(); hri.next()) {
      (*fphidata)[hri].shift_phase(-M_PI_2);
   }
}



void
molecule_class_info_t::save_previous_map_colour() {

   if (has_xmap() || has_nxmap())
      previous_map_colour = map_colour;
}


void
molecule_class_info_t::restore_previous_map_colour() {

   if (has_xmap() || has_nxmap())
	    map_colour = previous_map_colour;
   update_map();
}


void
molecule_class_info_t::set_initial_contour_level() {

   float level = 1.0;
   if (xmap_is_diff_map) {
      if (map_sigma_ > 0.05) {
	      level = nearest_step(map_mean_ +
			      graphics_info_t::default_sigma_level_for_fofc_map*map_sigma_, 0.01);
      } else {
	      level = 3.0*map_sigma_;
      }
   } else {
      if (map_sigma_ > 0.05) {
	      level = nearest_step(map_mean_ + graphics_info_t::default_sigma_level_for_map*map_sigma_, 0.01);
      } else {
	      level = graphics_info_t::default_sigma_level_for_map * map_sigma_;
      }
   }

   if (0)
      std::cout << "..... in set_initial_contour_level() xmap_is_diff_map is " << xmap_is_diff_map
		<< " and map_sigma_ is " << map_sigma_ << " and default sigma leve is "
		<< graphics_info_t::default_sigma_level_for_fofc_map << " and map_mean is "
		<< map_mean_ << std::endl;
   contour_level = level;
}


//
void
molecule_class_info_t::draw_skeleton(bool is_dark_background) {

   if (has_xmap()) {

      coot::CartesianPair pair;

      set_bond_colour(GREY_BOND);
      glLineWidth(2.0);

      if (greer_skeleton_draw_on == 1) {

	 //cout << "greer_skeleton_draw_on: "
	 //	   << greer_skel_box.bonds_[0].num_lines<< endl;

	 glBegin(GL_LINES);
	 for (int j=0; j<greer_skel_box.bonds_[0].num_lines; j++) {

            glVertex3f(greer_skel_box.bonds_[0].pair_list[j].positions.getStart().get_x(),
		       greer_skel_box.bonds_[0].pair_list[j].positions.getStart().get_y(),
		       greer_skel_box.bonds_[0].pair_list[j].positions.getStart().get_z());
            glVertex3f(greer_skel_box.bonds_[0].pair_list[j].positions.getFinish().get_x(),
		       greer_skel_box.bonds_[0].pair_list[j].positions.getFinish().get_y(),
		       greer_skel_box.bonds_[0].pair_list[j].positions.getFinish().get_z());
	 }
	 glEnd();
      }

      if (fc_skeleton_draw_on == 1) {

	 for (int l=0; l<fc_skel_box.num_colours; l++) {
 	    if (colour_skeleton_by_random) {
	       //  	       set_skeleton_bond_colour_random(l, colour_table);
	       set_skeleton_bond_colour(0.96);
 	    } else {
// 	       std::cout << "skel: " << l
// 			 << " of  " <<  fc_skel_box.num_colours <<  " "
// 			 << (float(l)/float(fc_skel_box.num_colours)+0.01)/1.011
// 			 << std::endl;
	       set_skeleton_bond_colour( (float(l)/float(fc_skel_box.num_colours)+0.01)/1.011 );
	    }


	    glBegin(GL_LINES);
	    for (int j=0; j<fc_skel_box.bonds_[l].num_lines; j++) {

	       glVertex3f(fc_skel_box.bonds_[l].pair_list[j].positions.getStart().get_x(),
			  fc_skel_box.bonds_[l].pair_list[j].positions.getStart().get_y(),
			  fc_skel_box.bonds_[l].pair_list[j].positions.getStart().get_z());
	       glVertex3f(fc_skel_box.bonds_[l].pair_list[j].positions.getFinish().get_x(),
			  fc_skel_box.bonds_[l].pair_list[j].positions.getFinish().get_y(),
			  fc_skel_box.bonds_[l].pair_list[j].positions.getFinish().get_z());
	    }
	    glEnd();
	 }
      }
   }
}

// Added rotate colour_map for EJD 5/5/2004.
void
molecule_class_info_t::set_skeleton_bond_colour(float f) {

   float rotation_size = float(imol_no) * 2.0*graphics_info_t::rotate_colour_map_on_read_pdb/360.0;
   while (rotation_size > 1.0) {
      rotation_size -= 1.0;
   }

   if (0) {
      std::vector<float> c(3);
      c[0] = 0.1+0.6*f*graphics_info_t::skeleton_colour[0];
      c[1] = 0.1+0.9*f*graphics_info_t::skeleton_colour[1];
      c[2] = 0.1+0.2*f*graphics_info_t::skeleton_colour[2];
      std::vector<float> rgb_new = rotate_rgb(c, rotation_size);
   }

   std::vector<float> rgb_new(3);
   for (int i=0; i<3; i++)
      rgb_new[i] = graphics_info_t::skeleton_colour[i];

   glColor3f(rgb_new[0], rgb_new[1], rgb_new[2]);
}



void
molecule_class_info_t::set_colour_skeleton_by_segment() { // use random colouring

   colour_skeleton_by_random = 1;
}

void
molecule_class_info_t::set_colour_skeleton_by_level() { // use random colouring

   colour_skeleton_by_random = 0;
}


//
void
molecule_class_info_t::draw_fc_skeleton() {

}

//
void
molecule_class_info_t::update_clipper_skeleton() {

   if (has_xmap()) {

      // Create map extents (the extents of the skeletonization)
      // from the current centre.

      if (xskel_is_filled == 1) {

	 // graphics_info_t g;

	 if (!xmap.is_null() && xmap_is_diff_map != 1) {
	    //
	    float skeleton_box_radius = graphics_info_t::skeleton_box_radius;

	    GraphicalSkel cowtan;

	    // fc_skel_box: class object type graphical_bonds_container
	    //
	    coot::Cartesian rc(graphics_info_t::RotationCentre_x(),
			       graphics_info_t::RotationCentre_y(),
			       graphics_info_t::RotationCentre_z());
	    fc_skel_box = cowtan.make_graphical_bonds(xmap,xskel_cowtan,
						      rc,
						      skeleton_box_radius,
						      graphics_info_t::skeleton_level);
	 }
      }
   }
}

void
molecule_class_info_t::unskeletonize_map() {

   fc_skeleton_draw_on = 0;
   xskel_is_filled = 0;
   clipper::Xmap<int> empty;
   xskel_cowtan = empty;
}

#include "coot-utils/slurp-map.hh"

// Return -1 on error
int
molecule_class_info_t::read_ccp4_map(std::string filename, int is_diff_map_flag,
				     const std::vector<std::string> &acceptable_extensions) {

   // For now, where we try to read in a map and we crash in importing
   // a file that is not a ccp4 map, lets do some checking: first that
   // the file exists and is not a directory, then that the file has
   // an extension of ".map" or ".ext".  If not, then complain and
   // return having done nothing.

   // stat filename
   struct stat s;
   int status = stat(filename.c_str(), &s);
   if (status != 0) {
      std::cout << "WARNING:: Error reading " << filename << std::endl;
      return -1;
   } else {
      if (!S_ISREG (s.st_mode)) {
	 if (S_ISDIR(s.st_mode)) {
	    std::cout << "WARNING:: " << filename << " is a directory." << endl;
	 } else {
	    std::cout << "WARNING:: " << filename << " not a regular file." << endl;
	 }
	 return -1;
      }
   }

   // was a regular file, let's check the extension:
   //
#ifdef WINDOWS_MINGW
   std::string::size_type islash = coot::util::intelligent_debackslash(filename).find_last_of("/");
#else
   std::string::size_type islash = filename.find_last_of("/");
#endif // MINGW
   std::string tstring;
   if (islash == std::string::npos) {
      // no slash found
      tstring = filename;
   } else {
      tstring = filename.substr(islash + 1);
   }

   bool good_extension_flag = 0;
   for (unsigned int iextension=0; iextension<acceptable_extensions.size(); iextension++) {
      std::string::size_type imap = tstring.rfind(acceptable_extensions[iextension]);
      if (imap != std::string::npos) {
	 good_extension_flag = 1;
	 break;
      }
   }

   // not really extension checking, just that it has it in the
   // filename:
   if (good_extension_flag == 0) {

      std::cout << "Filename for a CCP4 map must end in .map or .ext "
		<< "or some other approved extension - sorry\n";
      return -1;
      std::string ws = "The filename for a CCP4 map must\n";
      ws += "currently end in .map or .ext - sorry.\n\n";
      ws += "The map must be a CCP4 map or Badness Will Happen! :-)\n";
      GtkWidget *w = graphics_info_t::wrapped_nothing_bad_dialog(ws);
      gtk_widget_show(w);
   }

   // KDC: check map type
   enum MAP_FILE_TYPE { CCP4, CNS };
   MAP_FILE_TYPE map_file_type;
   {
     FILE* file = fopen( filename.c_str(), "r" );
     int c1, c2;
     c1 = c2 = 0;
     for ( int i = 0; i < 16; i++ ) {
       int c = getc(file);
       if ( c == EOF ) break;
       if ( c == 0 )                                 c1++;
       if ( std::isalpha(c) || std::isdigit(c) || std::isspace(c) ) c2++;
     }
     if ( c1 > c2 ) map_file_type = CCP4;
     else           map_file_type = CNS;
   }

   if (map_file_type == CCP4)
      std::cout << "INFO:: map file type was determined to be CCP4 type\n";
   if (map_file_type == CNS)
      std::cout << "INFO:: map file type was determined to be CNS type\n";

   bool bad_read = false; // so far
   bool em = false;

   if ( map_file_type == CCP4 ) {

      bool done = false;
      if (coot::util::is_basic_em_map_file(filename)) {
         em = true;
         // fill xmap
         // done = coot::util::slurp_fill_xmap_from_map_file(filename, &xmap);
      }

      if (! done) {
         std::cout << "INFO:: attempting to read CCP4 map: " << filename << std::endl;
         // clipper::CCP4MAPfile file;
         clipper_map_file_wrapper file;
         try {
            file.open_read(filename);

            em = set_is_em_map(file);

            bool use_xmap = true; // not an nxmap
            if (true) {

               clipper::Grid_sampling fgs = file.grid_sampling();

               std::cout << ".......................... grid sampling " << fgs.format() << std::endl;

               clipper::Cell fcell = file.cell();
               double vol = fcell.volume();
               if (vol < 1.0) {
                  std::cout << "WARNING:: non-sane unit cell volume " << vol << " - skip read"
                  << std::endl;
                  bad_read = true;
               } else {
                  try {
                     file.import_xmap(xmap);
                  }
                  catch (const clipper::Message_generic &exc) {
                     std::cout << "WARNING:: failed to read " << filename
                     << " Bad ASU (inconsistant gridding?)." << std::endl;
                     bad_read = true;
                  }
               }
            } else {

               // Should never happen.  Not yet.
               //
               std::cout << "=================== EM Map NXmap =================== " << std::endl;
               file.import_nxmap(nxmap);
               std::cout << "INFO:: created NX Map with grid " << nxmap.grid().format() << std::endl;
            }
         } catch (const clipper::Message_base &exc) {
            std::cout << "WARNING:: failed to open " << filename << std::endl;
            bad_read = true;
         }

         std::pair<bool, coot::Cartesian> new_centre(false, coot::Cartesian(0,0,0)); // used only for first EM map

         if (em) {

            // If this was the first map, recentre to the middle of the cell
            //
            if (imol_no == 0) {
               clipper::Cell c = file.cell();
               coot::Cartesian m(0.5*c.descr().a(),
               0.5*c.descr().b(),
               0.5*c.descr().c());
               new_centre.first = true;
               new_centre.second = m;
               std::cout << "INFO:: map appears to be EM map."<< std::endl;
            }
            std::cout << "INFO:: closing CCP4 map: " << filename << std::endl;
            file.close_read();

            if (new_centre.first) {
               graphics_info_t g;
               g.setRotationCentre(new_centre.second);
            }
         }
      }
   } else {
     std::cout << "INFO:: attempting to read CNS map: " << filename << std::endl;
     clipper::CNSMAPfile file;
     file.open_read(filename);
     try {
       file.import_xmap( xmap );
     }
     catch (const clipper::Message_base &exc) {
       std::cout << "WARNING:: failed to read " << filename << std::endl;
       bad_read = true;
     }
     file.close_read();
   }

   if (! bad_read) {

      bool is_anomalous_flag = false;
      initialize_map_things_on_read_molecule(filename, is_diff_map_flag, is_anomalous_flag,
					     graphics_info_t::swap_difference_map_colours);

      mean_and_variance<float> mv = map_density_distribution(xmap, 40, true, true);

      float mean = mv.mean;
      float var = mv.variance;

      xmap_is_diff_map = is_diff_map_flag; // but it may be...
      // fill class variables
      map_mean_ = mv.mean;
      map_sigma_ = sqrt(mv.variance);
      map_max_   = mv.max_density;
      map_min_   = mv.min_density;

      update_map_in_display_control_widget();
      contour_level    = nearest_step(mean + 1.5*sqrt(var), 0.05);

      bool em = is_EM_map();

      if (em) {
	 // make better defaults
	 contour_level = mean + nearest_step(mean + 5.0*sqrt(var), 0.2);
	 contour_sigma_step = 0.4;
      } else {
	 // "how it used to be" logic.  contour_level is set above and reset here
	 // Hmm.
	 set_initial_contour_level();
      }

      std::cout << "      Map extents: ..... "
		<< xmap.grid_sampling().nu() << " "
		<< xmap.grid_sampling().nv() << " "
		<< xmap.grid_sampling().nw() << " " << std::endl;
      std::cout << "      Map mean: ........ " << map_mean_ << std::endl;
      std::cout << "      Map rmsd: ........ " << map_sigma_ << std::endl;
      std::cout << "      Map maximum: ..... " << map_max_ << std::endl;
      std::cout << "      Map minimum: ..... " << map_min_ << std::endl;

      // save state strings
      save_state_command_strings_.push_back("handle-read-ccp4-map");
      save_state_command_strings_.push_back(single_quote(coot::util::intelligent_debackslash(filename)));
      save_state_command_strings_.push_back(graphics_info_t::int_to_string(is_diff_map_flag));

      update_map();
   }

   int stat = imol_no;
   if (bad_read)
      stat = -1;
   return stat;
}

// is the CCP4 map a EM map? (this is so that we can fill the
// NXmap, not the xmap)
//
bool
molecule_class_info_t::set_is_em_map(const clipper_map_file_wrapper &file) {

   // Even if mapdump says that the spacegroup is 0, file.spacegroup()
   // will be "P1".  So this returns true for maps with spacegroup 0
   // (and 90 degrees)

   if (file.spacegroup().num_symops() == 1) { // P1
      if (((file.cell().descr().alpha() - M_PI/2) <  0.0001) &&
	  ((file.cell().descr().alpha() - M_PI/2) > -0.0001) &&
	  ((file.cell().descr().beta()  - M_PI/2) > -0.0001) &&
	  ((file.cell().descr().beta()  - M_PI/2) <  0.0001) &&
	  ((file.cell().descr().gamma() - M_PI/2) > -0.0001) &&
	  ((file.cell().descr().gamma() - M_PI/2) <  0.0001)) {
	 if (file.starts_at_zero()) {
	    is_em_map_cached_flag = 1; // yes
	 } else {
	    is_em_map_cached_flag = 0;
	 }
      } else {
	 is_em_map_cached_flag = 0;
      }
   } else {
      is_em_map_cached_flag = 0;
   }
   return false; // not a useful return values, because flag can have 3 values
}

bool
molecule_class_info_t::is_EM_map() const {

   bool ret_is_em = false;

   if (has_xmap()) {
      if (is_em_map_cached_flag == 1) { // -1 means unset
	 ret_is_em = true;
      }
   }
   return ret_is_em;
}

short int
molecule_class_info_t::is_em_map_cached_state() {

   if (is_em_map_cached_flag == -1) {

      if (has_xmap()) { // FIXME - need to test for NXmap too.
	 bool is_em = is_EM_map();
	 is_em_map_cached_flag = is_em;
      }
   }
   return is_em_map_cached_flag;
}


void
molecule_class_info_t::install_new_map(const clipper::Xmap<float> &map_in, std::string name_in, bool is_em_map_flag_in) {

   xmap = map_in;
   if (is_em_map_flag_in)
      is_em_map_cached_flag = 1;
   // the map name is filled by using set_name(std::string)
   // sets name_ to name_in:
   initialize_map_things_on_read_molecule(name_in, false, false, false); // not a diff_map

   mean_and_variance<float> mv = map_density_distribution(xmap, 40, true);

   float mean = mv.mean;
   float var = mv.variance;
   contour_level  = nearest_step(mean + 1.5*sqrt(var), 0.05);
   update_map_in_display_control_widget();

   // fill class variables
   map_mean_ = mv.mean;
   map_sigma_ = sqrt(mv.variance);

   update_map();
}

void
molecule_class_info_t::set_mean_and_sigma() {

   mean_and_variance<float> mv = map_density_distribution(xmap, 40, true);
   map_mean_ = mv.mean;
   map_sigma_ = sqrt(mv.variance);

}

void
molecule_class_info_t::set_name(std::string name) {
   name_ = name;
   update_mol_in_display_control_widget();

}


int
molecule_class_info_t::make_map_from_phs(std::string pdb_filename,
                                         std::string phs_filename) {

   int iret = -1; // default error return status
   //
   std::cout << "INFO:: Make a map from " << phs_filename << " using "
	     << pdb_filename << " for the cell and symmetry information " << std::endl;

   atom_selection_container_t SelAtom = get_atom_selection(pdb_filename, true, true);

   if (SelAtom.read_success == 1) { // success
      try {
	 std::pair<clipper::Cell,clipper::Spacegroup> xtal =
	    coot::util::get_cell_symm( SelAtom.mol );
	 iret = make_map_from_phs(xtal.second, xtal.first, phs_filename);
      } catch (const std::runtime_error &except) {
	 cout << "!! get_cell_symm() fails in make_map_from_phs"
	      << endl;
      }
   }
   return iret;
}


int
molecule_class_info_t::make_map_from_phs_using_reso(std::string phs_filename,
						    const clipper::Spacegroup &sg,
						    const clipper::Cell &cell,
						    float reso_limit_low,
						    float reso_limit_high,
						    float map_sampling_rate) {

   clipper::PHSfile phs;

   phs.open_read(phs_filename);

   // std::cout << "creating resolution" << std::endl;
   clipper::Resolution resolution(reso_limit_high);

   clipper::HKL_info mydata(sg, cell, resolution);
   clipper::HKL_data<clipper::datatypes::F_sigF<float>  >  myfsig(mydata);
   clipper::HKL_data<clipper::datatypes::Phi_fom<float> >  myphwt(mydata);
   clipper::HKL_data<clipper::datatypes::F_phi<float>   >  fphidata(mydata);

   std::cout << "importing info" << std::endl;
   phs.import_hkl_info(mydata);
   std::cout << "importing data" << std::endl;
   phs.import_hkl_data(myfsig);
   phs.import_hkl_data(myphwt);

   phs.close_read();

   std::cout << "PHS file: Number of reflections: " << mydata.num_reflections() << "\n";

   fphidata.update();

   fphidata.compute(myfsig, myphwt,
 		    clipper::datatypes::Compute_fphi_from_fsigf_phifom<float>());

//    for (int i=0; i<10; i++) {
//       std::cout << "checking phi weight: " << i << " " << myphwt[i].phi() << "  " << myphwt[i].fom()
// 		<< std::endl;
//        std::cout << "checking f    sigf: " << i << " " << myfsig[i].f() << "   "
// 		 << myfsig[i].sigf() << std::endl;
//        std::cout << "checking missing: " << i << " " << myfsig[i].missing() << " "
// 		 << myphwt[i].missing() << " " << fphidata[i].missing() << std::endl;
//        // << " " << fphidata[i].phi() <<
//    }

  std::string mol_name = phs_filename;

  initialize_map_things_on_read_molecule(mol_name, false, false, false); // not diff map

  std::cout << "initializing map...";
  xmap.init(mydata.spacegroup(),
		    mydata.cell(),
		    clipper::Grid_sampling(mydata.spacegroup(),
					   mydata.cell(),
					   mydata.resolution(),
					   map_sampling_rate));
 std:cout << "done."<< std::endl;

//   cout << "Map Grid (from phs file)..."
//        << xmap.grid_sampling().format()
//        << endl;

  std::cout << "doing fft..." ;
  xmap.fft_from(fphidata);                  // generate map
  std::cout << "done." << std::endl;

  mean_and_variance<float> mv = map_density_distribution(xmap, 40, false);

  cout << "Mean and sigma of map from PHS file: " << mv.mean
       << " and " << sqrt(mv.variance) << endl;

  // fill class variables
  map_mean_ = mv.mean;
  map_sigma_ = sqrt(mv.variance);

  original_fphis_filled = 1;
  original_fphis.init(fphidata.spacegroup(),fphidata.cell(),fphidata.hkl_sampling());
  original_fphis = fphidata;


  xmap_is_diff_map = 0;
  update_map_in_display_control_widget();
  contour_level = nearest_step(mv.mean + 1.5*sqrt(mv.variance), 0.05);

  std::cout << "updating map..." << std::endl;
  update_map();
  std::cout << "done updating map..." << std::endl;

  // as for 'normal' maps
  std::string cwd = coot::util::current_working_dir();
  std::string f1  = coot::util::intelligent_debackslash(phs_filename);
  std::string f2  = coot::util::relativise_file_name(f1, cwd);
  // how do we restore this map?
  save_state_command_strings_.push_back("read-phs-and-make-map-using-cell-symm");
  save_state_command_strings_.push_back(single_quote(f2));
  save_state_command_strings_.push_back(single_quote(sg.symbol_hm()));
  save_state_command_strings_.push_back(coot::util::float_to_string(cell.descr().a()));
  save_state_command_strings_.push_back(coot::util::float_to_string(cell.descr().b()));
  save_state_command_strings_.push_back(coot::util::float_to_string(cell.descr().c()));
  save_state_command_strings_.push_back(coot::util::float_to_string(clipper::Util::rad2d(cell.descr().alpha())));
  save_state_command_strings_.push_back(coot::util::float_to_string(clipper::Util::rad2d(cell.descr().beta())));
  save_state_command_strings_.push_back(coot::util::float_to_string(clipper::Util::rad2d(cell.descr().gamma())));

  return imol_no;
}



// and the molecule number imol_coords where the coordinates are.
int
molecule_class_info_t::make_map_from_cif(int imol_no_in,
					 std::string cif_file_name,
					 int imol_coords) {

   graphics_info_t g;
   int r = -1;
   if (g.is_valid_model_molecule(imol_coords))
      r =  make_map_from_cif(imol_no_in, cif_file_name,
			     g.molecules[imol_coords].atom_sel);
   else
      std::cout << "WARNING:: " << imol_coords << " is not a valid model molecule" << std::endl;
   return r;
}

// and the molecule number imol_coords where the coordinates are.
int
molecule_class_info_t::make_map_from_cif_2fofc(int imol_no_in,
					       std::string cif_file_name,
					       int imol_coords) {

   graphics_info_t g;
   int r = -1;
   if (g.is_valid_model_molecule(imol_coords))
      r =  make_map_from_cif(imol_no_in, cif_file_name,
			     g.molecules[imol_coords].atom_sel);
   else
      std::cout << "WARNING:: " << imol_coords << " is not a valid model molecule" << std::endl;
   return r;
}

// and the molecule number imol_coords where the coordinates are.
int
molecule_class_info_t::make_map_from_cif_fofc(int imol_no_in,
					      std::string cif_file_name,
					      int imol_coords) {

   graphics_info_t g;
   int r = -1;
   if (g.is_valid_model_molecule(imol_coords))
      r = make_map_from_cif_generic(imol_no_in,
				    cif_file_name,
				    g.molecules[imol_coords].atom_sel,
				    2);  // 2 -> is Fo-Fc map
   else
      std::cout << "WARNING:: " << imol_coords << " is not a valid model molecule" << std::endl;
   return r;
}


int
molecule_class_info_t::make_map_from_cif(int imol_no_in,
					 std::string cif_file_name,
					 atom_selection_container_t SelAtom) {

   // 0 is not is_2fofc_type map (is sigmaa)
   return make_map_from_cif_generic(imol_no_in, cif_file_name, SelAtom, 0);

}

int
molecule_class_info_t::make_map_from_cif_2fofc(int imol_no_in,
					       std::string cif_file_name,
					       atom_selection_container_t SelAtom) {

   // 1 is is_2fofc_type map (not sigmaa)
   return make_map_from_cif_generic(imol_no_in, cif_file_name, SelAtom, 1);

}


int
molecule_class_info_t::make_map_from_cif_generic(int imol_in,
						 std::string cif_file_name,
						 atom_selection_container_t SelAtom,
						 short int is_2fofc_type) {

   clipper::HKL_data< clipper::datatypes::F_sigF<float> > myfsigf;
   clipper::CIFfile cif;
   cif.open_read ( cif_file_name );
   cif.import_hkl_data( myfsigf );
   cif.close_read();

   clipper::Spacegroup sg = myfsigf.spacegroup();
   if (! sg.is_null()) {
      std::cout << "DEBUG in make_map_from_cif_generic imol_in " << imol_in << std::endl;
      return calculate_sfs_and_make_map(imol_in, cif_file_name, myfsigf,
					SelAtom, is_2fofc_type);
   } else {
      std::cout << "ERROR:: null space group in make_map_from_cif_generic() " << std::endl;
      return -1;
   }
}


// fill original_fphis
void
molecule_class_info_t::save_original_fphis_from_map() {

   // clipper::HKL_data< clipper::datatypes::F_phi<float> > original_fphis;



}


int
molecule_class_info_t::calculate_sfs_and_make_map(int imol_no_in,
						  const std::string &mol_name,
						  const clipper::HKL_data< clipper::datatypes::F_sigF<float> > &myfsigf,
						  atom_selection_container_t SelAtom,
						  short int is_2fofc_type) {

   initialize_map_things_on_read_molecule(mol_name, false, false, false); // not diff map

   std::cout << "calculating structure factors..." << std::endl;

   // Fix up fphidata to contain the calculated structure factors

   // Calculated structure factors go here:
   const clipper::HKL_info& hkls = myfsigf.hkl_info();
   clipper::Spacegroup sg = myfsigf.spacegroup();
   if (sg.is_null()) {
      std::cout << "ERROR:: spacegroup from cif data is null" << std::endl;
      return -1;
   }


   clipper::HKL_data< clipper::datatypes::F_phi<float> > fphidata(sg, myfsigf.cell(),myfsigf.hkl_sampling());
   // map coefficients ((combined Fo and scaled Fc) and calc phi) go here:

   clipper::HKL_data< clipper::datatypes::F_phi<float> > map_fphidata(myfsigf.spacegroup(),myfsigf.cell(),myfsigf.hkl_sampling());

   // get a list of all the atoms
   clipper::MMDBAtom_list atoms(SelAtom.atom_selection, SelAtom.n_selected_atoms);

   std::cout << "isotropic fft of " << SelAtom.n_selected_atoms
	     << " atoms..." << std::endl;
   clipper::SFcalc_iso_fft<float>(fphidata, atoms);
   std::cout << "done iso fft..." << std::endl;

   // debug:: examine fphidata and myfsigf:
   std::cout << "INFO:: myfsigf  has " <<  myfsigf.data_size() << " data" << std::endl;
   std::cout << "INFO:: fphidata has " << fphidata.data_size() << " data" << std::endl;

   if (0) { // debug
      float sum_fo = 0;
      float sum_fc = 0;
      int n_fo = 0;
      int n_fc = 0;
      for (clipper::HKL_info::HKL_reference_index ih=myfsigf.first();
	   !ih.last(); ih.next()) {
	 if (!myfsigf[ih].missing()) {
	    n_fo++;
	    n_fc++;
	    sum_fo += myfsigf[ih].f();
	    sum_fc = fphidata[ih].f();
	 }
      }

      std::cout << "DEBUG:: fo: sum average: " << sum_fo << " " << sum_fo/float(n_fo)
		<< std::endl;
      std::cout << "DEBUG:: fc: sum average: " << sum_fc << " " << sum_fc/float(n_fc)
		<< std::endl;
      for (clipper::HKL_info::HKL_reference_index ih=myfsigf.first();
	   !ih.last(); ih.next())
	      std::cout << "DEBUG::  myfsigf " <<  " " <<  myfsigf[ih].f() << " "
		   << myfsigf[ih].sigf() << " " << myfsigf[ih].missing() << std::endl;
      for (int i=0; i<10; i++)
	 std::cout << "DEBUG:: fphidata " << i << " " << fphidata[i].f()
		   << " " << fphidata[i].phi() << std::endl;
   }

   int nprm = 10;
   std::vector<clipper::ftype> params_init( nprm, 1.0 );
   // clipper::BasisFn_spline basis_f1f2( mydata, nprm, 2.0 );
   //  target_f1f2( fc, fo );
   //clipper::TargetFn_scaleF1F2<clipper::datatypes::F_phi<float>,
   //                          clipper::datatypes::F_sigF<float> >
   // target_f1f2( fphidata, myfsigf );
   //clipper::ResolutionFn fscale( mydata, basis_f1f2,
   //                              target_f1f2, params_init );

   float r_top = 0.0, r_bot = 0.0;
   float sum_fo = 0.0, sum_fc = 0.0, sum_scale = 0.0;
   int n_data = 0;

   if (is_2fofc_type == molecule_map_type::TYPE_2FO_FC ||
       is_2fofc_type == molecule_map_type::TYPE_FO_FC) {

      if (is_2fofc_type == molecule_map_type::TYPE_2FO_FC)
	 std::cout << "INFO:: calculating 2fofc map..." << std::endl;
      if (is_2fofc_type == molecule_map_type::TYPE_FO_FC)
	 std::cout << "INFO:: calculating fofc map..." << std::endl;

      clipper::BasisFn_spline basis_f1f2( hkls, nprm, 2.0 );
      //  target_f1f2( fc, fo );
      clipper::TargetFn_scaleF1F2<clipper::datatypes::F_phi<float>,
	 clipper::datatypes::F_sigF<float> >
	 target_f1f2( fphidata, myfsigf );
      clipper::ResolutionFn fscale( hkls, basis_f1f2, target_f1f2, params_init );

      float multiplier = 2.0;
      if (is_2fofc_type == molecule_map_type::TYPE_FO_FC)
	 multiplier = 1.0;

      for ( clipper::HKL_info::HKL_reference_index ih=myfsigf.first();
	    !ih.last(); ih.next() ) {
	 map_fphidata[ih].phi() = fphidata[ih].phi();
	 if (!myfsigf[ih].missing()) {
	    map_fphidata[ih].f() = multiplier*myfsigf[ih].f() -
	       fphidata[ih].f()*sqrt(fscale.f(ih));
	    float top_tmp = fabs(myfsigf[ih].f() - fphidata[ih].f()*sqrt(fscale.f(ih)));
	    if (0) { // debug
	       std::cout << "debug:: fobs: " << myfsigf[ih].f() << " fcalc: "
			 << fphidata[ih].f() << " scale: " << fscale.f(ih)
			 << std::endl;
	    }

	    r_top += top_tmp;
	    r_bot += fabs(myfsigf[ih].f());
// 	    std::cout << "debug:: adding to top: " << top_tmp << " bot: "
// 		      << fabs(myfsigf[ih].f()) << std::endl;
	    sum_fo += myfsigf[ih].f();
	    sum_fc += fphidata[ih].f();
	    sum_scale += sqrt(fscale.f(ih));
	    n_data++;
	 } else {
	    map_fphidata[ih].f() = 0.0;
	 }
      }

   } else { // not 2fofc-style, i.e. is sigmaa style

      if (is_2fofc_type == molecule_map_type::TYPE_SIGMAA) {

	 std::cout << "sigmaa and scaling..." << std::endl;

	 // need an mmdb
	 mmdb::Manager *mmdb = SelAtom.mol;

	 // get a list of all the atoms
	 mmdb::PAtom *psel;
	 int hndl, nsel;
	 hndl = mmdb->NewSelection();
	 mmdb->SelectAtoms( hndl, 0, 0, mmdb::SKEY_NEW );
	 mmdb->GetSelIndex( hndl, psel, nsel );
	 clipper::MMDBAtom_list atoms( psel, nsel );
	 mmdb->DeleteSelection( hndl );

	 // calculate structure factors
	 const clipper::HKL_data<clipper::datatypes::F_sigF<float> >& fo = myfsigf;
	 clipper::HKL_data<clipper::datatypes::F_phi<float> > fc( hkls );
	 clipper::SFcalc_obs_bulk<float> sfcb;
	 sfcb( fc, fo, atoms );

	 // do anisotropic scaling
	 clipper::SFscale_aniso<float> sfscl;
	 // sfscl( fo, fc );  // scale Fobs
	 sfscl( fc, fo );  // scale Fcal

	 // now do sigmaa calc
	 clipper::HKL_data<clipper::datatypes::F_phi<float> >   fb( hkls ), fd( hkls );
	 clipper::HKL_data<clipper::datatypes::Phi_fom<float> > phiw( hkls );
	 clipper::HKL_data<clipper::datatypes::Flag>    flag( hkls );
	 typedef clipper::HKL_data_base::HKL_reference_index HRI;
	 // If no free flag is available, then use all reflections..
	 for (HRI ih = flag.first(); !ih.last(); ih.next() )
	    flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
	 /* This code uses free reflections only for sigmaa and scaling...
	 for (HRI ih = flag.first(); !ih.last(); ih.next() )
	    if ( !fo[ih].missing() && (free[ih].missing()||free[ih].flag()==0) )
	       flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
	    else
	       flag[ih].flag() = clipper::SFweight_spline<float>::NONE;
	 */

	 // do sigmaa calc

	 // Bug! 20071218 Kevin fixes it.
// 	 int n_refln = mydata.num_reflections();
	 int n_refln = 1000;
	 int n_param = 20;
	 clipper::SFweight_spline<float> sfw( n_refln, n_param );
	 sfw( fb, fd, phiw, fo, fc, flag );

	 // OK, so now fb and fd contain F_phis, one for "best"
	 // sigmaa, one for difference map.  Let's just use the "best"
	 // map for now.
	 map_fphidata = fb;

      }
   } // is 2fofc else sigmaa style check

   // std::cout << "DEBUG:: rdiffsum/rsum: " << r_top << "/" << r_bot << std::endl;
   if (is_2fofc_type != molecule_map_type::TYPE_SIGMAA) {
      if (r_bot>0.0) {
	 std::cout << "Isotropic R-factor: " << 100.0*r_top/r_bot << "%"
		   << " for " << n_data  << " reflections" <<  std::endl;
	 std::cout << "DEBUG:: sums: fo: " << sum_fo/float(n_data) << " fc: "
		   << sum_fc/n_data << " scale: " << sum_scale/n_data << " with "
		   << n_data << " data" << std::endl;
      } else {
	 std::cout << "Problem with structure factors, no structure factor sum!?"
		   << std::endl;
      }
   }
   std::cout << "Initializing map...";
   xmap.init(map_fphidata.spacegroup(), map_fphidata.cell(),
		     clipper::Grid_sampling(map_fphidata.spacegroup(),
					    map_fphidata.cell(),
					    map_fphidata.resolution()));
   cout << "done."<< endl;
   cout << "doing fft..." ;
   xmap.fft_from( map_fphidata ); // generate map
   cout << "done." << endl;

   mean_and_variance<float> mv = map_density_distribution(xmap, 40, false);

   cout << "Mean and sigma of map " << mol_name << " " << mv.mean
	<< " and " << sqrt(mv.variance) << endl;

   // fill class variables
   map_mean_ = mv.mean;
   map_sigma_ = sqrt(mv.variance);
   map_max_   = mv.max_density;
   map_min_   = mv.min_density;

   original_fphis.init(map_fphidata.spacegroup(),map_fphidata.cell(),map_fphidata.hkl_sampling());
   original_fphis = map_fphidata;

   xmap_is_diff_map = 0;
   update_map_in_display_control_widget();

   std::cout << "      Map mean: ........ " << map_mean_ << std::endl;
   std::cout << "      Map sigma: ....... " << map_sigma_ << std::endl;
   std::cout << "      Map maximum: ..... " << map_max_ << std::endl;
   std::cout << "      Map minimum: ..... " << map_min_ << std::endl;

  set_initial_contour_level();

   int imol = imol_no_in;
   update_map();
   return imol;
}

// This needs to be rationalized with the version that *does* pass the
// coordinates.
//
// This is the version that gets called when we use a file selector to
// get the file (i.e. it doesn't specify the coordinates (molecule))
// because there are calculated structure factors in the file.
//
// We make a Fc alpha-c map.  Which is not usually what we want.
//
int
molecule_class_info_t::make_map_from_cif(int imol_no_in,
					 std::string cif_file_name) {
   return make_map_from_cif_sigmaa(imol_no_in,
				   cif_file_name, molecule_map_type::TYPE_SIGMAA);
}

int
molecule_class_info_t::make_map_from_cif_diff_sigmaa(int imol_no_in,
						     std::string cif_file_name) {
   return make_map_from_cif_sigmaa(imol_no_in,
				   cif_file_name,
				   molecule_map_type::TYPE_DIFF_SIGMAA);
}

// SigmaA map type, either molecule_map_type::TYPE_SIGMAA or TYPE_DIFF_SIGMAA.
//
int
molecule_class_info_t::make_map_from_cif_sigmaa(int imol_no_in,
						std::string cif_file_name,
						int sigmaa_map_type) {

   imol_no = imol_no_in;
   clipper::HKL_info mydata;
   clipper::CIFfile cif;


   try {
      cif.open_read (cif_file_name);
      cif.import_hkl_info(mydata); // set spacegroup, cell and get hkl list.
      clipper::HKL_data< clipper::datatypes::F_sigF<float> > myfsigf(mydata); // Fobs
      clipper::HKL_data< clipper::datatypes::F_phi<float> > fc(mydata); // FC PHIC

      cif.import_hkl_data(myfsigf);
      cif.import_hkl_data(fc);

      cif.close_read();

      // std::cout << "DEBUG:: make_map_from_cif_sigmaa" << std::endl;
      std::cout << "Read " << mydata.num_reflections() << " from CIF file (sigmaa)."
		<< std::endl;

      if (mydata.num_reflections() == 0) {
	 return -1;
      } else {

	 // Are all the calculated sfs missing/zero?
	 //
	 int non_zero = 0;
	 for(int i=0; i< mydata.num_reflections(); i++) {
	    if (! fc[i].missing()) {
	       if (fc[i].f() > 0.0) {
		  non_zero = 1;
		  break;
	       }
	    }
	 }

	 if (non_zero == 0) {
	    std::cout << "WARNING:: Ooops - all the structure factor amplitudes "
		      << " appear to be zero - or missing.  " << std::endl;
	    std::cout << "WARNING:: Are you sure this file (" << cif_file_name
		      << ") contains calculated structure factors?" << std::endl;
	    std::cout << "WARNING:: No map calculated." << std::endl;
	    std::cout << "INFO:: if you want to calculate structure factors from a"
		      << " set of coordinates,  consider the function read_cif_data()"
		      << std::endl;
	 } else {

	    std::string mol_name = cif_file_name;
	    if (sigmaa_map_type == molecule_map_type::TYPE_SIGMAA)
	       mol_name += " SigmaA";
	    if (sigmaa_map_type == molecule_map_type::TYPE_DIFF_SIGMAA)
	       mol_name += " Difference SigmaA";

	    // new sigmaA code... needs to be updated to new Kevin
	    // code... but that is slightly tricky because here we have
	    // sfs, whereas KC code calculates them.

	    std::cout << "sigmaa and scaling..." << std::endl;

	    clipper::HKL_data< clipper::datatypes::F_phi<float> > map_fphidata(mydata);
	    clipper::HKL_data<clipper::datatypes::Phi_fom<float> > phifom(mydata);

	    clipper::Cell cxtl = myfsigf.hkl_info().cell();
	    // "Aliases" to fix Kevin's sigmaA code into mine
	    const clipper::HKL_data<clipper::datatypes::F_sigF<float> >& fo = myfsigf;
	    const clipper::HKL_info& hkls = mydata;

	    // now do sigmaa calc
	    clipper::HKL_data<clipper::datatypes::F_phi<float> >   fb(hkls, cxtl), fd(hkls, cxtl);
	    clipper::HKL_data<clipper::datatypes::Phi_fom<float> > phiw(hkls, cxtl);
	    clipper::HKL_data<clipper::datatypes::Flag>    flag(hkls, cxtl);
	    typedef clipper::HKL_data_base::HKL_reference_index HRI;
	    // If no free flag is available, then use all reflections..
	    for (HRI ih = flag.first(); !ih.last(); ih.next() )
	       flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;

	    /* This code uses free reflections only for sigmaa and scaling...
	       for (HRI ih = flag.first(); !ih.last(); ih.next() )
	       if ( !fo[ih].missing() && (free[ih].missing()||free[ih].flag()==0) )
	       flag[ih].flag() = clipper::SFweight_spline<float>::BOTH;
	       else
	       flag[ih].flag() = clipper::SFweight_spline<float>::NONE;
	    */

	    // do sigmaa calc
	    int n_refln = 1000;
	    int n_param = 20;
	    clipper::SFweight_spline<float> sfw(n_refln, n_param);
	    sfw( fb, fd, phiw, fo, fc, flag );
	    // fb is F+phi for "Best"
	    // fd is F+phi for difference map
	    short int is_diff = 0;
	    if (sigmaa_map_type == molecule_map_type::TYPE_DIFF_SIGMAA) {
	       map_fphidata = fd;
	       is_diff = 1;
	    } else {
	       map_fphidata = fb;
	    }


	    // 20091101 This fails to give a sensible cell, spacegroup
	    // and sampling for original_fphis.  Needs Kevin.
	    //
	    // original_fphis_filled = 1;
	    // original_fphis.init(map_fphidata.spacegroup(), map_fphidata.cell(), map_fphidata.hkl_sampling());
	    // original_fphis = map_fphidata;


	    // back to old code
	    //
	    cout << "initializing map...";
	    xmap.init(mydata.spacegroup(),
			      mydata.cell(),
			      clipper::Grid_sampling(mydata.spacegroup(),
						     mydata.cell(),
						     mydata.resolution(),
						     graphics_info_t::map_sampling_rate));
	    cout << "done."<< endl;

	    cout << "doing fft..." ;
	    // xmap.fft_from( fphidata );       // generate Fc alpha-c map
	    xmap.fft_from( map_fphidata );       // generate sigmaA map 20050804
	    cout << "done." << endl;
	    initialize_map_things_on_read_molecule(mol_name, is_diff, false, false);
	    // now need to fill contour_level, xmap_is_diff_map xmap_is_filled
	    if (is_diff)
	       xmap_is_diff_map = 1;
	    else
	       xmap_is_diff_map = 0;

	    mean_and_variance<float> mv = map_density_distribution(xmap, 40, false);

	    cout << "Mean and sigma of map from CIF file (make_map_from_cif): "
		 << mv.mean << " and " << sqrt(mv.variance) << endl;

	    update_map_in_display_control_widget();

	    map_mean_  = mv.mean;
	    map_sigma_ = sqrt(mv.variance);
	    map_max_   = mv.max_density;
	    map_min_   = mv.min_density;

	    set_initial_contour_level();

	    int imol = imol_no_in;
	    update_map();

	    if (sigmaa_map_type != molecule_map_type::TYPE_DIFF_SIGMAA) {
	       save_state_command_strings_.push_back("read-cif-data-with-phases-sigmaa");
	       save_state_command_strings_.push_back(single_quote(cif_file_name));
	    } else {
	       save_state_command_strings_.push_back("read-cif-data-with-phases-diff-sigmaa");
	       save_state_command_strings_.push_back(single_quote(cif_file_name));
	    }
	    return imol;
	 }
      }
   }
   catch (const clipper::Message_base &rte) {
      std::cout << "WARNING:: Problem reading " << cif_file_name << std::endl;
   }
   return -1;
}


//
// This is the version that gets called when we use a file selector to
// get the file (i.e. it doesn't specify the coordinates (molecule)) because
// this cif file has (or it is hoped that it has) calculated structure factors.
int
molecule_class_info_t::make_map_from_cif_nfofc(int imol_no_in,
					       std::string cif_file_name,
					       int map_type,
					       short int swap_difference_map_colours) {

   int ir = -1;
   imol_no = imol_no_in;

   clipper::HKL_info mydata;
   clipper::CIFfile cif;

   cif.open_read(cif_file_name);
   cif.import_hkl_info(mydata); // set spacegroup, cell and get hkl list.
   clipper::HKL_data< clipper::datatypes::F_sigF<float> >   fsigf(mydata);
   clipper::HKL_data< clipper::datatypes::F_phi<float> > fphidata(mydata);

   cif.import_hkl_data(fsigf);
   cif.import_hkl_data(fphidata);

   cif.close_read();

   std::cout << "Read " << mydata.num_reflections() << " from CIF file."
	     << std::endl;

   if (mydata.num_reflections() == 0) {
      return -1;
   } else {

      int non_zero = 0;
      for(int i=0; i< mydata.num_reflections(); i++) {
	 if (! fphidata[i].missing() ) {
	    if (fphidata[i].f() > 0.0) {
	       non_zero++;
	       break;
	    }
	 }
      }

      if (! non_zero) {
	 std::cout << "WARNING:: Ooops - all the calculated structure factor "
		   << "amplitudes appear"
		   << " to be zero - or missing.  " << std::endl;
	 std::cout << "WARNING:: Are you sure this file (" << cif_file_name
		   << ") contains calculated structure factors?" << std::endl;
	 std::cout << "WARNING:: No map calculated." << std::endl;
      } else {

	 std::string mol_name = cif_file_name;

	 int is_diff_map_flag = 0;
	 if (map_type == molecule_map_type::TYPE_FO_FC) {
	    is_diff_map_flag = 1;
	    mol_name += " Fo-Fc";
	 }
	 if (map_type == molecule_map_type::TYPE_2FO_FC) {
	    mol_name += " 2Fo-Fc";
	 }
	 if (map_type == molecule_map_type::TYPE_FO_ALPHA_CALC) {
	    mol_name += " Fo ac";
	 }

	 bool is_anomalous_flag = false;
	 initialize_map_things_on_read_molecule(mol_name, is_diff_map_flag, is_anomalous_flag,
						swap_difference_map_colours);

	 cout << "initializing map...";
	 xmap.init(mydata.spacegroup(),
			   mydata.cell(),
			   clipper::Grid_sampling(mydata.spacegroup(),
						  mydata.cell(),
						  mydata.resolution(),
						  graphics_info_t::map_sampling_rate));
	 std::cout << "done."<< std::endl;

	 // Here we need to fix up fphidata to be a combination
	 // of fsigf data and fphidata.
	 //
	 float fo_multiplier = 2.0;
	 float fc_multiplier = 1.0;
	 if (map_type == molecule_map_type::TYPE_FO_FC)
	    fo_multiplier = 1.0;
	 if (map_type == molecule_map_type::TYPE_FO_ALPHA_CALC) {
	    fo_multiplier = 1.0;
	    fc_multiplier = 0.0;
	 }

	 int nprm = 10;
	 std::vector<clipper::ftype> params_init(nprm, 1.0);
	 clipper::BasisFn_spline basis_f1f2( mydata, nprm, 2.0 );
	 clipper::TargetFn_scaleF1F2<clipper::datatypes::F_phi<float>,
	    clipper::datatypes::F_sigF<float> >
	    target_f1f2(fphidata, fsigf);
	 clipper::ResolutionFn fscale(mydata, basis_f1f2, target_f1f2, params_init);

	 int nrefl = 0;
	 int nmissing = 0;
	 for (clipper::HKL_info::HKL_reference_index ih=fsigf.first();
	      !ih.last(); ih.next()) {
	    nrefl++;
	    if (!fsigf[ih].missing()) {
	       fphidata[ih].f() = fo_multiplier * fsigf[ih].f() -
		  fc_multiplier * fphidata[ih].f() * sqrt(fscale.f(ih));
	       // std::cout << "scale: " << sqrt(fscale.f(ih)) << std::endl;
	    } else {
	       nmissing++;
	       // std::cout << "missing reflection: " << ih << std::endl;
	       fphidata[ih].f() = 0;
	    }
	 }
	 std::cout << "There were " << nrefl << " reflections of which "
		   << nmissing << " were missing\n";


	 std::cout << "doing fft..." ;
	 xmap.fft_from( fphidata );                  // generate map
	 std::cout << "done." << std::endl;

	 mean_and_variance<float> mv = map_density_distribution(xmap, 40, false);

	 std::cout << "Mean and sigma of map from CIF file (make_map_from_cif_nfofc): "
		   << mv.mean << " and " << sqrt(mv.variance) << std::endl;

	 if (is_diff_map_flag == 1) {
	    contour_level = nearest_step(mv.mean + 2.5*sqrt(mv.variance), 0.01);
	 } else {
	    contour_level = nearest_step(mv.mean + 1.5*sqrt(mv.variance), 0.05);
	 }

	 // fill class variables
	 map_mean_ = mv.mean;
	 map_sigma_ = sqrt(mv.variance);
	 xmap_is_diff_map = is_diff_map_flag;

	 int imol = imol_no_in;
	 update_map_in_display_control_widget();

	 update_map();

	 have_unsaved_changes_flag = 0;
	 std::vector<std::string> strings;
	 if (map_type == molecule_map_type::TYPE_FO_FC)
	    strings.push_back("read-cif-data-with-phases-fo-fc");
	 else
	    strings.push_back("read-cif-data-with-phases-2fo-fc");
	 strings.push_back(single_quote(cif_file_name));
	 save_state_command_strings_ = strings;

	 return imol;
      }
   }
   return ir;
}

int
molecule_class_info_t::make_map_from_mtz_by_calc_phases(int imol_no_in,
							const std::string &mtz_file_name,
							const std::string &f_col,
							const std::string &sigf_col,
							atom_selection_container_t SelAtom,
							short int is_2fofc_type) {

   clipper::CCP4MTZfile mtz;

   std::cout << "INFO:: reading mtz file..." << mtz_file_name << std::endl;
   mtz.open_read(mtz_file_name);

   // make the data names for import:
   std::pair<std::string, std::string> p = make_import_datanames(f_col, sigf_col, "", 0);
   clipper::HKL_data< clipper::datatypes::F_sigF<float> > myfsigf;
   mtz.import_hkl_data(myfsigf, p.first);
   mtz.close_read();

   return calculate_sfs_and_make_map(imol_no_in,
				     mtz_file_name, myfsigf,
				     SelAtom, is_2fofc_type);
}



// The rest was all interface fluff.  Here is where we do the real work
// (or get clipper to do it :).
//
int
molecule_class_info_t::make_map_from_phs(const clipper::Spacegroup &sg,
                                         const clipper::Cell &cell,
                                         std::string phs_filename) {

   // clipper::Resolution resolution(reso);  // no.

   // clipper::HKL_info mydata(sg, cell, resolution);

   clipper::PHSfile phs;

   if (! coot::file_exists(phs_filename)) {
      std::cout << "INFO:: file " << phs_filename << " does not exit " << std::endl;
      return -1;
   }

   try {
      std::cout << "INFO:: reading phs file: " << phs_filename << std::endl;
      phs.open_read(phs_filename);

      std::cout << "INFO:: phs: creating resolution" << std::endl;
      clipper::Resolution resolution = phs.resolution(cell);
      // mydata.init(sg, cell, resolution);

      std::cout << "PHS:: creating mydata" << std::endl;
      clipper::HKL_info mydata(sg, cell, resolution);
      clipper::HKL_data<clipper::datatypes::F_sigF<float>  >  myfsig(mydata);
      clipper::HKL_data<clipper::datatypes::Phi_fom<float> >  myphwt(mydata);
      clipper::HKL_data<clipper::datatypes::F_phi<float>   >  fphidata(mydata);

      std::cout << "INFO:: phs: importing info" << std::endl;
      phs.import_hkl_info(mydata);
      std::cout << "INFO:: phs: importing data" << std::endl;
      phs.import_hkl_data(myfsig);
      phs.import_hkl_data(myphwt);

      phs.close_read();

      std::cout << "INFO phs: using cell and symmetry: "
		<< cell.descr().a() << " "
		<< cell.descr().b() << " "
		<< cell.descr().c() << " "
		<< clipper::Util::rad2d(cell.descr().alpha()) << " "
		<< clipper::Util::rad2d(cell.descr().beta())  << " "
		<< clipper::Util::rad2d(cell.descr().gamma()) << " "
		<< single_quote(sg.symbol_hm()) << std::endl;

      std::cout << "INFO:: phs: number of reflections: " << mydata.num_reflections()
		<< "\n";

      fphidata.update();

      int ncount = 0;
      clipper::HKL_info::HKL_reference_index hri;
      //    for (hri=myfsig.first(); !hri.last(); hri.next()) {
      //       if (ncount < 300)
      // 	 std::cout << " PHS fsigf: " << hri.hkl().h() << " "
      // 		   << hri.hkl().k() << " "
      // 		   << hri.hkl().l() << " " << myfsig[hri].f() << " "
      // 		   << (myfsig[hri].sigf()) << std::endl;
      //       ncount++;
      //    }

      ncount = 0;
      //    for (hri=myphwt.first(); !hri.last(); hri.next()) {
      //       if (ncount < 300)
      // 	 std::cout << " PHS myphwt: " << hri.hkl().h() << " " << hri.hkl().k() << " "
      // 		   << hri.hkl().l() << " " << myphwt[hri].fom() << " "
      // 		   << clipper::Util::rad2d(myphwt[hri].phi()) << std::endl;
      //       ncount++;
      //    }

      fphidata.compute(myfsig, myphwt,
		       clipper::datatypes::Compute_fphi_from_fsigf_phifom<float>());

      //    for (int i=0; i<10; i++) {
      //       std::cout << "checking phi weight: " << i << " " << myphwt[i].phi() << "  "
      //              << myphwt[i].fom() << std::endl;
      //        std::cout << "checking f    sigf: " << i << " " << myfsig[i].f() << "   "
      // 		 << myfsig[i].sigf() << std::endl;
      //        std::cout << "checking missing: " << i << " " << myfsig[i].missing() << " "
      // 		 << myphwt[i].missing() << " " << fphidata[i].missing() << std::endl;
      //        // << " " << fphidata[i].phi() <<
      //    }

      std::string mol_name = phs_filename;

      initialize_map_things_on_read_molecule(mol_name, false, false, false); // not diff map

      std::cout << "initializing map...";
      xmap.init(mydata.spacegroup(),
                mydata.cell(),
                clipper::Grid_sampling(mydata.spacegroup(),
                                       mydata.cell(),
                                       mydata.resolution(),
                                       graphics_info_t::map_sampling_rate));
      std::cout << "done."<< std::endl;

      if (0) {
	 std::cout << "PHS:: debug:: " << mydata.spacegroup().symbol_hm() << " "
		   << mydata.cell().descr().a() << " "
		   << mydata.cell().descr().b() << " "
		   << mydata.cell().descr().c() << " "
		   << clipper::Util::rad2d(mydata.cell().descr().alpha()) << " "
		   << clipper::Util::rad2d(mydata.cell().descr().beta ()) << " "
		   << clipper::Util::rad2d(mydata.cell().descr().gamma()) << std::endl;
	 std::cout << "PHS:: debug:: n_reflections: " << mydata.num_reflections()
		   << std::endl;
      }

      ncount = 0;
      // clipper::HKL_info::HKL_reference_index hri;
      //   for (hri=fphidata.first(); !hri.last(); hri.next()) {
      //      if (ncount < 300)
      // 	std::cout << " PHS fphi: " << hri.hkl().h() << " " << hri.hkl().k() << " "
      // 		  << hri.hkl().l() << " " << fphidata[hri].f() << " "
      // 		  << clipper::Util::rad2d(fphidata[hri].phi()) << std::endl;
      //      ncount++;
      //   }


      //   cout << "Map Grid (from phs file)..."
      //        << xmap.grid_sampling().format()
      //        << endl;

      cout << "doing fft..." ;
      xmap.fft_from( fphidata );                  // generate map
      cout << "done." << endl;

      mean_and_variance<float> mv = map_density_distribution(xmap, 40, false);

      cout << "Mean and sigma of map from PHS file: " << mv.mean
	   << " and " << sqrt(mv.variance) << endl;

      // fill class variables
      map_mean_ = mv.mean;
      map_sigma_ = sqrt(mv.variance);

      xmap_is_diff_map = 0;
      contour_level = nearest_step(mv.mean + 1.5*sqrt(mv.variance), 0.05);
      update_map_in_display_control_widget();

      std::cout << "updating map..." << std::endl;
      update_map();
      std::cout << "done updating map..." << std::endl;
   }

   catch (...) {
      std::cout << "INFO:: problem reading phs file " << phs_filename << std::endl;
   }

  // as for 'normal' maps
  std::string cwd = coot::util::current_working_dir();
  std::string f1  = coot::util::intelligent_debackslash(phs_filename);
  std::string f2  = coot::util::relativise_file_name(f1, cwd);
  // how do we restore this map?
  save_state_command_strings_.push_back("read-phs-and-make-map-using-cell-symm");
  save_state_command_strings_.push_back(single_quote(f2));
  save_state_command_strings_.push_back(single_quote(sg.symbol_hm()));
  save_state_command_strings_.push_back(coot::util::float_to_string(cell.descr().a()));
  save_state_command_strings_.push_back(coot::util::float_to_string(cell.descr().b()));
  save_state_command_strings_.push_back(coot::util::float_to_string(cell.descr().c()));
  save_state_command_strings_.push_back(coot::util::float_to_string(clipper::Util::rad2d(cell.descr().alpha())));
  save_state_command_strings_.push_back(coot::util::float_to_string(clipper::Util::rad2d(cell.descr().beta())));
  save_state_command_strings_.push_back(coot::util::float_to_string(clipper::Util::rad2d(cell.descr().gamma())));

  return imol_no;
}

void
molecule_class_info_t::fill_skeleton_treenodemap() {

   // if we have a skeleton map but not treenodemap:
   //
   if (xskel_is_filled && !skeleton_treenodemap_is_filled) {

      // Chomp up the lovely memory! Yum!
      //
      skeleton_treenodemap.init(xskel_cowtan.spacegroup(),
				xskel_cowtan.cell(),
				xskel_cowtan.grid_sampling());
      clipper::Coord_grid c_g;
      clipper::Skeleton_basic::Neighbours skel_neighbs(xskel_cowtan);

//       std::cout << "Build tree: there are " << skel_neighbs.size() << " skel_neighbs"
// 		<< std::endl;  18, actually.

      clipper::Xmap_base::Map_reference_index ix;
      for (ix = xskel_cowtan.first(); !ix.last(); ix.next() ) {
	 if (xskel_cowtan[ix] > 0) {

	    coot::SkeletonTreeNode stn;

	    for(int i=0; i< skel_neighbs.size(); i++) {
	       c_g = ix.coord() + skel_neighbs[i];

	       if (xskel_cowtan.get_data(c_g) > 0 ) {

		  // OK, so this node has a neighbour:
		  //
		  stn.neighbs.push_back(c_g);
	       }
	    }
	    stn.near_grid_point = ix.coord();  // Strange but true!
	    //
	    // We do this because "out of cell" reference
	    // (e.g.  uvw = (  -1, -12, -19)) will get wrapped
	    // to some (hidden) value.  To get the wrapped
	    // value (i.e the grid), we look it up here.
	    // Cunning (if it works).
	    skeleton_treenodemap[ix] = stn;
	 }
      }
      // set the flag
      skeleton_treenodemap_is_filled = 1;
   }
}


float
molecule_class_info_t::density_at_point(const clipper::Coord_orth &co) const {

   if (xmap.is_null()) {
      std::cout << "WARNING:: null map. Returning bogus value from density_at_point()" << std::endl;
      return -1000.0;
   } else {

#ifdef HAVE_GSL
      float dv;
      clipper::Coord_frac af = co.coord_frac(xmap.cell());
      clipper::Coord_map  am = af.coord_map(xmap.grid_sampling());
      clipper::Interp_linear::interp(xmap, am, dv);
      return dv;
#else
      printf("no GSL so no density at point - remake \n");
      return -1000.0;
#endif // HAVE_GSL
   }
}
// Return status, was the contour level changed?  In that way, we
// don't try to recontour (which is a slow process) when the contour
// level has not been changed.
//
// We don't change the contour level if the contour level goes too
// low (typically below 0).
//
// We don't change the contour level if the contour level goes too
// high (above the maximum level of the map).
//
short int
molecule_class_info_t::change_contour(int direction) {

   short int istat = 0;
   // std::cout << "DEBUG:: contour_by_sigma_flag " << contour_by_sigma_flag << std::endl;
   // std::cout << "DEBUG:: adding " << contour_sigma_step << " * " << map_sigma_
   // << " to  " << contour_level << std::endl;
   if (has_xmap() || has_nxmap()) {

      float shift = graphics_info_t::diff_map_iso_level_increment;
      if (contour_by_sigma_flag) {
	 shift = contour_sigma_step * map_sigma_;
      } else {
	 if (xmap_is_diff_map) {
	    shift = graphics_info_t::diff_map_iso_level_increment;
	 } else {
	    shift = graphics_info_t::iso_level_increment;
	 }
      }

      if (xmap_is_diff_map) {
	 if (direction == -1) {
	    if (graphics_info_t::stop_scroll_diff_map_flag) {
	       if ((contour_level - shift) >
		   graphics_info_t::stop_scroll_diff_map_level) {
		  contour_level -= shift;
		  istat = 1;
	       }
	    } else {
	       contour_level -= shift;
	       istat = 1;
	    }
	 } else {
	    // add, but don't go past the top of the map or the bottom of the map
	    //
	    if (contour_level <= map_max_ || contour_level <= -map_min_) {
	       contour_level += shift;
	       istat = 1;
	    }
	 }
      } else {
	 // iso map

	 if (direction == -1) {
	    if (graphics_info_t::stop_scroll_iso_map_flag && ! is_patterson) {
	       if ((contour_level - shift) > graphics_info_t::stop_scroll_iso_map_level) {
		  contour_level -= shift;
		  istat = 1;
	       }
	    } else {
	       contour_level -= shift;
	       istat = 1;
	    }
	 } else {
	    if (contour_level <= map_max_) {
	       contour_level += shift;
	       istat = 1;
	    }
	 }
      }
   }
   return istat;
}

//
void
molecule_class_info_t::set_map_is_difference_map() {

   if (has_xmap() || has_nxmap()) {
      xmap_is_diff_map = 1;
      // we should update the contour level...
      set_initial_contour_level();
      // and set the right colors
      if (graphics_info_t::swap_difference_map_colours != 1) {
         map_colour.red   = 0.2;
         map_colour.green = 0.6;
         map_colour.blue  = 0.2;
      } else {
         map_colour.red   = 0.6;
         map_colour.green = 0.2;
         map_colour.blue  = 0.2;
      }
      update_map();
   }
}

short int
molecule_class_info_t::is_difference_map_p() const {

   short int istat = 0;
   if (has_xmap() || has_nxmap())
      if (xmap_is_diff_map)
	 istat = 1;
   return istat;
}


void
molecule_class_info_t::set_contour_by_sigma_step(float v, short int state) {
   contour_by_sigma_flag = state;
   if (state)
      contour_sigma_step = v;
}


// jiggle residue
float
molecule_class_info_t::fit_to_map_by_random_jiggle(coot::residue_spec_t &spec,
						   const clipper::Xmap<float> &xmap,
						   float map_sigma,
						   int n_trials,
						   float jiggle_scale_factor) {

   float v = -999.0;
   mmdb::Residue *residue_p = get_residue(spec);
   if (residue_p) {
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      bool use_biased_density_scoring = true;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      std::vector<mmdb::Chain *> chains; // empty - apply RTop to atoms of selection
      v = fit_to_map_by_random_jiggle(residue_atoms, n_residue_atoms, xmap, map_sigma,
                                      n_trials, jiggle_scale_factor, use_biased_density_scoring, chains);
   } else {
      std::cout << "WARNING:: residue " << spec << " not found" << std::endl;
   }
   return v;
}

// Sort so that the biggest numbers are at the top (lowest index) of the sorted list
bool
trial_results_comparer(const std::pair<clipper::RTop_orth, float> &a,
                       const std::pair<clipper::RTop_orth, float> &b) {

   return (b.second < a.second);

}

#ifdef HAVE_CXX_THREAD

float
molecule_class_info_t::fit_chain_to_map_by_random_jiggle(const std::string &chain_id, const clipper::Xmap<float> &xmap,
                                                         float map_sigma, int n_trials, float jiggle_scale_factor) {
   float r = 0;

   mmdb::PPAtom atom_selection = 0;
   int n_atoms;

   // If we have more than 20 residues, lets do an atom selection and use that
   // for fitting rather all the atoms.  No Side chains for protein,
   //
   std::pair<unsigned int, unsigned int> n_residues(0,0);
   mmdb::Chain *chain_p = 0;
   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p_this = model_p->GetChain(ichain);
         std::string chain_id_this(chain_p_this->GetChainID());
         if (chain_id_this == std::string(chain_id)) {
            chain_p = chain_p_this;
            break;
         }
      }
   }
   if (chain_p)
      n_residues = coot::util::get_number_of_protein_or_nucleotides(chain_p);

   int SelHnd = atom_sel.mol->NewSelection(); // d

   if (n_residues.first > 20) {
      atom_sel.mol->SelectAtoms(SelHnd, 0, chain_id.c_str(),
         mmdb::ANY_RES, "*",
         mmdb::ANY_RES, "*", "*",
         "CA,C,O,N","*","*",mmdb::SKEY_NEW);
   } else {
      if (n_residues.second > 20) {
         atom_sel.mol->SelectAtoms(SelHnd, 0, chain_id.c_str(),
               mmdb::ANY_RES, "*",
               mmdb::ANY_RES, "*", "*",
               "P,C1',N1,C2,N3,C4,N4,O2,C5,C6,O4,N9,C8,N7,N6","*","*",mmdb::SKEY_NEW);
      } else {
               atom_sel.mol->SelectAtoms(SelHnd, 0, chain_id.c_str(),
                  mmdb::ANY_RES, "*",
                  mmdb::ANY_RES, "*",
                  "*","*","*","*",mmdb::SKEY_NEW);
      }
   }

   atom_sel.mol->GetSelIndex(SelHnd, atom_selection, n_atoms);

   if (n_atoms) {
      bool use_biased_density_scoring = false; // not for all-molecule
      std::vector<mmdb::Chain *> chains;
      chains.push_back(chain_p);
      r = fit_to_map_by_random_jiggle(atom_selection, n_atoms,
                                      xmap, map_sigma,
                                      n_trials, jiggle_scale_factor,
                                      use_biased_density_scoring,
                                      chains);
   }
   atom_sel.mol->DeleteSelection(SelHnd);
   return r;
}

// static
void
molecule_class_info_t::test_jiggle_fit_func(unsigned int thread_index,
					    unsigned int i_trial,
					    unsigned int n_trials,
					    mmdb::PPAtom atom_selection,
					    int n_atoms,
					    const std::vector<mmdb::Atom *> &initial_atoms,
					    const clipper::Coord_orth &centre_pt,
					    const std::vector<std::pair<std::string, int> > &atom_numbers,
					    const clipper::Xmap<float> *xmap_masked,
					    float jiggle_scale_factor) {

}

// static
void
molecule_class_info_t::jiggle_fit_multi_thread_func_1(int thread_index,
						      unsigned int i_trial,
						      unsigned int n_trials,
						      mmdb::PPAtom atom_selection,
						      int n_atoms,
						      const std::vector<mmdb::Atom *> &initial_atoms,
						      const clipper::Coord_orth &centre_pt,
						      float jiggle_translation_scale_factor,
						      const std::vector<std::pair<std::string, int> > &atom_numbers,
						      const clipper::Xmap<float> *xmap_masked_p,
						      float (*density_scoring_function)(const coot::minimol::molecule &mol,
											const std::vector<std::pair<std::string, int> > &atom_number_list,
											const clipper::Xmap<float> &map),
						      std::pair<clipper::RTop_orth, float> *trial_results_p) {
   molecule_class_info_t m;
   float annealing_factor = 1.0 - static_cast<float>(i_trial)/static_cast<float>(n_trials);
   std::pair<clipper::RTop_orth, std::vector<mmdb::Atom> > jiggled_atoms =
      coot::util::jiggle_atoms(initial_atoms, centre_pt, jiggle_translation_scale_factor, annealing_factor);
   coot::minimol::molecule jiggled_mol(atom_selection, n_atoms, jiggled_atoms.second);
   float this_score = density_scoring_function(jiggled_mol, atom_numbers, std::cref(*xmap_masked_p));
   std::pair<clipper::RTop_orth, float> p(jiggled_atoms.first, this_score);
   *trial_results_p = p;
}
#endif // HAVE_CXX_THREAD

#ifdef HAVE_CXX_THREAD
// static
void
molecule_class_info_t::jiggle_fit_multi_thread_func_2(int thread_index,
						      const coot::minimol::molecule &direct_mol,
						      const clipper::Xmap<float> &xmap_masked,
						      float map_sigma,
						      const clipper::Coord_orth &centre_pt,
						      const std::vector<std::pair<std::string, int> > &atom_numbers,
						      float trial_results_pre_fit_score_for_trial,
						      float (*density_scoring_function)(const coot::minimol::molecule &mol,
											const std::vector<std::pair<std::string, int> > &atom_number_list,
											const clipper::Xmap<float> &map),
						      std::pair<clipper::RTop_orth, float> *post_fit_scores_p) {

   coot::minimol::molecule trial_mol = direct_mol;
   trial_mol.transform(post_fit_scores_p->first, centre_pt);
   float pre_score = density_scoring_function(trial_mol, atom_numbers, xmap_masked);
   molecule_class_info_t m;
   coot::minimol::molecule fitted_mol = m.rigid_body_fit(trial_mol, xmap_masked, map_sigma);
   // sorting and selection works by sorting the score of fitted_mols.
   float this_score = density_scoring_function(fitted_mol, atom_numbers, std::cref(xmap_masked));
   std::cout << "jiggle_fit_multi_thread_func_2() thread_index " << std::setw(2) << thread_index
	     << " pre-score " << std::setw(5) << pre_score
	     << " post-fit-score " << std::setw(5) << this_score << std::endl;
   post_fit_scores_p->second = this_score; // hand the score back
}
#endif // HAVE_CXX_THREAD


// called by above and split_water.
//
// chain_for_moving is default arg, with value 0.
//
// if chain_for_moving is not empty, apply the transformation
// the the atoms of chain_for_moving rather than to the atom of atom_selection
//
float
molecule_class_info_t::fit_to_map_by_random_jiggle(mmdb::PPAtom atom_selection,
                                                   int n_atoms,
                                                   const clipper::Xmap<float> &xmap,
                                                   float map_sigma,
                                                   int n_trials,
                                                   float jiggle_scale_factor,
                                                   bool use_biased_density_scoring,
                                                   std::vector<mmdb::Chain *> chains_for_moving) {
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
      std::pair<bool, clipper::Coord_orth> pt = residue_centre(res_ref);
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
#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
   int n_threads = coot::get_max_number_of_threads();
   if (n_threads > 0)
      do_multi_thread = true;
#endif

   if (do_multi_thread) {
#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY

      try {
         unsigned int n_threads = coot::get_max_number_of_threads();

         for (int itrial=0; itrial<n_trials; itrial++) {

            auto tp_1 = std::chrono::high_resolution_clock::now();

            // throw into the mix a model that has been only small rotated/translate
            // or maybe nothing at all.

            graphics_info_t::static_thread_pool.push(jiggle_fit_multi_thread_func_1, itrial, n_trials, atom_selection, n_atoms,
						     initial_atoms, centre_pt, jiggle_scale_factor, atom_numbers,
						     &xmap_masked, // pointer arg, no std::ref()
						     density_scoring_function, &trial_results[itrial]);

            auto tp_2 = std::chrono::high_resolution_clock::now();
            auto d21 = chrono::duration_cast<chrono::microseconds>(tp_2 - tp_1).count();
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
	    if (graphics_info_t::static_thread_pool.n_idle() == graphics_info_t::static_thread_pool.size())
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
      for (int i_trial=0; i_trial<n_for_rigid; i_trial++) {
         // does the fitting
         graphics_info_t::static_thread_pool.push(jiggle_fit_multi_thread_func_2, direct_mol,
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
         if (graphics_info_t::static_thread_pool.n_idle() == graphics_info_t::static_thread_pool.size())
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

   //
   if (bested) {
      make_backup();
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
                                                        COOT_LSQ_MAIN);
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

               // there were no chains for moving, so we apply the rtop to everything
               // with replace_coords().

               atom_selection_container_t asc_ligand = make_asc(mol);
               replace_coords(asc_ligand, false, true);
            }

            have_unsaved_changes_flag = 1;
            make_bonds_type_checked();
         } else {
            std::cout << "ERROR:: fit_to_map_by_random_jiggle(): mol is null! " << std::endl;
         }
         delete mol;
      } else {
         std::cout << "ERROR:: fit_to_map_by_random_jiggle(): best_molecule is empty!" << std::endl;
      }
   } else {
      std::cout << " nothing better found " << std::endl;
   }
   return v;
}

float
molecule_class_info_t::fit_molecule_to_map_by_random_jiggle(const clipper::Xmap<float> &xmap,
                                           float map_sigma, int n_trias, float jiggle_scale_factor) {
   float r = 0;

   // put the guts of fit_molecule_to_map_by_random_jiggle_and_blur() here.

   return r;
}


// return a fitted molecule
coot::minimol::molecule
molecule_class_info_t::rigid_body_fit(const coot::minimol::molecule &mol_in,
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

bool
molecule_class_info_t::map_is_too_blue_p() const {

   bool state = 0;

   if (has_xmap() || has_nxmap())
      if (! xmap_is_diff_map)
	      if (map_colour.red < 0.4)
	         if (map_colour.green < 0.4)
	            state = 1;

   std::cout << "Map is too blue: " << state << std::endl;
   return state;
}



map_statistics_t
molecule_class_info_t::map_statistics() const {

   double mean = 0;
   double sd   = 0;
   double skew     = 0;
   double kurtosis = 0;
   long   n = 0;

   double sum = 0;
   double sum_sq = 0;
   double sum_cubed = 0;
   double sum_4rd   = 0;

   // recall kurtosis, $k$ of $N$ observations:
   // k = \frac{\Sigma(x_i - \mu)^4} {N \sigma^4} - 3
   // (x_i - \mu)^4 = x_i^4 - 4x_i^3\mu + 6x_i^2\mu^2 - 4x_i\mu^3 + \mu^4


   clipper::Xmap_base::Map_reference_index ix;
   for (ix=xmap.first(); !ix.last(); ix.next()) {
      const float &rho = xmap[ix];
      if (! clipper::Util::is_nan(rho)) {
	 n++;
	 sum       += rho;
	 sum_sq    += rho * rho;
	 sum_cubed += rho * rho * rho;
	 sum_4rd   += rho * rho * rho * rho;
      }
   }

   if (n > 0) {
      double dn = double(n);
      mean = sum/dn;
      double v = sum_sq/dn - mean * mean;
      if (v < 0)
	 v = 0;
      sd = sqrt(v);
      skew = sum_cubed/dn - 3 * mean * v - mean * mean * mean;
      double kt =
	 sum_4rd
	 - 4 * sum_cubed * mean
	 + 6 * sum_sq    * mean * mean
	 - 4 * sum       * mean * mean * mean
	 + mean * mean * mean * mean * dn;
      kurtosis = kt/(dn * v * v);
   }
   map_statistics_t mp(mean, sd, skew, kurtosis);
   return mp;

}

std::vector<std::pair<clipper::Coord_orth, clipper::Coord_orth> >
molecule_class_info_t::get_contours(float contour_level,
                                    float radius,
                                    const coot::Cartesian &centre) const {

      // who calls this function?

   std::vector<std::pair<clipper::Coord_orth, clipper::Coord_orth> > r;

   bool is_em_map_local = false; // needs setting?

   int isample_step = 1;
   CIsoSurface<float> my_isosurface;
   // a pointer and a size
   coot::CartesianPairInfo v = my_isosurface.GenerateSurface_from_Xmap(xmap,
                                                                       contour_level,
                                                                       radius, centre,
                                                                       0,1,
                                                                       isample_step, is_em_map_local);
   if (v.data) {
      if (v.size > 0) {
         r.resize(v.size);
         for (int i=0; i<v.size; i++) {
            const coot::Cartesian &s = v.data[i].getStart();
            const coot::Cartesian &f = v.data[i].getFinish();
            clipper::Coord_orth p1(s.x(), s.y(), s.z());
            clipper::Coord_orth p2(f.x(), f.y(), f.z());
            r[i]= std::pair<clipper::Coord_orth, clipper::Coord_orth>(p1,p2);
         }
      }
   }
   return r;
}

// static
int
molecule_class_info_t::watch_mtz(gpointer data) {

   int status = 1; // continue

   updating_map_params_t *ump = static_cast<updating_map_params_t *>(data);
   const updating_map_params_t &rump = *ump;
   status = graphics_info_t::molecules[rump.imol].update_map_from_mtz_if_changed(rump);
   return status;
}

int
molecule_class_info_t::update_map_from_mtz_if_changed(const updating_map_params_t &ump_in) {

   int status = 1;
   if (continue_watching_mtz) {

      bool update_it = false;

      updating_map_params_t ump = ump_in;
      struct stat s;
      int status = stat(ump.mtz_file_name.c_str(), &s);
      if (status != 0) {
	 std::cout << "WARNING:: update_map_from_mtz_if_changed() Error reading "
		   << ump.mtz_file_name << std::endl;
      } else {
	 if (!S_ISREG (s.st_mode)) {
	    std::cout << "WARNING:: update_map_from_mtz_if_changed() not a reguular file: "
		      << ump.mtz_file_name << std::endl;
	    continue_watching_mtz = false;
	 } else {
	    // happy path
#ifndef _POSIX_SOURCE
	    ump.ctime = s.st_ctimespec; // mac?
#else
	    ump.ctime = s.st_ctim;
#endif
	 }
      }

      if (false)
	 std::cout << "#### ctime comparision: was "
		   << updating_map_previous.ctime.tv_sec << " " << updating_map_previous.ctime.tv_nsec
		   << " now " << ump.ctime.tv_sec << " " << ump.ctime.tv_nsec
		   << std::endl;

      if (ump.ctime.tv_sec > updating_map_previous.ctime.tv_sec) {
	 update_it = true;
      } else {
	 if (ump.ctime.tv_sec == updating_map_previous.ctime.tv_sec) {
	    if (ump.ctime.tv_nsec > updating_map_previous.ctime.tv_nsec) {
	       update_it = true;
	    }
	 }
      }
      if (update_it) {

	 // map_fill_from_mtz(ump) ?

	 // updating maps shouldn't update (add to) the Display Manager.

	 std::string cwd = coot::util::current_working_dir();
	 map_fill_from_mtz(ump.mtz_file_name,
			   cwd,
			   ump.f_col,
			   ump.phi_col,
			   ump.weight_col,
			   ump.use_weights,
			   ump.is_difference_map,
			   graphics_info_t::map_sampling_rate);
	 updating_map_previous = ump;
	 graphics_info_t::graphics_draw();
      }
   } else {
      status = 0;
   }
   return status;
}

#include "coot-utils/sfcalc-genmap.hh"

   // use this molecules mol and the passed data to make a map for some other
   // molecule
int
molecule_class_info_t::sfcalc_genmap(const clipper::HKL_data<clipper::data32::F_sigF> &fobs,
                                     const clipper::HKL_data<clipper::data32::Flag> &free,
                                     clipper::Xmap<float> *xmap_p) {
   coot::util::sfcalc_genmap(atom_sel.mol, fobs, free, xmap_p);
   return 0;
}

// Weird combination of passed variables and class variables.
GdkRGBA
molecule_class_info_t::radius_to_colour(float radius, float min_radius, float max_radius) {

   float f = 0.0;
   if (radius > min_radius) {
      if (radius > max_radius) {
         f = 1.0;
      } else {
         float range = max_radius - min_radius;
         f = (radius - min_radius)/range;
      }
   }
   if (radial_map_colour_invert_flag)
      f = 1.0 - f;
   return fraction_to_colour(f);
}

GdkRGBA
molecule_class_info_t::fraction_to_colour(float fraction) {
   GdkRGBA col;
   float sat = radial_map_colour_saturation;
   coot::colour_t cc(0.6+0.4*sat, 0.6-0.6*sat, 0.6-0.6*sat);
   // cc.rotate(1.05 * fraction); // blue end is a bit purple/indigo
   cc.rotate(0.66 * fraction);
   col.red   = cc.col[0];
   col.green = cc.col[1];
   col.blue  = cc.col[2];
   col.alpha = 1.0;

   return col;
}

void
molecule_class_info_t::colour_map_using_map(const clipper::Xmap<float> &xmap) {

  colour_map_using_other_map_flag = true;
  other_map_for_colouring_p = &xmap;
  update_map();

}


void
molecule_class_info_t::colour_map_using_map(const clipper::Xmap<float> &xmap, float table_bin_start, float table_bin_size,
                                            const std::vector<coot::colour_t> &colours) {

   if (colours.empty()) {
      std::cout << "WARNING:: no colours - no map colouring" << std::endl;
   } else {
      // we need to tell the glsl colour setup to use value_to_colour() (which uses our table)
      // (not done yet).
      colour_map_using_other_map_flag = true; // tell the triangle generator to use a function to get the colour
                                              // (position_to_colour_using_other_map(co))
      other_map_for_colouring_p = &xmap;
      other_map_for_colouring_min_value = table_bin_start;
      other_map_for_colouring_max_value = table_bin_start + colours.size() * table_bin_size;
      other_map_for_colouring_colour_table = colours;
      update_map();
   }

}

coot::density_contour_triangles_container_t
molecule_class_info_t::export_molecule_as_x3d() const {

   // return value should be a vector of vertices and triangle indices
   coot::density_contour_triangles_container_t tc;

   if (true) {

      unsigned int sum_tri_con_points = 0;
      unsigned int sum_tri_con_normals = 0;
      unsigned int sum_tri_con_triangles = 0;
      std::vector<coot::density_contour_triangles_container_t>::const_iterator it;
      for (it=draw_vector_sets.begin(); it!=draw_vector_sets.end(); it++) {
         const coot::density_contour_triangles_container_t &tri_con(*it);
         sum_tri_con_points    += tri_con.points.size();
         sum_tri_con_normals   += tri_con.normals.size();
         sum_tri_con_triangles += tri_con.point_indices.size();
      }
      if (xmap_is_diff_map) {
         for (it=draw_diff_map_vector_sets.begin(); it!=draw_diff_map_vector_sets.end(); it++) {
            const coot::density_contour_triangles_container_t &tri_con(*it);
            sum_tri_con_points    += tri_con.points.size();
            sum_tri_con_normals   += tri_con.normals.size();
            sum_tri_con_triangles += tri_con.point_indices.size();
         }
      }

      if (sum_tri_con_triangles > 0) {

         std::vector<int> idx_base_for_points_vec(   draw_vector_sets.size()); // and normals
         std::vector<int> idx_base_for_triangles_vec(draw_vector_sets.size());
         if (xmap_is_diff_map) {
            idx_base_for_points_vec.resize(draw_vector_sets.size() + draw_diff_map_vector_sets.size());
            idx_base_for_triangles_vec.resize(draw_vector_sets.size() + draw_diff_map_vector_sets.size());
         }
         for (unsigned int i=0; i<draw_vector_sets.size(); i++) {
            if (i==0) {
               idx_base_for_points_vec[i] = 0;
               idx_base_for_triangles_vec[i] = 0;
            } else {
               idx_base_for_points_vec[i]    = idx_base_for_points_vec[i-1]    + draw_vector_sets[i-1].points.size();
               idx_base_for_triangles_vec[i] = idx_base_for_triangles_vec[i-1] + draw_vector_sets[i-1].points.size();
            }
         }
         if (xmap_is_diff_map) {
            for (unsigned int i0=0; i0<draw_diff_map_vector_sets.size(); i0++) {
               unsigned int i = i0 + draw_vector_sets.size();
               const coot::density_contour_triangles_container_t &tri_con(draw_diff_map_vector_sets[i0]);
               try {
                  if (i0==0) {
                     idx_base_for_points_vec.at(i)    = idx_base_for_points_vec.at(i-1)    + draw_vector_sets.at(i-1).points.size();
                     idx_base_for_triangles_vec.at(i) = idx_base_for_triangles_vec.at(i-1) + draw_vector_sets.at(i-1).points.size();
                  } else {
                     idx_base_for_points_vec.at(i)    = idx_base_for_points_vec.at(i-1)    + draw_diff_map_vector_sets.at(i0-1).points.size();
                     idx_base_for_triangles_vec.at(i) = idx_base_for_triangles_vec.at(i-1) + draw_diff_map_vector_sets.at(i0-1).points.size();
                  }
               }
               catch (std::out_of_range &oor) {
                  std::cout << "ERROR:: caught out of range " << oor.what() << " at i " << i << std::endl;
               }
            }
         }

         std::cout << "resize points  " << sum_tri_con_points << std::endl;
         std::cout << "resize normals " << sum_tri_con_normals << " (should be the same)" << std::endl;
         tc.points.resize(sum_tri_con_points);
         tc.normals.resize(sum_tri_con_normals);
         tc.point_indices.resize(sum_tri_con_triangles);

         // Now transfer the points

         int idx_points = 0;
         for (unsigned int i=0; i<draw_vector_sets.size(); i++) {
            const coot::density_contour_triangles_container_t &tri_con(draw_vector_sets[i]);
            for (std::size_t j=0; j<tri_con.points.size(); j++) {
               tc.points[idx_points  ] = tri_con.points[j];
               idx_points++;
            }
         }
         if (xmap_is_diff_map) {
            for (unsigned int i0=0; i0<draw_diff_map_vector_sets.size(); i0++) {
               unsigned int i = i0 + draw_vector_sets.size();
               const coot::density_contour_triangles_container_t &tri_con(draw_diff_map_vector_sets[i0]);
               for (std::size_t j=0; j<tri_con.points.size(); j++) {
                  tc.points[idx_points] = tri_con.points[j];
                  idx_points++;
               }
            }
         }

         int n_indices_for_triangles    = 3 * sum_tri_con_triangles;

         int idx_for_triangles = 0;
         for (unsigned int i=0; i<draw_vector_sets.size(); i++) {
            const coot::density_contour_triangles_container_t &tri_con(draw_vector_sets[i]);
            int idx_base_for_triangles = idx_base_for_triangles_vec[i];
            for (std::size_t i=0; i<tri_con.point_indices.size(); i++) {
               // indices_for_triangles[3*idx_for_triangles  ] = idx_base_for_triangles + tri_con.point_indices[i].pointID[0];
               // indices_for_triangles[3*idx_for_triangles+1] = idx_base_for_triangles + tri_con.point_indices[i].pointID[1];
               // indices_for_triangles[3*idx_for_triangles+2] = idx_base_for_triangles + tri_con.point_indices[i].pointID[2];
               tc.point_indices[idx_for_triangles].pointID[0] = idx_base_for_triangles + tri_con.point_indices[i].pointID[0];
               tc.point_indices[idx_for_triangles].pointID[1] = idx_base_for_triangles + tri_con.point_indices[i].pointID[1];
               tc.point_indices[idx_for_triangles].pointID[2] = idx_base_for_triangles + tri_con.point_indices[i].pointID[2];
               idx_for_triangles++;
            }
         }
         if (xmap_is_diff_map) {
            for (unsigned int i0=0; i0<draw_diff_map_vector_sets.size(); i0++) {
               unsigned int i = i0 + draw_vector_sets.size();
               const coot::density_contour_triangles_container_t &tri_con(draw_diff_map_vector_sets[i0]);
               int idx_base_for_triangles = idx_base_for_triangles_vec[i];
               for (std::size_t i=0; i<tri_con.point_indices.size(); i++) {
                  // indices_for_triangles[3*idx_for_triangles  ] = idx_base_for_triangles + tri_con.point_indices[i].pointID[0];
                  // indices_for_triangles[3*idx_for_triangles+1] = idx_base_for_triangles + tri_con.point_indices[i].pointID[1];
                  // indices_for_triangles[3*idx_for_triangles+2] = idx_base_for_triangles + tri_con.point_indices[i].pointID[2];
                  tc.point_indices[idx_for_triangles].pointID[0] = idx_base_for_triangles + tri_con.point_indices[i].pointID[0];
                  tc.point_indices[idx_for_triangles].pointID[1] = idx_base_for_triangles + tri_con.point_indices[i].pointID[1];
                  tc.point_indices[idx_for_triangles].pointID[2] = idx_base_for_triangles + tri_con.point_indices[i].pointID[2];
                  idx_for_triangles++;
               }
            }
         }

         // each index has a normal

         int n_normals = sum_tri_con_normals;
         int idx_for_normals = 0;
         for (unsigned int i0=0; i0<draw_vector_sets.size(); i0++) {
            const coot::density_contour_triangles_container_t &tri_con(draw_vector_sets[i0]);
            for (std::size_t i=0; i<tri_con.normals.size(); i++) {
               // std::cout << "storing normal " << idx_for_normals << " " << tri_con.normals[i].format() << std::endl;
               tc.normals[idx_for_normals] = tri_con.normals[i];
               idx_for_normals++;
            }
         }
         if (xmap_is_diff_map) {
            for (unsigned int i0=0; i0<draw_diff_map_vector_sets.size(); i0++) {
               const coot::density_contour_triangles_container_t &tri_con(draw_diff_map_vector_sets[i0]);
               for (std::size_t i=0; i<tri_con.normals.size(); i++) {
                  tc.normals[idx_for_normals] = tri_con.normals[i];
                  idx_for_normals++;
               }
            }
         }
      }
   }
   return tc;
}
