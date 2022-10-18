
#include <thread>

#include <clipper/ccp4/ccp4_map_io.h>

#include "utils/coot-utils.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "density-contour/occlusion.hh"
#include "density-contour/transfer-occlusions.hh"
#include "coot_molecule.hh"

std::atomic<bool> coot::molecule_t::draw_vector_sets_lock(false);

int
coot::molecule_t::writeMap(const std::string &file_name) const {

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
               positions[i] = coot::augmented_position(pos, norm);
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

      coot::density_contour_triangles_container_t tri_con =
        my_isosurface.GenerateTriangles_from_Xmap(std::cref(*xmap_p),
                                                  contour_level, dy_radius, centre, isample_step,
                                                  iream_start, n_reams, is_em_map);

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
coot::molecule_t::get_map_contours_mesh(clipper::Coord_orth position, float radius, float contour_level) {


   std::cout << "!!! ######################################### get_map_contours_mesh() for imol " << imol_no << std::endl;

   auto coord_orth_to_glm = [] (const clipper::Coord_orth &co) {
      return glm::vec3(co.x(), co.y(), co.z());
   };

   coot::simple_mesh_t m;

   clipper::Coord_orth p(position.x(), position.y(), position.z());
   update_map_triangles(radius, p, contour_level);

   auto &vertices  = m.vertices;
   auto &triangles = m.triangles;

   // now convert the contents of the draw-vector sets to a simple_mesh_t.

   coot::colour_holder map_colour(0.4, 0.5, 0.8);

   std::vector<coot::density_contour_triangles_container_t>::const_iterator it;
   glm::vec4 col(map_colour.red, map_colour.green, map_colour.blue, 1.0f);
   for (it=draw_vector_sets.begin(); it!=draw_vector_sets.end(); ++it) {
      const coot::density_contour_triangles_container_t &tri_con(*it);
      unsigned int idx_base = vertices.size();
      for (unsigned int i=0; i<tri_con.points.size(); i++) {
         glm::vec3 pos    = coord_orth_to_glm(tri_con.points[i]);
         glm::vec3 normal = coord_orth_to_glm(tri_con.normals[i]);
         coot::api::vnc_vertex vert(pos, normal, col);
         vertices.push_back(vert);
      }
      for (unsigned int i=0; i<tri_con.point_indices.size(); i++) {
         g_triangle tri(tri_con.point_indices[i].pointID[0],
                        tri_con.point_indices[i].pointID[1],
                        tri_con.point_indices[i].pointID[2]);
         tri.rebase(idx_base);
         triangles.push_back(tri);

         // date removed map triangle centres block here.

      }
   }
   

   return m;
}

#include "coot-utils/peak-search.hh"

// the molecule is passed so that the peaks are placed around the protein
std::vector<coot::molecule_t::difference_map_peaks_info_t>
coot::molecule_t::difference_map_peaks(mmdb::Manager *mol, float n_rmsd) const {

   std::vector<coot::molecule_t::difference_map_peaks_info_t> v;
   if (mol) {
      coot::peak_search ps(xmap);
      std::vector<std::pair<clipper::Coord_orth, float> > peaks = ps.get_peaks(xmap, mol, n_rmsd, true, true, true);
      for (const auto &peak : peaks) {
         difference_map_peaks_info_t dmp(peak.first, peak.second);
         v.push_back(dmp);
      }
   } else {
      std::cout << "ERROR:: " << __FUNCTION__ << "() null mol" << std::endl;
   }
   return v;
}


#include "coot-utils/xmap-stats.hh"

// map functions, return -1.1 on not-a-map
float
coot::molecule_t::get_map_rmsd_approx() const {

   bool ignore_pseudo_zeros_for_map_stats = false; // set this to true for an EM map
   bool ipz = ignore_pseudo_zeros_for_map_stats;
   mean_and_variance<float> mv = map_density_distribution(xmap, 20, false, ipz);
   float rmsd = std::sqrt(mv.variance);
   return rmsd;

}

bool
coot::molecule_t::is_difference_map_p() const {

   short int istat = 0;
   if (is_valid_map_molecule())
      if (xmap_is_diff_map)
         istat = 1;
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

   // set original_fobs_sigfobs_filled when done

   bool have_sensible_refmac_params = true; // 20221016-PE need to be set properly!

   if (have_sensible_refmac_params) {

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
               std::cout << "INFO:: About to read " << Refmac_mtz_filename() << " with dataname " << dataname << std::endl;
               clipper::CCP4MTZfile *mtzin_rfree_p = new clipper::CCP4MTZfile;
               mtzin_rfree_p->open_read(Refmac_mtz_filename());
               mtzin_rfree_p->import_hkl_data(*original_r_free_flags_p, dataname);
               mtzin_rfree_p->close_read();

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
         std::cout << "Timings: read mtz file and store data " << d10 << " milliseconds" << std::endl;
      }
   } else {
      std::cout << "DEBUG:: fill_fobs_sigfobs() no Fobs parameters\n";
   }
}
