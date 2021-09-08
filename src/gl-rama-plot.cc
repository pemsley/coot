
#include <chrono>

#include "gl-rama-plot.hh"
#include "coot-utils/coot-coord-utils.hh"

void
gl_rama_plot_t::setup_from(int imol, mmdb::Manager *mol) {

   // auto tp_0 = std::chrono::high_resolution_clock::now();
   if (mol) {
      float position_hash_now = coot::get_position_hash(mol);
      // std::cout << "comparing hashes " << position_hash_now << " " << position_hash << std::endl;
      if (position_hash_now != position_hash) {
         std::map<coot::residue_spec_t, rama_plot::phi_psi_t> phi_psi_map = generate_phi_psis(imol, mol);
         update_hud_tmeshes(phi_psi_map);
         position_hash = position_hash_now;
      }
   }
   // auto tp_1 = std::chrono::high_resolution_clock::now();
   // auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
   // std::cout << "INFO:: setup_from() " << d10 << " milliseconds" << std::endl;
}

#include <glm/gtx/string_cast.hpp>

void
gl_rama_plot_t::update_hud_tmeshes(const std::map<coot::residue_spec_t, rama_plot::phi_psi_t> &phi_psi_map) {

   // auto tp_0 = std::chrono::high_resolution_clock::now();  ~10 microsceconds

   // hud_tmesh_for_other_normal;
   // hud_tmesh_for_other_outlier;
   // hud_tmesh_for_pro_normal;
   // hud_tmesh_for_pro_outlier;
   // hud_tmesh_for_gly_normal;
   // hud_tmesh_for_gly_outlier;

   float rama_plot_scale = 0.5; // pass this or make it a member item
   float sf = rama_plot_scale * 0.23; // 0.25 is too big, 0.23 is about right, 0.22 is a bit too small
   std::vector<glm::vec2> new_other_normal_positions;
   std::vector<glm::vec2> new_other_outlier_positions;
   std::vector<glm::vec2> new_pro_normal_positions;
   std::vector<glm::vec2> new_pro_outlier_positions;
   std::vector<glm::vec2> new_gly_normal_positions;
   std::vector<glm::vec2> new_gly_outlier_positions;
   std::map<coot::residue_spec_t, rama_plot::phi_psi_t>::const_iterator it;
   for (it=phi_psi_map.begin(); it!=phi_psi_map.end(); it++) {
      const auto &phi_psi = it->second; // phi and psi are in degrees
      double phi_r = clipper::Util::d2rad(phi_psi.phi);
      double psi_r = clipper::Util::d2rad(phi_psi.psi);
      glm::vec2 pos(sf * phi_psi.phi, sf * phi_psi.psi);
      // std::cout << "phi_psi pos: " << glm::to_string(pos) << std::endl;
      if (phi_psi.residue_name == "PRO") {
         // new_pro_normal_positions.push_back(pos);
      } else {
         if (phi_psi.residue_name == "GLY") {
            // new_gly_normal_positions.push_back(pos);
         } else {
            double probability = rama_non_gly_pro.probability(phi_r, psi_r);
            // bool is_outlier = (probability < rama_threshold_allowed);
            // std::cout << "debug:: probabilities " << probability << " " << rama_threshold_allowed << " " << is_outlier << std::endl;
            if (probability < rama_threshold_allowed) {
               new_other_outlier_positions.push_back(pos);
            } else {
               new_other_normal_positions.push_back(pos);
            }
         }
      }
   }

   // mark the corners of the plot for scaling and offset testing
   new_other_normal_positions.push_back(glm::vec2(sf * -180.0, sf * -180.0));
   new_other_normal_positions.push_back(glm::vec2(sf *  180.0, sf * -180.0));
   new_other_normal_positions.push_back(glm::vec2(sf * -180.0, sf *  180.0));
   new_other_normal_positions.push_back(glm::vec2(sf *  180.0, sf *  180.0));

   glm::vec2 plot_point_scales(0.012, 0.012); // this gets the marker point size right
   glm::vec2 offset(-5.66 * sf, -5.66 * sf); // -7 is too negative, -5.7 is too negative, -5.6 is too little

   hud_tmesh_for_other_normal.set_scales(plot_point_scales);
   hud_tmesh_for_other_normal.set_position(offset); // tweaked to match the axes mesh
   hud_tmesh_for_other_normal.update_instancing_buffer_data(new_other_normal_positions);

   hud_tmesh_for_pro_normal.set_scales(plot_point_scales);
   hud_tmesh_for_pro_normal.set_position(offset);
   hud_tmesh_for_pro_normal.update_instancing_buffer_data(new_pro_normal_positions);

   hud_tmesh_for_gly_normal.set_scales(plot_point_scales);
   hud_tmesh_for_gly_normal.set_position(offset);
   hud_tmesh_for_gly_normal.update_instancing_buffer_data(new_gly_normal_positions);

   hud_tmesh_for_other_outlier.set_scales(plot_point_scales);
   hud_tmesh_for_other_outlier.set_position(offset); // tweaked to match the axes mesh
   hud_tmesh_for_other_outlier.update_instancing_buffer_data(new_other_outlier_positions);

   hud_tmesh_for_pro_outlier.set_scales(plot_point_scales);
   hud_tmesh_for_pro_outlier.set_position(offset);
   hud_tmesh_for_pro_outlier.update_instancing_buffer_data(new_pro_outlier_positions);

   hud_tmesh_for_gly_outlier.set_scales(plot_point_scales);
   hud_tmesh_for_gly_outlier.set_position(offset);
   hud_tmesh_for_gly_outlier.update_instancing_buffer_data(new_gly_outlier_positions);

   // auto tp_1 = std::chrono::high_resolution_clock::now();
   // auto d10 = std::chrono::duration_cast<std::chrono::microseconds>(tp_1 - tp_0).count();
   // std::cout << "INFO:: update_hud_tmeshes() " << d10 << " microseconds" << std::endl;
}


void
gl_rama_plot_t::init() {

   current_background_type = clipper::Ramachandran::All;
   draw_outliers_only_flag = false;
   colour_for_OK_phi_psi     = "DodgerBlue";
   colour_for_OK_phi_psi_pro = "SkyBlue";
   colour_for_outlier_phi_psi = "red";

   rama_threshold_preferred = 0.02;
   rama_threshold_allowed   = 0.002;

   rama_all.init(clipper::Ramachandran::All);
   rama_gly.init(clipper::Ramachandran::Gly2);
   rama_pro.init(clipper::Ramachandran::Pro2);
   rama_non_gly_pro.init(clipper::Ramachandran::NonGlyPro);
   rama_ile_val.init(clipper::Ramachandran::IleVal2);
   rama_pre_pro.init(clipper::Ramachandran::PrePro2);
   rama_non_gly_pro_ile_val_or_pre_pro.init(clipper::Ramachandran::NoGPIVpreP2);

   rama_all.set_thresholds(            rama_threshold_preferred, rama_threshold_allowed);
   rama_gly.set_thresholds(            rama_threshold_preferred, rama_threshold_allowed);
   rama_pro.set_thresholds(            rama_threshold_preferred, rama_threshold_allowed);
   rama_non_gly_pro.set_thresholds(    rama_threshold_preferred, rama_threshold_allowed);
   rama_ile_val.set_thresholds(        rama_threshold_preferred, rama_threshold_allowed);
   rama_pre_pro.set_thresholds(        rama_threshold_preferred, rama_threshold_allowed);
   rama_non_gly_pro_ile_val_or_pre_pro.set_thresholds(rama_threshold_preferred, rama_threshold_allowed);

   position_hash = 0.0;
}

std::map<coot::residue_spec_t, rama_plot::phi_psi_t>
gl_rama_plot_t::generate_pseudo_phi_psis() {

   std::map<coot::residue_spec_t, rama_plot::phi_psi_t> m;

   const unsigned int n_spots = 100;
   for (unsigned int i=0; i<n_spots; i++) {
      coot::residue_spec_t rs("A", i+1, "");
      float f = static_cast<float>(i)/static_cast<float>(n_spots);
      float pi = 3.14159265;
      float phi = -60.0 * sinf(2.0 * pi * f) + 100.0 * sin(1000.0 * f);
      float psi = 180.0 * cosf(2.0 * pi * f);
      rama_plot::phi_psi_t pp(phi, psi);
      pp.label = std::string("A ") + std::to_string(i+1) + " GLN";
      m[rs] = pp;
   }
   return m;
}



std::map<coot::residue_spec_t, rama_plot::phi_psi_t>
gl_rama_plot_t::generate_phi_psis(int imol, mmdb::Manager *mol) {

   std::map<coot::residue_spec_t, rama_plot::phi_psi_t> r;
   int n_models = mol->GetNumberOfModels();
   // add a place-holder for the "0-th" model
   rama_plot::phi_psis_for_model_t empty(0);
   for (int imod=1; imod<=n_models; imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
	 int nchains = model_p->GetNumberOfChains();
	 for (int ichain=0; ichain<nchains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
	    int nres = chain_p->GetNumberOfResidues();
	    if (nres > 2) { 
	       for (int ires=1; ires<(nres-1); ires++) { 
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
		  mmdb::Residue *res_prev = coot::util::previous_residue(residue_p);
		  mmdb::Residue *res_next = coot::util::next_residue(residue_p);
		  if (res_prev && residue_p && res_next) {
		     try {
			// rama_plot::phi_psi_t constructor can throw an error
			// (e.g. bonding atoms too far apart).
			coot::residue_spec_t spec(residue_p);
			rama_plot::phi_psi_t pp(res_prev, residue_p, res_next);
                        pp.imol = imol;
                        r[spec] = pp;
		     }
		     catch (const std::runtime_error &rte) {
			// nothing too bad, just don't add that residue
			// to the plot
		     }
		  }
	       }
	    }
	 }
      }
   }
   return r;
}


#include "utils/coot-utils.hh"

void
gl_rama_plot_t::setup_buffers(float rama_plot_scale) {

   // I presume that we need to attach bufffers before this function is called?

   hud_tmesh_for_other_normal.setup_quad();
   hud_tmesh_for_other_outlier.setup_quad();
   hud_tmesh_for_pro_normal.setup_quad();
   hud_tmesh_for_pro_outlier.setup_quad();
   hud_tmesh_for_gly_normal.setup_quad();
   hud_tmesh_for_gly_outlier.setup_quad();

   // hud_tmesh_for_rama["other-normal"].setup_instancing_buffers()?
   hud_tmesh_for_other_normal.setup_instancing_buffers(10000);
   hud_tmesh_for_other_outlier.setup_instancing_buffers(10000);
   hud_tmesh_for_pro_normal.setup_instancing_buffers(10000);
   hud_tmesh_for_pro_outlier.setup_instancing_buffers(10000);
   hud_tmesh_for_gly_normal.setup_instancing_buffers(10000);
   hud_tmesh_for_gly_outlier.setup_instancing_buffers(10000);

   texture_for_other_normal.init( "rama-plot-other-normal.png");
   texture_for_other_outlier.init("rama-plot-other-outlier.png");
   texture_for_gly_normal.init(   "rama-plot-gly-normal.png");
   texture_for_gly_outlier.init(  "rama-plot-gly-outlier.png");
   texture_for_pro_normal.init(   "rama-plot-pro-normal.png");
   texture_for_pro_outlier.init(  "rama-plot-pro-outlier.png");

   // background image/texture
   hud_tmesh_for_global_distribution_non_gly_pro.setup_quad();
   texture_for_global_distribution_non_gly_pro.init("rama-standard-inverted.png");
   hud_tmesh_for_global_distribution_pro.setup_quad();
   texture_for_global_distribution_pro.init("test-image.png");
   hud_tmesh_for_global_distribution_gly.setup_quad();
   texture_for_global_distribution_gly.init("test-image.png");

   // add bar geometry to hud_mesh_for_axes_and_ticks

   hud_mesh_for_axes_and_ticks.set_name("hud_mesh_for_rama_axes_and_ticks");
   hud_mesh_for_axes_and_ticks.setup_simple_camera_facing_quad();
   hud_mesh_for_axes_and_ticks.setup_instancing_buffer(30, sizeof(HUD_bar_attribs_t));
   glm::vec4 col(0.8, 0.8, 0.8, 0.8); // white background plot will be different
   std::vector<HUD_bar_attribs_t> bars;
   float    thin = 0.008;
   float  v_thin = 0.007;
   float x_scale = rama_plot_scale;
   float y_scale = rama_plot_scale;
   bars.push_back(HUD_bar_attribs_t(col, glm::vec2(0.0,     0.0), x_scale, thin * y_scale));  // bottom bar
   bars.push_back(HUD_bar_attribs_t(col, glm::vec2(0.0, y_scale), x_scale, thin * y_scale));  // top bar
   bars.push_back(HUD_bar_attribs_t(col, glm::vec2(0.0,     0.0), thin * x_scale, y_scale));  // left bar
   bars.push_back(HUD_bar_attribs_t(col, glm::vec2(x_scale, 0.0), thin * x_scale, y_scale));  // right bar

   // I am not sure that I want ticks - I don't care about the values of phi,psi, just that they are
   // either good or outlier

   const unsigned int n_ticks = 6;
   for (unsigned int i=0; i<=n_ticks; i++) {
      float frac = static_cast<float>(i)/static_cast<float>(n_ticks);
      HUD_bar_attribs_t tick(col, glm::vec2(-0.04 * x_scale, frac * y_scale), 0.04 * x_scale, v_thin * y_scale);
      bars.push_back(tick);
   }
   for (unsigned int i=0; i<=n_ticks; i++) {
      float frac = static_cast<float>(i)/static_cast<float>(n_ticks);
      HUD_bar_attribs_t tick(col, glm::vec2(x_scale * frac, -0.04 * y_scale), v_thin * x_scale, 0.04 * y_scale);
      bars.push_back(tick);
   }
   // currently the graph is scaled to 40% with bottom left at the origin.
   for (auto &bar : bars) {
      bar.position_offset += glm::vec2(-0.9, -0.9); // move it to bottom left
   }

   hud_mesh_for_axes_and_ticks.update_instancing_buffer_data(bars);

}

std::pair<glm::vec2, glm::vec2>
gl_rama_plot_t::get_munged_offset_and_scale(screen_position_origins_t spo,
                                            const glm::vec2 &offset_natural,
                                            float scale_x_natural, float scale_y_natural,
                                            int glarea_width, int glarea_height) const {

   glm::vec2 offset_rel = glm::vec2(0,0);

   // glm::vec2 scales_new(scale_x_natural, scale_y_natural);

   // for 900 pixels and offset if 0.1 is 90 pixels.
   // 90 pixels in a 1000 pixels widths is 0.1/wr

   float wr = static_cast<float>(900)/static_cast<float>(glarea_width);
   float hr = static_cast<float>(900)/static_cast<float>(glarea_height);

   if (spo == TOP_LEFT)
      offset_rel = glm::vec2(-1.0 + offset_natural.x/wr, 1.0 + offset_natural.y/hr) - offset_natural;
   if (spo == BOTTOM_LEFT)
      offset_rel = glm::vec2(wr - 1.0, hr - 1.0) * offset_natural;
   if (spo == BOTTOM_RIGHT)
      offset_rel = glm::vec2(1.0 + offset_natural.x/wr, -1.0 + offset_natural.y/hr);
   if (spo == TOP_RIGHT)
      offset_rel = glm::vec2(1.0 + offset_natural.x/wr, 1.0 + offset_natural.y/hr);

   glm::vec2 scales_new(scale_x_natural * wr, scale_y_natural * hr);

   return std::pair<glm::vec2, glm::vec2>(offset_rel, scales_new);

}

   


void
gl_rama_plot_t::draw(Shader *shader_for_rama_plot_axes_and_ticks_p,
                     Shader *shader_for_rama_plot_phi_psis_markers_p,
                     Shader *shader_for_hud_image_textures_p,
                     int glarea_width, int glarea_height) {

   // draw() needs to:
   //
   //  1: draw the box outline
   //  2: draw the phi ticks
   //  3: draw the psi ticks
   //  4: draw the phi ticks
   //  5: draw the psi ticks
   //  6: draw/render the phi axis label
   //  7: draw/render the psi axis label
   //  8: draw/render the background texture
   //  9: draw/render mouse-over tooltips (just use text
   //                rendered to an opaque (black?) texture - it doesn't need a background)
   //                if I do that, how can I put a line around the outide? (that would be cool)
   // 10: draw the contour lines
   //
   //
   //
   // And how shall I draw the actual phi/psi points?
   // GLY, PRO, other?
   // With instanced textures I think. HUDTextureMesh_attribs_t gives us
   // vec2 positions - that's all we need
   // Use textures for Normal_other,  Normal_PRO,  Normal_GLY
   //                  Outlier_other, Outlier_PRO, Outlier_GLY
   // Path-tracker... track the path of the phi_psis with on-the-fly points
   // (that will need another shader, I think (as simple as can be)

   // Do I want to render this whole thing to a framebuffer?
   // No, but it's not much effort to do so if needed.
   // The scale and the position can be sent as a uniform.
   //
   // Turn off GL depth test and draw in order

   // Use shader_for_hud_geometry_tooltip_text (maybe rename this
   // "shader_for_hud_text" because it's general purpose)
   // for the tick labels and mouse-overp

   // Use a HUD mesh for each of the plot point types: normal, outlier, gly, pro, other
   // 6 draw calls - hmm.

   // shader_for_rama_plot_axes_and_ticks is like hud-bars.shader, but has some extra
   // scales and offsets to adjust the size and position of the plot (which "naturally"
   // is rendered at -1 to +1)

   // glDisable(GL_DEPTH_TEST); // interesting.
   // glDisable(GL_BLEND);

   glEnable(GL_BLEND);
   glm::vec2 offset_position_natural(0.1, -0.1);
   auto p_s = get_munged_offset_and_scale(BOTTOM_LEFT, offset_position_natural, 1.0, 1.0, glarea_width, glarea_height);
   glm::vec2 munged_position_offset = p_s.first;
   glm::vec2 munged_scales = p_s.second;

   hud_mesh_for_axes_and_ticks.set_scales(munged_scales);
   hud_mesh_for_axes_and_ticks.set_offset_positions(munged_position_offset * glm::vec2(10,10));
   hud_mesh_for_axes_and_ticks.draw(shader_for_rama_plot_axes_and_ticks_p); // draws instanced objects (bars)

   glEnable(GL_BLEND);
   glEnable(GL_DEPTH_TEST); // needed?
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   std::cout << "munged_scales " << glm::to_string(munged_scales) << " munged_offsets " << glm::to_string(munged_position_offset) << std::endl;

   GLenum err = glGetError(); if (err) std::cout << "GL ERROR:: gl_rama_plot_t::draw() A error " << err << std::endl;
   hud_tmesh_for_other_normal.set_window_resize_position_correction(munged_position_offset * glm::vec2(10,10));
   err = glGetError(); if (err) std::cout << "GL ERROR:: gl_rama_plot_t::draw() B error " << err << std::endl;
   hud_tmesh_for_other_normal.set_window_resize_scales_correction(munged_scales);
   err = glGetError(); if (err) std::cout << "GL ERROR:: gl_rama_plot_t::draw() C error " << err << std::endl;
   texture_for_other_normal.Bind(0);
   err = glGetError(); if (err) std::cout << "GL ERROR:: gl_rama_plot_t::draw() D error " << err << std::endl;
   hud_tmesh_for_other_normal.draw_instances(shader_for_rama_plot_phi_psis_markers_p);
   err = glGetError(); if (err) std::cout << "GL ERROR:: gl_rama_plot_t::draw() E error " << err << std::endl;

   texture_for_pro_normal.Bind(0);
   hud_tmesh_for_pro_normal.set_window_resize_position_correction(munged_position_offset * glm::vec2(10,10));
   hud_tmesh_for_pro_normal.set_window_resize_scales_correction(munged_scales);
   hud_tmesh_for_pro_normal.draw_instances(shader_for_rama_plot_phi_psis_markers_p);
   err = glGetError(); if (err) std::cout << "GL ERROR:: gl_rama_plot_t::draw() F error " << err << std::endl;

   texture_for_gly_normal.Bind(0);
   hud_tmesh_for_gly_normal.set_window_resize_position_correction(munged_position_offset * glm::vec2(10,10));
   hud_tmesh_for_gly_normal.set_window_resize_scales_correction(munged_scales);
   hud_tmesh_for_gly_normal.draw_instances(shader_for_rama_plot_phi_psis_markers_p);
   err = glGetError(); if (err) std::cout << "GL ERROR:: gl_rama_plot_t::draw() G error " << err << std::endl;

   texture_for_other_outlier.Bind(0);
   hud_tmesh_for_other_outlier.set_window_resize_position_correction(munged_position_offset * glm::vec2(10,10));
   hud_tmesh_for_other_outlier.set_window_resize_scales_correction(munged_scales);
   hud_tmesh_for_other_outlier.draw_instances(shader_for_rama_plot_phi_psis_markers_p);
   err = glGetError(); if (err) std::cout << "GL ERROR:: gl_rama_plot_t::draw() H error " << err << std::endl;

   texture_for_pro_outlier.Bind(0);
   hud_tmesh_for_pro_outlier.set_window_resize_position_correction(munged_position_offset * glm::vec2(10,10));
   hud_tmesh_for_pro_outlier.set_window_resize_scales_correction(munged_scales);
   hud_tmesh_for_pro_outlier.draw_instances(shader_for_rama_plot_phi_psis_markers_p);
   err = glGetError(); if (err) std::cout << "GL ERROR:: gl_rama_plot_t::draw() I error " << err << std::endl;

   texture_for_gly_outlier.Bind(0);
   hud_tmesh_for_gly_outlier.set_window_resize_position_correction(munged_position_offset * glm::vec2(10,10));
   hud_tmesh_for_gly_outlier.set_window_resize_scales_correction(munged_scales);
   hud_tmesh_for_gly_outlier.draw_instances(shader_for_rama_plot_phi_psis_markers_p);
   err = glGetError(); if (err) std::cout << "GL ERROR:: gl_rama_plot_t::draw() J error " << err << std::endl;

   if (true) {
      texture_for_global_distribution_non_gly_pro.Bind(0);
      // munged_position_offset = glm::vec2(0,0);
      // munged_scales = glm::vec2(1,1);
      hud_tmesh_for_global_distribution_non_gly_pro.set_scales(glm::vec2(0.25, 0.25));
      hud_tmesh_for_global_distribution_non_gly_pro.set_position(glm::vec2(-0.65, -0.65));
      hud_tmesh_for_global_distribution_non_gly_pro.set_window_resize_position_correction(munged_position_offset * glm::vec2(10,10));
      hud_tmesh_for_global_distribution_non_gly_pro.set_window_resize_scales_correction(munged_scales);
      hud_tmesh_for_global_distribution_non_gly_pro.draw(shader_for_hud_image_textures_p);
   }
}
