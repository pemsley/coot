
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
         phi_psi_map = generate_phi_psis(imol, mol);
         update_hud_tmeshes(phi_psi_map); // no need for attach_buffers() as this is instanced data.
         position_hash = position_hash_now;
      }
   }
   // auto tp_1 = std::chrono::high_resolution_clock::now();
   // auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
   // std::cout << "INFO:: setup_from() " << d10 << " milliseconds" << std::endl;
}

void
gl_rama_plot_t::clear() {

   phi_psi_map.clear();
   std::map<coot::residue_spec_t, rama_plot::phi_psi_t> empty;
   update_hud_tmeshes(empty);

   // I don't want to delete the VAO or the buffers for the instance object.

}

bool
gl_rama_plot_t::is_active() const {

   return (! phi_psi_map.empty());

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

   // float rama_plot_scale = 0.7; // pass this or make it a member item 20210908-PE  now it's a data member.
   float sf = rama_plot_scale * 0.23; // 0.25 is too big, 0.23 is about right, 0.22 is a bit too small
   std::vector<glm::vec2> new_other_normal_positions;
   std::vector<glm::vec2> new_other_outlier_positions;
   std::vector<glm::vec2> new_pro_normal_positions;
   std::vector<glm::vec2> new_pro_outlier_positions;
   std::vector<glm::vec2> new_gly_normal_positions;
   std::vector<glm::vec2> new_gly_outlier_positions;
   std::map<coot::residue_spec_t, rama_plot::phi_psi_t>::const_iterator it;
   for (it=phi_psi_map.begin(); it!=phi_psi_map.end(); ++it) {
      const auto &phi_psi = it->second; // phi and psi are in degrees
      double phi_r = clipper::Util::d2rad(phi_psi.phi);
      double psi_r = clipper::Util::d2rad(phi_psi.psi);
      glm::vec2 pos(sf * phi_psi.phi, sf * phi_psi.psi);
      // std::cout << "phi_psi pos: " << glm::to_string(pos) << std::endl;
      if (phi_psi.residue_name == "PRO") {
         double probability = rama_pro.probability(phi_r, psi_r);
         if (probability < rama_threshold_allowed) {
            new_pro_outlier_positions.push_back(pos);
         } else {
            new_pro_normal_positions.push_back(pos);
         }
      } else {
         if (phi_psi.residue_name == "GLY") {
            double probability = rama_gly.probability(phi_r, psi_r);
            if (probability < rama_threshold_allowed) {
               new_gly_outlier_positions.push_back(pos);
            } else {
               new_gly_normal_positions.push_back(pos);
            }
         } else {
            double probability = rama_non_gly_pro.probability(phi_r, psi_r);
            if (probability < rama_threshold_allowed) {
               new_other_outlier_positions.push_back(pos);
            } else {
               new_other_normal_positions.push_back(pos);
            }
         }
      }
   }

   if (false)
      std::cout << "debug:: update_hud_tmeshes() counts "
                << " " << new_other_normal_positions.size() << " " << new_other_outlier_positions.size()
                << " " << new_pro_normal_positions.size() << " " << new_pro_outlier_positions.size()
                << " " << new_gly_normal_positions.size() << " " << new_gly_outlier_positions.size()
                << std::endl;

   // mark the corners of the plot for scaling and offset testing

   if (false) {
      new_other_normal_positions.push_back(glm::vec2(sf * -180.0, sf * -180.0));
      new_other_normal_positions.push_back(glm::vec2(sf *  180.0, sf * -180.0));
      new_other_normal_positions.push_back(glm::vec2(sf * -180.0, sf *  180.0));
      new_other_normal_positions.push_back(glm::vec2(sf *  180.0, sf *  180.0));
   }

   if (false) // check the z-depth and/or draw order
      for (unsigned int i=0; i<new_other_normal_positions.size(); i++)
         std::cout << "new_other_normal_positions " << i << " " << glm::to_string(new_other_normal_positions[i])
                   << std::endl;


   // for the background
   // float ff = -0.5 * rama_plot_scale + 0.9;
   // hud_tmesh_for_global_distribution_non_gly_pro.set_position(glm::vec2(-ff, -ff));

   glm::vec2 plot_point_scales(0.012, 0.012); // this gets the marker point size right

   float ff =  -0.5 * rama_plot_scale + 0.9;

   glm::vec2 offset(-ff, -ff);

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

   if (false) {
      // add in some fake points for testing
      coot::residue_spec_t f1("Test", 1, "");
      coot::residue_spec_t f2("Test", 2, "");
      coot::residue_spec_t f3("Test", 3, "");
      coot::residue_spec_t f4("Test", 4, "");

      r[f1] = rama_plot::phi_psi_t(-180, -180, "Test-Res", "TR1", 1, "", "A", false);
      r[f2] = rama_plot::phi_psi_t( 180, -180, "Test-Res", "TR2", 2, "", "A", false);
      r[f3] = rama_plot::phi_psi_t( 180,  180, "Test-Res", "TR3", 3, "", "A", false);
      r[f4] = rama_plot::phi_psi_t(-180,  180, "Test-Res", "TR4", 4, "", "A", false);
   }

   return r;
}


#include "utils/coot-utils.hh"

void
gl_rama_plot_t::setup_buffers(float rama_plot_scale_in) {

   // I presume that we need to attach bufffers before this function is called?

   rama_plot_scale = rama_plot_scale_in;

   hud_tmesh_for_other_normal.set_name("hud_tmesh_for_other_normal");
   hud_tmesh_for_other_outlier.set_name("hud_tmesh_for_other_outlier");
   hud_tmesh_for_pro_normal.set_name("hud_tmesh_for_pro_normal");
   hud_tmesh_for_pro_outlier.set_name("hud_tmesh_for_pro_outlier");
   hud_tmesh_for_gly_normal.set_name("hud_tmesh_for_gly_normal");
   hud_tmesh_for_gly_outlier.set_name("hud_tmesh_for_gly_outlier");

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
   texture_for_global_distribution_pro.init("rama-standard-inverted.png"); // needs fixing
   hud_tmesh_for_global_distribution_gly.setup_quad();
   texture_for_global_distribution_gly.init("rama-standard-inverted.png"); // needs fixing

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

   float wr = static_cast<float>(700)/static_cast<float>(glarea_width);  // 20220328-PE was 900, 900
   float hr = static_cast<float>(700)/static_cast<float>(glarea_height);

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
                     int glarea_width_at_hud_start,
                     int glarea_heigth_at_hud_start,
                     int glarea_width, int glarea_height) {

   // draw() needs to:
   //
   //  1: draw the box outline
   //  2: draw the phi ticks
   //  3: draw the psi ticks
   //  4: draw the phi tick labels
   //  5: draw the psi tick labels
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

   glDisable(GL_BLEND);

   glm::vec2 offset_position_natural(0.1, 0.1);
   auto p_s = get_munged_offset_and_scale(BOTTOM_LEFT, offset_position_natural, 1.0, 1.0, glarea_width, glarea_height);
   glm::vec2 munged_position_offset = p_s.first;
   glm::vec2 munged_scales = p_s.second;

   if (true) {
      texture_for_global_distribution_non_gly_pro.Bind(0);
      // munged_position_offset = glm::vec2(0,0);
      // munged_scales = glm::vec2(1,1);

      // In on_glarea_realize() this is setup using:
      // g.gl_rama_plot.setup_buffers(0.6)

      // 20220329-PE use a reference to graphics_info_t::hud_tmesh_for_global_distribution, which can be reassigned,
      // depending on the moused-over residue type (that should happen in the moused-over-residue gtk callback).
      //
      hud_tmesh_for_global_distribution_non_gly_pro.set_scales(glm::vec2(rama_plot_scale * 0.5, rama_plot_scale * 0.5));
      float ff = -0.5 * rama_plot_scale + 0.9;
      hud_tmesh_for_global_distribution_non_gly_pro.set_position(glm::vec2(-ff, -ff));
      hud_tmesh_for_global_distribution_non_gly_pro.set_window_resize_position_correction(munged_position_offset * glm::vec2(10.0,10.0));
      hud_tmesh_for_global_distribution_non_gly_pro.set_window_resize_scales_correction(munged_scales);
      hud_tmesh_for_global_distribution_non_gly_pro.draw(shader_for_hud_image_textures_p);
   }

   glEnable(GL_BLEND);

   hud_mesh_for_axes_and_ticks.set_scales(munged_scales);
   hud_mesh_for_axes_and_ticks.set_offset_positions(munged_position_offset * glm::vec2(10,10));
   hud_mesh_for_axes_and_ticks.draw(shader_for_rama_plot_axes_and_ticks_p); // draws instanced objects (bars)

   glEnable(GL_BLEND);
   glEnable(GL_DEPTH_TEST); // needed?
   glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

   // std::cout << "munged_scales " << glm::to_string(munged_scales) << " munged_offsets "
   // << glm::to_string(munged_position_offset) << std::endl;
   GLenum err;

   err = glGetError(); if (err) std::cout << "GL ERROR:: gl_rama_plot_t::draw() A error " << err << std::endl;
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

   glDisable(GL_BLEND);
   // glDisable(GL_DEPTH_TEST);

}


// std::pair<bool, coot::residue_spec_t>
mouse_over_hit_t
gl_rama_plot_t::get_mouse_over_hit(double x_widget, double y_widget, int widget_width, int widget_height) const {

   bool residue_hit_status = false;
   bool plot_hit_status = false;
   coot::residue_spec_t hit_residue_spec;

   auto mouse_position_to_phi_psi = [x_widget,y_widget, widget_width, widget_height] (float rama_plot_scale,
                                                                                      const glm::vec2 &scales,
                                                                                      const glm::vec2 &position,
                                                                                      const glm::vec2 &window_resize_scales_correction,
                                                                                      const glm::vec2 &window_resize_position_correction) {

                                       float x_opengl =  2.0 * x_widget/static_cast<float>(widget_width)  - 1.0;
                                       float y_opengl = -2.0 * y_widget/static_cast<float>(widget_height) + 1.0;

                                       // reverse the position
                                       float sf = rama_plot_scale * 0.23; // as in update_hud_tmeshes()
                                       glm::vec2 oowrsc_v2(1.0/window_resize_scales_correction.x, 1.0/window_resize_scales_correction.y);
                                       // the sign of these offset will depend on the corner that they're offset from
                                       // this is lower-lefft of course.
                                       glm::vec2 wrpc(window_resize_position_correction.x, -window_resize_position_correction.y);
                                       glm::vec2 rpp_1 = glm::vec2(x_opengl, y_opengl) - wrpc; // window_resize_position_correction;
                                       glm::vec2 rpp_2 = rpp_1 * oowrsc_v2;
                                       glm::vec2 rpp_3 = rpp_2 - position;
                                       glm::vec2 rpp_4(rpp_3.x/scales.x, rpp_3.y/scales.y);
                                       glm::vec2 rpp_5(rpp_4.x/sf, rpp_4.y/sf);

                                       glm::vec2 rpp_1a = glm::vec2(x_opengl, y_opengl);
                                       glm::vec2 rpp_2a = rpp_1a * oowrsc_v2;
                                       glm::vec2 rpp_3a = rpp_2a - position;
                                       glm::vec2 rpp_4a(rpp_3a.x/scales.x, rpp_3a.y/scales.y);
                                       glm::vec2 rpp_5a(rpp_4a.x/sf, rpp_4a.y/sf);

                                       if (false)
                                          std::cout << "debug:: window-pos-corr " << glm::to_string(window_resize_position_correction)
                                                    << " in " << x_widget << " " << y_widget
                                                    << " opengl " << x_opengl << " " << y_opengl << " "
                                                    << "phi-psi-sans-win-res-pos-cor: " << rpp_5a.x << " " << rpp_5a.y << " "
                                                    << "cursor-phi-psi: " << rpp_5.x << " " << rpp_5.y <<std::endl;

                                       rama_plot::phi_psi_t phi_psi_mouse_pos(rpp_5.x, rpp_5.y);

                                       return phi_psi_mouse_pos;
                                    };

   if (! is_active()) {
      return mouse_over_hit_t(false, false, hit_residue_spec); // no.
   } else {

      // I want scales and position and window_resize_position_correction window_resize_scales_correction

      int glarea_width = widget_width;
      int glarea_height = widget_height;
      glm::vec2 offset_position_natural(0.1, -0.1);
      auto p_s = get_munged_offset_and_scale(BOTTOM_LEFT, offset_position_natural, 1.0, 1.0, glarea_width, glarea_height);
      glm::vec2 munged_position_offset = p_s.first;
      glm::vec2 munged_scales = p_s.second;
      glm::vec2 window_resize_scales_correction = munged_scales;
      glm::vec2 window_resize_position_correction = munged_position_offset * glm::vec2(10,10);
      glm::vec2 plot_point_scales(0.012, 0.012);
      glm::vec2 scales = plot_point_scales;
      float ff =  -0.5 * rama_plot_scale + 0.9;
      glm::vec2 offset(-ff, -ff);
      glm::vec2 position = offset;

      if (false)
         std::cout << "get_mouse_over_hit()"
                   << " scales " << glm::to_string(scales)
                   << " position " << glm::to_string(position)
                   << " window_resize_scales_correction "   << glm::to_string(window_resize_scales_correction)
                   << " window_resize_position_correction " << glm::to_string(window_resize_position_correction)
                   << std::endl;

      rama_plot::phi_psi_t phi_psi_mouse = mouse_position_to_phi_psi(rama_plot_scale, scales, position,
                                                                     window_resize_scales_correction,
                                                                     window_resize_position_correction);

      if (phi_psi_mouse.phi >= -180.5)
         if (phi_psi_mouse.phi <= 180.5)
            if (phi_psi_mouse.psi >= -180.5)
               if (phi_psi_mouse.psi <= 180.5)
                  plot_hit_status = true;

      // std::cout << "debug:: phi_psi_mouse: in " << x_widget << " " << y_widget
      // << " out: " << phi_psi_mouse.phi << " " << phi_psi_mouse.psi << std::endl;
      std::map<coot::residue_spec_t, rama_plot::phi_psi_t>::const_iterator it;
      float delta_r_sqrd_best = 99999999.9;
      bool best_has_been_found = false;
      rama_plot::phi_psi_t phi_psi_best;
      for (it=phi_psi_map.begin(); it!=phi_psi_map.end(); ++it) {
         const auto &phi_psi = it->second; // phi and psi are in degrees
         double delta_phi = fabs(phi_psi.phi - phi_psi_mouse.phi);
         double delta_psi = fabs(phi_psi.psi - phi_psi_mouse.psi);
         if (delta_phi < 10.0) {
            if (delta_psi < 10.0) {
               double delta_r_sqrd = delta_phi * delta_phi + delta_psi * delta_psi;
               if (delta_r_sqrd < delta_r_sqrd_best) {
                  phi_psi_best = phi_psi;
                  best_has_been_found = true;
               }
            }
         }
      }
      if (best_has_been_found) {
         if (phi_psi_best.residue_this) {
            hit_residue_spec = coot::residue_spec_t(phi_psi_best.residue_this);
            residue_hit_status = true;
         }
      }
   }
   //std::cout << "get_mouse_over_hit() returning " << status << " " << hit_residue_spec << std::endl;
   return mouse_over_hit_t(residue_hit_status, plot_hit_status, hit_residue_spec);
}
