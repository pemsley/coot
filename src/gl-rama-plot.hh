
#ifndef GL_RAMA_PLOT_HH
#define GL_RAMA_PLOT_HH

#include <clipper/core/ramachandran.h>
#include "rama-plot-phi-psi.hh"
#include "HUDMesh.hh"
#include "HUDTextureMesh.hh"
#include "Texture.hh"
#include "Shader.hh"

class gl_rama_plot_t {

   class canvas_tick_t {
   public:
      double x, y;
      short int axis;
      // i==0 for x axis, i==1 for y axis.
      canvas_tick_t(int i, double a, double b) {
         x = a;
         y = b;
         axis = i;
      }
      double start_x() {return x;};
      double start_y() {return y;};
      double end_x() { if (axis == 0) return x; else return x - 10;};
      double end_y() { if (axis == 1) return y; else return y + 10;}; 
   };
   clipper::Ramachandran rama_all;
   clipper::Ramachandran rama_gly;
   clipper::Ramachandran rama_pro;
   clipper::Ramachandran rama_non_gly_pro;
   clipper::Ramachandran rama_ile_val;
   clipper::Ramachandran rama_pre_pro;
   clipper::Ramachandran rama_non_gly_pro_ile_val_or_pre_pro;
   double rama_threshold_preferred;
   double rama_threshold_allowed;
   std::string colour_for_OK_phi_psi;
   std::string colour_for_OK_phi_psi_pro;
   std::string colour_for_outlier_phi_psi;
   std::map<coot::residue_spec_t, rama_plot::phi_psi_t> phi_psis;
   clipper::Ramachandran::TYPE current_background_type;
   std::map<coot::residue_spec_t, rama_plot::phi_psi_t> generate_pseudo_phi_psis();
   std::map<coot::residue_spec_t, rama_plot::phi_psi_t> generate_phi_psis(int imol, mmdb::Manager *mol_in);
   bool draw_outliers_only_flag; 
   void init(); // rama things
   HUDTextureMesh hud_tmesh_for_other_normal;
   HUDTextureMesh hud_tmesh_for_other_outlier;
   HUDTextureMesh hud_tmesh_for_pro_normal;
   HUDTextureMesh hud_tmesh_for_pro_outlier;
   HUDTextureMesh hud_tmesh_for_gly_normal;
   HUDTextureMesh hud_tmesh_for_gly_outlier;
   Texture texture_for_other_normal;
   Texture texture_for_other_outlier;
   Texture texture_for_pro_normal;
   Texture texture_for_pro_outlier;
   Texture texture_for_gly_normal;
   Texture texture_for_gly_outlier;
   HUDTextureMesh hud_tmesh_for_global_distribution_non_gly_pro;
   HUDTextureMesh hud_tmesh_for_global_distribution_pro;
   HUDTextureMesh hud_tmesh_for_global_distribution_gly;
   Texture texture_for_global_distribution_non_gly_pro;
   Texture texture_for_global_distribution_gly;
   Texture texture_for_global_distribution_pro;
   HUDMesh hud_mesh_for_axes_and_ticks;

   // probably these should be in a lower-lever library.
   enum screen_position_origins_t { TOP_LEFT, TOP_RIGHT, BOTTOM_LEFT, BOTTOM_RIGHT};

   std::pair<glm::vec2, glm::vec2>
    get_munged_offset_and_scale(screen_position_origins_t spo,
                                const glm::vec2 &offset_natural,
                                float scale_x_natural, float scale_y_natural,
                                int glarea_width, int glarea_height) const;
   
   void update_hud_tmeshes(const std::map<coot::residue_spec_t, rama_plot::phi_psi_t> &phi_psi_map);

public:
   gl_rama_plot_t() { init(); }
   void setup_buffers(float rama_plot_scale); // setup OpenGL things - must be done after OpenGL realize()
   void setup_from(int imol, mmdb::Manager *mol);
   void update_phi_psis_on_moved_atoms();
   // void background_to_type(GtkWidget *canvas, clipper::Ramachandran::TYPE); // don't change it if we are already there of course.
   void draw(Shader *shader_for_axes_and_tick,
             Shader *shader_for_rama_plot_phi_psis_markers_p,  // instanced
             Shader *shader_for_hud_textures_p,
             int glarea_width, int glarea_height);
   float position_hash; // updated and tested in setup_from() so that we don't recalculate if we don't need to
};

#endif // GL_RAMA_PLOT_HH
