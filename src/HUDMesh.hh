
#ifndef HUD_MESH_HH
#define HUD_MESH_HH

#include <string>
#include <vector>
#include <epoxy/gl.h>
#include <glm/glm.hpp>
#include "g_triangle.hh"
#include "Shader.hh"

// this is a bar_info
class HUD_bar_attribs_t {
public:
   glm::vec4 colour;
   glm::vec2 position_offset;
   // if this_is_a_non_moving_atoms_residue then we want to make 2 bars, an 80% thickness bar
   // and a small "grey" bar below. same x coords.
   float scale_x;
   float scale_y;
   HUD_bar_attribs_t(const glm::vec4 &c, const glm::vec2 &p, const float &l) :
      colour(c), position_offset(p), scale_x(l) {
      scale_y = 1.0;
   }
   HUD_bar_attribs_t(const glm::vec4 &c, const glm::vec2 &p, const float &ls, const float &ws) :
      colour(c), position_offset(p), scale_x(ls), scale_y(ws) {
   }
};

class HUD_button_limits_t {
public:
   float top;
   float bottom;
   float left;
   float right;
   HUD_button_limits_t(float t, float b, float l, float r) : top(t), bottom(b), left(l), right(r) {}
   bool is_hit(double x, double y) const {
      if (y >= bottom)
         if (y <= top)
            if (x >= left)
               if (x <= right)
                  return true;
      return false;
   }
};

// this is another type of bar_info
class HUD_button_info_t : public HUD_bar_attribs_t {
   static glm::vec2 calculate_position_offset(unsigned int button_index, int width, int height);
public:
   glm::vec4 colour_basic;
   glm::vec4 colour_pressed;
   glm::vec4 colour_highlighted;
   unsigned int position_offset_index; // used to determine the position of this bar when checking to see if this
                                       // button has been moused over

   bool is_drawn_flag;
   std::string button_label;
   static constexpr float button_width  = 0.30;
   static constexpr float button_height = 0.06;
   static glm::vec4 get_default_colour() { return glm::vec4(0.45f, 0.45f, 0.45f, 0.6f); }
   enum colour_highlight_type { BASIC, HIGHLIGHTED, PRESSED };
   bool (*callback_function)();
   void set_button_colour_for_mode(colour_highlight_type cht) {
      if (cht == BASIC)       colour = colour_basic;
      if (cht == HIGHLIGHTED) colour = colour_highlighted;
      if (cht == PRESSED)     colour = colour_pressed;
   }
   // we need the scale the button so that the button contains "Backrub Rotamer" at 900 pixels - this is
   // a terrible way of doing it.
   HUD_button_info_t() : HUD_bar_attribs_t(get_default_colour(), glm::vec2(0.0, 0.0), button_width, button_height),
                         colour_basic(get_default_colour()) {
      set_colours_from_basic();
      is_drawn_flag = true;
      callback_function = 0;
      position_offset_index = 0;
   }
   explicit HUD_button_info_t (const std::string &label) :
      HUD_bar_attribs_t(get_default_colour(), glm::vec2(0.0, 0.0), button_width, button_height),
      colour_basic(get_default_colour()), button_label(label) {
      set_colours_from_basic();
      is_drawn_flag = true;
      callback_function = 0;
      position_offset_index = 0;
   }
   void set_colours_from_basic() {
      colour_highlighted = colour_basic + glm::vec4( 0.07,  0.07,  0.07, colour_basic.w);
      colour_pressed     = colour_basic + glm::vec4(-0.2, -0.2, -0.2, 0.6);
   }
   void set_colour(const glm::vec4 &c) {
      colour_basic = c;
      colour = colour_basic;
      set_colours_from_basic();
   }
   void set_position_offset(unsigned int poi, const glm::vec2 &p) {
      position_offset_index = poi;
      position_offset = p;
   }
   void set_draw_flag(bool state) {
      is_drawn_flag = state;
   }
   void set_label(const std::string &l) {
      button_label = l;
   }
   void connect(bool (*cb_in)()) {
      callback_function = cb_in;
   }
   // counting from the bottom! (at the moment)
   void set_scales_and_position_offset(unsigned int button_index, int width, int height);
   // This doesn't change the button width
   void set_position_offset(unsigned int button_index, int width, int height);
   HUD_button_limits_t get_button_limits(int width, int height) const;
};

class HUDMesh {
   void setup_buffers();
   void init();
   bool first_time;
   unsigned int max_n_instances;
   unsigned int n_instances;
   unsigned int inst_hud_bar_attribs_buffer_id;
   bool use_shading_flag = true;
   bool scales_have_been_set;
   bool offset_position_has_been_set;
   glm::vec2 scales;
   glm::vec2 offset_position;
   glm::vec2 window_resize_scales_correction;
   glm::vec2 window_resize_position_correction;

public:
   GLuint vao;
   GLuint vertex_buffer_id;
   GLuint shades_buffer_id;
   GLuint index_buffer_id;
   bool this_mesh_is_closed;
   bool use_blending;
   std::vector<glm::vec2> vertices;
   std::vector<float> shades;
   std::vector<g_triangle> triangles;
   std::string name;
   HUDMesh() { init(); }
   HUDMesh(const std::string &n) : name(n) { init(); }
   void setup_simple_camera_facing_quad();
   void setup_camera_facing_quad_for_bar();
   void setup_vertices_and_triangles_for_button();
   void setup_vertices_and_triangles_for_tooltip_background();
   void set_name(const std::string &n) { name = n; }
   void set_scales(const glm::vec2 &s) { scales = s; scales_have_been_set = true; }
   void set_offset_positions(const glm::vec2 &p) { offset_position = p; offset_position_has_been_set = true; }
   void set_use_shading(bool state) {
      use_shading_flag = state;
   }
   // size_of_bar_info can be sizeof(HUD_bar_attribs_t) or sizeof(HUD_button_info_t).
   // all bar types should be derived from HUD_bar_attribs_t because the layout in
   // the shader data reflects the HUD_bar_attribs_t.
   void setup_instancing_buffer(unsigned int n_boxes, unsigned int size_of_bar_info);
   void update_instancing_buffer_data(const std::vector<HUD_bar_attribs_t> &new_bars);
   void update_instancing_buffer_data(const std::vector<HUD_button_info_t> &new_buttons);
   void set_window_resize_scales_correction(const glm::vec2 &v) {
      window_resize_scales_correction = v;
   }
   void set_window_resize_position_correction(const glm::vec2 &v) {
      window_resize_position_correction = v;
   }

   void draw(Shader *shader);
   void close() { this_mesh_is_closed = true; }

};

#endif // HUD_MESH_HH
