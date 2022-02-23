
#ifndef FRAMEBUFFER_HH
#define FRAMEBUFFER_HH

#include <vector>
#include <tuple>

#include <gtk/gtk.h>
#include <epoxy/gl.h>

std::tuple<unsigned int, unsigned int, unsigned int>
setup_frame_buffer_old(int width, int height);


class framebuffer {
   unsigned int fbo;
   unsigned int colour_texture;
   unsigned int depth_texture;
   std::vector<GLenum> drawbuffer; // texture attachments
   bool filled;

public:
   framebuffer();

   ~framebuffer(); // delete fbo and the textures

   std::string name;

   void init(int width, int height, unsigned int attachment_index_color_texture, const std::string &name_in="");
   void reset(int width, int height);
   void reset_test(int width, int height);
   void tear_down();

   void generate_colourtexture(unsigned int width, unsigned int height);
   void generate_depthtexture( unsigned int width, unsigned int height);
   void generate_framebuffer_object(unsigned int width, unsigned int height, unsigned int attachment_index_color_texture);
   void generate_framebuffer_object_test(unsigned int width, unsigned int height, unsigned int attachment_index_color_texture);

   void resize(int width, int height, unsigned int attachment_index_color_texture) { init(width, height, attachment_index_color_texture); }

   void bind();

   unsigned int get_texture_colour() const { return colour_texture; }
   unsigned int get_texture_depth()  const { return depth_texture; }

   unsigned int get_fbo() const { return fbo; }

   unsigned int gPosition;
   unsigned int gNormal;

   void do_gbuffer_stuff(int w, int h);

};

#endif // FRAMEBUFFER_HH
