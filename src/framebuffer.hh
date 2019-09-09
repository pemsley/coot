
#ifndef FRAMEBUFFER_HH
#define FRAMEBUFFER_HH

#include <tuple>

std::tuple<unsigned int, unsigned int, unsigned int>
setup_frame_buffer_old(int width, int height);


class framebuffer {
   unsigned int fbo;
   unsigned int texture_colour;
   unsigned int texture_depth;
   std::vector<GLenum> drawbuffer; // texture attachments
   bool filled;
   void tear_down();

public:
   framebuffer();

   ~framebuffer(); // delete fbo and the textures

   std::string name;

   void init(int width, int height, unsigned int attachment_index_color_texture, const std::string &name_in="");

   void generate_colourtexture(unsigned int width, unsigned int height);
   void generate_depthtexture( unsigned int width, unsigned int height);
   void generate_framebuffer_object(unsigned int width, unsigned int height, unsigned int attachment_index_color_texture);

   void resize(int width, int height, unsigned int attachment_index_color_texture) { init(width, height, attachment_index_color_texture); }

   void bind();

   unsigned int get_texture_colour() const { return texture_colour; }
   unsigned int get_texture_depth()  const { return texture_depth; }


};

#endif // FRAMEBUFFER_HH
