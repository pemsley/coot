
#ifndef HMT_TEXTURE_HH
#define HMT_TEXTURE_HH

#include <string>
#include <vector>
#include <epoxy/gl.h>

#include <glm/glm.hpp>

class Texture {

   std::string default_directory;
   unsigned int id; // not used
   int image_width;
   int image_height;
   void colour_bar(int image_width_in, int image_height_in, const std::vector<glm::vec4> &colours, unsigned int n_ticks);
   void add_tick_marks(unsigned int n_ticks, const glm::vec4 &tick_colour, unsigned char *image_data); // n_ticks should be at least 2

public:
   enum type_t { DIFFUSE, NORMAL, SPECULAR, ROUGHNESS, AMBIENT_OCCLUSION, AMBIENT_OCCLUSION_ROUGHNESS_METALICITY };
   Texture() {}
   explicit Texture(const std::string &file_name, type_t t, bool reversed_normals=false);
   // colour components are in the range 0 to 1.
   Texture(int image_width, int image_height, glm::vec4 colour); // create a solid colour texture
   Texture(int image_width, int image_height, const std::vector<glm::vec4> &colours); // colour bar
   Texture(int image_width, int image_height, const std::vector<glm::vec4> &colours, unsigned int n_ticks); // colour bar with n_ticks tick marks. n_ticks should be at least 2
   ~Texture(); // don't close
   type_t type;
   bool reversed_normals;
   GLuint m_texture_handle; // make this private after this testing              
   std::string file_name;
   // std::string type; 20211121-PE now we make it a enum
   void init(const std::string &file_name);
   void init(const std::string &local_file_name, const std::string &directory);
   void set_file_name(const std::string &fn) { file_name = fn; }
   void set_type(type_t t) { type = t; }
   void Bind(unsigned int unit);
   void set_default_directory(const std::string &dir);
   void close();
   std::pair<int, int> get_image_width_and_height() const;
   void handle_raw_image_data(const std::string &image_name, const std::vector<unsigned char> &image_data, int width, int height);

};

#endif // HMT_TEXTURE_HH

