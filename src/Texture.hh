
#ifndef HMT_TEXTURE_HH
#define HMT_TEXTURE_HH

#include <string>
#include <vector>
#include <epoxy/gl.h>

class Texture {

   std::string default_directory;
   unsigned int id; // not used
   int image_width;
   int image_height;

public:
   enum type_t { DIFFUSE, NORMAL, SPECULAR, ROUGHNESS, AMBIENT_OCCLUSION, AMBIENT_OCCLUSION_ROUGHNESS_METALICITY };
   Texture() {}
   explicit Texture(const std::string &file_name, type_t t, bool reversed_normals=false);
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

