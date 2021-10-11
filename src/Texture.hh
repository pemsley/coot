
#ifndef HMT_TEXTURE_HH
#define HMT_TEXTURE_HH

#include <string>
#include <vector>
#include <epoxy/gl.h>

class Texture {

   GLuint m_texture_handle;
   std::string default_directory;
   std::string file_name;
   std::string type;
   unsigned int id;
   int image_width;
   int image_height;

public:
   Texture() {}
   explicit Texture(const std::string &file_name);
   ~Texture(); // don't close

   void init(const std::string &file_name);
   void init(const std::string &local_file_name, const std::string &directory);
   void set_file_name(const std::string &fn) { file_name = fn; }
   void set_type(const std::string &t) { type = t; }
   void Bind(unsigned int unit);
   void set_default_directory(const std::string &dir);
   void close();
   std::pair<int, int> get_image_width_and_height() const;
   void handle_raw_image_data(const std::string &image_name, const std::vector<unsigned char> &image_data, int width, int height);

};

#endif // HMT_TEXTURE_HH

