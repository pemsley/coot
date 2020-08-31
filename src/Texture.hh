
#ifndef HMT_TEXTURE_HH
#define HMT_TEXTURE_HH

#include <string>
#include <epoxy/gl.h>

class Texture {

   GLuint m_texture_handle;
   std::string default_directory;

public:
   std::string id;
   Texture() {}
   Texture(const std::string &file_name);
   void init(const std::string &file_name);
   ~Texture(); // don't close
   void Bind(unsigned int unit);
   void set_default_directory(const std::string &dir);
   void close();

};

#endif // HMT_TEXTURE_HH

