
#ifndef HMT_TEXTURE_HH
#define HMT_TEXTURE_HH

#include <string>
#include <epoxy/gl.h>

class Texture {

   GLuint m_texture_handle;

public:
   std::string id;
   Texture() {}
   Texture(const std::string &file_name);
   void init(const std::string &file_name);
   ~Texture(); // don't close
   void Bind(unsigned int unit);
   void close();

};

#endif // HMT_TEXTURE_HH

