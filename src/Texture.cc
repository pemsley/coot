
#include <iostream>

#include "Texture.hh"
#include "stb_image.h"

#ifdef THIS_IS_HMT
#else
#include "utils/coot-utils.hh"
#endif

Texture::Texture(const std::string &file_name) {
   init(file_name);
}

void
Texture::init(const std::string &file_name_in) {

   file_name = file_name_in;

#ifdef THIS_IS_HMT
#else
   std::string pkg_data_dir = coot::package_data_dir();
   std::string default_directory = pkg_data_dir + "/textures";

   if (! coot::file_exists(file_name)) {
      file_name = default_directory + "/" + file_name;
   }

   if (! coot::file_exists(file_name)) {
      std::cout << "ERROR:: missing file " << file_name << std::endl;
      std::cout << "ERROR:: not in " << default_directory << std::endl;
      return;
   }
#endif // THIS_IS_HMT

   int width, height, num_components;
   stbi_uc* image_data = stbi_load(file_name.c_str(), &width, &height, &num_components, 4);
   id = 0;

   if (!image_data) {
      std::string s = stbi_failure_reason();
      std::cout << "Error loading image data from " << file_name << " : " << s << std::endl;
      return;
   }

   glGenTextures(1, &m_texture_handle);
   glBindTexture(GL_TEXTURE_2D, m_texture_handle);

   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

   glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
   glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image_data);

   stbi_image_free(image_data);
   // std::cout << "debug::  done Texture::init() " << file_name << std::endl;
}

std::pair<int, int>
Texture::get_image_width_and_height() const {
   return std::make_pair(image_width, image_height);
}


void
Texture::init(const std::string &file_name, const std::string &directory) {

   std::string full_file_name = directory + "/" + file_name;
   init(full_file_name);

}


void
Texture::set_default_directory(const std::string &dir) {
   default_directory = dir;
}


void
Texture::close() {

   std::cout << "INFO:: deleting texture with id: " << id << " handle " << m_texture_handle << std::endl;
   glDeleteTextures(1, &m_texture_handle);

}

Texture::~Texture() {
   // don't delete the texture here - otherwise we can't copy textures.
}

void
Texture::Bind(unsigned int unit) {

   // unit should be less than 32
   glActiveTexture(GL_TEXTURE0 + unit);
   glBindTexture(GL_TEXTURE_2D, m_texture_handle);
   GLenum err = glGetError();
   if (err)
      std::cout << "GL Error:: in Texture::Bind() image from file " << file_name
                << " unit " << unit << " err " << err << std::endl;
}
