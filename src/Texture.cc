/*
 * src/Texture.cc
 *
 * Copyright 2020 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#include <iostream>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>

#include "../coot-utils/stb_image.h"
#include "Texture.hh"

#ifdef THIS_IS_HMT
#else
#include "utils/coot-utils.hh"
#endif

Texture::Texture(const std::string &file_name, type_t t, bool reversed_normals_in) {
   type = t;
   init(file_name);
   reversed_normals = reversed_normals_in;
}

// create a solid colour texture
Texture::Texture(int image_width_in, int image_height_in, glm::vec4 colour) {

   image_width = image_width_in;
   image_height = image_height_in;

   // unsigned char image_data[image_height * image_width * 4];
   unsigned char *image_data = new unsigned char[image_height * image_width * 4];

   unsigned char r = colour[0] * 255;
   unsigned char g = colour[1] * 255;
   unsigned char b = colour[2] * 255;

   for (int i=0; i<image_width; i++) {
      for (int j=0; j<image_height; j++) {
         unsigned int idx = i * image_height + j;
         image_data[4 * idx + 0] = r;
         image_data[4 * idx + 1] = g;
         image_data[4 * idx + 2] = b;
         image_data[4 * idx + 3] = 255;
      }
   }

   glGenTextures(1, &m_texture_handle);
   glBindTexture(GL_TEXTURE_2D, m_texture_handle);

   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

   glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
   glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, image_width, image_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image_data);

   delete [] image_data;

}

// create a colour bar
Texture::Texture(int image_width_in, int image_height_in, const std::vector<glm::vec4> &colours) {

   colour_bar(image_width_in, image_height_in, colours, 0);
}

Texture::Texture(int image_width_in, int image_height_in, const std::vector<glm::vec4> &colours, unsigned int n_ticks) {

   colour_bar(image_width_in, image_height_in, colours, n_ticks);
}

Texture::Texture(const mini_texture_t &mt, const std::string &name_in) {

   auto _ = [] (int err) {
      std::string s = std::to_string(err);
      if (err == GL_INVALID_ENUM)      s = "GL_INVALID_ENUM";
      if (err == GL_INVALID_OPERATION) s = "GL_INVALID_OPERATION";
      if (err == GL_INVALID_VALUE)     s = "GL_INVALID_VALUE";
      return s;
   };

   file_name = name_in; // it's not a file

   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: Texture::Texture(mt) A " << _(err) << "\n";

   glGenTextures(1, &m_texture_handle);
   glBindTexture(GL_TEXTURE_2D, m_texture_handle);

   err = glGetError();
   if (err) std::cout << "GL ERROR:: Texture::Texture(mt) C " << _(err) << "\n";

   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

   err = glGetError();
   if (err) std::cout << "GL ERROR:: Texture::Texture(mt) E " << _(err) << "\n";

   if (true) {
      glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
   } else {
      glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
   }

   err = glGetError();
   if (err) std::cout << "GL ERROR:: Texture::Texture(mt) G " << _(err) << "\n";

   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, mt.width, mt.height, 0, GL_RGBA, GL_UNSIGNED_BYTE, mt.image_data);

   err = glGetError();
   if (err) std::cout << "GL ERROR:: Texture::Texture(mt) H " << _(err)
                      << " w and h " << mt.width << " " << mt.height << std::endl;

}

// Don't draw ticks if n_ticks < 2.
void
Texture::colour_bar(int image_width_in, int image_height_in, const std::vector<glm::vec4> &colours, unsigned int n_ticks) {

  if (colours.empty()) {
      std::cout << "ERROR:: failure to create Texture because colours was empty." << std::endl;
      return;
   }

   image_width = image_width_in;
   image_height = image_height_in;

   // unsigned char image_data[image_height * image_width * 4];
   unsigned char *image_data = new unsigned char[image_height * image_width * 4];
   float s = colours.size();

   // image stored in rows
   //
   for (int j=0; j<image_height; j++) {

      for (int i=0; i<image_width; i++) {

         float f = static_cast<float>(i)/static_cast<float>(image_width);
         unsigned int colours_idx = f * s;

         if (colours_idx > colours.size()) colours_idx = colours.size() -1;

         auto colour = colours[colours_idx];
         if (colour[0] > 1.0) colour[0] = 1.0;
         if (colour[1] > 1.0) colour[1] = 1.0;
         if (colour[2] > 1.0) colour[2] = 1.0;
         unsigned char r = colour[0] * 255;
         unsigned char g = colour[1] * 255;
         unsigned char b = colour[2] * 255;

         // std::cout << "colour_idx " << colours_idx << " colours size " << colours.size() << " " << glm::to_string(colour)
         // << " " << static_cast<int>(r) << " " << static_cast<int>(g) << " " << static_cast<int>(b) << std::endl;
         // if (i > 390) { r = 122; g = 0; b = 122; } testing

         unsigned int idx = j * image_width + i;
         image_data[4 * idx + 0] = r;
         image_data[4 * idx + 1] = g;
         image_data[4 * idx + 2] = b;
         image_data[4 * idx + 3] = 255;
      }
   }

   if (n_ticks > 1) {
      glm::vec4 tick_colour(0.155, 0.155, 0.155, 1.0);
      add_tick_marks(n_ticks, tick_colour, image_data);
   }

   glGenTextures(1, &m_texture_handle);
   glBindTexture(GL_TEXTURE_2D, m_texture_handle);

   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

   glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
   glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, image_width, image_height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image_data);

   delete [] image_data;
}


void
Texture::init(const std::string &file_name_in) {

   file_name = file_name_in;

#ifdef THIS_IS_HMT
#else

   if (false)
      std::cout << "Texture::init() was passed file_name_in " << file_name_in << std::endl;

   if (default_directory.empty()) {
      std::string pkg_data_dir = coot::package_data_dir();
      default_directory = coot::util::append_dir_dir(pkg_data_dir, "textures");
   }

   if (! coot::file_exists(file_name)) {
      file_name = coot::util::append_dir_file(default_directory, file_name);
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

void
Texture::handle_raw_image_data(const std::string &image_name,
                               const std::vector<unsigned char> &image_data,
                               int width, int height) {

   file_name = image_name;
   glGenTextures(1, &m_texture_handle);
   glBindTexture(GL_TEXTURE_2D, m_texture_handle);

   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

   glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
   glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image_data.data());

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
   if (false) // debug
      std::cout << "   Texture::Bind() " << file_name << " " << type << " unit " << unit
                << " m_texture_handle " << m_texture_handle << std::endl;

   glActiveTexture(GL_TEXTURE0 + unit);
   glBindTexture(GL_TEXTURE_2D, m_texture_handle);
   GLenum err = glGetError();
   if (err)
      std::cout << "GL Error:: in Texture::Bind() image from file \"" << file_name << "\""
                << " unit " << unit << " err " << err << std::endl;
}

void
Texture::add_tick_marks(unsigned int n_ticks, const glm::vec4 &tick_colour, unsigned char *image_data) {

   // strip off the top 2 rows of pixels and make them black and transparent.

   // blank the rows
   for (int irow=0; irow<40; irow++) {
      int j = irow;
      for (int i=0; i<image_width; i++) {
         unsigned int idx = j * image_width + i;
         image_data[4 * idx + 0] = 0;
         image_data[4 * idx + 1] = 0;
         image_data[4 * idx + 2] = 0;
         image_data[4 * idx + 3] = 0;
      }
   }

   // add ticks

   for (unsigned int i_tick=0; i_tick<n_ticks; i_tick++) {

      float f = static_cast<float>(i_tick)/static_cast<float>(n_ticks-1);
      int i = f * image_width;
      if (i >= image_width) i = image_width - 1;
      int tick_height = 100;
      if (tick_height > image_height) tick_height = image_height;
      for (int j=0; j<tick_height; j++) {

         int idx = j * image_width + i;

         if (idx >= image_width * image_height) {
            std::cout << "ERROR " << idx << std::endl;
         } else {
            unsigned char r = tick_colour[0] * 255;
            unsigned char g = tick_colour[1] * 255;
            unsigned char b = tick_colour[2] * 255;
            unsigned char a = tick_colour[3] * 255;

            if (j < 40) { r = 255; g = 255; b = 255; }

            image_data[4 * idx + 0] = r;
            image_data[4 * idx + 1] = g;
            image_data[4 * idx + 2] = b;
            image_data[4 * idx + 3] = a;
         }
      }
   }

}
