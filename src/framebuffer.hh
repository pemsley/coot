/*
 * src/framebuffer.hh
 *
 * Copyright 2019 by Medical Research Council
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
