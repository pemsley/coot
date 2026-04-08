/*
 * src/framebuffer.cc
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

#include <iostream>
#include <vector>

#include <gtk/gtk.h>
#include <epoxy/gl.h>

#include "framebuffer.hh"

framebuffer::framebuffer() {
   filled = false;
}

// caller does the scaling for width and heigh
void
framebuffer::init(int width, int height, unsigned int attachment_index_color_texture, const std::string &name_in) {

   name = name_in;

   if (width  == 0) std::cout << "ERROR:: in framebuffer::init() " << name << " width is 0" << std::endl;
   if (height == 0) std::cout << "ERROR:: in framebuffer::init() " << name << " height is 0" << std::endl;

   GLenum err = glGetError();
   if (err)
      std::cout << "--- start framebuffer " << name << " init() err is " << err << std::endl;

   if (false)
      std::cout << "debug:: framebuffer::init() calling generate_framebuffer_object() " << name
                << " with " << width << " " << height << std::endl;

   // std::cout << "in framebuffer::init() Here 1 " << name_in << std::endl;
   generate_framebuffer_object(width, height, attachment_index_color_texture);
   // std::cout << "in framebuffer::init() Here 2 " << name_in << std::endl;
   err = glGetError();
   if (err) std::cout << "GL ERROR:: finish framebuffer::init() with error " << err << std::endl;

}

// on window resize
void
framebuffer::reset(int width, int height) {

   GLenum err = glGetError();
   if (err)
      std::cout << "--- start framebuffer " << name << " init() err is " << err << std::endl;

   if (false)
      std::cout << "debug:: framebuffer::reset() calling generate_framebuffer_object() " << name
                << " with " << width << " " << height << std::endl;

   // 20220108-PE all my framebuffers only have 1 color attachment at the moment.

   unsigned int attachment_index_for_color_texture = 0;
   generate_framebuffer_object(width, height, attachment_index_for_color_texture);

   // // 20220108-PE alternatively I could run through all of the elements of drawbuffer
   // for (const auto &d : drawbuffer) {
   //    // GL_COLOR_ATTACHMENT0 is added to the attachment_index_for_color_texture when it gets added
   //    // to the drawbuffer vector
   //    generate_framebuffer_object(width, height, d - GL_COLOR_ATTACHMENT0);
   // }

   err = glGetError();
   if (err) std::cout << "done framebuffer::init() with error " << err << std::endl;

}

// on window resize
void
framebuffer::reset_test(int width, int height) {

   // this function is called in graphics-info-draw.cc and graphics-info-opengl.cc
   // it is not a test function and should be renamed

   GLenum err = glGetError();
   if (err)
      std::cout << "--- start framebuffer " << name << " init() err is " << err << std::endl;

   if (false)
      std::cout << "debug:: framebuffer::reset_test() calling generate_framebuffer_object() " << name
                << " with " << width << " " << height << std::endl;

   // 20220108-PE all my framebuffers only have 1 color attachment at the moment.

   // unsigned int attachment_index_for_color_texture = 0;
   // generate_framebuffer_object_test(width, height, attachment_index_for_color_texture);
   // attachment_index_for_color_texture = 1;
   //    generate_framebuffer_object_test(width, height, attachment_index_for_color_texture);

   // // 20220108-PE alternatively I could run through all of the elements of drawbuffer
   // for (const auto &d : drawbuffer) {
   //    // GL_COLOR_ATTACHMENT0 is added to the attachment_index_for_color_texture when it gets added
   //    // to the drawbuffer vector
   //    generate_framebuffer_object(width, height, d - GL_COLOR_ATTACHMENT0);
   // }

   err = glGetError();
   if (err) std::cout << "done framebuffer::init() with error " << err << std::endl;

   do_gbuffer_stuff(width, height);

}



framebuffer::~framebuffer() {

   // std::cout << "--- framebuffer destructor - deleting framebuffer and textures \"" << name << "\"" << std::endl;
   tear_down();
}

void
framebuffer::tear_down() {

   if (filled) {
      // std::cout << "framebuffer::tear_down()" << std::endl;
      glDeleteFramebuffers(1, &fbo); // 20250615-PE crash. // 20251123-PE crash. Hmm.
      glDeleteTextures(1, &colour_texture);
      glDeleteTextures(1, &depth_texture);
      drawbuffer.clear();
   }
}

void
framebuffer::bind() {

   GLint local_fbo;
   glGetIntegerv(GL_FRAMEBUFFER_BINDING, &local_fbo);
   // std::cout << "Here in framebuffer::bind() " << name << " with local_fbo pre-binding  " << local_fbo << std::endl;
   glBindFramebuffer(GL_FRAMEBUFFER, fbo);
   glGetIntegerv(GL_FRAMEBUFFER_BINDING, &local_fbo);
   // std::cout << "Here in framebuffer::bind() " << name << " with local_fbo post-binding " << local_fbo << std::endl;

   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: framebuffer::bind() " << name << " fbo " << fbo << " err is " << err << std::endl;
}

void
framebuffer::generate_framebuffer_object(unsigned int width, unsigned int height, unsigned int attachment_index_color_texture) {

   // when we run init for the second time (after a resize) then we only need do some of these things?
   // or maybe clean up before we start assigning things again

   if (filled) {
      // this happens when this gets called from a window resize callback
      tear_down();
      filled = false;
   }

   glGenFramebuffers(1, &fbo);
   GLenum err = glGetError();
   if (err) std::cout << "--- start generate_framebuffer_object() " << name
                      << " err is " << err << std::endl;

   glBindFramebuffer(GL_FRAMEBUFFER, fbo);
   err = glGetError();
   if (err) std::cout << "--- generate_framebuffer_object() A post glBindFramebuffer() "
		      << name << " err is " << err << std::endl;
   generate_colourtexture(width, height);
   err = glGetError();
   if (err) std::cout << "---- generate_framebuffer_object() post generate_colourtexture() "
		      << name << " err is " << err << std::endl;
   generate_depthtexture(width, height);
   err = glGetError();
   if (err) std::cout << "---- generate_framebuffer_object() post generate_depthtexture() "
		      << name << " err is " << err << std::endl;

   // unsigned int attachment_index_color_texture = 0;

   // bind textures to pipeline. texture_depth is optional
   // 0 is the mipmap level. 0 is the heightest

   glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0 + attachment_index_color_texture, colour_texture, 0);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: generate_framebuffer_object() C \"" << name << "\" err is " << err << std::endl;
   glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, depth_texture, 0);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: generate_framebuffer_object() D \"" << name << "\" err is " << err << std::endl;

   if (false)
      std::cout << "debug:: framebuffer::generate_framebuffer_object() currently drawbuffer.size() is " << drawbuffer.size()
                << " for " << name << std::endl;

   // add attachments
   drawbuffer.push_back(GL_COLOR_ATTACHMENT0 + attachment_index_color_texture);
   if (false)
      std::cout << "GL ERROR:: generate_framebuffer_object() \"" << name
                << "\" calling glDrawBuffers() with drawbuffer size "
                << drawbuffer.size() << " and contents of drawbuffer[0]: " << drawbuffer[0] << std::endl;
   glDrawBuffers(drawbuffer.size(), &drawbuffer[0]);
   err = glGetError(); if (err) std::cout << "GL ERROR:: generate_framebuffer_object() E \""
                                          << name << "\" err is " << err << std::endl;

   // GL_INVALID_OPERATION (1282) is generated if any of the entries in bufs (other than GL_NONE) indicates a color buffer that does not
   // exist in the current GL context.

   if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
      std::cout << "xxxxxxx GL ERROR:: framebuffer::generate_framebuffer_object() FrameBuffer \""
                << name << "\" width " << width << " height " << height << " is not complete" << std::endl;
   } else {
      // std::cout << "xxxxxxxxxxxxxxxxxxxxx FrameBuffer " << name << " is complete" << std::endl;
      filled = true;
   }

   err = glGetError(); if (err) std::cout << "--------------------- generate_framebuffer_object() " << name
                                          << " end err is " << err << std::endl;

   //std::cout << "in framebuffer::generate_framebuffer_object() done " << name << std::endl;

}

void
framebuffer::generate_framebuffer_object_test(unsigned int width, unsigned int height, unsigned int attachment_index_color_texture) {

   // when we run init for the second time (after a resize) then we only need do some of these things?
   // or maybe clean up before we start assigning things again

   std::cout << "debug:: in framebuffer::generate_framebuffer_object_test() Here 1 " << name << " with drawbuffer size "
             << drawbuffer.size() << " drawbuffer[0] " << drawbuffer[0] << std::endl;

   if (filled) {
      // this happens when this gets called from a window resize callback
      tear_down();
      filled = false;
   }


   glGenFramebuffers(1, &fbo);
   GLenum err = glGetError();
   if (err) std::cout << "--- start generate_framebuffer_object() " << name
                      << " err is " << err << std::endl;

   glBindFramebuffer(GL_FRAMEBUFFER, fbo);
   err = glGetError();
   if (err) std::cout << "--- generate_framebuffer_object() A post glBindFramebuffer() "
		      << name << " err is " << err << std::endl;
   generate_colourtexture(width, height);
   err = glGetError();
   if (err) std::cout << "---- generate_framebuffer_object() post generate_colourtexture() "
		      << name << " err is " << err << std::endl;
   generate_depthtexture(width, height);
   err = glGetError();
   if (err) std::cout << "---- generate_framebuffer_object() post generate_depthtexture() "
		      << name << " err is " << err << std::endl;

   // unsigned int attachment_index_color_texture = 0;

   // bind textures to pipeline. texture_depth is optional
   // 0 is the mipmap level. 0 is the heightest

   glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0 + attachment_index_color_texture, colour_texture, 0);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: generate_framebuffer_object() C \"" << name << "\" err is " << err << std::endl;
   glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, depth_texture, 0);
   err = glGetError();
   if (err) std::cout << "GL ERROR:: generate_framebuffer_object() D \"" << name << "\" err is " << err << std::endl;

   std::cout << "debug:: in framebuffer::generate_framebuffer_object() A currently drawbuffer.size() is "
             << drawbuffer.size() << " for " << name << std::endl;

   // add attachments
   drawbuffer.push_back(GL_COLOR_ATTACHMENT0 + attachment_index_color_texture);
   if (false)
      std::cout << "GL ERROR:: generate_framebuffer_object() \"" << name
                << "\" calling glDrawBuffers() with drawbuffer size "
                << drawbuffer.size() << " and contents of drawbuffer[0]: " << drawbuffer[0] << std::endl;
   glDrawBuffers(drawbuffer.size(), &drawbuffer[0]);
   err = glGetError(); if (err) std::cout << "GL ERROR:: generate_framebuffer_object() E \""
                                          << name << "\" err is " << err << std::endl;

   std::cout << "debug:: in framebuffer::generate_framebuffer_object() B currently drawbuffer.size() is "
             << drawbuffer.size() << " for " << name << std::endl;

   std::cout << "debug:: in framebuffer::generate_framebuffer_object_test() Here 2 " << name << " with drawbuffer size "
             << drawbuffer.size() << " drawbuffer[0] " << drawbuffer[0] << std::endl;

   // GL_INVALID_OPERATION (1282) is generated if any of the entries in bufs (other than GL_NONE) indicates a color buffer that does not
   // exist in the current GL context.

   if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
      std::cout << "xxxxxxxxxxxxxxxxxxxxx Error! FrameBuffer \"" << name << "\" is not complete" << std::endl;
   } else {
      // std::cout << "xxxxxxxxxxxxxxxxxxxxx FrameBuffer " << name << " is complete" << std::endl;
      filled = true;
   }

   err = glGetError(); if (err) std::cout << "--------------------- generate_framebuffer_object() " << name
                                          << " end err is " << err << std::endl;

   //std::cout << "in framebuffer::generate_framebuffer_object_test() done " << name << " filled " << filled << std::endl;
}

void
framebuffer::generate_colourtexture(unsigned int width, unsigned int height) {

   GLenum err = glGetError();
   if (err) std::cout << "ERROR generate_colourtexture() --start--  "
		      << name << " err is " << err << std::endl;
   glGenTextures(1, &colour_texture);
   err = glGetError();
   if (err) std::cout << "ERROR generate_colourtexture() A "
		      << name << " err is " << err << std::endl;
   glBindTexture(GL_TEXTURE_2D, colour_texture);
   err = glGetError();
   if (err) std::cout << "ERROR generate_colourtexture() B "
		      << name << " err is " << err << std::endl;
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
   err = glGetError();
   if (err) std::cout << "ERROR generate_colourtexture() C "
		      << name << " err is " << err << std::endl;
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
   err = glGetError();
   if (err) std::cout << "ERROR generate_colourtexture() D "
		      << name << " err is " << err << std::endl;
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
   err = glGetError();
   if (err) std::cout << "ERROR generate_colourtexture() E "
		      << name << " err is " << err << std::endl;
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
   err = glGetError();
   if (err) std::cout << "ERROR generate_colourtexture() F "
		      << name << " err is " << err << std::endl;
   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
   err = glGetError();
   if (err) std::cout << "ERROR generate_colourtexture() G "
		      << name << " err is " << err << std::endl;
   // std::cout << "in framebuffer::generate_colourtexture() with texture_colour "
   // << texture_colour << std::endl;

}

void
framebuffer::generate_depthtexture(unsigned int width, unsigned int height) {

   glGenTextures(1, &depth_texture);
   glBindTexture(GL_TEXTURE_2D, depth_texture);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
   glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, width, height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);
   // std::cout << "in framebuffer::generate_depthtexture() with texture_depth " << texture_depth << std::endl;
   
}


void
framebuffer::do_gbuffer_stuff(int w, int h) {

   glBindFramebuffer(GL_FRAMEBUFFER, fbo);

   glGenTextures(1, &gPosition);
   glBindTexture(GL_TEXTURE_2D, gPosition);
   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, w, h, 0, GL_RGBA, GL_FLOAT, NULL);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
   glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, gPosition, 0);
   // normal color buffer
   glGenTextures(1, &gNormal);
   glBindTexture(GL_TEXTURE_2D, gNormal);
   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, w, h, 0, GL_RGBA, GL_FLOAT, NULL);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
   glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, gNormal, 0);



}
