
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

   GLenum err = glGetError();
   if (err)
      std::cout << "--------------------- start screen_framebuffer " << name
                << " init() err is " << err << std::endl;

   if (false)
      std::cout << "debug:: calling generate_framebuffer_object() " << name
                << " with " << width << " " << height << std::endl;
   generate_framebuffer_object(width, height,
                               attachment_index_color_texture); // try 2 * width here for
                                                                // supersampling at some stage

   err = glGetError();
   if (err) std::cout << "done framebuffer::init() with error " << err << std::endl;

}

framebuffer::~framebuffer() {

   // std::cout << "---------------- framebuffer destuctor - deleting framebuffer and textures " << std::endl;
   tear_down();
}

void
framebuffer::tear_down() {

   if (filled) {
      // std::cout << "framebuffer::tear_down()" << std::endl;
      glDeleteFramebuffers(1, &fbo);
      glDeleteTextures(1, &texture_colour);
      glDeleteTextures(1, &texture_depth);
      drawbuffer.clear();
   }
}

void
framebuffer::bind() {
   // GLint local_fbo;
   // glGetIntegerv(GL_FRAMEBUFFER_BINDING, &local_fbo);
   // std::cout << "Here in framebuffer::bind() pre  with local_fbo binding " << local_fbo << std::endl;
   glBindFramebuffer(GL_FRAMEBUFFER, fbo);
   // glGetIntegerv(GL_FRAMEBUFFER_BINDING, &local_fbo);
   // std::cout << "Here in framebuffer::bind() post with local_fbo binding " << local_fbo << std::endl;
   GLenum err = glGetError();
   if (err) std::cout << "framebuffer::bind() " << name << " fbo " << fbo << " err is " << err << std::endl;
}

void
framebuffer::generate_framebuffer_object(unsigned int width, unsigned int height, unsigned int attachment_index_color_texture) {

   // when we run init for the second time (after a resize) then we only need do some of these things?
   // or maybe clean up before we start assigning things again

   if (filled) {
      tear_down();
      filled = false;
   }

   glGenFramebuffers(1, &fbo);
   GLenum err = glGetError();
   if (err) std::cout << "--------------------- start generate_framebuffer_object() " << name
                      << " err is " << err << std::endl;
   glBindFramebuffer(GL_FRAMEBUFFER, fbo);
   err = glGetError(); if (err) std::cout << "--------------------- generate_framebuffer_object() A  " << name
                                          << " err is " << err << std::endl;

   generate_colourtexture(width, height);
   generate_depthtexture( width, height);
   err = glGetError(); if (err) std::cout << "--------------------- generate_framebuffer_object() B  "
                                          << name << " err is " << err << std::endl;

   // unsigned int attachment_index_color_texture = 0;

   // bind textures to pipeline. texture_depth is optional
   // 0 is the mipmap level. 0 is the heightest

   glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0 + attachment_index_color_texture, texture_colour, 0);
   err = glGetError();
   if (err) std::cout << "--------------------- generate_framebuffer_object() C  " << name << " err is " << err << std::endl;
   glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, texture_depth, 0);
   err = glGetError();
   if (err) std::cout << "--------------------- generate_framebuffer_object() D  " << name << " err is " << err << std::endl;

   // add attachments
   drawbuffer.push_back(GL_COLOR_ATTACHMENT0 + attachment_index_color_texture);
   if (false)
      std::cout << "--------------------- generate_framebuffer_object() " << name
                << " calling glDrawBuffers() with drawbuffer size "
                << drawbuffer.size() << " and contents of drawbuffer[0]: " << drawbuffer[0] << std::endl;
   glDrawBuffers(drawbuffer.size(), &drawbuffer[0]);
   err = glGetError(); if (err) std::cout << "--------------------- generate_framebuffer_object() E  "
                                          << name << " err is " << err << std::endl;

   // GL_INVALID_OPERATION (1282) is generated if any of the entries in bufs (other than GL_NONE) indicates a color buffer that does not
   // exist in the current GL context.
   
   if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
      std::cout << "xxxxxxxxxxxxxxxxxxxxx Error! FrameBuffer " << name << " is not complete" << std::endl;
   } else {
      // std::cout << "xxxxxxxxxxxxxxxxxxxxx FrameBuffer " << name << " is complete" << std::endl;
      filled = true;
   }

   err = glGetError(); if (err) std::cout << "--------------------- generate_framebuffer_object() " << name
                                          << " end err is " << err << std::endl;

}

void
framebuffer::generate_colourtexture(unsigned int width, unsigned int height) {

   glGenTextures(1, &texture_colour);
   glBindTexture(GL_TEXTURE_2D, texture_colour);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
   // std::cout << "in framebuffer::generate_colourtexture() with texture_colour " << texture_colour << std::endl;

}

void
framebuffer::generate_depthtexture(unsigned int width, unsigned int height) {

   glGenTextures(1, &texture_depth);
   glBindTexture(GL_TEXTURE_2D, texture_depth);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
   glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, width, height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);
   // std::cout << "in framebuffer::generate_depthtexture() with texture_depth " << texture_depth << std::endl;
   
}


