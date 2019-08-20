
#include <iostream>
#include <vector>

#include <gtk/gtk.h>
#include <epoxy/gl.h>

#include "framebuffer.hh"

framebuffer::framebuffer() {
}

void
framebuffer::init(int width, int height) {

   GLenum err = glGetError(); std::cout << "--------------------- start screen_framebuffer init() err is " << err << std::endl;
   generate_framebuffer_object(width, height);

}

framebuffer::~framebuffer() {
   glDeleteFramebuffers(1, &fbo);
   glDeleteTextures(1, &texture_colour);
   glDeleteTextures(1, &texture_depth);

}

void
framebuffer::bind() {
   glBindFramebuffer(GL_FRAMEBUFFER, fbo);
}

void
framebuffer::generate_framebuffer_object(unsigned int width, unsigned int height) {

   glGenFramebuffers(1, &fbo);
   GLenum err = glGetError(); std::cout << "--------------------- start generate_framebuffer_object() err is " << err << std::endl;
   glBindFramebuffer(GL_FRAMEBUFFER, fbo);
   err = glGetError(); if (true) std::cout << "--------------------- generate_framebuffer_object() A err is " << err << std::endl;

   generate_colourtexture(width, height);
   generate_depthtexture( width, height);
   err = glGetError(); if (true) std::cout << "--------------------- generate_framebuffer_object() B err is " << err << std::endl;

   unsigned int attachment_index_color_texture = 0;

   // bind textures to pipeline. texture_depth is optional
   // 0 is the mipmap level. 0 is the heightest

   glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0 + attachment_index_color_texture, texture_colour, 0);
   glFramebufferTexture(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, texture_depth, 0);

   // add attachments
   drawbuffer.push_back(GL_COLOR_ATTACHMENT0 + attachment_index_color_texture);
   glDrawBuffers(drawbuffer.size(), &drawbuffer[0]);

   if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
      std::cout << "xxxxxxxxxxxxxxxxxxxxx Error! FrameBuffer is not complete" << std::endl;
   else
      std::cout << "xxxxxxxxxxxxxxxxxxxxx FrameBuffer is complete" << std::endl;

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

}




std::tuple<unsigned int, unsigned int, unsigned int>
setup_frame_buffer_old(int width, int height) {

   unsigned int fbo; // returned first
   unsigned int textureColorbuffer; // returned second
   unsigned int textureDepthbuffer; // returned last

   glGenFramebuffers(1, &fbo);
   GLenum err = glGetError(); std::cout << " start setup_frame_buffer() err is " << err << std::endl;
   glBindFramebuffer(GL_FRAMEBUFFER, fbo);
   err = glGetError(); if (err) std::cout << " setup_frame_buffer() A err is " << err << std::endl;

   glGenTextures(1, &textureColorbuffer);
   err = glGetError(); if (err) std::cout << " setup_frame_buffer() B err is " << err << std::endl;
   glBindTexture(GL_TEXTURE_2D, textureColorbuffer);
   err = glGetError(); if (err) std::cout << " setup_frame_buffer() C err is " << err << std::endl;

   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
   err = glGetError(); if (err) std::cout << " setup_frame_buffer() D err is " << err << std::endl;

   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
   err = glGetError(); if (err) std::cout << " setup_frame_buffer() E err is " << err << std::endl;
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
   err = glGetError(); if (err) std::cout << " setup_frame_buffer() F err is " << err << std::endl;
   glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, textureColorbuffer, 0);
   err = glGetError(); if (err) std::cout << " setup_frame_buffer() G err is " << err << std::endl;
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,  GL_TEXTURE_2D, textureColorbuffer, 0);
   err = glGetError(); if (err) std::cout << " setup_frame_buffer() H err is " << err << std::endl;

   unsigned int rbo;
   glGenRenderbuffers(1, &rbo);
   glBindRenderbuffer(GL_RENDERBUFFER, rbo);
   // glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, width, height);
   glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT, width, height);
   // glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, rbo);
   glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rbo);

   if (glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE) {
      std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Framebuffer Setup OK\n";
   } else {
      std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Framebuffer Setup incomplete...\n";
   }

   glBindFramebuffer(GL_FRAMEBUFFER, 0);

   return std::tuple<unsigned int, unsigned int, unsigned int> (fbo, textureColorbuffer, rbo);
}
