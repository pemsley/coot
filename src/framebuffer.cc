
#include <iostream>
#include <gtk/gtk.h>
#include <epoxy/gl.h>



void setup_frame_buffer() {
   unsigned int fbo;
   glGenFramebuffers(1, &fbo);
   GLenum err = glGetError(); std::cout << " start setup_frame_buffer() err is " << err << std::endl;
   glBindFramebuffer(GL_FRAMEBUFFER, fbo);
   err = glGetError(); if (err) std::cout << " setup_frame_buffer() A err is " << err << std::endl;

   unsigned int texture;
   glGenTextures(1, &texture);
   err = glGetError(); if (err) std::cout << " setup_frame_buffer() B err is " << err << std::endl;
   glBindTexture(GL_TEXTURE_2D, texture);
   err = glGetError(); if (err) std::cout << " setup_frame_buffer() C err is " << err << std::endl;

   glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 800, 600, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
   err = glGetError(); if (err) std::cout << " setup_frame_buffer() D err is " << err << std::endl;

   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
   err = glGetError(); if (err) std::cout << " setup_frame_buffer() E err is " << err << std::endl;
   glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
   err = glGetError(); if (err) std::cout << " setup_frame_buffer() F err is " << err << std::endl;
   glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texture, 0);
   err = glGetError(); if (err) std::cout << " setup_frame_buffer() G err is " << err << std::endl;
   glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT,  GL_TEXTURE_2D, texture, 0);
   err = glGetError(); if (err) std::cout << " setup_frame_buffer() H err is " << err << std::endl;

   if (glCheckFramebufferStatus(GL_FRAMEBUFFER) == GL_FRAMEBUFFER_COMPLETE) {
      std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Framebuffer Setup OK\n";
      glBindFramebuffer(GL_FRAMEBUFFER, 0);
   } else {
      std::cout << "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Framebuffer Setup incomplete...\n";
   }

}
