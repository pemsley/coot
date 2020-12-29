
#include <iostream>
#include <epoxy/gl.h>

#if THIS_IS_HMT
#include "display-info.hh"
#else
#ifdef USE_PYTHON // this test block should be in graphics-info.h - for another time.
#include "Python.h"
#endif
#include "graphics-info.h"
#endif

#include "screendump-tga.hh"

void set_framebuffer_scale_factor(unsigned int sf) {

   graphics_info_t g;
   gtk_gl_area_attach_buffers(GTK_GL_AREA(g.glareas[0]));
   gtk_gl_area_make_current(GTK_GL_AREA(g.glareas[0]));
   g.framebuffer_scale = sf;
   GtkAllocation allocation = g.get_glarea_allocation();
   unsigned int w = allocation.width;
   unsigned int h = allocation.height;
   g.reset_frame_buffers(w,h);
   g.graphics_draw();
}



void screendump_tga_internal(std::string tga_file,
                             int widget_width, int widget_height, int image_scale_factor,
                             unsigned int framebuffer_obj) { // care: pass by reference

   GLenum err = glGetError();
   if (err) std::cout << "error:: screendump_tga_internal() start " << err << std::endl;
   int w = widget_width;
   int h = widget_height;
   FILE *output_file = fopen(tga_file.c_str(), "w");
   short int sf = static_cast<short int>(image_scale_factor);
   unsigned char* pixel_data = new unsigned char[4 * sf * sf * w * h]; // 4 components? Probably not right.
   short int sfw = static_cast<short int>(sf * widget_width);
   short int sfh = static_cast<short int>(sf * widget_height);
   short int TGAhead[] = {0, 2, 0, 0, 0, 0, sfw, sfh, 24};

   std::cout << "screendump_tga application image: scaling " << sf << " " << w << " x " << h
             << " to " << tga_file << std::endl;

#if THIS_IS_HMT
   display_info_t di;
   di.attach_buffers();
#else
   gtk_gl_area_attach_buffers(GTK_GL_AREA(graphics_info_t::glareas[0]));
#endif
   err = glGetError();
   if (err) std::cout << "error:: screendump_tga_internal() post-attach " << err << std::endl;

   GLint local_fbo;
   glGetIntegerv(GL_FRAMEBUFFER_BINDING, &local_fbo);
   std::cout << "debug::  pre-bind with local_fbo binding " << local_fbo << std::endl;

   // framebuffer.bind();
   glBindFramebuffer(GL_FRAMEBUFFER, framebuffer_obj);

   err = glGetError(); if (err) std::cout << "error:: screendump_tga_internal() post-bind " << err << std::endl;

   glGetIntegerv(GL_FRAMEBUFFER_BINDING, &local_fbo);
   std::cout << "debug:: post-bind with local_fbo binding " << local_fbo << std::endl;

   std::cout << "debug:: Using framebuffer fbo " << framebuffer_obj << std::endl;

   glNamedFramebufferReadBuffer(framebuffer_obj, GL_BACK); // this often errors
   err = glGetError();
   if (err) std::cout << "error:: screendump_tga_internal() post-set glnamedreadbuffer "
                      << err << std::endl;

   glReadPixels(0, 0, sf * w, sf * h, GL_BGR, GL_UNSIGNED_BYTE, pixel_data);
   err = glGetError(); if (err) std::cout << "error:: screendump_tga_internal() post-glReadpixels "
                                          << err << std::endl;
   fwrite(&TGAhead, sizeof(TGAhead), 1, output_file);
   fwrite(pixel_data, 3 * sf * sf * w * h, 1, output_file); // 3 components
   fclose(output_file);
   delete [] pixel_data;

   std::cout << "INFO:: screendump_tga sf " << sf << " " << w << "x" << h
             << " wrote " << 3 * sf * sf * w * h << " bytes" << std::endl;

}

