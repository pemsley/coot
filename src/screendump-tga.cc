/*
 * src/screendump-tga.cc
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
#include "utils/logging.hh"

extern logging logger;

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
                             int widget_width, int widget_height, int framebuffer_scale_factor) {

   // The caller is responsible for binding the correct framebuffer before calling this function.
   // This function reads pixels from whatever framebuffer is currently bound.

#ifdef __APPLE__
   // framebuffer_scale_factor *= 2;
#endif

   int w = widget_width;
   int h = widget_height;
   short int sf = static_cast<short int>(framebuffer_scale_factor);
   int pixel_width  = sf * w;
   int pixel_height = sf * h;

   FILE *output_file = fopen(tga_file.c_str(), "wb");
   if (!output_file) {
      std::cout << "ERROR:: screendump_tga_internal() could not open " << tga_file << " for writing" << std::endl;
      return;
   }

   unsigned char* pixel_data = new unsigned char[4 * pixel_width * pixel_height];
   short int sfw = static_cast<short int>(pixel_width);
   short int sfh = static_cast<short int>(pixel_height);
   short int TGAhead[] = {0, 2, 0, 0, 0, 0, sfw, sfh, 32};

   logger.log(log_t::INFO, logging::ltw("screendump_tga_internal() "), logging::ltw(pixel_width), logging::ltw("x"), logging::ltw(pixel_height), logging::ltw(" (scale " + std::to_string(sf) + ") to "), logging::ltw(tga_file));
   // std::cout << "INFO:: screendump_tga_internal() " << pixel_width << "x" << pixel_height
   //           << " (scale " << sf << ") to " << tga_file << std::endl;

   glReadBuffer(GL_COLOR_ATTACHMENT0);
   GLenum err = glGetError();
   if (err) std::cout << "ERROR:: screendump_tga_internal() post-glReadBuffer " << err << std::endl;

   glFinish();
   glReadPixels(0, 0, pixel_width, pixel_height, GL_BGRA, GL_UNSIGNED_BYTE, pixel_data);
   err = glGetError();
   if (err) std::cout << "ERROR:: screendump_tga_internal() post-glReadPixels " << err << std::endl;

   fwrite(&TGAhead, sizeof(TGAhead), 1, output_file);
   fwrite(pixel_data, 4 * pixel_width * pixel_height, 1, output_file);
   fclose(output_file);
   delete [] pixel_data;

   logger.log(log_t::INFO, logging::ltw("screendump_tga_internal() wrote "), logging::ltw(tga_file), logging::ltw(" "), logging::ltw(pixel_width), logging::ltw("x"), logging::ltw(pixel_height));
   // std::cout << "INFO:: screendump_tga_internal() wrote " << tga_file << " "
   //           << pixel_width << "x" << pixel_height << std::endl;

}

