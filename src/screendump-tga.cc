
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

void screendump_tga_internal(std::string tga_file,
                             int widget_height, int widget_width, int image_scale_factor,
                             framebuffer &framebuffer_for_screen) {

   int w = widget_width;
   int h = widget_height;
   FILE *output_file = fopen(tga_file.c_str(), "w");
   short int sf = static_cast<short int>(image_scale_factor);
   unsigned char* pixel_data = new unsigned char[4 * sf * sf * w * h]; // 4 components?
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

   framebuffer_for_screen.bind();
   glNamedFramebufferReadBuffer(framebuffer_for_screen.get_fbo(), GL_FRONT);

   glReadPixels(0, 0, sf * w, sf * h, GL_BGR, GL_UNSIGNED_BYTE, pixel_data);
   fwrite(&TGAhead, sizeof(TGAhead), 1, output_file);
   fwrite(pixel_data, 3 * sf * sf * w * h, 1, output_file); // 4 components
   fclose(output_file);
   delete [] pixel_data;

   std::cout << "screendump_tga done, wrote " << 3 * sf * sf * w * h << " bytes" << std::endl;

}

