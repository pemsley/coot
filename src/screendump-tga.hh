
#ifndef SCREENDUMP_TGA_HH
#define SCREENDUMP_TGA_HH

#include <string>
#include <gtk/gtk.h>
#include "framebuffer.hh"

// There is an API function called screendump_tga now

void screendump_tga_internal(std::string tga_file,
                             int widget_height, int widget_width, int image_scale_factor,
                             unsigned int framebuffer_obj);

#endif // SCREENDUMP_TGA_HH
