
#include <iostream>
#include <fstream>
#include "utils/coot-utils.hh"
#include "pae.hh"
#include <cairo.h>
#include "json.hpp"
using json = nlohmann::json;

std::string
pae_t::file_to_string(const std::string &file_name) const {

   std::string s;
   std::string line;
   std::ifstream f(file_name.c_str());
   if (!f) {
      std::cout << "Failed to open " << file_name << std::endl;
   } else {
      while (std::getline(f, line)) {
         s += line;
         s += "\n";
      }
   }
   return s;
}


pae_t::pae_t(const std::string &file_name, int n_pixels_in) {

   auto file_to_string = [] (const std::string &file_name) {
      std::fstream f(file_name);
      std::string s;
      f.seekg(0, std::ios::end);
      s.reserve(f.tellg());
      f.seekg(0, std::ios::beg);
      s.assign((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
      return s;
   };

   n_pixels = n_pixels_in;
   if (coot::file_exists(file_name)) {
      std::string s = file_to_string(file_name);
      try {
         json j = json::parse(s);
         json item_1 = j[0];
         json ls = item_1["predicted_aligned_error"];
         std::vector<std::vector<int> > pae_vecs;
         for (json::iterator it_1=ls.begin(); it_1!=ls.end(); ++it_1) {
            json &item_2 = *it_1;
            // item_2 is a list with a number for each residue (self is 0)
            std::vector<int> arr = item_2;
            pae_vecs.push_back(arr);
         }
         bool is_square = true;
         unsigned int l = pae_vecs.size();
         for (auto &a : pae_vecs)
            if (a.size() != l) is_square = false;
         if (is_square) {
            image = make_image(pae_vecs);
         }

      }
      catch (const nlohmann::detail::type_error &e) {
         std::cout << "ERROR:: " << e.what() << std::endl;
      }
   }
}

std::string
pae_t::make_image(const std::vector<std::vector<int> > &pae_vecs) const {

   auto text_png_as_string_png_writer = [] (void *closure, const unsigned char *data, unsigned int length) {
      std::string *s_ptr = static_cast<std::string *>(closure);
      *s_ptr += std::string(reinterpret_cast<const char *>(data), length);
      return CAIRO_STATUS_SUCCESS;
   };

   std::string s;
   unsigned char *image_data = new unsigned char[n_pixels*n_pixels*4];
   // pre-colour with dark purple
   for (unsigned int i=0; i<n_pixels; i++) {
      for (unsigned int j=0; j<n_pixels; j++) {
         image_data[4*(i*n_pixels+j)]   = 70;
         image_data[4*(i*n_pixels+j)+1] = 0;
         image_data[4*(i*n_pixels+j)+2] = 70;
         image_data[4*(i*n_pixels+j)+3] = 255;
      }
   }

   unsigned int n_residues = pae_vecs.size();

   float max_value = 32.0;
   for (unsigned int i=0; i<n_pixels; i++) {
      for (unsigned int j=0; j<n_pixels; j++) {
         int idx = 4*(i*n_pixels+j);
         float f_x = static_cast<float>(i)/static_cast<float>(n_pixels);
         float f_y = static_cast<float>(j)/static_cast<float>(n_pixels);
         int i_res_index = static_cast<int>(f_x * static_cast<float>(n_residues));
         int j_res_index = static_cast<int>(f_y * static_cast<float>(n_residues));
         float f = static_cast<float>(pae_vecs[i_res_index][j_res_index])/max_value;
         unsigned char rb_col = static_cast<unsigned char>(255.0 * f);
         unsigned char  g_col = static_cast<unsigned char>(170.0 * f) + 85;
         if (rb_col < 0)   rb_col = 0;
         if (rb_col > 255) rb_col = 255;
         if (g_col < 0)   g_col = 0;
         if (g_col > 255) g_col = 255;
         image_data[idx]   = rb_col;
         image_data[idx+1] = g_col;
         image_data[idx+2] = rb_col;
      }
   }
   std::cout << "done setting image data" << std::endl;

   cairo_t *cr;
   cairo_surface_t *surface = cairo_image_surface_create_for_data(image_data, CAIRO_FORMAT_RGB24, n_pixels, n_pixels, n_pixels*4);
   cr = cairo_create(surface);
   cairo_set_source_rgb(cr, 0,0,0.5);

   if (cairo_surface_status(surface) == CAIRO_STATUS_SUCCESS) {
      std::cout << "########### cairo_surface_status() success " << std::endl;
      int w = cairo_image_surface_get_width(surface);
      int h = cairo_image_surface_get_height(surface);
      // std::cout << "on_draw_positron_plot() image_surface w h " << w << " " << h << std::endl;
      cairo_set_source_surface(cr, surface, 0,0);
      cairo_paint(cr);
      cairo_surface_write_to_png_stream(surface, text_png_as_string_png_writer,
                                        reinterpret_cast<void *>(&s));

   } else {
      std::cout << "########### cairo_surface_status() fail " << std::endl;
   }
   return s;

}