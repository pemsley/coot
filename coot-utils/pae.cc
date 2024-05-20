
#include <iostream>
#include <fstream>
#include "utils/coot-utils.hh"
#include "pae.hh"
#if RDKIT_HAS_CAIRO_SUPPORT // acutally just a test for CAIRO availability
#include <cairo.h>
#endif
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

float
pae_t::get_max_value(const std::vector<std::vector<int> > &pae_vecs) const {

   float max_value = 0.0;
   for (const auto &pae : pae_vecs) {
      for (const int &epe : pae) {
         if (epe > max_value) {
            max_value = epe;
         }
      }
   }
   return max_value;
}

std::string
pae_t::make_image(const std::vector<std::vector<int> > &pae_vecs) const {

#ifdef RDKIT_HAS_CAIRO_SUPPORT

   auto text_png_as_string_png_writer = [] (void *closure, const unsigned char *data, unsigned int length) {
      std::string *s_ptr = static_cast<std::string *>(closure);
      *s_ptr += std::string(reinterpret_cast<const char *>(data), length);
      return CAIRO_STATUS_SUCCESS;
   };

   auto value_to_colour = [] (float value, float max_value) {
      std::vector<unsigned char> col(3, 0);
      float f = value/max_value;
      unsigned char rb_col = static_cast<unsigned char>(255.0 * f);
      unsigned char  g_col = static_cast<unsigned char>(190.0 * f) + 65;
      if (rb_col > 255) rb_col = 255;
      if (g_col  > 255)  g_col = 255;
      col[0] = rb_col;
      col[1] =  g_col;
      col[2] = rb_col;
      return col;
   };

   int n_pixels_for_pae_image = n_pixels - 100;

   std::string s;
   unsigned char *image_data = new unsigned char[n_pixels_for_pae_image*n_pixels_for_pae_image*4];
   // pre-colour with dark purple
   for (int i=0; i<n_pixels_for_pae_image; i++) {
      for (int j=0; j<n_pixels_for_pae_image; j++) {
         image_data[4*(i*n_pixels_for_pae_image+j)]   = 70;
         image_data[4*(i*n_pixels_for_pae_image+j)+1] = 0;
         image_data[4*(i*n_pixels_for_pae_image+j)+2] = 70;
         image_data[4*(i*n_pixels_for_pae_image+j)+3] = 255;
      }
   }

   unsigned int n_residues = pae_vecs.size();

   float max_value = get_max_value(pae_vecs);
   for (int i=0; i<n_pixels_for_pae_image; i++) {
      for (int j=0; j<n_pixels_for_pae_image; j++) {
         int idx = 4*(i*n_pixels_for_pae_image+j);
         float f_x = static_cast<float>(i)/static_cast<float>(n_pixels_for_pae_image);
         float f_y = static_cast<float>(j)/static_cast<float>(n_pixels_for_pae_image);
         int i_res_index = static_cast<int>(f_x * static_cast<float>(n_residues));
         int j_res_index = static_cast<int>(f_y * static_cast<float>(n_residues));
         auto col = value_to_colour(static_cast<float>(pae_vecs[i_res_index][j_res_index]), max_value);
         image_data[idx]   = col[0];
         image_data[idx+1] = col[1];
         image_data[idx+2] = col[2];
      }
   }

   cairo_t *cr;
   cairo_surface_t *pae_img_surface = cairo_image_surface_create_for_data(image_data, CAIRO_FORMAT_RGB24,
                                                                   n_pixels_for_pae_image, n_pixels_for_pae_image,
                                                                   n_pixels_for_pae_image*4);
   // the image is taller than it is for for the legend to fit
   cairo_surface_t *base_surface = cairo_image_surface_create(CAIRO_FORMAT_RGB24, n_pixels, n_pixels+100);

   cr = cairo_create(base_surface);
   cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
   cairo_paint(cr);

   if (cairo_surface_status(base_surface) == CAIRO_STATUS_SUCCESS) {

      // move the PAE graph in the base image
      int offset_x = 80;
      int offset_y = 20;

      // std::cout << "### cairo_surface_status() success " << std::endl;
      cairo_set_source_surface(cr, pae_img_surface, 80, 20);
      cairo_paint(cr);

      // draw a box around the pae plot
      cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
      cairo_rectangle(cr, offset_x, offset_y, n_pixels_for_pae_image, n_pixels_for_pae_image); // x, y, width, height
      cairo_stroke(cr);

      cairo_set_font_size(cr, offset_y);

      // "scored residue" label
      cairo_move_to(cr, 260, 573);
      cairo_show_text(cr, "Scored Residue");

      // "aligned residue" label
      cairo_move_to(cr, 26, 360);
      cairo_save(cr);
      cairo_rotate(cr, - M_PI / 2.0); //
      cairo_show_text(cr, "Aligned Residue");
      cairo_restore(cr);

      // "Expected position error" label
      cairo_move_to(cr, 150, 680);
      cairo_show_text(cr, "Expected position error (Ångströms)");

      // actual tick marks would be good here.

      // tick labels - x axis
      unsigned int tick_res_no = 0;
      while (tick_res_no < n_residues) {
         float f = static_cast<float>(tick_res_no) / static_cast<float>(n_residues);
         float pixel_for_tick_res_no = static_cast<float>(n_pixels_for_pae_image) * f;
         int x = offset_x + static_cast<int>(pixel_for_tick_res_no);
         int x_tweak = -20;
         if (tick_res_no == 0) x_tweak = -5;
         x += x_tweak;
         cairo_move_to(cr, x, 545);
         std::string text = std::to_string(tick_res_no);
         cairo_show_text(cr, text.c_str());
         tick_res_no += 100;
      }

      // tick labels - y axis
      tick_res_no = 0;
      while (tick_res_no < n_residues) {
         float f = static_cast<float>(tick_res_no) / static_cast<float>(n_residues);
         float pixel_for_tick_res_no = static_cast<float>(n_pixels_for_pae_image) * f;
         int y = offset_y + static_cast<int>(pixel_for_tick_res_no);
         y += 8; // so that the middle of the label text is at the tick position (rather than the bottom)
         double text_width = 42.0;
         if (tick_res_no == 0) text_width = 18.0;
         cairo_move_to(cr, offset_x-text_width, y);
         std::string text = std::to_string(tick_res_no);
         cairo_show_text(cr, text.c_str());
         tick_res_no += 100;
      }

      // legend colour ramp
      //
      for(double epe=0.0; epe<max_value; epe += 0.5) {
         auto col = value_to_colour(epe, max_value);
         double r = static_cast<double>(col[0])/255.0;
         double g = static_cast<double>(col[1])/255.0;
         double b = static_cast<double>(col[2])/255.0;
         cairo_set_source_rgb(cr, r, g, b);
         double f = epe/max_value;
         double x = offset_x + f * static_cast<double>(n_pixels_for_pae_image);
         double y = offset_y + n_pixels_for_pae_image + 80.0;
         cairo_rectangle(cr, x, y, 9.0, 20.0); // x, y, width, height
         cairo_fill(cr);
      }

      // draw a box around the legend
      cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
      float offset_y_legend = offset_y + n_pixels_for_pae_image + 80.0;
      cairo_rectangle(cr, offset_x, offset_y_legend, n_pixels_for_pae_image, 20.0); // x, y, width, height
      cairo_stroke(cr);

      // tick labels for the legend
      int tick_pos_error = 0;
      while (tick_pos_error < max_value) {
         float f = static_cast<float>(tick_pos_error) / static_cast<float>(max_value);
         float pixel_for_tick = static_cast<float>(n_pixels_for_pae_image) * f; // the legend is the same length (width)
         int x = offset_x + static_cast<int>(pixel_for_tick);
         int x_tweak = -20;
         if (tick_pos_error < 10) x_tweak = -7;
         x += x_tweak;
         cairo_move_to(cr, x, 645);
         std::string text = std::to_string(tick_pos_error);
         cairo_show_text(cr, text.c_str());
         tick_pos_error += 5;
      }

      // output
      cairo_surface_write_to_png_stream(base_surface, text_png_as_string_png_writer,
                                        reinterpret_cast<void *>(&s));

   } else {
      std::cout << "########### cairo_surface_status() fail " << std::endl;
   }
   return s;

#else

   return "no-cairo";

#endif

}
