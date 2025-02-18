/*
 * src/cc-interface-png.cc
 *
 * Copyright 2017 by Medical Research Council
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

#ifdef USE_PYTHON
#include <Python.h>
#endif

#ifdef HAVE_CXX11
#include <cmath>
#else
#include <math.h>
#endif

#include <iostream>


#ifdef USE_GUILE
#include <libguile.h>
#endif

#include <cairo.h>

// All the above header includes (other than cairo) are needed to include this include
#include "cc-interface-image.hh" // include for text_png_as_string

#include "utils/coot-utils.hh" // colour_holder
#include "geometry/residue-and-atom-specs.hh"

// for spec conversion atom_spec_from_python_expression
#include "cc-interface.hh"
// for valide model molecule check (so that we can look up the residue type)
#include "graphics-info.h"

#ifdef USE_PYTHON


// cairo string helper function
// we don't need this in a header if it's abouve text_png_as_string
cairo_status_t
text_png_as_string_png_stream_writer(void *closure_in,
				     const unsigned char *data,
				     unsigned int length) {

   std::string *s_ptr = static_cast<std::string *>(closure_in);
   *s_ptr += std::string(reinterpret_cast<const char *>(data), length); // it's safe!
   return CAIRO_STATUS_SUCCESS;
}


// make a png given the info in text_info_dict
// fields: text, font, size, colour
std::string text_png_as_string(PyObject *text_info_dict_py) {

   std::string r;
   r.reserve(12000);

   std::pair<bool, long> npx(false, 300);
   std::pair<bool, long> npy(false, 100);
   coot::colour_holder bg_col(0.95, 0.95, 0.95);
   coot::colour_holder fg_col(0.2, 0.2, 0.2);
   std::string text = "Hello, World";
   std::string face = "sans";
   float font_size = 24;
   cairo_font_slant_t slant = CAIRO_FONT_SLANT_NORMAL;
   cairo_font_weight_t weight = CAIRO_FONT_WEIGHT_NORMAL;
   bool translucent = false;

   // extract values from the dictionary
   PyObject *key;
   PyObject *value;
   Py_ssize_t pos = 0;
   if (PyDict_Check(text_info_dict_py)) {
      while (PyDict_Next(text_info_dict_py, &pos, &key, &value)) {
        std::string key_string = PyBytes_AS_STRING(PyUnicode_AsUTF8String(key));
	 if (key_string == "text") {
	    if (PyUnicode_Check(value)) {
              std::string s = PyBytes_AS_STRING(PyUnicode_AsUTF8String(value));
	       text = s;
	    }
	 }
	 if (key_string == "face") {
	    if (PyUnicode_Check(value)) {
              std::string s = PyBytes_AS_STRING(PyUnicode_AsUTF8String(value));
	       face = s;
	    }
	 }
	 if (key_string == "font-size") {
	    if (PyFloat_Check(value) || PyLong_Check(value)) {
	       double fs = PyFloat_AsDouble(value);
	       font_size = fs;
	       std::cout << "debug:: font_size set to " << font_size << std::endl;
	    } else {
	       std::cout << "font-size: value is not an int or a float" << std::endl;
	    }
	 }
	 
	 if (key_string == "fg-colour") {
	    if (PyUnicode_Check(value)) {
              std::string s = PyBytes_AS_STRING(PyUnicode_AsUTF8String(value));
	       if (s.length() == 7)
		  if (s[0] == '#')
		     fg_col = coot::colour_holder(s);
	    }
	 }
	 if (key_string == "bg-colour") {
	    if (PyUnicode_Check(value)) {
              std::string s = PyBytes_AS_STRING(PyUnicode_AsUTF8String(value));
	       if (s.length() == 7)
		  if (s[0] == '#')
		     bg_col = coot::colour_holder(s);
	    }
	 }
	 if (key_string == "translucent") {
	    if (PyBool_Check(value)) {
	       // This is how you test python booleans
	       if (value == Py_False)
		  translucent = false; // default
	       if (value == Py_True)
		  translucent = true;
	    }
	 }
	 if (key_string == "npx") {
	    if (PyLong_Check(value)) {
	       npx.first = true;
	       npx.second = PyLong_AsLong(value);
	    }
	 }
	 if (key_string == "npy") {
	    if (PyLong_Check(value)) {
	       npy.first = true;
	       npy.second = PyLong_AsLong(value);
	    }
	 }
	 if (key_string == "atom-spec") {
	    coot::atom_spec_t atom_spec = atom_spec_from_python_expression(value);
	    int imol = atom_spec.int_user_data;
	    if (graphics_info_t::is_valid_model_molecule(imol)) {
	       coot::residue_spec_t residue_spec(atom_spec);
	       mmdb::Residue *residue_p =
		  graphics_info_t::molecules[imol].get_residue(residue_spec);
	       if (residue_p) {
		  std::string rn = residue_p->GetResName();
		  text = atom_spec.label(rn);
	       } else {
		  text = atom_spec.label();
	       }
	    } else {
	       text = atom_spec.label();
	    }
	 }
	 if (key_string == "residue-spec") {
	    coot::residue_spec_t spec = residue_spec_from_py(value);
	    text = spec.label();
	 }
	 if (key_string == "residue-spec") {
	    coot::residue_spec_t residue_spec = residue_spec_from_py(value);
	    int imol = residue_spec.int_user_data;
	    if (graphics_info_t::is_valid_model_molecule(imol)) {
	       mmdb::Residue *residue_p =
		  graphics_info_t::molecules[imol].get_residue(residue_spec);
	       if (residue_p) {
		  std::string rn = residue_p->GetResName();
		  text = residue_spec.label(rn);
	       } else {
		  text = residue_spec.label();
	       }
	    } else {
	       text = residue_spec.label();
	    }
	 }
      }
   }

   if (npx.first == false) {
      npx.second = text.length() * font_size * 0.7; // fudge factor
   }

   if (npy.first == false)
      npy.second = std::lround(std::floor(float(font_size)*1.3)) + 1; // fudge factors
   
   cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, npx.second, npy.second);
   cairo_t *cr = cairo_create(surface);

   // Drawing space is 0->1 with 0,0 top left
   // if cairo_image_surface_create called with 240,240, then cairo_scale called with
   // 120,120 puts the image, half-size in top left quadrant
   // cairo_scale(cr, npx, npx);

   cairo_select_font_face(cr, face.c_str(), slant, weight);
   cairo_set_font_size(cr, font_size);

   cairo_text_extents_t te;
   cairo_text_extents(cr, text.c_str(), &te);
   // can we use te.width and height to change the surface?
   // we can create a fake surface I suppose and a fake cairo_t cr_fake so that
   // we can ask what the text extents are - and *then* create surface and cr.
   // Hmm... Well, seems like a lot of work that may well not be needed.

   // background
   if (translucent) {
      cairo_save(cr);
      cairo_set_source_rgba(cr, 0.1, 0.1, 0.1, 0.0);
      cairo_set_operator(cr, CAIRO_OPERATOR_SOURCE);
      cairo_paint(cr);
      cairo_restore(cr);
   } else {
      cairo_set_source_rgb(cr, bg_col.red, bg_col.green, bg_col.blue);
      cairo_paint(cr);
   }

   cairo_move_to(cr, 0.0, font_size); // seems reaonable, puts text in top left of box
   cairo_set_source_rgb(cr, fg_col.red, fg_col.green, fg_col.blue);

   // in future, use pango-cairo:
   cairo_show_text(cr, text.c_str());

   cairo_surface_write_to_png_stream(surface, text_png_as_string_png_stream_writer,
				     reinterpret_cast<void *> (&r));
   cairo_destroy(cr);
   cairo_surface_destroy(surface);

   return r;
}

#endif // USE_PYTHON


void display_svg_from_file_in_a_dialog(const std::string &file_name) {

   auto file_to_string = [] (const std::string &file_name) {
      std::string s;
      std::string line;
      std::ifstream f(file_name);
      if (!f) {
         std::cout << "Failed to open " << file_name << std::endl;
      } else {
         while (std::getline(f, line)) {
            s += line;
            s += "\n";
         }
      }
      return s;
   };

   if (coot::file_exists(file_name)) {
      std::string s = file_to_string(file_name);
      display_svg_from_string_in_a_dialog(s);
   }

}


//  ------------------------- svg (not png) ------------------------------

#ifdef HAVE_RSVG
#include <librsvg/rsvg.h>
#endif

void display_svg_from_string_in_a_dialog(const std::string &string) {

   // the api function declaration exists

#ifdef HAVE_RSVG
   // Load the SVG file
   GError *error = NULL;
   RsvgHandle *handle = rsvg_handle_new_from_data((const unsigned char *)string.c_str(),
                                                  string.length(), &error);

   if (handle == NULL) {
      g_printerr("Error loading SVG: %s\n", error->message);
      g_error_free(error);
      return;
   } else {
      GtkWidget *window = gtk_window_new();
      gtk_window_set_title(GTK_WINDOW(window), "Coot: SVG Viewer");
      // Create a drawing area
      GtkWidget *drawing_area = gtk_drawing_area_new();
      gtk_widget_set_hexpand(drawing_area, TRUE); // Allow expansion in x direction
      gtk_widget_set_vexpand(drawing_area, TRUE); // Allow expansion in y direction
      gtk_window_set_child(GTK_WINDOW(window), drawing_area);

      // Connect the draw signal to our drawing function, passing the SVG handle

      // void rsvg_handle_get_intrinsic_dimensions (RsvgHandle* handle,
      //                                            gboolean* out_has_width,
      //                                            RsvgLength* out_width,
      //                                            gboolean* out_has_height,
      //                                            RsvgLength* out_height,
      //                                            gboolean* out_has_viewbox,
      //                                            RsvgRectangle* out_viewbox)


      auto draw_svg = +[] (GtkDrawingArea *da, cairo_t *cr, int w, int h, gpointer data) {

         GError *error = NULL;
         RsvgHandle *handle = static_cast<RsvgHandle *>(data);
         GdkRectangle allocation;
         gtk_widget_get_allocation(GTK_WIDGET(da), &allocation);
         // Scale the SVG to fit the widget (optional, adjust as needed)
         gboolean has_width, has_height, has_viewbox;
         RsvgLength width, height;
         RsvgRectangle viewbox;
         rsvg_handle_get_intrinsic_dimensions(handle, &has_width, &width, &has_height, &height,
                                              &has_viewbox, &viewbox);

         if (has_width && has_height) {

            double w_int = static_cast<double>(width.length);
            double h_int = static_cast<double>(height.length);
            std::cout << "intrinsic w unit " << width.unit << std::endl;
            double scale_x = static_cast<double>(allocation.width) / w_int;
            double scale_y = static_cast<double>(allocation.height) / h_int;
            double scale = std::min(scale_x, scale_y);
            std::cout << "intrinsic w and h " << w_int << " h_int " << h_int
                      << " window w and h " << allocation.width << " " << allocation.height
                      << std::endl;
            scale = static_cast<double>(allocation.width) / 60.0;
            cairo_scale(cr, scale, scale);
            cairo_translate(cr, 18.0, 18.0);
         }

         // Render the SVG
         rsvg_handle_render_document(handle, cr, &viewbox, &error);
      };

      gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(drawing_area), draw_svg, handle, NULL);
      gtk_widget_set_visible(window, TRUE);
   }
#endif // HAVE_LIBRSVG
}


