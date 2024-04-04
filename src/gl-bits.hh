/*
 * src/gl-bits.hh
 *
 * Copyright 2009 by University of York
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

#ifndef HAVE_GL_BITS_HH
#define HAVE_GL_BITS_HH

enum {GL_CONTEXT_MAIN = 0, GL_CONTEXT_SECONDARY = 1};

class gl_context_info_t {
public:
   GtkWidget *widget_1;
   GtkWidget *widget_2;

   gl_context_info_t(GtkWidget *widget_1_in, GtkWidget *widget_2_in) {
      widget_1 = widget_1_in;
      widget_2 = widget_2_in;
   } 
   gl_context_info_t() {
      widget_1 = 0;
      widget_2 = 0;
   }
};


#endif // HAVE_GL_BITS_HH
