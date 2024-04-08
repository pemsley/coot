/*
 * src/gtk-widget-conversion-utils.cc
 *
 * Copyright 2011 by University of York
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

#include <stdexcept>
#include <iostream>

#include <gtk/gtk.h>
#include "gtk-widget-conversion-utils.h"

#include "utils/coot-utils.hh"

struct entry_info_t coot_entry_to_val(GtkEntry *entry) { 

  struct entry_info_t ei;
  const gchar *text = gtk_editable_get_text(GTK_EDITABLE(entry));
  ei.float_is_set = 0;
  ei.val_as_float = 0;
  ei.val = 0;

  if (text) {
     ei.string_is_set = 1;
     ei.string = text;
     try {
	ei.val_as_float = coot::util::string_to_float(text);
     }
     catch (const std::runtime_error &rte) {
	std::cout << rte.what() << std::endl;
     } 
  } 

  return ei;

}

