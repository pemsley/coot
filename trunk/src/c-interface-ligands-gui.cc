/* src/c-interface-ligands.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Author: Paul Emsley
 * Copyright 2008, 2009 The University of Oxford
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */
#include <stdlib.h>
#include <iostream>
#include <stdexcept>

#if defined _MSC_VER
#include <windows.h>
#endif
 
#include "globjects.h" //includes gtk/gtk.h

#include "callbacks.h"
#include "interface.h" // now that we are moving callback
		       // functionality to the file, we need this
		       // header since some of the callbacks call
		       // fuctions built by glade.
#include <vector>
#include <string>

#include "mmdb_manager.h"
#include "mmdb-extras.h"
#include "mmdb.h"
#include "mmdb-crystal.h"

#include "graphics-info.h"
#include "c-interface.h"
#include "cc-interface.hh"
#include "coot-coord-utils.hh"
#include "peak-search.hh"

#ifdef USE_PYTHON
#include "Python.h"
#endif
#include "wligand.hh"

#include "guile-fixups.h"

/* in here we check if libcheck is available (if scripting is available) */
GtkWidget *wrapped_create_libcheck_monomer_dialog() {

   GtkWidget *w = create_libcheck_monomer_dialog();

#ifdef USE_GUILE

   std::string c = "(command-in-path? libcheck-exe)";
   SCM v = safe_scheme_command(c.c_str());

   if (!scm_is_true(v)) {
      GtkWidget *l = lookup_widget(w, "no_libcheck_frame");
      if (l) {
	 gtk_widget_show(l);
      }
   } 
   
#endif // USE_GUILE

#ifdef USE_PYTHON

   // something similar here.
   
#endif // USE_PYTHON   
   
   return w;
} 
