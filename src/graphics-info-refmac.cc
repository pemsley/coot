/* src/graphics-info.cc
 *
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by the University of York
 * Copyright 2007, 2008, 2009 by the University of Oxford
 * Copyright 2015 by Medical Research Council
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

// See also the below function, which should be used in future.
//
// c.f. the other function:
// graphics_info_t::fill_option_menu_with_map_options(GtkWidget *option_menu,
// 						   GtkSignalFunc signal_func,
//						   int imol_active_position).
//

// There is very little here worth saving.

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"


#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif

#include <gtk/gtk.h>  // must come after mmdb_manager on MacOS X Darwin
// #include <GL/glut.h>  // for some reason...  // Eh?

#include <iostream>
#include <dirent.h>   // for refmac dictionary files

#include <sys/types.h> // for stating
#include <sys/stat.h>

#if !defined _MSC_VER && !defined WINDOWS_MINGW
#include <unistd.h>
#else
//#include "coot-sysdep.h"
#endif

#include <mmdb2/mmdb_manager.h>
#include "coords/mmdb-extras.h"
#include "coords/mmdb.hh"
#include "coords/mmdb-crystal.h"
#include "coords/Cartesian.h"
#include "coords/Bond_lines.h"

#include "clipper/core/map_utils.h" // Map_stats
#include "skeleton/graphical_skel.h"

#include "graphics-info.h"
#include "interface.h"

#include "molecule-class-info.h"
#include "skeleton/BuildCas.h"

#include "gl-matrix.h" // for baton rotation
#include "trackball.h" // for baton rotation

#include "analysis/bfkurt.hh"

#include "globjects.h"
#include "ligand/ligand.hh"

#include "ligand/dunbrack.hh"

#include "utils/coot-utils.hh"

#include "cmtz-interface.hh"
#include "cmtz-interface-gui.hh"

#include "manipulation-modes.hh"

#include "guile-fixups.h"


int
graphics_info_t::fill_combobox_with_map_mtz_options(GtkWidget *combobox, GCallback signal_func,
						    int imol_active) {

   int imol = fill_combobox_with_map_options(combobox, signal_func, imol_active);
   return imol;
}

int
graphics_info_t::fill_combobox_with_map_options(GtkWidget *combobox,
						GCallback signal_func,
						int imol_preferred) {

   // delete this function on merge - hmm what did I mean by that?
   std::vector<int> maps_vec;
   for (int i=0; i<n_molecules(); i++)
      if (is_valid_map_molecule(i))
	 maps_vec.push_back(i);

   fill_combobox_with_molecule_options(combobox, signal_func, imol_preferred, maps_vec);

   // fill_combobox_with_molecule_options doesn't return a value (for the active molecule)
   // sadly, so set the presumed active imol if the imol_preferred is in the molecules in
   // maps_vec

   int imol = -1;
   if (std::find(maps_vec.begin(), maps_vec.end(), imol_preferred) != maps_vec.end())
      imol = imol_preferred;

   return imol;
}


void
graphics_info_t::fill_combobox_with_difference_map_options(GtkWidget *combobox,
							   GCallback signal_func,
							   int imol_active_position) {

   std::vector<int> maps_vec;
   for (int i=0; i<n_molecules(); i++) {
      if (molecules[i].is_difference_map_p())
	 maps_vec.push_back(i);
   }

   fill_combobox_with_molecule_options(combobox, signal_func, imol_active_position, maps_vec);

}



void
graphics_info_t::set_refmac_refinement_method(int method) {

  graphics_info_t g;

  switch (method) {

  case coot::refmac::RESTRAINED:
    g.refmac_refinement_method = coot::refmac::RESTRAINED;
    break;

  case coot::refmac::RIGID_BODY:
    g.refmac_refinement_method = coot::refmac::RIGID_BODY;
    break;

  case coot::refmac::RESTRAINED_TLS:
    g.refmac_refinement_method = coot::refmac::RESTRAINED_TLS;
    break;

  default:
    g.refmac_refinement_method = coot::refmac::RESTRAINED;
    break;
  }
}

void
graphics_info_t::refmac_change_refinement_method(GtkWidget *item, GtkPositionType pos) {

  set_refmac_refinement_method(pos);

}

