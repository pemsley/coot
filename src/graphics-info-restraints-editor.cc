/* src/graphics-info.h
 * -*-c++-*-
 *
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2007 by Paul Emsley
 * Copyright 2007, 2008 by The University of Oxford
 * Copyright 2016 by Medical Research Council
 *
 * Author: Paul Emsley
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
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

#include "graphics-info.h"

#include "utils/logging.hh"
extern logging logger;

coot::restraints_editor
graphics_info_t::get_restraints_editor(GtkWidget *button) {

   coot::restraints_editor r; // a null/unset restraints editor
   int found_index = -1;
   bool debug = true;

   // set in coot::restraints_editor::fill_dialog(const coot::dictionary_residue_restraints_t &restraints)
   //
   GtkWidget *dialog = GTK_WIDGET(g_object_get_data(G_OBJECT(button), "restraints_editor_dialog"));
   std::cout << ":::::: get_restraints_editor() says dialog is " << dialog << std::endl;

   if (debug) {
      for (unsigned int i=0; i<restraints_editors.size(); i++) {
	 if (restraints_editors[i].is_valid())
	    std::cout << " debug:: in get_restraints_editor() a stored restraints editor number "
		      << i << " of " << restraints_editors.size() << ": "
		      << restraints_editors[i].get_dialog()
		      << " (c.f. " << button << ")" << std::endl;
	 else
	    std::cout << " debug:: in get_restraints_editor() a stored restraints editor number "
		      << i << " of " << restraints_editors.size() << ": "
		      << "NULL" << std::endl;
      }
   }

   for (unsigned int i=0; i<restraints_editors.size(); i++) {
      if (restraints_editors[i].is_valid()) {
	 if (restraints_editors[i].matches_dialog(dialog)) {
	    found_index = i;
	    break;
	 }
      }
   }

   if (found_index != -1)
      r = restraints_editors[found_index];
   return r;
}

void
graphics_info_t::clear_restraints_editor_by_dialog(GtkWidget *dialog) {

   for (unsigned int i=0; i<restraints_editors.size(); i++) {
      if (restraints_editors[i].is_valid()) {
	 if (restraints_editors[i].matches_dialog(dialog)) {
	    coot::restraints_editor null_restraints;
	    restraints_editors[i] = null_restraints;
	 }
      }
   }
}
