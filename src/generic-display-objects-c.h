/* src/generic-display-objects-c.h
 *
 * Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Copyright 2007 by Paul Emsley
 * Copyright 2007, 2008, 2009, 2010, 2011, 2012 by The University of Oxford
 * Copyright 2014, 2015, 2016 by Medical Research Council
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

#ifndef GENERIC_DISPLAY_OBJECTS_C_H
#define GENERIC_DISPLAY_OBJECTS_C_H

#ifndef BEGIN_C_DECLS

#ifdef __cplusplus
#define BEGIN_C_DECLS extern "C" {
#define END_C_DECLS }

#else
#define BEGIN_C_DECLS
#define END_C_DECLS     
#endif
#endif /* BEGIN_C_DECLS */

BEGIN_C_DECLS


/*! \file
  \brief Coot Scripting Interface - Generic objects
*/

void clear_generic_objects_dialog_pointer();

/*! \brief display (1) or undisplay (0) all generic display objects */
void set_display_all_generic_objects(int state);

/*! \brief a kludgey thing, so that the generic objects gui can be
  called from a callback.  */
void generic_objects_gui_wrapper();

/*! \brief close all generic display objects  */
void close_all_generic_objects();

END_C_DECLS

#endif /* GENERIC_DISPLAY_OBJECTS_C_H */
