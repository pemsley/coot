/* src/c-interface-refmac.h
 * 
 * Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Copyright 2007 by Paul Emsley
 * Copyright 2007, 2008, 2009, 2010, 2011, 2012 by The University of Oxford
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

/* svn $Id: c-interface.h 1458 2007-01-26 20:20:18Z emsley $ */

/*! \file 
  \brief Coot Scripting Interface - Refmac interface

*/

#ifndef C_INTERFACE_REFMAC_H
#define C_INTERFACE_REFMAC_H

/*
  The following extern stuff here because we want to return the
  filename from the file entry box.  That code (e.g.) 
  on_ok_button_coordinates_clicked (callback.c), is written and
  compiled in c.
 
  But, we need that function to set the filename in mol_info, which 
  is a c++ class.
 
  So we need to have this function external for c++ linking.
 
*/

/* Francois says move this up here so that things don't get wrapped
   twice in C-declarations inside gmp library. Hmm! */
#ifdef __cplusplus
#ifdef USE_GUILE
#include <cstdio> /* for std::FILE in gmp.h for libguile.h */
#include <libguile.h>		/* for SCM type (returned by safe_scheme_command) */
#else
#include <string> /* for std::string; included (sic!) in above for guile */
#endif /*  USE_GUILE */
#endif /* c++ */

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

#ifndef COOT_SCHEME_DIR
#define COOT_SCHEME_DIR "COOT_SCHEME_DIR"
#endif
#ifndef COOT_PYTHON_DIR
#define COOT_PYTHON_DIR "COOT_PYTHON_DIR"
#endif

int get_refmac_refinement_method();
void set_refmac_refinement_method(int method);
int refmac_imol_coords(void);

END_C_DECLS

#endif /* C_INTERFACE_REFMAC_H */
