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
#define BEGIN_C_DECLS extern
#define END_C_DECLS
#endif
#endif /* BEGIN_C_DECLS */

BEGIN_C_DECLS

#define COOT_SCHEME_DIR "COOT_SCHEME_DIR"
#define COOT_PYTHON_DIR "COOT_PYTHON_DIR"

/*  this is the option menu callback - does nothing. */
int set_refmac_molecule(int imol); /* used by callback.c */

/* some methods to get refmac run parameters */
int get_refmac_refinement_method(void);
void set_refmac_refinement_method(int method);
int get_refmac_phase_input(void);
void set_refmac_phase_input(int phase_flag);
void set_refmac_use_tls(int state);
int refmac_use_tls_state(void);
void set_refmac_use_twin(int state);
int refmac_use_twin_state(void);
void set_refmac_use_sad(int state);
int refmac_use_sad_state(void);
int get_refmac_ncycles(void);
void set_refmac_ncycles(int no_cycles);
void add_refmac_ncycle_no(int cycle);
void set_refmac_use_ncs(int state);
int refmac_use_ncs_state(void);
void set_refmac_use_intensities(int state);
int refmac_use_intensities_state(void);
int refmac_imol_coords(void);
void add_refmac_sad_atom(const char *atom_name, float fp, float fpp, float lambda);
void add_refmac_sad_atom_fp(const char *atom_name, float fp, float fpp);
void add_refmac_sad_atom_lambda(const char *atom_name, float lambda);
void clear_refmac_sad_atoms();
short int get_refmac_used_mtz_file_state();
void set_refmac_used_mtz_file(int state);
const gchar *get_saved_refmac_file_filename(void);
void set_stored_refmac_file_mtz_filename(int imol, const char *mtz_filename);
void save_refmac_params_to_map(int imol_map,
			       const char *mtz_filename,
			       const char *fobs_col,
			       const char *sigfobs_col,
			       const char *r_free_col,
			       int r_free_flag_sensible);
void save_refmac_phase_params_to_map(int imol_map,
			     	     const char *phi,
				     const char *fom,
				     const char *hla,
				     const char *hlb,
				     const char *hlc,
				     const char *hld);

END_C_DECLS
#endif /* C_INTERFACE_REFMAC_H */
