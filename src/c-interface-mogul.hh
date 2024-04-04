/* src/c-interface-mogul.hh
 * 
 * Copyright 2011 by the University of Oxford
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

#ifndef C_INTERFACE_MOGUL_HH
#define C_INTERFACE_MOGUL_HH
#ifdef USE_GUILE
#include <cstdio> // for std::FILE in gmp.h for libguile.h
#include <libguile.h>		/* for SCM type (returned by safe_scheme_command) */
#endif /*  USE_GUILE */


void mogul_markup(int imol, const char *chain_id, int res_no, const char *ins_code,
		  const char *mogul_out_file_name);

int update_restraints_using_mogul(int imol, const char *chain_id, int res_no, const char *ins_code,
				  const char *residue_type, 
				  const char *mogul_out_file_name);

void set_mogul_max_badness(float b);
float get_mogul_max_badness();

#ifdef __cplusplus
#ifdef USE_GUILE
SCM mogul_results_scm(const char *mogul_out_file_name);

SCM mogul_results_process_many(const char *glob_str, const char *glob_dir);
#endif 

#ifdef USE_PYTHON
PyObject *mogul_results_py(const char *mogul_out_file_name);
#endif // USE_PYTHON
#endif 

#endif // C_INTERFACE_MOGUL_HH

