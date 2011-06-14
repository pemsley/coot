/* src/c-interface-mmdb.hh
 * 
 * Copyright 2007 The University of York
 * Author: Paul Emsley
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


// A misnamed file, perhaps.  This is for mmdb<->SCM interface.

#include "mmdb_manager.h"

#ifdef __cplusplus
#ifdef USE_GUILE
#include <libguile.h>	
#include <guile/gh.h> // needed for guile-1.6.x

// return 0 on failure.  An atom with an atom name of "" and element
// of "" is converted to a Ter atom.
// 
CMMDBManager *
mmdb_manager_from_scheme_expression(SCM molecule_expression);
SCM display_scm(SCM o);

#endif // USE_GUILE

#ifdef USE_PYTHON
#include "Python.h"
CMMDBManager * mmdb_manager_from_python_expression(PyObject *molecule_expression);
#endif // PYTHON
#endif 

// Delete the CONECT records from the MMDBManager? Depending on 
// write_conect_records_flag
void mmdb_manager_delete_conect(CMMDBManager *mol);
