/* src/probe-clash-score.cc
 * 
 * Copyright 2010 by the University of Oxford
 * Copyright 2012 by Medical Research Council
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

#ifdef USE_PYTHON
#include <Python.h>
#endif // USE_PYTHON

#include <stdio.h>
#include <string.h>

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <map>

#include "coot-utils/coot-coord-utils.hh"
#include "utils/coot-utils.hh"

#include "guile-fixups.h"
#include "c-interface-ligands-swig.hh"

#include "probe-clash-score.hh"

coot::probe_clash_score_t::probe_clash_score_t(const std::string &dots_file_name) {

   // consider ignoring water contacts and atoms with occ < 0.3

   float score = 0;
   bool success = false;

   filled = false;
   n_bad_overlaps   = 0;
   n_hydrogen_bonds = 0;
   n_small_overlaps = 0;
   n_close_contacts = 0;
   n_wide_contacts  = 0;


   FILE* dots = fopen(dots_file_name.c_str(), "r" );
   if ( dots == NULL ) { 
      std::cout << "probe_clash_score()  - Could not read: " << dots_file_name << std::endl;
   } else {
      int n_input_lines = 0; 
      int n_lines = 0;
      int n_points = 0;
      std::string current_colour = "blue"; // should be reset.

      int imol_source, imol_target;
      float gap1, factor3, gap2, factor4, factor5, factor6;
      char line[240];
      char c_type[20];
      char atom_id1[30];
      char atom_id2[30];
      char contact_type1[200];
      char contact_type2[3];
      float x1, x2, x3, x4, x5, x6;
      std::string current_contact_type = "none";
      int dot_size = 2; // 3 for hydrogen bonds

      one_way_probe_contact_container_t bad_overlaps_cont;
      one_way_probe_contact_container_t small_overlaps_cont;
      one_way_probe_contact_container_t hydrogen_bonds_cont;
      one_way_probe_contact_container_t wide_contact_cont;
      one_way_probe_contact_container_t close_contact_cont;

      // null string arrays
      for (int i=0; i<20; i++)   c_type[i] = 0;
      for (int i=0; i<30; i++) atom_id1[i] = 0;
      for (int i=0; i<30; i++) atom_id2[i] = 0;
      std::string current_useful_name;

      std::string scan_line = ":%d->%d:%2c:%16c:%16c:%f:%f:%f:%f:%f:%f:%f:%s";
	 
      while (fgets( line, 240, dots ) != NULL) {
	 n_input_lines++; 
	 for (int i=0; i<3; i++) contact_type1[i] = 0;
	 for (int i=0; i<3; i++) contact_type2[i] = 0;

	 if (sscanf(line,
		    scan_line.c_str(),
		    &imol_source, &imol_target, c_type, atom_id1, atom_id2,
		    &gap1, &gap2, &x1, &x2, &x3, &factor3, &factor4,
		    contact_type1)) {

	    if (strlen(contact_type1) > 2) {
	       int colon_count = 0;
	       int offset = 0;

	       // std::cout << "colon search |" << contact_type1 << std::endl;
	       for (int i=0; i<10; i++) {
		  if (contact_type1[i] == ':')
		     colon_count++;
		  if (colon_count == 2) {
		     offset = i + 1;
		     break;
		  }
	       }

	       if (1) {
		  if (offset > 0) { 
		     // std::cout << "scanning |" << contact_type1+offset << std::endl;
		     
		     if (sscanf((contact_type1+offset), "%f:%f:%f:%f:%f", 
				&x4, &x5, &x6, &factor5, &factor6)) {
			
			coot::probe_atom_spec_t as_1(atom_id1);
			coot::probe_atom_spec_t as_2(atom_id2);

			if (!as_1.empty() && !as_2.empty()) {
			   std::string contact_type(c_type);
			
			   std::pair<coot::probe_atom_spec_t, coot::probe_atom_spec_t> contact(as_1, as_2);
                           filled = true;

			   if (contact_type == "wc")
			      wide_contact_cont.add(as_1, as_2);
			   if (contact_type == "bo")
			      bad_overlaps_cont.add(as_1, as_2);
			   if (contact_type == "so")
			      small_overlaps_cont.add(as_1, as_2);
			   if (contact_type == "hb")
			      hydrogen_bonds_cont.add(as_1, as_2);
			   if (contact_type == "cc")
			      close_contact_cont.add(as_1, as_2);
			   
			} else {
			   std::cout << "atom spec unset..." << std::endl;
			}
		     }
		  }
	       }
	    }
	 }
      }

      // now remove from bad_overlaps, contacts that are in hydrogen_bonds
      
      
      std::cout << "found  " << bad_overlaps_cont.size()   << " bad overlaps "   << std::endl;
      std::cout << "found  " << hydrogen_bonds_cont.size() << " hydrogen bonds " << std::endl;
      std::cout << "found  " << small_overlaps_cont.size() << " small overlaps " << std::endl;
      std::cout << "found  " << close_contact_cont.size()  << " close contact "  << std::endl;
      std::cout << "found  " << wide_contact_cont.size()   << " wide contacts "  << std::endl;
      score =
	 float(hydrogen_bonds_cont.size()) - float(bad_overlaps_cont.size()) +
	 float(wide_contact_cont.size()) * 0.03 + float(close_contact_cont.size()) * 0.01;
	 
      n_bad_overlaps   = bad_overlaps_cont.size();
      n_hydrogen_bonds = hydrogen_bonds_cont.size();
      n_small_overlaps = small_overlaps_cont.size();
      n_close_contacts = close_contact_cont.size();
      n_wide_contacts  = wide_contact_cont.size();
      success = true;

   }

}

coot::probe_clash_score_t
probe_clash_score(const std::string &dots_file_name) {

   coot::probe_clash_score_t p(dots_file_name);
   return p;
}

#ifdef USE_GUILE
// SCM
SCM 
probe_clash_score_scm(const std::string &dots_file_name) { 

   coot::probe_clash_score_t p(dots_file_name);
   return probe_clash_score_as_scm(p);
}

SCM probe_clash_score_as_scm(const coot::probe_clash_score_t &p) {

   SCM r = SCM_BOOL_F;
   if (p.filled) {
      r = scm_list_5(SCM_MAKINUM(p.n_bad_overlaps),
		     SCM_MAKINUM(p.n_hydrogen_bonds),
		     SCM_MAKINUM(p.n_small_overlaps),
		     SCM_MAKINUM(p.n_close_contacts),
		     SCM_MAKINUM(p.n_wide_contacts));
   }
   return r;
}

// p is a list of 5 ints.  Convert that to a probe_clash_score_t
// 
coot::probe_clash_score_t
probe_clash_score_from_scm(SCM p) {

   coot::probe_clash_score_t pcs;
   std::cout << "debug:: probe_clash_score_from_scm() here 1 " << std::endl;
   if (scm_is_true(scm_list_p(p))) {
      SCM p_len_scm = scm_length(p);
      int p_len = scm_to_int(p_len_scm);
      std::cout << "debug:: probe_clash_score_from_scm() here 2 " << p_len << std::endl;
      if (p_len == 5) {
	 SCM n_bo_scm = scm_list_ref(p, SCM_MAKINUM(0));
	 SCM n_hb_scm = scm_list_ref(p, SCM_MAKINUM(1));
	 SCM n_so_scm = scm_list_ref(p, SCM_MAKINUM(2));
	 SCM n_cc_scm = scm_list_ref(p, SCM_MAKINUM(3));
	 SCM n_wc_scm = scm_list_ref(p, SCM_MAKINUM(4));
	 int n_bo = scm_to_int(n_bo_scm);
	 int n_hb = scm_to_int(n_hb_scm);
	 int n_so = scm_to_int(n_so_scm);
	 int n_cc = scm_to_int(n_cc_scm);
	 int n_wc = scm_to_int(n_cc_scm);
	 std::cout << "debug:: probe_clash_score_from_scm() here 3 " << n_bo << std::endl;
	 pcs = coot::probe_clash_score_t(n_bo, n_hb, n_so, n_cc, n_wc);
      }
   }
   return pcs;
}

#endif // USE_GUILE

#ifdef USE_PYTHON
// Python
PyObject * 
probe_clash_score_py(const std::string &dots_file_name) { 

   coot::probe_clash_score_t p(dots_file_name); 
   return probe_clash_score_as_py(p);
} 

PyObject *probe_clash_score_as_py(const coot::probe_clash_score_t &p) {

   PyObject *r = Py_False;
   if (p.filled) {
      r = PyList_New(5);
      PyList_SetItem(r, 0, PyLong_FromLong(p.n_bad_overlaps));
      PyList_SetItem(r, 1, PyLong_FromLong(p.n_hydrogen_bonds));
      PyList_SetItem(r, 2, PyLong_FromLong(p.n_small_overlaps));
      PyList_SetItem(r, 3, PyLong_FromLong(p.n_close_contacts));
      PyList_SetItem(r, 4, PyLong_FromLong(p.n_wide_contacts));
   }
   if (PyBool_Check(r)) {
      Py_XINCREF(r);
   }
   return r;
}

// p is a list of 5 ints.  Convert that to a probe_clash_score_t
// 
coot::probe_clash_score_t
probe_clash_score_from_py(PyObject *p) {

   coot::probe_clash_score_t pcs;
   std::cout << "debug:: probe_clash_score_from_py() here 1 " << std::endl;
   if (PyList_Check(p)) {
      Py_ssize_t p_len = PyList_Size(p);
      std::cout << "debug:: probe_clash_score_from_py() here 2 " << p_len << std::endl;
      if (p_len == 5) {
         PyObject *n_bo_py = PyList_GetItem(p, 0);
         PyObject *n_hb_py = PyList_GetItem(p, 1);
         PyObject *n_so_py = PyList_GetItem(p, 2);
         PyObject *n_cc_py = PyList_GetItem(p, 3);
         PyObject *n_wc_py = PyList_GetItem(p, 4);
         int n_bo = PyLong_AsLong(n_bo_py);
         int n_hb = PyLong_AsLong(n_hb_py);
         int n_so = PyLong_AsLong(n_so_py);
         int n_cc = PyLong_AsLong(n_cc_py);
         int n_wc = PyLong_AsLong(n_wc_py);
         std::cout << "debug:: probe_clash_score_from_py() here 3 " << n_bo << std::endl;
         pcs = coot::probe_clash_score_t(n_bo, n_hb, n_so, n_cc, n_wc);
      }
   }
   return pcs;
}

#endif // USE_PYTHON
