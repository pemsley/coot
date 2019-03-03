/* src/probe-clash-score.hh
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

#ifndef PROBE_CLASH_SCORE_HH
#define PROBE_CLASH_SCORE_HH

#include <algorithm> // for find

#ifdef USE_GUILE
#include <libguile.h>     
#endif // USE_GUILE

#include "utils/coot-utils.hh"
#include "coot-utils/coot-coord-utils.hh"

/*! \file
  \brief Coot Scripting Interface - Probe clash interface
*/

namespace coot {

   //! probe atom type
   //
   // " B  72 CYS  HG A: B  53 HIS  H   "
   // 
   class probe_atom_spec_t : public atom_spec_t {
   public: 
      probe_atom_spec_t(const std::string &s) : atom_spec_t() {
	 if (s.length() > 14) { 
	    std::string chain_local = s.substr(0,2);
	    std::string res_no_str = s.substr(2, 4);
	    std::string atom_name_local = s.substr(11, 4);
	    try {
	       int resno_local = util::string_to_int(res_no_str);
	       if (chain_local[0] == ' ') { 
		  if (chain_local.length() > 1)
		     chain_id = std::string(chain_local.substr(1));
	       } else { 
		  chain_id = chain_local;
	       }
	       res_no = resno_local;
	       atom_name = atom_name_local;
	    }
	    catch (const std::exception &e) {
	       std::cout << "WARNING:: " << e.what() << std::endl;
	    }
	 }
      }
      probe_atom_spec_t() : atom_spec_t() {}
   };

   //! one way probe contact
   class one_way_probe_contact_t {
   public:
      probe_atom_spec_t from_atom;
      std::vector<probe_atom_spec_t> to_atoms;
      one_way_probe_contact_t(const probe_atom_spec_t &spec) {
	 from_atom = spec;
      }
      unsigned int size() const { return to_atoms.size(); }
      void add(const probe_atom_spec_t &spec) {
	 std::vector<probe_atom_spec_t>::const_iterator it;
	 it = std::find(to_atoms.begin(), to_atoms.end(), spec);
	 if (it == to_atoms.end()) {
	    to_atoms.push_back(spec);
	 }
      }
   };

   //! one way probe contacts container
   class one_way_probe_contact_container_t {
   public:
      std::vector<one_way_probe_contact_t> contacts;
      void add(const probe_atom_spec_t &from_atom,
	       const probe_atom_spec_t &to_atom) {
	 bool found = false;
	 for (unsigned int i=0; i<contacts.size(); i++) { 
	    if (contacts[i].from_atom == from_atom) {
	       contacts[i].add(to_atom);
	       found = true;
	    }
	 }
	 if (! found) {
	    one_way_probe_contact_t new_contact(from_atom);
	    new_contact.add(to_atom);
	    contacts.push_back(new_contact);
	 }
      }
      unsigned int size() const {
	 unsigned int s = 0;
	 for (unsigned int i=0; i<contacts.size(); i++)
	    s += contacts[i].size();
	 return s;
      }
   };

   //! probe clash score
   class probe_clash_score_t {
   public:
      bool filled;
      int n_bad_overlaps;
      int n_hydrogen_bonds;
      int n_small_overlaps;
      int n_close_contacts;
      int n_wide_contacts;
      probe_clash_score_t() {
	 filled = false;
      }
      probe_clash_score_t(int n_bad_overlaps_in,
			  int n_hydrogen_bonds_in,
			  int n_small_overlaps_in,
			  int n_close_contacts_in,
			  int n_wide_contacts_in) {

	 n_bad_overlaps   = n_bad_overlaps_in;
	 n_hydrogen_bonds = n_hydrogen_bonds_in;
	 n_small_overlaps = n_small_overlaps_in;
	 n_close_contacts = n_close_contacts_in;
	 n_wide_contacts  = n_wide_contacts_in;
	 filled = true;
      } 
      probe_clash_score_t(const std::string &dots_file_name);
   }; 


   //! This doesn't work
   // couldn't get this to work - so I did by another method.
   // deleteable cruft.
   class spec_eraser {
   public:
      std::map<std::pair<probe_atom_spec_t, probe_atom_spec_t>, bool> ref_specs;
      spec_eraser(const std::map<std::pair<probe_atom_spec_t, probe_atom_spec_t>, bool> &ref_specs_in) {
	 ref_specs = ref_specs_in;
      }
      // bool operator() (const std::pair<probe_atom_spec_t, probe_atom_spec_t> &s) const {
      //    return true;
      // }
   };
}


#ifdef USE_GUILE
//! \brief return scheme false on failure or a scheme list
// (n_bad_overlaps n_hydrogen_bonds n_small_overlaps n_close_contacts
// n_wide_contacts)
// 
SCM probe_clash_score_scm(const std::string &dots_file_name);
SCM probe_clash_score_as_scm(const coot::probe_clash_score_t &p);
coot::probe_clash_score_t probe_clash_score_from_scm(SCM p);
#endif // USE_GUILE

#ifdef USE_PYTHON
//! \brief return scheme false on failure or a scheme list
// (n_bad_overlaps n_hydrogen_bonds n_small_overlaps n_close_contacts
// n_wide_contacts)
// 
PyObject *probe_clash_score_py(const std::string &dots_file_name);
PyObject *probe_clash_score_as_py(const coot::probe_clash_score_t &p);
coot::probe_clash_score_t probe_clash_score_from_py(PyObject *p);
#endif // USE_PYTHON


#endif // PROBE_CLASH_SCORE_HH
