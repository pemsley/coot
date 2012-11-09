
#include <stdio.h>
#include <string.h>

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <map>

#include "coot-coord-utils.hh"
#include "coot-utils.hh"

#include "guile-fixups.h"
#include "c-interface-ligands-swig.hh"

#include "probe-clash-score.hh"

coot::probe_clash_score_t::probe_clash_score_t(const std::string &dots_file_name) {

   // consider ignoring water contacts and atoms with occ < 0.3

   float score = 0;
   bool success = false;

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

      std::map<std::pair<coot::probe_atom_spec_t, coot::probe_atom_spec_t>, bool> bad_overlaps;
      std::map<std::pair<coot::probe_atom_spec_t, coot::probe_atom_spec_t>, bool> small_overlaps;
      std::map<std::pair<coot::probe_atom_spec_t, coot::probe_atom_spec_t>, bool> hydrogen_bonds;
      std::map<std::pair<coot::probe_atom_spec_t, coot::probe_atom_spec_t>, bool> wide_contact;
      std::map<std::pair<coot::probe_atom_spec_t, coot::probe_atom_spec_t>, bool> close_contact;
      

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

	       if (! c_type) {
		  std::cout << "null c_type!" << std::endl;
	       } else { 

		  if (offset > 0) { 
		     // std::cout << "scanning |" << contact_type1+offset << std::endl;
		     
		     if (sscanf((contact_type1+offset), "%f:%f:%f:%f:%f", 
				&x4, &x5, &x6, &factor5, &factor6)) {
			
			coot::probe_atom_spec_t as_1(atom_id1);
			coot::probe_atom_spec_t as_2(atom_id2);

			if (!as_1.empty() && !as_2.empty()) {
			   std::string contact_type(c_type);
			
			   std::pair<coot::probe_atom_spec_t, coot::probe_atom_spec_t> contact(as_1, as_2);

			   if (contact_type == "bo")
			      bad_overlaps[contact]   = true;
			   if (contact_type == "so")
			      small_overlaps[contact] = true;
			   if (contact_type == "hb")
			      hydrogen_bonds[contact] = true;
			   if (contact_type == "wc")
			      wide_contact[contact]   = true;
			   if (contact_type == "cc")
			      close_contact[contact]  = true;
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
      //
      // I don't understand how to fix this:
      // bad_overlaps.erase(std::remove_if(bad_overlaps.begin(), bad_overlaps.end(),
      // coot::spec_eraser(hydrogen_bonds)),
      // bad_overlaps.end());

      // do it long-hand:
      std::map<std::pair<coot::probe_atom_spec_t, coot::probe_atom_spec_t>, bool>
	 bad_overlaps_copy = bad_overlaps;
      
      std::map<std::pair<coot::probe_atom_spec_t, coot::probe_atom_spec_t>, bool>::const_iterator it;
      std::map<std::pair<coot::probe_atom_spec_t, coot::probe_atom_spec_t>, bool>::const_iterator ithb;
      for (it=bad_overlaps_copy.begin();
	   it!=bad_overlaps_copy.end();
	   it++) {
	 bool found = false;
	 for (ithb=hydrogen_bonds.begin();
	      ithb!=hydrogen_bonds.end();
	      ithb++) {
	    if (ithb->first == it->first) {
	       found = true;
	       break;
	    }
	 }
	 if (! found)
	    bad_overlaps[it->first] = true;
      }

      
      std::cout << "found " << bad_overlaps.size()   << " bad overlaps " << std::endl;
      std::cout << "found " << hydrogen_bonds.size() << " hydrogen bonds " << std::endl;
      std::cout << "found " << small_overlaps.size() << " small overlaps " << std::endl;
      std::cout << "found " << close_contact.size()  << " close contact " << std::endl;
      std::cout << "found " << wide_contact.size()   << " wide contacts " << std::endl;
      score =
	 float(hydrogen_bonds.size()) - float(bad_overlaps.size()) +
	 float(wide_contact.size()) * 0.03 + float(close_contact.size()) * 0.01;
	 
      success = true;
   }


}

#ifdef USE_GUILE
// SCM
coot::probe_clash_score_t
probe_clash_score(const std::string &dots_file_name) {

   SCM r = SCM_BOOL_F;
   coot::probe_clash_score_t p(dots_file_name);

   if (p.filled) {
      int n_bad_overlaps;
      int n_hydrogen_bonds;
      int n_small_overlaps;
      int n_close_contacts;
      int n_wide_contacts;

      r = scm_list_5(SCM_MAKINUM(p.n_bad_overlaps),
		     SCM_MAKINUM(p.n_hydrogen_bonds),
		     SCM_MAKINUM(p.n_small_overlaps),
		     SCM_MAKINUM(p.n_close_contacts),
		     SCM_MAKINUM(p.n_wide_contacts));
   }

   // return r;
   return p;
}

#endif // USE_GUILE
