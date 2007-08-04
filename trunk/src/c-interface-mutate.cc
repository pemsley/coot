/* src/c-interface-mutate.cc
 * 
 * Copyright 2007 by Paul Emsley
 * Copyright 2007 by The University of York
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
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

#include <stdlib.h>
#include <iostream>

#ifdef USE_GUILE
#include <guile/gh.h>

#if (SCM_MAJOR_VERSION > 1) || (SCM_MINOR_VERSION > 7)
// no fix up needed 
#else    
#define scm_to_int gh_scm2int
#define scm_to_locale_string SCM_STRING_CHARS
#define scm_to_double  gh_scm2double
#define  scm_is_true gh_scm2bool
#endif // SCM version

#endif // USE_GUILE

#ifdef USE_PYTHON
#include "Python.h"
#endif // USE_PYTHON

#include <gtk/gtk.h>

#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>
#include "cootaneer-sequence.h"

#include "c-interface.h"
#include "cc-interface.hh"
#include "coot-utils.hh"

#include "graphics-info.h"

#ifdef USE_GUILE
coot::atom_spec_t
atom_spec_from_scm_expression(SCM expr) {

   coot::atom_spec_t atom_spec;
   atom_spec.string_user_data = "Bad Spec";

   SCM len_expr = scm_length(expr);
   int len_view = scm_to_int(len_expr);
   if (len_view == 5) {
      
      SCM idx_scm = SCM_MAKINUM(0);
      SCM chain_id_scm = scm_list_ref(expr, idx_scm);
      std::string chain_id = scm_to_locale_string(chain_id_scm);

      idx_scm = SCM_MAKINUM(1);
      SCM resno_scm = scm_list_ref(expr, idx_scm);
      int resno = scm_to_int(resno_scm);

      idx_scm = SCM_MAKINUM(2);
      SCM ins_code_scm = scm_list_ref(expr, idx_scm);
      std::string ins_code = scm_to_locale_string(ins_code_scm);
      
      idx_scm = SCM_MAKINUM(3);
      SCM atom_name_scm = scm_list_ref(expr, idx_scm);
      std::string atom_name = scm_to_locale_string(atom_name_scm);
      
      idx_scm = SCM_MAKINUM(4);
      SCM alt_conf_scm = scm_list_ref(expr, idx_scm);
      std::string alt_conf = scm_to_locale_string(alt_conf_scm);

//       std::cout << "decoding spec :" << chain_id << ": " << resno << " :" << ins_code
// 		<< ": :" << atom_name << ": :" << alt_conf << ":" << std::endl;
      atom_spec = coot::atom_spec_t(chain_id, resno, ins_code, atom_name, alt_conf);
      // if all OK:
      atom_spec.string_user_data = "OK";
   }
   
   return atom_spec;
}
#endif // USE_GUILE

/*  ----------------------------------------------------------------------- */
/*                  Cootaneer                                               */
/*  ----------------------------------------------------------------------- */

#ifdef USE_GUILE
int cootaneer(int imol_map, int imol_model, SCM atom_in_fragment_atom_spec) {

   int istat = 0;
   coot::atom_spec_t atom_spec =
      atom_spec_from_scm_expression(atom_in_fragment_atom_spec);
   if (atom_spec.string_user_data == "Bad Spec") {
      std::cout << "Bad SCM expression for atom spec" << std::endl;
      return -1; 
   } else { 
      istat = cootaneer_internal(imol_map, imol_model, atom_spec);
      graphics_draw();
   }
   return istat;
}

int cootaneer_internal(int imol_map, int imol_model, coot::atom_spec_t &atom_spec) {
   int istat = 0;
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map)) {

	 std::string llkdfile = PKGDATADIR;
	 llkdfile += "/cootaneer-llk-2.40.dat";

	 // was that over-ridden?
	 const char *cp = getenv("COOT_PREFIX");
	 if (cp) {
	    llkdfile = cp;
	    llkdfile += "/share/coot/cootaneer-llk-2.40.dat";
	 }

	 if (!coot::file_exists(llkdfile)) {

	    std::cout << "Ooops! Can't find cootaneer likelihood data! - failure"
		      << std::endl;

	 } else { 

	    std::string chain_id = atom_spec.chain;

	    CMMDBManager *mol = graphics_info_t::molecules[imol_model].atom_sel.mol;
	    std::pair<CMMDBManager *, std::vector<coot::residue_spec_t> > mmdb_info =
	       coot::util::get_fragment_from_atom_spec(atom_spec, mol);

	    if (!mmdb_info.first) {
	       std::cout << "Bad - no fragment from atom spec" << std::endl;
	    } else {

	       // mmdb needs a spacegroup, it seems, or crash in sequence_chain().
	       // OK.  Let's bolt one in.
	       //
	       realtype a[6];
	       realtype vol;
	       int orthcode;
	       mol->GetCell(a[0], a[1], a[2], a[3], a[4], a[5], vol, orthcode);
	       char *sg = mol->GetSpaceGroup();
	       size_t l = strlen(sg+1);
	       char *sgc = new char[l];
	       strcpy(sgc, sg);
	       mmdb_info.first->SetCell(a[0], a[1], a[2], a[3], a[4], a[5]);
	       mmdb_info.first->SetSpaceGroup(sgc);

      
	       // create sequencer
	       Coot_sequence sequencer( llkdfile );
	       std::vector<std::pair<std::string, std::string> > seq =
		  graphics_info_t::molecules[imol_model].sequence_info();

	       if (seq.size() > 0) { 

		  // and apply
		  sequencer.sequence_chain(graphics_info_t::molecules[imol_map].xmap_list[0],
					   seq, *mmdb_info.first, chain_id);
		  std::string bestseq = sequencer.best_sequence();
		  std::string fullseq = sequencer.full_sequence();
		  double conf = sequencer.confidence();
		  int chnnum = sequencer.chain_number();
		  int chnoff = sequencer.chain_offset();
   
		  // write results
		  std::cout << "\nSequence: " << bestseq << "\nConfidence: " << conf << "\n";
		  if ( chnnum >= 0 ) { 
		     std::cout << "\nFrom    : " << fullseq << "\nChain id: "
			       << chnnum << "\tOffset: " << chnoff+1 << "\n";

		     if (conf > 0.9) {
			std::string chain_id = seq[chnnum].first;
			std::vector<coot::residue_spec_t> mmdb_residues = mmdb_info.second;
			int istat = 
			   graphics_info_t::molecules[imol_model].apply_sequence(imol_map,
										 mmdb_info.first, mmdb_residues, bestseq,
										 chain_id, chnoff+1);
		     }
		  }
	       } else {
		  std::string s = "Oops - no sequence information has been given to molecule\n";
		  s += "number ";
		  s += coot::util::int_to_string(imol_model);
		  info_dialog(s.c_str());
	       }
	       delete mmdb_info.first;
	    }
	 }
	 
      } else {
	 std::cout << "Not a valid map molecule " << imol_model << std::endl;
      }
   } else {
      std::cout << "Not a valid model molecule " << imol_model << std::endl;
   }

   return istat;
}
#endif // USE_GUILE

