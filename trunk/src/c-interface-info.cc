/* src/c-interface.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006 The University of York
 * Author: Paul Emsley
 * Copyright 2007 by Paul Emsley
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

#include <stdlib.h>
#include <iostream>
#include <fstream>

#ifdef USE_GUILE
#include <guile/gh.h>
#endif // USE_GUILE

#ifdef USE_PYTHON
#include "Python.h"
#endif // USE_PYTHON


#define HAVE_CIF  // will become unnessary at some stage.

#include <sys/types.h> // for stating
#include <sys/stat.h>

#if !defined _MSC_VER
#include <unistd.h>
#else
#include <windows.h>
#define snprintf _snprintf
#endif
 
#include "globjects.h" //includes gtk/gtk.h

#include "callbacks.h"
#include "cc-interface.hh"
#include "interface.h" // now that we are moving callback
		       // functionality to the file, we need this
		       // header since some of the callbacks call
		       // fuctions built by glade.

#include <vector>
#include <string>

#include "mmdb_manager.h"
#include "mmdb-extras.h"
#include "mmdb.h"
#include "mmdb-crystal.h"

#include "Cartesian.h"
#include "Bond_lines.h"

// #include "xmap-interface.h"
#include "graphics-info.h"

#include "atom-utils.h" // asc_to_graphics
// #include "db-main.h" not yet

#include "BuildCas.h"

#include "trackball.h" // adding exportable rotate interface

#include "c-interface.h"
#include "coot-database.hh"


/*  ------------------------------------------------------------------------ */
/*                         Molecule Functions       :                        */
/*  ------------------------------------------------------------------------ */
/* section Molecule Functions */
// return status, less than -9999 is for failure (eg. bad imol);
float molecule_centre_internal(int imol, int iaxis) {

   float fstat = -10000;

   if (is_valid_model_molecule(imol)) {
      if (iaxis >=0 && iaxis <=2) { 
	 coot::Cartesian c =
	    centre_of_molecule(graphics_info_t::molecules[imol].atom_sel);
	 if (iaxis == 0)
	    return c.x();
	 if (iaxis == 1)
	    return c.y();
	 if (iaxis == 2)
	    return c.z();
      }
   } else {
      std::cout << "WARNING: molecule " << imol
		<< " is not a valid model molecule number " << std::endl;
   }
   return fstat;
}

/*  ----------------------------------------------------------------------- */
/*               (Eleanor's) Residue info                                   */
/*  ----------------------------------------------------------------------- */
/* Similar to above, we need only one click though. */
void do_residue_info() {

   if (graphics_info_t::residue_info_edits->size() > 0) {

      std::string s =  "You have pending (un-Applied) residue edits\n";
      s += "Deal with them first.";
      GtkWidget *w = wrapped_nothing_bad_dialog(s);
      gtk_widget_show(w);
   } else { 
      std::cout << "Click on an atom..." << std::endl;
      graphics_info_t g;
      g.in_residue_info_define = 1;
      pick_cursor_maybe();
      graphics_info_t::pick_pending_flag = 1;
   }
} 

#ifdef USE_GUILE
SCM sequence_info(int imol) {

   SCM r = SCM_EOL;

   if (is_valid_model_molecule(imol)) { 
      std::vector<std::pair<std::string, std::string> > seq =
	 graphics_info_t::molecules[imol].sequence_info();
      
      if (seq.size() > 0) {
	 // unsigned int does't work here because then the termination
	 // condition never fails.
	 for (int iv=int(seq.size()-1); iv>=0; iv--) {
	    std::cout << "iv: " << iv << " seq.size: " << seq.size() << std::endl;
	    std::cout << "debug scming" << seq[iv].first.c_str()
		      << " and " << seq[iv].second.c_str() << std::endl;
	    SCM a = scm_makfrom0str(seq[iv].first.c_str());
	    SCM b = scm_makfrom0str(seq[iv].second.c_str());
	    SCM ls = scm_cons(a, b);
	    r = scm_cons(ls, r);
	 }
      }
   }
   return r;
} 
#endif // USE_GUILE



// Called from a graphics-info-defines routine, would you believe? :)
//
// This should be a graphics_info_t function. 
//
// The reader is graphics_info_t::apply_residue_info_changes(GtkWidget *dialog);
// 
void output_residue_info(int atom_index, int imol) {

   if (graphics_info_t::residue_info_edits->size() > 0) {

      std::string s =  "You have pending (un-Applied) residue edits.\n";
      s += "Deal with them first.";
      GtkWidget *w = wrapped_nothing_bad_dialog(s);
      gtk_widget_show(w);

   } else { 

      if (imol <graphics_info_t::n_molecules) {
	 if (graphics_info_t::molecules[imol].has_model()) {
	    if (atom_index < graphics_info_t::molecules[imol].atom_sel.n_selected_atoms) { 

	       graphics_info_t g;
	       output_residue_info_as_text(atom_index, imol);
	       PCAtom picked_atom = g.molecules[imol].atom_sel.atom_selection[atom_index];
   
	       PPCAtom atoms;
	       int n_atoms;

	       picked_atom->residue->GetAtomTable(atoms,n_atoms);
	       GtkWidget *widget = wrapped_create_residue_info_dialog();

	       CResidue *residue = picked_atom->residue; 
	       coot::residue_spec_t *res_spec_p =
		  new coot::residue_spec_t(residue->GetChainID(),
					   residue->GetSeqNum(),
					   residue->GetInsCode());

	       // fill the master atom
	       GtkWidget *master_occ_entry =
		  lookup_widget(widget, "residue_info_master_atom_occ_entry"); 
	       GtkWidget *master_b_factor_entry =
		  lookup_widget(widget, "residue_info_master_atom_b_factor_entry");

	       gtk_signal_connect (GTK_OBJECT (master_occ_entry), "changed",
				   GTK_SIGNAL_FUNC (graphics_info_t::on_residue_info_master_atom_occ_changed),
				   NULL);
	       gtk_signal_connect (GTK_OBJECT (master_b_factor_entry),
				   "changed",
				   GTK_SIGNAL_FUNC (graphics_info_t::on_residue_info_master_atom_b_factor_changed),
				   NULL);


	       gtk_entry_set_text(GTK_ENTRY(master_occ_entry), "1.00");
	       gtk_entry_set_text(GTK_ENTRY(master_b_factor_entry),
				  graphics_info_t::float_to_string(graphics_info_t::default_new_atoms_b_factor).c_str());
					   
	       gtk_object_set_user_data(GTK_OBJECT(widget), res_spec_p);
	       g.fill_output_residue_info_widget(widget, imol, atoms, n_atoms);
	       gtk_widget_show(widget);
	       g.reset_residue_info_edits();
	    }
	 }
      }
   }

   std::string cmd = "output-residue-info";
   std::vector<coot::command_arg_t> args;
   args.push_back(atom_index);
   args.push_back(imol);
   add_to_history_typed(cmd, args);
}

// This is a different action to wrapped_create_goto_atom_window,
// because we don't want to raise an already existing dialog, we need
// to destroy it.
// 
GtkWidget *wrapped_create_residue_info_dialog() { 

   GtkWidget *widget = graphics_info_t::residue_info_dialog;
   if (widget) {
      // raise/uniconify (or whatever) what we have:
      // 
// not this widget...
//       if (!GTK_WIDGET_MAPPED(widget))
// 	 gtk_widget_show(widget);
//       else
// 	 gdk_window_raise(widget->window);
      
      gtk_widget_destroy(widget);
      widget = create_residue_info_dialog();
      graphics_info_t::residue_info_dialog = widget;
   } else {

      // create (then store) a new one.
      
      widget = create_residue_info_dialog();
      graphics_info_t::residue_info_dialog = widget;
   }
   return widget; 

} 


// 23 Oct 2003: Why is this so difficult?  Because we want to attach
// atom info (what springs to mind is a pointer to the atom) for each
// entry, so that when the text in the entry is changed, we know to
// modify the atom.
//
// The problem with that is that behind our backs, that atom could
// disappear (e.g close molecule or delete residue, mutate or
// whatever), we are left with a valid looking (i.e. non-NULL)
// pointer, but the memory to which is points is invalid -> crash when
// we try to reference it.
//
// How shall we get round this?  refcounting?
//
// Instead, let's make a trivial class that contains the information
// we need to do a SelectAtoms to find the pointer to the atom, that
// class shall be called select_atom_info, it shall contain the
// molecule number, the chain id, the residue number, the insertion
// code, the atom name, the atom altconf.
// 


void output_residue_info_as_text(int atom_index, int imol) { 
   
   // It would be cool to flash the residue here.
   // (heh - it is).
   // 
   graphics_info_t g;
   PCAtom picked_atom = g.molecules[imol].atom_sel.atom_selection[atom_index];
   
   g.flash_selection(imol, 
		     picked_atom->residue->seqNum,
		     picked_atom->residue->seqNum,
		     picked_atom->altLoc,
		     picked_atom->residue->GetChainID());

   PPCAtom atoms;
   int n_atoms;

   picked_atom->residue->GetAtomTable(atoms,n_atoms);
   for (int i=0; i<n_atoms; i++) { 
      cout << "(" << imol << ") " 
	   << atoms[i]->name << "/"
	   << atoms[i]->GetModelNum()
	   << "/"
	   << atoms[i]->GetChainID()  << "/"
	   << atoms[i]->GetSeqNum()   << "/"
	   << atoms[i]->GetResName()
	   << ", occ: " 
	   << atoms[i]->occupancy 
	   << " with B-factor: "
	   << atoms[i]->tempFactor
	   << " element: \""
	   << atoms[i]->element
	   << "\""
	   << " at " << "("
	   << atoms[i]->x << "," << atoms[i]->y << "," 
	   << atoms[i]->z << ")" << endl;
   }
   std::string cmd = "output-residue-info-as-text";
   std::vector<coot::command_arg_t> args;
   args.push_back(atom_index);
   args.push_back(imol);
   add_to_history_typed(cmd, args);
}

// Actually I want to return a scheme object with occ, pos, b-factor info
// 
void
output_atom_info_as_text(int imol, const char *chain_id, int resno,
			 const char *inscode, const char *atname,
			 const char *altconf) {

   if (is_valid_model_molecule(imol)) {
      int index =
	 graphics_info_t::molecules[imol].full_atom_spec_to_atom_index(std::string(chain_id),
							resno,
							std::string(inscode),
							std::string(atname),
							std::string(altconf));
      
      CAtom *atom = graphics_info_t::molecules[imol].atom_sel.atom_selection[index];
      std::cout << "(" << imol << ") " 
		<< atom->name << "/"
		<< atom->GetModelNum()
		<< "/"
		<< atom->GetChainID()  << "/"
		<< atom->GetSeqNum()   << "/"
		<< atom->GetResName()
		<< ", occ: " 
		<< atom->occupancy 
		<< " with B-factor: "
		<< atom->tempFactor
		<< " element: \""
		<< atom->element
		<< "\""
		<< " at " << "("
		<< atom->x << "," << atom->y << "," 
		<< atom->z << ")" << std::endl;
   }
   std::string cmd = "output-atom-info-as-text";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(chain_id);
   args.push_back(resno);
   args.push_back(inscode);
   args.push_back(atname);
   args.push_back(altconf);
   add_to_history_typed(cmd, args);

}

#ifdef USE_GUILE
// (list occ temp-factor element x y z) or #f 
const char *atom_info_string(int imol, const char *chain_id, int resno,
		     const char *ins_code, const char *atname,
		     const char *altconf) {

   const char *r = 0;  // guile/SWIG sees this as #f
   if (is_valid_model_molecule(imol)) {
      int index =
	 graphics_info_t::molecules[imol].full_atom_spec_to_atom_index(std::string(chain_id),
								       resno,
								       std::string(ins_code),
								       std::string(atname),
								       std::string(altconf));
      if (index > -1) { 
	 CAtom *atom = graphics_info_t::molecules[imol].atom_sel.atom_selection[index];

	 // we need the ' because eval needs it.
	 std::string s = "(quote (";
	 s += coot::util::float_to_string(atom->occupancy);
	 s += " ";
	 s += coot::util::float_to_string(atom->tempFactor);
	 s += " ";
	 s += single_quote(atom->element);
	 s += " ";
	 s += coot::util::float_to_string(atom->x);
	 s += " ";
	 s += coot::util::float_to_string(atom->y);
	 s += " ";
	 s += coot::util::float_to_string(atom->z);
	 s += "))";
	 r = new char[s.length() + 1];
	 strcpy((char *)r, s.c_str());
      }
   }
   std::string cmd = "atom-info-string";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(chain_id);
   args.push_back(resno);
   args.push_back(ins_code);
   args.push_back(atname);
   args.push_back(altconf);
   add_to_history_typed(cmd, args);

   return r;
}
#endif // USE_GUILE

#ifdef USE_GUILE
// output is like this:
// (list
//    (list (list atom-name alt-conf)
//          (list occ temp-fact element)
//          (list x y z)))
// 
SCM residue_info(int imol, const char* chain_id, int resno, const char *ins_code) {

   SCM r = SCM_BOOL(0);
   if (is_valid_model_molecule(imol)) {
      CMMDBManager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      int imod = 1;
      
      CModel *model_p = mol->GetModel(imod);
      CChain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 std::string chain_id_mol(chain_p->GetChainID());
	 if (chain_id_mol == std::string(chain_id)) { 
	    int nres = chain_p->GetNumberOfResidues();
	    PCResidue residue_p;
	    CAtom *at;

	    // why use this bizarre contrivance to get a null list for
	    // starting? I must be missing something.
	    SCM all_atoms = SCM_CAR(scm_listofnull);
	    for (int ires=0; ires<nres; ires++) { 
	       residue_p = chain_p->GetResidue(ires);
	       std::string res_ins_code(residue_p->GetInsCode());
	       if (residue_p->GetSeqNum() == resno) { 
		  if (res_ins_code == std::string(ins_code)) {
		     int n_atoms = residue_p->GetNumberOfAtoms();
		     SCM at_info = SCM_BOOL(0);
		     SCM at_pos;
		     SCM at_occ, at_b, at_ele, at_name, at_altconf;
		     SCM at_x, at_y, at_z;
		     for (int iat=0; iat<n_atoms; iat++) {
			at = residue_p->GetAtom(iat);
			at_x  = scm_float2num(at->x);
			at_y  = scm_float2num(at->y);
			at_z  = scm_float2num(at->z);
			at_pos = scm_list_3(at_x, at_y, at_z);
			at_occ = scm_float2num(at->occupancy);
			at_b   = scm_float2num(at->tempFactor);
			at_ele = scm_makfrom0str(at->element);
			at_name = scm_makfrom0str(at->name);
			at_altconf = scm_makfrom0str(at->altLoc);
			SCM compound_name = scm_list_2(at_name, at_altconf);
			SCM compound_attrib = scm_list_3(at_occ, at_b, at_ele);
			at_info = scm_list_3(compound_name, compound_attrib, at_pos);
			all_atoms = scm_cons(at_info, all_atoms);
		     }
		  }
	       }
	    }
	    r = all_atoms;
	 }
      }
   }
   return r;
}
#endif // USE_GUILE


#ifdef USE_GUILE
SCM residue_name(int imol, const char* chain_id, int resno, const char *ins_code) {

   SCM r = SCM_BOOL(0);
   if (is_valid_model_molecule(imol)) {
      CMMDBManager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      int imod = 1;
      bool have_resname_flag = 0;
      
      CModel *model_p = mol->GetModel(imod);
      CChain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
	 chain_p = model_p->GetChain(ichain);
	 std::string chain_id_mol(chain_p->GetChainID());
	 if (chain_id_mol == std::string(chain_id)) { 
	    int nres = chain_p->GetNumberOfResidues();
	    PCResidue residue_p;
	    for (int ires=0; ires<nres; ires++) { 
	       residue_p = chain_p->GetResidue(ires);
	       if (residue_p->GetSeqNum() == resno) { 
		  std::string ins = residue_p->GetInsCode();
		  if (ins == ins_code) {
		     r = scm_makfrom0str(residue_p->GetResName());
		     have_resname_flag = 1;
		     break;
		  }
	       }
	    }
	 }
	 if (have_resname_flag)
	    break;
      }
   }
   return r;
}
#endif // USE_GUILE


// A C++ function interface:
// 
int set_go_to_atom_from_spec(const coot::atom_spec_t &atom_spec) {

   graphics_info_t g;

   g.set_go_to_atom_chain_residue_atom_name(atom_spec.chain.c_str(), 
					    atom_spec.resno,
					    atom_spec.insertion_code.c_str(), 
					    atom_spec.atom_name.c_str(),
					    atom_spec.alt_conf.c_str());

   int success = g.try_centre_from_new_go_to_atom(); 
   if (success)
      update_things_on_move_and_redraw(); 

   return success; 
}

std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom_spec() {

   coot::atom_spec_t spec;
   bool was_found_flag = 0;
   
   graphics_info_t g;
   float dist_best = 999999999.9;
   int imol_closest = -1;
   CAtom *at_close = 0;
   
   for (int imol=0; imol<graphics_info_t::n_molecules; imol++) {

      if (is_valid_model_molecule(imol)) {
	 if (graphics_info_t::molecules[imol].is_displayed_p()) { 
	    coot::at_dist_info_t at_info =
	       graphics_info_t::molecules[imol].closest_atom(g.RotationCentre());
	    if (at_info.atom) {
	       if (at_info.dist <= dist_best) {
		  dist_best = at_info.dist;
		  imol_closest = at_info.imol;
		  at_close = at_info.atom;
	       }
	    }
	 }
      }
   }
   if (at_close) {
      spec = coot::atom_spec_t(at_close);
      was_found_flag = 1;
   }

   std::pair<int, coot::atom_spec_t> p1(imol_closest, spec);
   return std::pair<bool, std::pair<int, coot::atom_spec_t> > (was_found_flag, p1);
} 

#ifdef USE_GUILE
//* \brief 
// Return a list of (list imol chain-id resno ins-code atom-name
// alt-conf) for atom that is closest to the screen centre.  If there
// are multiple models with the same coordinates at the screen centre,
// return the attributes of the atom in the highest number molecule
// number.
// 
SCM active_residue() {

   SCM s = SCM_BOOL(0);
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();

   if (pp.first) {
      s = SCM_EOL;
      s = scm_cons(scm_makfrom0str(pp.second.second.alt_conf.c_str()) , s);
      s = scm_cons(scm_makfrom0str(pp.second.second.atom_name.c_str()) , s);
      s = scm_cons(scm_makfrom0str(pp.second.second.insertion_code.c_str()) , s);
      s = scm_cons(scm_int2num(pp.second.second.resno) , s);
      s = scm_cons(scm_makfrom0str(pp.second.second.chain.c_str()) , s);
      s = scm_cons(scm_int2num(pp.second.first) ,s);
   } 
   return s;
}
#endif // USE_GUILE

#ifdef USE_GUILE
SCM closest_atom(int imol) {

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      coot::at_dist_info_t at_info =
	 graphics_info_t::molecules[imol].closest_atom(g.RotationCentre());
      if (at_info.atom) {
	 r = SCM_EOL;
	 r = scm_cons(scm_double2num(at_info.atom->z), r);
	 r = scm_cons(scm_double2num(at_info.atom->y), r);
	 r = scm_cons(scm_double2num(at_info.atom->x), r);
	 r = scm_cons(scm_makfrom0str(at_info.atom->altLoc), r);
	 r = scm_cons(scm_makfrom0str(at_info.atom->name), r);
	 r = scm_cons(scm_makfrom0str(at_info.atom->GetInsCode()), r);
	 r = scm_cons(scm_int2num(at_info.atom->GetSeqNum()), r);
	 r = scm_cons(scm_makfrom0str(at_info.atom->GetChainID()), r);
	 r = scm_cons(scm_int2num(imol), r);
      }
   }
   return r;
} 
#endif 


/*! \brief update the Go To Atom widget entries to atom closest to
  screen centre. */
void update_go_to_atom_from_current_position() {
   
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      set_go_to_atom_molecule(pp.second.first);
      set_go_to_atom_chain_residue_atom_name(pp.second.second.chain.c_str(),
					     pp.second.second.resno,
					     pp.second.second.atom_name.c_str());
      update_go_to_atom_window_on_other_molecule_chosen(pp.second.first);
   }
}


#ifdef USE_GUILE
SCM generic_string_vector_to_list_internal(const std::vector<std::string> &v) {

   SCM r = SCM_CAR(scm_listofnull);
   for (int i=v.size()-1; i>=0; i--) {
      r = scm_cons(scm_makfrom0str(v[i].c_str()), r);
   }
   return r; 
}
#endif // USE_GUILE

// and the reverse function:
#ifdef USE_GUILE
std::vector<std::string>
generic_list_to_string_vector_internal(SCM l) {
   std::vector<std::string> r;
   SCM l_length_scm = scm_length(l);

#if (SCM_MAJOR_VERSION > 1) || (SCM_MINOR_VERSION > 7)

   int l_length = scm_to_int(l_length_scm);
   for (int i=0; i<l_length; i++) {
      SCM le = scm_list_ref(l, SCM_MAKINUM(i));
      std::string s = scm_to_locale_string(le);
      r.push_back(s);
   } 
   
#else
   
   int l_length = gh_scm2int(l_length_scm);
   for (int i=0; i<l_length; i++) {
      SCM le = scm_list_ref(l, SCM_MAKINUM(i));
      std::string s =  SCM_STRING_CHARS(le);
      r.push_back(s);
   }

#endif

   return r;
}
#endif

#ifdef USE_GUILE
SCM generic_int_vector_to_list_internal(const std::vector<int> &v) {

   SCM r = SCM_EOL;
   for (int i=v.size()-1; i>=0; i--) {
      r = scm_cons(scm_int2num(v[i]), r);
   }
   return r; 
}
#endif // USE_GUILE




#ifdef USE_GUILE
SCM rtop_to_scm(const clipper::RTop_orth &rtop) {

   SCM r = SCM_EOL;

   SCM tr_list = SCM_EOL;
   SCM rot_list = SCM_EOL;

   clipper::Mat33<double>  mat = rtop.rot();
   clipper::Vec3<double> trans = rtop.trn();

   double x;
   tr_list = scm_cons(scm_double2num(trans[2]), tr_list);
   tr_list = scm_cons(scm_double2num(trans[1]), tr_list);
   tr_list = scm_cons(scm_double2num(trans[0]), tr_list);

   rot_list = scm_cons(scm_double2num(mat(2,2)), rot_list);
   rot_list = scm_cons(scm_double2num(mat(2,1)), rot_list);
   rot_list = scm_cons(scm_double2num(mat(2,0)), rot_list);
   rot_list = scm_cons(scm_double2num(mat(1,2)), rot_list);
   rot_list = scm_cons(scm_double2num(mat(1,1)), rot_list);
   rot_list = scm_cons(scm_double2num(mat(1,0)), rot_list);
   rot_list = scm_cons(scm_double2num(mat(0,2)), rot_list);
   rot_list = scm_cons(scm_double2num(mat(0,1)), rot_list);
   rot_list = scm_cons(scm_double2num(mat(0,0)), rot_list);

   r = scm_cons(tr_list, r);
   r = scm_cons(rot_list, r);
   return r;

}
#endif // USE_GUILE


#ifdef USE_GUILE
// get the symmetry operators strings for the given molecule
//
/*! \brief Return as a list of strings the symmetry operators of the
  given molecule. If imol is a not a valid molecule, return an empty
  list.*/
// 
SCM get_symmetry(int imol) {

   SCM r = SCM_CAR(scm_listofnull);
   if (is_valid_model_molecule(imol) ||
       is_valid_map_molecule(imol)) {
      std::vector<std::string> symop_list =
	 graphics_info_t::molecules[imol].get_symop_strings();
      r = generic_string_vector_to_list_internal(symop_list);
   }
   return r; 
} 
#endif 


void residue_info_apply_all_checkbutton_toggled() {

} 


void apply_residue_info_changes(GtkWidget *widget) {
   graphics_info_t g;
   g.apply_residue_info_changes(widget);
   graphics_draw();
} 

void do_distance_define() {

   std::cout << "Click on 2 atoms: " << std::endl;
   graphics_info_t g;
   g.pick_cursor_maybe();
   g.in_distance_define = 1;
   g.pick_pending_flag = 1;

} 

void do_angle_define() {

   std::cout << "Click on 3 atoms: " << std::endl;
   graphics_info_t g;
   g.pick_cursor_maybe();
   g.in_angle_define = 1;
   g.pick_pending_flag = 1;

} 

void do_torsion_define() {

   std::cout << "Click on 4 atoms: " << std::endl;
   graphics_info_t g;
   g.pick_cursor_maybe();
   g.in_torsion_define = 1;
   g.pick_pending_flag = 1;

} 

void clear_simple_distances() {
   graphics_info_t g;
   g.clear_simple_distances();
   g.normal_cursor();
   std::string cmd = "clear-simple-distances";
   std::vector<coot::command_arg_t> args;
   add_to_history_typed(cmd, args);
}

void clear_last_simple_distance() { 

   graphics_info_t g;
   g.clear_last_simple_distance();
   g.normal_cursor();
   std::string cmd = "clear-last-simple-distance";
   std::vector<coot::command_arg_t> args;
   add_to_history_typed(cmd, args);
} 

void store_geometry_dialog(GtkWidget *w) { 

   graphics_info_t g;
   g.geometry_dialog = w;
   if (w) { 
      gtk_window_set_transient_for(GTK_WINDOW(w),
				   GTK_WINDOW(lookup_widget(g.glarea, "window1")));
   }
}


void clear_residue_info_edit_list() {

   graphics_info_t::residue_info_edits->resize(0);
   std::string cmd = "clear-residue-info-edit-list";
   std::vector<coot::command_arg_t> args;
   add_to_history_typed(cmd, args);
}


/*  ----------------------------------------------------------------------- */
/*                  residue environment                                      */
/*  ----------------------------------------------------------------------- */
void fill_environment_widget(GtkWidget *widget) {

   GtkWidget *entry;
   char *text = (char *) malloc(100);
   graphics_info_t g;

   entry = lookup_widget(widget, "environment_distance_min_entry");
   snprintf(text, 99, "%-5.1f", g.environment_min_distance);
   gtk_entry_set_text(GTK_ENTRY(entry), text);
   
   entry = lookup_widget(widget, "environment_distance_max_entry");
   snprintf(text, 99, "%-5.1f" ,g.environment_max_distance);
   gtk_entry_set_text(GTK_ENTRY(entry), text);
   free(text);

   GtkWidget *toggle_button;
   toggle_button = lookup_widget(widget, "environment_distance_checkbutton");
   
   if (g.environment_show_distances == 1) {
      // we have to (temporarily) set the flag to 0 because the
      // set_active creates an event which causes
      // toggle_environment_show_distances to run (and thus turn off
      // distances if they were allowed to remain here at 1 (on).
      // Strange but true.
      g.environment_show_distances = 0;
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(toggle_button), 1);
      // std::cout << "filling: button is active" << std::endl;
   } else {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(toggle_button), 0);
      // std::cout << "filling: button is inactive" << std::endl;
   } 
}

// Called when the OK button of the environment distances dialog is clicked 
// (just before it is destroyed).
// 
void execute_environment_settings(GtkWidget *widget) {
   
   GtkWidget *entry;
   float val;
   graphics_info_t g;

   entry = lookup_widget(widget, "environment_distance_min_entry");
   const gchar *text = gtk_entry_get_text(GTK_ENTRY(entry));
   val = atof(text);
   if (val < 0 || val > 1000) {
      g.environment_min_distance = 2.2;
      std::cout <<  "nonsense value for limit using "
		<< g.environment_min_distance << std::endl;
   } else {
      g.environment_min_distance = val;
   } 

   entry = lookup_widget(widget, "environment_distance_max_entry");
   text = gtk_entry_get_text(GTK_ENTRY(entry));
   val = atof(text);
   if (val < 0 || val > 1000) {
      g.environment_max_distance = 3.2;
      std::cout <<  "nonsense value for limit using "
		<< g.environment_max_distance << std::endl;
   } else {
      g.environment_max_distance = val;
   }

   if (g.environment_max_distance < g.environment_min_distance) {
      float tmp = g.environment_max_distance;
      g.environment_max_distance = g.environment_min_distance;
      g.environment_min_distance = tmp;
   }

   GtkWidget *label_check_button;
   label_check_button = lookup_widget(widget, "environment_distance_label_atom_checkbutton");
   if (GTK_TOGGLE_BUTTON(label_check_button)->active) {
      g.environment_distance_label_atom = 1;
   }

   // not sure that this is necessary now that the toggle function is
   // active
   std::pair<int, int> r =  g.get_closest_atom();
   if (r.first >= 0) { 
      g.update_environment_distances_maybe(r.first, r.second);
   }
   graphics_draw();
}

void set_show_environment_distances(int state) {

   graphics_info_t g;
   g.environment_show_distances = state;
   if (state) {
      std::pair<int, int> r =  g.get_closest_atom();
      if (r.first >= 0) { 
	 g.update_environment_distances_maybe(r.first, r.second);
      }
   }
   graphics_draw();
}

int show_environment_distances_state() {
   return graphics_info_t::environment_show_distances;
} 


void toggle_environment_show_distances(GtkToggleButton *button) {

   graphics_info_t g;

   //    if (g.environment_show_distances == 0) {
   
   GtkWidget *hbox = lookup_widget(GTK_WIDGET(button),
				   "environment_distance_distances_frame");
   GtkWidget *label_atom_check_button =
      lookup_widget(GTK_WIDGET(button),
		    "environment_distance_label_atom_checkbutton");
   
   if (button->active) { 
      // std::cout << "toggled evironment distances on" << std::endl;
      g.environment_show_distances = 1;
      gtk_widget_set_sensitive(hbox, TRUE);
      gtk_widget_set_sensitive(label_atom_check_button, TRUE);

      // 
      std::pair<int, int> r =  g.get_closest_atom();
      std::cout << "DEBUG:: got close info: " 
		<< r.first << " " << r.second << std::endl;
      if (r.first >= 0) { 
	 g.update_environment_distances_maybe(r.first, r.second);
	 graphics_draw();
      }
      
   } else {
      // std::cout << "toggled evironment distances off" << std::endl;
      g.environment_show_distances = 0;
      gtk_widget_set_sensitive(hbox, FALSE);
      gtk_widget_set_sensitive(label_atom_check_button, FALSE);
   }
}

/* a graphics_info_t function wrapper: */
void residue_info_release_memory(GtkWidget *widget) { 

   graphics_info_t g;
   g.residue_info_release_memory(widget);

}

void unset_residue_info_widget() { 
   graphics_info_t g;
   g.residue_info_dialog = NULL; 
} 


/*  ----------------------------------------------------------------------- */
/*                  pointer distances                                      */
/*  ----------------------------------------------------------------------- */
void fill_pointer_distances_widget(GtkWidget *widget) {

   GtkWidget *min_entry   = lookup_widget(widget, "pointer_distances_min_dist_entry");
   GtkWidget *max_entry   = lookup_widget(widget, "pointer_distances_max_dist_entry");
   GtkWidget *checkbutton = lookup_widget(widget, "pointer_distances_checkbutton");
   GtkWidget *frame       = lookup_widget(widget, "pointer_distances_frame");

   float min_dist = graphics_info_t::pointer_min_dist;
   float max_dist = graphics_info_t::pointer_max_dist;
   
   gtk_entry_set_text(GTK_ENTRY(min_entry),
		      graphics_info_t::float_to_string(min_dist).c_str());
   gtk_entry_set_text(GTK_ENTRY(max_entry),
		      graphics_info_t::float_to_string(max_dist).c_str());

   if (graphics_info_t::show_pointer_distances_flag) {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), TRUE);
      gtk_widget_set_sensitive(frame, TRUE);
   } else {
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), FALSE);
      gtk_widget_set_sensitive(frame, FALSE);
   }

}

void execute_pointer_distances_settings(GtkWidget *widget) {

   GtkWidget *min_entry   = lookup_widget(widget, "pointer_distances_min_dist_entry");
   GtkWidget *max_entry   = lookup_widget(widget, "pointer_distances_max_dist_entry");
   // GtkWidget *checkbutton = lookup_widget(widget, "pointer_distances_checkbutton");

   float min_dist = 0.0;
   float max_dist = 0.0;

   float t;

   const gchar *tt = gtk_entry_get_text(GTK_ENTRY(min_entry));
   t = atof(tt);

   if (t >= 0.0 & t < 999.9)
      min_dist = t;

   tt = gtk_entry_get_text(GTK_ENTRY(max_entry));
   t = atof(tt);

   if (t >= 0.0 & t < 999.9)
      max_dist = t;

   graphics_info_t::pointer_max_dist = max_dist;
   graphics_info_t::pointer_min_dist = min_dist;

   
}


void toggle_pointer_distances_show_distances(GtkToggleButton *togglebutton) {

   GtkWidget *frame = lookup_widget(GTK_WIDGET(togglebutton),
				    "pointer_distances_frame");
   if (togglebutton->active) {
      set_show_pointer_distances(1);
      gtk_widget_set_sensitive(frame, TRUE);
   } else {
      set_show_pointer_distances(0);
      gtk_widget_set_sensitive(frame, FALSE);
   }
   
}


void set_show_pointer_distances(int istate) {

   // Use the graphics_info_t's pointer min and max dist

   std::cout << "in set_show_pointer_distances: state: "
	     << istate << std::endl;
   
   if (istate == 0) { 
      graphics_info_t::show_pointer_distances_flag = 0;
      graphics_info_t g;
      g.clear_pointer_distances();
   } else {
      graphics_info_t::show_pointer_distances_flag = 1;
      graphics_info_t g;
      // std::cout << "in set_show_pointer_distances: making distances.." << std::endl;
      g.make_pointer_distance_objects();
      // std::cout << "in set_show_pointer_distances: done making distances.." << std::endl;
   }
   graphics_draw();
   graphics_info_t::residue_info_edits->resize(0);
   std::string cmd = "set-show-pointer-distances";
   std::vector<coot::command_arg_t> args;
   args.push_back(istate);
   add_to_history_typed(cmd, args);
}


void fill_single_map_properties_dialog(GtkWidget *window, int imol) { 

   GtkWidget *cell_text = lookup_widget(window, "single_map_properties_cell_text");
   GtkWidget *spgr_text = lookup_widget(window, "single_map_properties_sg_text");

   std::string cell_text_string;
   std::string spgr_text_string;

   cell_text_string = graphics_info_t::molecules[imol].cell_text_with_embeded_newline();
   spgr_text_string = "   ";
   spgr_text_string += graphics_info_t::molecules[imol].xmap_list[0].spacegroup().descr().symbol_hm();
   spgr_text_string += "  [";
   spgr_text_string += graphics_info_t::molecules[imol].xmap_list[0].spacegroup().descr().symbol_hall();
   spgr_text_string += "]";


   gtk_label_set_text(GTK_LABEL(cell_text), cell_text_string.c_str());
   gtk_label_set_text(GTK_LABEL(spgr_text), spgr_text_string.c_str());
} 

/*  ----------------------------------------------------------------------- */
/*                  miguels orientation axes matrix                         */
/*  ----------------------------------------------------------------------- */

void
set_axis_orientation_matrix(float m11, float m12, float m13,
			    float m21, float m22, float m23,
			    float m31, float m32, float m33) {

   graphics_info_t::axes_orientation =
      GL_matrix(m11, m12, m13, m21, m22, m23, m31, m32, m33);

   std::string cmd = "set-axis-orientation-matrix";
   std::vector<coot::command_arg_t> args;
   args.push_back(m11);
   args.push_back(m12);
   args.push_back(m12);
   args.push_back(m21);
   args.push_back(m22);
   args.push_back(m23);
   args.push_back(m31);
   args.push_back(m32);
   args.push_back(m33);
   add_to_history_typed(cmd, args);
}

void
set_axis_orientation_matrix_usage(int state) {

   graphics_info_t::use_axes_orientation_matrix_flag = state;
   graphics_draw();
   std::string cmd = "set-axis-orientation-matrix-usage";
   std::vector<coot::command_arg_t> args;
   args.push_back(state);
   add_to_history_typed(cmd, args);

} 



/*  ----------------------------------------------------------------------- */
/*                  dynamic map                                             */
/*  ----------------------------------------------------------------------- */
void toggle_dynamic_map_sampling() {
   graphics_info_t g;
   // std::cout << "toggling from " << g.dynamic_map_resampling << std::endl;
   if (g.dynamic_map_resampling) {
      g.dynamic_map_resampling = 0;
   } else {
      g.dynamic_map_resampling = 1;
   }
   std::string cmd = "toggle-dynamic-map-sampling";
   std::vector<coot::command_arg_t> args;
   add_to_history_typed(cmd, args);
}


void toggle_dynamic_map_display_size() {
   graphics_info_t g;
   // std::cout << "toggling from " << g.dynamic_map_size_display << std::endl;
   if (g.dynamic_map_size_display) {
      g.dynamic_map_size_display = 0;
   } else {
      g.dynamic_map_size_display = 1;
   }
} 

void   set_map_dynamic_map_sampling_checkbutton(GtkWidget *checkbutton) {

   graphics_info_t g;
   if (g.dynamic_map_resampling) {
      g.dynamic_map_resampling = 0;
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), 1);
   }
}

void
set_map_dynamic_map_display_size_checkbutton(GtkWidget *checkbutton) {

   graphics_info_t g;
   if (g.dynamic_map_size_display) {
      g.dynamic_map_size_display = 0; 
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(checkbutton), 1);
   }
} 

void
set_dynamic_map_size_display_on() {
   graphics_info_t::dynamic_map_size_display = 1;
} 
void
set_dynamic_map_size_display_off() {
   graphics_info_t::dynamic_map_size_display = 0;
} 
void
set_dynamic_map_sampling_on() {
   graphics_info_t::dynamic_map_resampling = 1;
}
void
set_dynamic_map_sampling_off() {
   graphics_info_t::dynamic_map_resampling = 0;
}

void
set_dynamic_map_zoom_offset(int i) {
   graphics_info_t::dynamic_map_zoom_offset = i;
} 




/*  ------------------------------------------------------------------------ */
/*                         history                                           */
/*  ------------------------------------------------------------------------ */



void add_to_history(const std::vector<std::string> &command_strings) {
   // something
   graphics_info_t g;
   g.add_history_command(command_strings);

   if (g.console_display_commands) 
      std::cout << "INFO:: Command: "
		<< graphics_info_t::schemize_command_strings(command_strings)
		<< std::endl;

#ifdef USE_MYSQL_DATABASE

   add_to_database(command_strings);
#endif
   
}

void add_to_history_typed(const std::string &command,
			  const std::vector<coot::command_arg_t> &args) {

   std::vector<std::string> command_strings;

   command_strings.push_back(command);
   for (int i=0; i<args.size(); i++)
      command_strings.push_back(args[i].as_string());

   add_to_history(command_strings);
}

void add_to_history_simple(const std::string &s) {

   std::vector<std::string> command_strings;
   command_strings.push_back(s);
   add_to_history(command_strings);
} 

std::string
single_quote(const std::string &s) {
   std::string r("\"");
   r += s;
   r += "\"";
   return r;
} 

std::string pythonize_command_name(const std::string &s) {

   std::string ss;
   for (unsigned int i=0; i<s.length(); i++) {
      if (s[i] == '-') {
	 ss += '_';
      } else {
	 ss += s[i];	 
      }
   }
   return ss;
} 

std::string schemize_command_name(const std::string &s) {

   std::string ss;
   for (unsigned int i=0; i<s.length(); i++) {
      if (s[i] == '_') {
	 ss += '-';
      } else {
	 ss += s[i];	 
      }
   }
   return ss;
}

#ifdef USE_MYSQL_DATABASE
int db_query_insert(const std::string &insert_string) {

   int v = -1;
   if (graphics_info_t::mysql) {

      time_t timep = time(0);
      long int li = timep;
      
      clipper::String query("insert into session");
      query += " (userid, sessionid, command, commandnumber, timeasint)";
      query += " values ('";
      query += graphics_info_t::db_userid_username.first;
      query += "', '";
      query += graphics_info_t::sessionid;
      query += "', '";
      query += insert_string;
      query += "', ";
      query += graphics_info_t::int_to_string(graphics_info_t::query_number);
      query += ", ";
      query += coot::util::long_int_to_string(li);
      query += ") ; ";

//       query = "insert into session ";
//       query += "(userid, sessionid, command, commandnumber) values ";
//       query += "('pemsley', 'sesh', 'xxx', ";
//       query += graphics_info_t::int_to_string(graphics_info_t::query_number);
//       query += ") ;";
      
//       std::cout << "query: " << query << std::endl;
      unsigned long length = query.length();
      v = mysql_real_query(graphics_info_t::mysql, query.c_str(), length);
      if (v != 0) {
	 if (v == CR_COMMANDS_OUT_OF_SYNC)
	    std::cout << "WARNING:: MYSQL Commands executed in an"
		      << " improper order" << std::endl;
	 if (v == CR_SERVER_GONE_ERROR) 
	    std::cout << "WARNING:: MYSQL Server gone!"
		      << std::endl;
	 if (v == CR_SERVER_LOST) 
	    std::cout << "WARNING:: MYSQL Server lost during query!"
		      << std::endl;
	 if (v == CR_UNKNOWN_ERROR) 
	    std::cout << "WARNING:: MYSQL Server transaction had "
		      << "an uknown error!" << std::endl;
	 std::cout << "history: mysql_real_query returned " << v
		   << std::endl;
      }
      graphics_info_t::query_number++;
   }
   return v;
}
#endif // USE_MYSQL_DATABASE

void add_to_database(const std::vector<std::string> &command_strings) {

#ifdef USE_MYSQL_DATABASE
   std::string insert_string =
      graphics_info_t::schemize_command_strings(command_strings);
   db_query_insert(insert_string);
#endif // USE_MYSQL_DATABASE

}


#ifdef USE_MYSQL_DATABASE
// 
void db_finish_up() {

   db_query_insert(";#finish");

} 
#endif // USE_MYSQL_DATABASE


/*  ----------------------------------------------------------------------- */
/*                         History/scripting                                */
/*  ----------------------------------------------------------------------- */

void print_all_history_in_scheme() {

   graphics_info_t g;
   std::vector<std::vector<std::string> > ls = g.history_list.history_list();
   for (unsigned int i=0; i<ls.size(); i++)
      std::cout << i << "  " << graphics_info_t::schemize_command_strings(ls[i]) << "\n";

   add_to_history_simple("print-all-history-in-scheme");

}

void print_all_history_in_python() {

   graphics_info_t g;
   std::vector<std::vector<std::string> > ls = g.history_list.history_list();
   for (unsigned int i=0; i<ls.size(); i++)
      std::cout << i << "  " << graphics_info_t::pythonize_command_strings(ls[i]) << "\n";
   add_to_history_simple("print-all-history-in-python");
}

/*! \brief set a flag to show the text command equivalent of gui
  commands in the console as they happen. 

  1 for on, 0 for off. */
void set_console_display_commands_state(short int istate) {

   graphics_info_t::console_display_commands = istate;
} 


std::string languagize_command(const std::vector<std::string> &command_parts) {

   short int language = 0; 
#ifdef USE_PYTHON
#ifdef USE_GUILE
   language = coot::STATE_SCM;
#else   
   language = coot::STATE_PYTHON;
#endif
#endif

#ifdef USE_GUILE
   language = 1;
#endif

   std::string s;
   if (language) {
      if (language == coot::STATE_PYTHON)
	 s = graphics_info_t::pythonize_command_strings(command_parts);
      if (language == coot::STATE_SCM)
	 s = graphics_info_t::schemize_command_strings(command_parts);
   }
   return s; 
}


/*  ------------------------------------------------------------------------ */
/*                     state (a graphics_info thing)                         */
/*  ------------------------------------------------------------------------ */

void
save_state() {
   graphics_info_t g;
   g.save_state();
   add_to_history_simple("save-state");
}

void
save_state_file(const char *filename) {
   graphics_info_t g;
   g.save_state_file(std::string(filename));
   std::string cmd = "save-state-file";
   std::vector<coot::command_arg_t> args;
   args.push_back(filename);
   add_to_history_typed(cmd, args);
}


/*  ------------------------------------------------------------------------ */
/*                     resolution                                            */
/*  ------------------------------------------------------------------------ */

/* Return negative number on error, otherwise resolution in A (eg. 2.0) */
float data_resolution(int imol) {

   float r = -1;
   if (is_valid_map_molecule(imol)) {
      r = graphics_info_t::molecules[imol].data_resolution();
   }
   std::string cmd = "data-resolution";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   add_to_history_typed(cmd, args);
   return r;
}


/*  ------------------------------------------------------------------------ */
/*                     residue exists?                                       */
/*  ------------------------------------------------------------------------ */

int does_residue_exist_p(int imol, char *chain_id, int resno, char *inscode) {

   int istate = 0;
   if (is_valid_model_molecule(imol)) {
      istate = graphics_info_t::molecules[imol].does_residue_exist_p(std::string(chain_id),
								     resno,
								     std::string(inscode));
   }
   std::string cmd = "does-residue-exist-p";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(chain_id);
   args.push_back(resno);
   args.push_back(inscode);
   add_to_history_typed(cmd, args);
   return istate;
}

/*  ------------------------------------------------------------------------ */
/*                         Parameters from map:                              */
/*  ------------------------------------------------------------------------ */
/*! \brief return the mtz file that was use to generate the map

  return 0 when there is no mtz file associated with that map (it was
  generated from a CCP4 map file say). */
const char *mtz_hklin_for_map(int imol_map) {

   std::string mtz;

   if (is_valid_map_molecule(imol_map)) {
      mtz = graphics_info_t::molecules[imol_map].save_mtz_file_name;
   }
   std::string cmd = "mtz-hklin-for-map";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol_map);
   add_to_history_typed(cmd, args);
   return mtz.c_str();
}

/*! \brief return the FP column in the file that was use to generate
  the map

  return 0 when there is no mtz file associated with that map (it was
  generated from a CCP4 map file say). */
const char *mtz_fp_for_map(int imol_map) {

   std::string fp;
   if (is_valid_map_molecule(imol_map)) {
      std::cout << imol_map << " Is valid map molecule" << std::endl;
      fp = graphics_info_t::molecules[imol_map].save_f_col;
   }
   std::string cmd = "mtz-fp-for-map";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol_map);
   add_to_history_typed(cmd, args);
   return fp.c_str();
} 

/*! \brief return the phases column in mtz file that was use to generate
  the map

  return 0 when there is no mtz file associated with that map (it was
  generated from a CCP4 map file say). */
const char *mtz_phi_for_map(int imol_map) {

   std::string phi;
   if (is_valid_map_molecule(imol_map)) {
      phi = graphics_info_t::molecules[imol_map].save_phi_col;
   }
   std::string cmd = "mtz-phi-for-map";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol_map);
   add_to_history_typed(cmd, args);
   return phi.c_str();
}

/*! \brief return the weight column in the mtz file that was use to
  generate the map

  return 0 when there is no mtz file associated with that map (it was
  generated from a CCP4 map file say) or no weights were used. */
const char *mtz_weight_for_map(int imol_map) {

   std::string weight;
   if (is_valid_map_molecule(imol_map)) {
      weight = graphics_info_t::molecules[imol_map].save_weight_col;
   }
   std::string cmd = "mtz-weight-for-map";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol_map);
   add_to_history_typed(cmd, args);
   return weight.c_str();
}

/*! \brief return flag for whether weights were used that was use to
  generate the map

  return 0 when no weights were used or there is no mtz file
  associated with that map. */
short int mtz_use_weight_for_map(int imol_map) {

   short int i;
   if (is_valid_map_molecule(imol_map)) {
      i = graphics_info_t::molecules[imol_map].save_use_weights;
   }
   std::string cmd = "mtz-use-weight-for-map";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol_map);
   add_to_history_typed(cmd, args);
   return i;
} 


/*! \brief Put text at x,y,z  */
// use should be given access to colour and size.
int place_text(const char *text, float x, float y, float z, int size) {

   int handle = graphics_info_t::generic_texts_p->size();
   std::string s(text); 
   coot::generic_text_object_t o(s, handle, x, y, z); 
   graphics_info_t::generic_texts_p->push_back(o);
   //   return graphics_info_t::generic_text->size() -1; // the index of the
	  					    // thing we just
						    // pushed.
   std::string cmd = "place-text";
   std::vector<coot::command_arg_t> args;
   args.push_back(text);
   args.push_back(x);
   args.push_back(y);
   args.push_back(z);
   args.push_back(size);
   add_to_history_typed(cmd, args);

   return handle; // same value as above.
} 

void remove_text(int text_handle) {

   std::vector<coot::generic_text_object_t>::iterator it;
   for (it = graphics_info_t::generic_texts_p->begin();
	it != graphics_info_t::generic_texts_p->end();
	it++) {
      if (it->handle == text_handle) {
	 graphics_info_t::generic_texts_p->erase(it);
	 break;
      }
   }
   std::string cmd = "remove-text";
   std::vector<coot::command_arg_t> args;
   args.push_back(text_handle);
   add_to_history_typed(cmd, args);
   graphics_draw();
} 

/*  ----------------------------------------------------------------------- */
/*                         Dictionaries                                     */
/*  ----------------------------------------------------------------------- */
/*! \brief return a list of all the dictionaries read */

#ifdef USE_GUILE
SCM dictionaries_read() {

   return generic_string_vector_to_list_internal(*graphics_info_t::cif_dictionary_filename_vec);
}
#endif



/*  ----------------------------------------------------------------------- */
/*                  CCP4MG Interface                                        */
/*  ----------------------------------------------------------------------- */
void write_ccp4mg_picture_description(const char *filename) {

   std::ofstream mg_stream(filename);

   if (!mg_stream) {
      std::cout << "WARNING:: can't open file " << filename << std::endl;
   } else {
      mg_stream << "# CCP4mg picture definition file\n";
      mg_stream << "# See http://www.ysbl.york.ac.uk/~ccp4mg/ccp4mg_help/picture_definition.html\n";
      mg_stream << "# or your local copy of CCP4mg documentation.\n";
      mg_stream << "# created by Coot " << coot_version() << "\n";

      // View:
      graphics_info_t g;
      coot::Cartesian rc = g.RotationCentre();
      mg_stream << "view = View (\n";
      mg_stream << "    centre_xyz = [";
      mg_stream << -rc.x() << ", " << -rc.y() << ", " << -rc.z() << "],\n";
      mg_stream << "    radius = " << 0.75*graphics_info_t::zoom << ",\n";
      //       mg_stream << "    orientation = [ " << g.quat[0] << ", "
      // 		<< g.quat[1] << ", " << g.quat[2] << ", " << g.quat[3] << "]\n";
      // Stuart corrects the orientation specification:
      mg_stream << "    orientation = [ " << -g.quat[3] << ", "
		<< g.quat[0] << ", " << g.quat[1] << ", " << g.quat[2] << "]\n";
      mg_stream << ")\n";

      // Molecules:
      for (int imol=0; imol<graphics_info_t::n_molecules; imol++) {
	 if (is_valid_model_molecule(imol)) {
	    mg_stream << "MolData (\n";
	    mg_stream << "   filename = [\n";
	    mg_stream << "   'FULLPATH',\n";
	    mg_stream << "   " << single_quote(graphics_info_t::molecules[imol].name_) << ",\n";
	    mg_stream << "   " << single_quote(graphics_info_t::molecules[imol].name_) << "])\n";
	    mg_stream << "\n";
	    mg_stream << "MolDisp (\n";
	    mg_stream << "    selection = 'all',\n";
	    mg_stream << "    colour    = 'atomtype',\n"; 
	    mg_stream << "    style     = 'BONDS',\n";
	    mg_stream << "    selection_parameters = {\n";
	    mg_stream << "                'neighb_rad' : '5.0',\n";
	    mg_stream << "                'select' : 'all',\n";
	    mg_stream << "                'neighb_excl' : 0  },\n";
	    mg_stream << "    colour_parameters =  {\n";
	    mg_stream << "                'colour_mode' : 'atomtype',\n";
	    mg_stream << "                'user_scheme' : [ ]  })\n";
	    mg_stream << "\n";
	 }
	 if (is_valid_map_molecule(imol)) {
	    std::string phi = g.molecules[imol].save_phi_col;
	    std::string f   = g.molecules[imol].save_f_col;
	    std::string w   = g.molecules[imol].save_weight_col;
	    float lev       = g.molecules[imol].contour_level[0];
	    float r         = g.box_radius;
	    int use_weights_flag = g.molecules[imol].save_use_weights;
	    std::string name = single_quote(graphics_info_t::molecules[imol].save_mtz_file_name);
	    mg_stream << "MapData (\n";
	    mg_stream << "   filename = [\n";
	    mg_stream << "   'FULLPATH',\n";
	    mg_stream << "   " << name << ",\n";
	    mg_stream << "   " << name << "],\n";
	    mg_stream << "   phi = " << single_quote(phi) << ",\n";
	    mg_stream << "   f   = " << single_quote(f) << ",\n";
	    mg_stream << "   rate = 0.75\n";
	    mg_stream << ")\n";
	    mg_stream << "MapDisp (\n";
	    mg_stream << "    contour_level = " << lev << ",\n";
	    mg_stream << "    surface_style = 'chickenwire',\n";
	    mg_stream << "    radius = " << r << ",\n";
	    mg_stream << "    clip_mode = 'OFF')\n";
	    mg_stream << "\n";
	 }
      }
   }
} 
