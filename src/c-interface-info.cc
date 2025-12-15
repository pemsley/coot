
/* src/c-interface.cc
 *
 * Copyright 2002, 2003, 2004, 2005, 2006 The University of York
 * Author: Paul Emsley
 * Copyright 2007 by Paul Emsley
 * Copyright 2007, 2008, 2009 The University of Oxford
 * Copyright 2008 by The University of Oxford
 * Author: Paul Emsley
 * Copyright 2007, 2008 by Bernhard Lohkamp
 * Copyright 2007, 2008 The University of York
 * Copyright 2012, 2013, 2015, 2016 by Medical Research Council
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

#ifdef USE_PYTHON
#include <Python.h>  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#include "python-3-interface.hh"
#endif

#include "compat/coot-sysdep.h"


#include <stdlib.h>
#include <string.h> // strncmp
#include <iostream>
#include <fstream>
#include <stdexcept>


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

#include <vector>
#include <string>

#include <mmdb2/mmdb_manager.h>

#include "coords/mmdb-extras.hh"
#include "coords/mmdb.hh"
#include "coords/mmdb-crystal.hh"

#include "coot-utils/coot-map-utils.hh" // for make_rtop_orth_from()

#include "coords/Cartesian.hh"
#include "coords/Bond_lines.hh"

#include "graphics-info.h"

#include "skeleton/BuildCas.h"
#include "ligand/primitive-chi-angles.hh"

#include "c-interface.h"
#include "c-interface-gtk-widgets.h"
#include "coot-database.hh"

#include "widget-from-builder.hh"

#include "guile-fixups.h"

#include "cc-interface.hh"
#include "cc-interface-scripting.hh"
#include "c-interface-python.hh"

#include "widget-headers.hh"

/*  ----------------------------------------------------------------- */
/*                         Scripting:                                 */
/*  ----------------------------------------------------------------- */

#ifdef USE_GUILE
SCM coot_has_python_p() {

   SCM r = SCM_BOOL_F;
#ifdef USE_PYTHON
   r = SCM_BOOL_T;
#endif
   return r;
}
#endif

#ifdef USE_PYTHON
PyObject *coot_has_guile() {
   PyObject *r = Py_False;
#ifdef USE_GUILE
   r = Py_True;
#endif
   Py_INCREF(r);
   return r;
}
#endif

bool coot_can_do_lidia_p() {

   bool r = false;

#ifdef HAVE_GOOCANVAS
#if ( ( (GTK_MAJOR_VERSION == 2) && (GTK_MINOR_VERSION > 11) ) || GTK_MAJOR_VERSION > 2)
   r = true;
#endif
#endif

   return r;

}


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

int is_shelx_molecule(int imol) {

   int r=0;
   if (is_valid_model_molecule(imol)) {
      r = graphics_info_t::molecules[imol].is_from_shelx_ins();
   }
   return r;

}


/*  ----------------------------------------------------------------------- */
/*               (Eleanor's) Residue info                                   */
/*  ----------------------------------------------------------------------- */
/* Similar to above, we need only one click though. */
void do_residue_info_dialog() {

   // 20230515-PE It doesn't work like this now.
   // Delete this function FIXME

   if (graphics_info_t::residue_info_edits.size() > 0) {

      std::string s =  "WARNING:: You have pending (un-Applied) residue edits\n";
      s += "Deal with them first.";
      GtkWidget *w = wrapped_nothing_bad_dialog(s);
      gtk_widget_set_visible(w, TRUE);
   } else {
      std::cout << "INFO:: Click on an atom..." << std::endl;
      add_status_bar_text("Click on an atom");
      graphics_info_t g;
      g.in_residue_info_define = 1;
      pick_cursor_maybe();
      graphics_info_t::pick_pending_flag = 1;
   }
}

#ifdef USE_GUILE
SCM sequence_info(int imol) {

   SCM r = SCM_BOOL_F;

   if (is_valid_model_molecule(imol)) {
      std::vector<std::pair<std::string, std::string> > seq =
         graphics_info_t::molecules[imol].sequence_info();

      if (seq.size() > 0) {
         r = SCM_EOL;
         // unsigned int does't work here because then the termination
         // condition never fails.
         for (int iv=int(seq.size()-1); iv>=0; iv--) {
//             std::cout << "iv: " << iv << " seq.size: " << seq.size() << std::endl;
//             std::cout << "debug scming" << seq[iv].first.c_str()
//                       << " and " << seq[iv].second.c_str() << std::endl;
            SCM a = scm_from_locale_string(seq[iv].first.c_str());
            SCM b = scm_from_locale_string(seq[iv].second.c_str());
            SCM ls = scm_cons(a, b);
            r = scm_cons(ls, r);
         }
      }
   }
   return r;
}
#endif // USE_GUILE


#ifdef USE_PYTHON
PyObject *sequence_info_py(int imol) {

   PyObject *r = Py_False;

   if (is_valid_model_molecule(imol)) {

      std::vector<std::pair<std::string, std::string> > seq =
         graphics_info_t::molecules[imol].sequence_info();


      if (seq.size() > 0) {
         // unsigned int does't work here because then the termination
         // condition never fails.
         r = PyList_New(seq.size());
         PyObject *a;
         PyObject *b;
         PyObject *ls;
         for (int iv=int(seq.size()-1); iv>=0; iv--) {
            //std::cout << "iv: " << iv << " seq.size: " << seq.size() << std::endl;
            //std::cout << "debug pythoning " << seq[iv].first.c_str()
            //              << " and " << seq[iv].second.c_str() << std::endl;
            a = myPyString_FromString(seq[iv].first.c_str());
            b = myPyString_FromString(seq[iv].second.c_str());
            ls = PyList_New(2);
            PyList_SetItem(ls, 0, a);
            PyList_SetItem(ls, 1, b);
            PyList_SetItem(r, iv, ls);
         }
      }
   }
   if (PyBool_Check(r)) {
      Py_INCREF(r);
   }
   return r;
}
#endif // USE_PYTHON


// Called from a graphics-info-defines routine, would you believe? :)
//
// This should be a graphics_info_t function.
//
// The reader is graphics_info_t::apply_residue_info_changes(GtkWidget *dialog);
//
void output_residue_info_dialog(int imol, int atom_index) {

   graphics_info_t g;
   g.output_residue_info_dialog(imol, atom_index);
   std::string cmd = "output-residue-info";
   std::vector<coot::command_arg_t> args;
   args.push_back(atom_index);
   args.push_back(imol);
   add_to_history_typed(cmd, args);
}

void output_residue_info_dialog(int imol, const coot::residue_spec_t &rs) {

   graphics_info_t g;
   g.output_residue_info_dialog(imol, rs);
}

void
residue_info_dialog(int imol, const char *chain_id, int resno, const char *ins_code) {

   if (is_valid_model_molecule(imol)) {
      int atom_index = -1;
      mmdb::Residue *res = graphics_info_t::molecules[imol].residue_from_external(resno, ins_code, chain_id);
      mmdb::PPAtom residue_atoms;
      int n_residue_atoms;
      res->GetAtomTable(residue_atoms, n_residue_atoms);
      if (n_residue_atoms > 0) {
         mmdb::Atom *at = residue_atoms[0];
         int handle = graphics_info_t::molecules[imol].atom_sel.UDDAtomIndexHandle;
         int ierr = at->GetUDData(handle, atom_index);
         if (ierr == mmdb::UDDATA_Ok) {
            if (atom_index != -1) {
               output_residue_info_dialog(imol, atom_index);
            }
         }
      }
   }
}

// change this argument order one day. Care needed.
void
output_residue_info_as_text(int atom_index, int imol) {
   graphics_info_t g;
   g.output_residue_info_as_text(atom_index, imol);
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

      mmdb::Atom *atom = graphics_info_t::molecules[imol].atom_sel.atom_selection[index];
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
      try {
         // chi angles:
         coot::primitive_chi_angles chi_angles(atom->residue);
         std::vector<coot::alt_confed_chi_angles> chis = chi_angles.get_chi_angles();
         if (chis.size() > 0) {
            unsigned int i_chi_set = 0;
            std::cout << "   Chi Angles:" << std::endl;
            for (unsigned int ich=0; ich<chis[i_chi_set].chi_angles.size(); ich++) {
               std::cout << "     chi "<< chis[i_chi_set].chi_angles[ich].first << ": "
                         << chis[i_chi_set].chi_angles[ich].second
                         << " degrees" << std::endl;
            }
         } else {
            std::cout << "No Chi Angles for this residue" << std::endl;
         }
      }
      catch (const std::runtime_error &mess) {
         std::cout << mess.what() << std::endl;
      }
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

std::string atom_info_as_text_for_statusbar(int atom_index, int imol) {

  std::string ai;
  ai = "";
  if (is_valid_model_molecule(imol)) {
     if (atom_index >= 0) {
        graphics_info_t g;
        if (atom_index < g.molecules[imol].atom_sel.n_selected_atoms) {
           mmdb::Atom *at = g.molecules[imol].atom_sel.atom_selection[atom_index];
           std::string alt_conf_bit("");
           if (strncmp(at->altLoc, "", 1))
              alt_conf_bit=std::string(",") + std::string(at->altLoc);
           std::string mol_name = coot::util::file_name_non_directory(g.molecules[imol].get_name());
           ai += " [ ";
           ai += graphics_info_t::int_to_string(imol);
           ai += " \"";
           ai += mol_name;
           ai += "\"] ";
           ai += at->name;
           ai += alt_conf_bit;
           ai += "/";
           ai += graphics_info_t::int_to_string(at->GetModelNum());
           ai += "/";
           ai += at->GetChainID();
           ai += "/";
           ai += graphics_info_t::int_to_string(at->GetSeqNum());
           ai += at->GetInsCode();
           ai += " ";
           ai += at->GetResName();
           ai += " occupancy: ";
           ai += graphics_info_t::float_to_string(at->occupancy);
           ai += " b-factor: ";
           ai += graphics_info_t::float_to_string(at->tempFactor);
           ai += " ele: ";
           ai += at->element;
           ai += " pos: (";
           // using atom positions (ignoring symmetry etc)
           ai += graphics_info_t::float_to_string(at->x);
           ai += ",";
           ai += graphics_info_t::float_to_string(at->y);
           ai += ",";
           ai += graphics_info_t::float_to_string(at->z);
           ai += ")";
        }
     }
  }
  return ai;
}


std::string
atom_info_as_text_for_statusbar(int atom_index, int imol,
                                const std::pair<symm_trans_t, Cell_Translation> &sts) {

  std::string ai;
  ai = "";
  if (is_valid_model_molecule(imol)) {
    mmdb::Atom *at = graphics_info_t::molecules[imol].atom_sel.atom_selection[atom_index];
    std::string alt_conf_bit("");
    if (strncmp(at->altLoc, "", 1))
      alt_conf_bit=std::string(",") + std::string(at->altLoc);
    ai += "(mol. no: ";
    ai += graphics_info_t::int_to_string(imol);
    ai += ") ";
    ai += at->name;
    ai += alt_conf_bit;
    ai += "/";
    ai += graphics_info_t::int_to_string(at->GetModelNum());
    ai += "/";
    ai += at->GetChainID();
    ai += "/";
    ai += graphics_info_t::int_to_string(at->GetSeqNum());
    ai += at->GetInsCode();
    ai += " ";
    ai += at->GetResName();
    // ignoring symmetry?!, no
    ai += " ";
    ai += to_string(sts);
    ai += " occ: ";
    ai += graphics_info_t::float_to_string(at->occupancy);
    ai += " bf: ";
    ai += graphics_info_t::float_to_string(at->tempFactor);
    ai += " ele: ";
    ai += at->element;
    ai += " pos: (";
    // using atom positions (ignoring symmetry etc)
    ai += graphics_info_t::float_to_string(at->x);
    ai += ",";
    ai += graphics_info_t::float_to_string(at->y);
    ai += ",";
    ai += graphics_info_t::float_to_string(at->z);
    ai += ")";
  }

  return ai;

}


#ifdef USE_GUILE
// (list occ temp-factor element x y z) or #f
SCM atom_info_string_scm(int imol, const char *chain_id, int resno,
                         const char *ins_code, const char *atname,
                         const char *altconf) {

   SCM r = SCM_BOOL_F;

   if (is_valid_model_molecule(imol)) {
      int index =
         graphics_info_t::molecules[imol].full_atom_spec_to_atom_index(std::string(chain_id),
                                                                       resno,
                                                                       std::string(ins_code),
                                                                       std::string(atname),
                                                                       std::string(altconf));
      if (index > -1) {
         mmdb::Atom *atom = graphics_info_t::molecules[imol].atom_sel.atom_selection[index];

         r = SCM_EOL;

         r = scm_cons(scm_from_double(atom->occupancy), r);
         r = scm_cons(scm_from_double(atom->tempFactor), r);
         r = scm_cons(scm_from_locale_string(atom->element), r);
         r = scm_cons(scm_from_double(atom->x), r);
         r = scm_cons(scm_from_double(atom->y), r);
         r = scm_cons(scm_from_double(atom->z), r);
         r = scm_reverse(r);
      }
   }
   std::string cmd = "atom-info-string";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol);
   args.push_back(coot::util::single_quote(chain_id));
   args.push_back(resno);
   args.push_back(coot::util::single_quote(ins_code));
   args.push_back(coot::util::single_quote(atname));
   args.push_back(coot::util::single_quote(altconf));
   add_to_history_typed(cmd, args);

   return r;
}
#endif // USE_GUILE

#ifdef USE_GUILE
SCM molecule_to_pdb_string_scm(int imol) {
   SCM r = SCM_EOL;
   if (is_valid_model_molecule(imol)) {
      std::string s = graphics_info_t::molecules[imol].pdb_string();
      r = scm_from_locale_string(s.c_str());
   }
   return r;
}
#endif

#ifdef USE_PYTHON
PyObject *molecule_to_pdb_string_py(int imol) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      std::string s = graphics_info_t::molecules[imol].pdb_string();
      // std::cout << "s: " << s << std::endl;
      r = myPyString_FromString(s.c_str());
   }

   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;

}
#endif



// BL says:: we return a string in python list compatible format.
// to use it in python you need to eval the string!
#ifdef USE_PYTHON
// "[occ,temp_factor,element,x,y,z]" or 0
PyObject *atom_info_string_py(int imol, const char *chain_id, int resno,
                              const char *ins_code, const char *atname,
                              const char *altconf) {

   PyObject *r = Py_False;

   if (is_valid_model_molecule(imol)) {
      int index =
         graphics_info_t::molecules[imol].full_atom_spec_to_atom_index(std::string(chain_id),
                                                                       resno,
                                                                       std::string(ins_code),
                                                                       std::string(atname),
                                                                       std::string(altconf));
      if (index > -1) {
         mmdb::Atom *atom = graphics_info_t::molecules[imol].atom_sel.atom_selection[index];

         r = PyList_New(6);
         PyList_SetItem(r, 0, PyFloat_FromDouble(atom->occupancy));
         PyList_SetItem(r, 1, PyFloat_FromDouble(atom->tempFactor));
         PyList_SetItem(r, 2, myPyString_FromString(atom->element));
         PyList_SetItem(r, 3, PyFloat_FromDouble(atom->x));
         PyList_SetItem(r, 4, PyFloat_FromDouble(atom->y));
         PyList_SetItem(r, 5, PyFloat_FromDouble(atom->z));
      }
   }
   std::string cmd = "atom_info_string";
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
#endif // PYTHON



#ifdef USE_GUILE
#ifdef USE_PYTHON
PyObject *scm_to_py(SCM s) {

   PyObject *o = Py_None;
   if (scm_is_true(scm_list_p(s))) {
      SCM s_length_scm = scm_length(s);
      int s_length = scm_to_int(s_length_scm);
      o = PyList_New(s_length);
      for (int item=0; item<s_length; item++) {
         SCM item_scm = scm_list_ref(s, scm_from_int(item));
         PyList_SetItem(o, item, scm_to_py(item_scm));
      }
   } else {
      if (scm_is_true(scm_boolean_p(s))) {
         if (scm_is_true(s)) {
            o = Py_True;
         } else {
            o = Py_False;
         }
      } else {
         if (scm_is_true(scm_integer_p(s))) {
            int iscm = scm_to_int(s);
            o = PyLong_FromLong(iscm);
         } else {
            if (scm_is_true(scm_real_p(s))) {
               float f = scm_to_double(s);
               o = PyFloat_FromDouble(f);

            } else {
               if (scm_is_true(scm_string_p(s))) {
                  std::string str = scm_to_locale_string(s);
                  o = myPyString_FromString(str.c_str());
               }
            }
         }
      }
   }

   if (PyBool_Check(o) || o == Py_None) {
     Py_INCREF(o);
   }

   return o;
}
#endif // USE_GUILE
#endif // USE_PYTHON


#ifdef USE_GUILE
#ifdef USE_PYTHON
SCM py_to_scm(PyObject *o) {

   SCM s = SCM_BOOL_F;
   if (! o) return s;
   if (PyList_Check(o)) {
      int l = PyObject_Length(o);
      s = SCM_EOL;
      for (int item=0; item<l; item++) {
         PyObject *py_item = PyList_GetItem(o, item);
         if (py_item == NULL) {
           PyErr_Print();
         }
         SCM t = py_to_scm(py_item);
         s = scm_cons(t, s);
      }
      s = scm_reverse(s);
   } else {
      if (PyBool_Check(o)) {
         s = SCM_BOOL_F;
         int i = PyLong_AsLong(o);
         if (i)
            s = SCM_BOOL_T;
      } else {
         if (PyLong_Check(o)) {
            int i=PyLong_AsLong(o);
            s = scm_from_int(i);
         } else {
            if (PyFloat_Check(o)) {
               double f = PyFloat_AsDouble(o);
               s = scm_from_double(f);
            } else {
               if (PyUnicode_Check(o)) {
                  std::string str = PyBytes_AS_STRING(PyUnicode_AsUTF8String(o));
                  s = scm_from_locale_string(str.c_str());
               } else {
                  if (o == Py_None) {
                     s = SCM_UNSPECIFIED;
                  }
               }
            }
         }
      }
   }
   return s;
}

#endif // USE_GUILE
#endif // USE_PYTHON

// Return residue specs for residues that have atoms that are
// closer than radius Angstroems to any atom in the residue
// specified by res_in.
//
#ifdef USE_GUILE
SCM residues_near_residue(int imol, SCM residue_in, float radius) {

   SCM r = SCM_EOL;
   if (is_valid_model_molecule(imol)) {
      SCM chain_id_scm = scm_list_ref(residue_in, scm_from_int(0));
      SCM resno_scm    = scm_list_ref(residue_in, scm_from_int(1));
      SCM ins_code_scm = scm_list_ref(residue_in, scm_from_int(2));
      std::string chain_id = scm_to_locale_string(chain_id_scm);
      std::string ins_code = scm_to_locale_string(ins_code_scm);
      int resno            = scm_to_int(resno_scm);
      coot::residue_spec_t rspec(chain_id, resno, ins_code);
      std::vector<coot::residue_spec_t> v =
         graphics_info_t::molecules[imol].residues_near_residue(rspec, radius);
      for (unsigned int i=0; i<v.size(); i++) {
         SCM res_spec = SCM_EOL;
         res_spec = scm_cons(scm_from_locale_string(v[i].ins_code.c_str()), res_spec);
         res_spec = scm_cons(scm_from_int(v[i].res_no), res_spec);
         res_spec = scm_cons(scm_from_locale_string(v[i].chain_id.c_str()), res_spec);
         r = scm_cons(res_spec, r);
      }
   }
   return r;
}
#endif // USE_GUILE

//! \brief return residues near the given residues
//!
//! Return residue specs for residues that have atoms that are
//! closer than radius Angstroems to any atom in the residue
//! specified by res_in.
//!
#ifdef USE_GUILE
SCM residues_near_residues_scm(int imol, SCM residues_in_scm, float radius) {

   SCM r = SCM_EOL;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      if (scm_is_true(scm_list_p(residues_in_scm))) {
         SCM l_scm = scm_length(residues_in_scm);
         int l = scm_to_int(l_scm);
         std::vector<std::pair<bool, mmdb::Residue *> > res_vec;
         for (int i=0; i<l; i++) {
            SCM item_py = scm_list_ref(residues_in_scm, scm_from_int(i));
            coot::residue_spec_t spec = residue_spec_from_scm(item_py);
            mmdb::Residue *residue_p = graphics_info_t::molecules[imol].get_residue(spec);
            if (residue_p) {
               std::pair<bool, mmdb::Residue *> p(true, residue_p);
               res_vec.push_back(p);
            }
         }
         std::map<mmdb::Residue *, std::set<mmdb::Residue *> > rnrs =
            coot::residues_near_residues(res_vec, mol, radius);
         std::map<mmdb::Residue *, std::set<mmdb::Residue *> >::const_iterator it;
         for(it=rnrs.begin(); it!=rnrs.end(); it++) {
            mmdb::Residue *res_key = it->first;
            SCM key_scm = residue_spec_to_scm(coot::residue_spec_t(res_key));
            SCM value_scm = SCM_EOL;
            const std::set<mmdb::Residue *> &s = it->second;
            std::set<mmdb::Residue *>::const_iterator it_s;
            for (it_s=s.begin(); it_s!=s.end(); it_s++) {
               mmdb::Residue *r = *it_s;
               coot::residue_spec_t r_spec(r);
               SCM r_spec_scm = residue_spec_to_scm(r_spec);
               value_scm = scm_cons(r_spec_scm, value_scm);
            }
            SCM result_item_key_value_pair_scm = scm_list_2(key_scm, value_scm);
            r = scm_cons(result_item_key_value_pair_scm, r);
         }
      }
   }
   return r;
}
#endif // USE_GUILE


//! \brief
// Return residue specs for residues that have atoms that are
// closer than radius Angstroems to any atom in the residues
// specified by residues_in.
//
#ifdef USE_PYTHON
PyObject *residues_near_residues_py(int imol, PyObject *residues_in, float radius) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      if (PyList_Check(residues_in)) {
         int l = PyList_Size(residues_in);
         r = PyList_New(l);
         std::vector<std::pair<bool, mmdb::Residue *> > res_vec;
         for (int i=0; i<l; i++) {
            PyObject *item_py = PyList_GetItem(residues_in, i);
            coot::residue_spec_t spec = residue_spec_from_py(item_py);
            mmdb::Residue *residue_p = graphics_info_t::molecules[imol].get_residue(spec);
            if (residue_p) {
               std::pair<bool, mmdb::Residue *> p(true, residue_p);
               res_vec.push_back(p);
            }
         }
         std::map<mmdb::Residue *, std::set<mmdb::Residue *> > rnrs =
            coot::residues_near_residues(res_vec, mol, radius);
         std::map<mmdb::Residue *, std::set<mmdb::Residue *> >::const_iterator it;
         int map_idx = 0;
         for(it=rnrs.begin(); it!=rnrs.end(); it++) {
            mmdb::Residue *res_key = it->first;
            const std::set<mmdb::Residue *> &s = it->second;
            // FIXME
            // residue_spec_to_py adds a True prefix for a 4-item list
            // for historical reasons - get rid of that one day.
            PyObject *res_spec_py = residue_spec_to_py(coot::residue_spec_t(res_key));
            PyObject *key_py = residue_spec_make_triple_py(res_spec_py);
            PyObject *val_py = PyList_New(s.size());
            std::set<mmdb::Residue *>::const_iterator it_s;
            int idx = 0;
            for (it_s=s.begin(); it_s!=s.end(); it_s++) {
               mmdb::Residue *r = *it_s;
               coot::residue_spec_t r_spec(r);
               // like above, fiddle with the residue spec
               PyObject *r_4_py = residue_spec_to_py(r_spec);
               PyObject *r_3_py = residue_spec_make_triple_py(r_4_py);
               PyList_SetItem(val_py, idx, r_3_py);
               idx++;
            }
            PyObject *item_pair_py = PyList_New(2);
            PyList_SetItem(item_pair_py, 0, key_py);
            PyList_SetItem(item_pair_py, 1, val_py);
            PyList_SetItem(r, map_idx, item_pair_py);
            map_idx++;
         }
      }
   }

   if (PyBool_Check(r))
     Py_INCREF(r);

   return r;
}
#endif // USE_PYTHON


// Return residue specs for residues that have atoms that are
// closer than radius Angstroems to any atom in the residue
// specified by res_in.
//
// This will return a list of residue specifiers as triples
//
#ifdef USE_PYTHON
PyObject *residues_near_residue_py(int imol, PyObject *residue_spec_in, float radius) {

   PyObject *r = PyList_New(0);
   if (is_valid_model_molecule(imol)) {
      if (PyList_Check(residue_spec_in)) {
         std::pair<bool, coot::residue_spec_t> rspec =
            make_residue_spec_py(residue_spec_in);

         if (rspec.first) {

            std::vector<coot::residue_spec_t> v =
               graphics_info_t::molecules[imol].residues_near_residue(rspec.second, radius);
            for (unsigned int i=0; i<v.size(); i++) {
               PyObject *res_spec_py = residue_spec_to_py(v[i]);
               PyObject *res_spec_triple_py = residue_spec_make_triple_py(res_spec_py);
               PyList_Append(r, res_spec_triple_py);
               // Py_XDECREF(res_spec); - what did this do before I removed it?
            }
         } else {
            std::cout << "ERROR:: residues_near_residue_py() failed to construct "
                      << "residue spec" << std::endl;
         }
      } else {
         std::cout << "ERROR:: residues_near_residue_py() res_spec_in not a list"
                   << std::endl;
      }
   }
   return r;
}
#endif // USE_PYTHON

#ifdef USE_GUILE
SCM residues_near_position_scm(int imol, SCM pt_in_scm, float radius) {

   SCM r = SCM_EOL;

   if (is_valid_model_molecule(imol)) {

      SCM pt_in_length_scm = scm_length(pt_in_scm);
      int pt_in_length = scm_to_int(pt_in_length_scm);
      if (pt_in_length != 3) {
         std::cout << "WARNING:: input pt is not a list of 3 elements"
                   << std::endl;
      } else {

         double x = scm_to_double(scm_list_ref(pt_in_scm, scm_from_int(0)));
         double y = scm_to_double(scm_list_ref(pt_in_scm, scm_from_int(1)));
         double z = scm_to_double(scm_list_ref(pt_in_scm, scm_from_int(2)));

         clipper::Coord_orth pt(x,y,z);

         mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
         std::vector<mmdb::Residue *> v = coot::residues_near_position(pt, mol, radius);
         for (unsigned int i=0; i<v.size(); i++) {
            SCM r_scm = residue_spec_to_scm(coot::residue_spec_t(v[i]));
            r = scm_cons(r_scm, r);
         }
      }
   }
   return r;
}
#endif

#ifdef USE_PYTHON
PyObject *residues_near_position_py(int imol, PyObject *pt_in_py, float radius) {

   PyObject *r = PyList_New(0);

   if (is_valid_model_molecule(imol)) {

      int pt_in_length = PyObject_Length(pt_in_py);
      if (pt_in_length != 3) {
         std::cout << "WARNING:: input pt is not a list of 3 elements"
                   << std::endl;
      } else {

        double x = PyFloat_AsDouble(PyList_GetItem(pt_in_py, 0));
        double y = PyFloat_AsDouble(PyList_GetItem(pt_in_py, 1));
        double z = PyFloat_AsDouble(PyList_GetItem(pt_in_py, 2));

        clipper::Coord_orth pt(x,y,z);

        mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
        std::vector<mmdb::Residue *> v = coot::residues_near_position(pt, mol, radius);
        for (unsigned int i=0; i<v.size(); i++) {
          PyObject *r_py = residue_spec_to_py(coot::residue_spec_t(v[i]));
          PyList_Append(r, r_py);
          Py_XDECREF(r_py);
        }
      }
   }
   return r;
}
#endif

/*! \brief Label the atoms in the residues around the central residue */
void label_neighbours() {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      float radius = 4.0;
      int imol = pp.second.first;
      coot::residue_spec_t central_residue(pp.second.second);
      graphics_info_t g;
      g.molecules[imol].label_closest_atoms_in_neighbour_atoms(central_residue, radius);
      graphics_draw();
   }
}

/*! \brief Label the atoms in the central residue */
void label_atoms_in_residue() {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      graphics_info_t g;
      coot::residue_spec_t residue_spec(pp.second.second);
      mmdb::Residue *residue_p = g.molecules[imol].get_residue(residue_spec);
      if (residue_p) {
         g.molecules[imol].add_atom_labels_for_residue(residue_p);
         graphics_draw();
      }
   }
}

/*! \brief Label the atoms with their B-factors */
void set_show_local_b_factors(short int state) {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      graphics_info_t g;
      coot::Cartesian screen_centre = g.RotationCentre();
      g.molecules[imol].local_b_factor_display(state, screen_centre);
      graphics_draw();
   }
}



#include "c-interface-scm.hh"
#include "c-interface-python.hh"

#ifdef USE_GUILE
void label_closest_atoms_in_neighbour_residues_scm(int imol, SCM residue_spec_scm, float radius) {

   if (is_valid_model_molecule(imol)) {
      std::pair<bool, coot::residue_spec_t> res_spec = make_residue_spec(residue_spec_scm);
      if (res_spec.first) {
         graphics_info_t g;
         g.molecules[imol].label_closest_atoms_in_neighbour_atoms(res_spec.second, radius);
         graphics_draw();
      } else {
         std::cout << "WARNING:: bad spec " << std::endl;
      }
   }

}
#endif

#ifdef USE_PYTHON
void label_closest_atoms_in_neighbour_residues_py(int imol, PyObject *res_spec_py, float radius) {

   if (is_valid_model_molecule(imol)) {
      std::pair<bool, coot::residue_spec_t> res_spec = make_residue_spec_py(res_spec_py);
      if (res_spec.first) {
         graphics_info_t g;
         g.molecules[imol].label_closest_atoms_in_neighbour_atoms(res_spec.second, radius);
         graphics_draw();
      } else {
         std::cout << "WARNING:: bad spec " << std::endl;
      }
   }
}
#endif


//! find the active residue, find the near residues (within radius)
//! create a new molecule, run reduce on that, import hydrogens from
//! the result and apply them to the molecule of the active residue.
void hydrogenate_region(float radius) {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      coot::residue_spec_t central_residue(pp.second.second);
      std::cout << "----------- hydrogenating " << central_residue << " in " << imol << std::endl;
      coot::residue_spec_t res_spec(pp.second.second);
      std::vector<coot::residue_spec_t> v =
         graphics_info_t::molecules[imol].residues_near_residue(res_spec, radius);
      v.push_back(central_residue);
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      mmdb::Manager *new_mol = coot::util::create_mmdbmanager_from_residue_specs(v, mol);
      if (new_mol) {

         coot::util::create_directory("coot-molprobity"); // exists already maybe? Handled.
         std::string name_part = graphics_info_t::molecules[imol].Refmac_name_stub() + ".pdb";

         std::string pdb_in_file_name  = "hydrogenate-region-in-"  + name_part;
         std::string pdb_out_file_name = "hydrogenate-region-out-" + name_part;

         std::string pdb_in =  coot::util::append_dir_file("coot-molprobity", pdb_in_file_name);
         std::string pdb_out = coot::util::append_dir_file("coot-molprobity", pdb_out_file_name);

         new_mol->WritePDBASCII(pdb_in.c_str());

         if (graphics_info_t::prefer_python) {
#ifdef USE_PYTHON

            graphics_info_t g;
            short int lang = coot::STATE_PYTHON;
            std::string module = "generic_objects";
            std::string function = "reduce_on_pdb_file_no_flip";
            std::vector<coot::command_arg_t> args = {
               coot::command_arg_t(imol), pdb_in, pdb_out };
            std::string sc = g.state_command(module, function, args, lang);
            safe_python_command("import generic_objects");
            PyObject *r = safe_python_command_with_return(sc);
            std::cout << "::: A safe_python_command_with_return() returned " << r << std::endl;
            if (r)
               std::cout << "::: B safe_python_command_with_return() returned "
                         << PyBytes_AS_STRING(PyUnicode_AsUTF8String(display_python(r))) << std::endl;
            // 20230605-PE frustratingly the return value is None, even though I expect it to
            // be true. So just ignore this test for now.
            // if (r == Py_True) {
            if (true) {
               if (coot::file_exists(pdb_out)) {
                  std::cout << "DEBUG:: calling add_hydrogens_from_file() with pdb_out "
                            << pdb_out << std::endl;
                  graphics_info_t::molecules[imol].add_hydrogens_from_file(pdb_out);
               } else {
                  std::cout << "WARNING:: file does not exist " << pdb_out << std::endl;
               }
            }
            Py_XDECREF(r);

#endif // PYTHON
         } else {
#ifdef USE_GUILE
            // write a PDB file and run reduce, read it in
            //
            std::string scheme_command = "(reduce-on-pdb-file-no-flip ";
            scheme_command += coot::util::int_to_string(imol);
            scheme_command += " ";
            scheme_command += single_quote(pdb_in);
            scheme_command += " ";
            scheme_command += single_quote(pdb_out);
            scheme_command += ")";

            SCM r = safe_scheme_command(scheme_command);
            if (scm_is_true(r)) {
               graphics_info_t::molecules[imol].add_hydrogens_from_file(pdb_out);
            }
#endif
         }

         graphics_draw();
         delete new_mol;

      }
   }
}



#ifdef USE_GUILE
coot::residue_spec_t residue_spec_from_scm(SCM residue_in) {

   if (scm_is_true(scm_list_p(residue_in))) {
      SCM len_scm = scm_length(residue_in);
      int len = scm_to_int(len_scm);
      int offset = 0;
      if (len == 4)
         offset = 1;
      SCM chain_id_scm = scm_list_ref(residue_in, scm_from_int(0+offset));
      SCM resno_scm    = scm_list_ref(residue_in, scm_from_int(1+offset));
      SCM ins_code_scm = scm_list_ref(residue_in, scm_from_int(2+offset));
      std::string chain_id = scm_to_locale_string(chain_id_scm);
      std::string ins_code = scm_to_locale_string(ins_code_scm);
      int resno            = scm_to_int(resno_scm);
      coot::residue_spec_t rspec(chain_id, resno, ins_code);
      return rspec;
   } else {
      return coot::residue_spec_t();
   }
}
#endif // USE_GUILE

#ifdef USE_PYTHON
coot::residue_spec_t residue_spec_from_py(PyObject *residue_in) {

   // What about make_residue_spec_py()?

   coot::residue_spec_t rspec; // default
   int offset = 0;

   if (PyList_Check(residue_in)) {
      long len = PyList_Size(residue_in);
      if (len == 4)
         offset = 1;
      PyObject *chain_id_py = PyList_GetItem(residue_in, 0+offset);
      PyObject *resno_py    = PyList_GetItem(residue_in, 1+offset);
      PyObject *ins_code_py = PyList_GetItem(residue_in, 2+offset);
      if (PyUnicode_Check(chain_id_py)) {
         if (PyUnicode_Check(ins_code_py)) {
            if (PyLong_Check(resno_py)) {
               std::string chain_id  = PyBytes_AS_STRING(PyUnicode_AsUTF8String(chain_id_py));
               std::string ins_code  = PyBytes_AS_STRING(PyUnicode_AsUTF8String(ins_code_py));
               long resno            = PyLong_AsLong(resno_py);
               rspec = coot::residue_spec_t(chain_id, resno, ins_code);
               if (len == 4) {
                  PyObject *o = PyList_GetItem(residue_in, 0);
                  if (PyLong_Check(o)) {
                     long imol = PyLong_AsLong(o);
                     rspec.int_user_data = imol;
                  }
               }
               return rspec;
            }
         }
      }
   }
   return coot::residue_spec_t();

}
#endif // USE_PYTHON

#ifdef USE_GUILE
// output is like this:
// (list
//    (list (list atom-name alt-conf)
//          (list occ temp-fact element)
//          (list x y z)))
//
SCM residue_info(int imol, const char* chain_id, int resno, const char *ins_code) {

  SCM r = SCM_BOOL_F;
  if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      int imod = 1;

      mmdb::Model *model_p = mol->GetModel(imod);
      mmdb::Chain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
         chain_p = model_p->GetChain(ichain);
         std::string chain_id_mol(chain_p->GetChainID());
         if (chain_id_mol == std::string(chain_id)) {
            int nres = chain_p->GetNumberOfResidues();
            mmdb::PResidue residue_p;
            mmdb::Atom *at;

            for (int ires=0; ires<nres; ires++) {
               residue_p = chain_p->GetResidue(ires);
               std::string res_ins_code(residue_p->GetInsCode());
               if (residue_p->GetSeqNum() == resno) {
                  if (res_ins_code == std::string(ins_code)) {
                     int n_atoms = residue_p->GetNumberOfAtoms();
                     SCM at_info = SCM_BOOL(0);
                     SCM at_pos;
                     SCM at_occ, at_biso, at_ele, at_name, at_altconf;
                     SCM at_segid;
                     SCM at_x, at_y, at_z;
                     SCM all_atoms = SCM_EOL;
                     for (int iat=0; iat<n_atoms; iat++) {
                        at = residue_p->GetAtom(iat);
                        if (at->Ter) continue; //ignore TER records

                        at_x  = scm_from_double(at->x);
                        at_y  = scm_from_double(at->y);
                        at_z  = scm_from_double(at->z);
                        at_pos = scm_list_3(at_x, at_y, at_z);
                        at_occ = scm_from_double(at->occupancy);
                        at_biso= scm_from_double(at->tempFactor);
                        at_ele = scm_from_locale_string(at->element);
                        at_name = scm_from_locale_string(at->name);
                        at_segid = scm_from_locale_string(at->segID);
                        at_altconf = scm_from_locale_string(at->altLoc);
                        SCM at_b = at_biso;
                        if (at->WhatIsSet & mmdb::ASET_Anis_tFac) {
                           at_b = SCM_EOL;
                           at_b = scm_cons(at_biso, at_b);
                           at_b = scm_cons(scm_from_double(at->u11), at_b);
                           at_b = scm_cons(scm_from_double(at->u22), at_b);
                           at_b = scm_cons(scm_from_double(at->u33), at_b);
                           at_b = scm_cons(scm_from_double(at->u12), at_b);
                           at_b = scm_cons(scm_from_double(at->u13), at_b);
                           at_b = scm_cons(scm_from_double(at->u23), at_b);
                           at_b = scm_reverse(at_b);
                        }
                        SCM compound_name = scm_list_2(at_name, at_altconf);
                        SCM compound_attrib = scm_list_4(at_occ, at_b, at_ele, at_segid);
                        at_info = scm_list_3(compound_name, compound_attrib, at_pos);
                        all_atoms = scm_cons(at_info, all_atoms);
                     }
                     r = scm_reverse(all_atoms);
                  }
               }
            }
         }
      }
   }
   return r;
}
#endif // USE_GUILE

#ifdef USE_PYTHON
PyObject *residue_info_py(int imol, const char* chain_id, int resno, const char *ins_code) {

   PyObject *r = Py_False;
   PyObject *all_atoms;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      int udd_handle     = graphics_info_t::molecules[imol].atom_sel.UDDAtomIndexHandle;
      int imod = 1;

      mmdb::Model *model_p = mol->GetModel(imod);
      mmdb::Chain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
         chain_p = model_p->GetChain(ichain);
         std::string chain_id_mol(chain_p->GetChainID());
         if (chain_id_mol == std::string(chain_id)) {
            int nres = chain_p->GetNumberOfResidues();
            mmdb::PResidue residue_p;
            mmdb::Atom *at;

            // why use this bizarre contrivance to get a null list for
            // starting? I must be missing something.
            for (int ires=0; ires<nres; ires++) {
               residue_p = chain_p->GetResidue(ires);
               std::string res_ins_code(residue_p->GetInsCode());
               if (residue_p->GetSeqNum() == resno) {
                  if (res_ins_code == std::string(ins_code)) {
                     int n_atoms = residue_p->GetNumberOfAtoms();
                     PyObject *at_info = Py_False;
                     PyObject *at_pos;
                     PyObject *at_occ, *at_b, *at_biso, *at_ele, *at_name, *at_altconf;
                     PyObject *at_segid;
                     PyObject *at_x, *at_y, *at_z;
                     PyObject *compound_name;
                     PyObject *compound_attrib;
                     all_atoms = PyList_New(0);

                     for (int iat=0; iat<n_atoms; iat++) {

                        at = residue_p->GetAtom(iat);
                        if (at->Ter) continue; //ignore TER records

                        int idx = -1;
                        int ierr = at->GetUDData(udd_handle, idx);
                        if (ierr != mmdb::UDDATA_Ok) {
                           std::cout << "WARNING:: residue_info_py(): error getting uddata for atom "
                                     << at << std::endl;
                           idx = -1; // maybe not needed
                        }
                        PyObject *atom_idx_py = PyLong_FromLong(idx);

                        at_x  = PyFloat_FromDouble(at->x);
                        at_y  = PyFloat_FromDouble(at->y);
                        at_z  = PyFloat_FromDouble(at->z);
                        at_pos = PyList_New(3);
                        PyList_SetItem(at_pos, 0, at_x);
                        PyList_SetItem(at_pos, 1, at_y);
                        PyList_SetItem(at_pos, 2, at_z);

                        at_occ = PyFloat_FromDouble(at->occupancy);
                        at_biso= PyFloat_FromDouble(at->tempFactor);
                        at_ele = myPyString_FromString(at->element);
                        at_name = myPyString_FromString(at->name);
                        at_segid = myPyString_FromString(at->segID);
                        at_altconf = myPyString_FromString(at->altLoc);

                        at_b = at_biso;
                        if (at->WhatIsSet & mmdb::ASET_Anis_tFac) {
                           at_b = PyList_New(7);
                           PyList_SetItem(at_b, 0, at_biso);
                           PyList_SetItem(at_b, 1, PyFloat_FromDouble(at->u11));
                           PyList_SetItem(at_b, 2, PyFloat_FromDouble(at->u22));
                           PyList_SetItem(at_b, 3, PyFloat_FromDouble(at->u33));
                           PyList_SetItem(at_b, 4, PyFloat_FromDouble(at->u12));
                           PyList_SetItem(at_b, 5, PyFloat_FromDouble(at->u13));
                           PyList_SetItem(at_b, 6, PyFloat_FromDouble(at->u23));
                        }

                        compound_name = PyList_New(2);
                        PyList_SetItem(compound_name, 0 ,at_name);
                        PyList_SetItem(compound_name, 1 ,at_altconf);

                        compound_attrib = PyList_New(4);
                        PyList_SetItem(compound_attrib, 0, at_occ);
                        PyList_SetItem(compound_attrib, 1, at_b);
                        PyList_SetItem(compound_attrib, 2, at_ele);
                        PyList_SetItem(compound_attrib, 3, at_segid);

                        at_info = PyList_New(4);
                        PyList_SetItem(at_info, 0, compound_name);
                        PyList_SetItem(at_info, 1, compound_attrib);
                        PyList_SetItem(at_info, 2, at_pos);
                        PyList_SetItem(at_info, 3, atom_idx_py);

                        PyList_Append(all_atoms, at_info);
                     }
                     r = all_atoms;
                  }
               }
            }
         }
      }
   }

   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }

   return r;

}

#endif //PYTHON

std::string residue_name(int imol, const std::string &chain_id, int resno,
                         const std::string &ins_code) {

   std::string r = "";
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
         bool have_resname_flag = 0;

         mmdb::Model *model_p = mol->GetModel(imod);

         if (model_p) {
            // run over chains of the existing mol
            int nchains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<nchains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               std::string chain_id_mol(chain_p->GetChainID());
               if (chain_id_mol == std::string(chain_id)) {
                  int nres = chain_p->GetNumberOfResidues();
                  for (int ires=0; ires<nres; ires++) {
                     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                     if (residue_p->GetSeqNum() == resno) {
                        std::string ins = residue_p->GetInsCode();
                        if (ins == ins_code) {
                           r = residue_p->GetResName();
                           have_resname_flag = 1;
                           break;
                        }
                     }
                  }
               }
            }
            if (have_resname_flag)
               break;
         }
         if (have_resname_flag)
            break;
      }
   }
   return r;
}


#ifdef USE_GUILE
SCM residue_name_scm(int imol, const char* chain_id, int resno, const char *ins_code) {

   SCM r = SCM_BOOL(0);
   std::string res_name = residue_name(imol, chain_id, resno, ins_code);
   if (res_name.size() > 0) {
      r = scm_from_locale_string(res_name.c_str());
   }
   return r;
}
#endif // USE_GUILE

#ifdef USE_PYTHON
PyObject *residue_name_py(int imol, const char* chain_id, int resno, const char *ins_code) {

   PyObject *r;
   r = Py_False;
   std::string res_name = residue_name(imol, chain_id, resno, ins_code);
   if (res_name.size() > 0) {
      r = myPyString_FromString(res_name.c_str());
   }

   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif // USE_PYTHON


#ifdef USE_GUILE
SCM chain_fragments_scm(int imol, short int screen_output_also) {

   SCM r = SCM_BOOL_F;

   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      std::vector<coot::fragment_info_t> f = g.molecules[imol].get_fragment_info(screen_output_also);
   }
   return r;
}
#endif // USE_GUILE

#ifdef USE_PYTHON
PyObject *chain_fragments_py(int imol, short int screen_output_also) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      std::vector<coot::fragment_info_t> f = g.molecules[imol].get_fragment_info(screen_output_also);
   }

   if (PyBool_Check(r)) {
      Py_INCREF(r);
   }

   return r;
}
#endif // USE_PYTHON


void set_rotation_centre(const clipper::Coord_orth &pos) {
   graphics_info_t g;
   g.setRotationCentre(coot::Cartesian(pos));
}


#ifdef USE_GUILE
// Bernie, no need to pythonize this, it's just to test the return
// values on pressing "next residue" and "previous residue" (you can
// if you wish of course).
//
// Pass the current values, return new values
//
// Hmm maybe these function should pass the atom name too?  Yes they should
SCM goto_next_atom_maybe_scm(const char *chain_id, int resno, const char *ins_code,
                             const char *atom_name) {

   SCM r = SCM_BOOL_F;

   int imol = go_to_atom_molecule_number();
   if (is_valid_model_molecule(imol)) {

      graphics_info_t g;
      coot::Cartesian rc = g.RotationCentre();

      int atom_index =
         graphics_info_t::molecules[imol].intelligent_next_atom(chain_id, resno,
                                                                atom_name, ins_code, rc);

      if (atom_index != -1) {
         mmdb::Atom *next_atom = graphics_info_t::molecules[imol].atom_sel.atom_selection[atom_index];

         std::string next_chain_id  = next_atom->GetChainID();
         std::string next_atom_name = next_atom->name;
         int next_residue_number    = next_atom->GetSeqNum();
         std::string next_ins_code  = next_atom->GetInsCode();

         r = SCM_EOL;
         r = scm_cons(scm_from_locale_string(next_atom_name.c_str()), r);
         r = scm_cons(scm_from_locale_string(next_ins_code.c_str()), r);
         r = scm_cons(scm_from_int(next_residue_number) ,r);
         r = scm_cons(scm_from_locale_string(next_chain_id.c_str()), r);
      }
   }
   return r;
}
#endif

#ifdef USE_GUILE
SCM goto_prev_atom_maybe_scm(const char *chain_id, int resno, const char *ins_code,
                             const char *atom_name) {

   SCM r = SCM_BOOL_F;
   int imol = go_to_atom_molecule_number();
   if (is_valid_model_molecule(imol)) {

      graphics_info_t g;
      coot::Cartesian rc = g.RotationCentre();
      int atom_index =
         graphics_info_t::molecules[imol].intelligent_previous_atom(chain_id, resno,
                                                                    atom_name, ins_code, rc);

      if (atom_index != -1) {
         mmdb::Atom *next_atom = graphics_info_t::molecules[imol].atom_sel.atom_selection[atom_index];

         std::string next_chain_id  = next_atom->GetChainID();
         std::string next_atom_name = next_atom->name;
         int next_residue_number    = next_atom->GetSeqNum();
         std::string next_ins_code  = next_atom->GetInsCode();

         r = SCM_EOL;
         r = scm_cons(scm_from_locale_string(next_atom_name.c_str()), r);
         r = scm_cons(scm_from_locale_string(next_ins_code.c_str()), r);
         r = scm_cons(scm_from_int(next_residue_number) ,r);
         r = scm_cons(scm_from_locale_string(next_chain_id.c_str()), r);
      }
   }
   return r;
}
#endif

#ifdef USE_PYTHON
// Pass the current values, return new values
//
// Hmm maybe these function should pass the atom name too?  Yes they should
PyObject *goto_next_atom_maybe_py(const char *chain_id, int resno, const char *ins_code,
                                  const char *atom_name) {

   PyObject *r = Py_False;

   int imol = go_to_atom_molecule_number();
   if (is_valid_model_molecule(imol)) {

      graphics_info_t g;
      coot::Cartesian rc = g.RotationCentre();
      int atom_index =
         graphics_info_t::molecules[imol].intelligent_next_atom(chain_id, resno,
                                                                    atom_name, ins_code, rc);
      if (atom_index != -1) {
         mmdb::Atom *next_atom = graphics_info_t::molecules[imol].atom_sel.atom_selection[atom_index];

         std::string next_chain_id  = next_atom->GetChainID();
         std::string next_atom_name = next_atom->name;
         int next_residue_number    = next_atom->GetSeqNum();
         std::string next_ins_code  = next_atom->GetInsCode();

         r = PyList_New(4);
         PyList_SetItem(r, 0, myPyString_FromString(next_chain_id.c_str()));
         PyList_SetItem(r, 1, PyLong_FromLong(next_residue_number));
         PyList_SetItem(r, 2, myPyString_FromString(next_ins_code.c_str()));
         PyList_SetItem(r, 3, myPyString_FromString(next_atom_name.c_str()));
      }
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif // USE_PYTHON

#ifdef USE_PYTHON
PyObject *goto_prev_atom_maybe_py(const char *chain_id, int resno, const char *ins_code,
                                  const char *atom_name) {

   PyObject *r = Py_False;
   int imol = go_to_atom_molecule_number();
   if (is_valid_model_molecule(imol)) {

      graphics_info_t g;
      coot::Cartesian rc = g.RotationCentre();
      int atom_index =
         graphics_info_t::molecules[imol].intelligent_previous_atom(chain_id, resno,
                                                                    atom_name, ins_code, rc);

      if (atom_index != -1) {
         mmdb::Atom *next_atom = graphics_info_t::molecules[imol].atom_sel.atom_selection[atom_index];

         std::string next_chain_id  = next_atom->GetChainID();
         std::string next_atom_name = next_atom->name;
         int next_residue_number    = next_atom->GetSeqNum();
         std::string next_ins_code  = next_atom->GetInsCode();

         r = PyList_New(4);
         PyList_SetItem(r, 0, myPyString_FromString(next_chain_id.c_str()));
         PyList_SetItem(r, 1, PyLong_FromLong(next_residue_number));
         PyList_SetItem(r, 2, myPyString_FromString(next_ins_code.c_str()));
         PyList_SetItem(r, 3, myPyString_FromString(next_atom_name.c_str()));
      }
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif // USE_PYTHON

// A C++ function interface:
//
int set_go_to_atom_from_spec(const coot::atom_spec_t &atom_spec) {

   int success = 0;
   graphics_info_t g;

   if (! atom_spec.empty()) {
      g.set_go_to_atom_chain_residue_atom_name(atom_spec.chain_id.c_str(),
                                               atom_spec.res_no,
                                               atom_spec.ins_code.c_str(),
                                               atom_spec.atom_name.c_str(),
                                               atom_spec.alt_conf.c_str());

      success = g.try_centre_from_new_go_to_atom();
      if (success)
         g.update_things_on_move_and_redraw();
   }
   return success;
}

int set_go_to_atom_from_res_spec(const coot::residue_spec_t &spec) {

   int success = 0;
   graphics_info_t g;
   int imol = g.go_to_atom_molecule();

   if (is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec = g.molecules[imol].intelligent_this_residue_atom(spec);
      if (! atom_spec.empty()) {
         success = set_go_to_atom_from_spec(atom_spec);
      }
   }
   return success;
}

#ifdef USE_GUILE
int set_go_to_atom_from_res_spec_scm(SCM residue_spec_scm) {

   coot::residue_spec_t spec = residue_spec_from_scm(residue_spec_scm);
   return set_go_to_atom_from_res_spec(spec);
}

#ifdef USE_GUILE
int set_go_to_atom_from_atom_spec_scm(SCM atom_spec_scm) {
   coot::atom_spec_t spec = atom_spec_from_scm_expression(atom_spec_scm);
   return set_go_to_atom_from_spec(spec);
}
#endif // USE_GUILE


#endif

#ifdef USE_PYTHON
int set_go_to_atom_from_res_spec_py(PyObject *residue_spec_py) {

   std::pair<bool, coot::residue_spec_t> spec = make_residue_spec_py(residue_spec_py);
   if (spec.first) {
      return set_go_to_atom_from_res_spec(spec.second);
   } else {
      return -1;
   }
}
#endif


#ifdef USE_PYTHON
int set_go_to_atom_from_atom_spec_py(PyObject *atom_spec_py) {

   coot::atom_spec_t spec = atom_spec_from_python_expression(atom_spec_py);
   return set_go_to_atom_from_spec(spec);
}
#endif


// (is-it-valid? (active-molecule-number spec))
std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom_spec() {

   return graphics_info_t::active_atom_spec();
}

#ifdef USE_PYTHON
// return a tuple of (Py_Bool (number, atom_spec))
PyObject *active_atom_spec_py() {

   PyObject *rv = PyTuple_New(2);

   std::pair<bool, std::pair<int, coot::atom_spec_t> > r = active_atom_spec();
      PyObject *state_py = Py_True;
      PyObject *mol_no = PyLong_FromLong(r.second.first);
      PyObject *spec = residue_spec_to_py(coot::residue_spec_t(r.second.second));
      PyObject *tuple_inner = PyTuple_New(2);
   if (! r.first) {
      state_py = Py_False;
   }
   Py_INCREF(state_py);

   PyTuple_SetItem(tuple_inner, 0, mol_no);
   PyTuple_SetItem(tuple_inner, 1, spec);
   PyTuple_SetItem(rv, 0, state_py);
   PyTuple_SetItem(rv, 1, tuple_inner);
   return rv;
}
#endif // USE_PYTHON


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
      s = scm_cons(scm_from_locale_string(pp.second.second.alt_conf.c_str()) , s);
      s = scm_cons(scm_from_locale_string(pp.second.second.atom_name.c_str()), s);
      s = scm_cons(scm_from_locale_string(pp.second.second.ins_code.c_str()) , s);
      s = scm_cons(scm_from_int(pp.second.second.res_no) , s);
      s = scm_cons(scm_from_locale_string(pp.second.second.chain_id.c_str()) , s);
      s = scm_cons(scm_from_int(pp.second.first) ,s);
   }
   return s;
}
#endif // USE_GUILE

#ifdef USE_GUILE
//! \brief return the specs of the closest displayed atom
//!
//! Return a list of (list imol chain-id resno ins-code atom-name
//! alt-conf (list x y z)) for atom that is closest to the screen
//! centre in the given molecule (unlike active-residue, potentila CA
//! substition is not performed).  If there is no atom, or if imol is
//! not a valid model molecule, return scheme false.
//!
SCM closest_atom_simple_scm() {

   SCM s = SCM_BOOL(0);
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = graphics_info_t::active_atom_spec_simple();

   if (pp.first) {
      s = SCM_EOL;
      s = scm_cons(scm_from_locale_string(pp.second.second.alt_conf.c_str()), s);
      s = scm_cons(scm_from_locale_string(pp.second.second.atom_name.c_str()), s);
      s = scm_cons(scm_from_locale_string(pp.second.second.ins_code.c_str()), s);
      s = scm_cons(scm_from_int(pp.second.second.res_no), s);
      s = scm_cons(scm_from_locale_string(pp.second.second.chain_id.c_str()), s);
      s = scm_cons(scm_from_int(pp.second.first), s);
   }
   return s;
}
#endif // USE_GUILE


#ifdef USE_PYTHON
PyObject *active_residue_py() {

   PyObject *s;
   s = Py_False;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();

   if (pp.first) {
      s = PyList_New(6);
      PyList_SetItem(s, 0, PyLong_FromLong(pp.second.first));
      PyList_SetItem(s, 1, myPyString_FromString(pp.second.second.chain_id.c_str()));
      PyList_SetItem(s, 2, PyLong_FromLong(pp.second.second.res_no));
      PyList_SetItem(s, 3, myPyString_FromString(pp.second.second.ins_code.c_str()));
      PyList_SetItem(s, 4, myPyString_FromString(pp.second.second.atom_name.c_str()));
      PyList_SetItem(s, 5, myPyString_FromString(pp.second.second.alt_conf.c_str()));
   }
   if (PyBool_Check(s)) {
     Py_INCREF(s);
   }
   return s;
}
#endif // PYTHON

#ifdef USE_PYTHON
//! \brief return the specs of the closest displayed atom
//!
//! Return a list of (list imol chain-id resno ins-code atom-name
//! alt-conf (list x y z)) for atom that is closest to the screen
//! centre in the given molecule (unlike active-residue, potentila CA
//! substition is not performed).  If there is no atom, or if imol is
//! not a valid model molecule, return scheme false.
//!
PyObject *closest_atom_simple_py() {

   PyObject *s = Py_False;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = graphics_info_t::active_atom_spec_simple();
   if (pp.first) {
      s = PyList_New(6);
      PyList_SetItem(s, 0, PyLong_FromLong(pp.second.first));
      PyList_SetItem(s, 1, myPyString_FromString(pp.second.second.chain_id.c_str()));
      PyList_SetItem(s, 2, PyLong_FromLong(pp.second.second.res_no));
      PyList_SetItem(s, 3, myPyString_FromString(pp.second.second.ins_code.c_str()));
      PyList_SetItem(s, 4, myPyString_FromString(pp.second.second.atom_name.c_str()));
      PyList_SetItem(s, 5, myPyString_FromString(pp.second.second.alt_conf.c_str()));
   }
   if (PyBool_Check(s)) {
     Py_INCREF(s);
   }
   return s;
}
#endif // USE_PYTHON


#ifdef USE_GUILE
SCM closest_atom(int imol) {

   SCM r = SCM_BOOL_F;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      coot::at_dist_info_t at_info =
         graphics_info_t::molecules[imol].closest_atom(g.RotationCentre());
      if (at_info.atom) {
         r = SCM_EOL;
         r = scm_cons(scm_from_double(at_info.atom->z), r);
         r = scm_cons(scm_from_double(at_info.atom->y), r);
         r = scm_cons(scm_from_double(at_info.atom->x), r);
         r = scm_cons(scm_from_locale_string(at_info.atom->altLoc), r);
         r = scm_cons(scm_from_locale_string(at_info.atom->name), r);
         r = scm_cons(scm_from_locale_string(at_info.atom->GetInsCode()), r);
         r = scm_cons(scm_from_int(at_info.atom->GetSeqNum()), r);
         r = scm_cons(scm_from_locale_string(at_info.atom->GetChainID()), r);
         r = scm_cons(scm_from_int(imol), r);
      }
   }
   return r;
}
#endif

#ifdef USE_PYTHON
PyObject *closest_atom_py(int imol) {

   PyObject *r;
   r = Py_False;
   if (is_valid_model_molecule(imol)) {
      graphics_info_t g;
      coot::at_dist_info_t at_info =
         graphics_info_t::molecules[imol].closest_atom(g.RotationCentre());
      if (at_info.atom) {
         r = PyList_New(9);
         PyList_SetItem(r, 0, PyLong_FromLong(imol));
         PyList_SetItem(r, 1, myPyString_FromString(at_info.atom->GetChainID()));
         PyList_SetItem(r, 2, PyLong_FromLong(at_info.atom->GetSeqNum()));
         PyList_SetItem(r, 3, myPyString_FromString(at_info.atom->GetInsCode()));
         PyList_SetItem(r, 4, myPyString_FromString(at_info.atom->name));
         PyList_SetItem(r, 5, myPyString_FromString(at_info.atom->altLoc));
         PyList_SetItem(r, 6, PyFloat_FromDouble(at_info.atom->x));
         PyList_SetItem(r, 7, PyFloat_FromDouble(at_info.atom->y));
         PyList_SetItem(r, 8, PyFloat_FromDouble(at_info.atom->z));
      }
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif // USE_PYTHON

#ifdef USE_GUILE
SCM closest_atom_raw_scm() {

   SCM r = SCM_BOOL_F;
   graphics_info_t g;
   // index, imol
   std::pair<int, int> ca_ii = g.get_closest_atom();

   int imol = ca_ii.second;
   if (is_valid_model_molecule(imol)) {
      mmdb::Atom *at = g.molecules[imol].get_atom(ca_ii.first);
      if (at) {
         r = SCM_EOL;
         r = scm_cons(scm_from_double(at->z), r);
         r = scm_cons(scm_from_double(at->y), r);
         r = scm_cons(scm_from_double(at->x), r);
         r = scm_cons(scm_from_locale_string(at->altLoc), r);
         r = scm_cons(scm_from_locale_string(at->name), r);
         r = scm_cons(scm_from_locale_string(at->GetInsCode()), r);
         r = scm_cons(scm_from_int(at->GetSeqNum()), r);
         r = scm_cons(scm_from_locale_string(at->GetChainID()), r);
         r = scm_cons(scm_from_int(imol), r);
      }
   }
   return r;
}
#endif

#ifdef USE_PYTHON
PyObject *closest_atom_raw_py() {

   PyObject *r;
   r = Py_False;
   graphics_info_t g;
   // index, imol
   std::pair<int, int> ca_ii = g.get_closest_atom();

   int imol = ca_ii.second;
   if (is_valid_model_molecule(imol)) {

      mmdb::Atom *at = g.molecules[imol].get_atom(ca_ii.first);
      if (at) {
         r = PyList_New(9);
         PyList_SetItem(r, 0, PyLong_FromLong(imol));
         PyList_SetItem(r, 1, myPyString_FromString(at->GetChainID()));
         PyList_SetItem(r, 2, PyLong_FromLong(at->GetSeqNum()));
         PyList_SetItem(r, 3, myPyString_FromString(at->GetInsCode()));
         PyList_SetItem(r, 4, myPyString_FromString(at->name));
         PyList_SetItem(r, 5, myPyString_FromString(at->altLoc));
         PyList_SetItem(r, 6, PyFloat_FromDouble(at->x));
         PyList_SetItem(r, 7, PyFloat_FromDouble(at->y));
         PyList_SetItem(r, 8, PyFloat_FromDouble(at->z));
      }
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif

/*! \brief update the Go To Atom widget entries to atom closest to
  screen centre. */
void update_go_to_atom_from_current_position() {

   // move this function into graphics_info_t I think

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
   if (pp.first) {
      int imol = pp.second.first;
      set_go_to_atom_molecule(imol);
      const coot::atom_spec_t &atom_spec = pp.second.second;
      set_go_to_atom_chain_residue_atom_name(atom_spec.chain_id.c_str(),
                                             atom_spec.res_no,
                                             atom_spec.atom_name.c_str());
      update_go_to_atom_window_on_other_molecule_chosen(pp.second.first);

      graphics_info_t g;
      std::cout << "if sequence view is displayed update highlighted position here A " << std::endl;
      // now run a graphics_info_t function.

   }
}

void
update_sequence_view_current_position_highlight_from_active_atom() {

   std::cout << "if sequence view is displayed update highlighted position here B " << std::endl;
   // now run a graphics_info_t function.
}


#ifdef USE_GUILE
SCM generic_string_vector_to_list_internal(const std::vector<std::string> &v) {

   SCM r = SCM_EOL;
   for (int i=v.size()-1; i>=0; i--) {
      r = scm_cons(scm_from_locale_string(v[i].c_str()), r);
   }
   return r;
}
#endif // USE_GUILE

// BL says:: python version
#ifdef USE_PYTHON
PyObject *generic_string_vector_to_list_internal_py(const std::vector<std::string> &v) {

  PyObject *r;

   r = PyList_New(v.size());
   for (int i=v.size()-1; i>=0; i--) {
      PyList_SetItem(r, i, myPyString_FromString(v[i].c_str()));
   }
   return r;
}
#endif // PYTHON

// and the reverse function:
#ifdef USE_GUILE
std::vector<std::string>
generic_list_to_string_vector_internal(SCM l) {
   std::vector<std::string> r;
   SCM l_length_scm = scm_length(l);

#if (SCM_MAJOR_VERSION > 1) || (SCM_MINOR_VERSION > 7)

   int l_length = scm_to_int(l_length_scm);
   for (int i=0; i<l_length; i++) {
      SCM le = scm_list_ref(l, scm_from_int(i));
      std::string s = scm_to_locale_string(le);
      r.push_back(s);
   }

#else

   int l_length = gh_scm2int(l_length_scm);
   for (int i=0; i<l_length; i++) {
      SCM le = scm_list_ref(l, scm_from_int(i));
      std::string s = SCM_STRING_CHARS(le);
      r.push_back(s);
   }

#endif

   return r;
}
#endif

#ifdef USE_PYTHON
std::vector<std::string>
generic_list_to_string_vector_internal_py(PyObject *l) {
   std::vector<std::string> r;

   int l_length = PyObject_Length(l);
   for (int i=0; i<l_length; i++) {
      PyObject *le = PyList_GetItem(l, i);
      std::string s = PyBytes_AS_STRING(PyUnicode_AsUTF8String(le));
      r.push_back(s);
   }

   return r;
}

#endif // USE_PYTHON

#ifdef USE_GUILE
SCM generic_int_vector_to_list_internal(const std::vector<int> &v) {

   SCM r = SCM_EOL;
   for (int i=v.size()-1; i>=0; i--) {
      r = scm_cons(scm_from_int(v[i]), r);
   }
   return r;
}
#endif // USE_GUILE

#ifdef USE_PYTHON
PyObject *generic_int_vector_to_list_internal_py(const std::vector<int> &v) {

   PyObject *r;
   r = PyList_New(v.size());
   for (int i=v.size()-1; i>=0; i--) {
      PyList_SetItem(r, i, PyLong_FromLong(v[i]));
   }
   return r;
}
#endif // USE_PYTHON

#ifdef USE_GUILE
SCM rtop_to_scm(const clipper::RTop_orth &rtop) {

   SCM r = SCM_EOL;

   SCM tr_list = SCM_EOL;
   SCM rot_list = SCM_EOL;

   clipper::Mat33<double>  mat = rtop.rot();
   clipper::Vec3<double> trans = rtop.trn();

   tr_list = scm_cons(scm_from_double(trans[2]), tr_list);
   tr_list = scm_cons(scm_from_double(trans[1]), tr_list);
   tr_list = scm_cons(scm_from_double(trans[0]), tr_list);

   rot_list = scm_cons(scm_from_double(mat(2,2)), rot_list);
   rot_list = scm_cons(scm_from_double(mat(2,1)), rot_list);
   rot_list = scm_cons(scm_from_double(mat(2,0)), rot_list);
   rot_list = scm_cons(scm_from_double(mat(1,2)), rot_list);
   rot_list = scm_cons(scm_from_double(mat(1,1)), rot_list);
   rot_list = scm_cons(scm_from_double(mat(1,0)), rot_list);
   rot_list = scm_cons(scm_from_double(mat(0,2)), rot_list);
   rot_list = scm_cons(scm_from_double(mat(0,1)), rot_list);
   rot_list = scm_cons(scm_from_double(mat(0,0)), rot_list);

   r = scm_cons(tr_list, r);
   r = scm_cons(rot_list, r);
   return r;

}
#endif // USE_GUILE

#ifdef  USE_GUILE
SCM inverse_rtop_scm(SCM rtop_scm) {

   clipper::RTop_orth r;
   SCM rtop_len_scm = scm_length(rtop_scm);
   int rtop_len = scm_to_int(rtop_len_scm);
   if (rtop_len == 2) {
      SCM rot_scm = scm_list_ref(rtop_scm, scm_from_int(0));
      SCM rot_length_scm = scm_length(rot_scm);
      int rot_length = scm_to_int(rot_length_scm);
      if (rot_length == 9) {
         SCM trn_scm = scm_list_ref(rtop_scm, scm_from_int(1));
         SCM trn_length_scm = scm_length(trn_scm);
         int trn_length = scm_to_int(trn_length_scm);
         double rot_arr[9];
         double trn_arr[3];
         if (trn_length == 3) {
            for (int i=0; i<9; i++) {
               rot_arr[i] = scm_to_double(scm_list_ref(rot_scm, scm_from_int(i)));
            }
            for (int i=0; i<3; i++) {
               trn_arr[i] = scm_to_double(scm_list_ref(trn_scm, scm_from_int(i)));
            }
            clipper::Mat33<double> rot(rot_arr[0], rot_arr[1], rot_arr[2],
                                       rot_arr[3], rot_arr[4], rot_arr[5],
                                       rot_arr[6], rot_arr[7], rot_arr[8]);
            clipper::Coord_orth trn(trn_arr[0], trn_arr[1], trn_arr[2]);
            clipper::RTop_orth rtop_in(rot, trn);
            r = rtop_in.inverse();
         }
      }
   }
   return rtop_to_scm(r);

}
#endif // USE_GUILE


#ifdef USE_PYTHON
PyObject *rtop_to_python(const clipper::RTop_orth &rtop) {

   PyObject *r;
   PyObject *tr_list;
   PyObject *rot_list;

   r = PyList_New(2);
   tr_list = PyList_New(3);
   rot_list = PyList_New(9);

   clipper::Mat33<double>  mat = rtop.rot();
   clipper::Vec3<double> trans = rtop.trn();

   PyList_SetItem(tr_list, 0, PyFloat_FromDouble(trans[0]));
   PyList_SetItem(tr_list, 1, PyFloat_FromDouble(trans[1]));
   PyList_SetItem(tr_list, 2, PyFloat_FromDouble(trans[2]));

   PyList_SetItem(rot_list, 0, PyFloat_FromDouble(mat(0,0)));
   PyList_SetItem(rot_list, 1, PyFloat_FromDouble(mat(0,1)));
   PyList_SetItem(rot_list, 2, PyFloat_FromDouble(mat(0,2)));
   PyList_SetItem(rot_list, 3, PyFloat_FromDouble(mat(1,0)));
   PyList_SetItem(rot_list, 4, PyFloat_FromDouble(mat(1,1)));
   PyList_SetItem(rot_list, 5, PyFloat_FromDouble(mat(1,2)));
   PyList_SetItem(rot_list, 6, PyFloat_FromDouble(mat(2,0)));
   PyList_SetItem(rot_list, 7, PyFloat_FromDouble(mat(2,1)));
   PyList_SetItem(rot_list, 8, PyFloat_FromDouble(mat(2,2)));

// BL says:: or maybe the other way round
   PyList_SetItem(r, 0, rot_list);
   PyList_SetItem(r, 1, tr_list);
   return r;

}

PyObject *inverse_rtop_py(PyObject *rtop_py) {

   clipper::RTop_orth r;
   int rtop_len = PyList_Size(rtop_py);
   if (rtop_len == 2) {
      PyObject *rot_py = PyList_GetItem(rtop_py, 0);
      int rot_length = PyList_Size(rot_py);
      if (rot_length == 9) {
         PyObject *trn_py = PyList_GetItem(rtop_py, 1);
         int trn_length = PyList_Size(trn_py);
         double rot_arr[9];
         double trn_arr[3];
         if (trn_length == 3) {
            for (int i=0; i<9; i++) {
              rot_arr[i] = PyFloat_AsDouble(PyList_GetItem(rot_py, i));
            }
            for (int i=0; i<3; i++) {
              trn_arr[i] = PyFloat_AsDouble(PyList_GetItem(trn_py, i));
            }
            clipper::Mat33<double> rot(rot_arr[0], rot_arr[1], rot_arr[2],
                                       rot_arr[3], rot_arr[4], rot_arr[5],
                                       rot_arr[6], rot_arr[7], rot_arr[8]);
            clipper::Coord_orth trn(trn_arr[0], trn_arr[1], trn_arr[2]);
            clipper::RTop_orth rtop_in(rot, trn);
            r = rtop_in.inverse();
         }
      }
   }
   return rtop_to_python(r);

}
#endif // USE_PYTHON

#ifdef USE_GUILE
// get the symmetry operators strings for the given molecule
//
/*! \brief Return as a list of strings the symmetry operators of the
  given molecule. If imol is a not a valid molecule, return an empty
  list.*/
//
SCM get_symmetry(int imol) {

   SCM r = SCM_EOL;
   if (is_valid_model_molecule(imol) ||
       is_valid_map_molecule(imol)) {
      std::vector<std::string> symop_list =
         graphics_info_t::molecules[imol].get_symop_strings();
      r = generic_string_vector_to_list_internal(symop_list);
   }
   return r;
}
#endif

// BL says:: python's get_symmetry:
#ifdef USE_PYTHON
PyObject *get_symmetry_py(int imol) {

   PyObject *r = PyList_New(0);
   if (is_valid_model_molecule(imol) ||
       is_valid_map_molecule(imol)) {
      std::vector<std::string> symop_list =
         graphics_info_t::molecules[imol].get_symop_strings();
      r = generic_string_vector_to_list_internal_py(symop_list);
   }
   return r;
}
#endif //PYTHON


/*! \brief Undo symmetry view. Translate back to main molecule from
  this symmetry position.  */
int undo_symmetry_view() {

   int r=0;

   int imol = first_molecule_with_symmetry_displayed();

   if (is_valid_model_molecule(imol)) {

      graphics_info_t g;
      atom_selection_container_t atom_sel = g.molecules[imol].atom_sel;
      mmdb::Manager *mol = atom_sel.mol;
      float symmetry_search_radius = 1;
      coot::Cartesian screen_centre = g.RotationCentre();
      molecule_extents_t mol_extents(atom_sel, symmetry_search_radius);
      std::vector<std::pair<symm_trans_t, Cell_Translation> > boxes =
         mol_extents.which_boxes(screen_centre, atom_sel);
      if (boxes.size() > 0) {
         std::vector<std::pair<clipper::RTop_orth, clipper::Coord_orth> > symm_mat_and_pre_shift_vec;
         for (unsigned int ibox=0; ibox<boxes.size(); ibox++) {
            symm_trans_t st = boxes[ibox].first;
            Cell_Translation pre_shift = boxes[ibox].second;
            mmdb::mat44 my_matt;
            int err = atom_sel.mol->GetTMatrix(my_matt, st.isym(), st.x(), st.y(), st.z());
            if (err == mmdb::SYMOP_Ok) {
               clipper::RTop_orth rtop_symm = coot::util::make_rtop_orth_from(my_matt);
               // and we also need an RTop for the pre-shift
               clipper::Coord_frac pre_shift_cf(pre_shift.us, pre_shift.vs, pre_shift.ws);
               std::pair<clipper::Cell, clipper::Spacegroup> cs = coot::util::get_cell_symm(mol);
               clipper::Coord_orth pre_shift_co = pre_shift_cf.coord_orth(cs.first);
               std::pair<const clipper::RTop_orth, clipper::Coord_orth> p(rtop_symm, pre_shift_co);
               symm_mat_and_pre_shift_vec.push_back(p);
            }
         }
         // so we have a set of matrices and origins shifts, find the
         // one that brings us closest to an atom in imol
         //
         g.unapply_symmetry_to_view(imol, symm_mat_and_pre_shift_vec);
      }
   } else {
      std::cout << "WARNING:: No molecule found that was displaying symmetry"
                << std::endl;
   }
   return r;
}


int first_molecule_with_symmetry_displayed() {

   int imol = -1; // unset
   int n = graphics_n_molecules();
   graphics_info_t g;
   for (int i=0; i<n; i++) {
      if (is_valid_model_molecule(i)) {
         std::pair<std::vector<float>, std::string> cv =
            g.molecules[i].get_cell_and_symm();
         if (cv.first.size() == 6) {
            if (g.molecules[i].show_symmetry) {
               imol = i;
               break;
            }
         }
      }
   }
   return imol;
}


void residue_info_apply_all_checkbutton_toggled() {

}


void apply_residue_info_changes() {
   graphics_info_t g;
   g.apply_residue_info_changes();
   graphics_draw();
}

void do_distance_define() {

   // std::cout << "Click on 2 atoms: " << std::endl;
   graphics_info_t g;
   g.pick_cursor_maybe();
   g.in_distance_define = 1;
   g.pick_pending_flag = 1;

}

void do_angle_define() {

   // std::cout << "Click on 3 atoms: " << std::endl;
   graphics_info_t g;
   g.pick_cursor_maybe();
   g.in_angle_define = 1;
   g.pick_pending_flag = 1;

}

void do_torsion_define() {

   // std::cout << "Click on 4 atoms: " << std::endl;
   graphics_info_t g;
   g.pick_cursor_maybe();
   g.in_torsion_define = 1;
   g.pick_pending_flag = 1;

}

void clear_measure_distances() {
   graphics_info_t g;
   g.clear_measure_distances();
   g.normal_cursor();
   std::string cmd = "clear-simple-distances";
   std::vector<coot::command_arg_t> args;
   add_to_history_typed(cmd, args);
}

void clear_last_measure_distance() {

   graphics_info_t g;
   g.clear_last_measure_distance();
   g.normal_cursor();
   std::string cmd = "clear-last-simple-distance";
   std::vector<coot::command_arg_t> args;
   add_to_history_typed(cmd, args);
}


void clear_residue_info_edit_list() {

   graphics_info_t g;
   g.reset_residue_info_edits();
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

   entry = widget_from_builder("environment_distance_min_entry");
   snprintf(text, 99, "%-5.1f", g.environment_min_distance);
   gtk_editable_set_text(GTK_EDITABLE(entry), text);

   entry = widget_from_builder("environment_distance_max_entry");
   snprintf(text, 99, "%-5.1f" ,g.environment_max_distance);
   gtk_editable_set_text(GTK_EDITABLE(entry), text);
   free(text);

   GtkWidget *check_button = widget_from_builder("environment_distance_checkbutton");

   if (g.environment_show_distances == 1) {
      // we have to (temporarily) set the flag to 0 because the
      // set_active creates an event which causes
      // toggle_environment_show_distances to run (and thus turn off
      // distances if they were allowed to remain here at 1 (on).
      // Strange but true.
      g.environment_show_distances = 0;
      gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(check_button), 1);
      // std::cout << "filling: button is active" << std::endl;
   } else {
      gtk_check_button_set_active(GTK_CHECK_BUTTON(check_button), 0);
      // std::cout << "filling: button is inactive" << std::endl;
   }
   // set the label button
   check_button = widget_from_builder("environment_distance_label_atom_checkbutton");
   if (g.environment_distance_label_atom) {
     gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(check_button), 1);
   } else {
     gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(check_button), 0);
   }
}

// Called when the OK button of the environment distances dialog is clicked
// (just before it is destroyed).
//
void execute_environment_settings(GtkWidget *widget) {

   GtkWidget *entry;
   float val;
   graphics_info_t g;

   entry = widget_from_builder("environment_distance_min_entry");
   const gchar *text = gtk_editable_get_text(GTK_EDITABLE(entry));
   val = atof(text);
   if (val < 0 || val > 1000) {
      g.environment_min_distance = 2.2;
      std::cout <<  "nonsense value for limit using "
                << g.environment_min_distance << std::endl;
   } else {
      g.environment_min_distance = val;
   }

   entry = widget_from_builder("environment_distance_max_entry");
   text = gtk_editable_get_text(GTK_EDITABLE(entry));
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
   label_check_button = widget_from_builder("environment_distance_label_atom_checkbutton");
   if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(label_check_button))) {
      g.environment_distance_label_atom = 1;
   }

   // not sure that this is necessary now that the toggle function is
   // active
   std::pair<int, int> r =  g.get_closest_atom();
   if (r.first >= 0) {
      g.mol_no_for_environment_distances = r.second;
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
         g.mol_no_for_environment_distances = r.second;
         g.update_environment_distances_maybe(r.first, r.second);
      }
   }
   graphics_draw();
}

void set_show_environment_distances_bumps(int state) {
   graphics_info_t g;
   std::pair<int, int> r =  g.get_closest_atom();
   g.environment_distances_show_bumps = state;
   g.update_environment_distances_maybe(r.first, r.second);
   graphics_draw();
}

void set_show_environment_distances_h_bonds(int state) {
   graphics_info_t g;
   std::pair<int, int> r =  g.get_closest_atom();
   g.environment_distances_show_h_bonds = state;
   g.update_environment_distances_maybe(r.first, r.second);
   graphics_draw();
}


int show_environment_distances_state() {
   return graphics_info_t::environment_show_distances;
}

/*! \brief min and max distances for the environment distances */
void set_environment_distances_distance_limits(float min_dist, float max_dist) {

   graphics_info_t::environment_min_distance = min_dist;
   graphics_info_t::environment_max_distance = max_dist;
}

void set_show_environment_distances_as_solid(int state) {
   graphics_info_t::display_environment_graphics_object_as_solid_flag = state;
}

void set_environment_distances_label_atom(int state) {
  graphics_info_t::environment_distance_label_atom = state;
}

double
add_geometry_distance(int imol_1, float x_1, float y_1, float z_1, int imol_2, float x_2, float y_2, float z_2) {

   graphics_info_t g;
   double d = g.add_measure_distance(coot::Cartesian(x_1, y_1, z_1), coot::Cartesian(x_2, y_2, z_2));
   return d;
}

#ifdef USE_GUILE
double
add_atom_geometry_distance_scm(int imol_1, SCM atom_spec_1, int imol_2, SCM atom_spec_2) {

   double d = -1;
   if (is_valid_model_molecule(imol_1)) {
      if (is_valid_model_molecule(imol_2)) {

         graphics_info_t g;
         coot::atom_spec_t spec_1 = atom_spec_from_scm_expression(atom_spec_1);
         coot::atom_spec_t spec_2 = atom_spec_from_scm_expression(atom_spec_2);
         mmdb::Atom *at_1 = g.molecules[imol_1].get_atom(spec_1);
         mmdb::Atom *at_2 = g.molecules[imol_2].get_atom(spec_2);
         if (! at_1) {
            std::cout << "WARNING:: atom not found from spec " << spec_1 << std::endl;
         } else {
            if (! at_2) {
               std::cout << "WARNING:: atom not found from spec " << spec_2 << std::endl;
            } else {
               // happy path
               coot::Cartesian pos_1(at_1->x, at_1->y, at_1->z);
               coot::Cartesian pos_2(at_2->x, at_2->y, at_2->z);
               // d = g.display_geometry_distance(imol_1, pos_1, imol_2, pos_2); 20211006-PE old function name
               d = g.add_measure_distance(pos_1, pos_2);
               std::cout << "Distance: " << spec_1 << " to " << spec_2 << " is " << d << " A" << std::endl;
            }
         }
      }
   }
   return d;
}
#endif

#ifdef USE_PYTHON
double add_atom_geometry_distance_py(int imol_1, PyObject *atom_spec_1, int imol_2, PyObject *atom_spec_2) {

   double d = -1;
   if (is_valid_model_molecule(imol_1)) {
      if (is_valid_model_molecule(imol_2)) {

         graphics_info_t g;
         coot::atom_spec_t spec_1 = atom_spec_from_python_expression(atom_spec_1);
         coot::atom_spec_t spec_2 = atom_spec_from_python_expression(atom_spec_2);
         mmdb::Atom *at_1 = g.molecules[imol_1].get_atom(spec_1);
         mmdb::Atom *at_2 = g.molecules[imol_2].get_atom(spec_2);
         if (! at_1) {
            std::cout << "WARNING:: atom not found from spec " << spec_1 << std::endl;
         } else {
            if (! at_2) {
               std::cout << "WARNING:: atom not found from spec " << spec_2 << std::endl;
            } else {
               // happy path
               coot::Cartesian pos_1(at_1->x, at_1->y, at_1->z);
               coot::Cartesian pos_2(at_2->x, at_2->y, at_2->z);
               // d = g.display_geometry_distance(imol_1, pos_1, imol_2, pos_2); 20211006-PE old function name
               d = g.add_measure_distance(pos_1, pos_2);
               std::cout << "Distance: " << spec_1 << " to " << spec_2 << " is " << d << " A" << std::endl;
            }
         }
      }
   }
   return d;
}
#endif


/*  ----------------------------------------------------------------------- */
/*                  pointer position                                        */
/*  ----------------------------------------------------------------------- */
/* section Pointer Position Function */
/*! \name Pointer Position Function */
/* \{ */
/*! \brief return the [x,y] position of the pointer in fractional coordinates.

    may return false if pointer is not available */
#ifdef USE_PYTHON
PyObject *get_pointer_position_frac_py() {

   PyObject *r = Py_False;

   if (graphics_info_t::use_graphics_interface_flag) {

      graphics_info_t g;
      std::pair<double, double> xy = g.get_pointer_position_frac();
      r = PyList_New(2);
      PyList_SetItem(r, 0, PyFloat_FromDouble(xy.first));
      PyList_SetItem(r, 1, PyFloat_FromDouble(xy.second));

   }
   if (PyBool_Check(r))
     Py_INCREF(r);
   return r;
}
#endif // USE_PYTHON



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
   graphics_info_t g;
   g.reset_residue_info_edits(); // 20230515-PE why?
   std::string cmd = "set-show-pointer-distances";
   std::vector<coot::command_arg_t> args;
   args.push_back(istate);
   add_to_history_typed(cmd, args);
}


/*! \brief show the state of display of the  pointer distances  */
int  show_pointer_distances_state() {

   return graphics_info_t::show_pointer_distances_flag;
}

#ifdef HAVE_GOOCANVAS
#include "goograph/goograph.hh"
#endif

#include "coot-utils/xmap-stats.hh"


void
fill_map_histogram_widget(int imol, GtkWidget *map_contour_frame) {

#ifdef HAVE_GOOCANVAS

   if (is_valid_map_molecule(imol)) {
      // set_and_get_histogram_values(); surely?

      unsigned int n_bins = 1000;
      bool ipz = graphics_info_t::ignore_pseudo_zeros_for_map_stats;
      mean_and_variance<float> mv = graphics_info_t::molecules[imol].set_and_get_histogram_values(n_bins, ipz);

      unsigned int n = mv.size();

      if (n == 1) {
         // pass, previous fail
      } else {

         if (true) {

            float rmsd = sqrt(mv.variance);

            if (false) {
               std::cout << "mv: mean: " << mv.mean << std::endl;
               std::cout << "mv: var: " << mv.variance << std::endl;
               std::cout << "mv: sd: " << sqrt(mv.variance) << std::endl;
            }

            if (mv.bins.size() > 0) {
               std::vector<std::pair<double, double> > data(mv.bins.size());
               for (unsigned int ibin=0; ibin<mv.bins.size(); ibin++) {
                  double x = (ibin+0.5)*mv.bin_width + mv.min_density;
                  double y = mv.bins[ibin];
                  data[ibin] = std::pair<double, double> (x, y);
               }

               int graph_x_n_pixels = 300;
               int graph_y_n_pixels =  64;
               coot::goograph* g = new coot::goograph(graph_x_n_pixels, graph_y_n_pixels);
               int trace = g->trace_new();

               g->set_plot_title("");
               g->set_data(trace, data);
               // g->set_axis_label(coot::goograph::X_AXIS, "Density Value");
               // g->set_axis_label(coot::goograph::Y_AXIS, "Counts");
               g->set_trace_type(trace, coot::graph_trace_info_t::PLOT_TYPE_BAR);
               g->set_trace_colour(trace, "#111111");
               if (true) {
                  if (data.size() == 0) {
                  } else {
                     // find y_max ignoring the peak
                     double y_max           = -1e100;
                     double y_max_secondary = -1e100;
                     unsigned int idata_peak = 0;
                     for (unsigned int idata=0; idata<data.size(); idata++) {
                        if (data[idata].second > y_max) {
                           y_max = data[idata].second;
                           idata_peak = idata;
                        }
                     }
                     for (unsigned int idata=0; idata<data.size(); idata++) {
                        if (idata != idata_peak)
                           if (data[idata].second > y_max_secondary)
                              y_max_secondary = data[idata].second;
                     }

                     // std::cout << ":::::::::: y_max_secondary " << y_max_secondary << std::endl;

                     g->set_extents(coot::goograph::X_AXIS,
                                    mv.mean-2*sqrt(mv.variance),
                                    mv.mean+2*sqrt(mv.variance)
                                    );
                     // the bar width depends on the X extents (for aesthetics)
                     double contour_level_bar_width = sqrt(mv.variance) * 0.1;
                     if (false)
                        std::cout << "::::: set_extents() X: "
                                  << mv.mean-2*sqrt(mv.variance) << " "
                                  << mv.mean+2*sqrt(mv.variance) << "\n";
                     if (y_max_secondary > 0) {
                        double y_max_graph = y_max_secondary * 1.3;
                        g->set_extents(coot::goograph::Y_AXIS,
                                       0,
                                       0.7 * y_max_graph); // calls set_data_scales() internall
                        if (false)
                           std::cout << "::::: set_extents() Y: "
                                     << 0 << " " << y_max_graph << std::endl;
                        // draw x-axis ticks only
                        g->set_draw_axis(coot::goograph::Y_AXIS, false);
                        g->set_draw_axis(coot::goograph::X_AXIS, false);
                        g->set_draw_ticks(coot::goograph::Y_AXIS, false);

                        // draw the contour level bar
                        float cl = graphics_info_t::molecules[imol].get_contour_level();

                        // std::pair<GdkRGBA, GdkRGBA> map_colours() const;
                        std::pair<GdkRGBA, GdkRGBA> map_colours =
                           graphics_info_t::molecules[imol].get_map_colours();

#if 0
                        if (true) { // this test is needed?
                           coot::colour_holder ch(map_colours);
                           void (*func)(int, float) = set_contour_level_absolute;
                           GtkWidget *canvas = g->get_canvas();
                           // moving the contour level bar redraws the graph and calls func
                           g->add_contour_level_box(cl, "#111111", 1.4, ch.hex(), imol, rmsd, func);
                           // ticks fall off the graph, sigh - add an offset
                           gtk_widget_set_size_request(canvas, graph_x_n_pixels, graph_y_n_pixels+10);
                           g->draw_graph();

                           gtk_widget_set_visible(canvas, TRUE);
                           gtk_container_add(GTK_CONTAINER(map_contour_frame), canvas);
                        }
#endif
                     }
                  }
               }
            }
         }
      }
   }
#else
   gtk_widget_set_visible(map_contour_frame, FALSE);
#endif
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
int
get_dynamic_map_size_display(){
   int ret = graphics_info_t::dynamic_map_size_display;
   return ret;
}
void
set_dynamic_map_sampling_on() {
   graphics_info_t::dynamic_map_resampling = 1;
}
void
set_dynamic_map_sampling_off() {
   graphics_info_t::dynamic_map_resampling = 0;
}
int
get_dynamic_map_sampling(){
   int ret = graphics_info_t::dynamic_map_resampling;
   return ret;
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

   if (g.console_display_commands.display_commands_flag) {

      char esc = 27;
      // std::string esc = "esc";
      if (g.console_display_commands.hilight_flag) {
        // std::cout << esc << "[34m";
#ifdef WINDOWS_MINGW
        // use the console cursor infot to distinguish between DOS and MSYS
        // shell
        CONSOLE_CURSOR_INFO ConCurInfo;
        if (GetConsoleCursorInfo(GetStdHandle(STD_OUTPUT_HANDLE), &ConCurInfo)) {
          // we have a DOS shell
          SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
                                  FOREGROUND_RED |
                                  FOREGROUND_GREEN |
                                  FOREGROUND_BLUE |
                                  FOREGROUND_INTENSITY);
        } else {
          // we have MSYS (or whatever else shell)
          std::cout << esc << "[1m";
        }
#else
         std::cout << esc << "[1m";
#endif // MINGW
      } else {
         std::cout << "INFO:: Command: ";
      }

      // Make it colourful?
      if (g.console_display_commands.hilight_colour_flag) {
#ifdef WINDOWS_MINGW
        CONSOLE_CURSOR_INFO ConCurInfo;
        if (GetConsoleCursorInfo(GetStdHandle(STD_OUTPUT_HANDLE), &ConCurInfo)) {
          // we have a DOS shell
          switch (g.console_display_commands.colour_prefix) {
          case(1):
            // red
            SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
                                    FOREGROUND_RED);
            break;
          case(2):
            // green
            SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
                                    FOREGROUND_GREEN);
            break;
          case(3):
            // yellow
            SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
                                    FOREGROUND_RED |
                                    FOREGROUND_GREEN);
            break;
          case(4):
            // blue
            SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
                                    FOREGROUND_BLUE);
            break;
          case(5):
            // magenta
            SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
                                    FOREGROUND_RED |
                                    FOREGROUND_BLUE);
            break;
          case(6):
            // cyan
            SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
                                    FOREGROUND_GREEN |
                                    FOREGROUND_BLUE);
            break;
          default:
            //white
            SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
                                    FOREGROUND_RED |
                                    FOREGROUND_GREEN |
                                    FOREGROUND_BLUE);
          }
        } else {
          // MSYS shell
         std::cout << esc << "[3"
                   << g.console_display_commands.colour_prefix << "m";
        }
#else
         std::cout << esc << "[3"
                   << g.console_display_commands.colour_prefix << "m";
#endif // MINGW
      }

#if defined USE_GUILE && !defined WINDOWS_MINGW
      std::cout << graphics_info_t::schemize_command_strings(command_strings);
#else
#ifdef USE_PYTHON
      std::cout << graphics_info_t::pythonize_command_strings(command_strings);
#endif // USE_PYTHON
#endif // USE_GUILE/MINGW

      if (g.console_display_commands.hilight_flag) {// hilight off
#ifdef WINDOWS_MINGW
        CONSOLE_CURSOR_INFO ConCurInfo;
        if (GetConsoleCursorInfo(GetStdHandle(STD_OUTPUT_HANDLE), &ConCurInfo)) {
          // we have a DOS shell (reset to white)
          SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),
                                  FOREGROUND_RED |
                                  FOREGROUND_GREEN |
                                  FOREGROUND_BLUE);
        } else {
          // MSYS shell
         std::cout << esc << "[0m"; // reset
        }
#else
         std::cout << esc << "[0m"; // reset
         //std::cout << esc; // reset
#endif // MINGW
      }
      std::cout << std::endl;
   }

#ifdef USE_MYSQL_DATABASE

   add_to_database(command_strings);
#endif
}

void add_to_history_typed(const std::string &command,
                          const std::vector<coot::command_arg_t> &args) {

   std::vector<std::string> command_strings;

   command_strings.push_back(command);
   for (unsigned int i=0; i<args.size(); i++)
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

   if (s == "run-refmac-by-filename") {
      ss = "refmac.run_refmac_by_filename";
   } else  {

      // 20231208-PE
      // std::cout << "************************* namespace this python command? " << ss << std::endl;
      // ss = "coot." + ss; // hideous hack. presumes that all state script commands are in coot module
                            // (most of them are, of course).
                            // 20231208-PE Are they though? If the function is in the coot namespace,
                            // why not just call the function directly?
   }
   return ss;
}

std::string schemize_command_name(const std::string &s_in) {

   std::string ss;
   std::string s(s_in);
   if (s.length() > 5)
      if (s.substr(0,5) == ".coot")
         s.erase(0,5);
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

   graphics_info_t::console_display_commands.display_commands_flag = istate;
}

void set_console_display_commands_hilights(short int bold_flag, short int colour_flag, int colour_index) {

   graphics_info_t g;
   g.console_display_commands.hilight_flag = bold_flag;
   g.console_display_commands.hilight_colour_flag = colour_flag;
   g.console_display_commands.colour_prefix = colour_index;
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
   args.push_back(single_quote(filename));
   add_to_history_typed(cmd, args);
}

/*! \brief save the current state to file filename */
void save_state_file_py(const char *filename) {
   graphics_info_t g;
   g.save_state_file(std::string(filename), coot::PYTHON_SCRIPT);
   std::string cmd = "save-state-file";
   std::vector<coot::command_arg_t> args;
   args.push_back(single_quote(filename));
   add_to_history_typed(cmd, args);
}



#ifdef USE_GUILE
/* Return the default file name suggestion (that would come up in the
   save coordinates dialog) or scheme faalse if imol is not a valid
   model molecule. */
SCM save_coords_name_suggestion_scm(int imol) {

   SCM r = SCM_BOOL_F;

   if (is_valid_model_molecule(imol)) {
      std::string s = graphics_info_t::molecules[imol].stripped_save_name_suggestion();
      r = scm_from_locale_string(s.c_str());
   }
   return r;
}
#endif /*  USE_GUILE */


#ifdef USE_PYTHON
PyObject *save_coords_name_suggestion_py(int imol) {

   PyObject *r = Py_False;
   if (is_valid_model_molecule(imol)) {
      std::string s = graphics_info_t::molecules[imol].stripped_save_name_suggestion();
      r = myPyString_FromString(s.c_str());
   }

   return r;
}
#endif /*  USE_PYTHON */



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

/*! \brief return the resolution set in the header of the
  model/coordinates file.  If this number is not available, return a
  number less than 0.  */
float model_resolution(int imol) {

   float r = -1;
   if (is_valid_model_molecule(imol)) {
      r = graphics_info_t::molecules[imol].atom_sel.mol->GetResolution();
   }
   return r;
}



/*  ------------------------------------------------------------------------ */
/*                     residue exists?                                       */
/*  ------------------------------------------------------------------------ */

int does_residue_exist_p(int imol, const char *chain_id, int resno, const char *inscode) {

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
// Caller should dispose of returned pointer.
const char *mtz_hklin_for_map(int imol_map) {

   std::string mtz;

   if (is_valid_map_molecule(imol_map)) {
      mtz = graphics_info_t::molecules[imol_map].save_mtz_file_name;
   }
   std::string cmd = "mtz-hklin-for-map";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol_map);
   add_to_history_typed(cmd, args);
   const char *s = strdup(mtz.c_str());
   return s;
}

/*! \brief return the FP column in the file that was use to generate
  the map

  return 0 when there is no mtz file associated with that map (it was
  generated from a CCP4 map file say). */
// Caller should dispose of returned pointer.
const char *mtz_fp_for_map(int imol_map) {

   std::string fp;
   if (is_valid_map_molecule(imol_map)) {
      fp = graphics_info_t::molecules[imol_map].save_f_col;
   }
   std::string cmd = "mtz-fp-for-map";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol_map);
   add_to_history_typed(cmd, args);
   const char *s = strdup(fp.c_str());
   return s;
}

/*! \brief return the phases column in mtz file that was use to generate
  the map

  return 0 when there is no mtz file associated with that map (it was
  generated from a CCP4 map file say). */
// Caller should dispose of returned pointer.
const char *mtz_phi_for_map(int imol_map) {

   std::string phi;
   if (is_valid_map_molecule(imol_map)) {
      phi = graphics_info_t::molecules[imol_map].save_phi_col;
   }
   std::string cmd = "mtz-phi-for-map";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol_map);
   add_to_history_typed(cmd, args);
   const char *s = strdup(phi.c_str());
   return s;

}

/*! \brief return the weight column in the mtz file that was use to
  generate the map

  return 0 when there is no mtz file associated with that map (it was
  generated from a CCP4 map file say) or no weights were used. */
// Caller should dispose of returned pointer.
const char *mtz_weight_for_map(int imol_map) {

   std::string weight;
   if (is_valid_map_molecule(imol_map)) {
      weight = graphics_info_t::molecules[imol_map].save_weight_col;
   }
   std::string cmd = "mtz-weight-for-map";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol_map);
   add_to_history_typed(cmd, args);
   const char *s = strdup(weight.c_str());
   return s;
}

/*! \brief return flag for whether weights were used that was use to
  generate the map

  return 0 when no weights were used or there is no mtz file
  associated with that map. */
short int mtz_use_weight_for_map(int imol_map) {

   short int i = 0;
   if (is_valid_map_molecule(imol_map)) {
      i = graphics_info_t::molecules[imol_map].save_use_weights;
   }
   std::string cmd = "mtz-use-weight-for-map";
   std::vector<coot::command_arg_t> args;
   args.push_back(imol_map);
   add_to_history_typed(cmd, args);
   return i;
}

// return #f or ("xxx.mtz" "FPH" "PHWT" "" #f)
//
#ifdef USE_GUILE
SCM map_parameters_scm(int imol) {

   SCM r = SCM_BOOL_F;
   if (is_valid_map_molecule(imol)) {
      r = SCM_EOL;
      if (graphics_info_t::molecules[imol].save_use_weights)
         r = scm_cons(SCM_BOOL_T, r);
      else
         r = scm_cons(SCM_BOOL_F, r);
      r = scm_cons(scm_from_locale_string(graphics_info_t::molecules[imol].save_weight_col.c_str()), r);
      r = scm_cons(scm_from_locale_string(graphics_info_t::molecules[imol].save_phi_col.c_str()), r);
      r = scm_cons(scm_from_locale_string(graphics_info_t::molecules[imol].save_f_col.c_str()), r);
      r = scm_cons(scm_from_locale_string(graphics_info_t::molecules[imol].save_mtz_file_name.c_str()), r);
   }
   return r;
}
#endif // USE_GUILE

// return False or ["xxx.mtz", "FPH", "PHWT", "", False]
//
#ifdef USE_PYTHON
PyObject *map_parameters_py(int imol) {

   PyObject *r = Py_False;
   if (is_valid_map_molecule(imol)) {
      r = PyList_New(5);
      PyList_SetItem(r, 0, myPyString_FromString(graphics_info_t::molecules[imol].save_mtz_file_name.c_str()));
      PyList_SetItem(r, 1, myPyString_FromString(graphics_info_t::molecules[imol].save_f_col.c_str()));
      PyList_SetItem(r, 2, myPyString_FromString(graphics_info_t::molecules[imol].save_phi_col.c_str()));
      PyList_SetItem(r, 3, myPyString_FromString(graphics_info_t::molecules[imol].save_weight_col.c_str()));
      if (graphics_info_t::molecules[imol].save_use_weights) {
         // Py_INCREF(Py_True);
         PyObject *o_py = PyBool_FromLong(true);
         PyList_SetItem(r, 4, o_py);
      } else {
         // Py_INCREF(Py_False);
         PyObject *o_py = PyBool_FromLong(false);
         PyList_SetItem(r, 4, o_py);
      }
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif // USE_PYTHON

#ifdef USE_GUILE
// return #f or (45 46 47 90 90 120), angles in degrees.
SCM cell_scm(int imol) {
   SCM r = SCM_BOOL_F;
   if (is_valid_map_molecule(imol) || (is_valid_model_molecule(imol))) {
      std::pair<bool, clipper::Cell> cell = graphics_info_t::molecules[imol].cell();
      if (cell.first) {
         r = SCM_EOL;
         r = scm_cons(scm_from_double(clipper::Util::rad2d(cell.second.descr().gamma())), r);
         r = scm_cons(scm_from_double(clipper::Util::rad2d(cell.second.descr().beta() )), r);
         r = scm_cons(scm_from_double(clipper::Util::rad2d(cell.second.descr().alpha())), r);
         r = scm_cons(scm_from_double(cell.second.descr().c()), r);
         r = scm_cons(scm_from_double(cell.second.descr().b()), r);
         r = scm_cons(scm_from_double(cell.second.descr().a()), r);
      }
   }
   return r;
}
#endif // USE_GUILE



#ifdef USE_PYTHON
// return False or [45, 46, 47, 90, 90, 120), angles in degrees.
PyObject *cell_py(int imol) {
   PyObject *r = Py_False;
   if (is_valid_map_molecule(imol) || (is_valid_model_molecule(imol))) {
      std::pair<bool, clipper::Cell> cell = graphics_info_t::molecules[imol].cell();
      if (cell.first) {
         r = PyList_New(6);
         PyList_SetItem(r, 0, PyFloat_FromDouble(cell.second.descr().a()));
         PyList_SetItem(r, 1, PyFloat_FromDouble(cell.second.descr().b()));
         PyList_SetItem(r, 2, PyFloat_FromDouble(cell.second.descr().c()));
         PyList_SetItem(r, 3, PyFloat_FromDouble(clipper::Util::rad2d(cell.second.descr().alpha())));
         PyList_SetItem(r, 4, PyFloat_FromDouble(clipper::Util::rad2d(cell.second.descr().beta())));
         PyList_SetItem(r, 5, PyFloat_FromDouble(clipper::Util::rad2d(cell.second.descr().gamma())));
      }
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif // USE_PYTHON



/*! \brief Put text at x,y,z  */
// use should be given access to colour and size.
int place_text(const char *text, float x, float y, float z, int size) {

   int handle = graphics_info_t::generic_texts.size();
   std::string s(text);
   coot::generic_text_object_t o(s, handle, x, y, z);
   graphics_info_t::generic_texts.push_back(o);

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
   graphics_draw();

   return handle; // same value as above.
}

void remove_text(int text_handle) {

   std::vector<coot::generic_text_object_t>::iterator it;
   for (it = graphics_info_t::generic_texts.begin();
        it != graphics_info_t::generic_texts.end();
        it++) {
      if (it->handle == text_handle) {
         graphics_info_t::generic_texts.erase(it);
         break;
      }
   }
   std::string cmd = "remove-text";
   std::vector<coot::command_arg_t> args;
   args.push_back(text_handle);
   add_to_history_typed(cmd, args);
   graphics_draw();
}


void edit_text(int text_handle, const char *str) {

   graphics_info_t g;
   if (str) {
      if (text_handle >= 0) {
         unsigned int ui_text_handle = text_handle;
         if (ui_text_handle < g.generic_texts.size()) {
            g.generic_texts[ui_text_handle].s = str;
         }
      }
   }
   std::string cmd = "edit-text";
   std::vector<coot::command_arg_t> args;
   args.push_back(text_handle);
   args.push_back(str);
   add_to_history_typed(cmd, args);
   graphics_draw();
}

/*! \brief return the closest text that is with r A of the given
  position.  If no text item is close, then return -1 */
int text_index_near_position(float x, float y, float z, float rad) {

   int r = -1;
   graphics_info_t g;
   double best_dist = 999999999.9; // not (long) integer, conversion to double problems in GCC 4.1.2

   std::cout << "size: " << g.generic_texts.size() << std::endl;

   for (unsigned int i=0; i<g.generic_texts.size(); i++) {
      std::cout << "i " << i << std::endl;
      clipper::Coord_orth p1(x,y,z);
      clipper::Coord_orth p2(g.generic_texts[i].x,
                             g.generic_texts[i].y,
                             g.generic_texts[i].z);
      double d = (p1-p2).lengthsq();
      std::cout << "   d " << d  << std::endl;
      if (d < rad*rad) {
         if (d < best_dist) {
            best_dist = d;
            r = i;
         }
      }
   }
   return r;
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

// BL says:: python's func.
#ifdef USE_PYTHON
PyObject *dictionaries_read_py() {

   return generic_string_vector_to_list_internal_py(*graphics_info_t::cif_dictionary_filename_vec);
}
#endif // PYTHON

std::vector<std::string>
dictionary_entries() {
   graphics_info_t g;
   return g.Geom_p()->monomer_restraints_comp_ids();
}

void
debug_dictionary() {
   graphics_info_t g;
   g.Geom_p()->debug();
}

#ifdef USE_GUILE
SCM dictionary_entries_scm() {
   std::vector<std::string> comp_ids = dictionary_entries();
   return generic_string_vector_to_list_internal(comp_ids);
}
#endif // USE_GUILE

#ifdef USE_PYTHON
PyObject *dictionary_entries_py() {

   std::vector<std::string> comp_ids = dictionary_entries();
   return generic_string_vector_to_list_internal_py(comp_ids);
}
#endif // USE_PYTHON





#ifdef USE_GUILE
SCM cif_file_for_comp_id_scm(const std::string &comp_id) {

   graphics_info_t g;
   int imol = 0; // dummy. Hmm.
   std::string f = g.Geom_p()->get_cif_file_name(comp_id, imol);
   return scm_from_locale_string(f.c_str());
}
#endif // GUILE


#ifdef USE_PYTHON
PyObject *cif_file_for_comp_id_py(const std::string &comp_id) {

   int imol = 0; // dummy. Hmm.
   graphics_info_t g;
   return myPyString_FromString(g.Geom_p()->get_cif_file_name(comp_id, imol).c_str());
}
#endif // PYTHON

std::string SMILES_for_comp_id(const std::string &comp_id) {

   int imol_enc = coot::protein_geometry::IMOL_ENC_ANY; // pass this?

   graphics_info_t g;
   std::string s;
   try {
      s = g.Geom_p()->Get_SMILES_for_comp_id(comp_id, imol_enc); // can throw
   }
   catch (const std::runtime_error &e) {
      std::cout << "WARNING::" << e.what() << std::endl;
   }
   catch (...) {
      std::cout << "SMILES_for_comp_id() caught generic throw" << std::endl;
   }
   return s;
}

#ifdef USE_GUILE
SCM SMILES_for_comp_id_scm(const std::string &comp_id) {

   SCM r = SCM_BOOL_F;
   try {
      std::string s = SMILES_for_comp_id(comp_id);
      r = scm_from_locale_string(s.c_str());
   }
   catch (const std::runtime_error &rte) {
      std::cout << "WARNING:: " << rte.what() << std::endl;
   }
   return r;
}
#endif


#ifdef USE_PYTHON
PyObject *SMILES_for_comp_id_py(const std::string &comp_id) {
   PyObject *r = Py_False;
   try {
      std::string s = SMILES_for_comp_id(comp_id);
      r = myPyString_FromString(s.c_str());
   }
   catch (const std::runtime_error &rte) {
      std::cout << "WARNING:: " << rte.what() << std::endl;
   }
   if (PyBool_Check(r))
      Py_INCREF(r);
   return r;
}
#endif // USE_PYTHON


/*  ----------------------------------------------------------------------- */
/*                         Restraints                                       */
/*  ----------------------------------------------------------------------- */
#ifdef USE_GUILE
SCM monomer_restraints(const char *monomer_type) {

   int imol = 0; // dummy.  This should be passed?

   SCM r = SCM_BOOL_F;

   graphics_info_t g;
   // this forces try_dynamic_add()
   g.Geom_p()->have_dictionary_for_residue_type(monomer_type, imol, ++g.cif_dictionary_read_number);
   // this doesn't force try_dynamic_add().  Hmmm.. FIXME (or think about it).
   std::pair<short int, coot::dictionary_residue_restraints_t> p =
      g.Geom_p()->get_monomer_restraints(monomer_type, imol);
   if (p.first) {

      coot::dictionary_residue_restraints_t restraints = p.second;
      r = SCM_EOL;

      // ------------------ chem_comp -------------------------
      coot::dict_chem_comp_t info = restraints.residue_info;
      SCM chem_comp_scm = SCM_EOL;
      chem_comp_scm = scm_cons(scm_from_locale_string(info.comp_id.c_str()),           chem_comp_scm);
      chem_comp_scm = scm_cons(scm_from_locale_string(info.three_letter_code.c_str()), chem_comp_scm);
      chem_comp_scm = scm_cons(scm_from_locale_string(info.name.c_str()),              chem_comp_scm);
      chem_comp_scm = scm_cons(scm_from_locale_string(info.group.c_str()),             chem_comp_scm);
      chem_comp_scm = scm_cons(scm_from_int(info.number_atoms_all),              chem_comp_scm);
      chem_comp_scm = scm_cons(scm_from_int(info.number_atoms_nh),               chem_comp_scm);
      chem_comp_scm = scm_cons(scm_from_locale_string(info.description_level.c_str()), chem_comp_scm);
      chem_comp_scm = scm_reverse(chem_comp_scm);
      SCM chem_comp_container = SCM_EOL;
      // chem_comp_container = scm_cons(chem_comp_scm, chem_comp_container);
      chem_comp_container = scm_cons(scm_from_locale_string("_chem_comp"), chem_comp_scm);

      // ------------------ chem_comp_atom -------------------------
      std::vector<coot::dict_atom> atom_info = restraints.atom_info;
      int n_atoms = atom_info.size();
      SCM atom_info_list = SCM_EOL;
      for (int iat=0; iat<n_atoms; iat++) {
         SCM atom_attributes_list = SCM_EOL;
         atom_attributes_list = scm_cons(scm_from_locale_string(atom_info[iat].atom_id_4c.c_str()),   atom_attributes_list);
         atom_attributes_list = scm_cons(scm_from_locale_string(atom_info[iat].type_symbol.c_str()),  atom_attributes_list);
         atom_attributes_list = scm_cons(scm_from_locale_string(atom_info[iat].type_energy.c_str()),  atom_attributes_list);
         atom_attributes_list = scm_cons(scm_from_double(atom_info[iat].partial_charge.second), atom_attributes_list);
         SCM partial_flag = SCM_BOOL_F;
         if (atom_info[iat].partial_charge.first)
            partial_flag = SCM_BOOL_T;
         atom_attributes_list = scm_cons(partial_flag, atom_attributes_list);
         atom_attributes_list = scm_reverse(atom_attributes_list);
         atom_info_list = scm_cons(atom_attributes_list, atom_info_list);
      }
      atom_info_list = scm_reverse(atom_info_list);
      SCM atom_info_list_container = SCM_EOL;
      // atom_info_list_container = scm_cons(atom_info_list, atom_info_list_container);
      atom_info_list_container = scm_cons(scm_from_locale_string("_chem_comp_atom"), atom_info_list);


      // ------------------ Bonds -------------------------
      SCM bond_restraint_list = SCM_EOL;

      for (unsigned int ibond=0; ibond<restraints.bond_restraint.size(); ibond++) {
         const coot::dict_bond_restraint_t &bond_restraint = restraints.bond_restraint[ibond];
         std::string a1 = bond_restraint.atom_id_1_4c();
         std::string a2 = bond_restraint.atom_id_2_4c();
         std::string type = bond_restraint.type();
         SCM bond_restraint_scm = SCM_EOL;
         SCM esd_scm = SCM_BOOL_F;
         SCM d_scm   = SCM_BOOL_F;
         try {
            double esd = bond_restraint.value_esd();
            double d   = bond_restraint.value_dist();
            esd_scm = scm_from_double(esd);
            d_scm   = scm_from_double(d);
         }
         catch (const std::runtime_error &rte) {
            // we use the default values of #f, if the esd or dist is not set.
         }
         bond_restraint_scm = scm_cons(esd_scm, bond_restraint_scm);
         bond_restraint_scm = scm_cons(d_scm,   bond_restraint_scm);
         bond_restraint_scm = scm_cons(scm_from_locale_string(type.c_str()), bond_restraint_scm);
         bond_restraint_scm = scm_cons(scm_from_locale_string(a2.c_str()),   bond_restraint_scm);
         bond_restraint_scm = scm_cons(scm_from_locale_string(a1.c_str()),   bond_restraint_scm);
         bond_restraint_list = scm_cons(bond_restraint_scm, bond_restraint_list);
      }
      SCM bond_restraints_container = SCM_EOL;
      // bond_restraints_container = scm_cons(bond_restraint_list, bond_restraints_container);
      bond_restraints_container = scm_cons(scm_from_locale_string("_chem_comp_bond"), bond_restraint_list);

      // ------------------ Angles -------------------------
      SCM angle_restraint_list = SCM_EOL;
      for (unsigned int iangle=0; iangle<restraints.angle_restraint.size(); iangle++) {
         coot::dict_angle_restraint_t angle_restraint = restraints.angle_restraint[iangle];
         std::string a1 = angle_restraint.atom_id_1_4c();
         std::string a2 = angle_restraint.atom_id_2_4c();
         std::string a3 = angle_restraint.atom_id_3_4c();
         double d   = angle_restraint.angle();
         double esd = angle_restraint.esd();
         SCM angle_restraint_scm = SCM_EOL;
         angle_restraint_scm = scm_cons(scm_from_double(esd), angle_restraint_scm);
         angle_restraint_scm = scm_cons(scm_from_double(d),   angle_restraint_scm);
         angle_restraint_scm = scm_cons(scm_from_locale_string(a3.c_str()),   angle_restraint_scm);
         angle_restraint_scm = scm_cons(scm_from_locale_string(a2.c_str()),   angle_restraint_scm);
         angle_restraint_scm = scm_cons(scm_from_locale_string(a1.c_str()),   angle_restraint_scm);
         angle_restraint_list = scm_cons(angle_restraint_scm, angle_restraint_list);
      }
      SCM angle_restraints_container = SCM_EOL;
      // angle_restraints_container = scm_cons(angle_restraint_list, angle_restraints_container);
      angle_restraints_container = scm_cons(scm_from_locale_string("_chem_comp_angle"), angle_restraint_list);

      // ------------------ Torsions -------------------------
      SCM torsion_restraint_list = SCM_EOL;
      for (unsigned int itorsion=0; itorsion<restraints.torsion_restraint.size(); itorsion++) {
         coot::dict_torsion_restraint_t torsion_restraint = restraints.torsion_restraint[itorsion];
         std::string id = torsion_restraint.id();
         std::string a1 = torsion_restraint.atom_id_1_4c();
         std::string a2 = torsion_restraint.atom_id_2_4c();
         std::string a3 = torsion_restraint.atom_id_3_4c();
         std::string a4 = torsion_restraint.atom_id_4_4c();
         double tor  = torsion_restraint.angle();
         double esd = torsion_restraint.esd();
         int period = torsion_restraint.periodicity();
         SCM torsion_restraint_scm = SCM_EOL;
         torsion_restraint_scm = scm_cons(scm_from_int(period), torsion_restraint_scm);
         torsion_restraint_scm = scm_cons(scm_from_double(esd), torsion_restraint_scm);
         torsion_restraint_scm = scm_cons(scm_from_double(tor), torsion_restraint_scm);
         torsion_restraint_scm = scm_cons(scm_from_locale_string(a4.c_str()),   torsion_restraint_scm);
         torsion_restraint_scm = scm_cons(scm_from_locale_string(a3.c_str()),   torsion_restraint_scm);
         torsion_restraint_scm = scm_cons(scm_from_locale_string(a2.c_str()),   torsion_restraint_scm);
         torsion_restraint_scm = scm_cons(scm_from_locale_string(a1.c_str()),   torsion_restraint_scm);
         torsion_restraint_scm = scm_cons(scm_from_locale_string(id.c_str()),   torsion_restraint_scm);
         torsion_restraint_list = scm_cons(torsion_restraint_scm, torsion_restraint_list);
      }
      SCM torsion_restraints_container = SCM_EOL;
      // torsion_restraints_container = scm_cons(torsion_restraint_list, torsion_restraints_container);
      torsion_restraints_container = scm_cons(scm_from_locale_string("_chem_comp_tor"), torsion_restraint_list);


      // ------------------ Planes -------------------------
      SCM plane_restraint_list = SCM_EOL;
      for (unsigned int iplane=0; iplane<restraints.plane_restraint.size(); iplane++) {
         coot::dict_plane_restraint_t plane_restraint = restraints.plane_restraint[iplane];
         SCM atom_list = SCM_EOL;
         for (int iat=0; iat<plane_restraint.n_atoms(); iat++) {

            std::string at = plane_restraint[iat].first;
            atom_list = scm_cons(scm_from_locale_string(at.c_str()), atom_list);
         }
         atom_list = scm_reverse(atom_list);

         double esd = plane_restraint.dist_esd(0); // fixme
         SCM plane_id_scm = scm_from_locale_string(plane_restraint.plane_id.c_str());

         SCM plane_restraint_scm = SCM_EOL;
         plane_restraint_scm = scm_cons(scm_from_double(esd), plane_restraint_scm);
         plane_restraint_scm = scm_cons(atom_list, plane_restraint_scm);
         plane_restraint_scm = scm_cons(plane_id_scm, plane_restraint_scm);
         plane_restraint_list = scm_cons(plane_restraint_scm, plane_restraint_list);
      }
      SCM plane_restraints_container = SCM_EOL;
      // plane_restraints_container = scm_cons(plane_restraint_list, plane_restraints_container);
      plane_restraints_container = scm_cons(scm_from_locale_string("_chem_comp_plane_atom"),
                                            plane_restraint_list);


      // ------------------ Chirals -------------------------
      SCM chiral_restraint_list = SCM_EOL;
      for (unsigned int ichiral=0; ichiral<restraints.chiral_restraint.size(); ichiral++) {
         coot::dict_chiral_restraint_t chiral_restraint = restraints.chiral_restraint[ichiral];

         std::string a1 = chiral_restraint.atom_id_1_4c();
         std::string a2 = chiral_restraint.atom_id_2_4c();
         std::string a3 = chiral_restraint.atom_id_3_4c();
         std::string ac = chiral_restraint.atom_id_c_4c();
         std::string chiral_id = chiral_restraint.Chiral_Id();
         int vol_sign = chiral_restraint.volume_sign;

         double esd = chiral_restraint.volume_sigma();
         // int volume_sign = chiral_restraint.volume_sign;
         SCM chiral_restraint_scm = SCM_EOL;
         chiral_restraint_scm = scm_cons(scm_from_double(esd), chiral_restraint_scm);
         chiral_restraint_scm = scm_cons(scm_from_int(vol_sign), chiral_restraint_scm);
         chiral_restraint_scm = scm_cons(scm_from_locale_string(a3.c_str()), chiral_restraint_scm);
         chiral_restraint_scm = scm_cons(scm_from_locale_string(a2.c_str()), chiral_restraint_scm);
         chiral_restraint_scm = scm_cons(scm_from_locale_string(a1.c_str()), chiral_restraint_scm);
         chiral_restraint_scm = scm_cons(scm_from_locale_string(ac.c_str()), chiral_restraint_scm);
         chiral_restraint_scm = scm_cons(scm_from_locale_string(chiral_id.c_str()), chiral_restraint_scm);
         chiral_restraint_list = scm_cons(chiral_restraint_scm, chiral_restraint_list);
      }
      SCM chiral_restraints_container = SCM_EOL;
      // chiral_restraints_container = scm_cons(chiral_restraint_list, chiral_restraints_container);
      chiral_restraints_container = scm_cons(scm_from_locale_string("_chem_comp_chir"), chiral_restraint_list);


      r = scm_cons( chiral_restraints_container, r);
      r = scm_cons(  plane_restraints_container, r);
      r = scm_cons(torsion_restraints_container, r);
      r = scm_cons(  angle_restraints_container, r);
      r = scm_cons(   bond_restraints_container, r);
      r = scm_cons(    atom_info_list_container, r);
      r = scm_cons(         chem_comp_container, r);

   }
   return r;
}
#endif // USE_GUILE

#ifdef USE_PYTHON
PyObject *monomer_restraints_py(std::string monomer_type) {
   return monomer_restraints_for_molecule_py(monomer_type, coot::protein_geometry::IMOL_ENC_ANY);
}
#endif // USE_PYTHON

#ifdef USE_PYTHON
PyObject *monomer_restraints_for_molecule_py(std::string monomer_type, int imol) {

   PyObject *r;
   r = Py_False;

   graphics_info_t g;
   // this forces try_dynamic_add()
   g.Geom_p()->have_dictionary_for_residue_type(monomer_type, imol, ++g.cif_dictionary_read_number);
   // this doesn't force try_dynamic_add().  Hmmm.. FIXME (or think about it).
   std::pair<short int, coot::dictionary_residue_restraints_t> p =
      g.Geom_p()->get_monomer_restraints(monomer_type, imol);
   if (!p.first) {
      std::cout << "WARNING:: can't find " << monomer_type << " in monomer dictionary"
                << std::endl;
   } else {

      r = PyDict_New();

      coot::dictionary_residue_restraints_t restraints = p.second;

      // ------------------ chem_comp -------------------------
      coot::dict_chem_comp_t info = restraints.residue_info;

      PyObject *chem_comp_py = PyList_New(7);
      PyList_SetItem(chem_comp_py, 0, myPyString_FromString(info.comp_id.c_str()));
      PyList_SetItem(chem_comp_py, 1, myPyString_FromString(info.three_letter_code.c_str()));
      PyList_SetItem(chem_comp_py, 2, myPyString_FromString(info.name.c_str()));
      PyList_SetItem(chem_comp_py, 3, myPyString_FromString(info.group.c_str()));
      PyList_SetItem(chem_comp_py, 4, PyLong_FromLong(info.number_atoms_all));
      PyList_SetItem(chem_comp_py, 5, PyLong_FromLong(info.number_atoms_nh));
      PyList_SetItem(chem_comp_py, 6, myPyString_FromString(info.description_level.c_str()));

      // Put chem_comp_py into a dictionary?
      PyDict_SetItem(r, myPyString_FromString("_chem_comp"), chem_comp_py);


      // ------------------ chem_comp_atom -------------------------
      std::vector<coot::dict_atom> atom_info = restraints.atom_info;
      int n_atoms = atom_info.size();
      PyObject *atom_info_list = PyList_New(n_atoms);
      for (int iat=0; iat<n_atoms; iat++) {
         PyObject *atom_attributes_list = PyList_New(5);
         PyList_SetItem(atom_attributes_list, 0, myPyString_FromString(atom_info[iat].atom_id_4c.c_str()));
         PyList_SetItem(atom_attributes_list, 1, myPyString_FromString(atom_info[iat].type_symbol.c_str()));
         PyList_SetItem(atom_attributes_list, 2, myPyString_FromString(atom_info[iat].type_energy.c_str()));
         PyList_SetItem(atom_attributes_list, 3, PyFloat_FromDouble(atom_info[iat].partial_charge.second));
         PyObject *flag = Py_False;
         if (atom_info[iat].partial_charge.first)
            flag = Py_True;
     Py_INCREF(flag);
         PyList_SetItem(atom_attributes_list, 4, flag);
         PyList_SetItem(atom_info_list, iat, atom_attributes_list);
      }

      PyDict_SetItem(r, myPyString_FromString("_chem_comp_atom"), atom_info_list);

      // ------------------ Bonds -------------------------
      PyObject *bond_restraint_list = PyList_New(restraints.bond_restraint.size());
      for (unsigned int ibond=0; ibond<restraints.bond_restraint.size(); ibond++) {
         std::string a1   = restraints.bond_restraint[ibond].atom_id_1_4c();
         std::string a2   = restraints.bond_restraint[ibond].atom_id_2_4c();
         std::string type = restraints.bond_restraint[ibond].type();

         PyObject *py_value_dist = Py_False;
         PyObject *py_value_esd = Py_False;

         try {
            double d   = restraints.bond_restraint[ibond].value_dist();
            double esd = restraints.bond_restraint[ibond].value_esd();
            py_value_dist = PyFloat_FromDouble(d);
            py_value_esd  = PyFloat_FromDouble(esd);
         }
         catch (const std::runtime_error &rte) {

            // Use default false values.
            // So I suppose that I need to do this then:
            if (PyBool_Check(py_value_dist))
               Py_INCREF(py_value_dist);
            if (PyBool_Check(py_value_esd))
               Py_INCREF(py_value_esd);
         }

         PyObject *bond_restraint = PyList_New(5);
         PyList_SetItem(bond_restraint, 0, myPyString_FromString(a1.c_str()));
         PyList_SetItem(bond_restraint, 1, myPyString_FromString(a2.c_str()));
         PyList_SetItem(bond_restraint, 2, myPyString_FromString(type.c_str()));
         PyList_SetItem(bond_restraint, 3, py_value_dist);
         PyList_SetItem(bond_restraint, 4, py_value_esd);
         PyList_SetItem(bond_restraint_list, ibond, bond_restraint);
      }

      PyDict_SetItem(r, myPyString_FromString("_chem_comp_bond"), bond_restraint_list);


      // ------------------ Angles -------------------------
      PyObject *angle_restraint_list = PyList_New(restraints.angle_restraint.size());
      for (unsigned int iangle=0; iangle<restraints.angle_restraint.size(); iangle++) {
         std::string a1 = restraints.angle_restraint[iangle].atom_id_1_4c();
         std::string a2 = restraints.angle_restraint[iangle].atom_id_2_4c();
         std::string a3 = restraints.angle_restraint[iangle].atom_id_3_4c();
         double d   = restraints.angle_restraint[iangle].angle();
         double esd = restraints.angle_restraint[iangle].esd();
         PyObject *angle_restraint = PyList_New(5);
         PyList_SetItem(angle_restraint, 0, myPyString_FromString(a1.c_str()));
         PyList_SetItem(angle_restraint, 1, myPyString_FromString(a2.c_str()));
         PyList_SetItem(angle_restraint, 2, myPyString_FromString(a3.c_str()));
         PyList_SetItem(angle_restraint, 3, PyFloat_FromDouble(d));
         PyList_SetItem(angle_restraint, 4, PyFloat_FromDouble(esd));
         PyList_SetItem(angle_restraint_list, iangle, angle_restraint);
      }

      PyDict_SetItem(r, myPyString_FromString("_chem_comp_angle"), angle_restraint_list);


      // ------------------ Torsions -------------------------
      PyObject *torsion_restraint_list = PyList_New(restraints.torsion_restraint.size());
      for (unsigned int itorsion=0; itorsion<restraints.torsion_restraint.size(); itorsion++) {
         std::string id = restraints.torsion_restraint[itorsion].id();
         std::string a1 = restraints.torsion_restraint[itorsion].atom_id_1_4c();
         std::string a2 = restraints.torsion_restraint[itorsion].atom_id_2_4c();
         std::string a3 = restraints.torsion_restraint[itorsion].atom_id_3_4c();
         std::string a4 = restraints.torsion_restraint[itorsion].atom_id_4_4c();
         double tor  = restraints.torsion_restraint[itorsion].angle();
         double esd = restraints.torsion_restraint[itorsion].esd();
         int period = restraints.torsion_restraint[itorsion].periodicity();
         PyObject *torsion_restraint = PyList_New(8);
         PyList_SetItem(torsion_restraint, 0, myPyString_FromString(id.c_str()));
         PyList_SetItem(torsion_restraint, 1, myPyString_FromString(a1.c_str()));
         PyList_SetItem(torsion_restraint, 2, myPyString_FromString(a2.c_str()));
         PyList_SetItem(torsion_restraint, 3, myPyString_FromString(a3.c_str()));
         PyList_SetItem(torsion_restraint, 4, myPyString_FromString(a4.c_str()));
         PyList_SetItem(torsion_restraint, 5, PyFloat_FromDouble(tor));
         PyList_SetItem(torsion_restraint, 6, PyFloat_FromDouble(esd));
         PyList_SetItem(torsion_restraint, 7, PyLong_FromLong(period));
         PyList_SetItem(torsion_restraint_list, itorsion, torsion_restraint);
      }

      PyDict_SetItem(r, myPyString_FromString("_chem_comp_tor"), torsion_restraint_list);

      // ------------------ Planes -------------------------
      PyObject *plane_restraints_list = PyList_New(restraints.plane_restraint.size());
      for (unsigned int iplane=0; iplane<restraints.plane_restraint.size(); iplane++) {
         PyObject *atom_list = PyList_New(restraints.plane_restraint[iplane].n_atoms());
         for (int iat=0; iat<restraints.plane_restraint[iplane].n_atoms(); iat++) {
            std::string at = restraints.plane_restraint[iplane][iat].first;
            PyList_SetItem(atom_list, iat, myPyString_FromString(at.c_str()));
         }
         double esd = restraints.plane_restraint[iplane].dist_esd(0);
         PyObject *plane_restraint = PyList_New(3);
         PyList_SetItem(plane_restraint, 0, myPyString_FromString(restraints.plane_restraint[iplane].plane_id.c_str()));
         PyList_SetItem(plane_restraint, 1, atom_list);
         PyList_SetItem(plane_restraint, 2, PyFloat_FromDouble(esd));
         PyList_SetItem(plane_restraints_list, iplane, plane_restraint);
      }

      PyDict_SetItem(r, myPyString_FromString("_chem_comp_plane_atom"), plane_restraints_list);

      // ------------------ Chirals -------------------------
      PyObject *chiral_restraint_list = PyList_New(restraints.chiral_restraint.size());
      for (unsigned int ichiral=0; ichiral<restraints.chiral_restraint.size(); ichiral++) {

         std::string a1 = restraints.chiral_restraint[ichiral].atom_id_1_4c();
         std::string a2 = restraints.chiral_restraint[ichiral].atom_id_2_4c();
         std::string a3 = restraints.chiral_restraint[ichiral].atom_id_3_4c();
         std::string ac = restraints.chiral_restraint[ichiral].atom_id_c_4c();
         std::string chiral_id = restraints.chiral_restraint[ichiral].Chiral_Id();

         double esd = restraints.chiral_restraint[ichiral].volume_sigma();
         int volume_sign = restraints.chiral_restraint[ichiral].volume_sign;
         PyObject *chiral_restraint = PyList_New(7);
         PyList_SetItem(chiral_restraint, 0, myPyString_FromString(chiral_id.c_str()));
         PyList_SetItem(chiral_restraint, 1, myPyString_FromString(ac.c_str()));
         PyList_SetItem(chiral_restraint, 2, myPyString_FromString(a1.c_str()));
         PyList_SetItem(chiral_restraint, 3, myPyString_FromString(a2.c_str()));
         PyList_SetItem(chiral_restraint, 4, myPyString_FromString(a3.c_str()));
         PyList_SetItem(chiral_restraint, 5, PyLong_FromLong(volume_sign));
         PyList_SetItem(chiral_restraint, 6, PyFloat_FromDouble(esd));
         PyList_SetItem(chiral_restraint_list, ichiral, chiral_restraint);
      }

      PyDict_SetItem(r, myPyString_FromString("_chem_comp_chir"), chiral_restraint_list);
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}
#endif // USE_PYTHON

#ifdef USE_GUILE
SCM set_monomer_restraints(const char *monomer_type, SCM restraints) {

   SCM retval = SCM_BOOL_F;

   if (!scm_is_true(scm_list_p(restraints))) {
      std::cout << " Failed to read restraints - not a list" << std::endl;
   } else {

      std::vector<coot::dict_bond_restraint_t> bond_restraints;
      std::vector<coot::dict_angle_restraint_t> angle_restraints;
      std::vector<coot::dict_torsion_restraint_t> torsion_restraints;
      std::vector<coot::dict_chiral_restraint_t> chiral_restraints;
      std::vector<coot::dict_plane_restraint_t> plane_restraints;
      std::vector<coot::dict_atom> atoms;
      coot::dict_chem_comp_t residue_info;

      SCM restraints_length_scm = scm_length(restraints);
      int restraints_length = scm_to_int(restraints_length_scm);
      if (restraints_length > 0) {
         for (int i_rest_type=0; i_rest_type<restraints_length; i_rest_type++) {
            SCM rest_container = scm_list_ref(restraints, scm_from_int(i_rest_type));
            if (scm_is_true(scm_list_p(rest_container))) {
               SCM rest_container_length_scm = scm_length(rest_container);
               int rest_container_length = scm_to_int(rest_container_length_scm);
               if (rest_container_length > 1) {
                  SCM restraints_type_scm = SCM_CAR(rest_container);
                  if (scm_string_p(restraints_type_scm)) {
                     std::string restraints_type = scm_to_locale_string(restraints_type_scm);

                     if (restraints_type == "_chem_comp") {
                        SCM chem_comp_info_scm = SCM_CDR(rest_container);
                        SCM chem_comp_info_length_scm = scm_length(chem_comp_info_scm);
                        int chem_comp_info_length = scm_to_int(chem_comp_info_length_scm);

                        if (chem_comp_info_length != 7) {
                           std::cout << "WARNING:: chem_comp_info length " << chem_comp_info_length
                                     << " should be " << 7 << std::endl;
                        } else {
                           SCM  comp_id_scm = scm_list_ref(chem_comp_info_scm, scm_from_int(0));
                           SCM      tlc_scm = scm_list_ref(chem_comp_info_scm, scm_from_int(1));
                           SCM     name_scm = scm_list_ref(chem_comp_info_scm, scm_from_int(2));
                           SCM    group_scm = scm_list_ref(chem_comp_info_scm, scm_from_int(3));
                           SCM      noa_scm = scm_list_ref(chem_comp_info_scm, scm_from_int(4));
                           SCM    nonha_scm = scm_list_ref(chem_comp_info_scm, scm_from_int(5));
                           SCM desc_lev_scm = scm_list_ref(chem_comp_info_scm, scm_from_int(6));
                           if (scm_is_true(scm_string_p(comp_id_scm)) &&
                               scm_is_true(scm_string_p(tlc_scm)) &&
                               scm_is_true(scm_string_p(name_scm)) &&
                               scm_is_true(scm_string_p(group_scm)) &&
                               scm_is_true(scm_number_p(noa_scm)) &&
                               scm_is_true(scm_number_p(nonha_scm)) &&
                               scm_is_true(scm_string_p(desc_lev_scm))) {
                              std::string comp_id = scm_to_locale_string(comp_id_scm);
                              std::string     tlc = scm_to_locale_string(tlc_scm);
                              std::string    name = scm_to_locale_string(name_scm);
                              std::string   group = scm_to_locale_string(group_scm);
                              std::string des_lev = scm_to_locale_string(desc_lev_scm);
                              int no_of_atoms = scm_to_int(noa_scm);
                              int no_of_non_H_atoms = scm_to_int(nonha_scm);
                              coot::dict_chem_comp_t n(comp_id, tlc, name, group,
                                                       no_of_atoms, no_of_non_H_atoms,
                                                       des_lev);
                              residue_info = n;
                           }
                        }
                     }

                     if (restraints_type == "_chem_comp_atom") {
                        SCM chem_comp_atoms = SCM_CDR(rest_container);
                        SCM chem_comp_atoms_length_scm = scm_length(chem_comp_atoms);
                        int chem_comp_atoms_length = scm_to_int(chem_comp_atoms_length_scm);

                        for (int iat=0; iat<chem_comp_atoms_length; iat++) {
                           SCM chem_comp_atom_scm = scm_list_ref(chem_comp_atoms, scm_from_int(iat));
                           SCM chem_comp_atom_length_scm = scm_length(chem_comp_atom_scm);
                           int chem_comp_atom_length = scm_to_int(chem_comp_atom_length_scm);

                           if (chem_comp_atom_length != 5) {
                              std::cout << "WARNING:: chem_comp_atom length " << chem_comp_atom_length
                                        << " should be " << 5 << std::endl;
                           } else {
                              SCM atom_id_scm  = scm_list_ref(chem_comp_atom_scm, scm_from_int(0));
                              SCM element_scm  = scm_list_ref(chem_comp_atom_scm, scm_from_int(1));
                              SCM energy_scm   = scm_list_ref(chem_comp_atom_scm, scm_from_int(2));
                              SCM partial_charge_scm = scm_list_ref(chem_comp_atom_scm, scm_from_int(3));
                              SCM valid_pc_scm = scm_list_ref(chem_comp_atom_scm, scm_from_int(4));

                              if (scm_string_p(atom_id_scm) && scm_string_p(element_scm) &&
                                  scm_number_p(partial_charge_scm)) {
                                 std::string atom_id(scm_to_locale_string(atom_id_scm));
                                 std::string element(scm_to_locale_string(element_scm));
                                 std::string energy(scm_to_locale_string(energy_scm));
                                 float partial_charge = scm_to_double(partial_charge_scm);
                                 short int valid_partial_charge = 1;
                                 if SCM_FALSEP(valid_pc_scm)
                                    valid_partial_charge = 0;
                                 coot::dict_atom at(atom_id, atom_id, element, energy,
                                                    std::pair<bool, float>(valid_partial_charge,
                                                                           partial_charge));

                                 atoms.push_back(at);
                              }
                           }
                        }

                     }


                     if (restraints_type == "_chem_comp_bond") {
                        SCM bond_restraints_list_scm = SCM_CDR(rest_container);
                        SCM bond_restraints_list_length_scm = scm_length(bond_restraints_list_scm);
                        int bond_restraints_list_length = scm_to_int(bond_restraints_list_length_scm);

                        for (int ibr=0; ibr<bond_restraints_list_length; ibr++) {
                           SCM bond_restraint = scm_list_ref(bond_restraints_list_scm, scm_from_int(ibr));
                           SCM bond_restraint_length_scm = scm_length(bond_restraint);
                           int bond_restraint_length = scm_to_int(bond_restraint_length_scm);

                           if (bond_restraint_length != 5) {
                              std::cout << "WARNING:: bond_restraint_length " << bond_restraint_length
                                        << " should be " << 5 << std::endl;
                           } else {
                              SCM atom_1_scm = scm_list_ref(bond_restraint, scm_from_int(0));
                              SCM atom_2_scm = scm_list_ref(bond_restraint, scm_from_int(1));
                              SCM type_scm   = scm_list_ref(bond_restraint, scm_from_int(2));
                              SCM dist_scm   = scm_list_ref(bond_restraint, scm_from_int(3));
                              SCM esd_scm    = scm_list_ref(bond_restraint, scm_from_int(4));
                              if (scm_string_p(atom_1_scm) && scm_string_p(atom_2_scm) &&
                                  scm_number_p(dist_scm) && scm_number_p(esd_scm)) {
                                 std::string atom_1 = scm_to_locale_string(atom_1_scm);
                                 std::string atom_2 = scm_to_locale_string(atom_2_scm);
                                 std::string type   = scm_to_locale_string(type_scm);
                                 double dist        = scm_to_double(dist_scm);
                                 double esd         = scm_to_double(esd_scm);
                                 coot::dict_bond_restraint_t rest(atom_1, atom_2, type, dist, esd, 0.0, 0.0, false);
                                 bond_restraints.push_back(rest);
                              }
                           }
                        }
                     }

                     if (restraints_type == "_chem_comp_angle") {
                        SCM angle_restraints_list = SCM_CDR(rest_container);
                        SCM angle_restraints_list_length_scm = scm_length(angle_restraints_list);
                        int angle_restraints_list_length = scm_to_int(angle_restraints_list_length_scm);

                        for (int iar=0; iar<angle_restraints_list_length; iar++) {
                           SCM angle_restraint = scm_list_ref(angle_restraints_list, scm_from_int(iar));
                           SCM angle_restraint_length_scm = scm_length(angle_restraint);
                           int angle_restraint_length = scm_to_int(angle_restraint_length_scm);

                           if (angle_restraint_length != 5) {
                              std::cout << "WARNING:: angle_restraint_length length "
                                        << angle_restraint_length << " should be " << 5 << std::endl;
                           } else {
                              SCM atom_1_scm = scm_list_ref(angle_restraint, scm_from_int(0));
                              SCM atom_2_scm = scm_list_ref(angle_restraint, scm_from_int(1));
                              SCM atom_3_scm = scm_list_ref(angle_restraint, scm_from_int(2));
                              SCM angle_scm  = scm_list_ref(angle_restraint, scm_from_int(3));
                              SCM esd_scm    = scm_list_ref(angle_restraint, scm_from_int(4));
                              if (scm_string_p(atom_1_scm) && scm_string_p(atom_2_scm) &&
                                  scm_string_p(atom_3_scm) &&
                                  scm_number_p(angle_scm) && scm_number_p(esd_scm)) {
                                 std::string atom_1 = scm_to_locale_string(atom_1_scm);
                                 std::string atom_2 = scm_to_locale_string(atom_2_scm);
                                 std::string atom_3 = scm_to_locale_string(atom_3_scm);
                                 double angle       = scm_to_double(angle_scm);
                                 double esd         = scm_to_double(esd_scm);
                                 coot::dict_angle_restraint_t rest(atom_1, atom_2, atom_3, angle, esd);
                                 angle_restraints.push_back(rest);
                              }
                           }
                        }
                     }


                     if (restraints_type == "_chem_comp_tor") {
                        SCM torsion_restraints_list = SCM_CDR(rest_container);
                        SCM torsion_restraints_list_length_scm = scm_length(torsion_restraints_list);
                        int torsion_restraints_list_length = scm_to_int(torsion_restraints_list_length_scm);

                        for (int itr=0; itr<torsion_restraints_list_length; itr++) {
                           SCM torsion_restraint = scm_list_ref(torsion_restraints_list, scm_from_int(itr));
                           SCM torsion_restraint_length_scm = scm_length(torsion_restraint);
                           int torsion_restraint_length = scm_to_int(torsion_restraint_length_scm);

                           if (torsion_restraint_length == 8) {
                              SCM torsion_id_scm = scm_list_ref(torsion_restraint, scm_from_int(0));
                              SCM atom_1_scm     = scm_list_ref(torsion_restraint, scm_from_int(1));
                              SCM atom_2_scm     = scm_list_ref(torsion_restraint, scm_from_int(2));
                              SCM atom_3_scm     = scm_list_ref(torsion_restraint, scm_from_int(3));
                              SCM atom_4_scm     = scm_list_ref(torsion_restraint, scm_from_int(4));
                              SCM torsion_scm    = scm_list_ref(torsion_restraint, scm_from_int(5));
                              SCM esd_scm        = scm_list_ref(torsion_restraint, scm_from_int(6));
                              SCM period_scm     = scm_list_ref(torsion_restraint, scm_from_int(7));
                              if (scm_is_true(scm_string_p(atom_1_scm)) &&
                                  scm_is_true(scm_string_p(atom_2_scm)) &&
                                  scm_is_true(scm_string_p(atom_3_scm)) &&
                                  scm_is_true(scm_string_p(atom_4_scm)) &&
                                  scm_is_true(scm_number_p(torsion_scm)) &&
                                  scm_is_true(scm_number_p(esd_scm)) &&
                                  scm_is_true(scm_number_p(period_scm))) {
                                 std::string torsion_id = scm_to_locale_string(torsion_id_scm);
                                 std::string atom_1     = scm_to_locale_string(atom_1_scm);
                                 std::string atom_2     = scm_to_locale_string(atom_2_scm);
                                 std::string atom_3     = scm_to_locale_string(atom_3_scm);
                                 std::string atom_4     = scm_to_locale_string(atom_4_scm);
                                 double torsion         = scm_to_double(torsion_scm);
                                 double esd             = scm_to_double(esd_scm);
                                 int period             = scm_to_int(period_scm);
                                 coot::dict_torsion_restraint_t rest(torsion_id,
                                                                     atom_1, atom_2, atom_3, atom_4,
                                                                     torsion, esd, period);
                                 torsion_restraints.push_back(rest);
                              }
                           }
                        }
                     }

                     if (restraints_type == "_chem_comp_plane_atom") {
                        SCM plane_restraints_list = SCM_CDR(rest_container);
                        SCM plane_restraints_list_length_scm = scm_length(plane_restraints_list);
                        int plane_restraints_list_length = scm_to_int(plane_restraints_list_length_scm);

                        for (int ipr=0; ipr<plane_restraints_list_length; ipr++) {
                           SCM plane_restraint = scm_list_ref(plane_restraints_list, scm_from_int(ipr));
                           SCM plane_restraint_length_scm = scm_length(plane_restraint);
                           int plane_restraint_length = scm_to_int(plane_restraint_length_scm);

                           if (plane_restraint_length == 3) {

                              std::vector<SCM> plane_atoms;
                              SCM plane_id_scm   = scm_list_ref(plane_restraint, scm_from_int(0));
                              SCM esd_scm        = scm_list_ref(plane_restraint, scm_from_int(2));
                              SCM atom_list_scm  = scm_list_ref(plane_restraint, scm_from_int(1));
                              SCM atom_list_length_scm = scm_length(atom_list_scm);
                              int atom_list_length = scm_to_int(atom_list_length_scm);
                              bool atoms_pass = 1;
                              for (int iat=0; iat<atom_list_length; iat++) {
                                 SCM atom_scm   = scm_list_ref(atom_list_scm, scm_from_int(iat));
                                 plane_atoms.push_back(atom_scm);
                                 if (!scm_string_p(atom_scm))
                                    atoms_pass = 0;
                              }

                              if (atoms_pass && scm_string_p(plane_id_scm) &&  scm_number_p(esd_scm)) {
                                 std::vector<std::string> atom_names;
                                 for (unsigned int i=0; i<plane_atoms.size(); i++)
                                    atom_names.push_back(std::string(scm_to_locale_string(plane_atoms[i])));

                                 std::string plane_id = scm_to_locale_string(plane_id_scm);
                                 double esd           = scm_to_double(esd_scm);
                                 if (atom_names.size() > 0) {
                                    coot::dict_plane_restraint_t rest(plane_id, atom_names[0], esd);

                                    for (unsigned int i=1; i<atom_names.size(); i++) {
                                       rest.push_back_atom(atom_names[i], esd);
                                    }
                                    plane_restraints.push_back(rest);
                                    // std::cout << "plane restraint: " << rest << std::endl;
                                 }
                              }
                           }
                        }
                     }


                     if (restraints_type == "_chem_comp_chir") {
                        SCM chiral_restraints_list = SCM_CDR(rest_container);
                        SCM chiral_restraints_list_length_scm = scm_length(chiral_restraints_list);
                        int chiral_restraints_list_length = scm_to_int(chiral_restraints_list_length_scm);

                        for (int icr=0; icr<chiral_restraints_list_length; icr++) {
                           SCM chiral_restraint = scm_list_ref(chiral_restraints_list, scm_from_int(icr));
                           SCM chiral_restraint_length_scm = scm_length(chiral_restraint);
                           int chiral_restraint_length = scm_to_int(chiral_restraint_length_scm);

                           if (chiral_restraint_length == 7) {
                              SCM chiral_id_scm= scm_list_ref(chiral_restraint, scm_from_int(0));
                              SCM atom_c_scm   = scm_list_ref(chiral_restraint, scm_from_int(1));
                              SCM atom_1_scm   = scm_list_ref(chiral_restraint, scm_from_int(2));
                              SCM atom_2_scm   = scm_list_ref(chiral_restraint, scm_from_int(3));
                              SCM atom_3_scm   = scm_list_ref(chiral_restraint, scm_from_int(4));
                              SCM chiral_vol_sign_scm = scm_list_ref(chiral_restraint, scm_from_int(5));
                              // SCM esd_scm      = scm_list_ref(chiral_restraint, scm_from_int(6));
                              if (scm_string_p(atom_1_scm) && scm_string_p(atom_2_scm) &&
                                  scm_string_p(atom_3_scm) && scm_string_p(atom_c_scm)) {
                                 std::string chiral_id = scm_to_locale_string(chiral_id_scm);
                                 std::string atom_1 = scm_to_locale_string(atom_1_scm);
                                 std::string atom_2 = scm_to_locale_string(atom_2_scm);
                                 std::string atom_3 = scm_to_locale_string(atom_3_scm);
                                 std::string atom_c = scm_to_locale_string(atom_c_scm);
                                 // double esd         = scm_to_double(esd_scm);
                                 int chiral_vol_sign= scm_to_int(chiral_vol_sign_scm);
                                 coot::dict_chiral_restraint_t rest(chiral_id,
                                                                    atom_c, atom_1, atom_2, atom_3,
                                                                    chiral_vol_sign);

                                 chiral_restraints.push_back(rest);
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }

//       std::cout << "Found " <<    bond_restraints.size() << "   bond  restraints" << std::endl;
//       std::cout << "Found " <<   angle_restraints.size() << "   angle restraints" << std::endl;
//       std::cout << "Found " << torsion_restraints.size() << " torsion restraints" << std::endl;
//       std::cout << "Found " <<   plane_restraints.size() << "   plane restraints" << std::endl;
//       std::cout << "Found " <<  chiral_restraints.size() << "  chiral restraints" << std::endl;

      graphics_info_t g;
      coot::dictionary_residue_restraints_t monomer_restraints(monomer_type,
                                                               g.cif_dictionary_read_number++);
      monomer_restraints.bond_restraint    = bond_restraints;
      monomer_restraints.angle_restraint   = angle_restraints;
      monomer_restraints.torsion_restraint = torsion_restraints;
      monomer_restraints.chiral_restraint  = chiral_restraints;
      monomer_restraints.plane_restraint   = plane_restraints;
      monomer_restraints.residue_info      = residue_info;
      monomer_restraints.atom_info         = atoms;

      // the imol is not specified, we want to replace the restraints that
      // can be used by any molecule
      int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
      bool s = g.Geom_p()->replace_monomer_restraints(monomer_type, imol_enc, monomer_restraints);
      if (s)
         retval = SCM_BOOL_T;
   }

   return retval;
}
#endif // USE_GUILE

#ifdef USE_PYTHON
PyObject *set_monomer_restraints_py(const char *monomer_type, PyObject *restraints) {

   PyObject *retval = Py_False;

   if (!PyDict_Check(restraints)) {
      std::cout << " Failed to read restraints - not a list" << std::endl;
   } else {

      std::vector<coot::dict_bond_restraint_t> bond_restraints;
      std::vector<coot::dict_angle_restraint_t> angle_restraints;
      std::vector<coot::dict_torsion_restraint_t> torsion_restraints;
      std::vector<coot::dict_chiral_restraint_t> chiral_restraints;
      std::vector<coot::dict_plane_restraint_t> plane_restraints;
      std::vector<coot::dict_atom> atoms;
      coot::dict_chem_comp_t residue_info;

      PyObject *key;
      PyObject *value;
      Py_ssize_t pos = 0;

      std::cout << "looping over restraint" << std::endl;
      while (PyDict_Next(restraints, &pos, &key, &value)) {
         // std::cout << ":::::::key: " << PyUnicode_AsUTF8String(key) << std::endl;

         std::string key_string = PyBytes_AS_STRING(PyUnicode_AsUTF8String(key));
         if (key_string == "_chem_comp") {
            PyObject *chem_comp_list = value;
            if (PyList_Check(chem_comp_list)) {
               if (PyObject_Length(chem_comp_list) == 7) {
                  std::string comp_id  = PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyList_GetItem(chem_comp_list, 0)));
                  std::string tlc      = PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyList_GetItem(chem_comp_list, 1)));
                  std::string name     = PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyList_GetItem(chem_comp_list, 2)));
                  std::string group    = PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyList_GetItem(chem_comp_list, 3)));
                  int n_atoms_all      = PyLong_AsLong(PyList_GetItem(chem_comp_list, 4));
                  int n_atoms_nh       = PyLong_AsLong(PyList_GetItem(chem_comp_list, 5));
                  std::string desc_lev = PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyList_GetItem(chem_comp_list, 6)));

                  coot::dict_chem_comp_t n(comp_id, tlc, name, group,
                                           n_atoms_all, n_atoms_nh, desc_lev);
                  residue_info = n;
               }
            }
         }


         if (key_string == "_chem_comp_atom") {
            PyObject *chem_comp_atom_list = value;
            if (PyList_Check(chem_comp_atom_list)) {
               int n_atoms = PyObject_Length(chem_comp_atom_list);
               for (int iat=0; iat<n_atoms; iat++) {
                  PyObject *chem_comp_atom = PyList_GetItem(chem_comp_atom_list, iat);
                  if (PyObject_Length(chem_comp_atom) == 5) {
                     std::string atom_id  = PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyList_GetItem(chem_comp_atom, 0)));
                     std::string element  = PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyList_GetItem(chem_comp_atom, 1)));
                     std::string energy_t = PyBytes_AS_STRING(PyUnicode_AsUTF8String(PyList_GetItem(chem_comp_atom, 2)));
                     float part_chr        = PyFloat_AsDouble(PyList_GetItem(chem_comp_atom, 3));
                     bool flag = 0;
                     if (PyLong_AsLong(PyList_GetItem(chem_comp_atom, 4))) {
                        flag = 1;
                     }
                     std::pair<bool, float> part_charge_info(flag, part_chr);
                     coot::dict_atom at(atom_id, atom_id, element, energy_t, part_charge_info);
                     atoms.push_back(at);
                  }
               }
            }
         }

         if (key_string == "_chem_comp_bond") {
            PyObject *bond_restraint_list = value;
            if (PyList_Check(bond_restraint_list)) {
               int n_bonds = PyObject_Length(bond_restraint_list);
               for (int i_bond=0; i_bond<n_bonds; i_bond++) {
                  PyObject *bond_restraint = PyList_GetItem(bond_restraint_list, i_bond);
                  if (PyObject_Length(bond_restraint) == 5) {
                     PyObject *atom_1_py = PyList_GetItem(bond_restraint, 0);
                     PyObject *atom_2_py = PyList_GetItem(bond_restraint, 1);
                     PyObject *type_py   = PyList_GetItem(bond_restraint, 2);
                     PyObject *dist_py   = PyList_GetItem(bond_restraint, 3);
                     PyObject *esd_py    = PyList_GetItem(bond_restraint, 4);

                     if (PyUnicode_Check(atom_1_py) &&
                         PyUnicode_Check(atom_2_py) &&
                         PyUnicode_Check(type_py) &&
                         PyFloat_Check(dist_py) &&
                         PyFloat_Check(esd_py)) {
                        std::string atom_1 = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_1_py));
                        std::string atom_2 = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_2_py));
                        std::string type   = PyBytes_AS_STRING(PyUnicode_AsUTF8String(type_py));
                        float  dist = PyFloat_AsDouble(dist_py);
                        float  esd  = PyFloat_AsDouble(esd_py);
                        coot::dict_bond_restraint_t rest(atom_1, atom_2, type, dist, esd, 0.0, 0.0, false);
                        bond_restraints.push_back(rest);
                     }
                  }
               }
            }
         }


         if (key_string == "_chem_comp_angle") {
            PyObject *angle_restraint_list = value;
            if (PyList_Check(angle_restraint_list)) {
               int n_angles = PyObject_Length(angle_restraint_list);
               for (int i_angle=0; i_angle<n_angles; i_angle++) {
                  PyObject *angle_restraint = PyList_GetItem(angle_restraint_list, i_angle);
                  if (PyObject_Length(angle_restraint) == 5) {
                     PyObject *atom_1_py = PyList_GetItem(angle_restraint, 0);
                     PyObject *atom_2_py = PyList_GetItem(angle_restraint, 1);
                     PyObject *atom_3_py = PyList_GetItem(angle_restraint, 2);
                     PyObject *angle_py  = PyList_GetItem(angle_restraint, 3);
                     PyObject *esd_py    = PyList_GetItem(angle_restraint, 4);

                     if (PyUnicode_Check(atom_1_py) &&
                         PyUnicode_Check(atom_2_py) &&
                         PyUnicode_Check(atom_3_py) &&
                         PyFloat_Check(angle_py) &&
                         PyFloat_Check(esd_py)) {
                        std::string atom_1 = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_1_py));
                        std::string atom_2 = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_2_py));
                        std::string atom_3 = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_3_py));
                        float  angle = PyFloat_AsDouble(angle_py);
                        float  esd   = PyFloat_AsDouble(esd_py);
                        coot::dict_angle_restraint_t rest(atom_1, atom_2, atom_3, angle, esd);
                        angle_restraints.push_back(rest);
                     }
                  }
               }
            }
         }

         if (key_string == "_chem_comp_tor") {
            PyObject *torsion_restraint_list = value;
            if (PyList_Check(torsion_restraint_list)) {
               int n_torsions = PyObject_Length(torsion_restraint_list);
               for (int i_torsion=0; i_torsion<n_torsions; i_torsion++) {
                  PyObject *torsion_restraint = PyList_GetItem(torsion_restraint_list, i_torsion);
                  if (PyObject_Length(torsion_restraint) == 7) { // info for Nigel.
                     std::cout << "torsions now have 8 elements starting with the torsion id\n";
                  }
                  if (PyObject_Length(torsion_restraint) == 8) {
                     PyObject *id_py     = PyList_GetItem(torsion_restraint, 0);
                     PyObject *atom_1_py = PyList_GetItem(torsion_restraint, 1);
                     PyObject *atom_2_py = PyList_GetItem(torsion_restraint, 2);
                     PyObject *atom_3_py = PyList_GetItem(torsion_restraint, 3);
                     PyObject *atom_4_py = PyList_GetItem(torsion_restraint, 4);
                     PyObject *torsion_py= PyList_GetItem(torsion_restraint, 5);
                     PyObject *esd_py    = PyList_GetItem(torsion_restraint, 6);
                     PyObject *period_py = PyList_GetItem(torsion_restraint, 7);

                     if (PyUnicode_Check(atom_1_py) &&
                         PyUnicode_Check(atom_2_py) &&
                         PyUnicode_Check(atom_3_py) &&
                         PyUnicode_Check(atom_4_py) &&
                         PyFloat_Check(torsion_py) &&
                         PyFloat_Check(esd_py)    &&
                         PyLong_Check(period_py)) {
                        std::string id     = PyBytes_AS_STRING(PyUnicode_AsUTF8String(id_py));
                        std::string atom_1 = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_1_py));
                        std::string atom_2 = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_2_py));
                        std::string atom_3 = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_3_py));
                        std::string atom_4 = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_4_py));
                        float  torsion = PyFloat_AsDouble(torsion_py);
                        float  esd     = PyFloat_AsDouble(esd_py);
                        int  period    = PyLong_AsLong(period_py);
                        coot::dict_torsion_restraint_t rest(id, atom_1, atom_2, atom_3, atom_4,
                                                            torsion, esd, period);
                        torsion_restraints.push_back(rest);
                     }
                  }
               }
            }
         }

         if (key_string == "_chem_comp_chir") {
            PyObject *chiral_restraint_list = value;
            if (PyList_Check(chiral_restraint_list)) {
               int n_chirals = PyObject_Length(chiral_restraint_list);
               for (int i_chiral=0; i_chiral<n_chirals; i_chiral++) {
                  PyObject *chiral_restraint = PyList_GetItem(chiral_restraint_list, i_chiral);
                  if (PyObject_Length(chiral_restraint) == 7) {
                     PyObject *chiral_id_py= PyList_GetItem(chiral_restraint, 0);
                     PyObject *atom_c_py   = PyList_GetItem(chiral_restraint, 1);
                     PyObject *atom_1_py   = PyList_GetItem(chiral_restraint, 2);
                     PyObject *atom_2_py   = PyList_GetItem(chiral_restraint, 3);
                     PyObject *atom_3_py   = PyList_GetItem(chiral_restraint, 4);
                     PyObject *vol_sign_py = PyList_GetItem(chiral_restraint, 5);
                     PyObject *esd_py      = PyList_GetItem(chiral_restraint, 6);

                     if (PyUnicode_Check(atom_1_py) &&
                         PyUnicode_Check(atom_2_py) &&
                         PyUnicode_Check(atom_3_py) &&
                         PyUnicode_Check(atom_c_py) &&
                         PyUnicode_Check(chiral_id_py) &&
                         PyFloat_Check(esd_py)    &&
                         PyLong_Check(vol_sign_py)) {
                        std::string chiral_id = PyBytes_AS_STRING(PyUnicode_AsUTF8String(chiral_id_py));
                        std::string atom_c    = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_c_py));
                        std::string atom_1    = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_1_py));
                        std::string atom_2    = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_2_py));
                        std::string atom_3    = PyBytes_AS_STRING(PyUnicode_AsUTF8String(atom_3_py));
                        float  esd            = PyFloat_AsDouble(esd_py);
                        int  vol_sign         = PyLong_AsLong(vol_sign_py);
                        coot::dict_chiral_restraint_t rest(chiral_id,
                                                           atom_c, atom_1, atom_2, atom_3,
                                                           vol_sign);
                        chiral_restraints.push_back(rest);
                     }
                  }
               }
            }
         }


    if (key_string == "_chem_comp_plane_atom") {
       PyObject *plane_restraint_list = value;
       if (PyList_Check(plane_restraint_list)) {
          int n_planes = PyObject_Length(plane_restraint_list);
          for (int i_plane=0; i_plane<n_planes; i_plane++) {
             PyObject *plane_restraint = PyList_GetItem(plane_restraint_list, i_plane);
             if (PyObject_Length(plane_restraint) == 3) {
                std::vector<std::string> atoms;
                PyObject *plane_id_py = PyList_GetItem(plane_restraint, 0);
                PyObject *esd_py      = PyList_GetItem(plane_restraint, 2);
                PyObject *py_atoms_py = PyList_GetItem(plane_restraint, 1);

                bool atoms_pass = 1;
                if (PyList_Check(py_atoms_py)) {
                   int n_atoms = PyObject_Length(py_atoms_py);
                   for (int iat=0; iat<n_atoms; iat++) {
                      PyObject *at_py = PyList_GetItem(py_atoms_py, iat);
                      if (PyUnicode_Check(at_py)) {
                         atoms.push_back(PyBytes_AS_STRING(PyUnicode_AsUTF8String(at_py)));
                      } else {
                         atoms_pass = 0;
                      }
                   }
                   if (atoms_pass) {
                      if (PyUnicode_Check(plane_id_py)) {
                         if (PyFloat_Check(esd_py)) {
                            std::string plane_id = PyBytes_AS_STRING(PyUnicode_AsUTF8String(plane_id_py));
                            float esd = PyFloat_AsDouble(esd_py);
                            if (atoms.size() > 0) {
                               coot::dict_plane_restraint_t rest(plane_id, atoms[0], esd);
                               for (unsigned int i=1; i<atoms.size(); i++) {
                                  rest.push_back_atom(atoms[i], esd);
                               }
                               plane_restraints.push_back(rest);
                               std::cout << "plane restraint: " << rest <<std::endl;
                            }
                         }
                      }
                   }
                }
             }
          }
       }
    }
      }

      coot::dictionary_residue_restraints_t monomer_restraints(monomer_type, 1);
      monomer_restraints.bond_restraint    = bond_restraints;
      monomer_restraints.angle_restraint   = angle_restraints;
      monomer_restraints.torsion_restraint = torsion_restraints;
      monomer_restraints.chiral_restraint  = chiral_restraints;
      monomer_restraints.plane_restraint   = plane_restraints;
      monomer_restraints.residue_info      = residue_info;
      monomer_restraints.atom_info         = atoms;

      // the imol is not specified, we want to replace the restraints that
      // can be used by any molecule
      int imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
      graphics_info_t g;
      bool s = g.Geom_p()->replace_monomer_restraints(monomer_type, imol_enc, monomer_restraints);
      if (s)
         retval = Py_True;

   }
   if (PyBool_Check(retval)) {
     Py_INCREF(retval);
   }
   return retval;
}
#endif // USE_PYTHON



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
      mg_stream << "# created by Coot " << VERSION << "\n";

      // View:
      graphics_info_t g;
      coot::Cartesian rc = g.RotationCentre();
      mg_stream << "view = View (\n";
      mg_stream << "    centre_xyz = [";
      mg_stream << -rc.x() << ", " << -rc.y() << ", " << -rc.z() << "],\n";
      mg_stream << "    radius = " << 0.75*graphics_info_t::zoom << ",\n";
      //       mg_stream << "    orientation = [ " << g.quat[0] << ", "
      //                 << g.quat[1] << ", " << g.quat[2] << ", " << g.quat[3] << "]\n";
      // Stuart corrects the orientation specification:
      // mg_stream << "    orientation = [ " << -g.quat[3] << ", "
      // << g.quat[0] << ", " << g.quat[1] << ", " << g.quat[2] << "]\n";
      //       mg_stream << ")\n";

      // Parameters (maybe further down?)
      // GUI Parameters (for bg colour e.g.)
      mg_stream << "\n";
      mg_stream << "ParamsManager (\n";
      mg_stream << "   name = 'gui_params',\n";
      mg_stream << "   background_colour = [";
      mg_stream << (int)graphics_info_t::background_colour[0] * 255 << ", ";
      mg_stream << (int)graphics_info_t::background_colour[1] * 255 << ", ";
      mg_stream << (int)graphics_info_t::background_colour[2] * 255 << "]\n";
      mg_stream << ")\n";
      mg_stream << "\n";
      // Display Params (not working I think; default bond width for all)
      mg_stream << "\n";
      mg_stream << "ParamsManager (\n";
      mg_stream << "   name = 'model_drawing_style',\n";
      mg_stream << "   bond_width = ";
      mg_stream << graphics_info_t::default_bond_width << "         )\n";

      // Molecules:
      std::ostringstream map_colour_stream;
      for (int imol=0; imol<graphics_info_t::n_molecules(); imol++) {
         if (is_valid_model_molecule(imol)) {
            mg_stream << "MolData (\n";
            mg_stream << "   filename = [\n";
            mg_stream << "   'FULLPATH',\n";
            mg_stream << "   " << single_quote(coot::util::absolutise_file_name(graphics_info_t::molecules[imol].name_)) << ",\n";
            mg_stream << "   " << single_quote(coot::util::absolutise_file_name(graphics_info_t::molecules[imol].name_)) << "])\n";
            mg_stream << "\n";
            mg_stream << "MolDisp (\n";
            mg_stream << "    select    = 'all',\n";
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
            float lev       = g.molecules[imol].contour_level;
            float r         = g.box_radius_xray;
            // first create a colour
            map_colour_stream << "   coot_mapcolour" << imol << " = ["<<
               graphics_info_t::molecules[imol].map_colour.red   << ", " <<
               graphics_info_t::molecules[imol].map_colour.green << ", " <<
               graphics_info_t::molecules[imol].map_colour.blue  << ", " <<
               "1.0],\n";
            //colour_definitions.push_back(map_colour);
            // int use_weights_flag = g.molecules[imol].save_use_weights;
            std::string name = single_quote(coot::util::absolutise_file_name(graphics_info_t::molecules[imol].save_mtz_file_name));
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
            mg_stream << "    colour = 'coot_mapcolour" << imol << "',\n";
            mg_stream << "    clip_mode = 'OFF')\n";
            mg_stream << "\n";
            if (map_is_difference_map(imol)) {
                   // make a second contour level
                   mg_stream << "MapDisp (\n";
                   mg_stream << "    contour_level = -" << lev << ",\n";
                   mg_stream << "    surface_style = 'chickenwire',\n";
                   mg_stream << "    radius = " << r << ",\n";
                   mg_stream << "    colour = 'coot_mapcolour" << imol << "_2',\n";
                   mg_stream << "    clip_mode = 'OFF')\n";
                   mg_stream << "\n";
                   // and a color
                   map_colour_stream << "   coot_mapcolour" << imol << "_2 = ["<<
                      graphics_info_t::molecules[imol].map_colour.red   << ", " <<
                      graphics_info_t::molecules[imol].map_colour.green << ", " <<
                      graphics_info_t::molecules[imol].map_colour.blue  << ", " <<
                      "1.0],\n";
                }
         }
      }
      // now define map colours
      mg_stream << "Colours (\n";
      mg_stream << map_colour_stream.str();
      mg_stream << ")\n";
      mg_stream << "\n";

   }
}


/*  ----------------------------------------------------------------------- */
/*                  Restratints                                             */
/*  ----------------------------------------------------------------------- */

void write_restraints_cif_dictionary(std::string monomer_type, std::string file_name) {

   int imol = 0;
   graphics_info_t g;
   std::pair<short int, coot::dictionary_residue_restraints_t> r =
      g.Geom_p()->get_monomer_restraints(monomer_type, imol);
   if (!r.first) {
      std::string s =  "Failed to find ";
      s += monomer_type;
      s += " in dictionary";
      add_status_bar_text(s.c_str());
      std::cout << s << std::endl;
   } else {
      r.second.write_cif(file_name);
   }
}

/*  ----------------------------------------------------------------------- */
/*                  PKGDATDDIR                                              */
/*  ----------------------------------------------------------------------- */
#ifdef USE_PYTHON
PyObject *get_pkgdatadir_py() {

   // std::string pkgdatadir = PKGDATADIR;
   std::string pkgdatadir = coot::package_data_dir();
   return myPyString_FromString(pkgdatadir.c_str());
}
#endif

#ifdef USE_GUILE
SCM get_pkgdatadir_scm() {

   // std::string pkgdatadir = PKGDATADIR;
   std::string pkgdatadir = coot::package_data_dir();
   return scm_from_locale_string(pkgdatadir.c_str());
}
#endif


/* refmac version */
/* returns:
   1 if 5.4 or newer
   2 if 5.5 or newer */

int refmac_runs_with_nolabels() {

  int ret = 0;


#ifdef USE_GUILE
  SCM refmac_version = safe_scheme_command("(get-refmac-version)");
  if (scm_is_true(scm_list_p(refmac_version))) {
     int major = scm_to_int(scm_list_ref(refmac_version, scm_from_int(0)));
     int minor = scm_to_int(scm_list_ref(refmac_version, scm_from_int(1)));
     if ((major == 5 && minor >= 4) || (major > 5)) {
        ret = 1;
        if (minor >= 5 || major > 5) {
          // for SAD and TWIN
          ret = 2;
        }
     }
  }

#else
#ifdef USE_PYTHON
  PyObject *refmac_version;
  refmac_version = safe_python_command_with_return("get_refmac_version()");
  if (refmac_version) {
     int major = PyLong_AsLong(PyList_GetItem(refmac_version, 0));
     int minor = PyLong_AsLong(PyList_GetItem(refmac_version, 1));
     if ((major == 5 && minor >= 4) || (major > 5)) {
        ret = 1;
        if (minor >= 5 || major > 5) {
          // for SAD and TWIN
          ret = 2;
        }
     }
  }
  Py_XDECREF(refmac_version);
#endif
#endif


  return ret;
}


/*! \brief allow the user to not add ccp4i directories to the file choosers

use state=0 to turn it off */
void set_add_ccp4i_projects_to_file_dialogs(short int state) {

   graphics_info_t::add_ccp4i_projects_to_optionmenu_flag = state;

}


#ifdef USE_GUILE
SCM ccp4i_projects_scm() {
   SCM r = SCM_EOL;
   std::string ccp4_defs_file_name = graphics_info_t::ccp4_defs_file_name();
   std::vector<std::pair<std::string, std::string> > project_pairs =
      parse_ccp4i_defs(ccp4_defs_file_name);
   for (unsigned int i=0; i<project_pairs.size(); i++) {
      SCM p = SCM_EOL;
      p = scm_cons(scm_from_locale_string(project_pairs[i].second.c_str()), p);
      p = scm_cons(scm_from_locale_string(project_pairs[i].first.c_str()),  p);
      r = scm_cons(p, r);
   }
   r = scm_reverse(r);
   return r;
}
#endif // USE_GUILE

#ifdef USE_PYTHON
PyObject *ccp4i_projects_py() {
   PyObject *r = PyList_New(0);
   std::string ccp4_defs_file_name = graphics_info_t::ccp4_defs_file_name();
   std::vector<std::pair<std::string, std::string> > project_pairs =
      parse_ccp4i_defs(ccp4_defs_file_name);
   for (unsigned int i=0; i<project_pairs.size(); i++) {
      PyObject *p = PyList_New(2);
      PyList_SetItem(p, 0, myPyString_FromString(project_pairs[i].first.c_str()));
      PyList_SetItem(p, 1, myPyString_FromString(project_pairs[i].second.c_str()));
      PyList_Append(r, p);
      Py_XDECREF(p);
   }
   return r;
}
#endif // USE_PYTHON



#ifdef USE_GUILE
SCM remarks_scm(int imol) {

   SCM r = SCM_EOL;

   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      mmdb::TitleContainer *tc_p = mol->GetRemarks();
      int l = tc_p->Length();
      for (int i=0; i<l; i++) {
         mmdb::Remark *cr = static_cast<mmdb::Remark *> (tc_p->GetContainerClass(i));
         SCM a_scm = scm_from_int(cr->remarkNum);
         SCM b_scm = scm_from_locale_string(cr->remark);

         // #define SCM_LIST2(e0, e1) scm_cons2 ((e0), (e1), SCM_EOL)
         // SCM l2 = SCM_LIST2(a_scm, b_scm);
         SCM l2 = scm_cons2(a_scm, b_scm, SCM_EOL);
         r = scm_cons(l2, r);
      }
      r = scm_reverse(r); // undo schemey backwardsness
   }
   return r;
}

#endif


#ifdef USE_PYTHON
PyObject *remarks_py(int imol) {

   PyObject *o = Py_False;

   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      mmdb::TitleContainer *tc_p = mol->GetRemarks();
      int n_records = tc_p->Length();
      o = PyList_New(n_records);
      for (int i=0; i<n_records; i++) {
         mmdb::Remark *cr = static_cast<mmdb::Remark *> (tc_p->GetContainerClass(i));
         PyObject *l = PyList_New(2);
         PyList_SetItem(l, 0, PyLong_FromLong(cr->remarkNum));
         PyList_SetItem(l, 1, myPyString_FromString(cr->remark));
         PyList_SetItem(o, i, l);
      }
   }
   if (PyBool_Check(o))
      Py_INCREF(o);
   return o;
}

#endif


#ifdef USE_GUILE
SCM residue_centre_scm(int imol, const char *chain_id, int resno, const char *ins_code) {

   SCM r = SCM_BOOL_F;

   if (is_valid_model_molecule(imol)) {
      std::pair<bool, clipper::Coord_orth> rr =
         graphics_info_t::molecules[imol].residue_centre(chain_id, resno, ins_code);
      if (rr.first) {
         r = scm_list_3(scm_from_double(rr.second.x()),
                        scm_from_double(rr.second.y()),
                        scm_from_double(rr.second.z()));
      }
   }
   return r;
}

#endif


#ifdef USE_PYTHON
PyObject *residue_centre_py(int imol, const char *chain_id, int resno, const char *ins_code) {

   PyObject *r = Py_False;

   if (is_valid_model_molecule(imol)) {
      std::pair<bool, clipper::Coord_orth> rr =
         graphics_info_t::molecules[imol].residue_centre(chain_id, resno, ins_code);
      if (rr.first) {
         r = PyList_New(3);
         PyList_SetItem(r, 0, PyFloat_FromDouble(rr.second.x()));
         PyList_SetItem(r, 1, PyFloat_FromDouble(rr.second.y()));
         PyList_SetItem(r, 2, PyFloat_FromDouble(rr.second.z()));
      }
   }
   if (PyBool_Check(r)) {
      Py_INCREF(r);
   }
   return r;
}

#endif

#ifdef USE_PYTHON
// the expanded form of this is in c-interface.h
PyObject *residue_centre_from_spec_py(int imol,
                                      PyObject *spec_py) {

   PyObject *r = Py_False;
   std::pair<bool, coot::residue_spec_t> spec = make_residue_spec_py(spec_py);

   if (spec.first) {
      if (is_valid_model_molecule(imol)) {
         std::pair<bool, clipper::Coord_orth> rr =
            graphics_info_t::molecules[imol].residue_centre(spec.second);
         if (rr.first) {
            r = PyList_New(3);
            PyList_SetItem(r, 0, PyFloat_FromDouble(rr.second.x()));
            PyList_SetItem(r, 1, PyFloat_FromDouble(rr.second.y()));
            PyList_SetItem(r, 2, PyFloat_FromDouble(rr.second.z()));
         }
      }
   }
   if (PyBool_Check(r))
      Py_INCREF(r);

   return r;
}
#endif // USE_PYTHON



#ifdef USE_GUILE
SCM link_info_scm(int imol) {

   SCM r = SCM_EOL;
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      if (mol) {

         for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {

            mmdb::Model *model_p = mol->GetModel(imod);
            int n_links = model_p->GetNumberOfLinks();
            if (n_links > 0) {
               for (int i_link=1; i_link<=n_links; i_link++) {
                  mmdb::PLink link = model_p->GetLink(i_link);

                  std::pair<coot::atom_spec_t, coot::atom_spec_t> atoms = coot::link_atoms(link, model_p);
                  SCM l = scm_list_3(scm_from_int(imod),
                                     atom_spec_to_scm(atoms.first),
                                     atom_spec_to_scm(atoms.second));
                  r = scm_cons(l,r);
               }
            }
         }
         r = scm_reverse(r);
      }
   }
   return r;
}

#endif




#ifdef USE_PYTHON
PyObject *link_info_py(int imol) {

   PyObject *r = PyList_New(0);
   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
      if (mol) {
         for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
            mmdb::Model *model_p = mol->GetModel(imod);
            int n_links = model_p->GetNumberOfLinks();
            if (n_links > 0) {
               for (int i_link=1; i_link<=n_links; i_link++) {
                  mmdb::PLink link = model_p->GetLink(i_link);
                  std::pair<coot::atom_spec_t, coot::atom_spec_t> atoms = coot::link_atoms(link, model_p);
                  PyObject *l = PyList_New(3);
                  PyList_SetItem(l, 0, PyLong_FromLong(imod));
                  PyList_SetItem(l, 1, atom_spec_to_py(atoms.first));
                  PyList_SetItem(l, 2, atom_spec_to_py(atoms.second));
                  PyList_Append(r, l);
               }
            }
         }
      }
   }
   return r;
}

#endif
