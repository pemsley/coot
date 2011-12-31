/* src/c-interface.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Copyright 2008, 2009 by The University of Oxford
 * Author: Paul Emsley, Bernhard Lohkamp
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 3 of the
 * License, or (at your option) any later version.
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

// $Id: c-interface.cc 1458 2007-01-26 20:20:18Z emsley $
// $LastChangedDate: 2007-01-26 20:20:18 +0000 (Fri, 26 Jan 2007) $
// $Rev: 1458 $
 
// Load the head if it hasn't been included.
#ifdef USE_PYTHON
#ifndef PYTHONH
#define PYTHONH
#include <Python.h>
#endif
#endif

#include <stdlib.h>
#include <iostream>
#include <fstream>

#if !defined(_MSC_VER)
#include <glob.h> // for globbing.  Needed here?
#endif

#ifdef USE_GUILE
#include <guile/gh.h>
#include "c-interface-scm.hh"
#include "guile-fixups.h"
#endif // USE_GUILE

#ifdef USE_PYTHON
#include "c-interface-python.hh"
#endif // USE_PYTHON

#if defined (WINDOWS_MINGW)
#ifdef DATADIR
#undef DATADIR
#endif // DATADIR
#endif /* MINGW */
#include "sleep-fixups.h"

// Here we used to define GTK_ENABLE_BROKEN if defined(WINDOWS_MINGW)
// Now we don't want to enable broken stuff.  That is not the way.

#define HAVE_CIF  // will become unnessary at some stage.

#include <sys/types.h> // for stating
#include <sys/stat.h>
#if !defined _MSC_VER
#include <unistd.h>
#else
#define S_IRUSR S_IREAD
#define S_IWUSR S_IWRITE
#define S_IXUSR S_IEXEC
#define S_ISDIR(m) (((m) & S_IFMT) == S_IFDIR)
#define S_ISREG(m) (((m) & S_IFMT) == S_IFREG)
#define snprintf _snprintf
#include <windows.h>
#include <direct.h>
#endif // _MSC_VER

#include "clipper/ccp4/ccp4_map_io.h"
 
#include "globjects.h" //includes gtk/gtk.h

#include "callbacks.h"
#include "interface.h" // now that we are moving callback
		       // functionality to the file, we need this
		       // header since some of the callbacks call
		       // fuctions built by glade.

#include <vector>
#include <string>

#include <mmdb/mmdb_manager.h>
#include "mmdb-extras.h"
#include "mmdb.h"
#include "mmdb-crystal.h"

#include "Cartesian.h"
#include "Bond_lines.h"
#include "coot-utils.hh"
#include "coot-map-utils.hh"
#include "coot-database.hh"
#include "coot-fileselections.h"

// #include "xmap-interface.h"
#include "graphics-info.h"

#include "BuildCas.h"

#include "trackball.h" // adding exportable rotate interface


#include "c-interface.h"
#include "cc-interface.hh"
#include "c-interface-ligands.hh"

#include "nsv.hh"

#include "testing.hh"

#include "positioned-widgets.h"

// moving column_label selection to c-interface from mtz bits.
#include "cmtz-interface.hh"
// #include "mtz-bits.h" stuff from here moved to cmtz-interface

// Including python needs to come after graphics-info.h, because
// something in Python.h (2.4 - chihiro) is redefining FF1 (in
// ssm_superpose.h) to be 0x00004000 (Grrr).
//
#ifdef USE_PYTHON
#include "Python.h"
#endif // USE_PYTHON

#ifdef HAVE_GOOCANVAS
#include "goocanvas.h"
#include "wmolecule.hh"
#include "goograph.hh"
#endif

#include "simple-restraint.hh"  // for multi-residue torsion map fitting.


int test_function(int i, int j) {

   graphics_info_t g;

   // Is this the function you are really looking for (these days)?

   if (1) {
   } 
   

   if (0) {

      if (is_valid_model_molecule(i)) { 
	 if (is_valid_map_molecule(j)) { 
	    const clipper::Xmap<float> &xmap = g.molecules[j].xmap_list[0];
	    CMMDBManager *mol = g.molecules[i].atom_sel.mol;
	    std::vector<coot::residue_spec_t> v;
	    v.push_back(coot::residue_spec_t("G", 160, ""));
	    v.push_back(coot::residue_spec_t("G", 847, ""));

	    int n_rounds = 10;
	    for (unsigned int iround=0; iround<n_rounds; iround++) { 
   
	       CMMDBManager *moving_mol = coot::util::create_mmdbmanager_from_residue_specs(v, mol);

	       // do we need to send over the base atom too?  Or just say
	       // that it's the first atom in moving_mol?
	       // 
	       coot::multi_residue_torsion_fit_map(moving_mol, xmap, g.Geom_p());

	       atom_selection_container_t moving_atoms_asc = make_asc(moving_mol);

	       std::pair<CMMDBManager *, int> new_mol =
		  coot::util::create_mmdbmanager_from_mmdbmanager(moving_mol);
	       atom_selection_container_t asc_new = make_asc(new_mol.first);
	       std::string name = "test-" + coot::util::int_to_string(iround);
	       bool shelx_flag = 0;
	       int imol_new = g.create_molecule();
	       g.molecules[imol_new].install_model(imol_new, asc_new, name, 1, shelx_flag);

	       // Don't update - not at the moment at least.
	       // 
	       // g.molecules[i].replace_coords(moving_atoms_asc, 1, 1);
	       
	       delete moving_mol;
	       graphics_draw();
	    }
	 }
      }
   }



   if (0) {
      // coot::atom_spec_t spec_1("A", 41, "", " OE1", "");
      // coot::atom_spec_t spec_2("A", 39, "", " N  ", "");
//       coot::atom_spec_t spec_1("B", 48, "", " OG ", "");
//       coot::atom_spec_t spec_2("A", 48, "", " N  ", "");
//       graphics_info_t::molecules[0].add_animated_ligand_interaction(spec_1, spec_2);
//       coot::atom_spec_t spec_3("B", 47, "", " O  ", "");
//       coot::atom_spec_t spec_4("B", 48, "", " N  ", "");
//       graphics_info_t::molecules[0].add_animated_ligand_interaction(spec_3, spec_4);
//       coot::atom_spec_t spec_5("B", 48, "", " O  ", "");
//       coot::atom_spec_t spec_6("B", 49, "", " N  ", "");
//       graphics_info_t::molecules[0].add_animated_ligand_interaction(spec_5, spec_6);
//       coot::atom_spec_t spec_7("B", 49, "", " O  ", "");
//       coot::atom_spec_t spec_8("B", 50, "", " N  ", "");
//       graphics_info_t::molecules[0].add_animated_ligand_interaction(spec_7, spec_8);
   } 

       
#ifdef HAVE_GOOCANVAS      
   if (0) {
      std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
      graphics_info_t g;
      if (pp.first) {
	 CResidue *residue = g.molecules[pp.second.first].get_residue(pp.second.second);
	 CMMDBManager *mol = g.molecules[pp.second.first].atom_sel.mol;
	 if (residue) {
	    std::pair<bool, coot::dictionary_residue_restraints_t> restraints =
	       g.Geom_p()->get_monomer_restraints(residue->GetResName());
	    lig_build::molfile_molecule_t mm(residue, restraints.second);
	    widgeted_molecule_t wm(mm, mol);
	    topological_equivalence_t top_eq(wm.atoms, wm.bonds);
	 } 
      }
   }
#endif // HAVE_GOOCANVAS
   
   if (0) {

      if (is_valid_model_molecule(0)) {
	 CMMDBManager *mol = graphics_info_t::molecules[0].atom_sel.mol;
	 std::vector<std::string> h;
	 CTitleContainer *tc_p = mol->GetRemarks();
	 int l = tc_p->Length();
	 for (unsigned int i=0; i<l; i++) { 
	    CRemark *cr = static_cast<CRemark *> (tc_p->GetContainerClass(i));
	    std::cout << "container: " << cr->Remark << std::endl;
	 }
      }
   }


   if (0) {
      std::vector<std::pair<std::string, int> > h = 
	 coot::get_prodrg_hybridizations("coot-ccp4/tmp-prodrg-flat.log");
      
   } 

   if (0) {
      // atom_selection_container_t asc = get_atom_selection("double.pdb");
      atom_selection_container_t asc = get_atom_selection("test-frag.pdb", 1);
      coot::dots_representation_info_t dots;
      int sel_hnd = asc.SelectionHandle;
      std::vector<std::pair<CAtom *, float> > v =
	 dots.solvent_exposure(sel_hnd, asc.mol);

   } 

   if (0) {

      std::cout << "sizeof(int): " << sizeof(int) << std::endl;
      
      if (graphics_info_t::use_graphics_interface_flag) { 
	 GtkWidget *w = lookup_widget(graphics_info_t::glarea,
				      "main_window_model_fit_dialog_frame");
	 if (!w) {
	    std::cout << "failed to lookup toolbar" << std::endl;
	 } else {
	    gtk_widget_hide(w);
	 }
      }
   }

   if (0) {
      graphics_info_t::molecules[i].test_function();
   } 

   if (0) {
      GtkWidget *w = wrapped_create_add_additional_representation_gui();
      gtk_widget_show(w);
   } 

   if (0) {
      coot::util::quaternion::test_quaternion();
   }
   

   if (0) {
      graphics_info_t g;
      g.Geom_p()->hydrogens_connect_file("THH", "thh_connect.txt");
   }

   if (0) {

      // GTK2 GTkGLExt code
//       GdkGLContext *glcontext = gtk_widget_get_gl_context(graphics_info_t::glarea);
//       GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable(graphics_info_t::glarea);
//       GdkGLConfig *glconfig = gtk_widget_get_gl_config(graphics_info_t::glarea);
//       Display *dpy = gdk_x11_gl_config_get_xdisplay (glconfig);
      // Bool glXMakeCurrent(Display * dpy,
      //                     GLXDrawable  Drawable,
      //                     GLXContext  Context)
      // gdk_gl_glXMakeContextCurrent(dpy, gldrawable, glcontext);

      // bwah!
      // glXMakeCurrent(dpy, gldrawable, glcontext);

      // another way?
//       GtkWidget *w = graphics_info_t::glarea;
//       GdkGLContext *glcontext = gtk_widget_get_gl_context (w);
//       GdkGLDrawable *gldrawable = gtk_widget_get_gl_drawable (w);
//       int i = gdk_gl_drawable_gl_begin (gldrawable, glcontext);
//       std::cout << "DEBUG gdk_gl_drawable_gl_begin returns state: "
// 		<< i << std::endl;
//       return i;
   } 

   if (0) {
      int imol = i;
      if (is_valid_model_molecule(imol)) { 
	 const coot::residue_spec_t clicked_residue("A", 1);
	 short int is_n_term_addition = 1;
	 CAtom *at = graphics_info_t::molecules[imol].atom_sel.atom_selection[10];
	 CChain *chain_p = at->GetChain();
	 std::pair<bool, std::string> p = 
	    graphics_info_t::molecules[imol].residue_type_next_residue_by_alignment(clicked_residue, chain_p, is_n_term_addition, graphics_info_t::alignment_wgap, graphics_info_t::alignment_wspace);
	 if (p.first == 1) { 
	    std::cout << "next residue: " << p.second << std::endl;
	 } else {
	    std::cout << "no next residue found." << std::endl;
	 }
      }
   } 


   if (0) { 
      GtkWidget *w = wrapped_create_least_squares_dialog();
      gtk_widget_show(w);
   }
      

   if (0) { 
      std::vector<std::string> s;
      s.push_back("");
      s.push_back("123");
      s.push_back("123/456");
      s.push_back("123/456/");
      
      for (unsigned int i=0; i<s.size(); i++) { 
	 std::pair<std::string, std::string> p = coot::util::split_string_on_last_slash(s[i]);
	 std::cout << "For string :" << s[i] << ": split is :"
		   << p.first << ": :" << p.second << ":" << std::endl; 
      }

      std::string t = "/my/thing/int.mtz data/crystal/FWT data/crystal/PHWT";
      std::vector<std::string> v = coot::util::split_string(t, " ");

      std::cout << "splitting :" << t << ": on " << " " << std::endl;
      for (unsigned int i=0; i<v.size(); i++) {
	 std::cout << "split " << i << " :" << v[i] << ":\n";
      }
   }
   return 0;
}


#ifdef __cplusplus
#ifdef USE_GUILE
SCM test_function_scm(SCM i_scm, SCM j_scm) {

   graphics_info_t g;
   SCM r = SCM_BOOL_F;

   if (1) {
      coot::goograph g;
      g.show_dialog();
   } 

   if (0) {
      int i = scm_to_int(i_scm); // map molecule
      int j = scm_to_int(j_scm);

      // was_found, imol, atom_spec
      std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = active_atom_spec();
      if (active_atom.first) {

	 int imol = active_atom.second.first;
	 coot::atom_spec_t spec = active_atom.second.second;
	 if (! is_valid_map_molecule(i)) {
	    std::cout << "Not valid map " << i << std::endl;
	 } else { 
	    std::vector<coot::residue_spec_t> v;
	    v.push_back(coot::residue_spec_t(spec));
	    int n_rounds = 10;
	    CMMDBManager *mol = g.molecules[imol].atom_sel.mol;
	    const clipper::Xmap<float> &xmap = g.molecules[j].xmap_list[0];
	    for (unsigned int iround=0; iround<n_rounds; iround++) {
	       std::cout << "round " << iround << std::endl;
	       CMMDBManager *moving_mol = coot::util::create_mmdbmanager_from_residue_specs(v, mol);
	       
	       coot::multi_residue_torsion_fit_map(moving_mol, xmap, g.Geom_p());
	       atom_selection_container_t moving_atoms_asc = make_asc(moving_mol);
	       std::pair<CMMDBManager *, int> new_mol =
		  coot::util::create_mmdbmanager_from_mmdbmanager(moving_mol);
	       atom_selection_container_t asc_new = make_asc(new_mol.first);
	       std::string name = "test-" + coot::util::int_to_string(iround);
	       bool shelx_flag = 0;
	       int imol_new = g.create_molecule();
	       g.molecules[imol_new].install_model(imol_new, asc_new, name, 1, shelx_flag);
	       add_linked_residue(imol_new,
				  active_atom.second.second.chain.c_str(),
				  active_atom.second.second.resno,
				  active_atom.second.second.insertion_code.c_str(),
				  "NAG", "ASN-NAG");
	       
	       delete moving_mol;
	       graphics_draw();
	    }
	 } 
      }
   }
   return r;
}

#endif

#ifdef USE_PYTHON
PyObject *test_function_py(PyObject *i_py, PyObject *j_py) {

   graphics_info_t g;
   PyObject *r = Py_False;

   if (1) {
     int i = PyInt_AsLong(i_py); // map molecule
     int j = PyInt_AsLong(j_py);

     // was_found, imol, atom_spec
     std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = active_atom_spec();
     if (active_atom.first) {

       int imol = active_atom.second.first;
       coot::atom_spec_t spec = active_atom.second.second;
       if (! is_valid_map_molecule(i)) {
         std::cout << "Not valid map " << i << std::endl;
       } else { 
         std::vector<coot::residue_spec_t> v;
         v.push_back(coot::residue_spec_t(spec));
         int n_rounds = 10;
         CMMDBManager *mol = g.molecules[imol].atom_sel.mol;
         const clipper::Xmap<float> &xmap = g.molecules[j].xmap_list[0];
         for (unsigned int iround=0; iround<n_rounds; iround++) {
	       std::cout << "round " << iround << std::endl;
	       CMMDBManager *moving_mol = coot::util::create_mmdbmanager_from_residue_specs(v, mol);
	       
	       coot::multi_residue_torsion_fit_map(moving_mol, xmap, g.Geom_p());
	       atom_selection_container_t moving_atoms_asc = make_asc(moving_mol);
	       std::pair<CMMDBManager *, int> new_mol =
             coot::util::create_mmdbmanager_from_mmdbmanager(moving_mol);
	       atom_selection_container_t asc_new = make_asc(new_mol.first);
	       std::string name = "test-" + coot::util::int_to_string(iround);
	       bool shelx_flag = 0;
	       int imol_new = g.create_molecule();
	       g.molecules[imol_new].install_model(imol_new, asc_new, name, 1, shelx_flag);
	       add_linked_residue(imol_new,
                              active_atom.second.second.chain.c_str(),
                              active_atom.second.second.resno,
                              active_atom.second.second.insertion_code.c_str(),
                              "NAG", "ASN-NAG");
	       
	       delete moving_mol;
	       graphics_draw();
         }
       } 
     }
   }
   if (PyBool_Check(r)) {
     Py_INCREF(r);
   }
   return r;
}

#endif // PYTHON
#endif //C_PLUS_PLUS


/* glyco tools test  */
void glyco_tree_test() {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();

   if (pp.first) {
      int imol = pp.second.first;
      graphics_info_t g;
      CResidue *residue_p = g.molecules[imol].get_residue(pp.second.second);
      CMMDBManager *mol = g.molecules[imol].atom_sel.mol;

      std::vector<std::string> types_with_no_dictionary =
	 g.molecules[imol].no_dictionary_for_residue_type_as_yet(*g.Geom_p());
      std::cout << "glyco-test found " << types_with_no_dictionary.size() << " types with no dictionary"
		<< std::endl;
      for (unsigned int i=0; i<types_with_no_dictionary.size(); i++) {
	 std::cout << "trying to dynamic add: " << types_with_no_dictionary[i] << std::endl;
	 g.Geom_p()->try_dynamic_add(types_with_no_dictionary[i], 41);
      }
      
      coot::glyco_tree_t(residue_p, mol, g.Geom_p());
   } 

}


