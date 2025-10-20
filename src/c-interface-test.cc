/* src/c-interface-test.cc
 *
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Copyright 2008, 2009 by The University of Oxford
 * Copyright 2014 by Medical Research Council
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
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

// $Id: c-interface.cc 1458 2007-01-26 20:20:18Z emsley $
// $LastChangedDate: 2007-01-26 20:20:18 +0000 (Fri, 26 Jan 2007) $
// $Rev: 1458 $

// Load the head if it hasn't been included.

#include "gtk/gtk.h"
#include "lidia-core/svg-container.hh"
#include <filesystem>
#include <unordered_map>
#ifdef USE_PYTHON
#ifndef PYTHONH
#define PYTHONH
#include <Python.h>
#endif
#endif

#include "compat/coot-sysdep.h"

#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>

#if !defined(_MSC_VER)
#include <glob.h> // for globbing.  Needed here?
#endif

#ifdef USE_GUILE
#include <libguile.h>
#include "c-interface-scm.hh"
#include "guile-fixups.h"
#endif // USE_GUILE

#ifdef USE_PYTHON
#include "c-interface-python.hh"
#endif // USE_PYTHON

#include "compat/sleep-fixups.h"

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

#include <clipper/ccp4/ccp4_map_io.h>

#include "globjects.h" //includes gtk/gtk.h

#include <vector>
#include <string>

#include <mmdb2/mmdb_manager.h>
#include "coords/mmdb-extras.hh"
#include "coords/mmdb.hh"
#include "coords/mmdb-crystal.hh"
#include "coords/Cartesian.hh"
#include "coords/Bond_lines.hh"

#include "utils/coot-utils.hh"
#include "coot-utils/read-sm-cif.hh"
#include "coot-utils/coot-map-utils.hh"
#include "coot-utils/glyco-tree.hh"
#include "coot-database.hh"
#include "coot-fileselections.h"

// #include "xmap-interface.h"
#include "graphics-info.h"

#include "skeleton/BuildCas.h"

#include "c-interface.h"
#include "c-interface-generic-objects.h"
#include "c-interface-gtk-widgets.h"
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
#include <goocanvas.h>
#include "lbg/wmolecule.hh"
#include "goograph/goograph.hh"
#endif

#include "curl-utils.hh"

#include "ideal/simple-restraint.hh"  // for multi-residue torsion map fitting.
#include "ideal/torsion-bonds.hh"     // for multi-residue torsion map fitting.

#include "cc-interface-network.hh"
#include "c-interface-ligands-swig.hh"

#include "coot-utils/emma.hh"

#include "c-interface-widgets.hh" // for wrapped_create_generic_objects_dialog();

#ifdef USE_ASSIMP
#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing flags
#endif

#include "read-molecule.hh" // now with std::string args

#include "draw-2.hh" // for glarea_tick_func

int test_function(int i, int j) {

   graphics_info_t g;

   // Is this the function you are really looking for (these days)?

   if (true) {
      if (is_valid_model_molecule(i)) {
         if (is_valid_map_molecule(j)) {
            const clipper::Xmap<float> &xmap(g.molecules[j].xmap);
            float scale_factor = 4;
            float offset = 3;
            g.molecules[i].recolour_ribbon_by_map(xmap, scale_factor, offset);
            graphics_draw();
         }
      }
   }

   if (false) {
      g.setup_draw_for_happy_face_residue_markers();
   }

   if (false) {

      std::cout << "Hydrogen bonds mesh test" << std::endl;
      Material material;
      material.shininess = 10.0;
      material.specular_strength = 0.02;
      Mesh mesh("Test cyclinders");
      glm::vec3 p1(42.08, 9.67, 14.42);
      glm::vec3 p2(40.59, 5.68, 13.24);
      glm::vec3 p3(44.88, 12.95, 8.76);
      glm::vec3 p4(46.13, 10.59, 9.97);
      graphics_info_t::hydrogen_bonds_atom_position_pairs.push_back(std::pair<glm::vec3, glm::vec3>(p1, p2));
      graphics_info_t::hydrogen_bonds_atom_position_pairs.push_back(std::pair<glm::vec3, glm::vec3>(p3, p4));
      gtk_gl_area_attach_buffers(GTK_GL_AREA(graphics_info_t::glareas[0]));
      graphics_info_t::mesh_for_hydrogen_bonds = mesh;
      Shader &shader = graphics_info_t::shader_for_instanced_objects;
      graphics_info_t::mesh_for_hydrogen_bonds.setup_hydrogen_bond_cyclinders(&shader, material);
      graphics_info_t::do_tick_hydrogen_bonds_mesh = true;
      int new_tick_id = gtk_widget_add_tick_callback(graphics_info_t::glareas[0], g.glarea_tick_func, 0, 0);
   }

#ifdef USE_ASSIMP
   if (false) {
      std::string file_name = "cube.obj";
      file_name = "cessna.obj";

      Assimp::Importer importer;

      // And have it read the given file with some example postprocessing
      // Usually - if speed is not the most important aspect for you - you'll
      // probably to request more postprocessing than we do in this example.
      const aiScene* scene = importer.ReadFile(file_name,
                                                aiProcess_CalcTangentSpace       |
                                                aiProcess_Triangulate            |
                                                aiProcess_JoinIdenticalVertices  |
                                                aiProcess_SortByPType);

      // If the import failed, report it
      if( !scene) {
         std::cout << "Error in read of " <<  file_name << " " << importer.GetErrorString() << std::endl;
      } else {
         std::cout << "------------ scene read OK from " << file_name << std::endl;
      }
   }
#endif

   if (0) {

      if (is_valid_model_molecule(i)) {
         if (is_valid_map_molecule(j)) {
            const clipper::Xmap<float> &xmap = g.molecules[j].xmap;
            mmdb::Manager *mol = g.molecules[i].atom_sel.mol;
            int imol = 0; // dummy
            std::vector<coot::residue_spec_t> v;
            v.push_back(coot::residue_spec_t("G", 160, ""));
            v.push_back(coot::residue_spec_t("G", 847, ""));

            unsigned int n_rounds = 10;
            for (unsigned int iround=0; iround<n_rounds; iround++) {

               mmdb::Manager *moving_mol = coot::util::create_mmdbmanager_from_residue_specs(v, mol);

               std::vector<std::pair<bool, clipper::Coord_orth> > avoid_these_atoms;

               // do we need to send over the base atom too?  Or just say
               // that it's the first atom in moving_mol?
               //
               coot::multi_residue_torsion_fit_map(imol, moving_mol, xmap, avoid_these_atoms, 400, g.Geom_p());

               atom_selection_container_t moving_atoms_asc = make_asc(moving_mol);

               std::pair<mmdb::Manager *, int> new_mol =
                  coot::util::create_mmdbmanager_from_mmdbmanager(moving_mol);
               atom_selection_container_t asc_new = make_asc(new_mol.first);
               std::string name = "test-" + coot::util::int_to_string(iround);
               bool shelx_flag = 0;
               int imol_new = g.create_molecule();
               g.molecules[imol_new].install_model(imol_new, asc_new, g.Geom_p(), name, 1, shelx_flag);

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
         coot::residue_spec_t rs(coot::atom_spec_t(pp.second.second));
	 mmdb::Residue *residue = g.molecules[pp.second.first].get_residue(rs);
	 mmdb::Manager *mol = g.molecules[pp.second.first].atom_sel.mol;
	 if (residue) {
	    int imol = 0;
	    std::pair<bool, coot::dictionary_residue_restraints_t> restraints =
	       g.Geom_p()->get_monomer_restraints(residue->GetResName(), imol);
	    lig_build::molfile_molecule_t mm(residue, restraints.second);
	    widgeted_molecule_t wm(mm, mol);
	    topological_equivalence_t top_eq(wm.atoms, wm.bonds);
	 }
      }
   }
#endif // HAVE_GOOCANVAS

   if (0) {

      if (is_valid_model_molecule(0)) {
	 mmdb::Manager *mol = graphics_info_t::molecules[0].atom_sel.mol;
	 std::vector<std::string> h;
	 mmdb::TitleContainer *tc_p = mol->GetRemarks();
	 unsigned int l = tc_p->Length();
	 for (unsigned int i=0; i<l; i++) {
	    mmdb::Remark *cr = static_cast<mmdb::Remark *> (tc_p->GetContainerClass(i));
	    std::cout << "container: " << cr->remark << std::endl;
	 }
      }
   }


   if (0) {
      std::vector<std::pair<std::string, int> > h =
	 coot::get_prodrg_hybridizations("coot-ccp4/tmp-prodrg-flat.log");

   }

   if (0) {
      // atom_selection_container_t asc = get_atom_selection("double.pdb");
      atom_selection_container_t asc = get_atom_selection("test-frag.pdb", false, true, false);
      pli::dots_representation_info_t dots;
      int sel_hnd = asc.SelectionHandle;
      std::vector<std::pair<mmdb::Atom *, float> > v = dots.solvent_exposure(sel_hnd, asc.mol);

   }

   if (0) {
      graphics_info_t::molecules[i].test_function();
   }

   if (0) {
      GtkWidget *w = wrapped_create_add_additional_representation_gui();
      gtk_widget_set_visible(w, TRUE);
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
	 mmdb::Atom *at = graphics_info_t::molecules[imol].atom_sel.atom_selection[10];
	 mmdb::Chain *chain_p = at->GetChain();
	 std::pair<bool, std::string> p =
	    graphics_info_t::molecules[imol].residue_type_next_residue_by_alignment(clicked_residue, chain_p, is_n_term_addition, graphics_info_t::alignment_wgap, graphics_info_t::alignment_wspace);
	 if (p.first == 1) {
	    std::cout << "next residue: " << p.second << std::endl;
	 } else {
	    std::cout << "no next residue found." << std::endl;
	 }
      }
   }


   if (false) {
      GtkWidget *w = wrapped_create_least_squares_dialog();
      gtk_widget_set_visible(w, TRUE);
   }


   if (false) {
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


#include "analysis/kolmogorov.hh"
#include "analysis/stats.hh"

#ifdef USE_MOLECULES_TO_TRIANGLES

// Martin's MoleculeToTriangles
//
//
// #include <MoleculesToTriangles/CXXClasses/RendererGL.h>
#include <MoleculesToTriangles/CXXClasses/Light.h>
#include <MoleculesToTriangles/CXXClasses/Camera.h>
// #include <CXXClasses/CameraPort.h>
#include <MoleculesToTriangles/CXXClasses/SceneSetup.h>
#include <MoleculesToTriangles/CXXClasses/ColorScheme.h>
#include <MoleculesToTriangles/CXXClasses/MyMolecule.h>
#include <MoleculesToTriangles/CXXClasses/RepresentationInstance.h>
#include <MoleculesToTriangles/CXXClasses/MolecularRepresentationInstance.h>
#endif // USE_MOLECULES_TO_TRIANGLES

#include "coot-utils/c-beta-deviations.hh"

#include "ligand/richardson-rotamer.hh"

#include "coot-utils/cablam-markup.hh"
#include "coot-utils/pepflip-using-difference-map.hh"

#include "coot-utils/atom-tools.hh"

#include "ligand/libres-tracer.hh"

void get_mol_edit_lock(std::atomic<bool> &mol_edit_lock);
void release_mol_edit_lock(std::atomic<bool> &mol_edit_lock);

void testing_get_mol_edit_lock(std::atomic<bool> &mol_edit_lock) {
   // std::cout << "debug:: test_function_scm() trying to get the lock with mol_edit_lock " << mol_edit_lock << std::endl;
   bool unlocked = false;
   while (! mol_edit_lock.compare_exchange_weak(unlocked, true)) {
      // std::cout << "test_function_scm() failed to get the mol_edit_lock" << std::endl;
      std::this_thread::sleep_for(std::chrono::microseconds(100));
      unlocked = false;
   }
   // std::cout << "debug:: test_function_scm() got the lock" << std::endl;
}

void testing_release_mol_edit_lock(std::atomic<bool> &mol_edit_lock) {
   mol_edit_lock = false;
   // std::cout << "debug:: test_function_scm() released the lock" << std::endl;
};

#ifdef USE_GUILE
SCM test_function_scm(SCM i_scm, SCM j_scm) {

   graphics_info_t g;
   SCM r = SCM_BOOL_F;

   if (true) {

      auto get_alphafold_model_via_uniprot = [] (const std::string &uniprot_id) {

         std::string fn_tail = std::string("AF-") + uniprot_id + std::string("-F1-model_v3.pdb");
         std::string fn = coot::util::append_dir_file("coot-download", fn_tail);
         // make coot-download if needed
         std::string url = std::string("https://alphafold.ebi.ac.uk/files/") + fn_tail;
         bool needs_downloading = true;
         if (coot::file_exists_and_non_empty(fn))
            needs_downloading = false;
         if (needs_downloading) {
            coot_get_url(url.c_str(), fn.c_str());
         }
         int imol = read_pdb(fn);
         return imol;
      };

      int imol = scm_to_int(i_scm);
      if (is_valid_model_molecule(imol)) {
         mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
         if (mol) {
            int imod = 1;
            mmdb::Model *model_p = mol->GetModel(imod);
            if (model_p) {
               int n_chains = model_p->GetNumberOfChains();
               for (int ichain=0; ichain<n_chains; ichain++) {
                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                  int n_refs = chain_p->GetNumberOfDBRefs();
                  std::string chain_id = chain_p->GetChainID();
                  for (int ref_no=0; ref_no<n_refs; ref_no++) {
                     mmdb::DBReference  *ref = chain_p->GetDBRef(ref_no);  // 0..nDBRefs-1
                     std::string db = ref->database;
                     std::string db_accession = ref->dbAccession;
                     std::cout << "Chain " << chain_id << " " << db << " " << db_accession << std::endl;
                     if (db == "UNP") {  // uniprot

                        int imol_af = get_alphafold_model_via_uniprot(db_accession);
                        if (is_valid_model_molecule(imol_af)) {
                           move_molecule_to_screen_centre_internal(imol);
                           superpose_with_chain_selection(imol, imol_af, chain_id.c_str(), "A", 1, 0, 0);
                        }
                     }
                  }
               }
            }
         }
      }
   }

   if (false) {

      std::string hklin_file_name = "coot-download/1gwd_map.mtz";
      std::string f_col_label   = "FWT";
      std::string phi_col_label = "PHWT";
      std::string pir_file_name = "1gwd.pir";
      
      std::cout << "Read mtz file " << hklin_file_name << " " << f_col_label << " " << phi_col_label << std::endl;
      bool use_weights = false;
      bool is_diff_map = false;
      coot::fasta_multi fam;
      fam.read(pir_file_name);
      clipper::Xmap<float> xmap;
      double variation = 0.4; // speed
      unsigned int n_top_spin_pairs = 500; // Use for tracing at most this many spin score pairs (which have been sorted).
      // This and variation affect the run-time (and results?)
      // n_top_spin_pairs = 1000; // was 1000

      unsigned int n_top_fragments = 1000; // was 4000 // The top 1000 fragments at least are all the same trace for no-side-chain lyso test
      float flood_atom_mask_radius = 1.0; // was 0.6 for emdb
      unsigned int n_phi_psi_trials = 40000; // was 5000
      float weight = 8.0f; // calculate this (using rmsd)
      bool with_ncs = false;
      float rmsd_cuffoff = 2.3;

      mmdb::Manager *working_mol = new mmdb::Manager;

      int imol_new = graphics_info_t::create_molecule();
      atom_selection_container_t asc = make_asc(working_mol);
      std::string label = "Building Molecule";
      const std::vector<coot::ghost_molecule_display_t> ghosts;
      bool shelx_flag = false;
      g.molecules[imol_new].install_model_with_ghosts(imol_new, asc, g.Geom_p(), label, 1, ghosts,
                                                      shelx_flag, false, false);
      update_go_to_atom_window_on_new_mol();

      coot::util::map_fill_from_mtz(&xmap, hklin_file_name, f_col_label, phi_col_label, "", use_weights, is_diff_map);
      int imol_new_map = g.create_molecule();
      label = "Map";
      bool is_em_map_flag = false;
      g.molecules[imol_new_map].install_new_map(xmap, label, is_em_map_flag);
      float xmap_rmsd = g.molecules[imol_new_map].map_sigma();
      g.graphics_draw();

      watch_res_tracer_data_t *watch_data_p = new watch_res_tracer_data_t(working_mol, imol_new);
      std::cout << "post-constructor with mol_edit_lock: " << watch_data_p->mol_edit_lock << std::endl;

      // pass geom to this too.
      std::thread t(res_tracer_proc, xmap, xmap_rmsd, fam, variation, n_top_spin_pairs, n_top_fragments, rmsd_cuffoff, flood_atom_mask_radius,
                    weight, n_phi_psi_trials, with_ncs, watch_data_p);

      auto watching_timeout_func = [] (gpointer data) {
                                      watch_res_tracer_data_t *watch_data_p = static_cast<watch_res_tracer_data_t *>(data);
                                      if (false)
                                         std::cout << "debug:: watching_timeout_func runs... finished: " << watch_data_p->finished
                                                   << " lock: " << watch_data_p->mol_edit_lock
                                                   << " update_flag: " << watch_data_p->update_flag << std::endl;
                                      if (watch_data_p->update_flag) {
                                         watch_data_p->update_flag = false;
                                         graphics_info_t g;
                                         testing_get_mol_edit_lock(watch_data_p->mol_edit_lock);
                                         atom_selection_container_t asc_new = make_asc(watch_data_p->working_mol);
                                         g.molecules[watch_data_p->imol_new].atom_sel = asc_new;
                                         g.molecules[watch_data_p->imol_new].make_bonds_type_checked();
                                         testing_release_mol_edit_lock(watch_data_p->mol_edit_lock);
                                         if (watch_data_p->update_count == 1) {
                                            auto rc = g.molecules[watch_data_p->imol_new].centre_of_molecule();
                                            g.setRotationCentreSimple(rc);
                                            update_maps();
                                         }
                                         g.graphics_draw();
                                      }
                                      if (watch_data_p->finished) {
                                         std::cout << "Final update of working_mol..." << std::endl;
                                         testing_get_mol_edit_lock(watch_data_p->mol_edit_lock);
                                         atom_selection_container_t asc_new = make_asc(watch_data_p->working_mol);
                                         graphics_info_t g;
                                         g.molecules[watch_data_p->imol_new].atom_sel = asc_new;
                                         g.molecules[watch_data_p->imol_new].make_bonds_type_checked();
                                         testing_release_mol_edit_lock(watch_data_p->mol_edit_lock);
                                         g.graphics_draw();
                                      }
                                      int return_status = TRUE;
                                      if (watch_data_p->finished)
                                         return_status = FALSE; // don't continue
                                      return return_status;
                                   };

      g_timeout_add(500, watching_timeout_func, watch_data_p);

      t.detach();

   }

   if (false) {

      // test local map sharpening

      // I had intended to use this for interactive local map sharpening.
      // i.e. move a small map fragment to the origin,
      //      sharpen/blur that map
      //      transfer those grid coordinates back to where they came from
      // but the function crashes.

      // multisharpen_box(int imol_map, imol_map_ref, position, box_radius, float b_factor)

      // first create a target molecule for the interactively sharpened map
      int imol_map = scm_to_int(i_scm);
      int imol_map_updating = copy_molecule(imol_map);
      float cl = get_contour_level_absolute(imol_map);
      set_contour_level_absolute(imol_map_updating, cl);

      // now sharpen/blur
      float box_radius = 10.0;
      float b_factor = 12.0; // blur
      graphics_info_t g;
      coot::Cartesian pos = g.RotationCentre();
      clipper::Coord_orth centre(pos.x(), pos.y(), pos.z());

      if (is_valid_map_molecule(imol_map)) {
         const clipper::Xmap<float>    &xmap = g.molecules[imol_map].xmap;
         clipper::Xmap<float> &xmap_updating = g.molecules[imol_map_updating].xmap;

         coot::util::map_fragment_info_t mfi(xmap, centre, box_radius, true);
         mfi.simple_origin_shift(xmap, centre, box_radius);

         // is the map in mfi sensible?
         if (true) {
            int imol_new_map = graphics_info_t::create_molecule();
            std::cout << "INFO:: Creating origin map in molecule number " << imol_new_map << std::endl;
            bool is_em_map_flag = graphics_info_t::molecules[imol_map].is_EM_map();
            is_em_map_flag = false;
            std::string old_name = graphics_info_t::molecules[imol_map].get_name();
            std::string new_name = "Origin Fragment Map from " + old_name;
            g.molecules[imol_new_map].install_new_map(mfi.xmap, new_name, is_em_map_flag);
         }

         // coot::util::sharpen_blur_map(&mfi.xmap, b_factor); // put this inside mfi

         // is the blur map in mfi sensible?
         if (false) {
            int imol_new_map = graphics_info_t::create_molecule();
            std::cout << "INFO:: Creating blurred origin map in molecule number " << imol_new_map << std::endl;
            std::string old_name = graphics_info_t::molecules[imol_map].get_name();
            std::string new_name = "Blur Origin Fragment Map from " + old_name;
            bool is_em_map_flag = false;
            g.molecules[imol_new_map].install_new_map(mfi.xmap, new_name, is_em_map_flag);
         }

         clipper::Xmap_base::Map_reference_coord ix(xmap_updating);

         // now transfer, with shift, mfi.xmap back into xmap_updating
         mfi.unshift(&xmap_updating, centre);

         // recontour molecule imol_map_updating
         g.molecules[imol_map_updating].update_map(true);

      }
      graphics_draw();
   }


   if (false) {

      int imol_1 = read_pdb("good-test-for-out-of-register-errors-em-tutorial-partial.pdb");
      int imol_2 = read_pdb("EMD-3908/fittedModels/PDB/6eoj.ent");

      // imol_1 = read_pdb("6eoj-fragment-RSR-in-coot.pdb");

      if (is_valid_model_molecule(imol_1)) {
         if (is_valid_model_molecule(imol_2)) {
            mmdb::Manager *mol_1 = graphics_info_t::molecules[imol_1].atom_sel.mol;
            mmdb::Manager *mol_2 = graphics_info_t::molecules[imol_2].atom_sel.mol;
            coot::find_out_of_register_errors(mol_1, mol_2);
         }
      }
   }

   if (false) {
      std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
      if (pp.first) {
         int imol = pp.second.first;
	 int imol_map = scm_to_int(i_scm);
	 if (is_valid_map_molecule(imol_map)) {
	    graphics_info_t g;
	    if (g.molecules[imol_map].is_difference_map_p()) {
	       const clipper::Xmap<float> &diff_xmap = g.molecules[imol_map].xmap;
	       mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
	       coot::pepflip_using_difference_map pf(mol, diff_xmap);
	       float n_sigma = 4.0;
	       std::vector<coot::residue_spec_t> flips = pf.get_suggested_flips(n_sigma);
	       for (std::size_t i=0; i<flips.size(); i++) {
		  std::cout << i << " " << flips[i] << std::endl;
	       }
	    }
	 }
      }
   }

   if (false) {
      graphics_info_t g;
      std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
      if (pp.first) {
         int imol = pp.second.first;
         bool show_stub_flag = true;
         GtkWidget *widget = g.wrapped_create_residue_type_chooser_window(show_stub_flag);
         gtk_widget_set_visible(widget, TRUE);
         g.in_mutate_auto_fit_define = 0;
         g.residue_type_chooser_auto_fit_flag = 1;
         g.pick_pending_flag = 0;
      }
   }


   if (false) {
      int imol_model     = scm_to_int(i_scm);
      int imol_with_data = scm_to_int(j_scm);
      if (is_valid_model_molecule(imol_model)) {
         graphics_info_t g;

         std::string cablam_log_file_name = "6kzp-cablam.log";
         atom_selection_container_t asc = g.molecules[imol_model].atom_sel;
         if (asc.read_success) {
            // parse this log file and call the above function for each cablam outlier residue
            std::vector<coot::cablam_markup_t> v =
            coot::make_cablam_markups(asc.mol, cablam_log_file_name);

            std::cout << "Made " << v.size() << " cablam markups " << std::endl;
            std::vector<coot::cablam_markup_t>::const_iterator it;
            int idx_cablam = new_generic_object_number("Cablam");
            set_display_generic_object(idx_cablam, 1);
            for (it=v.begin(); it!=v.end(); it++) {
               const coot::cablam_markup_t &cm(*it);
               to_generic_object_add_point(idx_cablam, "pink", 14, cm.O_prev_pos.x(), cm.O_prev_pos.y(), cm.O_prev_pos.z());
               to_generic_object_add_point(idx_cablam, "pink", 14, cm.O_this_pos.x(), cm.O_this_pos.y(), cm.O_this_pos.z());
               to_generic_object_add_point(idx_cablam, "pink", 14, cm.O_next_pos.x(), cm.O_next_pos.y(), cm.O_next_pos.z());

               std::cout << "line 1: " << cm.O_this_pos.format() << " to " << cm.CA_proj_point_this.format() << std::endl;
               std::cout << "line 2: " << cm.O_prev_pos.format() << " to " << cm.CA_proj_point_prev.format() << std::endl;
               std::cout << "line 3: " << cm.O_next_pos.format() << " to " << cm.CA_proj_point_next.format() << std::endl;

               to_generic_object_add_line(idx_cablam, "pink", 4,
                                          cm.O_this_pos.x(), cm.O_this_pos.y(), cm.O_this_pos.z(),
                                          cm.CA_proj_point_this.x(), cm.CA_proj_point_this.y(), cm.CA_proj_point_this.z());

               to_generic_object_add_line(idx_cablam, "pink", 4,
                                          cm.O_prev_pos.x(), cm.O_prev_pos.y(), cm.O_prev_pos.z(),
                                          cm.CA_proj_point_prev.x(), cm.CA_proj_point_prev.y(), cm.CA_proj_point_prev.z());

               to_generic_object_add_line(idx_cablam, "pink", 4,
                                          cm.O_next_pos.x(), cm.O_next_pos.y(), cm.O_next_pos.z(),
                                          cm.CA_proj_point_next.x(), cm.CA_proj_point_next.y(), cm.CA_proj_point_next.z());

               to_generic_object_add_line(idx_cablam, "pink", 4,
                                          cm.CA_proj_point_this.x(), cm.CA_proj_point_this.y(), cm.CA_proj_point_this.z(),
                                          cm.CA_proj_point_prev.x(), cm.CA_proj_point_prev.y(), cm.CA_proj_point_prev.z());

               to_generic_object_add_line(idx_cablam, "pink", 4,
                                          cm.CA_proj_point_this.x(), cm.CA_proj_point_this.y(), cm.CA_proj_point_this.z(),
                                          cm.CA_proj_point_next.x(), cm.CA_proj_point_next.y(), cm.CA_proj_point_next.z());

            }
         }
      }
   }
   if (false) {
      // this has a bonefide interface now
      int imol_model     = scm_to_int(i_scm);
      int imol_with_data = scm_to_int(j_scm);
      if (is_valid_model_molecule(imol_model)) {
         graphics_info_t g;
         clipper::Xmap<float> *xmap_p = &g.molecules[imol_with_data].xmap;
         try {
            g.molecules[imol_with_data].fill_fobs_sigfobs();
            const clipper::HKL_data<clipper::data32::F_sigF> *fobs_data = g.molecules[imol_with_data].get_original_fobs_sigfobs();
            const clipper::HKL_data<clipper::data32::Flag> *free_flag = g.molecules[imol_with_data] .get_original_rfree_flags();
            if (fobs_data && free_flag) {
               g.molecules[imol_model].sfcalc_genmap(*fobs_data, *free_flag, xmap_p);
               g.molecules[imol_with_data].update_map(true);
               graphics_draw();
            }
         }
         catch (const std::runtime_error &rte) {
            std::cout << rte.what() << std::endl;
         }
      }
   }

   if (false) {
      int imol   = scm_to_int(i_scm);
      int res_no = scm_to_int(j_scm);
      if (is_valid_model_molecule(imol)) {
         std::string chain_id = "A";
         coot::residue_spec_t spec("A", res_no, "");
         mmdb::Residue *r = graphics_info_t::molecules[imol].get_residue(spec);
         if (r) {
            mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
            std::string alt_conf("");
            coot::richardson_rotamer d(r, alt_conf, mol, 0.0, 1);
            coot::rotamer_probability_info_t prob = d.probability_of_this_rotamer();
            std::string rn = residue_name(imol, chain_id, res_no, "");
            std::cout << "INFO:: " << coot::residue_spec_t(r) << " " << rn << " "
                      << prob << " with rotamer name \"" << prob.rotamer_name << "\"" << std::endl;
            // ------------------------
            short int add_extra_PHE_and_TYR_rotamers_flag = 1; // true
            coot::rotamer rotamer(r, alt_conf, add_extra_PHE_and_TYR_rotamers_flag);
            std::vector<std::pair<int,float> > chi_angles = rotamer.get_chi_angles();
            std::cout << "current chi angles ";
            for (unsigned int i=0; i<chi_angles.size(); i++)
               std::cout << " " << chi_angles[i].first << " " << chi_angles[i].second << " ";
            std::cout << " " << std::endl;
            coot::closest_rotamer_info_t closest_rotamer = rotamer.get_closest_rotamer(rn);
            std::cout << " drive to " << closest_rotamer.rotamer_probability_info.rotamer_name << " ";
            for (unsigned int i=0; i<closest_rotamer.residue_chi_angles.size(); i++)
               std::cout << " " << closest_rotamer.residue_chi_angles[i].first << " "
                         << closest_rotamer.residue_chi_angles[i].second << " ";
            std::cout << std::endl;
         }
      }
   }

   if (false) {
      mmdb::Manager *mol = new mmdb::Manager;
      mol->ReadPDBASCII("test.pdb");
      coot::get_c_beta_deviations(mol);
   }

   if (false) {
      int imol_1 = scm_to_int(i_scm); // from
      int imol_2 = scm_to_int(j_scm); // to

      if (is_valid_model_molecule(imol_1)) {
	 if (is_valid_model_molecule(imol_2)) {
	    coot::util::copy_headers(g.molecules[imol_1].atom_sel.mol,
				     g.molecules[imol_2].atom_sel.mol,
				     false); // no crystal info
	    write_pdb_file(imol_2, "copied-here.pdb");
	 }
      }
   }

   if (false) {
      dodec d;
      d.test("dodec.xyz");
   }

   if (false) {
      std::string file_name = scm_to_locale_string(i_scm);
   }

   if (false) {

      std::ifstream f("normal-data-0-1.tab");
      std::vector<double> data;
      if (f) {
	 double v;
	 std::string line;
	 while (std::getline(f, line)) {
	    try {
	       v = coot::util::string_to_float(line);
	       data.push_back(v);
	    }
	    catch (const std::exception &e) {
	       std::cout << e.what() << " failed to read " << line << std::endl;
	    }
	 }
      }
      double mean = 0;
      double variance = 1;
      double d = coot::stats::get_kolmogorov_smirnov_vs_normal(data, mean, variance);
      std::cout << "D: " << d << std::endl;

      // Compare my (Taylor expansion of the error function) with the GSL version.
      if (false) {
	 coot::stats::pnorm pn;
	 for (double v = -5; v<=5; v += 0.2) {
	    double v1 = pn.erf(v);
	    double v2 = gsl_sf_erf(v);
	    std::cout << "v " << std::setw(4) << v << " " << std::setw(10) << v1 << " " << v2 << std::endl;
	    // OK, so they are pretty close
	 }
      }
   }

   if (false) {

      std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
      if (pp.first) {
	 if (is_valid_model_molecule(pp.second.first)) {

	    coot::residue_spec_t residue_spec(pp.second.second);
	    int imol = pp.second.first;
	    g.perpendicular_ligand_view(imol, residue_spec);
	 }
      }
   }


   // ------------------------ spherical density overlap -------------------------
   //
   if (false) {
      int imol = scm_to_int(i_scm); // map molecule
      int imol_map = scm_to_int(j_scm); // map molecule

      if (is_valid_model_molecule(imol)) {
	 if (is_valid_map_molecule(imol_map)) {

	    const clipper::Xmap<float> &m = g.molecules[imol_map].xmap;
	    clipper::Coord_orth c(0,0,0); // (set-rotation-centre -15 -4 21)
	    coot::util::map_fragment_info_t mf(m, c, 50, true);

	    if (mf.xmap.is_null()) {
	       std::cout << "null map fragment xmap " << std::endl;
	    } else {
	       clipper::CCP4MAPfile mapout;
	       mapout.open_write("map-fragment-at-origin.map");
	       mapout.export_xmap(mf.xmap);
	       mapout.close_write();

	       mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
	       coot::util::emma sphd(mol, 5); // 5 is border
	       sphd.overlap_simple(mf.xmap);
	       // sphd.overlap(mf.xmap);
	    }
	 }
      }
   }

   if (0) {

      for (unsigned int io=0; io<20; io++) {
	 std::string name = "Test " + coot::util::int_to_string(io);
	 int n = new_generic_object_number(name.c_str());
	 to_generic_object_add_line(n, "green", 2+io, 1+io, 2, 3, 4, 5, 6);
	 set_display_generic_object(n, 1);
      }

      GtkWidget *w = wrapped_create_generic_objects_dialog();
      gtk_widget_set_visible(w, TRUE);
   }

   if (0) {

      coot::minimol::molecule m;
      m.read_file("oliver-clarke/oliver-clarke-test.pdb");
      m.write_file("mini-mol-out.pdb", 20);
   }

   if (0) {

      for (unsigned int io=0; io<20; io++) {
	 std::string name = "Test " + coot::util::int_to_string(io);
	 int n = new_generic_object_number(name.c_str());
	 to_generic_object_add_line(n, "green", 2+io, 1+io, 2, 3, 4, 5, 6);
	 set_display_generic_object(n, 1);
      }

      GtkWidget *w = wrapped_create_generic_objects_dialog();
      gtk_widget_set_visible(w, TRUE);
   }

   if (0) {
      int imol = 0;
      std::string file_name = "with-mtrix.pdb";
      // file_name = "coot-download/pdb1qex.ent";
      std::vector<clipper::RTop_orth> mv = coot::mtrix_info(file_name);
      for (unsigned int i=0; i<mv.size(); i++) {
	 const clipper::RTop_orth &rt = mv[i];
	 add_strict_ncs_matrix(imol, "A", "A",
			       rt.rot()(0,0), rt.rot()(0,1), rt.rot()(0,2),
			       rt.rot()(1,0), rt.rot()(1,1), rt.rot()(1,2),
			       rt.rot()(2,0), rt.rot()(2,1), rt.rot()(2,2),
			       rt.trn()[0],   rt.trn()[1],   rt.trn()[2]);
      }
   }

   if (0) {

      int imol = scm_to_int(i_scm); // map molecule

      print_residue_distortions(imol, "A", 1, "");

      // now test making a dictionary
      int imol_2 = scm_to_int(j_scm);
      if (is_valid_model_molecule(imol_2)) {
	 mmdb::Residue *residue_2_p = g.molecules[imol_2].get_residue("A", 304, "");
	 if (! residue_2_p) {
	    std::cout << " residue not found " << std::endl;
	 } else {
	    mmdb::Manager *mol = coot::util::create_mmdbmanager_from_residue(residue_2_p);
	    if (mol) {
	       coot::dictionary_residue_restraints_t rest(mol);
	       rest.write_cif("testing.cif");
	    }
	 }
      }
   }

   if (0) {
      std::cout << "size of a molecule " << sizeof(molecule_class_info_t) << std::endl;
   }

#ifdef USE_LIBCURL

   if (0) {
      curl_test_make_a_post();
   }

#endif

   if (0) {
      coot::smcif s;
      s.read_data_sm_cif("hof.fcf");
   }

   if (0) {
      std::cout << "======== n monomers in dictionary: " << g.Geom_p()->size() << std::endl;
      for (unsigned int irest=0; irest<g.Geom_p()->size(); irest++) {
	 std::cout << "   " << irest << "  " << (*g.Geom_p())[irest].first
		   << " " << (*g.Geom_p())[irest].second.residue_info.comp_id << std::endl;
      }
   }

   if (0) {
#if HAVE_GOOCANVAS
      coot::goograph *g = new coot::goograph;
      std::vector<std::pair<double, double> > data;
      data.push_back(std::pair<double, double> ( 104.5,  4));
      data.push_back(std::pair<double, double> ( 104.75, 1));
      data.push_back(std::pair<double, double> ( 105,    2));
      data.push_back(std::pair<double, double> ( 105.25, 1));
      data.push_back(std::pair<double, double> ( 105.5,  0));
      data.push_back(std::pair<double, double> ( 105.75, 0));
      data.push_back(std::pair<double, double> ( 106,    0));
      data.push_back(std::pair<double, double> ( 106.25, 0));
      data.push_back(std::pair<double, double> ( 106.5,  0));
      data.push_back(std::pair<double, double> ( 106.75, 0));
      data.push_back(std::pair<double, double> ( 107,    0));
      data.push_back(std::pair<double, double> ( 107.25, 3));
      data.push_back(std::pair<double, double> ( 107.5,  6));
      data.push_back(std::pair<double, double> ( 107.75, 9));
      int trace = g->trace_new();
      g->set_plot_title("Density Histogram");
      g->set_data(trace, data);
      g->show_dialog();
#endif
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
	    unsigned int n_rounds = 10;
	    mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
	    const clipper::Xmap<float> &xmap = g.molecules[j].xmap;
	    for (unsigned int iround=0; iround<n_rounds; iround++) {
	       std::cout << "round " << iround << std::endl;
	       mmdb::Manager *moving_mol = coot::util::create_mmdbmanager_from_residue_specs(v, mol);

	       std::vector<std::pair<bool, clipper::Coord_orth> > avoid_these_atoms;
	       coot::multi_residue_torsion_fit_map(imol, moving_mol, xmap, avoid_these_atoms, 400, g.Geom_p());
	       atom_selection_container_t moving_atoms_asc = make_asc(moving_mol);
	       std::pair<mmdb::Manager *, int> new_mol =
		  coot::util::create_mmdbmanager_from_mmdbmanager(moving_mol);
	       atom_selection_container_t asc_new = make_asc(new_mol.first);
	       std::string name = "test-" + coot::util::int_to_string(iround);
	       bool shelx_flag = 0;
	       int imol_new = g.create_molecule();
	       g.molecules[imol_new].install_model(imol_new, asc_new, g.Geom_p(), name, 1, shelx_flag);
	       add_linked_residue(imol_new,
				  active_atom.second.second.chain_id.c_str(),
				  active_atom.second.second.res_no,
				  active_atom.second.second.ins_code.c_str(),
				  "NAG", "ASN-NAG", 400);

	       delete moving_mol;
	       graphics_draw();
	    }
	 }
      }
   }

   if (false) {
      int i = scm_to_int(i_scm); // map molecule
      int j = scm_to_int(j_scm);

      if (is_valid_model_molecule(i)) {
         if (is_valid_map_molecule(j)) {
            const clipper::Xmap<float> &xmap = g.molecules[j].xmap;
            g.molecules[i].em_ringer(xmap);
         }
      }
   }

   if (false) {
      int i = scm_to_int(i_scm); // map molecule
      int j = scm_to_int(j_scm);

      if (is_valid_model_molecule(i)) {
         if (is_valid_map_molecule(j)) {
            const clipper::Xmap<float> &xmap = g.molecules[j].xmap;

	    int icount = 0;
	    int imod = 1;
	    mmdb::Model *model_p = graphics_info_t::molecules[i].atom_sel.mol->GetModel(imod);
	    if (model_p) {
	       int n_chains = model_p->GetNumberOfChains();
	       for (int ichain=0; ichain<n_chains; ichain++) {
		  mmdb::Chain *chain_p = model_p->GetChain(ichain);
		  int nres = chain_p->GetNumberOfResidues();
		  for (int ires=0; ires<nres; ires++) {
		     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
		     int n_atoms = residue_p->GetNumberOfAtoms();
		     for (int iat=0; iat<n_atoms; iat++) {
			mmdb::Atom *at = residue_p->GetAtom(iat);
			clipper::Coord_orth co = coot::co(at);
			float dv = coot::util::density_at_point_by_linear_interpolation(xmap, co);
			clipper::Grad_orth<double> grad = coot::util::gradient_at_point(xmap, co);

			std::cout << icount << " " << coot::atom_spec_t(at) << " dv "
				  << dv << " grad " << grad.format() << std::endl;
			icount++;
		     }
		  }
	       }
	    }
	 }
      }
   }


   if (false) {
      int imol_map   = scm_to_int(i_scm);
      float b_factor = scm_to_double(j_scm);
      graphics_info_t g;
      int imol_new = graphics_info_t::create_molecule();
      const clipper::Xmap<float> &xmap_orig = g.molecules[imol_map].xmap;
      clipper::Xmap<float> xmap_new = coot::util::sharpen_blur_map(xmap_orig, b_factor);
      bool is_em_flag = false;
      g.molecules[imol_new].install_new_map(xmap_new, "Blur map map", is_em_flag);
      float contour_level = graphics_info_t::molecules[imol_map].get_contour_level();
      graphics_info_t::molecules[imol_new].set_contour_level(contour_level);
      graphics_draw();
   }


   if (false) {

      try {
	 int imol_map   = scm_to_int(i_scm);
	 graphics_info_t g;
	 const clipper::Xmap<float> &xmap_orig = g.molecules[imol_map].xmap;
	 std::vector<float> b_factors = {-100, -50, 5, 50, 100, 160};
	 std::vector<clipper::Xmap<float> > xmaps(b_factors.size());
	 coot::util::multi_sharpen_blur_map(xmap_orig, b_factors, &xmaps);
	 float contour_level = g.molecules[imol_map].get_contour_level();
	 for (std::size_t i=0; i<b_factors.size(); i++) {
	    const clipper::Xmap<float> &xmap_new = xmaps[i];
	    float b_factor = b_factors[i];
	    int imol_new = graphics_info_t::create_molecule();
	    std::string map_name = "Map";
	    if (b_factor < 0)
	       map_name += " Sharpen ";
	    else
	       map_name += " Blur ";
	    map_name += coot::util::float_to_string(b_factor);
	    bool is_em_flag = false;
	    g.molecules[imol_new].install_new_map(xmap_new, map_name, is_em_flag);
	    graphics_info_t::molecules[imol_new].set_contour_level(contour_level*exp(-0.01*b_factor));
	 }
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "ERROR:: " << rte.what() << std::endl;
      }
   }

#if 0 // what is this doing here? Added to the wrong file?
   if (true) {
      graphics_info_t::draw_the_other_things = true;
      gtk_widget_queue_draw(glarea);
   }
#endif

   return r;
}

#endif


#include "coot-utils/ortep.hh" // ortep 20230108-PE // 20241202-PE
#include "utils/dodec.hh"
#include "widget-from-builder.hh"

#include "density-contour/gaussian-surface.hh"

#include "pli/sdf-interface-for-export.hh"
#include "get-monomer.hh"

#include <boost/python.hpp> // because we try to call python functions from c++ using an rdkit molecule

void pandda();

#ifdef USE_PYTHON
static PyObject *rdkit_chem_qed_func        = NULL;
static PyObject *rdkit_chem_properties_func = NULL;
PyObject *test_function_py(PyObject *i_py, PyObject *j_py) {

   std::cout << "-------------------------- test_function_py() " << std::endl;

   if (true) {
      pandda();
   }

   if (false)   {
      int imol = 0;
      graphics_info_t::molecules[imol].debug_ghosts();
   }

   if (false) {
      std::string d = coot::prefix_dir();
      std::cout << "--------- prefix_dir " << d << std::endl;
   }


#if 0
   if (true) {
      auto qed_stuff = [] (RDKit::RWMol &rdkm) {
         PyObject *pName = PyUnicode_FromString("rdkit.Chem"); // 20230513-PE merge: was just "silicos_it" - hmm.
         // Load the module object
         PyObject *pModule = PyImport_Import(pName);
         if (pModule == NULL) {
            std::cout << "OK... pModule was null " << pModule << std::endl;
         } else {
            std::cout << "OK... pModule was not null " << pModule << std::endl;
            PyObject *pName_qed = PyUnicode_FromString("rdkit.Chem.QED");
            // pDict is a borrowed reference
            PyObject *pModule_qed = PyImport_Import(pName_qed);
            PyObject *pDict_qed = PyModule_GetDict(pModule_qed);
            if (! PyDict_Check(pDict_qed)) {
               std::cout << "pDict is not a dict"<< std::endl;
            } else {
               std::cout << "OK... found the QED dictionary" << std::endl;

               PyObject *pFunc = PyDict_GetItemString(pDict_qed, "default");
               if (! pFunc) {
                  std::cout << "BAD:: pFunc default is NULL" << std::endl;
               } else {
                  if (PyCallable_Check(pFunc)) {
                     std::cout << "OK... Yeay - storing rdkit chem qed default_func" << std::endl;
                     rdkit_chem_qed_func = pFunc;
                  } else {
                     std::cout << "BAD:: default() function is not callable"  << std::endl;
                  }
               }

               pFunc = PyDict_GetItemString(pDict_qed, "properties");
               if (! pFunc) {
                  std::cout << "BAD:: pFunc properties is NULL" << std::endl;
               } else {
                  if (PyCallable_Check(pFunc)) {
                     std::cout << "OK... Yeay - storing rdkit chem qed properties_func" << std::endl;
                     rdkit_chem_properties_func = pFunc;
                  } else {
                     std::cout << "BAD:: properties() function is not callable"  << std::endl;
                  }
               }

               // pads2 has been renamed to adsParameters
               PyObject *rdkit_chem_qed_pads = PyDict_GetItemString(pDict_qed, "adsParameters");
            }
         }
      };

      graphics_info_t g;
      std::string ligand_res_name("801");
      int imol = get_monomer(ligand_res_name);
      std::cout << "imol: " << imol << std::endl;
      std::pair<bool, coot::dictionary_residue_restraints_t> p = 
         g.Geom_p()->get_monomer_restraints_at_least_minimal(ligand_res_name, imol);
      if (! p.first) {
         std::cout << "Failed to get dictionary: " << ligand_res_name << std::endl;
      } else {
         mmdb::Residue  *res_ref = g.molecules[imol].get_residue("A", 1, "");
         mmdb::Manager *mol_for_res_ref = g.molecules[imol].atom_sel.mol;
         if (res_ref) {
	    try {

	       // can throw a std::runtime_error
	       RDKit::RWMol rdkm = coot::rdkit_mol(res_ref, imol, *g.Geom_p());

	       // assign atom names
	       if (int(rdkm.getNumAtoms()) < res_ref->GetNumberOfAtoms()) {
		  std::cout << "WARNING:: failure to construct rdkit molecule " << rdkm.getNumAtoms() << " vs " << res_ref->GetNumberOfAtoms()
	                    << std::endl;
	       } else {
		  mmdb::PPAtom residue_atoms = 0;
		  int n_residue_atoms;
		  res_ref->GetAtomTable(residue_atoms, n_residue_atoms);

		  // polar Hs only, that is - need new function here.
		  // (can throw a std::exception)
		  coot::undelocalise(&rdkm);
		  coot::assign_formal_charges(&rdkm);
		  // coot::remove_non_polar_Hs(&rdkm);

		  // we need to sanitizeMol() after remove_non_polar_Hs, and have
		  // a kekulized representation.
		  // we need to sanitize to get ring info,
		  // then we need to kekulize because we are making a 2D chemical diagram
		  //
		  // failed_op sets set by sanitizeMol() - perhaps we shouldn't ignore it.
		  //
		  unsigned int failed_op_1 = 0;
		  unsigned int failed_op_2 = 0;
		  RDKit::MolOps::sanitizeMol(rdkm, failed_op_1, RDKit::MolOps::SANITIZE_ALL);
		  RDKit::MolOps::sanitizeMol(rdkm, failed_op_2, RDKit::MolOps::SANITIZE_KEKULIZE);

		  std::cout << "DEBUG:: sanitizeMol() returned with failed_op: "
			    << failed_op_1 << " " << failed_op_2
			    << " (note 'no-failure' is value 0)." << std::endl;

		  try {
		     RDKit::RingInfo *ri = rdkm.getRingInfo();
		     unsigned int n_rings = ri->numRings();
                     std::cout << "----------- OK good so far, n_rings " << n_rings << std::endl;
                     qed_stuff(rdkm);

                     if (rdkit_chem_qed_func) {
                        // yay.
                        double r = -1.0;
                        PyObject *arg_list = PyTuple_New(1);
                        RDKit::ROMol *mol_copy_p = new RDKit::ROMol(rdkm);
                        boost::shared_ptr<RDKit::ROMol> xx(mol_copy_p);
                        boost::python::object obj(xx);
                        PyObject *rdkit_mol_py = obj.ptr();
                        PyTuple_SetItem(arg_list, 0, rdkit_mol_py);
                        PyObject *result_default = PyEval_CallObject(rdkit_chem_qed_func, arg_list);
                        if (result_default) {
                           if (PyFloat_Check(result_default)) {
                              r = PyFloat_AsDouble(result_default);
                              std::cout << " ...................... r " << r << std::endl;
                           } else {
                              std::cout << "BAD:: result default was not a float " << std::endl;
                           }
                        } else {
                           std::cout << "BAD:: null result_default " << std::endl;
                        }
                     }
                  }
		  catch (const std::runtime_error &rte) {
		     std::vector<std::vector<int> > ring_info;
		     RDKit::MolOps::findSSSR(rdkm, ring_info);
		  }
               }
            }
            catch (const std::runtime_error &rte) {
               std::cout << "ERROR:: (runtime error) in fle_view_with_rdkit(): "
                         << rte.what() << std::endl;
            }
            catch (const std::exception &e) {
               std::cout << "ERROR (exception) in fle_view_with_rdkit(): " << e.what() << std::endl;
            }
         }
      }
   }
#endif

#if 0
   // There are lots of items in the list, but they don't seem to be toplevels - hmm.

   GList* top_level_list = gtk_window_list_toplevels();
   while (top_level_list) {
      GdkToplevel* toplevel = GDK_TOPLEVEL(top_level_list->data);
      if (GDK_IS_TOPLEVEL(toplevel)) {
         GdkToplevelState state = gdk_toplevel_get_state(toplevel);
         std::cout << "Debug:: toplevel state " << " " << state << std::endl;
      } else {
         std::cout << "that was not a toplevel " << toplevel << std::endl;
      }
      top_level_list = top_level_list->next;
   };
#endif


   graphics_info_t g;
   PyObject *r = Py_False;

   auto vnc_vertex_to_generic_vertex = [] (const coot::api::vnc_vertex &v) {
      return s_generic_vertex(v.pos, v.normal, v.color);
   };

   auto vnc_vertex_vector_to_generic_vertex_vector = [vnc_vertex_to_generic_vertex] (const std::vector<coot::api::vnc_vertex> &vv) {
      std::vector<s_generic_vertex> vo(vv.size());
      for (unsigned int i=0; i<vv.size(); i++)
         vo[i] = vnc_vertex_to_generic_vertex(vv[i]);
      return vo;
   };


   if (false) {
      int i = PyLong_AsLong(i_py);
      int j = PyLong_AsLong(j_py);

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

      int imol = i;
      coot::residue_spec_t res_spec("A", 1, "");
      mmdb::Residue *residue_p = g.get_residue(imol, res_spec);
      std::vector<coot::simple_mesh_t> meshes = chemical_features::generate_meshes(imol, residue_p, *g.Geom_p());

      std::cout << "generate_meshes() made " << meshes.size() << " meshes" << std::endl;

      for (unsigned int ii=0; ii<meshes.size(); ii++) {
         const auto &mesh = meshes[ii];
         std::cout << "mesh " << mesh.name << " " << mesh.vandt() << std::endl;
         g.attach_buffers();

         std::string object_name("Chemical Feature for #");
         object_name += std::to_string(imol);
         object_name += " ";
         object_name += mesh.name;
         int obj_mesh = new_generic_object_number(object_name);
         meshed_generic_display_object &obj = g.generic_display_objects[obj_mesh];
         obj.mesh.name = object_name;
         obj.mesh.set_draw_mesh_state(true);
         std::vector<s_generic_vertex> converted_vertices = vnc_vertex_vector_to_generic_vertex_vector(mesh.vertices);
         obj.mesh.import(converted_vertices, mesh.triangles);
         obj.mesh.set_material_specularity(0.7, 128);
         obj.mesh.setup_buffers();
      }
      graphics_draw();
#else
      std::cout << "No MAKE_ENHANCED_LIGAND_TOOLS (sad-face)" << std::endl;
#endif // MAKE_ENHANCED_LIGAND_TOOLS
   }

   if (false) {
      int i = PyLong_AsLong(i_py);
      int j = PyLong_AsLong(j_py);

      auto make_a_surface = [] (int imol, mmdb::Manager *mol,
                                unsigned int i_ch, const std::string &chain_id,
                                const std::vector<coot::ghost_molecule_display_t> &gi,
                                bool colour_by_ncs_ghost,
                                const std::map<std::string, int> &chain_id_map) {

         graphics_info_t g;
         coot::colour_holder ch(0.8, 0.4, 0.4);

         if (colour_by_ncs_ghost) {
            for (const auto &ghost : gi) {
               if (ghost.chain_id == chain_id) {
                  const std::string &target_chain_id = ghost.target_chain_id;
                  std::map<std::string, int>::const_iterator it = chain_id_map.find(target_chain_id);
                  if (it != chain_id_map.end()) {
                     int i_ch_for_target = it->second;
                     ch.rotate_by(0.22 * i_ch_for_target);
                  }
               }
            }
         } else {
            ch.rotate_by(0.22 * i_ch);
         }
         glm::vec4 col = colour_holder_to_glm(ch);

         coot::gaussian_surface_t gauss_surf(mol, chain_id);
         coot::simple_mesh_t smesh = gauss_surf.get_surface();
         std::vector<s_generic_vertex> vertices(smesh.vertices.size());
         for (unsigned int i = 0; i < smesh.vertices.size(); i++) {
            vertices[i] = s_generic_vertex(smesh.vertices[i].pos,
                                          smesh.vertices[i].normal,
                                          smesh.vertices[i].color);
            vertices[i].color = col;
            // std::cout << i << " " << glm::to_string(vertices[i].pos) << "\n";
         }

         g.attach_buffers();

         std::string object_name("Gaussian Surface #");
         object_name += std::to_string(imol);
         object_name += " Chain ";
         object_name += chain_id;
         int obj_mesh = new_generic_object_number(object_name);
         meshed_generic_display_object &obj = g.generic_display_objects[obj_mesh];
         obj.mesh.name = object_name;
         obj.mesh.set_draw_mesh_state(true);
         obj.mesh.import(vertices, smesh.triangles);
         obj.mesh.set_material_specularity(0.7, 128);
         obj.mesh.setup_buffers();
         g.graphics_draw();
      };


      if (is_valid_model_molecule(i)) {

         std::vector<drawn_ghost_molecule_display_t> gi = g.molecules[i].NCS_ghosts();

         mmdb::Manager *mol = g.molecules[i].atom_sel.mol;
         std::vector<std::string> chain_ids = g.molecules[i].get_chain_ids();
         std::map<std::string, int> chain_id_map;
         for (unsigned int i_ch=0; i_ch<chain_ids.size(); i_ch++) {
            auto chain_id = chain_ids[i_ch];
            chain_id_map[chain_id] = i_ch;
         }

         if (false) {
            std::cout << ":::::::::::::::::::::::: ghosts " << gi.size() << std::endl;
            for (unsigned int ighost=0; ighost<gi.size(); ighost++) {
               const auto &ghost = gi[ighost];
               std::cout << "   " << ighost << " " << ghost.name << " " << ghost.chain_id << " " << ghost.target_chain_id
                         << "\n" << ghost.rtop.format() << std::endl;
            }
         }

         for (unsigned int i_ch=0; i_ch<chain_ids.size(); i_ch++) {
            auto chain_id = chain_ids[i_ch];
            // gi is used if colour_by_ncs_ghost is true
            bool colour_by_ncs_ghost = false;
            // this doesn't compile
            // make_a_surface(i, mol, i_ch, chain_id, gi, colour_by_ncs_ghost, chain_id_map);
         }

      }
   }

   if (false) {
     int i = PyLong_AsLong(i_py);
     int j = PyLong_AsLong(j_py);

     GtkWidget *gl_area = graphics_info_t::glareas[0];

     GtkAllocation allocation;
     gtk_widget_get_allocation(graphics_info_t::glareas[0], &allocation);
     int w = allocation.width;
     int h = allocation.height;
     int w_new = w + 8;
     if (i < 0) w_new = w - 8;

     GtkWidget *main_window = graphics_info_t::get_main_window();

     gtk_window_set_resizable(GTK_WINDOW(main_window), TRUE);
     gtk_widget_set_size_request(GTK_WIDGET(gl_area), w_new, h);
     gtk_window_set_resizable(GTK_WINDOW(main_window), FALSE);
   }

   if (false) {

      g.attach_buffers();

      // unsigned int num_subdivisions = 3; // tessellate_sphere_sans_octant is built based on 3
      ortep_t o = tessellate_sphere_sans_octant();
      glm::vec4 col(0.5, 0.5, 0.98, 1.0);

      // ------------------------------------ solid -----------------------------------------

      {
         std::string object_name("Ortep Testing Solid");
         int obj_index = new_generic_object_number(object_name);
         meshed_generic_display_object &obj = g.generic_display_objects[obj_index];

         std::vector<s_generic_vertex> vertices(o.vertices.size());
         for (unsigned int i=0; i<vertices.size(); i++) {
            auto &v = vertices[i];
            v.pos    = o.vertices[i];
            v.normal = o.normals[i];
            v.color  = col;
         }
         obj.mesh.import(vertices, o.triangles);
         obj.mesh.set_draw_mesh_state(true);
         obj.mesh.setup_buffers();
      }

      // ------------------------------------ lines -----------------------------------------

      {
         g.attach_buffers();
         std::string object_name = std::string("Ortep Testing Lines");
         int obj_index = new_generic_object_number(object_name);
         meshed_generic_display_object &obj = g.generic_display_objects[obj_index];

         std::vector<s_generic_vertex> vertices = std::vector<s_generic_vertex>(o.vertices_for_lines.size());
         glm::vec4 black(0.1, 0.1, 0.1, 1.0);
         for (unsigned int i=0; i<vertices.size(); i++) {
            auto &v  = vertices[i];
            v.pos    = o.vertices_for_lines[i];
            v.normal = o.vertices_for_lines[i]; // not used
            v.color  = black;
         }
         std::vector<unsigned int> lvi(o.indices_for_lines.size() * 2);
         for (unsigned int i=0; i<o.indices_for_lines.size(); i++) {
            lvi[2 * i    ] = o.indices_for_lines[i].first;
            lvi[2 * i + 1] = o.indices_for_lines[i].second;
         }
         obj.wireframe_mode = true;
         obj.mesh.set_draw_mesh_state(true);
         obj.mesh.import_lines(vertices, lvi);
         obj.mesh.setup_buffers();
      }

      // ----------------- now draw it --------------------------

      graphics_draw();
   }

   if (false) {
#ifdef USE_MOLECULES_TO_TRIANGLES
      Mesh mesh("test-molecule-to-triangle");
      // mesh.load_from_glTF("blue-eyelashes-1.glb");
      // mesh.load_from_glTF("augmented-box.glb");
      // mesh.load_from_glTF("Triangle.gltf");
      //nnmesh.draw_this_mesh = true;
      // int imol = 0;
      // meshed_generic_display_object mgdo(mesh);
      // g.generic_display_objects.push_back(mgdo);
      // int obj_index = new_generic_object_number_for_molecule(object_name, imol);
      // set_display_generic_object_simple(obj_index, 1);

      std::string file_name = PyBytes_AS_STRING(PyUnicode_AsUTF8String(i_py));

      gtk_gl_area_attach_buffers(GTK_GL_AREA(graphics_info_t::glareas[0]));

      if (true) {
         std::string object_name("Imported Object");
         int obj_mesh = new_generic_object_number(object_name);
         meshed_generic_display_object &obj = g.generic_display_objects[obj_mesh];
         obj.mesh.name = object_name;
         Material material;
         material.shininess = 30.0;
         material.specular_strength = 0.3;

         // obj.mesh.load_from_glTF("augmented-box.glb");
         bool status = obj.mesh.load_from_glTF(file_name);

         std::cout << "debug:: obj mesh load_from_glTF returned status " << status << std::endl;
// not always available
#ifdef USE_MOLECULES_TO_TRIANGLES
         if (status) {

            // obj.mesh.setup(material); // calls setup_buffers() (again)
            // set_display_generic_object(obj_mesh, 1);

            // can I get a ribbon diagram out?
            int imol = read_pdb("A-chain.pdb");
            std::string ass = "//A";
            std::string cs = "colorRampChainsScheme";
            std::string style = "Ribbon";
            int secondary_structure_usage_flag = CALC_SECONDARY_STRUCTURE;
            graphics_info_t::molecules[imol].add_molecular_representation(ass, cs, style, secondary_structure_usage_flag);
            std::cout << "debug:: meshes size " << graphics_info_t::molecules[imol].meshes.size() << std::endl;
            auto &mesh = graphics_info_t::molecules[imol].meshes[0];
            bool is_binary_format = true; // no need for this here. binary type is tested internally now.
            mesh.export_to_glTF("exported-to-gltf.glb", is_binary_format);
         }
#endif // USE_MOLECULES_TO_TRIANGLES
      }

      if (false) {
         TextureMesh tm;
         g.texture_meshes.push_back(tm);
         g.texture_meshes.back().load_from_glTF(file_name);
      }

      GtkWidget *overlay = widget_from_builder("main_window_graphics_overlay");
      GtkWidget *label = gtk_label_new("Ramachandran Restraints Enabled ");
      gtk_widget_set_visible(label, TRUE);
      gtk_overlay_add_overlay(GTK_OVERLAY(overlay), label);

      // 20220807-PE not declared in this scope - hmm
      // gtk_overlay_set_overlay_pass_through(GTK_OVERLAY(overlay), label, TRUE);

      // top right
      gtk_widget_set_halign(label, GTK_ALIGN_END);
      gtk_widget_set_valign(label, GTK_ALIGN_START);
#endif
   }


   if (false) {
      GtkWidget *glarea = g.glareas[0];
      std::cout << "DEBUG:: sending " << glarea << " size request 900 900 " << std::endl;
      gtk_widget_set_size_request(glarea, 900, 900);
      g.graphics_draw();
   }

   if (false) {
      int imol = 0;
      if (is_valid_model_molecule(imol)) {
            molecule_class_info_t &m = g.molecules[imol];
            gl_rama_plot_t rama;
            gl_rama_plot_t::draw_mode_t draw_mode = gl_rama_plot_t::draw_mode_t::CHECK_IF_PICKED;
            rama.setup_from(imol, m.atom_sel.mol, "//", draw_mode);
      }
   }

   if (false) {
      // g.show_test_buttons(); // it got renamed show_accept_reject_hud_buttons()
      g.show_accept_reject_hud_buttons();
   }

   if (false) {
      int imol = g.create_molecule();

      std::string mtz_file_name = "coot-download/r1ucssf.mtz";
      mtz_file_name = "coot-download/r2xirsf.mtz";
      std::string i_col = "I.I_sigI.I";
      std::string sigi_col = "I.I_sigI.sigI";
      int status = g.molecules[imol].make_patterson_using_intensities(mtz_file_name,
								      i_col, sigi_col,
								      g.map_sampling_rate);
   }

   if (false) {
      dodec d;
      d.test("dodec.xyz");
      d.face(0);
   }

   if (false) {
     int i = PyLong_AsLong(i_py); // map molecule
     int j = PyLong_AsLong(j_py);

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
         unsigned int n_rounds = 10;
         mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
         const clipper::Xmap<float> &xmap = g.molecules[j].xmap;
         for (unsigned int iround=0; iround<n_rounds; iround++) {
	       std::cout << "round " << iround << std::endl;
	       mmdb::Manager *moving_mol = coot::util::create_mmdbmanager_from_residue_specs(v, mol);

	       std::vector<std::pair<bool, clipper::Coord_orth> > avoid_these_atoms;
	       coot::multi_residue_torsion_fit_map(imol, moving_mol, xmap, avoid_these_atoms, 400, g.Geom_p());
	       atom_selection_container_t moving_atoms_asc = make_asc(moving_mol);
	       std::pair<mmdb::Manager *, int> new_mol =
             coot::util::create_mmdbmanager_from_mmdbmanager(moving_mol);
	       atom_selection_container_t asc_new = make_asc(new_mol.first);
	       std::string name = "test-" + coot::util::int_to_string(iround);
	       bool shelx_flag = 0;
	       int imol_new = g.create_molecule();
	       g.molecules[imol_new].install_model(imol_new, asc_new, g.Geom_p(), name, 1, shelx_flag);
	       add_linked_residue(imol_new,
				  active_atom.second.second.chain_id.c_str(),
				  active_atom.second.second.res_no,
				  active_atom.second.second.ins_code.c_str(),
				  "NAG", "ASN-NAG", 400);

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


/* glyco tools test  */
void glyco_tree_test() {

   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();

   if (pp.first) {
      int imol = pp.second.first;
      graphics_info_t g;
      coot::residue_spec_t rs(coot::atom_spec_t(pp.second.second));
      mmdb::Residue *residue_p = g.molecules[imol].get_residue(rs);
      mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;

      std::vector<std::string> types_with_no_dictionary =
	 g.molecules[imol].no_dictionary_for_residue_type_as_yet(*g.Geom_p());
      std::cout << "glyco-test found " << types_with_no_dictionary.size()
		<< " types with no dictionary" << std::endl;
      for (unsigned int i=0; i<types_with_no_dictionary.size(); i++) {
	 std::cout << "trying to dynamic add: " << types_with_no_dictionary[i] << std::endl;
	 g.Geom_p()->try_dynamic_add(types_with_no_dictionary[i], 41);
      }

      coot::glyco_tree_t(residue_p, mol, g.Geom_p());
   }

}
