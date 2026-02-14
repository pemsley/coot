/* src/molecule-class-info.cc
 *
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2007 by Paul Emsley
 * Copyright 2007, 2008, 2009 by The University of Oxford
 * Copyright 2013, 2014, 2015, 2016 by Medical Research Council
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

#include "molecule-class-info.h"
#include "stereo-eye.hh"
#include <filesystem>
#ifdef USE_PYTHON
#include <Python.h> // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#ifndef EMSCRIPTEN
#include <epoxy/gl.h>
#endif

#include "compat/coot-sysdep.h"

#include <stdlib.h>

#ifndef _MSC_VER
#include <unistd.h>
#else
#include <windows.h>
#include <direct.h>
#endif

#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>

// For stat, mkdir:
#include <sys/types.h>
#include <sys/stat.h>

#define _USE_MATH_DEFINES
#include <cmath>
const double pi = M_PI;

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/rotate_vector.hpp>
#include <glm/gtx/string_cast.hpp>  // to_string()

#include <string.h> // strcmp

#include <mmdb2/mmdb_manager.h>
#include <mmdb2/mmdb_tables.h>

#include <clipper/ccp4/ccp4_mtz_io.h>
#include <clipper/ccp4/ccp4_map_io.h>
#include <clipper/core/xmap.h>
#include <clipper/cns/cns_map_io.h>
#include <clipper/core/hkl_compute.h>
#include <clipper/core/map_utils.h> // Map_stats
#include <clipper/core/resol_basisfn.h>
#include <clipper/core/resol_targetfn.h>
#include <clipper/mmdb/clipper_mmdb.h>
#include <clipper/clipper-phs.h>
#include <clipper/contrib/sfcalc_obs.h>
#include <clipper/contrib/sfscale.h>
#include <clipper/contrib/sfweight.h>

#include "coords/mmdb-extras.hh"
#include "coords/mmdb.hh"
#include "coords/mmdb-crystal.hh"
#include "gtk-manual.hh"

#include "coords/Bond_lines.hh"

#include "coot-utils/gl-matrix.h"
#include "graphics-info.h"

#include "coords/Bond_lines_ext.hh"

#include "coot-utils/coot-coord-utils.hh"
#include "utils/coot-utils.hh"

// #include <GL/glut.h> // needed (only?) for wirecube

#ifndef CLIPPER_MAP_INTERP
#include "clipper/core/map_interp.h"
#endif

#include "ligand/ligand.hh"
#include "ligand/residue_by_phi_psi.hh"
#include "mini-mol/mini-mol-utils.hh"

// #include "cylinder-with-rotation-translation.hh" // for bonds

// for debugging
#include "c-interface.h"

#include "molecular-mesh-generator.hh"
#include "make-a-dodec.hh"

#include "api/coot-molecule.hh" // the integration with api begins...

#include "widget-from-builder.hh"

#include "utils/logging.hh"
extern logging logger;


glm::vec3
cartesian_to_glm(const coot::Cartesian &c) {
   return glm::vec3(c.x(), c.y(), c.z());
}

void
molecule_class_info_t::setup_internal() { // init

   atom_sel.atom_selection = NULL;
   atom_sel.n_selected_atoms = 0;
   atom_sel.mol = NULL;

   // while zero maps, don't need to intialise the arrays (xmap_is_filled)
   is_patterson = 0;
   // draw_vectors = NULL;
   // diff_map_draw_vectors = NULL;

   //
   xskel_is_filled = 0; // not filled.
   skeleton_treenodemap_is_filled = 0;

   this_molecule_has_crystallographic_symmetry = false;
   draw_hydrogens_flag = 1;
   bond_width = 3.0;
   // By default atoms have the same radius as bonds. For ball and stick
   // we want the radius to be bigger. So allow the user to control that.
   atom_radius_scale_factor = 1.0; // used in making balls for atoms

   show_atoms_as_aniso_flag = false;
   show_aniso_atoms_as_ortep_flag = false;
   show_aniso_atoms_as_empty_flag = false;

   ghost_bond_width = 2.0;

   // initial bonds type (checked and reset in handle_read_draw_molecule)
   bonds_box_type = coot::UNSET_TYPE;
   bonds_rotate_colour_map_flag = 0;

   model_representation_mode = Mesh::representation_mode_t::BALL_AND_STICK;
   save_time_string = "";

   pickable_atom_selection = 1;

   is_intermediate_atoms_molecule = false;

   // refmac stuff
   //
   refmac_count = 0;
   have_sensible_refmac_params = 0;  // initially no refmac params.
   have_refmac_phase_params    = 0;  // initially no refmac params.

   // history stuff
   //
   history_index = 0;
   max_history_index = 0;

   have_unsaved_changes_flag = 0; // no unsaved changes initially.
   show_unit_cell_flag = 0;
   fc_skeleton_draw_on = 0;
   greer_skeleton_draw_on = 0;

   // Don't draw it initially.
   draw_it = 0;
   draw_model_molecule_as_lines = false;
   draw_it_for_map = 0;
   draw_it_for_map_standard_lines = 1;
   draw_it_for_extra_restraints = true;
   extra_restraints_representation_for_bonds_go_to_CA = false;

   // backup on by default, turned off for dummy atoms (baton building)
   backup_this_molecule = 1;

   // Map stuff
   map_max_ = 100.0;
   map_min_ = -100.0;
   sharpen_b_factor_ = 0.0;
   sharpen_b_factor_kurtosis_optimised_ = -999999.0;
   pending_contour_level_change_count = 0;
   data_resolution_ = -1; // unset

   use_vertex_gradients_for_map_normals_flag = false;

   // fourier (for phase recombination (potentially) in refmac:
   fourier_weight_label = ""; // unset initially.

   // HL coeff and phi (for phase recombination (potentially) in refmac:
   // should be enough to unset the first one for testing for HL
   refmac_phi_col = ""; // unset initially.
   refmac_hla_col = ""; // unset initially.

   //
   colour_skeleton_by_random = 0;

   // original Fs saved? (could be from map)
   original_fphis_filled = 0;
   original_fobs_sigfobs_filled = false;
   original_fobs_sigfobs_fill_tried_and_failed = false;
   original_r_free_flags_p = 0; // deleted on close_yourself()

   //  bond width (now changeable).
   bond_width = 3.0;
   display_stick_mode_atoms_flag = true;

   // bespoke colouring
   use_bespoke_grey_colour_for_carbon_atoms = false;
   bespoke_carbon_atoms_colour = coot::colour_t(0.4, 0.4, 0.4);

   //
   rotate_colour_map_for_difference_map = 240.0; // degrees

   // save index
   coot_save_index = 0; // how is this related to the backup history_index?

   other_molecule_backup_index = -1; // unset

   // ligand flipping. save the index
   ligand_flip_number = 0; // or 1 2 3 for the 4 different
   // orientations round the eigen vectors.
   // (For ligand molecules).  In future,
   // this may need to be keyed on the
   // resno and chain_id for multiple
   // ligands in a molecule.

   show_symmetry = 1; // individually on by default.

   // Draw NCS ghosts?
   show_ghosts_flag = 0;
   is_dynamically_transformed_map_flag = 0;
   ncs_ghosts_have_rtops_flag = 0;

   // contour by sigma
   contour_by_sigma_flag = 1;
   contour_sigma_step = 0.1;

   //
   // cootsurface = NULL; // no surface initial, updated by make_surface() removed for now
   theSurface = 0;
   transparent_molecular_surface_flag = 0;

   m_VertexArrayID_for_map = VAO_NOT_SET;
   //
   theMapContours.first = 0;
   theMapContours.second = 0;
   is_em_map_cached_flag = -1; // unset
   n_vertices_for_map_VertexArray = 0;
   n_vertices_for_model_VertexArray = 0;
   n_indices_for_triangles = 0;
   n_indices_for_lines = 0;
   // Assigning to GLuint. Hmm
   m_VertexArrayID_for_map  = -1;
   m_VertexArrayID_for_map_cap  = -1;
   m_VertexBufferID = -1;
   m_IndexBuffer_for_map_lines_ID  = -1;
   m_IndexBuffer_for_map_triangles_ID = -1;
   // m_ColourBufferID = -1; // 20220211-PE now a map is a Mesh.
   m_VertexArray_for_model_ID = -1;
   m_VertexBuffer_for_model_ID = -1;
   m_IndexBuffer_for_model_ID = -1;
   n_indices_for_model_triangles = 0;
   n_vertices_for_map_cap = 0;
   shader_shininess = 6.0;
   shader_specular_strength = 0.5;

   map_contours_outdated = false;
   map_mesh_first_time = true;
   model_mesh_first_time = true;

   material_for_maps.do_specularity = false;
   material_for_maps.specular_strength = 0.5; // non-shiny maps by default.

   material_for_models.do_specularity = true;
   material_for_models.specular_strength = 1.0;

   map_as_mesh.set_name("empty map molecule mesh");
   model_molecule_meshes.set_name("empty model molecule mesh");

   // molecule_as_mesh_atoms_1   = Mesh("molecule_as_mesh_atoms_1");
   // molecule_as_mesh_atoms_2   = Mesh("molecule_as_mesh_atoms_2");
   // molecule_as_mesh_bonds_c00 = Mesh("molecule_as_mesh_bonds_c00");
   // molecule_as_mesh_bonds_c01 = Mesh("molecule_as_mesh_bonds_c01");
   // molecule_as_mesh_bonds_c10 = Mesh("molecule_as_mesh_bonds_c10");
   // molecule_as_mesh_bonds_c11 = Mesh("molecule_as_mesh_bonds_c11");

   // draw vectors
   draw_vector_sets.reserve(120); // more than enough
   draw_vector_sets.resize(120);
   draw_diff_map_vector_sets.resize(120);

   write_model_vertices_and_triangles_to_file_mode = false;

   // don't show strict ncs unless it's turned on.
   show_strict_ncs_flag = 1;

   // shelx stuff
   is_from_shelx_ins_flag = 0;

   // symmetry, 20060202 valgrind found that these was not set.
   symmetry_rotate_colour_map_flag = 0;
   symmetry_as_calphas = 0;
   symmetry_whole_chain_flag = 0;
   symmetry_colour_by_symop_flag = 0;

   // dots colour can now be set from outside.
   //
   dots_colour = coot::colour_t(0.3,0.4,0.5);
   dots_colour_set = false; // so use atom-element colouring
   dots_colour_set = true; // don't use atom-element colouring

   // solid surface density representation
   //
   // draw_it_for_solid_density_surface = 0;
   density_surface_opacity = 1.0;

   // other map
   colour_map_using_other_map_flag = false;
   other_map_for_colouring_p = 0;

   // animated ligand interaction representation
   draw_animated_ligand_interactions_flag = 0;

   // single model view
   single_model_view_current_model_number = 0; // all models

   // mtz updating
   continue_watching_mtz = false;
   continue_watching_coordinates_file = false;

   previous_eye_position = clipper::Coord_orth(-999, -999, -999);

   radial_map_colour_saturation = 0.5;
   radial_map_colour_invert_flag = false;
   radial_map_colour_radius_min =  5.0;
   radial_map_colour_radius_max = 65.0;
   draw_it_for_parallel_plane_restraints = false;
   bonds_rotate_colour_map_flag = false;
   bonds_colour_map_rotation = 0.0;

   contour_level = 0.25;
   xmap_is_diff_map = false;
   have_unit_cell = false;
   save_use_reso_limits = false;
   save_low_reso_limit = 9999.9;
   save_high_reso_limit = 2.0;
   save_is_diff_map_flag = false;
   save_is_anomalous_map_flag = false;
   save_use_weights = false;
   refmac_r_free_flag_sensible = false;
   manual_bond_colour = false;
   map_mean_ = 0.0;
   map_sigma_ = 1.0;

   draw_chiral_volume_outlier_markers_flag = false;
}

std::string
coot::backup_file_info_t::get_timespec_string() const {

   char buffer[80];
   struct tm* timeinfo = localtime(&ctime.tv_sec);
   strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", timeinfo);
   std::ostringstream oss;
   oss << buffer << "." << std::setfill('0') << std::setw(3) << (ctime.tv_nsec / 1000000);
   return oss.str();
}


int
molecule_class_info_t::update_molecule(std::string file_name, std::string cwd) {

   return handle_read_draw_molecule(imol_no,
                                    file_name, cwd,
                                    graphics_info_t::Geom_p(),
                                    0, 0, true, false, bond_width, bonds_box_type, false);
}


void
molecule_class_info_t::set_draw_model_molecule_as_lines(bool state) {

   if (state != draw_model_molecule_as_lines) {
      draw_model_molecule_as_lines = state;
      make_bonds_type_checked();
   }
}

// Return the molecule number of the molecule that we just filled.
// Return -1 if there was a failure.
//
int
molecule_class_info_t::handle_read_draw_molecule(int imol_no_in,
                                                 std::string filename,
                                                 std::string cwd,
                                                 coot::protein_geometry *geom_p,
                                                 short int reset_rotation_centre,
                                                 short int is_undo_or_redo,
                                                 bool allow_duplseqnum,
                                                 bool convert_to_v2_atom_names_flag,
                                                 float bond_width_in,
                                                 int bonds_box_type_in,
                                                 bool warn_about_missing_symmetry_flag) {

   //
   graphics_info_t g;
   imol_no = imol_no_in;

   if (! is_undo_or_redo) {
      bond_width = bond_width_in;
      bonds_box_type = bonds_box_type_in;
      if (g.draw_stick_mode_atoms_default == false) {
         display_stick_mode_atoms_flag = false;
      }
   }

   // std::cout << "DEBUG:: ---- imol_no is now " << imol_no << std::endl;

   // need to check that filename exists and is a file.
   //
   struct stat s;
   int status = stat(filename.c_str(), &s);

   // stat check the link targets not the link itself, lstat stats the
   // link itself.
   //
   if (status != 0 || !S_ISREG (s.st_mode)) {
      std::cout << "WARNING:: Error reading " << filename << std::endl;
      if (S_ISDIR(s.st_mode)) {
         std::cout << filename << " is a directory." << std::endl;
      }
      return -1; // which is status in an error
   }


   // Read in pdb, [shelx files use the read_shelx_ins_file method]
   //
   bool verbose = true;

   input_molecule_was_in_mmcif = false;
   if (coot::is_mmcif_filename(filename))
      input_molecule_was_in_mmcif = true;

   bool use_gemmi = graphics_info_t::use_gemmi;
   atom_sel = get_atom_selection(filename, use_gemmi, allow_duplseqnum, verbose);

   if (atom_sel.read_success == 1) {

      // update the geometry as needed
      if (geom_p) {
         geom_p->read_extra_dictionaries_for_molecule(atom_sel.mol, imol_no,
                                                      &graphics_info_t::cif_dictionary_read_number);
      } else {
         std::cout << "ERROR:: mci::handle_read_draw_molecule(): geom_p is null" << std::endl;
      }

      // LINK info:
      int n_models = atom_sel.mol->GetNumberOfModels();
      // std::cout << "INFO:: Found " << n_models << " models\n";
      logger.log(log_t::INFO, "Found", n_models, " models");
      for (int imod=1; imod<=n_models; imod++) {
         mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
         if (model_p) {
            int n_links = model_p->GetNumberOfLinks();
            // std::cout << "   Model "  << imod << " had " << n_links << " links\n";}
            logger.log(log_t::INFO, "Model", imod, "had", n_links, "links");
         }
      }

      //
      // and move mol_class_info to indexed molecule[n_molecules];
      // note indexing difficulties/considerations.

      // save it in the static molecule_class_info_t
      //

      // mmdb::CMMDBCryst *cryst_p =  (atom_sel.mol)->get_cell_p();

      mmdb::mat44 my_matt;

      //
      //
      int err = atom_sel.mol->GetTMatrix(my_matt, 0, 0, 0, 0);
      if (warn_about_missing_symmetry_flag) {
         if (err != mmdb::SYMOP_Ok) {
            std::cout << "WARNING:: No symmetry available for this molecule"
                      << std::endl;
         }
      }

      // initialize some things.
      //
//       std::cout << "initialize_coordinate_things_on_read_molecule_internal for mol no "
//                 << imol_no << std::endl;
      initialize_coordinate_things_on_read_molecule_internal(filename, is_undo_or_redo);

      set_have_unit_cell_flag_maybe(warn_about_missing_symmetry_flag);

      add_molecular_symmetry_matrices(); // process REMARK 350 BIOMT[123]

      if (molecule_is_all_c_alphas()) {
         bool force_rebonding = true;
         ca_representation(force_rebonding);
      } else {

         if (! is_undo_or_redo) {
            //          std::cout << "DEBUG:: filling ghost info in
            //          handle_read_draw_molecule" << std::endl;
            short int do_rtops_flag = 0;
            // 0.7 is not used (I think) if do_rtops_flag is 0.
            // hack to fix Mac bug/strangeness

            // It only makes sense to do NCS chain searching on
            // crystallographic models.  Which generally only have one
            // model per manager.  This may change in future...
            int nmodels = atom_sel.mol->GetNumberOfModels();
            if (nmodels == 1) {
               fill_ghost_info(do_rtops_flag, 0.7); // returns nghosts
               // std::cout << "INFO:: found " << nghosts << " ghosts\n";
            }
         } else {
            update_mols_in_additional_representations(); //uses atom_sel.mol
         }
         // Generate bonds and save them in the graphical_bonds_container
         // which has static data members.
         //
         if (bonds_box_type == coot::UNSET_TYPE)
            bonds_box_type = coot::NORMAL_BONDS;
         // std::cout << "debug:: ---- handle_read_draw_molecule() calls make_bonds_type_checked()" << std::endl;
         make_bonds_type_checked(__FUNCTION__);
      }

      draw_it = 1;
      if (g.show_symmetry == 1) {
         if (show_symmetry) {  // internal
            update_symmetry();
         }
      }

      // Now, we have no map assocaited with this molecule,
      // so we set the draw_vects to zero.
      //
      // However, in future, we will have maps associated with
      // coordinates, so we should put a test here first before
      // we:
      //

      // debug();

      //
      if (g.recentre_on_read_pdb || imol_no_in == 0)  // n_molecules is updated
                                                      // in c-interface.cc
         if (reset_rotation_centre)
            g.setRotationCentre(::centre_of_molecule(atom_sel));

      // update the maps so that they appear around the new centre.
      //
      if (reset_rotation_centre) {
         for (int ii=0; ii<g.n_molecules(); ii++) {
            g.molecules[ii].update_map(graphics_info_t::auto_recontour_map_flag);
         }
      }

      g.run_post_read_model_hook(imol_no);

      // save state strings
      // 20231228-PE hack in the "coot." for now and remove it later when
      // writing a scheme script (not done ATM as far as I know).
      save_state_command_strings_.push_back("coot.handle-read-draw-molecule");
      std::string f1 = coot::util::intelligent_debackslash(filename);
      std::string f2 = coot::util::relativise_file_name(f1, cwd);
      save_state_command_strings_.push_back(single_quote(f2));

      return 1;

   } else {
      std::cout << "There was a coordinates read error\n";
      return -1;
   }
}

// Cleaner interface to molecule's attributes:
std::pair<bool, clipper::Spacegroup>
molecule_class_info_t::space_group() const {

   clipper::Spacegroup sg;
   std::pair<bool, clipper::Spacegroup> p(0, sg);

   if (has_model()) {
      // I want just the symmetry

      try { // for now
         std::pair<clipper::Cell, clipper::Spacegroup> cell_sg =
            coot::util::get_cell_symm(atom_sel.mol);
         if (!cell_sg.second.is_null()) {
            p.first = 1;
            p.second = cell_sg.second;
         }
      }
      catch (const std::runtime_error &rte) {
         std::cout << "ERROR:: " << rte.what() << std::endl;
      }
   }
   return p;
}

std::pair<bool, clipper::Cell>
molecule_class_info_t::cell() const {

   // NXMAP-FIXME

   clipper::Cell cell;
   std::pair<bool, clipper::Cell> p(0, cell);
   if (has_xmap()) {
      p = std::pair<bool, clipper::Cell> (1, xmap.cell());
   }

   if (has_model()) {
      mmdb::realtype a[6];
      mmdb::realtype vol;
      int orthcode;
      atom_sel.mol->GetCell(a[0], a[1], a[2], a[3], a[4], a[5], vol, orthcode);
      clipper::Cell_descr cdr(a[0], a[1], a[2],
                              clipper::Util::d2rad(a[3]),
                              clipper::Util::d2rad(a[4]),
                              clipper::Util::d2rad(a[5]));
      p.first = 1;
      p.second = clipper::Cell(cdr);
   }
   return p;
}


coot::Cartesian
molecule_class_info_t::centre_of_molecule() const {

   double xs=0, ys=0, zs=0;
   int n_atoms = 0;
   for (int i=0; i<atom_sel.n_selected_atoms; i++) {
      mmdb::realtype x = atom_sel.atom_selection[i]->x;
      mmdb::realtype y = atom_sel.atom_selection[i]->y;
      mmdb::realtype z = atom_sel.atom_selection[i]->z;
      if (x > -9999.9)
         if (x < 9999.9)
            if (y > -9999.9)
               if (y < 9999.9)
                  if (z > -9999.9)
                     if (z < 9999.9) {
                        xs += x;
                        ys += y;
                        zs += z;
                        n_atoms++;
                     }
   }
   if (n_atoms > 0) {
      xs /= double(n_atoms);
      ys /= double(n_atoms);
      zs /= double(n_atoms);
   }
   return coot::Cartesian(xs, ys, zs);
}

float
molecule_class_info_t::size_of_molecule() const {

   double d2_sum = 0;
   float r = 0;
   coot::Cartesian centre = centre_of_molecule();
   if (atom_sel.n_selected_atoms > 0) {
      for (int i=0; i<atom_sel.n_selected_atoms; i++) {
         double d =
            (atom_sel.atom_selection[i]->x - centre.x()) *
            (atom_sel.atom_selection[i]->x - centre.x()) +
            (atom_sel.atom_selection[i]->y - centre.y()) *
            (atom_sel.atom_selection[i]->y - centre.y()) +
            (atom_sel.atom_selection[i]->z - centre.z()) *
            (atom_sel.atom_selection[i]->z - centre.z());
         // std::cout << i << " adding in d = " << d << std::endl;
         d2_sum += d;
      }
      double msd = d2_sum/double(atom_sel.n_selected_atoms);
      r = sqrt(msd);
//       std::cout << " for imol = " << imol_no << " "
//                << r << " = sqrt(" << msd << ") " << " = sqrt("
//                << d2_sum << "/" << atom_sel.n_selected_atoms << ")" << std::endl;
   }
   return r;
}


std::string
molecule_class_info_t::show_spacegroup() const {

   std::string s("No spacegroup");

   if (has_model()) {
      char *st = atom_sel.mol->GetSpaceGroup();
      if (st)
         s = st;
   }

   if (has_xmap())
      s = xmap.spacegroup().symbol_hm();

   return s;
}


coot::at_dist_info_t
molecule_class_info_t::closest_atom(const coot::Cartesian &pt) const {

   coot::at_dist_info_t at_info = closest_atom(pt, true);
   if (at_info.atom)
      return at_info;
   else
      return closest_atom(pt, false);

}

coot::at_dist_info_t
molecule_class_info_t::closest_atom(const coot::Cartesian &pt, bool ca_check_flag) const {
   return closest_atom(pt, ca_check_flag, "", false); // don't use chain-id
}

coot::at_dist_info_t
molecule_class_info_t::closest_atom(const coot::Cartesian &pt, bool ca_check_flag,
                                    const std::string &chain_id,
                                    bool use_this_chain_id) const {

   coot::at_dist_info_t at_info(0,0,0);
   mmdb::Atom *at_best = 0;
   float dist_best = 99999999999.9;

   for (int iat=0; iat<atom_sel.n_selected_atoms; iat++) {
      mmdb::Atom *at = atom_sel.atom_selection[iat];
      if (! at->isTer()) {
         std::string chain_id_from_at(at->GetChainID());
         if ((chain_id_from_at == chain_id) || !use_this_chain_id) {
            float d2 = (at->x - pt.x()) * (at->x - pt.x());
            d2 += (at->y - pt.y()) * (at->y - pt.y());
            d2 += (at->z - pt.z()) * (at->z - pt.z());
            if (d2 < dist_best) {
               dist_best = d2;
               at_best = at;
               // Now, does this at belong to a residue that has a CA?  If
               // it does, reset at_best to be the CA of the residue, but
               // keep dist_best as it was, of course.
               if (ca_check_flag) {
                  mmdb::Residue *res = at->residue;
                  if (res) {
                     int natoms = 0;
                     mmdb::PPAtom residue_atoms = 0;
                     res->GetAtomTable(residue_atoms, natoms);
                     for (int iatom=0; iatom<natoms; iatom++) {
                        mmdb::Atom *r_at = residue_atoms[iatom];
                        if (! r_at->isTer()) {
                           if (! strcmp(r_at->name, " CA ")) {
                              if (! strcmp(r_at->altLoc, at->altLoc)) {
                                 at_best = r_at;
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
   if (at_best) {
      at_info.dist = sqrt(dist_best);
      at_info.atom = at_best;
      at_info.imol = imol_no;
   }
   return at_info;
}




std::string
molecule_class_info_t::single_quote(const std::string &s) const {
   std::string r("\"");
   r += s;
   r += "\"";
   return r;
}

// needs show_symmetry flag also
//
void
molecule_class_info_t::install_model(int imol_no_in,
                                     atom_selection_container_t asc,
                                     const coot::protein_geometry *geom_p,
                                     const std::string &name,
                                     short int display_in_display_control_widget_status,
                                     bool is_from_shelx_ins, // default false
                                     bool warn_about_missing_symmetry_flag // default false
                                     ) {

   bool generate_ghost_info = true;
   std::vector<coot::ghost_molecule_display_t> dummy_ghosts;
   install_model_with_ghosts(imol_no_in, asc, geom_p, name,
                             display_in_display_control_widget_status,
                             dummy_ghosts,
                             is_from_shelx_ins,
                             warn_about_missing_symmetry_flag,
                             generate_ghost_info);
}

// if generate_ghost_info is false, copy the ghost info, don't regnerate it.
// otherwise, behave as you used to.
void
molecule_class_info_t::install_model_with_ghosts(int imol_no_in,
                                                 atom_selection_container_t asc,
                                                 const coot::protein_geometry *geom_p,
                                                 const std::string &name,
                                                 short int display_in_display_control_widget_status,
                                                 const std::vector<coot::ghost_molecule_display_t> &ncs_ghosts_in,
                                                 bool is_from_shelx_ins, // default false
                                                 bool warn_about_missing_symmetry_flag, // default false
                                                 bool generate_ghost_info
                                                 ) {

   imol_no = imol_no_in;
   graphics_info_t g;  // pass g.Geom_p()
   bond_width = g.default_bond_width; // bleugh, perhaps this should
                                      // be a passed parameter?
   is_from_shelx_ins_flag = is_from_shelx_ins;

   atom_sel = asc;

   mmdb::mat44 my_matt;

   int err = asc.mol->GetTMatrix(my_matt, 0, 0, 0, 0);
   if (warn_about_missing_symmetry_flag)
      if (err != 0)
         std::cout << "WARNING:: No symmetry available for this molecule" << std::endl;
   set_have_unit_cell_flag_maybe(warn_about_missing_symmetry_flag);

   std::set<int> dummy;
   makebonds(geom_p, dummy);
   if (g.show_symmetry == 1)
      if (show_symmetry)
         update_symmetry();

   have_unsaved_changes_flag = 1;

   short int is_undo_or_redo = 0;

   if (display_in_display_control_widget_status == 0) {
      // treat as undo then (e.g. "terminal residue") made by add_cb_to_terminal_res().
      is_undo_or_redo = 1;
   } else {
      pickable_atom_selection = 1;
   }

   if (generate_ghost_info) {
      bool do_rtops = true;
      fill_ghost_info(do_rtops, graphics_info_t::ncs_homology_level);
   } else {
      // ncs_ghosts = ncs_ghosts_in;
      ncs_ghosts.clear();
      for (unsigned int i=0; i<ncs_ghosts_in.size(); i++)
         ncs_ghosts.push_back(drawn_ghost_molecule_display_t(ncs_ghosts_in[i]));
   }
   initialize_coordinate_things_on_read_molecule_internal(name, is_undo_or_redo);
}

void
molecule_class_info_t::install_model(int imol_no_in,
                                     mmdb::Manager *mol,
                                     const coot::protein_geometry *geom_p,
                                     const std::string &mol_name,
                                     short int display_in_display_control_widget_status,
                                     bool is_from_shelx_ins, // default false
                                     bool warn_about_missing_symmetry_flag // default false
                                     ) {

   atom_selection_container_t asc = make_asc(mol);
   install_model(imol_no_in, asc, geom_p, mol_name, display_in_display_control_widget_status, is_from_shelx_ins, warn_about_missing_symmetry_flag);
}



void
molecule_class_info_t::draw_atom_labels(int brief_atom_labels_flag,
                                        short int seg_ids_in_atom_labels_flag,
                                        const glm::vec4 &atom_label_colour,
                                        stereo_eye_t eye,
                                        const glm::mat4 &mvp,
                                        const glm::mat4 &view_rotation) {


   if (draw_it) {
      if (has_model()) {

         // keep a list of atoms that have labels (either in graphics_info
         // or mol_class_info) and loop over them calling label_atom(i)
         // which labels the i'th atom of the atom selection in
         // mol_class_info.
         //
         //
         //
         int n_atoms_to_label = labelled_atom_index_list.size();

         // also remove labels from atom indexes list of over the end.
         for (int ii=0; ii<n_atoms_to_label ; ii++)
            draw_atom_label(labelled_atom_index_list[ii], brief_atom_labels_flag,
                            seg_ids_in_atom_labels_flag, atom_label_colour,
                            eye, mvp, view_rotation);

         unsigned int n_symm_atoms_to_label = labelled_symm_atom_index_list.size();

         glm::vec4 blueish(0.7, 0.7, 1.0, 1.0); // this could/should be user-settable.
         for (unsigned int ii=0; ii<n_symm_atoms_to_label ; ii++) {
            const auto &st = labelled_symm_atom_symm_trans_[ii];
            draw_symm_atom_label(labelled_symm_atom_index_list[ii], st,
                            blueish, mvp, view_rotation);

         }
      }
   }
}

void
molecule_class_info_t::trim_atom_label_table() {

   int new_max_atom_index = atom_sel.n_selected_atoms;

   // 20120721 Yes, it was.  There is a better way, using a functor.
   //
//    // Man, this is a clumsy way of removing specific ints from a vector<int>.
//    // It's needed because as we do an erase, the labelled_atom_index_list changes
//    // so that we were double erasing an element:
//    // erasing *it: 141 limit: 135
//    // erasing *it: 141 limit: 135
//    //    -> crash.
//    // So do one at a time.
//    //
//    // Perhaps it would be better to use a queue or construct a new
//    // vector<int> and reassign labelled_atom_index_list at the end.
//    //
//    // Anyway, this seems to work as it is, so I'll leave it for now.
//    //

   // modern
   //
   labelled_atom_index_list.erase(std::remove_if(labelled_atom_index_list.begin(),
                                                 labelled_atom_index_list.end(),
                                                 labelled_atom_remover(new_max_atom_index)),
                                  labelled_atom_index_list.end());

   labelled_symm_atom_index_list.erase(std::remove_if(labelled_symm_atom_index_list.begin(),
                                                      labelled_symm_atom_index_list.end(),
                                                      labelled_atom_remover(new_max_atom_index)),
                                       labelled_symm_atom_index_list.end());

   // nice macro?
   // erase_items(labelled_atom_index_list, labelled_atom_remover_test(new_max_atom_index));
}


void
molecule_class_info_t::old_draw_anisotropic_atoms() {
#if 0
   int c; // atom colour

   if (has_model()) {
      graphics_info_t g; // bleugh
      bool is_bb = g.background_is_black_p(); //

      if (draw_it) {
         if (g.show_aniso_atoms_flag == 1 ) {
            glPushMatrix();

            float rx = g.X();
            float ry = g.Y();
            float rz = g.Z();

            float x1, y1, z1;
            float x_diff, y_diff, z_diff;
            float d2, mc_r2 = g.show_aniso_atoms_radius*g.show_aniso_atoms_radius;
            float r;

            for (int i=0; i<atom_sel.n_selected_atoms; i++) {

               // put a wiresphere at the atom positions

               if (atom_sel.atom_selection[i]->u11 > 0) {

                  mmdb::Atom *at = atom_sel.atom_selection[i];
                  std::string ele(at->element);

                  // if (draw_hydrogens_flag || ! mmdb_utils::is_hydrogen(ele))
                  if (draw_hydrogens_flag || ele != " H") {

                     glLineWidth(1.0);
                     glPushMatrix();

                     x1 = atom_sel.atom_selection[i]->x;
                     y1 = atom_sel.atom_selection[i]->y;
                     z1 = atom_sel.atom_selection[i]->z;

                     x_diff = x1 - rx;
                     y_diff = y1 - ry;
                     z_diff = z1 - rz;

                     d2 = x_diff*x_diff + y_diff*y_diff + z_diff*z_diff;

                     // are we either inside the distance or there is no distance set?
                     //
                     if ( (d2 <= mc_r2) || (g.show_aniso_atoms_radius_flag == 0) ) {

                        c = get_atom_colour_from_element(ele);
                        set_bond_colour_by_mol_no(c, is_bb);

                        GL_matrix mat(atom_sel.atom_selection[i]->u11,
                                      atom_sel.atom_selection[i]->u12,
                                      atom_sel.atom_selection[i]->u13,
                                      atom_sel.atom_selection[i]->u12,
                                      atom_sel.atom_selection[i]->u22,
                                      atom_sel.atom_selection[i]->u23,
                                      atom_sel.atom_selection[i]->u13,
                                      atom_sel.atom_selection[i]->u23,
                                      atom_sel.atom_selection[i]->u33);

                        glTranslatef(x1, y1, z1);
                        // glMultMatrixf(mat.get());
                        // std::cout << "Atom Us " << std::endl;
                        // mat.print_matrix();
                        // std::cout << "Choleskied: " << std::endl;
                        // mat.cholesky().print_matrix();
                        // std::pair<bool,GL_matrix> chol_pair = mat.cholesky();
                        std::pair<bool,GL_matrix> chol_pair = mat.eigensystem();
                        if (chol_pair.first) {
                           glMultMatrixf(chol_pair.second.get());
                           r = prob_to_radius(g.show_aniso_atoms_probability);
                           // note: g.show_aniso_atoms_probability is in the range
                           // 0.0 -> 100.0
                           // glutWireSphere(r, 10, 10);
                           std::cout << "Fix wire sphere\n";
                        } else {
                           std::cout << "Bad Anistropic Us for " << atom_sel.atom_selection[i]
                                     << std::endl;
                        }
                     }
                     glPopMatrix();
                  }
               }
            }
            glPopMatrix();
         }
      }
   }
#endif
}

coot::colour_t
molecule_class_info_t::get_bond_colour_basic(int colour_index, bool against_a_dark_background) const {

   // this is used for intermediate atoms

   coot::colour_t col(0.5, 0.5, 0.5);

   if (true) {  // background check
      switch (colour_index) {
      case CARBON_BOND:
         col = coot::colour_t (0.2, 0.7, 0.1);
         break;
      case GREEN_BOND:
         col = coot::colour_t (0.0, 0.7, 0.0);
         break;
      case BLUE_BOND:
         col = coot::colour_t (0.2, 0.2, 0.8);
         break;
      case RED_BOND:
         col = coot::colour_t (0.8, 0.1, 0.1);
         break;
      case YELLOW_BOND:
         col = coot::colour_t (0.7, 0.7, 0.0);
         break;
      case GREY_BOND:
         col = coot::colour_t (0.5, 0.5, 0.5);
         break;
      case HYDROGEN_GREY_BOND:
         col = coot::colour_t (0.7, 0.7, 0.7);
         break;
      case DEUTERIUM_PINK:
         col = coot::colour_t (0.8, 0.6, 0.64);
         break;
      case MAGENTA_BOND:
         col = coot::colour_t (0.8, 0.1, 0.8);
         break;
      case DARK_GREEN_BOND:
         col = coot::colour_t (0.05, 0.69, 0.05);
         break;
      case DARK_ORANGE_BOND:
         col = coot::colour_t (0.7, 0.7, 0.05);
         break;
      case DARK_BROWN_BOND:
         col = coot::colour_t (0.5, 0.5, 0.1);
         break;
      case VIOLET:
         col = coot::colour_t(0.93, 0.51, 0.93);
         break;
      case DARK_VIOLET:
         col = coot::colour_t(0.58, 0.0, 0.83);
         break;
      case BORON_PINK:
         col = coot::colour_t(0.98, 0.78, 0.69); // 0.98 0.72 0.63
         break;
      default:
         col = coot::colour_t (0.7, 0.8, 0.8);
      }
   }
   return col;
}

coot::colour_t
molecule_class_info_t::get_bond_colour_by_mol_no(int colour_index, bool against_a_dark_background) const {

   // No OpenGL here now.

   //std::cout << "get_bond_colour_by_mol_no() " << colour_index << " " << against_a_dark_background << std::endl;
   // GLenum err = glGetError(); if (err) std::cout << "GL status in get_bond_colour_by_mol_no() --start-- " << err << std::endl;

   coot::colour_t rgb;

   if (bonds_rotate_colour_map_flag == 0) {
      rgb = get_bond_colour_basic(colour_index, against_a_dark_background);
   } else {

      float rotation_size = bonds_colour_map_rotation/360.0;

      // rotation_size typically then: 2*32/360 = 0.178

      // This is for colour-by-chain-carbons-only carbon colour

      if (colour_index >= 50) {
         float ii_f = colour_index - 50;
         ii_f += 1.2 * static_cast<float>(imol_no);
         if (against_a_dark_background) {
            rgb[0] = 0.75; rgb[1] = 0.55; rgb[2] = 0.45; // pale/cream
         } else {
            rgb[0] = 0.55; rgb[1] = 0.4; rgb[2] = 0.3;
         }
         float ra = ii_f*79.0/360.0;
         ra += rotation_size;
         while (ra > 1.0) ra -= 1.0;
         if (ra > 0) {
            rgb.rotate(ra);
         }
         // std::cout << "get_bond_colour_by_mol_no() get chain colour for colour_index "
         // << colour_index << " " << rgb << std::endl;
      } else {

         while (rotation_size > 1.0) { // no more black bonds?
            rotation_size -= 1.0;
         }

         if (against_a_dark_background) {

            if (false)
               std::cout << "get_bond_colour_by_mol_no() against_a_dark_background==true, idx: " << colour_index << " vs "
                         << " green "   << GREEN_BOND << " "
                         << " blue "    << BLUE_BOND << " "
                         << " red "     << RED_BOND << " "
                         << " yellow "  << YELLOW_BOND << " "
                         << " grey "    << GREY_BOND << " "
                         << " H-grey "  << HYDROGEN_GREY_BOND << " "
                         << " magenta " << MAGENTA_BOND << " "
                         << std::endl;

            switch (colour_index) {
            case CARBON_BOND:
               if (use_bespoke_grey_colour_for_carbon_atoms) {
                  rgb = bespoke_carbon_atoms_colour;
               } else {
                  rgb[0] = 0.7; rgb[1] =  0.7; rgb[2] =  0.2;
               }
               break;
            case YELLOW_BOND:
               rgb[0] = 0.6; rgb[1] = 0.98; rgb[2] =  0.2;
               break;
            case BLUE_BOND:
               rgb[0] = 0.25; rgb[1] =  0.25; rgb[2] = 1.0;
               break;
            case RED_BOND:
               rgb[0] = 0.9; rgb[1] =  0.0; rgb[2] =  0.0;
               break;
            case GREEN_BOND:
               rgb[0] = 0.2; rgb[1] =  0.9; rgb[2] =  0.2;
               break;
            case GREY_BOND:
               rgb[0] = 0.6; rgb[1] =  0.6; rgb[2] =  0.6;
               break;
            case HYDROGEN_GREY_BOND:
               rgb[0] = 0.7; rgb[1] =  0.7; rgb[2] =  0.7;
               break;
            case DEUTERIUM_PINK:
               rgb[0] = 0.8; rgb[1] =  0.6; rgb[2] =  0.64;
               break;
               // replaced in mmdb-extras.h
               //       case white:
               //          rgb[0] = 0.99; rgb[1] =  0.99; rgb[2] = 0.99;
               //          break;
            case MAGENTA_BOND:
               rgb[0] = 0.99; rgb[1] =  0.2; rgb[2] = 0.99;
               break;
            case ORANGE_BOND:
               rgb[0] = 0.89; rgb[1] =  0.89; rgb[2] = 0.1;
               break;
            case CYAN_BOND:
               rgb[0] = 0.1; rgb[1] =  0.89; rgb[2] = 0.89;
               break;
            case DARK_GREEN_BOND:
               rgb[0] = 0.05; rgb[1] =  0.5; rgb[2] =  0.05;
               break;
            case DARK_ORANGE_BOND:
               rgb[0] = 0.7; rgb[1] =  0.7; rgb[2] = 0.05;
               break;
            case DARK_BROWN_BOND:
               rgb[0] = 0.5; rgb[1] =  0.5; rgb[2] = 0.1;
               break;
            case VIOLET:
               rgb[0] = 0.93; rgb[1] = 0.51; rgb[2] = 0.93;
               break;
            case DARK_VIOLET:
               rgb[0] = 0.58; rgb[1] = 0.0; rgb[2] = 0.83;
               break;
            case BORON_PINK:
               rgb[0] = 0.98; rgb[1] = 0.78; rgb[2] = 0.69;
               break;
            default:
               rgb[0] = 0.8; rgb[1] =  0.2; rgb[2] =  0.2;
               rgb.rotate(colour_index*26.0/360.0);
            }

         } else {

            // against a white background.  Less pale (more saturated) and darker.

            switch (colour_index) {
            case CARBON_BOND:
               if (use_bespoke_grey_colour_for_carbon_atoms) {
                  rgb = bespoke_carbon_atoms_colour;
               } else {
                  rgb[0] = 0.6; rgb[1] =  0.7; rgb[2] =  0.05;
               }
               break;
            case YELLOW_BOND:
               rgb[0] = 0.7; rgb[1] =  0.7; rgb[2] =  0.0;
               break;
            case BLUE_BOND:
               rgb[0] = 0.1; rgb[1] =  0.1; rgb[2] =  0.6;
               break;
            case RED_BOND:
               rgb[0] = 0.6; rgb[1] =  0.1; rgb[2] =  0.075; // more tomatoey.
               break;
            case GREEN_BOND:
               rgb[0] = 0.05; rgb[1] =  0.6; rgb[2] =  0.05;
               break;
            case GREY_BOND:
               rgb[0] = 0.5; rgb[1] =  0.5; rgb[2] =  0.5;
               break;
            case HYDROGEN_GREY_BOND:
               rgb[0] = 0.6; rgb[1] =  0.6; rgb[2] =  0.6;
               break;
            case DEUTERIUM_PINK:
               rgb[0] = 0.8; rgb[1] =  0.6; rgb[2] =  0.64;
               break;
            case MAGENTA_BOND:
               rgb[0] = 0.5; rgb[1] =  0.1; rgb[2] = 0.5;
               break;
            case ORANGE_BOND:
               rgb[0] = 0.5; rgb[1] =  0.5; rgb[2] = 0.1;
               break;
            case CYAN_BOND:
               rgb[0] = 0.1; rgb[1] =  0.5; rgb[2] = 0.5;
               break;
            case DARK_GREEN_BOND:
               rgb[0] = 0.05; rgb[1] =  0.69; rgb[2] =  0.05;
               break;
            case DARK_ORANGE_BOND:
               rgb[0] = 0.7; rgb[1] =  0.7; rgb[2] = 0.05;
               break;
            case DARK_BROWN_BOND:
               rgb[0] = 0.5; rgb[1] =  0.5; rgb[2] = 0.1;
               break;
            case VIOLET:
               rgb[0] = 0.93; rgb[1] = 0.51; rgb[2] = 0.93;
               break;
            case DARK_VIOLET:
               rgb[0] = 0.58; rgb[1] = 0.0; rgb[2] = 0.83;
               break;
            case BORON_PINK:
               rgb[0] = 0.98; rgb[1] = 0.78; rgb[2] = 0.69;
               break;

            default:
               rgb[0] = 0.5; rgb[1] =  0.1; rgb[2] =  0.1;
               // rgb = rotate_rgb(rgb, float(i*26.0/360.0));
               rgb.rotate(colour_index*26.0/360.0);
            }
         }

         // "correct" for the +1 added in the calculation of the rotation
         // size.
         // 21. is the default colour map rotation

         rgb.rotate(float(1.0 - 21.0/360.0));

         if (graphics_info_t::rotate_colour_map_on_read_pdb_c_only_flag) {

            if (colour_index == CARBON_BOND) {
               if (use_bespoke_grey_colour_for_carbon_atoms) {
                  rgb = bespoke_carbon_atoms_colour;
               } else {
                  rgb.rotate(rotation_size);
               }
            }
         } else {
            rgb.rotate(rotation_size);
         }
      }
   }
   if (false)
      std::cout << "DEBUG:: returning colour " << rgb << std::endl;
   return rgb;
}


// fix the name to something involving rotation perhaps?
//
// not const because bond_colour_internal is set.
void
molecule_class_info_t::set_bond_colour_by_mol_no(int colour_index, bool against_a_dark_background) {

   coot::colour_t col = get_bond_colour_by_mol_no(colour_index, against_a_dark_background);
   glColor3f(col.col[0], col.col[1], col.col[2]);
   // std::vector<float> bond_colour_internal;
   bond_colour_internal = {col.col[0], col.col[1], col.col[2]};
}

void
molecule_class_info_t::set_bond_colour_for_goodsell_mode(int icol, bool against_a_dark_background) {

   // 20221114-PE this function seems not to be used, in fact!

   bool is_C = !(icol %2);
   int n_steps = icol/2;

   coot::colour_t col(0.9, 0.52, 0.52);
   if (is_C) col = coot::colour_t(0.82, 0.6, 0.6); // more pastel

   col.rotate(0.06 * n_steps);

   glColor3f(col[0], col[1], col[2]);

}


#include "colour-functions.hh"

// aka rainbow - or maybe b factor, occupancy or user defined colour index
//
void
molecule_class_info_t::set_bond_colour_by_colour_wheel_position(int i, int bonds_box_type) {

   std::vector<float> rgb(3);
   rgb[0] = 0.2f; rgb[1] =  0.2f; rgb[2] =  0.8f; // blue

   bool done = false;
   int offset = 0; // blue starts at 0

   if (bonds_box_type == coot::COLOUR_BY_USER_DEFINED_COLOURS____BONDS || bonds_box_type == coot::COLOUR_BY_USER_DEFINED_COLOURS_CA_BONDS) {
      if (i == 0) {
         rgb[0] = 0.8f; rgb[1] =  0.8f; rgb[2] =  0.8f; // white
         done = true;
      }
      if (i == 1) {
         rgb[0] = 0.3f; rgb[1] =  0.3f; rgb[2] =  0.3f; // dark-grey
         done = true;
      }
      offset=2; // blue starts at 2.
   }

   if (false)
      std::cout << "debug set_bond_colour_by_colour_wheel_position() " << i
                << " box_type " << bonds_box_type << " vs " << coot::COLOUR_BY_USER_DEFINED_COLOURS____BONDS
                << std::endl;

   if (bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS_B_FACTOR_COLOUR || bonds_box_type == coot::COLOUR_BY_B_FACTOR_BONDS) {
      rgb[0] = 0.3f; rgb[1] =  0.3f; rgb[2] =  0.95f;
      const unsigned int n_b_factor_colours = 48; // matches index_for_b_factor() in my_atom_colour_map_t
      float f = static_cast<float>(i)/static_cast<float>(n_b_factor_colours);
      // f is in the range 0 to 1
      const float pi = 3.1415926535;
      float rotation_size = -0.11 * f * 2.0  * pi;
      if (rotation_size < -0.6666) rotation_size = -0.66666; // otherwise black bonds
      // std::cout << "rotation_size: " << rotation_size << std::endl;
      rgb = rotate_rgb(rgb, rotation_size);
      done = true;
   }
   if (! done) {
      float max_colour = 30;

      // 30 is the size of rainbow colours, 0 -> 1.0 is the range of rainbow colours

      float rotation_size = 1.0 - float(i-offset) * 0.7/max_colour + bonds_colour_map_rotation/360.0;
      rgb = rotate_rgb(rgb, rotation_size);
   }

   // rotation_size size has useful colours between
   // 1.0 (or higher?) and 0.0.
   // Below 0.0, to -0.65 (more than -0.68) (dependant on starting rgb
   // values I guess) there is amusing colour continuation of the colour
   // wheel (red through purple to blue and further).
   //
   if (false)
      std::cout << "set_bond_colour_by_colour_wheel_position "  << i << " " << " "
                << rgb[0] << " " << rgb[1] << " " << rgb[2] << " " << std::endl;
   bond_colour_internal = rgb;
   // glColor3f(rgb[0], rgb[1], rgb[2]); old.
}

// 20220214-PE modern graphics
glm::vec4
molecule_class_info_t::get_bond_colour_by_colour_wheel_position(int icol, int bonds_box_type) const {

   std::vector<float> rgb(3);
   rgb[0] = 0.2f; rgb[1] =  0.2f; rgb[2] =  0.8f; // blue

   bool done = false;
   int offset = 0; // blue starts at 0

   if (bonds_box_type == coot::COLOUR_BY_USER_DEFINED_COLOURS____BONDS) {
      if (icol == 0) {
         rgb[0] = 0.8f; rgb[1] =  0.8f; rgb[2] =  0.8f; // white
         done = true;
      }
      if (icol == 1) {
         rgb[0] = 0.3f; rgb[1] =  0.3f; rgb[2] =  0.3f; // dark-grey
         done = true;
      }
      offset=2; // blue starts at 2.
   }

   if (false)
      std::cout << "debug set_bond_colour_by_colour_wheel_position() " << icol
                << " box_type " << bonds_box_type << " vs " << coot::COLOUR_BY_USER_DEFINED_COLOURS____BONDS
                << std::endl;

   if (bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS_B_FACTOR_COLOUR || bonds_box_type == coot::COLOUR_BY_B_FACTOR_BONDS) {
      rgb[0] = 0.3f; rgb[1] =  0.3f; rgb[2] =  0.95f;
      const unsigned int n_b_factor_colours = 48; // matches index_for_b_factor() in my_atom_colour_map_t
      float f = static_cast<float>(icol)/static_cast<float>(n_b_factor_colours);
      // f is in the range 0 to 1
      const float pi = 3.1415926535;
      float rotation_size = -0.11 * f * 2.0  * pi;
      if (rotation_size < -0.6666) rotation_size = -0.66666; // otherwise black bonds
      // std::cout << "rotation_size: " << rotation_size << std::endl;
      rgb = rotate_rgb(rgb, rotation_size);
      done = true;
   }
   if (! done) {
      float max_colour = 30;

      // 30 is the size of rainbow colours, 0 -> 1.0 is the range of rainbow colours

      if (bonds_box_type == coot::COLOUR_BY_RAINBOW_BONDS) {
         // rotation_size is a fraction of a circle
         float rotation_size = 1.0 - float(icol-offset) * 0.7/max_colour;
         rgb = rotate_rgb(rgb, rotation_size);
         if (false)
            std::cout << "icol " << std::setw(2) <<  icol << " rotation size " << std::setw(7) << rotation_size
                      << "  rgb " << std::setw(6) << rgb[0] << " " << std::setw(6) << rgb[1] << " "
                      << std::setw(6) << rgb[2] << std::endl;
      } else {
	 // std::cout << "this fallback block...." << std::endl;
         float rotation_size = 1.0 - float(icol-offset) * 0.7/max_colour + bonds_colour_map_rotation/360.0;
         rgb = rotate_rgb(rgb, rotation_size);
      }
   }

   // rotation_size size has useful colours between
   // 1.0 (or higher?) and 0.0.
   // Below 0.0, to -0.65 (more than -0.68) (dependant on starting rgb
   // values I guess) there is amusing colour continuation of the colour
   // wheel (red through purple to blue and further).
   //
   if (false)
      std::cout << "set_bond_colour_by_colour_wheel_position "  << icol << " " << " "
                << rgb[0] << " " << rgb[1] << " " << rgb[2] << " " << std::endl;


   glm::vec4 colour(rgb[0], rgb[1], rgb[2], 1.0f);
   return colour;
}


// We find a box (symm: 2 trans: 0 0 0), but we don't find any atoms
// in it, but we still get though to bonds,
// bonds.make_graphical_symmetry_bonds
//
// But there are no bonds generated.
//
void
molecule_class_info_t::update_symmetry() {

   graphics_info_t g;

   // a bit of a hack...
   int shift_search_size = g.symmetry_shift_search_size;

   // std::cout << "DEBUG:: ---- update_symmetry start ----- " << std::endl;

   if ((graphics_info_t::show_symmetry == 1) && (show_symmetry == 1)) {

      // don't do stuff until we have read in a model molecule.
      //
      if (draw_it == 1) {

         molecule_extents_t extents(atom_sel, g.symmetry_search_radius);
         coot::Cartesian point = g.RotationCentre();

         // cout << "extents " << extents << endl;
         // cout << "point:  " << point << endl;
         std::vector<std::pair<symm_trans_t, Cell_Translation> > symm_trans_boxes =
            extents.which_boxes(point, atom_sel, shift_search_size);

         if (false)  {
            std::cout << "DEBUG:: imol_no " << imol_no << " symm_trans_boxes.size() is "
                      << symm_trans_boxes.size() << std::endl;
            std::cout << "Here are the symms we should check:" << std::endl;
            for(unsigned int ii=0; ii<symm_trans_boxes.size(); ii++)
               std::cout << ii << " " << symm_trans_boxes[ii].first << " "
                        << symm_trans_boxes[ii].second << std::endl;
         }

         if (symm_trans_boxes.size() > 0) {

            // when bonds goes out of scope (i.e. immediate after
            // this) then the class data member vector "bonds" of the
            // Bond_lines_container gets given back.
            //
            // It is with the "new"ly allocated graphical_symmetry_bonds
            // that we need to concern ourselves.
            //
            Bond_lines_container bonds;

            //
            // delete the old symmetry_bonds_box
            //
            // symmetry_bonds_box.clear_up();
            clear_up_all_symmetry();
            symmetry_bonds_box.clear();

            //     for (unsigned int ibox=0; ibox<symm_trans_boxes.size(); ibox++)
            //        std::cout << "box " << ibox << "/" << symm_trans_boxes.size()
            //        << " " << symm_trans_boxes[ibox] << "\n";

            bool do_intermolecular_symmetry_bonds = false; // for now

            symmetry_bonds_box =
               bonds.addSymmetry_vector_symms(atom_sel, imol_no,
                                              point,
                                              graphics_info_t::symmetry_search_radius,
                                              symm_trans_boxes,
                                              symmetry_as_calphas,
                                              symmetry_whole_chain_flag,
                                              draw_hydrogens_flag,
                                              do_intermolecular_symmetry_bonds);

            make_glsl_symmetry_bonds();
            this_molecule_has_crystallographic_symmetry = true;

         } else {
            Bond_lines_container bonds(NO_SYMMETRY_BONDS);
         }

         if (false) // come back and debug crystallographic strict NCS one day!
            std::cout << "Here in imol " << imol_no << " with show_strict_ncs_flag " << show_strict_ncs_flag
                      << " and  strict_ncs_matrices size " << strict_ncs_matrices.size() << std::endl;
         if (show_strict_ncs_flag == 1) {
            if (! strict_ncs_matrices.empty()) {
               update_strict_ncs_symmetry(point, extents);
            }
         }

      } else {
         // cout << "update_symmetry: no molecule yet" << endl;
      }
   }
}

void
molecule_class_info_t::draw_extra_restraints_representation() {

   std::cout << "old code in draw_extra_restraints_representation() " << std::endl;
#if 0
   if (draw_it) {
      if (draw_it_for_extra_restraints) {
         if (extra_restraints_representation.bonds.size() > 0) {
            glLineWidth(1.0);
            if (extra_restraints_representation_for_bonds_go_to_CA) {
               glLineWidth(3.0);
            } else {
            }
            glColor3f(0.6, 0.6, 0.8);

            glBegin(GL_LINES);
            for (unsigned int ib=0; ib<extra_restraints_representation.bonds.size(); ib++) {

               const coot::extra_restraints_representation_t::extra_bond_restraints_respresentation_t &res =
                  extra_restraints_representation.bonds[ib];

               // red if actual distance is greater than target
               //
               double d_sqd = (res.second - res.first).clipper::Coord_orth::lengthsq();

               if (res.esd > 0) {
                  double nz = (sqrt(d_sqd) - res.target_dist)/res.esd;

                  /*
                    std::cout << "debug:: nz " << nz << " target " << res.target_dist << " model "
                              << sqrt(d_sqd) << " esd " << res.esd << std::endl;
                  */

                  // we want to make short be green and long be purple
                  float b_2 = 0.05 * nz;
                  if (b_2 >  0.4999) b_2 =  0.4999;
                  if (b_2 < -0.4999) b_2 = -0.4999;
                  // b_2 is now between -0.5 and +0.5
                  float r = 0.5 - b_2;
                  float g = 0.5 + b_2;
                  float b = 0.5 - b_2;
                  glColor3f(r, g, b);
               }
               glVertex3f(res.first.x(), res.first.y(), res.first.z());
               glVertex3f(res.second.x(), res.second.y(), res.second.z());
            }
            glEnd();
         }
      }
   }

   draw_parallel_plane_restraints_representation();
#endif
}

void
molecule_class_info_t::draw_parallel_plane_restraints_representation() {

   std::cout << "old code in draw_parallel_plane_restraints_representation() " << std::endl;

#if 0
   if (draw_it) {
      if (draw_it_for_extra_restraints) {
         if (extra_restraints_representation.parallel_planes.size() > 0) {
            glLineWidth(2.0);
            glColor3f(0.55, 0.55, 0.3);
            glBegin(GL_LINES);
            for (unsigned int i=0; i<extra_restraints_representation.parallel_planes.size(); i++) {
               const coot::extra_restraints_representation_t::extra_parallel_planes_restraints_representation_t &r =
                  extra_restraints_representation.parallel_planes[i];

               clipper::Coord_orth arb(0.2, 0.8, 0.1);
               clipper::Coord_orth cr(clipper::Coord_orth::cross(r.normal, arb).unit());
               clipper::Coord_orth first_pt = r.ring_centre + r.ring_radius * cr;
               clipper::Coord_orth first_pt_pp = r.plane_projection_point + r.pp_radius * cr;
               // std::cout << i << " r.plane_projection_point: " << r.plane_projection_point.format() << std::endl;

               unsigned int n_steps = 32;
               double step_frac = 1/double(n_steps);
               clipper::Coord_orth pt_1;
               clipper::Coord_orth pt_2;
               for (unsigned int istep=0; istep<n_steps; istep++) {
                  double angle_1 = step_frac * 2.0 * M_PI * istep;
                  double angle_2 = step_frac * 2.0 * M_PI * (istep + 1);
                  pt_1 = coot::util::rotate_around_vector(r.normal, first_pt, r.ring_centre, angle_1);
                  pt_2 = coot::util::rotate_around_vector(r.normal, first_pt, r.ring_centre, angle_2);
                  glVertex3f(pt_1.x(), pt_1.y(), pt_1.z());
                  glVertex3f(pt_2.x(), pt_2.y(), pt_2.z());
               }

               n_steps = 16;
               step_frac = 1/double(n_steps);
               for (unsigned int istep=0; istep<n_steps; istep++) {
                  double angle_1 = step_frac * 2.0 * M_PI * istep;
                  double angle_2 = step_frac * 2.0 * M_PI * (istep + 1);
                  pt_1 = coot::util::rotate_around_vector(r.normal, first_pt_pp, r.plane_projection_point, angle_1);
                  pt_2 = coot::util::rotate_around_vector(r.normal, first_pt_pp, r.plane_projection_point, angle_2);
                  glVertex3f(pt_1.x(), pt_1.y(), pt_1.z());
                  glVertex3f(pt_2.x(), pt_2.y(), pt_2.z());
               }

               // now the lines between planes
               glVertex3d(r.ring_centre.x(), r.ring_centre.y(), r.ring_centre.z());
               glVertex3d(r.plane_projection_point.x(), r.plane_projection_point.y(), r.plane_projection_point.z());
            }
            glEnd();
         }

         // points
         float zsc = graphics_info_t::zoom;
         glPointSize(120.0/zsc);
         glBegin(GL_POINTS);
         for (unsigned int i=0; i<extra_restraints_representation.parallel_planes.size(); i++) {
            const coot::extra_restraints_representation_t::extra_parallel_planes_restraints_representation_t &r =
               extra_restraints_representation.parallel_planes[i];
            glVertex3d(r.plane_projection_point.x(), r.plane_projection_point.y(), r.plane_projection_point.z());
         }
         glEnd();
      }
   }
#endif
}


void
molecule_class_info_t::set_show_unit_cell(bool state) {

#ifndef EMSCRIPTEN
   if (state)
      setup_unit_cell();
#endif
   show_unit_cell_flag = state;

}

#ifndef EMSCRIPTEN
void
molecule_class_info_t::setup_unit_cell() {

   // modify the reference
   auto setup = [] (LinesMesh &lines_mesh_for_cell,
                    const clipper::Cell &cell) {
                   lines_mesh_for_cell = LinesMesh(cell);
                   lines_mesh_for_cell.setup();
                };

   if (lines_mesh_for_cell.empty()) {

      gtk_gl_area_attach_buffers(GTK_GL_AREA(graphics_info_t::glareas[0]));

      if (atom_sel.mol) {
         mmdb::realtype mmdb_cell[6];
         mmdb::realtype vol;
         int orthcode;
         atom_sel.mol->GetCell(mmdb_cell[0], mmdb_cell[1], mmdb_cell[2],
                               mmdb_cell[3], mmdb_cell[4], mmdb_cell[5],
                               vol, orthcode);
         clipper::Cell cell(clipper::Cell_descr(mmdb_cell[0], mmdb_cell[1], mmdb_cell[2],
                                                clipper::Util::d2rad(mmdb_cell[3]),
                                                clipper::Util::d2rad(mmdb_cell[4]),
                                                clipper::Util::d2rad(mmdb_cell[5])));
         setup(lines_mesh_for_cell, cell);
      }

      if (! xmap.is_null()) {
         setup(lines_mesh_for_cell, xmap.cell());
      }
   }

}
#endif

#ifndef EMSCRIPTEN
void
molecule_class_info_t::draw_unit_cell(Shader *shader_p,
                                      const glm::mat4 &mvp) {

   // 20220404-PE I can't use graphics_draw() like that - it puts coot into continuous-draw mode
   //             setup and draw need to be untangled. Another time.

   if (draw_it || draw_it_for_map) {
      if (show_unit_cell_flag) { // should be draw_it_for_unit_cell
         // if (lines_mesh_for_cell.empty())
         // setup_unit_cell(shader_p);
         glm::mat4 dummy(1.0f);
         lines_mesh_for_cell.draw(shader_p, mvp, dummy);
      }
   }

   // 20220320-PE
   // this is needed beacuse the unit cell is setup *during* (the first time) a draw call.
   // That's not right - it should be setup before this draw call.
   // setup and draw need to be untangled - but for now let's add an extra draw.
   //
   // graphics_info_t::graphics_draw();


}
#endif

// --------------------------------------------------------------------
//   Conversion functions
// --------------------------------------------------------------------
//
void
molecule_class_info_t::initialize_coordinate_things_on_read_molecule(std::string molecule_name) {

   // presume not an undo/redo by default:
   initialize_coordinate_things_on_read_molecule_internal(molecule_name, 0);
}

// If we are a redo/undo, then we don't want to update (add a) mol in
// display control widget
//
// Or a non graphics_info_t::molecules[] usage of this class.
//
void
molecule_class_info_t::initialize_coordinate_things_on_read_molecule_internal(std::string molecule_name,
                                                                              short int is_undo_or_redo) {

   //
   name_ = molecule_name;

   //
   draw_it = 1; // by default, display it, we change change this later, if we want.

   //
   if (! is_undo_or_redo) {
      bonds_colour_map_rotation = (imol_no + 1) * graphics_info_t::rotate_colour_map_on_read_pdb;
      while (bonds_colour_map_rotation > 360.0)
         bonds_colour_map_rotation -= 360.0;
      bonds_rotate_colour_map_flag = graphics_info_t::rotate_colour_map_on_read_pdb_flag;
//       std::cout << "::::::: in initialization setting bonds_colour_map_rotation "
//                 << bonds_colour_map_rotation << " for imol no " << imol_no << std::endl;
   }

   graphics_info_t g;
   if (g.use_graphics_interface_flag) {

      if (! is_undo_or_redo) {
         // std::cout << "DEBUG:: not an undo/redo!\n";
         // std::cout << "----------------------- initialize_coordinate_things_on_read_molecule_internal() calls "
         //           << "new_coords_mol_in_display_control_widget() " << std::endl;
         new_coords_mol_in_display_control_widget(); // uses draw_it
      }
      graphics_info_t::refresh_validation_graph_model_list();
      graphics_info_t::refresh_ramachandran_plot_model_list();
   }
}

void
molecule_class_info_t::set_symm_bond_colour_mol(int icol) {

   switch (icol) {
      case GREEN_BOND:
         glColor3f (combine_colour(0.1,0),
                    combine_colour(0.8,1),
                    combine_colour(0.1,2));
         break;
      case BLUE_BOND:
         glColor3f (combine_colour(0.2,0),
                    combine_colour(0.2,1),
                    combine_colour(0.8,2));
         break;
      case RED_BOND:
         glColor3f (combine_colour(0.8,0),
                    combine_colour(0.1,1),
                    combine_colour(0.1,2));
         break;
      case YELLOW_BOND:
         glColor3f (combine_colour(0.7,0),
                    combine_colour(0.7,1),
                    combine_colour(0.0,2));
         break;

      default:
         glColor3f (combine_colour(0.7, 0),
                    combine_colour(0.8, 1),
                    combine_colour(0.8, 2));
   }
}

void
molecule_class_info_t::set_symm_bond_colour_mol_and_symop(int icol, int isymop) {

//    std::cout << "in set_symm_bond_colour_mol_and_symop " << imol_no << " " << icol << " "
//              << isymop << " symmetry_rotate_colour_map_flag: "
//              << symmetry_rotate_colour_map_flag << "\n";

   if (symmetry_rotate_colour_map_flag) {
      if (symmetry_colour_by_symop_flag) {
         set_symm_bond_colour_mol_rotate_colour_map(icol, isymop);
      } else {
         set_symm_bond_colour_mol_rotate_colour_map(icol, 0);
      }
   } else {
      set_symm_bond_colour_mol(icol);
   }

}

void
molecule_class_info_t::set_symm_bond_colour_mol_rotate_colour_map(int icol, int isymop) {

   float step = graphics_info_t::rotate_colour_map_on_read_pdb/360.0;
   float rotation_size = float(icol + isymop) * step;
   std::vector<float> orig_colours(3);
   std::vector<float> rgb_new(3);
   std::vector<float> t_colours(3);

   switch (icol) {
   case GREEN_BOND:
      t_colours[0] = combine_colour(0.1, 0);
      t_colours[1] = combine_colour(0.8, 1);
      t_colours[2] = combine_colour(0.1, 2);
      rgb_new = rotate_rgb(t_colours, rotation_size);
      glColor3f (rgb_new[0], rgb_new[1], rgb_new[2]);
      break;
   case BLUE_BOND:
      t_colours[0] = combine_colour(0.2, 0);
      t_colours[1] = combine_colour(0.2, 1);
      t_colours[2] = combine_colour(0.8, 2);
      rgb_new = rotate_rgb(t_colours, rotation_size);
      glColor3f (rgb_new[0], rgb_new[1], rgb_new[2]);
      break;
   case RED_BOND:
      t_colours[0] = combine_colour(0.8, 0);
      t_colours[1] = combine_colour(0.1, 1);
      t_colours[2] = combine_colour(0.1, 2);
      rgb_new = rotate_rgb(t_colours, rotation_size);
      glColor3f (rgb_new[0], rgb_new[1], rgb_new[2]);
      break;
   case YELLOW_BOND:
      t_colours[0] = combine_colour(0.7, 0);
      t_colours[1] = combine_colour(0.7, 1);
      t_colours[2] = combine_colour(0.0, 2);
      rgb_new = rotate_rgb(t_colours, rotation_size);
      glColor3f (rgb_new[0], rgb_new[1], rgb_new[2]);
      break;

   default:
      t_colours[0] = combine_colour(0.6, 0);
      t_colours[1] = combine_colour(0.7, 1);
      t_colours[2] = combine_colour(0.7, 2);
      rgb_new = rotate_rgb(t_colours, rotation_size);
      glColor3f (rgb_new[0], rgb_new[1], rgb_new[2]);
   }
}



float
molecule_class_info_t::combine_colour(float v, int col_part_index) {

   // col_part_index is 0,1,2 for red gree blue components of the colour
   double w = graphics_info_t::symmetry_colour_merge_weight;
   return w*graphics_info_t::symmetry_colour[col_part_index] + v*(1.0-w);
}

// amount is not in degrees, it is in fractions of a circle, e.g. 10/360.
//
void
molecule_class_info_t::rotate_rgb_in_place(float *rgb, const float &amount) const {

   float hsv[3];
   convert_rgb_to_hsv_in_place(rgb, hsv);
   hsv[0] += amount;
   if (hsv[0] > 1.0) hsv[0] -= 1.0;
   convert_hsv_to_rgb_in_place(hsv, rgb);

}

// This allocated memory for xmap_is_diff_map, xmap_is_filled and
// contour_level, but *does not filll them!*
// So they need to be filled after calling this function.
void
molecule_class_info_t::initialize_map_things_on_read_molecule(std::string molecule_name,
                                                              bool is_diff_map,
                                                              bool is_anomalous_map,
                                                              bool swap_difference_map_colours) {

   if (false)
      std::cout << "------------------- initialize_map_things_on_read_molecule() "
                << " imol_no " << imol_no << " is_anomalous_map: " << is_anomalous_map
                << " is difference map " << is_diff_map << " swapcol: " << swap_difference_map_colours
                << std::endl;

   // unset coordinates, this is not a set of coordinates:
   atom_sel.n_selected_atoms = 0;
   atom_sel.mol = 0;  // tested (in set_undo_molecule()) to see if this
                      // was a coordinates molecule.  So maps have to
                      // set this to NULL.

   // Map initialization:
   // n_draw_vectors = 0;
   // draw_vectors = NULL;

   // n_diff_map_draw_vectors = 0;
   // diff_map_draw_vectors = NULL;

   xmap_is_diff_map = is_diff_map;

   show_unit_cell_flag = 0;
   have_unit_cell      = 0; // hmmm - CHECKME.

   if (is_diff_map) {
      if (! swap_difference_map_colours) {
         if (! is_anomalous_map) {
            map_colour.red   = 0.2;
            map_colour.green = 0.6;
            map_colour.blue  = 0.2;
         } else {
            map_colour.red   = 0.6;
            map_colour.green = 0.65;
            map_colour.blue  = 0.4;
         }
      } else {
         map_colour.red   = 0.6;
         map_colour.green = 0.2;
         map_colour.blue  = 0.2;
      }
   } else {
      std::vector<float> orig_colours(3); // convert this to using GdkRGBA
      orig_colours[0] =  0.3;
      orig_colours[1] =  0.62;
      orig_colours[2] =  0.8;
      float rotation_size = float(imol_no) * graphics_info_t::rotate_colour_map_for_map/360.0;
      // std::cout << "rotating map colour by " << rotation_size * 360.0 << std::endl;
      std::vector<float> rgb_new = rotate_rgb(orig_colours, rotation_size);
      map_colour.red   = rgb_new[0];
      map_colour.green = rgb_new[1];
      map_colour.blue  = rgb_new[2];
   }

   // negative contour level
   //
   if (! swap_difference_map_colours) {
      if (! is_anomalous_map) {
         map_colour_negative_level.red   = 0.6;
         map_colour_negative_level.green = 0.2;
         map_colour_negative_level.blue  = 0.2;
      } else {
         map_colour_negative_level.red   = 0.55;
         map_colour_negative_level.green = 0.25;
         map_colour_negative_level.blue  = 0.45;
      }
   } else {
      map_colour_negative_level.red   = 0.2;
      map_colour_negative_level.green = 0.6;
      map_colour_negative_level.blue  = 0.2;
   }
   name_ = molecule_name;

   clipper::Coord_orth cen(xmap.cell().a() * 0.5,
                           xmap.cell().b() * 0.5,
                           xmap.cell().c() * 0.5);
   float cell_a = xmap.cell().a();
   radial_map_colouring_do_radial_colouring = false;
   radial_map_colour_centre = cen;
   radial_map_colour_radius_min = 0.0;
   radial_map_colour_radius_max = 0.3 * cell_a;
   radial_map_colour_invert_flag = false;
   radial_map_colour_saturation = 0.5;

   colour_map_using_other_map_flag = false;

   if (graphics_info_t::use_graphics_interface_flag) {
      draw_it_for_map = 1;
   } else {
      draw_it_for_map = 0;
   }
   draw_it_for_map_standard_lines = 1; // display the map initially, by default

   // We can't call this untill xmap_is_filled[0] has been assigned,
   // and here we only make room for it.
   //
   // update_map_in_display_control_widget();

}


void
molecule_class_info_t::update_mol_in_display_control_widget() const {

   graphics_info_t g;

   // we don't want to add a display control hbox if we are simply
   // doing an undo: This is now deal with by the calling function.
   //
//    std::cout << "update_mol_in_display_control_widget() now" << std::endl;
//    std::cout << "update_mol_in_display_control_widget() passed derefrerence imol_no_ptr: "
//              << *imol_no_ptr << std::endl;
   std::string dmn = name_for_display_manager();

   update_name_in_display_control_molecule_combo_box(imol_no, dmn.c_str()); // because it's in gtk-manual.h - Fix later.
                                                              // note that display_control_map_combo_box() uses a string
}

void
molecule_class_info_t::new_coords_mol_in_display_control_widget() const {

   graphics_info_t g;

   // we don't want to add a display control hbox if we are simply
   // doing an undo: This is now handled by the calling function.
   //
   bool show_add_reps_flag = 0;
   if (add_reps.size() > 0)
      show_add_reps_flag = 1;

   std::string dmn = name_for_display_manager();
   display_control_molecule_combo_box(dmn.c_str(), imol_no, show_add_reps_flag);
   if (add_reps.size() > 0) {
      GtkWidget *vbox = display_control_add_reps_container(g.display_control_window(), imol_no);
      for (unsigned int iar=0; iar<add_reps.size(); iar++) {
         std::string name = coot::util::int_to_string(iar);
         name += " ";
         name += add_reps[iar].info_string();
         display_control_add_reps(vbox, imol_no, iar, add_reps[iar].show_it,
                                  add_reps[iar].bonds_box_type, name);
      }
   }

}

std::string
molecule_class_info_t::name_for_display_manager() const {

   std::string s("");
   if (graphics_info_t::show_paths_in_display_manager_flag) {
      s = name_;
   } else {
      if (has_model()) {
#ifdef WINDOWS_MINGW
         std::string::size_type islash = coot::util::intelligent_debackslash(name_).find_last_of("/");
#else
         std::string::size_type islash = name_.find_last_of("/");
#endif // MINGW
         if (islash == std::string::npos) {
            s = name_;
         } else {
            s = name_.substr(islash+1, name_.length());
         }
      } else {
         // This is a map, so we want to strip of xxx/ from each of
         // the (space separated) strings.
         // e.g.:
         // thing/other.mtz -> other.mtz
         // but
         // Averaged -> Averaged

         std::vector<std::string> v = coot::util::split_string(name_, " ");
         for (unsigned int i=0; i<v.size(); i++) {
            if (i > 0)
               s += " ";
            std::pair<std::string, std::string> p = coot::util::split_string_on_last_slash(v[i]);
            if (p.second == "")
               s += v[i];
            else
               s += p.second;
         }
      }
   }
   return s;
}

std::string
molecule_class_info_t::dotted_chopped_name() const {

   std::string ss = coot::util::int_to_string(imol_no);
   ss += " " ;
   int ilen = name_.length();
   int left_size = ilen-graphics_info_t::go_to_atom_menu_label_n_chars_max;
   if (left_size <= 0) {
      // no chop
      left_size = 0;
   } else {
      // chop
      ss += "...";
   }
   ss += name_.substr(left_size, ilen);
   return ss;
}



void
molecule_class_info_t::add_to_labelled_atom_list(int atom_index) {

   // this is quite passive. Caller does the atom label drawing.

   // note initialization n_labelled_atoms is 0;
   //
   if (is_in_labelled_list(atom_index) == true) {
      unlabel_atom(atom_index);
   } else {
      labelled_atom_index_list.push_back(atom_index);
   }
}


// or as we would say in lisp: rember
void
molecule_class_info_t::unlabel_atom(int atom_index) {

   //
   // Remove atom_index from the list of atoms to be labelled.
   //
   std::vector<int>::iterator it;
   for (it = labelled_atom_index_list.begin(); it != labelled_atom_index_list.end(); it++) {
      if ( *it == atom_index) {
         labelled_atom_index_list.erase(it);
         break;
      }
   }
}

void
molecule_class_info_t::unlabel_last_atom() {
   // remove the last atom from the list (if
   // there *are* atoms in the list, else do
   // nothing).
   unsigned int las = labelled_atom_index_list.size();
   if (las > 0) {
      int atom_index = labelled_atom_index_list[las-1];
      // std::cout << "unlabelling atom index" << atom_index << std::endl;
      unlabel_atom(atom_index);
   }
}

// or as we would say in lisp: member?
bool
molecule_class_info_t::is_in_labelled_list(int i) {

   // is the i'th atom in the list of atoms to be labelled?

   for (unsigned int ii=0; ii<labelled_atom_index_list.size(); ii++) {
      if (labelled_atom_index_list[ii] == i) {
         return 1;
      }
   }
   return 0;
}

// ------------------------- residue exists? -------------------------

int
molecule_class_info_t::does_residue_exist_p(const std::string &chain_id,
                                            int resno,
                                            const std::string &inscode) const {
   int state = 0;
   if (atom_sel.n_selected_atoms > 0) {
      int n_models = atom_sel.mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) {

         mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
         mmdb::Chain *chain_p;
         // run over chains of the existing mol
         int nchains = model_p->GetNumberOfChains();
         if (nchains <= 0) {
            std::cout << "ERROR:: bad nchains in molecule " << nchains
                      << std::endl;
         } else {
            for (int ichain=0; ichain<nchains; ichain++) {
               chain_p = model_p->GetChain(ichain);
               if (chain_p == NULL) {
                  // This should not be necessary. It seem to be a
                  // result of mmdb corruption elsewhere - possibly
                  // DeleteChain in update_molecule_to().
                  std::cout << "NULL chain in ... " << std::endl;
               } else {
                  mmdb::PResidue residue_p;
                  if (chain_id == chain_p->GetChainID()) {
                     int nres = chain_p->GetNumberOfResidues();
                     for (int ires=0; ires<nres; ires++) {
                        residue_p = chain_p->GetResidue(ires);
                        if (resno == residue_p->seqNum) {
                           if (inscode == residue_p->GetInsCode()) {
                              state = 1;
                              break;
                           }
                        }
                     }
                  }
               }
            }
         }
         if (state)
            break;
      }
   }
   return state;
}



// ------------------------- symmmetry atom labels -------------------------

void
molecule_class_info_t::add_atom_to_labelled_symm_atom_list(int atom_index,
                                                           const symm_trans_t &symm_trans,
                                                           const Cell_Translation &pre_shift_cell_trans) {

   if ( is_in_labelled_symm_list(atom_index) == 1 ) {
      unlabel_symm_atom(atom_index);
   } else {
      labelled_symm_atom_index_list.push_back(atom_index);
      std::pair<symm_trans_t, Cell_Translation> p(symm_trans, pre_shift_cell_trans);
      labelled_symm_atom_symm_trans_.push_back(p);
   }
}

// no need for this now we are using vectors.
// int
// molecule_class_info_t::labelled_symm_atom(int i) {
//    return labelled_symm_atom_index_list[i];
// }



std::pair<symm_trans_t, Cell_Translation>
molecule_class_info_t::labelled_symm_atom_symm_trans(int i) {

   return labelled_symm_atom_symm_trans_[i];
}


// old syle pointer using function
// int
// molecule_class_info_t::max_labelled_symm_atom() {

//    return n_labelled_symm_atoms;
// }

void
molecule_class_info_t::unlabel_symm_atom(int atom_index) {

   std::vector<int>::iterator it;
   for (it = labelled_symm_atom_index_list.begin();
        it != labelled_symm_atom_index_list.end();
        it++) {
      if ( *it == atom_index) {
         labelled_symm_atom_index_list.erase(it);
         break;
      }
   }
}

// shall we pass the symm_trans too?  Ideally we should, I think.
//
bool
molecule_class_info_t::is_in_labelled_symm_list(int i) {

   // is the i'th atom in the list of atoms to be labelled?

   for (unsigned int ii=0; ii<labelled_symm_atom_index_list.size(); ii++) {
      if (labelled_symm_atom_index_list[ii] == i) {
         return 1;
      }
   }
   return 0;
}


int molecule_class_info_t::add_atom_label(const char *chain_id, int iresno, const char *atom_id) {

   // int i = atom_index(chain_id, iresno, atom_id);
   int i = atom_spec_to_atom_index(std::string(chain_id),
                                   iresno,
                                   std::string(atom_id));
   if (i >= 0) // thanks Gabriele Balducci
      add_to_labelled_atom_list(i);
   else
      std::cout << atom_id << "/" << iresno << "/" << chain_id
                << " is not found in this molecule: (" <<  imol_no << ") "
                << name_ << std::endl;

   return i;
}

int
molecule_class_info_t::add_atom_labels_for_residue(mmdb::Residue *residue_p) {

   int n_atoms = 0;
   if (residue_p) {
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for(int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            int i = get_atom_index(at);
            add_to_labelled_atom_list(i);
            n_atoms++;
         }
      }
   }
   return n_atoms;
}

void
molecule_class_info_t::add_labels_for_all_CAs() {

   int imod = 1;
   if (! atom_sel.mol) return;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int n_res = chain_p->GetNumberOfResidues();
         for (int ires=0; ires<n_res; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            if (residue_p) {
               int n_atoms = residue_p->GetNumberOfAtoms();
               for (int iat=0; iat<n_atoms; iat++) {
                  mmdb::Atom *at = residue_p->GetAtom(iat);
                  if (! at->isTer()) {
                     std::string atom_name(at->name);
                     if (atom_name == " CA ") { // PDBv3 FIXME
                        int i = get_atom_index(at);
                        add_to_labelled_atom_list(i);
                     }
                  }
               }
            }
         }
      }
   }
}

void
molecule_class_info_t::local_b_factor_display(bool state,
                                              const coot::Cartesian &screen_centre) {

   float close_dist = 8.0;
   float close_dist_sqrd = close_dist * close_dist;
   if (state) {
      int imod = 1;
      if (! atom_sel.mol) return;
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      if (model_p) {
         std::vector<coot::generic_text_object_t> text_objects;
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  int n_atoms = residue_p->GetNumberOfAtoms();
                  for (int iat=0; iat<n_atoms; iat++) {
                     mmdb::Atom *at = residue_p->GetAtom(iat);
                     if (! at->isTer()) {
                        std::string ele = at->element;
                        if (ele != " H") { // don't show B factors of H atoms - crowded and not useful
                           float dx = at->x - screen_centre.x();
                           float dy = at->y - screen_centre.y();
                           float dz = at->z - screen_centre.z();
                           float dd = dx * dx + dy * dy + dz * dz;
                           if (dd < close_dist_sqrd) {
                              int handle = -1; // not set
                              std::string label = coot::util::float_to_string_using_dec_pl(at->tempFactor, 1);
                              coot::generic_text_object_t gto(label, handle, at->x + 0.2, at->y, at->z);
                              text_objects.push_back(gto);
                           }
                        }
                     }
                  }
               }
            }
         }
         if (! text_objects.empty())
            graphics_info_t::generic_texts = text_objects;
      }
   } else {
      graphics_info_t::generic_texts.clear();
   }
}



int molecule_class_info_t::remove_atom_label(const char *chain_id, int iresno, const char *atom_id) {

   int i = atom_index(chain_id, iresno, atom_id);
   if (i > 0)
      unlabel_atom(i);
   return i;
}


void
molecule_class_info_t::draw_molecule(short int do_zero_occ_spots,
                                     bool against_a_dark_background,
                                     bool show_cis_peptide_markups) {

   // show_cis_peptide_markups gets turned off by caller when there are intermediate atoms
   // displayed.
   if (has_model()) {
      if (draw_it == 1) {
         if (true) {
#ifdef USE_MOLECULES_TO_TRIANGLES
            if (! molrepinsts.size()) {
#endif
          deuterium_spots();
          if (do_zero_occ_spots)
             zero_occupancy_spots();
          display_bonds(against_a_dark_background);
          draw_fixed_atom_positions();
          if (show_ghosts_flag) {
             if (ncs_ghosts.size() > 0) {
                for (unsigned int ighost=0; ighost<ncs_ghosts.size(); ighost++) {
                   draw_ghost_bonds(ighost);
                }
             }
          }
          if (show_cis_peptide_markups)
             draw_cis_peptide_markups();
          draw_bad_CA_CA_dist_spots();
#ifdef USE_MOLECULES_TO_TRIANGLES
            }
#endif
         }
      }
   }
}


void
molecule_class_info_t::zero_occupancy_spots() const {

   if (bonds_box.n_zero_occ_spots > 0) {

      const std::pair<bool, float> &use_radius_limit = graphics_info_t::model_display_radius;

      glColor3f(0.8, 0.7, 0.7);
      float zsc = graphics_info_t::zoom;
      //glPointSize(145.0/zsc);
      // scale the pointer with the bond width
      glPointSize(30.0/std::min(zsc,(float)20)*std::max(bond_width, (float)4));
      glBegin(GL_POINTS);
      for (int i=0; i<bonds_box.n_zero_occ_spots; i++) {

         if ((use_radius_limit.first == false) ||
             (graphics_info_t::is_within_display_radius(bonds_box.zero_occ_spots_ptr[i]))) {
            glVertex3f(bonds_box.zero_occ_spots_ptr[i].x(),
                       bonds_box.zero_occ_spots_ptr[i].y(),
                       bonds_box.zero_occ_spots_ptr[i].z());
         }
      }
      glEnd();
   }
}

void
molecule_class_info_t::deuterium_spots() const {

   if (bonds_box.n_deuterium_spots > 0) {

      glColor3f(1.0, 0.2, 0.4);
      float zsc = graphics_info_t::zoom;
      glPointSize(165.0/zsc);
      glBegin(GL_POINTS);
      for (int i=0; i<bonds_box.n_deuterium_spots; i++) {
         glVertex3f(bonds_box.deuterium_spots_ptr[i].x(),
                    bonds_box.deuterium_spots_ptr[i].y(),
                    bonds_box.deuterium_spots_ptr[i].z());
      }
      glEnd();
   }
}

void
molecule_class_info_t::draw_cis_peptide_markups() const {

   const std::pair<bool, float> &use_radius_limit = graphics_info_t::model_display_radius;
   if (bonds_box.n_cis_peptide_markups > 0) {
      for (int i=0; i<bonds_box.n_cis_peptide_markups; i++) {
         const graphical_bonds_cis_peptide_markup &m = bonds_box.cis_peptide_markups[i];

         if ((single_model_view_current_model_number == 0) ||
             (single_model_view_current_model_number == m.model_number)) {

            if (! m.is_pre_pro_cis_peptide) {
               if (m.is_twisted) {
                  glColor3f(0.7, 0.6, 0.1);
               } else {
                  glColor3f(0.7, 0.2, 0.2);
               }
            } else {
               glColor3f(0.2, 0.7, 0.2);
            }

            coot::Cartesian fan_centre = m.pt_ca_1.mid_point(m.pt_ca_2);

            if ((! use_radius_limit.first)
                 || graphics_info_t::is_within_display_radius(fan_centre)) {

               coot::Cartesian v1 = fan_centre - m.pt_ca_1;
               coot::Cartesian v2 = fan_centre - m.pt_c_1;
               coot::Cartesian v3 = fan_centre - m.pt_n_2;
               coot::Cartesian v4 = fan_centre - m.pt_ca_2;

               coot::Cartesian pt_ca_1 = m.pt_ca_1 + v1 * 0.15;
               coot::Cartesian pt_c_1  = m.pt_c_1  + v2 * 0.15;
               coot::Cartesian pt_n_2  = m.pt_n_2  + v3 * 0.15;
               coot::Cartesian pt_ca_2 = m.pt_ca_2 + v4 * 0.15;

               glBegin(GL_TRIANGLE_FAN);

               glVertex3f(fan_centre.x(), fan_centre.y(), fan_centre.z());
               glVertex3f(pt_ca_1.x(), pt_ca_1.y(), pt_ca_1.z());
               glVertex3f(pt_c_1.x(),  pt_c_1.y(),  pt_c_1.z());
               glVertex3f(pt_n_2.x(),  pt_n_2.y(),  pt_n_2.z());
               glVertex3f(pt_ca_2.x(), pt_ca_2.y(), pt_ca_2.z());

               glEnd();
            }
         }
      }
   }
}

void
molecule_class_info_t::draw_bad_CA_CA_dist_spots() const {

   if (bonds_box.n_bad_CA_CA_dist_spots > 0) {

      glColor3f(0.9, 0.6, 0.3);
      const float &z = graphics_info_t::zoom;
      glPointSize(200.0 / z);
      glBegin(GL_POINTS);
      for (int i=0; i<bonds_box.n_bad_CA_CA_dist_spots; i++) {
        glVertex3f(bonds_box.bad_CA_CA_dist_spots_ptr[i].x(),
                   bonds_box.bad_CA_CA_dist_spots_ptr[i].y(),
                   bonds_box.bad_CA_CA_dist_spots_ptr[i].z());
      }
      glEnd();
   }
}

void
molecule_class_info_t::draw_fixed_atom_positions() const {

   if (fixed_atom_positions.size() > 0) {

      glColor3f(0.6, 0.95, 0.6);
      glPointSize(10.5);
      glBegin(GL_POINTS);
      for (unsigned int i=0; i<fixed_atom_positions.size(); i++) {
         glVertex3f(fixed_atom_positions[i].x(),
                    fixed_atom_positions[i].y(),
                    fixed_atom_positions[i].z());
      }
      glEnd();
   }
}

void
molecule_class_info_t::draw_ghost_bonds(int ighost) {

   stereo_eye_t eye = stereo_eye_t::MONO; // PASS THIS

#if 0 // olden code
   // hack in a value
   bool against_a_dark_background = true;

   if (ighost<int(ncs_ghosts.size())) {
      if (ncs_ghosts[ighost].display_it_flag) {
         glLineWidth(ghost_bond_width);
         for (int i=0; i<ncs_ghosts[ighost].bonds_box.num_colours; i++) {
            mmdb::Atom *at = atom_sel.atom_selection[i];
            std::string ele(at->element);
            int c = get_atom_colour_from_element(ele);
            if (ncs_ghosts[ighost].bonds_box.bonds_[i].num_lines > 0)
               set_bond_colour_by_mol_no(ighost, against_a_dark_background);
            glBegin(GL_LINES);
            for (int j=0; j< ncs_ghosts[ighost].bonds_box.bonds_[i].num_lines; j++) {
               glVertex3f(ncs_ghosts[ighost].bonds_box.bonds_[i].pair_list[j].positions.getStart().get_x(),
                          ncs_ghosts[ighost].bonds_box.bonds_[i].pair_list[j].positions.getStart().get_y(),
                          ncs_ghosts[ighost].bonds_box.bonds_[i].pair_list[j].positions.getStart().get_z());
               glVertex3f(ncs_ghosts[ighost].bonds_box.bonds_[i].pair_list[j].positions.getFinish().get_x(),
                          ncs_ghosts[ighost].bonds_box.bonds_[i].pair_list[j].positions.getFinish().get_y(),
                          ncs_ghosts[ighost].bonds_box.bonds_[i].pair_list[j].positions.getFinish().get_z());
            }
            glEnd();
         }
      }
   }
#endif

   if (ighost<int(ncs_ghosts.size())) {
      if (ncs_ghosts[ighost].display_it_flag) {
         Shader *shader_p = &graphics_info_t::shader_for_meshes_with_shadows;
         glm::mat4 mvp = graphics_info_t::get_molecule_mvp(eye);
         glm::mat4 model_rotation_matrix = graphics_info_t::get_model_rotation();
         glm::vec4 background_colour = graphics_info_t::get_background_colour();
         const auto &lights = graphics_info_t::lights;
         const auto &eye_position = graphics_info_t::eye_position;
         // mabye draw_with_shadows() should be used?

         // 20230826-PE ghost-molecule-display.hh was moved to the api directory some time ago.
         // That contains the draw method, but is if defed out (of course).
         // To bring back this draw() method, we need to derive a class, say draw_ncs_ghost_t,
         // that adds the draw method. For now, I will comment it out, but I should come back to it.
         //
         // ncs_ghosts[ighost].draw(shader_p, mvp, model_rotation_matrix, lights, eye_position, background_colour);
      }
   }

}



void
molecule_class_info_t::display_bonds(bool against_a_dark_background) {

}


// no symmetry here - for that see display_bonds(bool)
//
void
molecule_class_info_t::display_bonds(const graphical_bonds_container &bonds_box,
                                     float p_bond_width,
                                     bool against_a_dark_background) {

   // goodbye my old friend
}


void
molecule_class_info_t::display_bonds_stick_mode_atoms(const graphical_bonds_container &bonds_box,
                                                      const coot::Cartesian &front,
                                                      const coot::Cartesian &back,
                                                      bool against_a_dark_background) {

   // goodbye
}

coot::Cartesian
molecule_class_info_t::get_vector_pependicular_to_screen_z(const coot::Cartesian &front,
                                                           const coot::Cartesian &back,
                                                           const coot::Cartesian &bond_dir,
                                                           float zoom,
                                                           float p_bond_width)  const {

   coot::Cartesian bmf = back - front;
   // coot::Cartesian arb(0, 0.1, 0.9);
   coot::Cartesian p1 = coot::Cartesian::CrossProduct(bmf, bond_dir);

//    std::cout << "   crossproduct: " << p1 << " from "
//     << bmf << " and " << arb << "\n";

   p1.unit_vector_yourself();
   p1 *= zoom * 0.0004 * p_bond_width;

   return p1;
}


void
molecule_class_info_t::display_symmetry_bonds() {

   std::cout << "old code FIXME in display_symmetry_bonds() " << std::endl;

#if 0
   // We may come here after having done additional_representations -
   // which would change the line width.
   //
   glLineWidth(bond_width);

   auto rtop_to_opengl_matrix = [] (const clipper::RTop_orth &rtop, float *mat_for_glm) {

      // not used. Maybe useful one day (but not sure that it works as intended). c.f. gl_matrix class also.
      for (unsigned int i=0; i<3; i++)
         for (unsigned int j=0; j<3; j++)
            mat_for_glm[i*4 + j] = rtop.rot()(i,j);
      for (unsigned int i=0; i<3; i++)
         mat_for_glm[i*4 + 3 ] = 0.0;
      for (unsigned int i=0; i<3; i++)
         mat_for_glm[i + 12] = rtop.trn()[i];
      mat_for_glm[15] = 1.0f;
   };

   auto display_molecular_symmetry = [] (const graphical_bonds_container &bonds_box,
                                         const std::vector<std::pair<clipper::Mat33<double>, clipper::Coord_orth> > &molecular_symmetry_matrices) {

                                        // each symmetry-related molecule is drawn in it's own colour (no colour by atom (although it could be
                                        // easily changed to be so)).

                                        for (unsigned int i_ms=0; i_ms<molecular_symmetry_matrices.size(); i_ms++) {
                                           coot::colour_t rgb(0.66, 0.66, 0.3);
                                           float rotation_size = static_cast<float>(i_ms) * 0.28450f;
                                           while (rotation_size > 1.0f) rotation_size -= 1.0f;
                                           rgb.rotate(rotation_size);
                                           glColor3f(rgb[0], rgb[1], rgb[2]);
                                           const clipper::Mat33<double> &mat(molecular_symmetry_matrices[i_ms].first);
                                           const clipper::Coord_orth  &trans(molecular_symmetry_matrices[i_ms].second);
                                           const clipper::Coord_orth  minus_trans = -1.0 * trans;
                                           // std::cout << i_ms  << " mat\n" << mat.format() << "\n trans " << trans.format() << std::endl;
                                           clipper::Coord_orth zero_trans(0.0, 0.0, 0.0);
                                           clipper::RTop_orth rtop_symm_rot(mat, zero_trans);
                                           clipper::Mat33<double> identity_m(clipper::Mat33<double>::identity());
                                           clipper::RTop_orth origin_shift(identity_m, trans);
                                           clipper::RTop_orth minus_origin_shift(identity_m, minus_trans);
                                           clipper::RTop_orth back_origin_shift(identity_m, trans);
                                           clipper::RTop_orth rtop_A(minus_origin_shift);
                                           clipper::RTop_orth rtop_B(rtop_symm_rot * rtop_A);
                                           clipper::RTop_orth rtop_C(rtop_B * back_origin_shift);
                                           clipper::RTop_orth rtop(rtop_B);

                                           glBegin(GL_LINES);
                                           for (int icol=0; icol<bonds_box.num_colours; icol++) {
                                              graphical_bonds_lines_list<graphics_line_t> &ll = bonds_box.bonds_[icol];
                                              for (int j=0; j<ll.num_lines; j++) {
                                                 const coot::CartesianPair &pospair = ll.pair_list[j].positions;
                                                 const coot::Cartesian &p_1 = pospair.getStart();
                                                 const coot::Cartesian &p_2 = pospair.getFinish();

                                                 clipper::Coord_orth pt_A_1(p_1.x(), p_1.y(), p_1.z());
                                                 clipper::Coord_orth pt_A_2(p_2.x(), p_2.y(), p_2.z());

                                                 clipper::Coord_orth pt_B_1 = pt_A_1.transform(rtop_A);
                                                 clipper::Coord_orth pt_B_2 = pt_A_2.transform(rtop_A);

                                                 clipper::Coord_orth pt_C_1 = pt_B_1.transform(rtop_symm_rot);
                                                 clipper::Coord_orth pt_C_2 = pt_B_2.transform(rtop_symm_rot);

                                                 clipper::Coord_orth pt_D_1 = pt_C_1 + trans;
                                                 clipper::Coord_orth pt_D_2 = pt_C_2 + trans;

                                                 if (false) {
                                                    std::cout << i_ms << std::endl;
                                                    std::cout << "pt_A_1 " << pt_A_1.format() << std::endl;
                                                    std::cout << "pt_B_1 " << pt_B_1.format() << std::endl;
                                                    std::cout << "rtop A " << std::endl;
                                                    std::cout << rtop_A.format() << std::endl;
                                                    std::cout << "rtop_rot " << std::endl;
                                                    std::cout << rtop_B.format() << std::endl;
                                                 }

                                                 glVertex3d(pt_D_1.x(), pt_D_1.y(), pt_D_1.z());
                                                 glVertex3d(pt_D_2.x(), pt_D_2.y(), pt_D_2.z());
                                              }
                                           }
                                           glEnd();
                                           glPopMatrix();
                                        }
                                     };

   if ((show_symmetry == 1) && (graphics_info_t::show_symmetry == 1)) {

      if (! molecular_symmetry_matrices.empty()) {
         // display_molecular_symmetry(bonds_box, molecular_symmetry_matrices);
      } else {
         // std::cout << "molecular_symmetry_matrices empty\n";
      }

      for (unsigned int isym=0; isym<symmetry_bonds_box.size(); isym++) {
         int isymop = symmetry_bonds_box[isym].second.first.isym();
         if (symmetry_bonds_box[isym].first.symmetry_has_been_created == 1) {

            for (int icol=0; icol<symmetry_bonds_box[isym].first.num_colours; icol++) {

               set_symm_bond_colour_mol_and_symop(icol, isymop);
               int linesdrawn = 0;

               graphical_bonds_lines_list<graphics_line_t> &ll = symmetry_bonds_box[isym].first.symmetry_bonds_[icol];

               glBegin(GL_LINES);
               for (int j=0; j< symmetry_bonds_box[isym].first.symmetry_bonds_[icol].num_lines; j++) {

                  glVertex3f(ll.pair_list[j].positions.getStart().get_x(),
                             ll.pair_list[j].positions.getStart().get_y(),
                             ll.pair_list[j].positions.getStart().get_z());
                  glVertex3f(ll.pair_list[j].positions.getFinish().get_x(),
                             ll.pair_list[j].positions.getFinish().get_y(),
                             ll.pair_list[j].positions.getFinish().get_z());
                  if ( (++linesdrawn & 60023) == 0) {
                     glEnd();
                     glBegin(GL_LINES);
                     linesdrawn = 0;
                  }
               }
               glEnd();
            }
         }
      }

      if (show_strict_ncs_flag == 1) {
         // isn -> i_strict_ncs
         for (unsigned int isn=0; isn<strict_ncs_bonds_box.size(); isn++) {

            const graphical_bonds_container &gbc = strict_ncs_bonds_box[isn].first;

            if (0)
               std::cout << "here 3: isn "
                         << isn << " created_flag: "
                         << gbc.symmetry_has_been_created << " "
                         << "\n" ;

            // std::cout << "display_symmetry_bonds() here 5A " << gbc.symmetry_has_been_created
            // << std::endl;

            if (gbc.symmetry_has_been_created == 1) {

               // std::cout << "display_symmetry_bonds() here 6 " << std::endl;

               if (false)
                  std::cout << "num_colours: " << gbc.num_colours
                            << std::endl;

               for (int icol=0; icol<gbc.num_colours; icol++) {

                  if (0)
                     std::cout << "here 4 - isn: " << isn << " "
                               << "icol: " << icol << " num lines: "
                               << gbc.symmetry_bonds_[icol].num_lines
                               << "\n" ;

                  set_symm_bond_colour_mol_and_symop(icol, isn);
                  int linesdrawn = 0;

                  graphical_bonds_lines_list<graphics_line_t> &ll = gbc.symmetry_bonds_[icol];

                  glBegin(GL_LINES);
                  for (int j=0; j< gbc.symmetry_bonds_[icol].num_lines; j++) {

                     // pair = ll.pair_list[j];

                     glVertex3f(ll.pair_list[j].positions.getStart().get_x(),
                                ll.pair_list[j].positions.getStart().get_y(),
                                ll.pair_list[j].positions.getStart().get_z());
                     glVertex3f(ll.pair_list[j].positions.getFinish().get_x(),
                                ll.pair_list[j].positions.getFinish().get_y(),
                                ll.pair_list[j].positions.getFinish().get_z());
                     if ( (++linesdrawn & 60023) == 0) {
                        glEnd();
                        glBegin(GL_LINES);
                        linesdrawn = 0;
                     }
                  }
                  glEnd();
               }
            }
         }
      }
   }
#endif
}

// publically accessible
std::pair<coot::dipole, int>
molecule_class_info_t::add_dipole(const std::vector<coot::residue_spec_t> &res_specs,
                                  const coot::protein_geometry &geom) {

   int id = -1;
   coot::dipole d;
   std::vector<std::pair<coot::dictionary_residue_restraints_t, mmdb::Residue *> > pairs;

   for (unsigned int ires=0; ires<res_specs.size(); ires++) {
      mmdb::Residue *residue_p = get_residue(res_specs[ires]);

      if (residue_p) {
         try {
            std::string res_type = residue_p->GetResName();
            std::pair<short int, coot::dictionary_residue_restraints_t> rp =
               geom.get_monomer_restraints(res_type, imol_no);
            if (rp.first) {
               std::pair<coot::dictionary_residue_restraints_t, mmdb::Residue *> p(rp.second, residue_p);
               pairs.push_back(p);
            } else {
               std::cout << "INFO:: no monomer restraints found for "
                         << coot::residue_spec_t(residue_p) << " type: " << res_type << std::endl;
            }
         }
         catch (const std::runtime_error &mess) {
            std::cout << mess.what() << std::endl;
         }
      } else {
         std::cout << "   add_dipole: trapped null residue" << std::endl;
      }
   }

   if (pairs.size() > 0) {
      try {
         coot::dipole dl(pairs);
         dipoles.push_back(dl);
         id = dipoles.size() -1;
         d = dl;
      }
      catch (const std::runtime_error &mess) {
            std::cout << mess.what() << std::endl;
      }
   }
   return std::pair<coot::dipole, int> (d,id);
}

void
molecule_class_info_t::delete_dipole(int dipole_number) {

   if (dipole_number < int(dipoles.size())) {
      std::vector<coot::dipole>::iterator it;
      int n=0;
      for (it=dipoles.begin(); it!=dipoles.end(); ++it) {
         if (n == dipole_number) {
            dipoles.erase(it);
            break;
         }
         n++;
      }
   }
}


void
molecule_class_info_t::draw_dipoles() const {

   std::cout << "old code in draw_dipoles() " << std::endl;

#if 0
   if (! draw_it)
      return;

   if (dipoles.size() > 0) {
      glPushMatrix();
      glLineWidth(2.0);
      std::vector<clipper::Coord_orth> arrow_points;
      arrow_points.push_back(clipper::Coord_orth(0,0,0));
      arrow_points.push_back(clipper::Coord_orth(0.13,0,0));
      arrow_points.push_back(clipper::Coord_orth(0.13,0,0.66));
      arrow_points.push_back(clipper::Coord_orth(0.3,0,0.66));
      arrow_points.push_back(clipper::Coord_orth(0,0,1));
      // shift so that the centre of the arrow is at the origin
      for (unsigned int i=0; i<arrow_points.size(); i++)
         arrow_points[i] -= clipper::Coord_orth(0,0,0.4); // not 0.5,
                                                          // esthetics

      // Guide-line object
      if (0) {
         glColor3f(0.8,0.6,0.4);
         for (unsigned int i=0; i<dipoles.size(); i++) {
            clipper::Coord_orth pt = dipoles[i].position();
            clipper::Coord_orth d = dipoles[i].get_dipole();
            double sc = 3.0;
            glBegin(GL_LINES);
            glVertex3d(pt.x(), pt.y(), pt.z());
            glVertex3d(pt.x() + d.x() * sc,
                       pt.y() + d.y() * sc,
                       pt.z() + d.z() * sc);
            glVertex3d(pt.x(), pt.y(), pt.z());
            glVertex3d(pt.x() - d.x() * sc,
                       pt.y() - d.y() * sc,
                       pt.z() - d.z() * sc);
            glEnd();
         }
      }

      glColor3f(0.9, 0.6, 0.8);
      for (unsigned int i=0; i<dipoles.size(); i++) {
         clipper::Coord_orth pt = dipoles[i].position();
         clipper::Coord_orth  d = dipoles[i].get_dipole();
         clipper::Coord_orth d_unit = dipoles[i].get_unit_dipole();

         // make an arbitrary vector not parallel to d_unit.
         //
         clipper::Coord_orth arb(0,0.1,0.9);
         if (d_unit.y() < d_unit.z())
            arb = clipper::Coord_orth(0.0, 0.9, 0.1);
         if (d_unit.x() < d_unit.y())
            arb = clipper::Coord_orth(0.9, 0.0, 0.1);

         clipper::Coord_orth p1(clipper::Coord_orth::cross(arb, d_unit).unit());
         clipper::Coord_orth p2(clipper::Coord_orth::cross( p1, d_unit).unit());
         clipper::Coord_orth p3 = d_unit;

         GL_matrix m(p1.x(), p1.y(), p1.z(),
                     p2.x(), p2.y(), p2.z(),
                     p3.x(), p3.y(), p3.z());

         glPushMatrix();
         glTranslated(pt.x(), pt.y(), pt.z());
         glMultMatrixf(m.get());
         glBegin(GL_LINES);

         // scale the dipole
         //
         double ds = sqrt(d.lengthsq());
         std::vector<clipper::Coord_orth> local_arrow_points = arrow_points;
         for (unsigned int i=0; i<arrow_points.size(); i++)
            local_arrow_points[i] = 0.181818181818181818181818181818181818 * ds * arrow_points[i];

         for (unsigned int i=0; i<local_arrow_points.size()-1; i++) {
            glVertex3d(local_arrow_points[i].x(),
                       local_arrow_points[i].y(),
                       local_arrow_points[i].z());
            glVertex3d(local_arrow_points[i+1].x(),
                       local_arrow_points[i+1].y(),
                       local_arrow_points[i+1].z());
            glVertex3d(-local_arrow_points[i].x(),
                        local_arrow_points[i].y(),
                        local_arrow_points[i].z());
            glVertex3d(-local_arrow_points[i+1].x(),
                        local_arrow_points[i+1].y(),
                        local_arrow_points[i+1].z());
         }
         glEnd();
         glPopMatrix();
      }
      glPopMatrix();
   }
#endif
}

std::string
coot::atom_selection_info_t::name () const {

   std::string s = "An Add Rep atom sel string";
   return s;
}


// return the atom selection and the number of atoms
int
coot::atom_selection_info_t::select_atoms(mmdb::Manager *mol) const {

   int SelHnd = -1;
   const char *alt_conf_str = "*";
   if (alt_conf_is_set)
      alt_conf_str = altconf.c_str();
   if (type == BY_ATTRIBUTES) {
      SelHnd = mol->NewSelection();
      mol->SelectAtoms(SelHnd, 0, chain_id.c_str(),
                       resno_start, // starting resno, an int
                       ins_code.c_str(), // any insertion code
                       resno_start, // ending resno
                       ins_code.c_str(), // ending insertion code
                       "*", // any residue name
                       "*", // atom name
                       "*", // elements
                       alt_conf_str  // alt loc.
                       );
   }
   if (type == BY_STRING) {
      SelHnd = mol->NewSelection();
      mol->Select(SelHnd, mmdb::STYPE_ATOM, atom_selection_str.c_str(), mmdb::SKEY_NEW);
   }
   return SelHnd;
}


std::string
coot::atom_selection_info_t::mmdb_string() const {

   std::string s = atom_selection_str;
   if (type == BY_ATTRIBUTES) {
      s = "//";
      s += chain_id;
      s += "/";
      s += coot::util::int_to_string(resno_start);
      if (resno_end != resno_start) {
         s += "-";
         s += coot::util::int_to_string(resno_end);
      } else {
         if (!ins_code.empty()) {
            s += ".";
            s += ins_code;
         }
      }
   }
   return s;
}


void
coot::additional_representations_t::fill_bonds_box() {

   if (representation_type != coot::BALL_AND_STICK) {
      atom_selection_container_t atom_sel;

      atom_sel.mol = mol;
      atom_sel.SelectionHandle = mol->NewSelection();

      if (atom_sel_info.type == coot::atom_selection_info_t::BY_ATTRIBUTES) {

         mol->SelectAtoms(atom_sel.SelectionHandle,
                          0, atom_sel_info.chain_id.c_str(),
                          atom_sel_info.resno_start, atom_sel_info.ins_code.c_str(),
                          atom_sel_info.resno_end,   atom_sel_info.ins_code.c_str(),
                          "*", "*", "*", "*");
      }
      if (atom_sel_info.type == coot::atom_selection_info_t::BY_STRING) {
         mol->Select(atom_sel.SelectionHandle, mmdb::STYPE_ATOM,
                     atom_sel_info.atom_selection_str.c_str(),
                     mmdb::SKEY_NEW);

      }
      mol->GetSelIndex(atom_sel.SelectionHandle,
                       atom_sel.atom_selection,
                       atom_sel.n_selected_atoms);

      if (bonds_box_type == coot::NORMAL_BONDS) {
         Bond_lines_container bonds(atom_sel, 1, draw_hydrogens_flag);
         bonds_box.clear_up();
         bonds_box = bonds.make_graphical_bonds();
      }
      mol->DeleteSelection(atom_sel.SelectionHandle);
   }
}

std::string
coot::additional_representations_t::info_string() const {

   std::string s("Fat Bonds: ");

   if (representation_type == coot::BALL_AND_STICK) {
     s = "Ball and Stick: ";
   }

   if (representation_type == coot::STICKS) {
     s = "Sticks: ";
   }

   if (atom_sel_info.type == coot::atom_selection_info_t::BY_STRING)
      s += atom_sel_info.atom_selection_str;
   if (atom_sel_info.type == coot::atom_selection_info_t::BY_ATTRIBUTES) {
      s += atom_sel_info.chain_id;
      s += " ";
      s += coot::util::int_to_string(atom_sel_info.resno_start);
      if (atom_sel_info.resno_end != atom_sel_info.resno_start) {
         s += "-";
         s += coot::util::int_to_string(atom_sel_info.resno_end);
      }
      s += atom_sel_info.ins_code;
   }
   return s;
}

#ifndef EMSCRIPTEN
int
molecule_class_info_t::add_additional_representation(int representation_type,
                                                     const int &bonds_box_type,
                                                     float bonds_width,
                                                     bool draw_hydrogens_flag,
                                                     const coot::atom_selection_info_t &info,
                                                     GtkWidget *display_control_window,
                                                     const gl_context_info_t &glci,
                                                     const coot::protein_geometry *geom) {

/*   representation_types:

   enum { SIMPLE_LINES, STICKS, BALL_AND_STICK, LIQUORICE, SURFACE };

  bonds_box_type:
  enum {  UNSET_TYPE = -1, NORMAL_BONDS=1, CA_BONDS=2, COLOUR_BY_CHAIN_BONDS=3,
          CA_BONDS_PLUS_LIGANDS=4, BONDS_NO_WATERS=5, BONDS_SEC_STRUCT_COLOUR=6,
          BONDS_NO_HYDROGENS=15,
          CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR=7,
          CA_BONDS_PLUS_LIGANDS_B_FACTOR_COLOUR=14,
          COLOUR_BY_MOLECULE_BONDS=8,
          COLOUR_BY_RAINBOW_BONDS=9, COLOUR_BY_B_FACTOR_BONDS=10,
          COLOUR_BY_OCCUPANCY_BONDS=11};

*/


   float sphere_size    = bonds_width;
   bool  do_spheres     = true;

   coot::additional_representations_t rep(atom_sel.mol,
                                          representation_type,
                                          bonds_box_type,
                                          bonds_width, sphere_size, do_spheres,
                                          draw_hydrogens_flag, info);

   add_reps.push_back(rep);
   int n_rep = add_reps.size() -1;
   std::string name = rep.info_string();
   GtkWidget *vbox = display_control_add_reps_container(display_control_window, imol_no);
   display_control_add_reps(vbox, imol_no, n_rep, rep.show_it, rep.bonds_box_type, name);
   if (representation_type == coot::BALL_AND_STICK) {
      int display_list_handle_index = make_ball_and_stick(info.mmdb_string(),
                                                          bonds_width, sphere_size, do_spheres,
                                                          glci, geom);
      int n_display_list_tags = display_list_tags.size();
      if ((display_list_handle_index >= 0) &&
          (display_list_handle_index < n_display_list_tags)) {
         add_reps[n_rep].add_display_list_handle(display_list_handle_index);
      }
   }

   return n_rep;
}
#endif


// representation_number should be an unsigned int.
int
molecule_class_info_t::adjust_additional_representation(int representation_number,
                                                        const int &bonds_box_type_in,
                                                        float bonds_width,
                                                        bool draw_hydrogens_flag,
                                                        const coot::atom_selection_info_t &info,
                                                        bool show_it_flag_in) {
   return -1;
}


void
molecule_class_info_t::clear_additional_representation(int representation_number) {

   int n_add_reps = add_reps.size();
   if (n_add_reps > representation_number) {
      if (representation_number >= 0) {
         add_reps[representation_number].clear();
      }
   }
}

void
molecule_class_info_t::set_show_additional_representation(int representation_number,
                                                          bool on_off_flag) {

   int n_add_reps = add_reps.size();
   if (n_add_reps > representation_number) {
      if (representation_number >= 0) {
         add_reps[representation_number].show_it = on_off_flag;
         if (add_reps[representation_number].representation_type == coot::BALL_AND_STICK ||
             add_reps[representation_number].representation_type == coot::STICKS) {
           int dl_index = add_reps[representation_number].display_list_handle;
           // std::cout << "Ball and stick add rep toggled to " << on_off_flag << std::endl;
           display_list_tags[dl_index].display_it = on_off_flag;
         }
      }
   }
}

void
molecule_class_info_t::set_show_all_additional_representations(bool on_off_flag) {
   int n_reps = add_reps.size();
   for (int i=0; i<n_reps; i++)
      set_show_additional_representation(i, on_off_flag);
}

void
molecule_class_info_t::all_additional_representations_off_except(int rep_no,
                                                                 bool ball_and_sticks_off_too_flag) {

   int n_reps = add_reps.size();
   for (int i=0; i<n_reps; i++)
      if (i != rep_no)
         if (ball_and_sticks_off_too_flag ||
             add_reps[i].representation_type != coot::BALL_AND_STICK)
         set_show_additional_representation(i, 0);
}




// Return a pair.
//
// If first string of length 0 on error to construct dataname(s).
std::pair<std::string, std::string>
molecule_class_info_t::make_import_datanames(const std::string &f_col_in,
                                             const std::string &phi_col_in,
                                             const std::string &weight_col_in,
                                             int use_weights) const {

   // If use_weights return 2 strings, else set something useful only for pair.first

   std::string f_col = f_col_in;
   std::string phi_col = phi_col_in;
   std::string weight_col = weight_col_in;

#ifdef WINDOWS_MINGW
   std::string::size_type islash_f   = coot::util::intelligent_debackslash(  f_col).find_last_of("/");
   std::string::size_type islash_phi = coot::util::intelligent_debackslash(phi_col).find_last_of("/");
#else
   std::string::size_type islash_f   =      f_col.find_last_of("/");
   std::string::size_type islash_phi =    phi_col.find_last_of("/");
#endif // MINGW

   short int label_error = 0;

   if (islash_f != std::string::npos) {
      // f_col is of form e.g. xxx/yyy/FWT
      if (f_col.length() > islash_f)
         f_col = f_col.substr(islash_f+1);
      else
         label_error = 1;
   }

   if (islash_phi != std::string::npos) {
      // phi_col is of form e.g. xxx/yyy/PHWT
      if (phi_col.length() > islash_phi)
         phi_col = phi_col.substr(islash_phi+1);
      else
         label_error = 1;
   }

   if (use_weights) {
#ifdef WINDOWS_MINGW
      std::string::size_type islash_fom = coot::util::intelligent_debackslash(weight_col).find_last_of("/");
#else
      std::string::size_type islash_fom = weight_col.find_last_of("/");
#endif
      if (islash_fom != std::string::npos) {
         // weight_col is of form e.g. xxx/yyy/WT
         if (weight_col.length() > islash_fom)
            weight_col = weight_col.substr(islash_fom+1);
         else
            label_error = 1;
      }
   }


   std::pair<std::string, std::string> p("", "");

   if (!label_error) {
      std::string no_xtal_dataset_prefix= "/*/*/";
      if (use_weights) {
         p.first  = no_xtal_dataset_prefix + "[" +   f_col + " " +      f_col + "]";
         p.second = no_xtal_dataset_prefix + "[" + phi_col + " " + weight_col + "]";
      } else {
         p.first  = no_xtal_dataset_prefix + "[" +   f_col + " " + phi_col + "]";
      }
   }
   return p;
}



void
molecule_class_info_t::filter_by_resolution(clipper::HKL_data< clipper::datatypes::F_phi<float> > *fphidata,
                                            const float &reso_low,
                                            const float &reso_high) const {

   float inv_low  = 1.0/(reso_low*reso_low);
   float inv_high = 1.0/(reso_high*reso_high);
   int n_data = 0;
   int n_reset = 0;


   for (clipper::HKL_info::HKL_reference_index hri = fphidata->first(); !hri.last(); hri.next()) {
//        std::cout << "high: " << inv_high << " low: " << inv_low
//                  << " data: " << hri.invresolsq() << std::endl;
      n_data++;

      if ( hri.invresolsq() > inv_low &&
           hri.invresolsq() < inv_high) {
      } else {
         (*fphidata)[hri].f() = 0.0;
         n_reset++;
      }
   }
   if (n_data > 0) {
      float f = static_cast<float>(n_reset)/static_cast<float>(n_data);
      std::cout << "INFO:: Chopped " << n_reset << " data out of " << n_data << " (" << f << "%)" << std::endl;
   } else {
      std::cout << "INFO:: Chopped " << n_reset << " data out of " << n_data << std::endl;
   }
}


void
molecule_class_info_t::label_symmetry_atom(int i) {
    //

    // same test as has_model():
    if (has_model()) {

       unsigned int i_unsigned(i);

       if (i_unsigned < labelled_symm_atom_index_list.size()) {

          int iatom_index = labelled_symm_atom_index_list[i];

          if (iatom_index < atom_sel.n_selected_atoms) {

             // look at translate_atom_with_pre_shift(), it translate
             // the negative of the passed translation, so print the
             // negative.
             //
             std::pair <symm_trans_t, Cell_Translation> st = labelled_symm_atom_symm_trans_[i];
             std::pair <symm_trans_t, Cell_Translation> st_inv(st.first, st.second.inv());
             std::string label = make_symm_atom_label_string(atom_sel.atom_selection[iatom_index], st_inv);

             GLfloat blueish[3] = { 0.7, 0.7, 1.0 };
             glColor3fv(blueish);
             coot::Cartesian symm_point = translate_atom_with_pre_shift(atom_sel, iatom_index, st);

             // 	    glRasterPos3f(symm_point.get_x(),
             // 			  symm_point.get_y()+0.02,
             // 			  symm_point.get_z()+0.02);

             graphics_info_t::printString(label,
                                          symm_point.get_x(),
                                          symm_point.get_y()+0.02,
                                          symm_point.get_z()+0.02);
          }
       }
    }
 }

std::pair<std::string, clipper::Coord_orth>
molecule_class_info_t::make_atom_label_string(unsigned int ith_labelled_atom,
                                              int brief_atom_labels_flag,
                                              short int seg_ids_in_atom_labels_flag) const {

   mmdb::Atom *at = atom_sel.atom_selection[labelled_atom_index_list[ith_labelled_atom]];
   std::string label = make_atom_label_string(at, brief_atom_labels_flag, seg_ids_in_atom_labels_flag);
   clipper::Coord_orth p = coot::co(at);
   p += clipper::Coord_orth(0.02, 0.02, 0.02);

   return std::pair<std::string, clipper::Coord_orth>(label, p);
}



// Put a label at the ith atom of mol_class_info::atom_selection.
//
void
molecule_class_info_t::draw_atom_label(int atom_index,
                                       int brief_atom_labels_flag,
                                       short int seg_ids_in_atom_labels_flag,
                                       const glm::vec4 &atom_label_colour,
                                       stereo_eye_t eye,
                                       const glm::mat4 &mvp,
                                       const glm::mat4 &view_rotation) {

   if (has_model()) {
      if (atom_index < atom_sel.n_selected_atoms) {
         mmdb::Atom *atom = atom_sel.atom_selection[atom_index];
         if (atom) {
            std::string label;
            glm::vec3 position;

            position = glm::vec3(atom->x, atom->y, atom->z);
            label = make_atom_label_string(atom, brief_atom_labels_flag, seg_ids_in_atom_labels_flag);

            graphics_info_t g;
            g.tmesh_for_labels.draw_atom_label(label, position, atom_label_colour,
                                               &g.shader_for_atom_labels, eye, mvp, view_rotation,
                                               glm::vec4(g.background_colour, 1.0),
                                               g.shader_do_depth_fog_flag,
                                               g.perspective_projection_flag);

         }
      } else {
         std::cout << "ERROR:: draw_atom_label() trying to label atom out of range: "
                   << atom_index << " " << atom_sel.n_selected_atoms
                   << " Removing label\n";
         unlabel_atom(atom_index);
      }
   }

}

// Put a label at the ith atom of mol_class_info::atom_selection.
//
void
molecule_class_info_t::draw_symm_atom_label(int atom_index,
                                            const std::pair <symm_trans_t, Cell_Translation> &st,
                                            const glm::vec4 &atom_label_colour,
                                            const glm::mat4 &mvp,
                                            const glm::mat4 &view_rotation) {

   stereo_eye_t eye = stereo_eye_t::MONO;

   if (has_model()) {
      if (atom_index < atom_sel.n_selected_atoms) {
         mmdb::Atom *atom = atom_sel.atom_selection[atom_index];
         if (atom) {

            std::pair <symm_trans_t, Cell_Translation> st_inv(st.first, st.second.inv());
            std::string label = make_symm_atom_label_string(atom_sel.atom_selection[atom_index], st_inv);
            coot::Cartesian symm_point = translate_atom_with_pre_shift(atom_sel, atom_index, st);
            glm::vec3 position = glm::vec3(symm_point.x(), symm_point.y(), symm_point.z());

            graphics_info_t g;
            g.tmesh_for_labels.draw_atom_label(label, position, atom_label_colour,
                                               &g.shader_for_atom_labels, eye, mvp, view_rotation,
                                               glm::vec4(g.background_colour, 1.0),
                                               g.shader_do_depth_fog_flag,
                                               g.perspective_projection_flag);

         }
      } else {
         std::cout << "ERROR:: draw_atom_label() trying to label atom out of range: "
                   << atom_index << " " << atom_sel.n_selected_atoms
                   << " Removing label\n";
         unlabel_atom(atom_index);
      }
   }

}



void
molecule_class_info_t::set_have_unit_cell_flag_maybe(bool warn_about_missing_symmetry_flag) {

    // mmdb::CMMDBCryst *cryst_p = atom_sel.mol->get_cell_p();

    mmdb::mat44 my_matt;

    int err = atom_sel.mol->GetTMatrix(my_matt, 0, 0, 0, 0);

    if (err != 0) {
       have_unit_cell = 0;
       if (warn_about_missing_symmetry_flag)
          std::cout << "WARNING:: No Symmetry for this model" << std::endl;
    } else {
       have_unit_cell = 1;
    }
 }

void
   molecule_class_info_t::update_bonds_colour_using_map_rotation(float f) {

   bonds_colour_map_rotation = f;
   make_glsl_bonds_type_checked(__FUNCTION__);
}


// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------
// add_residue_indices argument is no longer used - the atom indices are added
// ------------------------------------------------------------------------------
// ------------------------------------------------------------------------------

void
molecule_class_info_t::makebonds(float min_dist, float max_dist, const coot::protein_geometry *geom_p) {

   // debug_atom_selection_container(atom_sel);
   // std::cout << "---------------------------------- makebonds() A " << std::endl;

   Bond_lines_container bonds(atom_sel, min_dist, max_dist);
   bonds_box.clear_up();
   bonds_box = bonds.make_graphical_bonds();
   bonds_box_type = coot::NORMAL_BONDS;
   if (! draw_hydrogens_flag)
      bonds_box_type = coot::BONDS_NO_HYDROGENS;

   make_glsl_bonds_type_checked(__FUNCTION__);

}

void
molecule_class_info_t::makebonds(float max_dist, const coot::protein_geometry *geom_p) {

   // std::cout << "---------------------------------- makebonds() B " << std::endl;
   Bond_lines_container bonds(atom_sel, max_dist, graphics_info_t::draw_missing_loops_flag);

   bonds_box.clear_up();
   bonds_box = bonds.make_graphical_bonds();
   make_glsl_bonds_type_checked(__FUNCTION__);
}

// we remove the argument for add_residue_indices because they are no longer useful.
// bond descriptions now have atom indices.
//
void
molecule_class_info_t::makebonds(const coot::protein_geometry *geom_p,
                                 const std::set<int> &no_bonds_to_these_atoms) {

   if (false)
      std::cout << "---------------------------------- makebonds() --- start --- "
                << "with is_intermediate_atoms_molecule " << is_intermediate_atoms_molecule
                << " draw_it: " << draw_it << std::endl;

   // Don't try to use OpenGL if we don't have graphics
   if (! graphics_info_t::use_graphics_interface_flag) return;

   // come back to this

   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: in makebonds() C -- start A --\n";
   err = glGetError();
   if (err) std::cout << "GL ERROR:: in makebonds() C -- start B --\n";

   std::set<int>::const_iterator it;
   if (false) { // debug no_bonds_to_these_atoms
      for (it=no_bonds_to_these_atoms.begin(); it!=no_bonds_to_these_atoms.end(); ++it) {
         int idx = *it;
         mmdb::Atom *at = atom_sel.atom_selection[idx];
         std::cout << "   makebonds() C: No bond to " << idx << " " << coot::atom_spec_t(at) << std::endl;
      }
   }

   int do_disulphide_flag = 1;
   int model_number = 0; // flag for all models
   bool do_sticks_for_waters = false;

   if (single_model_view_current_model_number != 0)
      model_number = single_model_view_current_model_number;

   Bond_lines_container bonds(atom_sel, imol_no, no_bonds_to_these_atoms,
                              geom_p, do_disulphide_flag, draw_hydrogens_flag,
                              graphics_info_t::draw_missing_loops_flag,
                              model_number, "dummy", false, false, do_sticks_for_waters);
   bonds_box.clear_up();
   bonds_box = bonds.make_graphical_bonds();

   if (false) { // debug
      int n_atoms = bonds_box.n_atoms();
      int n_bonds = bonds_box.n_bonds();
      std::cout << "########### makebonds() C: bonds_box from make_graphical_bonds() contains "
                << n_bonds << " bonds" << " and " << n_atoms << " atoms " <<  std::endl;
   }
   bonds_box_type = coot::NORMAL_BONDS;
   if (! draw_hydrogens_flag)
      bonds_box_type = coot::BONDS_NO_HYDROGENS;

   if (false)
      std::cout << "   makebonds() C calls make_glsl_bonds_type_checked() imol "
                << imol_no << " " << name_
                << " intermediate-atoms: " << is_intermediate_atoms_molecule << std::endl;

   make_glsl_bonds_type_checked(__FUNCTION__);

   if (false)
      std::cout << "---------------------------------- makebonds() ---  end  --- "
                << "with is_intermediate_atoms_molecule " << is_intermediate_atoms_molecule
                << " draw_it: " << draw_it << std::endl;

}

void
molecule_class_info_t::make_ca_bonds(float min_dist, float max_dist) {

   std::set<int> no_bonds_to_these_atom_indices;
   Bond_lines_container bonds(graphics_info_t::Geom_p(), "dummy-CA-mode", no_bonds_to_these_atom_indices, false);
   bonds.do_Ca_bonds(atom_sel, min_dist, max_dist, graphics_info_t::draw_missing_loops_flag);
   bonds_box = bonds.make_graphical_bonds_no_thinning();
   bonds_box_type = coot::CA_BONDS;
   // std::cout << "DEBUG()::"  << __FUNCTION__ << "() ca: bonds_box_type is now "
   // << bonds_box_type << std::endl;

   make_glsl_bonds_type_checked(__FUNCTION__);

}

void

molecule_class_info_t::make_ca_bonds(float min_dist, float max_dist, const std::set<int> &no_bonds_to_these_atom_indices) {

   Bond_lines_container bonds(graphics_info_t::Geom_p(), no_bonds_to_these_atom_indices);
   bonds.do_Ca_bonds(atom_sel, min_dist, max_dist, graphics_info_t::draw_missing_loops_flag);
   bonds_box = bonds.make_graphical_bonds_no_thinning();
   bonds_box_type = coot::CA_BONDS;
   make_glsl_bonds_type_checked(__FUNCTION__);

}



void
molecule_class_info_t::make_ca_bonds() {
   make_ca_bonds(2.4, 4.7);
}

void
molecule_class_info_t::make_ca_plus_ligands_bonds(coot::protein_geometry *geom_p) {

   std::set<int> no_bonds_to_these_atom_indices;
   Bond_lines_container bonds(geom_p, "dummy-CA-mode", no_bonds_to_these_atom_indices, false);
   bonds.do_Ca_plus_ligands_bonds(atom_sel, imol_no, geom_p, 2.4, 4.7, draw_hydrogens_flag,
                                  graphics_info_t::draw_missing_loops_flag);

   // 20250124-PE is this a hostage to fortune? We don't want ligands with fat bonds to hydrogen
   // atoms - but does that mean that one of the Protein chains will be thin?
   // bonds_box = bonds.make_graphical_bonds_no_thinning();
   // withe the thinning flag on, bonds with colour HYDROGEN_GREY_BOND are drawn thin.
   bonds_box = bonds.make_graphical_bonds();

   bonds_box_type = coot::CA_BONDS_PLUS_LIGANDS;
   make_glsl_bonds_type_checked(__FUNCTION__);

   // std::cout << "ca: bonds_box_type is now " << bonds_box_type << std::endl;
}

void
molecule_class_info_t::make_ca_plus_ligands_and_sidechains_bonds(coot::protein_geometry *geom_p) {

   std::set<int> no_bonds_to_these_atom_indices;
   Bond_lines_container bonds(geom_p, "dummy-CA-mode", no_bonds_to_these_atom_indices, false);
   bonds.do_Ca_plus_ligands_and_sidechains_bonds(atom_sel, imol_no, geom_p, 2.4, 4.7,
                                                 0.01, 1.9, draw_hydrogens_flag,
                                                 graphics_info_t::draw_missing_loops_flag);
   bonds_box = bonds.make_graphical_bonds_no_thinning();
   bonds_box_type = coot::CA_BONDS_PLUS_LIGANDS_AND_SIDECHAINS;
   make_glsl_bonds_type_checked(__FUNCTION__);

   // std::cout << "ca: bonds_box_type is now " << bonds_box_type << std::endl;
}

void
molecule_class_info_t::make_colour_by_chain_bonds(bool force_rebond) {

   std::set<int> no_bonds_to_these_atom_indices;
   // c-only-flag, goodsell-flag
   make_colour_by_chain_bonds(no_bonds_to_these_atom_indices, true, false, force_rebond);
}

void
molecule_class_info_t::make_colour_by_chain_bonds(const std::set<int> &no_bonds_to_these_atoms,
                                                  bool change_c_only_flag,
                                                  bool goodsell_mode,
                                                  bool force_rebonding) {

   // this function is called as a result of chaning the bonding mode in the Display Manager.
   // or by even opening the Display Manager.
   // We don't want to rebond if we don't have to (i.e the mode requested is the current mode)
   // so check the previous value of bonds_box_type so that we can know if it can be skipped.

   Bond_lines_container bonds(graphics_info_t::Geom_p(), no_bonds_to_these_atoms, draw_hydrogens_flag);

   bool do_rama_markup = false; // should we be more clever?
   bonds.do_colour_by_chain_bonds(atom_sel, false, imol_no, draw_hydrogens_flag,
                                  graphics_info_t::draw_missing_loops_flag,
                                  change_c_only_flag, goodsell_mode, do_rama_markup);
   bonds_box = bonds.make_graphical_bonds_no_thinning(); // make_graphical_bonds() is pretty
                                                         // stupid when it comes to thining.

   bonds_box = bonds.make_graphical_bonds(); // make_graphical_bonds() is pretty
                                             // stupid when it comes to thining.

   // testing previous values of bonds_box_type
   if (bonds_box_type != coot::COLOUR_BY_CHAIN_BONDS)
      force_rebonding = true;

   // maybe I always want to rebond in goodsell mode?
   if (goodsell_mode)
      // if (bonds_box_type != coot::COLOUR_BY_CHAIN_GOODSELL)
         force_rebonding = true;

   bonds_box_type = coot::COLOUR_BY_CHAIN_BONDS;

   if (goodsell_mode)
      bonds_box_type = coot::COLOUR_BY_CHAIN_GOODSELL;

   if (force_rebonding)
      make_glsl_bonds_type_checked(__FUNCTION__);

   // I don't think that this should be here - it should be in caller function
   //
   // OK (I agree)
   // if (graphics_info_t::glareas[0])
   //    graphics_info_t::graphics_draw();
}

void
molecule_class_info_t::make_colour_by_ncs_related_chains(bool goodsell_mode) {

   std::set<int> no_bonds_to_these_atoms; // should this be passed?
   Bond_lines_container bonds(graphics_info_t::Geom_p(), no_bonds_to_these_atoms, draw_hydrogens_flag);

   int model_number = 1;
   std::vector<std::vector<mmdb::Chain *> > ncs_related_chains = coot::ncs_related_chains(atom_sel.mol, model_number);
   bool change_c_only_flag = false;
   int draw_mode = 1; // not used currently (for atoms-and-bond or just atoms only)
   bonds.do_colour_by_ncs_related_chain_bonds(atom_sel, imol_no, ncs_related_chains, draw_mode, change_c_only_flag, goodsell_mode);
   bonds_box = bonds.make_graphical_bonds();
   bonds_box_type = coot::COLOUR_BY_CHAIN_GOODSELL;
   model_representation_mode = Mesh::representation_mode_t::BALLS_NOT_BONDS;
   make_glsl_bonds_type_checked(__FUNCTION__);

}

void
molecule_class_info_t::make_colour_by_molecule_bonds(bool force_rebonding) {

   //bool force_rebonding?

   Bond_lines_container bonds;
   // the imol_no is passed because we search the dictionary to find if
   // the residue of the atoms has a dictionary (if not, they are drawn large).
   bonds.do_colour_by_molecule_bonds(atom_sel, imol_no, draw_hydrogens_flag);
   bonds_box = bonds.make_graphical_bonds();
   bonds_box_type = coot::COLOUR_BY_MOLECULE_BONDS;
   make_glsl_bonds_type_checked(__FUNCTION__);

   // Put this in the caller
   // if (graphics_info_t::glarea)
   //    graphics_info_t::graphics_draw();

}




// caller is an optional argument
void
molecule_class_info_t::make_bonds_type_checked(const char *caller) {

   bool debug = false;

   // Note caller can be 0 (e.g. with clang) - so be aware of that when debugging.

   if (bonds_box_type == coot::UNSET_TYPE) bonds_box_type = coot::NORMAL_BONDS;

   std::string caller_s("NULL");
   if (caller) caller_s = std::string(caller);

   if (debug)
      std::cout << "debug:: plain make_bonds_type_checked() --------start--------- called by "
                << caller_s << "() with is_intermediate_atoms_molecule: " << is_intermediate_atoms_molecule
                << std::endl;
   if (debug)
      std::cout << "--------- make_bonds_type_checked() called with bonds_box_type "
                << bonds_box_type << " vs "
                << "NORMAL_BONDS " << coot::NORMAL_BONDS << " "
                << "BONDS_NO_HYDROGENS " << coot::BONDS_NO_HYDROGENS << " "
                << "COLOUR_BY_CHAIN_BONDS " << coot::COLOUR_BY_CHAIN_BONDS << " "
                << "COLOUR_BY_MOLECULE_BONDS " << coot::COLOUR_BY_MOLECULE_BONDS << " "
                << "CA_BONDS " << coot::CA_BONDS << " "
                << "CA_BONDS_PLUS_LIGANDS " << coot::CA_BONDS_PLUS_LIGANDS << " "
                << "COLOUR_BY_USER_DEFINED_COLOURS___BONDS " << coot::COLOUR_BY_USER_DEFINED_COLOURS____BONDS << " "
                << std::endl;

   // Delete this in due course
   graphics_info_t g; // urgh!  (But the best solution?)

   bool force_rebonding = true; // if we get here, this must be true (?)

   if (! g.use_graphics_interface_flag) return;

   coot::protein_geometry *geom_p = g.Geom_p();

   std::set<int> dummy;

   // std::cout << "bonds_box_type " << bonds_box_type << std::endl;

   if (bonds_box_type == coot::NORMAL_BONDS) {
      if (debug)
         std::cout << "debug:: plain make_bonds_type_checked() calls makebonds() with geom_p " << geom_p << std::endl;
      makebonds(geom_p, dummy);
   }
   if (bonds_box_type == coot::BONDS_NO_HYDROGENS)
      makebonds(geom_p, dummy);
   if (bonds_box_type == coot::CA_BONDS)
      make_ca_bonds();
   if (bonds_box_type == coot::COLOUR_BY_CHAIN_BONDS || bonds_box_type == coot::COLOUR_BY_CHAIN_GOODSELL) {
      // Baah, we have to use the static in graphics_info_t here as it
      // is not a per-molecule property.
      std::set<int> s;
      bool goodsell_mode = false;
      if (bonds_box_type == coot::COLOUR_BY_CHAIN_GOODSELL)
         goodsell_mode = true;
      make_colour_by_chain_bonds(s, g.rotate_colour_map_on_read_pdb_c_only_flag, goodsell_mode, force_rebonding);
   }
   if (bonds_box_type == coot::COLOUR_BY_MOLECULE_BONDS)
      make_colour_by_molecule_bonds(force_rebonding);
   if (bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS)
      make_ca_plus_ligands_bonds(g.Geom_p());
   if (bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS_AND_SIDECHAINS)
      make_ca_plus_ligands_and_sidechains_bonds(g.Geom_p());
   if (bonds_box_type == coot::BONDS_NO_WATERS)
      bonds_no_waters_representation();
   if (bonds_box_type == coot::BONDS_SEC_STRUCT_COLOUR)
      bonds_sec_struct_representation();
   if (bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR)
      ca_plus_ligands_sec_struct_representation(g.Geom_p());
   if (bonds_box_type == coot::COLOUR_BY_RAINBOW_BONDS)
      ca_plus_ligands_rainbow_representation(g.Geom_p());
   if (bonds_box_type == coot::COLOUR_BY_OCCUPANCY_BONDS)
      occupancy_representation();
   if (bonds_box_type == coot::COLOUR_BY_B_FACTOR_BONDS)
      b_factor_representation();
   if (bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS_B_FACTOR_COLOUR)
      b_factor_representation_as_cas();
   if (bonds_box_type == coot::COLOUR_BY_USER_DEFINED_COLOURS____BONDS)
      user_defined_colours_representation(g.Geom_p(), true, g.draw_missing_loops_flag); // hack,
                                                             // because we need to remeber somehow
                                                             // if this was called with all-atom or CA-only.
                                                             // See c-interface.cc
                                                             // graphics_to_user_defined_atom_colours_representation()
                                                             // Perhaps we need two functions
                                                             // user_defined_colours_representation_all()
                                                             // user_defined_colours_representation_Calpha() [+ ligands]

   if (bonds_box_type == coot::COLOUR_BY_USER_DEFINED_COLOURS_CA_BONDS)
      user_defined_colours_representation(g.Geom_p(), false, g.draw_missing_loops_flag); // hack,

   if (debug) bonds_box.debug();

   // bleugh. But if we don't do this here, where *do* we do it?
   // Should the glci be passed to make_bonds_type_checked()?  Urgh.
   // That is called from many places....
   //

   gl_context_info_t glci = graphics_info_t::get_gl_context_info();

   // make glsl triangles
   glUseProgram(graphics_info_t::shader_for_models.get_program_id());
   // std::cout << "make_bonds_type_checked() using model shader program_id is "
   //           << graphics_info_t::shader_for_models.get_program_id() << std::endl;
   GLenum err = glGetError();
   if (err) std::cout << "Error in glUseProgram() in make_bonds_type_checked() " << err << "\n";

   GLint current_program;
   glGetIntegerv(GL_CURRENT_PROGRAM, &current_program);
   // std::cout << "INFO:: make_bonds_type_checked() current program " << current_program << std::endl;

   // all these will need to be changed or removed
   update_additional_representations(glci, g.Geom_p());
   update_fixed_atom_positions();
   update_ghosts();
   update_extra_restraints_representation();

   if (debug) {
      std::cout << "debug:: -------------- make_bonds_type_checked() done " << draw_it << std::endl;
   }
}

void
molecule_class_info_t::set_atom_radius_scale_factor(float sf) {

   atom_radius_scale_factor = sf;
   make_glsl_bonds_type_checked(__FUNCTION__);
}

std::vector<glm::vec4>
molecule_class_info_t::make_colour_table() const {

   if (false)
      std::cout << ":::::::::::: in make_colour_table() imol: " << imol_no
                << " bonds_box_type is " << bonds_box_type << " vs "
                << coot::COLOUR_BY_B_FACTOR_BONDS << std::endl;

   graphics_info_t g; // Hmm..

   bool debug_colour_table = false;

   float goodselliness = graphics_info_t::goodselliness;

   // 20220214-PE does this matter (is it useful?) now with modern graphics?
   // 20221125-PE it sure does (yes it is).
   bool dark_bg_flag = graphics_info_t::background_is_black_p();

   float gcwrs = graphics_info_t::goodsell_chain_colour_wheel_rotation_step;

   std::vector<glm::vec4> colour_table(bonds_box.num_colours, glm::vec4(0.5f, 0.5f, 0.5f, 1.0f));

   if (debug_colour_table) {
      std::cout << "make_colour_table(): A bonds_box.n_consolidated_atom_centres "
                << bonds_box.n_consolidated_atom_centres << std::endl;
      std::cout << "make_colour_table(): A bonds_box has " << bonds_box.num_colours << " colours" << std::endl;
      std::cout << "make_colour_table(): A colour_table has size " << colour_table.size() << std::endl;
   }

   for (int icol=0; icol<bonds_box.num_colours; icol++) {
      if (bonds_box_type == coot::COLOUR_BY_RAINBOW_BONDS) {
         glm::vec4 col = get_bond_colour_by_colour_wheel_position(icol, coot::COLOUR_BY_RAINBOW_BONDS);
         // std::cout << "rainbow " << icol << " " << glm::to_string(col) << std::endl;
         colour_table[icol] = col;
      } else {
         // this is the old way of doing user-defined colours. Now we use
         // set_user_defined_atom_colour_by_selection()
         if (bonds_box_type == coot::COLOUR_BY_USER_DEFINED_COLOURS_CA_BONDS) {
            if (! graphics_info_t::user_defined_colours.empty()) {
               int n_ud_colours = graphics_info_t::user_defined_colours.size();
               if (icol < n_ud_colours) {
                  unsigned int idx               = graphics_info_t::user_defined_colours[icol].first;
                  const coot::colour_holder &col = graphics_info_t::user_defined_colours[icol].second;
                  glm::vec4 ud_col(col.red, col.green, col.blue, 1.0);
                  if (idx < colour_table.size())
                     colour_table[idx] = ud_col;
                  else
                     std::cout << "ERROR:: in make_colour_table() trapped bad index " << std::endl;
               } else {
                  std::cout << "WARNING:: in make_colour_table() out of index colour COLOUR_BY_USER_DEFINED_COLOURS_CA_BONDS "
                            << icol << " " << graphics_info_t::user_defined_colours.size() << std::endl;
               }
            } else {
               std::cout << "WARNING:: in make_colour_table() user_defined_colours was empty " << std::endl;
            }
         } else {
            if (bonds_box_type == coot::COLOUR_BY_CHAIN_GOODSELL) {

               // goodsell colours start at 100. There are 2 colours per chain, so for A and be chains the
               // colour indices are 100, 101, 102, 103.
               if (debug_colour_table)
                   std::cout << "make_colour_table(): goodsell mode icol: " << icol << " n_bonds: "
                   << bonds_box.bonds_[icol].num_lines << std::endl;
               // if (bonds_box.bonds_[icol].num_lines > 0) { // we need this to work if there are only atoms
               if (true) {
                  // the first 100 colours don't count in Goodsell mode.
                  // coot::colour_holder ch(0.8, 0.5, 0.6);
                  // coot::colour_holder ch(0.8, 0.2, 0.1);
                  coot::colour_holder ch(0.8, 0.5, 0.6);
                  int ic = icol - 100;
                  bool is_C = !(ic %2);
                  int chain_index = ic/2;
                  if (debug_colour_table)
                     std::cout << "  icol " << icol << " ic " << ic << " is_C: " << is_C
                               << " chain_index " << chain_index << std::endl;
                  float rotation_amount = gcwrs * static_cast<float>(chain_index);
                  if (is_C)
                     ch.make_pale(goodselliness);
                  ch.rotate_by(rotation_amount);
                  colour_table[icol] = colour_holder_to_glm(ch);
                  if (debug_colour_table)
                     std::cout << "make_colour_table(): col " << icol << " "
                               << glm::to_string(colour_table[icol]) << std::endl;
               }
            } else {
               if (bonds_box_type == coot::COLOUR_BY_B_FACTOR_BONDS ||
                   bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS_B_FACTOR_COLOUR) {
                  glm::vec4 col = get_bond_colour_by_colour_wheel_position(icol, bonds_box_type);
                  colour_table[icol] = col;
               } else {

                  coot::colour_t cc = get_bond_colour_by_mol_no(icol, dark_bg_flag);
                  cc.brighter(0.8); // calm down - now that we are using the instanced-object.shader - the molecule is too bright.
                  colour_table[icol] = cc.to_glm();
                  // std::cout << "..................... this path icol " << icol << " cc " << cc << std::endl;
               }
            }
         }
      }

      // was there a graphics_info_t user-defined bond colour that superceeds the colour for this icol?

      if (! g.user_defined_colours.empty()) {
         int colour_table_size = g.user_defined_colours.size();
         if (debug_colour_table) {
            std::cout << "checking" << std::endl;
            g.print_user_defined_colour_table();
         }
         // replace this colour table with a colour from the user_defined_colours?
         for (unsigned int i=0; i<g.user_defined_colours.size(); i++) {
            int idx = g.user_defined_colours[i].first;
            if (idx == icol) {
               colour_table[icol] = colour_holder_to_glm(g.user_defined_colours[i].second);
            }
         }
      }
   }

   if (debug_colour_table)
      std::cout << "------------ make_colour_table(): B colour_table has size " << colour_table.size()
                << std::endl;

   // 20220303-PE why does this happen? (it happens when refining the newly imported 3GP ligand)
   // I guess the bonds_box for the remaining atoms (there are none of them) is incorrectly constructed.
   // FIXME later.
   // Note: we were called fromm this function:
   // void
   // molecule_class_info_t::makebonds(const coot::protein_geometry *geom_p,
   //                              const std::set<int> &no_bonds_to_these_atoms)
   //
   if (bonds_box.n_consolidated_atom_centres > bonds_box.num_colours) {
      std::cout << "WARNING:: make_colour_table() n_consolidated_atom_centres is " << bonds_box.n_consolidated_atom_centres
                << " therefore resizing the colour table " << std::endl;
      colour_table = std::vector<glm::vec4>(bonds_box.n_consolidated_atom_centres, glm::vec4(0.6f, 0.0f, 0.6f, 1.0f));
   }

   auto pastelize = [] (glm::vec4 &col, float degree) {
                       for (unsigned int i=0; i<3; i++) {
                          const float &cc = col[i];
                          float r = 1.0f - cc;
                          col[i] += r * degree;
                          col[i] *= (1.0f - 0.5f * degree); // I don't want bright pastel
                       }
                    };

   if (debug_colour_table)
      std::cout << "------------ make_colour_table(): C colour_table has size " << colour_table.size()
                << std::endl;

   if (is_intermediate_atoms_molecule) {
      // pastelize the colour table
      float degree = 0.5f;
      for (auto &col : colour_table) {
         pastelize(col, degree); //ref
      }
   }

   if (debug_colour_table) {
      std::cout << "------------ make_colour_table(): D colour_table has size " << colour_table.size()
                << std::endl;
      std::cout << "------------ make_colour_table(): E colour table for bonds_box_type " << bonds_box_type
                << " ---------" << std::endl;
      std::cout << "------------ make_colour_table(): E colour_table has size " << colour_table.size()
                << std::endl;
      for (unsigned int icol=0; icol<colour_table.size(); icol++) {
         graphical_bonds_lines_list<graphics_line_t> &ll = bonds_box.bonds_[icol];
         int n_bonds = ll.num_lines;
         float s = colour_table[icol][0] + colour_table[icol][1] + colour_table[icol][2];
         std::cout << "make_colour_table(): imol " << imol_no << " colour-table index "
                   << std::setw(2) << icol << " n-bonds: " << std::setw(4) << n_bonds << " "
                   << glm::to_string(colour_table[icol]) << " br: " << s << std::endl;
      }
   }

   if (false)
      std::cout << ":::::::::::: in make_colour_table() done." << std::endl;

   return colour_table;
}


void
molecule_class_info_t::make_mesh_from_bonds_box() {

   // it's all instanced now.
   std::cout << "don't use make_mesh_from_bonds_box() - it's all instanced now " << std::endl;
}

//! user-defined atom selection to colour index
void
molecule_class_info_t::set_user_defined_atom_colour_by_selection(const std::vector<std::pair<std::string, unsigned int> > &indexed_residues_cids,
                                                                 bool apply_colour_to_non_carbon_atoms_also) {

   // Fill user_defined_bond_colours

   if (! atom_sel.mol) return;

   int udd_handle = atom_sel.mol->GetUDDHandle(mmdb::UDR_ATOM, "user-defined-atom-colour-index");
   if (udd_handle == 0)
      udd_handle = atom_sel.mol->RegisterUDInteger(mmdb::UDR_ATOM, "user-defined-atom-colour-index");

   for (unsigned int i=0; i<indexed_residues_cids.size(); i++) {
      const auto &rc = indexed_residues_cids[i];
      const std::string &cid = rc.first;
      int colour_index = rc.second; // change type
      int selHnd = atom_sel.mol->NewSelection(); // d

      mmdb::Atom **SelAtoms = nullptr;
      int nSelAtoms = 0;
      atom_sel.mol->Select(selHnd, mmdb::STYPE_ATOM, cid.c_str(), mmdb::SKEY_NEW);
      atom_sel.mol->GetSelIndex(selHnd, SelAtoms, nSelAtoms);
      if (nSelAtoms > 0) {
         for (int iat=0; iat<nSelAtoms; iat++) {
            mmdb::Atom *at = SelAtoms[iat];
            std::string ele(at->element);
            if (apply_colour_to_non_carbon_atoms_also || ele == " C") {
               int ierr = at->PutUDData(udd_handle, colour_index);
               if (ierr != mmdb::UDDATA_Ok) {
                  std::cout << "WARNING:: in set_user_defined_atom_colour_by_residue() problem setting udd on atom "
                            << coot::atom_spec_t(at) << std::endl;
               }
            }
         }
      }
      atom_sel.mol->DeleteSelection(selHnd);

#if 0 // 20230919-PE old - select residues - this is wrong.
      mmdb::Residue **SelResidues;
      int nSelResidues = 0;
      atom_sel.mol->Select(selHnd, mmdb::STYPE_RESIDUE, cid.c_str(), mmdb::SKEY_NEW);
      atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
      if (nSelResidues > 0) {
         for(int ires=0; ires<nSelResidues; ires++) {
            mmdb::Residue *residue_p = SelResidues[ires];

	    mmdb::Atom **residue_atoms = 0;
	    int n_residue_atoms;
	    residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
	    for (int iat=0; iat<n_residue_atoms; iat++) {
	       mmdb::Atom *at = residue_atoms[iat];
	       int ierr = at->PutUDData(udd_handle, colour_index);
               std::string ele(at->element);
               if (apply_colour_to_non_carbon_atoms_also || ele == " C") {
                  if (ierr != mmdb::UDDATA_Ok) {
                     std::cout << "WARNING:: in set_user_defined_atom_colour_by_residue() problem setting udd on atom "
                               << coot::atom_spec_t(at) << std::endl;
                  } else {
                     if (false)
                        std::cout << "set user-defined colour index for " << at << " to " << colour_index << std::endl;
                  }
               }
            }
         }
      }
      atom_sel.mol->DeleteSelection(selHnd);
#endif


   }
}


void
molecule_class_info_t::make_meshes_from_bonds_box_instanced_version() {

   // this function presumes that bonds_box has been set before this function is called.

   // what is the api_bond_colour_t for the given bbt?
   auto convert_box_box_type = [] (int bbt) {
      coot::api_bond_colour_t abbt(coot::api_bond_colour_t::NORMAL_BONDS);

      if (bbt == coot::UNSET_TYPE) abbt = coot::api_bond_colour_t::UNSET_TYPE;
      if (bbt == coot::NORMAL_BONDS) abbt = coot::api_bond_colour_t::NORMAL_BONDS;
      if (bbt == coot::CA_BONDS) abbt = coot::api_bond_colour_t::CA_BONDS;
      if (bbt == coot::COLOUR_BY_CHAIN_BONDS) abbt = coot::api_bond_colour_t::COLOUR_BY_CHAIN_BONDS;
      if (bbt == coot::CA_BONDS_PLUS_LIGANDS) abbt = coot::api_bond_colour_t::CA_BONDS_PLUS_LIGANDS;
      if (bbt == coot::BONDS_NO_WATERS) abbt = coot::api_bond_colour_t::BONDS_NO_WATERS;
      if (bbt == coot::BONDS_SEC_STRUCT_COLOUR) abbt = coot::api_bond_colour_t::BONDS_SEC_STRUCT_COLOUR;
      if (bbt == coot::BONDS_NO_HYDROGENS) abbt = coot::api_bond_colour_t::BONDS_NO_HYDROGENS;
      if (bbt == coot::CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR) abbt = coot::api_bond_colour_t::CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR;
      if (bbt == coot::CA_BONDS_PLUS_LIGANDS_B_FACTOR_COLOUR) abbt = coot::api_bond_colour_t::CA_BONDS_PLUS_LIGANDS_B_FACTOR_COLOUR;
      if (bbt == coot::CA_BONDS_PLUS_LIGANDS_AND_SIDECHAINS) abbt = coot::api_bond_colour_t::CA_BONDS_PLUS_LIGANDS_AND_SIDECHAINS;
      if (bbt == coot::COLOUR_BY_MOLECULE_BONDS) abbt = coot::api_bond_colour_t::COLOUR_BY_MOLECULE_BONDS;
      if (bbt == coot::COLOUR_BY_RAINBOW_BONDS) abbt = coot::api_bond_colour_t::COLOUR_BY_RAINBOW_BONDS;
      if (bbt == coot::COLOUR_BY_B_FACTOR_BONDS) abbt = coot::api_bond_colour_t::COLOUR_BY_B_FACTOR_BONDS;
      if (bbt == coot::COLOUR_BY_OCCUPANCY_BONDS) abbt = coot::api_bond_colour_t::COLOUR_BY_OCCUPANCY_BONDS;
      if (bbt == coot::COLOUR_BY_USER_DEFINED_COLOURS____BONDS) abbt = coot::api_bond_colour_t::COLOUR_BY_USER_DEFINED_COLOURS____BONDS;
      if (bbt == coot::COLOUR_BY_USER_DEFINED_COLOURS_CA_BONDS) abbt = coot::api_bond_colour_t::COLOUR_BY_USER_DEFINED_COLOURS_CA_BONDS;
      if (bbt == coot::COLOUR_BY_CHAIN_GOODSELL) abbt = coot::api_bond_colour_t::COLOUR_BY_CHAIN_GOODSELL;

      return abbt;
   };

   auto print_colour_table = [this] (const std::string &l) {

      std::vector<glm::vec4> colour_table = this->make_colour_table();
      std::cout << "----------- Here is the colour table for imol " << imol_no << " : " << l << " -------" << std::endl;
      for (unsigned int i=0; i<colour_table.size(); i++) {
         std::cout << "    " << std::setw(2) << i << " " << glm::to_string(colour_table[i]) << std::endl;
      }
   };

   if (false)
      std::cout << "debug:: make_meshes_from_bonds_box_instanced_version() --- start --- " << std::endl;

   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: in make_glsl_bonds_type_checked() --- start ---\n";

   if (atom_sel.mol) {

      unsigned int num_subdivisions = 2;
      unsigned int n_slices = 8;
      unsigned int n_stacks = 2;
      // do smooth
      if (graphics_info_t::bond_smoothness_factor == 1) {
         num_subdivisions = 1;
         n_slices = 8;
      }
      if (graphics_info_t::bond_smoothness_factor == 2) {
         num_subdivisions = 2;
         n_slices = 16;
      }
      if (graphics_info_t::bond_smoothness_factor == 3) {
         num_subdivisions = 3;
         n_slices = 32;
      }
      if (graphics_info_t::bond_smoothness_factor == 4) {
         num_subdivisions = 4;
         n_slices = 64;
      }
      float bond_radius = 0.026 * bond_width;
      float atom_radius = bond_radius * atom_radius_scale_factor;

      // something like this
      // float radius_scale = 0.2 * bond_width; // arbs
      // if (is_intermediate_atoms_molecule) radius_scale *= 1.8f;
      // radius_scale *= atom_radius_scale_factor;

      if (false) {
         std::cout << "DEBUG:: ********** model_representation_mode: BALL_AND_STICK " << int(Mesh::representation_mode_t::BALL_AND_STICK) << std::endl;
         std::cout << "DEBUG:: ********** model_representation_mode: BALLS_NOT_BONDS " << int(Mesh::representation_mode_t::BALLS_NOT_BONDS) << std::endl;
         std::cout << "DEBUG:: ********** model_representation_mode: VDW_BALLS " << int(Mesh::representation_mode_t::VDW_BALLS) << std::endl;
         std::cout << "DEBUG:: ********** model_representation_mode: " << int(model_representation_mode) << std::endl;
      }

      if (model_representation_mode == Mesh::representation_mode_t::BALLS_NOT_BONDS) {
         atom_radius = 1.67; // 20220226-PE  compromise between C, N, O. Actually we should of course get
                             // the radius of each atom from its type when model_representation_mode == Mesh::BALLS_NOT_BONDS.
                             // That's for another day.
      }

      std::vector<glm::vec4> colour_table = make_colour_table();
      if (false)
         print_colour_table(" ");

      err = glGetError();
      if (err) std::cout << "error in make_glsl_bonds_type_checked() pre molecules_as_mesh\n";
      float aniso_probability = graphics_info_t::show_aniso_atoms_probability;

      model_molecule_meshes.make_graphical_bonds(imol_no, bonds_box, atom_radius, bond_radius,
                                                 show_atoms_as_aniso_flag, // class member - user setable
                                                 aniso_probability,
                                                 show_aniso_atoms_as_ortep_flag, // ditto
                                                 show_aniso_atoms_as_empty_flag,
                                                 num_subdivisions, n_slices, n_stacks, colour_table);

      // Restore the user's material settings after recreating the bonds mesh
      model_molecule_meshes.set_material(material_for_models);

      // 2025-07-28 10:14 I don't want to set this here, surely.
      // There should be some other control.
      // I want to be able to update the mesh without seeing the bonds (Ctrl F)
      // if (true) // test that model_molecule_meshes is not empty()
      //    draw_it = 1;

      err = glGetError();
      if (err) std::cout << "error in make_glsl_bonds_type_checked() post molecules_as_mesh\n";
   } else {
      std::cout << "ERROR:: Null mol in make_glsl_bonds_type_checked() " << std::endl;
   }
}



 // static
glm::vec4
molecule_class_info_t::get_glm_colour_func(int idx_col, int bonds_box_type) {

   // 20220208-PE Future Paul, when you come to fix this, use vector<glm::vec4> index_to_colour rather than calling this function
   // because it's a tangle when you make it call other functions which  are not yet static (get_bond_colour_by_mol_no()).
   // i.e. don't use glm::vec4 (*get_glm_colour_for_bonds) (int, int)) in the function call.

   // 20230828-PE OK Past-Paul, I did it the way that you suggested.

   glm::vec4 col(0.7, 0.65, 0.4, 1.0);
   if (idx_col == 1) col = glm::vec4(0.7, 0.7, 0.2, 1.0);
   if (idx_col == 2) col = glm::vec4(0.7, 0.2, 0.2, 1.0);
   if (idx_col == 3) col = glm::vec4(0.2, 0.3, 0.8, 1.0);

   if (idx_col == 24) col = glm::vec4(0.9, 0.3, 0.2, 1.0);
   if (idx_col == 25) col = glm::vec4(0.2, 0.9, 0.0, 1.0);
   if (idx_col == 26) col = glm::vec4(0.2, 0.2, 0.8, 1.0);
   if (idx_col == 27) col = glm::vec4(0.5, 0.0, 0.8, 1.0);
   if (idx_col == 28) col = glm::vec4(0.5, 0.7, 0.2, 1.0);

   // bool dark_bg_flag = true; // use against_a_dark_background()?
   // coot::colour_t cc = get_bond_colour_by_mol_no(idx_col, dark_bg_flag);
   //glm::vec4 col = cc.to_glm();

   return col;
}

void molecule_class_info_t::make_glsl_bonds_type_checked(const char *caller) {

   if (false)
      std::cout << "debug:: make_glsl_bonds_type_checked() called by " << caller << "()"
                << " with is_intermediate_atoms_molecule " << is_intermediate_atoms_molecule
                << std::endl;

   if (false)
      std::cout << "debug:: ---- in make_glsl_bonds_type_checked() --- start ---" << std::endl;

   // Add no-graphics protection
   if (!graphics_info_t::use_graphics_interface_flag) return;

   GLenum err = glGetError();
   if (err) std::cout << "GL ERROR:: in make_glsl_bonds_type_checked() -- start A --\n";

   graphics_info_t::attach_buffers(); // needed.

   make_meshes_from_bonds_box_instanced_version();

   // make_mesh_from_bonds_box(); // non-instanced version - add lots of vectors of vertices and triangles
                               //
                               // needs:
                               // cis peptides,
                               // missing residue loops
                               // and rama balls if intermediate atoms.

}

void
molecule_class_info_t::make_glsl_symmetry_bonds() {

   auto pastelize_colour_table = [] (const std::vector<glm::vec4> &colour_table) {
      glm::vec4 grey(0.5, 0.5, 0.5, 1.0);
      std::vector<glm::vec4> new_colour_table = colour_table;
      for (unsigned int i=0; i<colour_table.size(); i++)
         new_colour_table[i] = (colour_table[i] + grey * 2.0f) * 0.33f;
      return new_colour_table;
   };

   // do things with symmetry_bonds_box;
   // std::vector<std::pair<graphical_bonds_container, std::pair<symm_trans_t, Cell_Translation> > > symmetry_bonds_box;
   graphics_info_t::attach_buffers();

#if 0
   mesh_for_symmetry_atoms.make_symmetry_atoms_bond_lines(symmetry_bonds_box, // boxes
                                                          graphics_info_t::symmetry_colour,
                                                          graphics_info_t::symmetry_colour_merge_weight);
#endif

   float atom_radius = 0.1;
   float bond_radius = 0.1;
   int num_subdivisions = 2;
   int n_slices = 8;
   int n_stacks = 2;
   std::vector<glm::vec4> colour_table = make_colour_table();

   std::vector<glm::vec4> new_colour_table = pastelize_colour_table(colour_table);

   meshes_for_symmetry_atoms.make_symmetry_bonds(imol_no, symmetry_bonds_box,
                                                 atom_radius, bond_radius,
                                                 num_subdivisions, n_slices, n_stacks,
                                                 new_colour_table);
}

// either we have licorice/ball-and-stick (licorice is a form of ball-and-stick) or big-ball-no-bonds
//enum { BALL_AND_STICK, BALLS_NOT_BONDS };
void
molecule_class_info_t::set_model_molecule_representation_style(unsigned int mode) {

   // we should use goodsell colouring by default here

   if (int(mode) == int(Mesh::representation_mode_t::BALL_AND_STICK)) {
      if (model_representation_mode != Mesh::representation_mode_t::BALL_AND_STICK) {
         model_representation_mode = Mesh::representation_mode_t::BALL_AND_STICK;
         make_glsl_bonds_type_checked(__FUNCTION__);
      }
   }
   if (int(mode) == int(Mesh::representation_mode_t::BALLS_NOT_BONDS)) {
      if (model_representation_mode != Mesh::representation_mode_t::BALLS_NOT_BONDS) {
         model_representation_mode = Mesh::representation_mode_t::BALLS_NOT_BONDS;
         make_glsl_bonds_type_checked(__FUNCTION__);
      }
   }
   if (int(mode) == int(Mesh::representation_mode_t::VDW_BALLS)) {
      if (model_representation_mode != Mesh::representation_mode_t::VDW_BALLS) {
         model_representation_mode = Mesh::representation_mode_t::VDW_BALLS;
         make_glsl_bonds_type_checked(__FUNCTION__);
      }
   }

}




// draw molecule as instanced meshes.
void
molecule_class_info_t::draw_molecule_as_meshes(Shader *shader_p,
                                               stereo_eye_t eye,
                                               const glm::mat4 &mvp,
                                               const glm::mat4 &view_rotation_matrix,
                                               const std::map<unsigned int, lights_info_t> &lights,
                                               const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                                               const glm::vec4 &background_colour,
                                               bool do_depth_fog) {


   if (false) {
      std::cout << "draw_molecule_as_meshes() shader " << shader_p->name << " " << std::endl;
      std::cout << "   mvp                  " << glm::to_string(mvp) << std::endl;
      std::cout << "   view_rotation_matrix " << glm::to_string(view_rotation_matrix) << std::endl;
      std::cout << "   eye pos              " << glm::to_string(eye_position) << std::endl;
   }

   Shader *shader_for_simple_mesh = &graphics_info_t::shader_for_model_as_meshes; // is this right?
   Shader *shader_for_instances = shader_p;
   float opacity = 1.0f;
   bool gl_lines_mode = false;
   bool show_just_shadows = false;
   model_molecule_meshes.draw(shader_for_simple_mesh, shader_for_instances, eye, mvp, view_rotation_matrix, lights, eye_position, opacity, background_colour,
                              gl_lines_mode, do_depth_fog, show_just_shadows);

}


// draw molecule as instanced meshes.
void
molecule_class_info_t::draw_molecule_as_meshes_for_ssao(Shader *shader_for_meshes_for_ssao,
                                                        Shader *shader_for_instanced_meshes_for_ssao,
                                                        const glm::mat4 &model_matrix,
                                                        const glm::mat4 &view_matrix,
                                                        const glm::mat4 &projection_matrix) {

   if (false) {
      std::cout << "draw_molecule_as_meshes_for_ssao() shader " << shader_for_meshes_for_ssao->name << " and "
                << shader_for_meshes_for_ssao->name << std::endl;
      std::cout << "   model_matrix " << glm::to_string(model_matrix) << std::endl;
      std::cout << "   view_matrix " << glm::to_string(view_matrix) << std::endl;
      std::cout << "   proj_matrix " << glm::to_string(projection_matrix) << std::endl;
   }

   model_molecule_meshes.draw_for_ssao(shader_for_meshes_for_ssao, shader_for_instanced_meshes_for_ssao,
                                       model_matrix, view_matrix, projection_matrix);

}

// instanced models
void
molecule_class_info_t::draw_molecule_as_meshes_with_shadows(Shader *shader,
                                                            const glm::mat4 &mvp,
                                                            const glm::mat4 &model_rotation_matrix,
                                                            const std::map<unsigned int, lights_info_t> &lights,
                                                            const glm::vec3 &eye_position, // eye position in view space (not molecule space)
                                                            float opacity,
                                                            const glm::vec4 &background_colour,
                                                            bool do_depth_fog,
                                                            const glm::mat4 &light_view_mvp,
                                                            unsigned int shadow_depthMap,
                                                            float shadow_strength,
                                                            unsigned int shadow_softness, // 1, 2 or 3.
                                                            bool show_just_shadows) {

   std::cout << "draw_molecule_as_meshes_with_shadows() replacement code needed here" << std::endl;
   // 20230909-PE if this function is bypassed, then remove this fun
}



void
molecule_class_info_t::draw_symmetry(Shader *shader_p,
                                     const glm::mat4 &mvp,
                                     const glm::mat4 &view_rotation,
                                     const std::map<unsigned int, lights_info_t> &lights,
                                     const glm::vec3 &eye_position,
                                     const glm::vec4 &background_colour,
                                     bool do_depth_fog) {

   stereo_eye_t eye = stereo_eye_t::MONO;

   if (draw_it) {
      if (show_symmetry) {
         if (this_molecule_has_crystallographic_symmetry) {

#if 0 // 20250312-PE old line symmetry
            mesh_for_symmetry_atoms.draw_symmetry(shader_p, mvp, view_rotation, lights,
                                                  eye_position, background_colour, do_depth_fog);
#endif

            Shader *shader_for_simple_mesh = &graphics_info_t::shader_for_model_as_meshes;
            Shader *shader_for_instances = &graphics_info_t::shader_for_instanced_objects;
            float opacity = 1.0;
            bool gl_lines_mode = false;
            bool show_just_shadows = false;
            meshes_for_symmetry_atoms.draw(shader_for_simple_mesh, shader_for_instances, eye,
                                           mvp, view_rotation, lights, eye_position, opacity,
                                           background_colour, gl_lines_mode, do_depth_fog,
                                           show_just_shadows);

         }
      }
   }
}



void
molecule_class_info_t::export_these_as_3d_object(const std::vector<coot::api::vertex_with_rotation_translation> &vertices,
                                                 const std::vector<g_triangle> &triangles) {

   export_vertices_and_triangles_func(vertices, triangles);
}


void
molecule_class_info_t::setup_glsl_bonds_buffers(const std::vector<coot::api::vertex_with_rotation_translation> &vertices,
                                                const std::vector<g_triangle> &triangles) {

    if (false)
       std::cout << "debug:: in setup_glsl_bonds_buffers() with vertices size " << vertices.size()
                 << " and triangles size " << triangles.size() << std::endl;

    // Nope. We need to clear the old molecule representation with nothing sometimes 20210619-PE
    // if (triangles.empty()) return;
    // if (vertices.empty()) return;

   graphics_info_t::shader_for_models.Use(); // is this called from the outside before we get here?

   n_vertices_for_model_VertexArray = vertices.size(); // the "signal" to the draw function to draw this model

   GLenum err = glGetError(); if (err) std::cout << "GL error in setup_glsl_bonds_buffers() -- start --\n";

   if (model_mesh_first_time) {
      glGenVertexArrays(1, &m_VertexArray_for_model_ID);
      err = glGetError();
      if (err) std::cout << "GL error in setup_glsl_bonds_buffers() 1\n";
   }

   glBindVertexArray(m_VertexArray_for_model_ID);
   err = glGetError();
   if (err) std::cout << "GL error in molecule_class_info_t::setup_glsl_bonds_buffers()"
                      << " glBindVertexArray() " << m_VertexArray_for_model_ID
                      << " model_mesh_first_time " << model_mesh_first_time
                      << "\n";

   if (model_mesh_first_time) {
      glGenBuffers(1, &m_VertexBuffer_for_model_ID);
      glBindBuffer(GL_ARRAY_BUFFER, m_VertexBuffer_for_model_ID);
      err = glGetError(); if (err) std::cout << "GL error in setup_glsl_bonds_buffers() 3\n";
      unsigned int n_vertices = vertices.size();
      if (is_intermediate_atoms_molecule)
         glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(vertices[0]), &(vertices[0]), GL_DYNAMIC_DRAW);
      else
         glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(vertices[0]), &(vertices[0]), GL_STATIC_DRAW);
      err = glGetError(); if (err) std::cout << "GL error in setup_glsl_bonds_buffers()  5\n";
   } else {
      glDeleteBuffers(1, &m_VertexBuffer_for_model_ID);
      glGenBuffers(1, &m_VertexBuffer_for_model_ID);
      glBindBuffer(GL_ARRAY_BUFFER, m_VertexBuffer_for_model_ID);
      err = glGetError(); if (err) std::cout << "GL error in setup_glsl_bonds_buffers() 3\n";
      unsigned int n_vertices = vertices.size();
      if (is_intermediate_atoms_molecule)
         glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(vertices[0]), &(vertices[0]), GL_DYNAMIC_DRAW);
      else
         glBufferData(GL_ARRAY_BUFFER, n_vertices * sizeof(vertices[0]), &(vertices[0]), GL_STATIC_DRAW);
      err = glGetError(); if (err) std::cout << "GL error in setup_glsl_bonds_buffers()  5\n";
   }

   // "from-origin" model matrix (orientation)
   glEnableVertexAttribArray(0);
   glEnableVertexAttribArray(1);
   glEnableVertexAttribArray(2);
   glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(coot::api::vertex_with_rotation_translation), reinterpret_cast<void *>(0 * sizeof(glm::vec3)));
   glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(coot::api::vertex_with_rotation_translation), reinterpret_cast<void *>(1 * sizeof(glm::vec3)));
   glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(coot::api::vertex_with_rotation_translation), reinterpret_cast<void *>(2 * sizeof(glm::vec3)));

   // "from origin" translate position, 3, size 3 floats
   glEnableVertexAttribArray(3);
   glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(coot::api::vertex_with_rotation_translation), reinterpret_cast<void *>(3 * sizeof(glm::vec3)));

   // positions, 4, size 3 floats
   glEnableVertexAttribArray(4);
   err = glGetError(); if (err) std::cout << "GL error in setup_glsl_bonds_buffers() 6\n";
   glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, sizeof(coot::api::vertex_with_rotation_translation), reinterpret_cast<void *>(4 * sizeof(glm::vec3)));
   err = glGetError(); if (err) std::cout << "GL error bonds 7\n";

   //  normals, 5, size 3 floats
   glEnableVertexAttribArray(5);
   err = glGetError(); if (err) std::cout << "GL error in setup_glsl_bonds_buffers() 11\n";
   glVertexAttribPointer(5, 3, GL_FLOAT, GL_FALSE, sizeof(coot::api::vertex_with_rotation_translation), reinterpret_cast<void *>(5 * sizeof(glm::vec3)));
   err = glGetError(); if (err) std::cout << "GL error in setup_glsl_bonds_buffers() 12\n";

   //  colours, 6, size 4 floats
   glEnableVertexAttribArray(6);
   err = glGetError(); if (err) std::cout << "GL error setup_glsl_bonds_buffers()  16\n";
   glVertexAttribPointer(6, 4, GL_FLOAT, GL_FALSE, sizeof(coot::api::vertex_with_rotation_translation), reinterpret_cast<void *>(6 * sizeof(glm::vec3)));
   err = glGetError(); if (err) std::cout << "GL error bonds 17\n";

   // Indices
   n_indices_for_model_triangles = triangles.size() * 3;
   // if (! is_intermediate_atoms_molecule)
   // std::cout << "DEBUG:: n_triangles in model: " << triangles.size() << std::endl;
   unsigned int n_bytes = triangles.size() * 3 * sizeof(unsigned int);
   if (model_mesh_first_time) {
      glGenBuffers(1, &m_IndexBuffer_for_model_ID);
      err = glGetError(); if (err) std::cout << "GL error bonds setup_glsl_bonds_buffers() 18\n";
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_IndexBuffer_for_model_ID);
      err = glGetError(); if (err) std::cout << "GL error bonds setup_glsl_bonds_buffers() 19\n";
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_bytes, &triangles[0], GL_STATIC_DRAW);
   } else {
      glDeleteBuffers(1, &m_IndexBuffer_for_model_ID);
      glGenBuffers(1, &m_IndexBuffer_for_model_ID);
      err = glGetError(); if (err) std::cout << "GL error bonds setup_glsl_bonds_buffers() 18\n";
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_IndexBuffer_for_model_ID);
      err = glGetError(); if (err) std::cout << "GL error bonds setup_glsl_bonds_buffers() 19\n";
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, n_bytes, &triangles[0], GL_STATIC_DRAW);
   }
   err = glGetError(); if (err) std::cout << "GL error bonds --- end ---\n";
   model_mesh_first_time = false;

}

// caller is an optional argument
void
molecule_class_info_t::make_bonds_type_checked(const std::set<int> &no_bonds_to_these_atom_indices,
                                               const char *caller) {

   if (true)
      std::cout << "debug:: ---- in make_bonds_type_checked(2args) --- start ---" << std::endl;

   if (false) {
      std::string caller_s = "NULL";
      if (caller)
         caller_s = std::string(caller);
      std::cout << "debug::make_bonds_type_checked() no-bonds-to-these_atoms-set size "
                << no_bonds_to_these_atom_indices.size() << " called by "
                << caller_s << std::endl;
   }

   if (false) {
      for (auto it=no_bonds_to_these_atom_indices.begin();
           it!=no_bonds_to_these_atom_indices.end(); ++it) {
         mmdb::Atom *at = atom_sel.atom_selection[*it];
         std::cout << "   make_bonds_type_checked(): No bond to atom " << coot::atom_spec_t(at) << std::endl;
      }
   }

   bool force_rebonding = true; // if we are here then we need to do rebonding (I think)

   graphics_info_t g; // urgh!  (But the best solution?)
   coot::protein_geometry *geom_p = g.Geom_p();
   if (bonds_box_type == coot::NORMAL_BONDS) {
      makebonds(geom_p, no_bonds_to_these_atom_indices);
   } else {
      if (bonds_box_type == coot::COLOUR_BY_CHAIN_BONDS) {
         bool goodsell_mode = false;
         make_colour_by_chain_bonds(no_bonds_to_these_atom_indices,
                                    g.rotate_colour_map_on_read_pdb_c_only_flag,
                                    goodsell_mode, force_rebonding);
      } else {
         if (bonds_box_type == coot::COLOUR_BY_CHAIN_GOODSELL) {
            bool goodsell_mode = true;
            make_colour_by_chain_bonds(no_bonds_to_these_atom_indices,
                                       g.rotate_colour_map_on_read_pdb_c_only_flag,
                                       goodsell_mode, force_rebonding);
         } else {
            if (bonds_box_type == coot::CA_BONDS) {
               make_ca_bonds(2.4, 4.7, no_bonds_to_these_atom_indices);
            } else {
               if (bonds_box_type == coot::CA_BONDS_PLUS_LIGANDS) {
                  make_ca_bonds(2.4, 4.7, no_bonds_to_these_atom_indices);
               } else {
                  make_bonds_type_checked(__FUNCTION__); // function above
               }
            }
         }
      }
   }
}


void
molecule_class_info_t::update_bonds_using_phenix_geo(const coot::phenix_geo::phenix_geometry &phenix_geo) {

   Bond_lines_container bonds(atom_sel.mol, phenix_geo);
   bonds_box.clear_up();
   bonds_box = bonds.make_graphical_bonds();
   make_glsl_bonds_type_checked(__FUNCTION__);

}


// n-models
int
molecule_class_info_t::n_models() const {
   int r = -1;
   if (has_model()) {
      r = atom_sel.mol->GetNumberOfModels();
   }
   return r;
}

// single model view
void
molecule_class_info_t::single_model_view_model_number(int imodel) {
   if (has_model()) {
      single_model_view_current_model_number = imodel;
      // make_bonds_type_checked();
   }
}

int
molecule_class_info_t::single_model_view_this_model_number() const {
   int model_no = 0;
   if (has_model()) {
      model_no = single_model_view_current_model_number;
   }
   return model_no;
}

int
molecule_class_info_t::single_model_view_prev_model_number() {
   int model_no = 0;
   if (has_model()) {
      int n = n_models();
      if (n > 1) {
         int prev = single_model_view_current_model_number - 1;
         if (prev >= 1) {
            // OK
         } else {
            prev = n;
         }
         mmdb::Model *model = atom_sel.mol->GetModel(prev);
         if (model) {
            model_no = prev;
         }
      }
   }
   single_model_view_model_number(model_no);
   return model_no;
}

int
molecule_class_info_t::single_model_view_next_model_number() {
   int model_no = 0;
   if (has_model()) {
      int n = n_models();
      if (n > 1) {
         int next = single_model_view_current_model_number + 1;
         if (next <= n) {
            // OK
         } else {
            next = 1;
         }
         mmdb::Model *model = atom_sel.mol->GetModel(next);
         if (model) {
            model_no = next;
         }
      }
   }
   single_model_view_model_number(model_no);
   return model_no;
}

#ifndef EMSCRIPTEN
void
molecule_class_info_t::update_additional_representations(const gl_context_info_t &gl_info,
                                                         const coot::protein_geometry *geom) {

   for (unsigned int i=0; i<add_reps.size(); i++) {
      // make_ball_and_stick is not available from inside an add_rep,
      // so we do it outside.
      if (add_reps[i].representation_type != coot::BALL_AND_STICK) {
         add_reps[i].update_self();
      } else {

         // let's remove the old/current dlo from the tags list before we add a new one.

         int old_handle = add_reps[i].display_list_handle;
         remove_display_list_object_with_handle(old_handle);

         int handle = make_ball_and_stick(add_reps[i].atom_sel_info.mmdb_string(),
                                          add_reps[i].bond_width,
                                          add_reps[i].sphere_radius,
                                          add_reps[i].draw_atom_spheres_flag,
                                          gl_info, geom);
         int n_display_list_tags = display_list_tags.size();
         if ((handle >= 0) && (handle < n_display_list_tags))
            add_reps[i].update_self_display_list_entity(handle);
         display_list_tags[handle].display_it = add_reps[i].show_it;
      }
   }
}
#endif


void
molecule_class_info_t::remove_display_list_object_with_handle(int handle_index) {

//    for (unsigned int i=0; i<display_list_tags.size(); i++) {
//       std::cout << "   display_list_tags: index " << i << " tag " << display_list_tags[i].tag;
//       if (display_list_tags[i].is_closed)
//          std::cout << " closed";
//       std::cout << std::endl;
//    }

//    std::cout << "closing display_list_tags number..." << handle_index
//              << " which has GL tag " << display_list_tags[handle_index].tag << std::endl;
   display_list_tags[handle_index].close_yourself();
}

void
molecule_class_info_t::update_mols_in_additional_representations() {

   // std::cout << "........................... update_mols_in_additional_representations() called..."
   // << std::endl;

   for (unsigned int i=0; i<add_reps.size(); i++) {
      add_reps[i].change_mol(atom_sel.mol);
   }
}


void
molecule_class_info_t::update_fixed_atom_positions() {

   fixed_atom_positions.clear();
   bool found_match = 0;
   for(unsigned int i=0; i<fixed_atom_specs.size(); i++) {
      int ifast_index = fixed_atom_specs[i].int_user_data;
      if (ifast_index != -1) {
         if (ifast_index < atom_sel.n_selected_atoms) {
            mmdb::Atom *at = atom_sel.atom_selection[ifast_index];
            if (fixed_atom_specs[i].matches_spec(at)) {
               found_match = 1;
               coot::Cartesian pos(at->x, at->y, at->z);
               fixed_atom_positions.push_back(pos);
            }
         }
      }
      if (! found_match) {
         // use a slower method to find atom
         int idx = full_atom_spec_to_atom_index(fixed_atom_specs[i]);
         if (idx != -1) {
            mmdb::Atom *at = atom_sel.atom_selection[idx];
            if (fixed_atom_specs[i].matches_spec(at)) {
               coot::Cartesian pos(at->x, at->y, at->z);
               fixed_atom_positions.push_back(pos);
            }
         }
      }
   }
}

std::vector<coot::atom_spec_t>
molecule_class_info_t::get_fixed_atoms() const {
   return fixed_atom_specs;
}

void
molecule_class_info_t::update_extra_restraints_representation() {

   if (false)
      std::cout << "here we are in update_extra_restraints_representation() "
                << extra_restraints.bond_restraints.size() << " "
                << extra_restraints.geman_mcclure_restraints.size() << std::endl;

   extra_restraints_representation.clear();
   update_extra_restraints_representation_bonds();
   update_extra_restraints_representation_geman_mcclure();
   update_extra_restraints_representation_parallel_planes();

}

void
molecule_class_info_t::update_extra_restraints_representation_bonds() {

   // extra_restraints_representation.clear() should be called before calling this function.

   // make things redraw fast - this is a hack for morph-and-refine.

   if (false)
      std::cout << "here with extra_restraints_representation.bond_restraints size "
                << extra_restraints.bond_restraints.size() << " " << draw_it_for_extra_restraints
                << std::endl;

   // I want  to update them even if they are not drawn, I think.
   // if (! draw_it_for_extra_restraints || ! draw_it)
   //      return;

   for (unsigned int i=0; i<extra_restraints.bond_restraints.size(); i++) {
      mmdb::Atom *at_1 = NULL;
      mmdb::Atom *at_2 = NULL;
      clipper::Coord_orth p1(0,0,0);
      clipper::Coord_orth p2(0,0,0);
      bool ifound_1 = false;
      bool ifound_2 = false;
      int ifast_index_1 = extra_restraints.bond_restraints[i].atom_1.int_user_data;
      int ifast_index_2 = extra_restraints.bond_restraints[i].atom_2.int_user_data;
      const coot::extra_restraints_t::extra_bond_restraint_t &res = extra_restraints.bond_restraints[i];

      // set p1 from ifast_index_1 (if possible)
      //
      if (ifast_index_1 != -1) {
         if (ifast_index_1 < atom_sel.n_selected_atoms) {
            at_1 = atom_sel.atom_selection[ifast_index_1];
            if (extra_restraints.bond_restraints[i].atom_1.matches_spec(at_1)) {
               p1 = clipper::Coord_orth(at_1->x, at_1->y, at_1->z);
               ifound_1 = true;
            }
         }
      }
      if (! ifound_1) {
         int idx = full_atom_spec_to_atom_index(extra_restraints.bond_restraints[i].atom_1);
         if (idx != -1) {
            at_1 = atom_sel.atom_selection[idx];
            if (extra_restraints.bond_restraints[i].atom_1.matches_spec(at_1)) {
               p1 = clipper::Coord_orth(at_1->x, at_1->y, at_1->z);
               ifound_1 = true;
            }
         }
      }
   }

}

void
molecule_class_info_t::update_extra_restraints_representation_geman_mcclure() {

   for (unsigned int i=0; i<extra_restraints.geman_mcclure_restraints.size(); i++) {
      const coot::extra_restraints_t::extra_geman_mcclure_restraint_t &rest = extra_restraints.geman_mcclure_restraints[i];
      mmdb::Atom *at_1 = NULL;
      mmdb::Atom *at_2 = NULL;
      clipper::Coord_orth p1(0,0,0);
      clipper::Coord_orth p2(0,0,0);
      bool ifound_1 = false;
      bool ifound_2 = false;
      int ifast_index_1 = rest.atom_1.int_user_data;
      int ifast_index_2 = rest.atom_2.int_user_data;

      if (ifast_index_1 != -1) {
         if (ifast_index_1 < atom_sel.n_selected_atoms) {
            at_1 = atom_sel.atom_selection[ifast_index_1];
            if (rest.atom_1.matches_spec(at_1)) {
               p1 = clipper::Coord_orth(at_1->x, at_1->y, at_1->z);
               ifound_1 = true;
            }
         }
      }
      if (! ifound_1) {
         int idx = full_atom_spec_to_atom_index(rest.atom_1);
         if (idx != -1) {
            if (idx < atom_sel.n_selected_atoms) {
               at_1 = atom_sel.atom_selection[idx];
               if (rest.atom_1.matches_spec(at_1)) {
                  p1 = clipper::Coord_orth(at_1->x, at_1->y, at_1->z);
                  ifound_1 = true;
               }
            }
         }
      }
      if (ifast_index_2 != -1) {
         if (ifast_index_2 < atom_sel.n_selected_atoms) {
            at_2 = atom_sel.atom_selection[ifast_index_2];
            if (rest.atom_2.matches_spec(at_2)) {
               p2 = clipper::Coord_orth(at_2->x, at_2->y, at_2->z);
               ifound_2 = true;
            }
         }
      }
      if (! ifound_2) {
         int idx = full_atom_spec_to_atom_index(rest.atom_1);
         if (idx != -1) {
            if (idx < atom_sel.n_selected_atoms) {
               at_1 = atom_sel.atom_selection[idx];
               if (rest.atom_2.matches_spec(at_2)) {
                  p2 = clipper::Coord_orth(at_2->x, at_2->y, at_2->z);
                  ifound_2 = true;
               }
            }
         }
      }

      if (ifound_1 && ifound_2) {

         // if the distance (actually, n-sigma) is within limits, draw it.
         //
         double dist_sq = (p1-p2).lengthsq();
         double dist = sqrt(dist_sq);
         double this_n_sigma = (dist - rest.bond_dist)/rest.esd;

         // std::cout << "dist " << dist << " rest.bond_dist " << rest.bond_dist << std::endl;

         if (false)
            std::cout << "comparing this_n_sigma " << this_n_sigma << " with "
                      << extra_restraints_representation.prosmart_restraint_display_limit_low << " "
                      << extra_restraints_representation.prosmart_restraint_display_limit_high << " and to-CA-mode "
                      << extra_restraints_representation_for_bonds_go_to_CA
                      << "\n";

         if (this_n_sigma >= extra_restraints_representation.prosmart_restraint_display_limit_high ||
             this_n_sigma <= extra_restraints_representation.prosmart_restraint_display_limit_low) {

            extra_restraints_representation.add_bond(p1, p2, rest.bond_dist, rest.esd);
         }
      }
   }
}


void
molecule_class_info_t::update_extra_restraints_representation_bonds_internal(const coot::extra_restraints_t::extra_bond_restraint_t &res) {

   mmdb::Atom *at_1 = NULL;
   mmdb::Atom *at_2 = NULL;
   clipper::Coord_orth p1(0,0,0);
   clipper::Coord_orth p2(0,0,0);
   bool ifound_1 = false;
   bool ifound_2 = false;
   int ifast_index_1 = res.atom_1.int_user_data;
   int ifast_index_2 = res.atom_2.int_user_data;

   // set p1 from ifast_index_1 (if possible)
   //
   if (ifast_index_1 != -1) {
      if (ifast_index_1 < atom_sel.n_selected_atoms) {
         at_1 = atom_sel.atom_selection[ifast_index_1];
         if (res.atom_1.matches_spec(at_1)) {
            p1 = clipper::Coord_orth(at_1->x, at_1->y, at_1->z);
            ifound_1 = true;
         }
      }
   }
   if (! ifound_1) {
      int idx = full_atom_spec_to_atom_index(res.atom_1);
      if (idx != -1) {
         at_1 = atom_sel.atom_selection[idx];
         if (res.atom_1.matches_spec(at_1)) {
            p1 = clipper::Coord_orth(at_1->x, at_1->y, at_1->z);
            ifound_1 = true;
         }
      }
   }

   // set p2 from ifast_index_2 (if possible)
   //
   if (ifast_index_2 != -1) {
      if (ifast_index_2 < atom_sel.n_selected_atoms) {
         at_2 = atom_sel.atom_selection[ifast_index_2];
         if (res.atom_2.matches_spec(at_2)) {
            p2 = clipper::Coord_orth(at_2->x, at_2->y, at_2->z);
            ifound_2 = true;
         }
      }
   }
   if (! ifound_2) {
      int idx = full_atom_spec_to_atom_index(res.atom_2);
      if (idx != -1) {
         at_2 = atom_sel.atom_selection[idx];
         if (res.atom_2.matches_spec(at_2)) {
            p2 = clipper::Coord_orth(at_2->x, at_2->y, at_2->z);
            ifound_2 = 1;
         }
      }
   }

   // std::cout << "post: debug ifound_1 and ifound_2 " << ifound_1 << " " << ifound_2 << std::endl;

   if (! ifound_1) {
      std::cout << "no spec for " << res.atom_1 << std::endl;
   }

   if (ifound_1 && ifound_2) {

      // if the distance (actually, n-sigma) is within limits, draw it.
      //
      double dist_sq = (p1-p2).lengthsq();
      double dist = sqrt(dist_sq);
      double this_n_sigma = (dist - res.bond_dist)/res.esd;

      if (0)
         std::cout << "comparing this_n_sigma " << this_n_sigma << " with "
                   << extra_restraints_representation.prosmart_restraint_display_limit_low << " "
                   << extra_restraints_representation.prosmart_restraint_display_limit_high << " and to-CA-mode "
                   << extra_restraints_representation_for_bonds_go_to_CA
                   << "\n";

      if (this_n_sigma >= extra_restraints_representation.prosmart_restraint_display_limit_high ||
          this_n_sigma <= extra_restraints_representation.prosmart_restraint_display_limit_low) {
         if (extra_restraints_representation_for_bonds_go_to_CA) {
            if (at_1->residue != at_2->residue) {
               int idx_1 = intelligent_this_residue_atom(at_1->residue);
               int idx_2 = intelligent_this_residue_atom(at_2->residue);
               if (idx_1 >= 0 && idx_2 >= 0) {
                  clipper::Coord_orth ca_p1 = coot::co(atom_sel.atom_selection[idx_1]);
                  clipper::Coord_orth ca_p2 = coot::co(atom_sel.atom_selection[idx_2]);
                  // hack a distance.
                  double d = sqrt((ca_p2-ca_p1).lengthsq());

                  // undashed:
                  // extra_restraints_representation.add_bond(ca_p1, ca_p2, d, res.esd);

                  // dashed:
                  //
                  double dash_density = 4.0;
                  int n_dashes = int(dash_density * d);
                  bool visible = true;
                  for (int idash=0; idash<(n_dashes-1); idash++) {
                     if (0)
                        std::cout << "idash " << idash << " n_dashes " << n_dashes << " visible "
                                  << visible << std::endl;
                     if (visible) {
                        double frac_s = double(idash  )/double(n_dashes);
                        double frac_e = double(idash+1)/double(n_dashes);
                        clipper::Coord_orth dash_pos_1 = ca_p1 + frac_s * (ca_p2 - ca_p1);
                        clipper::Coord_orth dash_pos_2 = ca_p1 + frac_e * (ca_p2 - ca_p1);
                        std::cout << "   " << dash_pos_1.format() << " " << dash_pos_2.format() << " "
                                  << d << " " << res.esd << std::endl;
                        double fake_d = d/double(n_dashes);
                        extra_restraints_representation.add_bond(dash_pos_1, dash_pos_2, fake_d, res.esd);
                     }
                     visible = !visible;
                  }
               }
            }
         } else {
            // Normal case - not CA exception
            extra_restraints_representation.add_bond(p1, p2, res.bond_dist, res.esd);
         }
      }
   }
}

void
molecule_class_info_t::update_extra_restraints_representation_parallel_planes() {

   std::string s = "Mol " + coot::util::int_to_string(imol_no) + " Parallel-Plane-Restraints";
   coot::old_generic_display_object_t rest_rep(s);

   // needed for parallel plane restraints atom name lookup.
   const coot::protein_geometry &geom = *graphics_info_t::Geom_p();  // pass this?

   for (unsigned int i=0; i<extra_restraints.parallel_plane_restraints.size(); i++) {

      const coot::parallel_planes_t &pp = extra_restraints.parallel_plane_restraints[i];
      mmdb::Residue *r_1 = get_residue(pp.plane_1_atoms.res_spec);
      mmdb::Residue *r_2 = get_residue(pp.plane_2_atoms.res_spec);

      if (r_1 && r_2) {

    std::string res_type_1 = r_1->GetResName();
    std::string res_type_2 = r_2->GetResName();
    std::pair<bool, coot::dictionary_residue_restraints_t> dri_1 = geom.get_monomer_restraints(res_type_1, imol_no);
    std::pair<bool, coot::dictionary_residue_restraints_t> dri_2 = geom.get_monomer_restraints(res_type_2, imol_no);

    std::vector<clipper::Coord_orth> p_1_positions;
    std::vector<clipper::Coord_orth> p_2_positions;

    mmdb::PPAtom residue_atoms = 0;
    int n_residue_atoms;
    r_1->GetAtomTable(residue_atoms, n_residue_atoms);

    // first plane
    for (unsigned int i_rest_at=0; i_rest_at<pp.plane_1_atoms.atom_names.size(); i_rest_at++) {
       std::string plane_atom_expanded_name =
          dri_1.second.atom_name_for_tree_4c(pp.plane_1_atoms.atom_names[i_rest_at]);
       for (int iat=0; iat<n_residue_atoms; iat++) {
          mmdb::Atom *at = residue_atoms[iat];
          std::string atom_name(at->name);
          std::string alt_conf(at->altLoc);
          if (plane_atom_expanded_name == atom_name) {
             if (pp.plane_1_atoms.alt_conf == alt_conf) {
                p_1_positions.push_back(clipper::Coord_orth(at->x, at->y, at->z));
             }
          }
       }
    }
    // second plane
   residue_atoms = 0;
   r_2->GetAtomTable(residue_atoms, n_residue_atoms);
   std::cout << "second plane " << coot::residue_spec_t(r_2) << std::endl;
   for (unsigned int i_rest_at=0; i_rest_at<pp.plane_2_atoms.atom_names.size(); i_rest_at++) {
      std::string plane_atom_expanded_name =
          dri_2.second.atom_name_for_tree_4c(pp.plane_2_atoms.atom_names[i_rest_at]);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         std::string atom_name(at->name);
         std::string alt_conf(at->altLoc);
         if (false)
            std::cout << "comparing ref-plane-name \"" << plane_atom_expanded_name << "\" and molecule atom name \""
                      << atom_name << "\"" << std::endl;
         if (plane_atom_expanded_name == atom_name) {
            if (pp.plane_2_atoms.alt_conf == alt_conf) {
               p_2_positions.push_back(clipper::Coord_orth(at->x, at->y, at->z));
            }
         }
       }
    }

    if (p_1_positions.size() > 2) {
       if (p_2_positions.size() > 2) {
          coot::lsq_plane_info_t pi_1(p_1_positions);
          coot::lsq_plane_info_t pi_2(p_2_positions);
          extra_restraints_representation.add_parallel_plane(pi_1, pi_2);
       }
    }
      }
   }
}


// export the molecule in atom_selection_container_t atom_sel;
//
int
molecule_class_info_t::export_coordinates(std::string filename) const {

   //
   int err = atom_sel.mol->WritePDBASCII(filename.c_str());

   if (err) {
      std::cout << "WARNING:: export coords: There was an error in writing "
                << filename << std::endl;
      std::cout << mmdb::GetErrorDescription(mmdb::ERROR_CODE(err)) << std::endl;
      graphics_info_t g;
      std::string s = "ERROR:: writing coordinates file ";
      s += filename;
      g.add_status_bar_text(s);
   } else {
      std::string s = "INFO:: coordinates file ";
      s += filename;
      s += " saved successfully";
      graphics_info_t g;
      g.add_status_bar_text(s);
   }
   return err;
}

// Perhaps this should be a util function?
mmdb::Manager *
molecule_class_info_t::get_residue_range_as_mol(const std::string &chain_id,
                                                int resno_start,
                                                int resno_end) const {

   int imod = 1;

   mmdb::Manager *mol_new = new mmdb::Manager;
   mmdb::Model *model_new = new mmdb::Model;
   mmdb::Chain *chain_new = new mmdb::Chain;

   mmdb::realtype cell[6];
   mmdb::realtype vol;
   int orthcode;
   char *spacegroup_str = atom_sel.mol->GetSpaceGroup();
   atom_sel.mol->GetCell(cell[0], cell[1], cell[2],
                         cell[3], cell[4], cell[5],
                         vol, orthcode);
   mol_new->SetCell(cell[0], cell[1], cell[2],
                    cell[3], cell[4], cell[5], orthcode);
   mol_new->SetSpaceGroup(spacegroup_str);

   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   mmdb::Chain *chain_p;
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      if (std::string(chain_p->GetChainID()) == chain_id) {
         int nres = chain_p->GetNumberOfResidues();
         mmdb::PResidue residue_p;
         for (int ires=0; ires<nres; ires++) {
            residue_p = chain_p->GetResidue(ires);
            if (residue_p->GetSeqNum() >= resno_start) {
               if (residue_p->GetSeqNum() <= resno_end) {
                  mmdb::Residue *res_new =
                     coot::util::deep_copy_this_residue(residue_p);
                  chain_new->AddResidue(res_new);
               }
            }
         }
      }
   }

   chain_new->SetChainID(chain_id.c_str());
   model_new->AddChain(chain_new);
   mol_new->AddModel(model_new);
   mol_new->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
   mol_new->FinishStructEdit();
   return mol_new;
}


std::string
molecule_class_info_t::make_symm_atom_label_string(mmdb::PAtom atom,
                                                   const std::pair <symm_trans_t, Cell_Translation> &symm_trans) const {

   std::string s = make_atom_label_string(atom, 0, 0);
   s += ": ";
   s += to_string(symm_trans);
   return s;
}

// This took one minute to write.
//
// And then another 2 to add the alt location code (and test it).
//
std::string
molecule_class_info_t::make_atom_label_string(mmdb::PAtom atom,
                                              int brief_atom_labels_flag,
                                              short int seg_ids_in_atom_labels_flag) const {

   char *chain_id  = atom->GetChainID();
   char *res_name  = atom->GetResName();
   int   res_no    = atom->GetSeqNum();
   char *atom_name = atom->name;
   char *ins_code  = atom->GetInsCode();

   // format: atom_name/res_no res_name/chain_id
   // new format: atom_name,alt_conf/res_no res_name/chain_id if altconf != ""
   graphics_info_t g;

   std::string s(atom_name);
   std::string alt_loc(atom->altLoc);
   if (alt_loc != "") {
      int slen = s.length();
      if (slen > 0) {
         if (s[slen-1] == ' ') {
            s = s.substr(0,slen-1) + ",";
         } else {
            s += ",";
         }
      } else {
         s += ",";
      }
      s += alt_loc;
   }

   if (brief_atom_labels_flag) {
      s += g.int_to_string(res_no);
      if (strlen(ins_code) > 0) {
         s += ins_code;
         s += " ";
      }
      s += chain_id;
   } else {
      s += "/";
      s += g.int_to_string(res_no);
      s += ins_code;
      s += " ";
      s += res_name;
      s += "/";
      s += chain_id;
      if (seg_ids_in_atom_labels_flag) {
         std::string seg_id(atom->segID);
         // std::cout << "Here with seg_id :" << seg_id << ":" << std::endl;
         if (! seg_id.empty()) {
            s += " ";
            s += seg_id;
         }
      }
   }

   return s;
}

// Don't use this function.
//
// For labelling from guile and replacing coordinates and others.
//
// Return -1 on failure to find match
int
molecule_class_info_t::atom_spec_to_atom_index(std::string chain, int resno,
                                               std::string atom_name) const {

   int iatom_index = -1;
   int selHnd = atom_sel.mol->NewSelection();

   atom_sel.mol->SelectAtoms(selHnd, 0, chain.c_str(),
                             resno, "*", // start, insertion code
                             resno, "*", // end, insertion code
                             "*", // residue name
                             atom_name.c_str(),
                             "*", // elements
                             "*"); // alt locs

   int nSelAtoms;
   mmdb::PPAtom local_SelAtom;
   atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);

   if (nSelAtoms == 0) {
      std::cout << "Sorry (atom_spec_to_atom_index): Could not find " << atom_name << "/"
                << resno << "/" << chain << " in this molecule: ("
                <<  imol_no << ") " << name_ << std::endl;

      // debug:
      selHnd = atom_sel.mol->NewSelection();

      atom_sel.mol->SelectAtoms(selHnd, 0, "*",
                                mmdb::ANY_RES, "*", // start, insertion code
                                mmdb::ANY_RES, "*", // end, insertion code
                                "*", // residue name
                                atom_name.c_str(),
                                "*", // elements
                                "*"); // alt locs

      atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);

      std::cout << "There were " << nSelAtoms << " atoms with resno "
                << resno << std::endl;


   } else {
      // compare pointers
      for (int i=0; i<atom_sel.n_selected_atoms; i++) {
         if (atom_sel.atom_selection[i] == local_SelAtom[0]) {
            iatom_index = i;
            break;
         }
      }
   }
   return iatom_index;
}

int
molecule_class_info_t::atom_to_atom_index(mmdb::Atom *at) const {

   int iatom_index_udd = -1;
   int ic;
   if (at->GetUDData(atom_sel.UDDAtomIndexHandle, ic) == mmdb::UDDATA_Ok) {
      iatom_index_udd = ic;
   }
   if (iatom_index_udd == -1)
      iatom_index_udd=full_atom_spec_to_atom_index(coot::atom_spec_t(at));
   return iatom_index_udd;
}

int
molecule_class_info_t::full_atom_spec_to_atom_index(const coot::atom_spec_t &atom_spec) const {

   return full_atom_spec_to_atom_index(atom_spec.chain_id,
                                       atom_spec.res_no,
                                       atom_spec.ins_code,
                                       atom_spec.atom_name,
                                       atom_spec.alt_conf);

}


// return -1 on no atom found.
int
molecule_class_info_t::full_atom_spec_to_atom_index(const std::string &chain,
                                                    int resno,
                                                    const std::string &insertion_code,
                                                    const std::string &atom_name,
                                                    const std::string &alt_conf) const {

   int iatom_index = -1;

   // some protection for null molecule.
   if (! atom_sel.mol) {
      std::cout << "ERROR:: null molecule for molecule number "
                << imol_no << " pointer: " << atom_sel.mol
                << " (in full_atom_spec_to_atom_index)" << std::endl;
      return -1;
   }

   int selHnd = atom_sel.mol->NewSelection();
   int idx = 0;

   atom_sel.mol->SelectAtoms(selHnd, 0, chain.c_str(),
                            resno, insertion_code.c_str(), // start, insertion code
                            resno, insertion_code.c_str(), // end, insertion code
                            "*", // residue name
                            atom_name.c_str(),
                            "*", // elements
                            alt_conf.c_str()); // alt locs

   int nSelAtoms;
   mmdb::PPAtom local_SelAtom;
   atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);

   if (false)
      std::cout << "DEBUG:: full_atom_spec_to_atom_index() for :" << chain << ": "
                << resno << " :" << insertion_code << ": :"
                << atom_name << ": :" << alt_conf << ": finds " << nSelAtoms <<  " atoms\n";

   if (nSelAtoms == 0) {

      std::cout << "WARNING:: full_atom_spec_to_atom_index() Could not find "
                << "\"" << atom_name << "\"," << "\"" << alt_conf  << "\"" << "/"
                << resno << insertion_code << "/" << chain << " in this molecule: ("
                <<  imol_no << ") " << name_ << std::endl;

      int selHnd2 = atom_sel.mol->NewSelection(); // d

      atom_sel.mol->SelectAtoms(selHnd2, 0,
                                chain.c_str(),
                                resno, "*", // start, insertion code
                                resno, "*", // end, insertion code
                                "*", // residue name
                                "*", // atom name
                                "*", // elements
                                "*"); // alt locs

      atom_sel.mol->GetSelIndex(selHnd2, local_SelAtom, nSelAtoms);

      if (false) { // debugging.
         std::cout << "There were " << nSelAtoms << " atoms in that residue:\n";
         std::cout << "debgu:: full_atom_spec_to_atom_index() resno " << resno
                   << " (cf MinInt4) " << mmdb::MinInt4 << "\n";
         if (resno == mmdb::MinInt4) {
            std::cout << "      residue with resno MinInt4\n";
         } else {
            for (int i=0; i<nSelAtoms; i++) {
               std::cout << "      " << local_SelAtom[i] << "\n";
            }
         }
      }

      atom_sel.mol->DeleteSelection(selHnd2);

   } else {

      if (nSelAtoms != 1) {
         // the wildcard atom selection case "*HO2"
         short int found = 0;
         for (int i=0; i<nSelAtoms; i++) {
            if (std::string(local_SelAtom[i]->GetChainID()) == chain) {
               if (local_SelAtom[i]->residue->seqNum == resno) {
                  if (std::string(local_SelAtom[i]->GetInsCode()) == insertion_code) {
                     if (std::string(local_SelAtom[i]->name) == atom_name) {
                        if (std::string(local_SelAtom[i]->altLoc) == alt_conf) {
                           found = 0;
                           idx = i;
                           break;
                        }
                     }
                  }
               }
            }
         }
      }

      int iatom_index_udd = -1;
      int ic;
      if (local_SelAtom[idx]->GetUDData(atom_sel.UDDAtomIndexHandle, ic) == mmdb::UDDATA_Ok) {
         iatom_index_udd = ic;
      }
      iatom_index = iatom_index_udd;
   }
   atom_sel.mol->DeleteSelection(selHnd); // Oh dear, this should have
                                          // been in place for years
                                          // (shouldn't it?) 20071121
   return iatom_index;
}
//       // compare pointers
//       for (int i=0; i<atom_sel.n_selected_atoms; i++) {
//          if (atom_sel.atom_selection[i] == local_SelAtom[0]) {
//             iatom_index = i;
//             break;
//          }
//       }
         //          if (iatom_index != iatom_index_udd) {
         //             std::cout << "ERROR: atom indexes (UDD test) dont match "
         //                       << iatom_index << " " << iatom_index_udd << std::endl;
         //          }


// Does atom at from moving atoms match atom_sel.atom_selection[this_mol_index_maybe]?
// or has atom_sel changed in the mean time?
bool
molecule_class_info_t::moving_atom_matches(mmdb::Atom *at, int this_mol_index_maybe) const {

   bool matches = false;
   if (atom_sel.n_selected_atoms > 0) {
      if (this_mol_index_maybe >= atom_sel.n_selected_atoms) {
         return false;
      } else {
         std::string atom_name_mov = at->name;
         std::string ins_code_mov  = at->GetInsCode();
         std::string alt_conf_mov  = at->altLoc;
         std::string chain_id_mov  = at->GetChainID();
         int resno_mov = at->GetSeqNum();

         std::string atom_name_ref = atom_sel.atom_selection[this_mol_index_maybe]->name;
         std::string ins_code_ref  = atom_sel.atom_selection[this_mol_index_maybe]->GetInsCode();
         std::string alt_conf_ref  = atom_sel.atom_selection[this_mol_index_maybe]->altLoc;
         std::string chain_id_ref  = atom_sel.atom_selection[this_mol_index_maybe]->GetChainID();
         int resno_ref = atom_sel.atom_selection[this_mol_index_maybe]->GetSeqNum();

         if (atom_name_ref == atom_name_mov) {
            if (ins_code_ref == ins_code_mov) {
               if (resno_ref == resno_mov) {
                  if (alt_conf_ref == alt_conf_mov) {
                     if (chain_id_mov == chain_id_ref) { // 20170612 extra condition added, Oliver Clarke bug
                        matches = true;
                     }
                  }
               }
            }
         }
      }
   }
   return matches;
}



// Another attempt at that:
//
// find "by hand" the atom with the given characteristics in
// the atom selection.
//
// return -1 if atom not found.
// Note we have to search for " CA " etc
//
int molecule_class_info_t::atom_index(const char *chain_id, int iresno, const char *atom_id) {

   int n = atom_sel.n_selected_atoms;
   for (int i=0; i<n; i++) {
      if ( ( ! strcmp(atom_id,atom_sel.atom_selection[i]->name) ) &&
           (atom_sel.atom_selection[i]->residue->seqNum == iresno)  &&
           ( ! strcmp(chain_id,atom_sel.atom_selection[i]->residue->GetChainID()) )
           ) {
         return i;
      }
   }

   return -1;
}
int
molecule_class_info_t::atom_index_first_atom_in_residue(const std::string &chain_id,
                                                          int iresno,
                                                        const std::string &ins_code) const {

   bool tacf = 0; // test_alt_conf_flag
   return atom_index_first_atom_in_residue_internal(chain_id, iresno, ins_code, "", tacf);
}

int
molecule_class_info_t::atom_index_first_atom_in_residue(const std::string &chain_id,
                                                          int iresno,
                                                        const std::string &ins_code,
                                                        const std::string &alt_conf) const {

   bool tacf = 1; // test_alt_conf_flag
   return atom_index_first_atom_in_residue_internal(chain_id, iresno, ins_code, alt_conf, tacf);

}

int
molecule_class_info_t::atom_index_first_atom_in_residue_internal(const std::string &chain_id,
                                                                 int iresno,
                                                                 const std::string &ins_code,
                                                                 const std::string &alt_conf,
                                                                 bool test_alt_conf_flag) const {

   int index = -1; // failure
   int selHnd = atom_sel.mol->NewSelection();
   int nSelResidues;
   mmdb::PPResidue SelResidues;
   atom_sel.mol->Select(selHnd, mmdb::STYPE_RESIDUE, 1,
                        chain_id.c_str(),
                        iresno, ins_code.c_str(),
                        iresno, ins_code.c_str(),
                        "*",  // residue name
                        "*",  // Residue must contain this atom name?
                        "*",  // Residue must contain this Element?
                        "*",  // altLocs
                        mmdb::SKEY_NEW // selection key
                        );
   atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
   if (nSelResidues > 0) {
      int ic = -1;
      for (int ires=0; ires<nSelResidues; ires++) {
         int natoms;
         mmdb::PPAtom residue_atoms;
         SelResidues[ires]->GetAtomTable(residue_atoms, natoms);
         for (int iatom=0; iatom<natoms; iatom++) {
            if (test_alt_conf_flag == 0
                || alt_conf == residue_atoms[iatom]->altLoc) {
               if (residue_atoms[iatom]->GetUDData(atom_sel.UDDAtomIndexHandle, ic) == mmdb::UDDATA_Ok) {
                  index = ic;
                  break;
               }
            }
         }
         if (index > -1)
            break;
      }
   }
   atom_sel.mol->DeleteSelection(selHnd);
   // std::cout << "DEBUG:: atom_index_first_atom_in_residue returns " << index << std::endl;
   return index;

}

// Put the regularization results back into the molecule:
//
//// Recall that regularized_asc contains an atom_selection_container_t
// with the new coordinates in.  the mol contains all the molecule and
// the atom_selection contains just the moving parts.
//
void
molecule_class_info_t::replace_coords(const atom_selection_container_t &asc,
                                      bool change_altconf_occs_flag,
                                      bool replace_coords_with_zero_occ_flag) {


   int n_atom = 0;
   int tmp_index;
   bool debug = false;

   if (!asc.mol) {
      std::cout << "ERROR:: unset moving_atoms_asc mol in replace_coords()" << std::endl;
      return;
   }

   make_backup("replace_coords()");

   // debug::
   if (debug) {
      std::cout << "DEBUG:: --------------- replace_coords replacing "
                << asc.n_selected_atoms << " atoms " << std::endl;
      for (int i=0; i<asc.n_selected_atoms; i++) {
         mmdb::Atom *atom = asc.atom_selection[i];
         bool is_ter_state = atom->isTer();
         std::cout << "DEBUG:: in replace_coords, intermediate atom: chain-id :"
                   << atom->residue->GetChainID() <<  ": "
                   << atom->residue->seqNum << " inscode \""
                   << atom->GetInsCode() << "\" name \""
                   << atom->name << "\" altloc \""
                   << atom->altLoc << "\" occupancy: "
                   << atom->occupancy << " :"
                   << " TER state: " << is_ter_state << std::endl;
      }
   }

   // For each atom in the new set of atoms:
   //
   for (int i=0; i<asc.n_selected_atoms; i++) {
      int idx = -1;
      mmdb::Atom *atom = asc.atom_selection[i];
      if (! atom->isTer()) {

         if (debug) { // debug
            std::cout << "considering replacement for selected atom " << coot::atom_spec_t(atom) << std::endl;

            //
            // idx = atom_spec_to_atom_index(std::string(atom->residue->GetChainID()),
            // atom->residue->seqNum, std::string(atom->name));

         }

         if (asc.UDDOldAtomIndexHandle >= 0) { // OK for fast atom indexing

            if (debug)
               std::cout << "... OK for fast atom indexing, asc.UDDOldAtomIndexHandle: "
                         << asc.UDDOldAtomIndexHandle
                         << " for atom " << coot::atom_spec_t(atom) << std::endl;

            if (atom->GetUDData(asc.UDDOldAtomIndexHandle, tmp_index) == mmdb::UDDATA_Ok) {

               if (debug)
                  std::cout << "OK, good GetUDData() for atom " << coot::atom_spec_t(atom) << std::endl;
               if (tmp_index >= 0) {
                  if (moving_atom_matches(atom, tmp_index)) {
                     // std::cout << "      DEBUG:: successfully found old atom index" << std::endl;
                     idx = tmp_index;
                  } else {
                     // std::cout << "DEBUG:: atom index mismatch" << std::endl;
                     idx = full_atom_spec_to_atom_index(std::string(atom->residue->GetChainID()),
                                                        atom->residue->seqNum,
                                                        std::string(atom->GetInsCode()),
                                                        std::string(atom->name),
                                                        std::string(atom->altLoc));
                     // std::cout << "DEBUG:: full_atom_spec_to_atom_index gives index: " << idx << std::endl;
                  }
               } else {
                  // This shouldn't happen.
                  std::cout << "Good Handle, bad index found for old atom: specing" << std::endl;
                  idx = full_atom_spec_to_atom_index(std::string(atom->residue->GetChainID()),
                                                     atom->residue->seqNum,
                                                     std::string(atom->GetInsCode()),
                                                     std::string(atom->name),
                                                     std::string(atom->altLoc));
               }
            } else {

               std::cout << "ERROR:: non-bad handle (" << asc.UDDOldAtomIndexHandle
                         <<  "), but bad GetUDData() for atom " << coot::atom_spec_t(atom) << std::endl;

            }
         } else {

            if (debug)
               std::cout << "DEBUG:: asc.UDDOldAtomIndexHandle is "
                         << asc.UDDOldAtomIndexHandle << " using full atom spec to atom index..."
                         << std::endl;

            idx = full_atom_spec_to_atom_index(std::string(atom->residue->GetChainID()),
                                               atom->residue->seqNum,
                                               std::string(atom->GetInsCode()),
                                               std::string(atom->name),
                                               std::string(atom->altLoc));
            if (idx == -1) {
               std::cout << "DEBUG:: idx: " << idx << "\n";
               std::cout << "ERROR:: failed to find atom in molecule: chain-id :"
                         << std::string(atom->residue->GetChainID()) <<  ": res_no "
                         << atom->residue->seqNum << " inscode :"
                         << std::string(atom->GetInsCode()) << ": name :"
                         << std::string(atom->name) << ": altloc :"
                         << std::string(atom->altLoc) << ":" << std::endl;
            }
         }

         if (change_altconf_occs_flag) {
            if (idx >= 0) {
               n_atom++;
               mmdb::Atom *mol_atom = atom_sel.atom_selection[idx];
               float atom_occ = atom->occupancy;
               // if this is a shelx molecule, then we don't change
               // occupancies this way.  We do it by changing the FVAR
               if (is_from_shelx_ins_flag) {
                  atom_occ = mol_atom->occupancy;

                  // OK, one more go.  We have an occupancy of 31 or -31
                  // say.  Now, the alt conf atoms has been immmediately
                  // added with the old occupancy for the actual FVAR number
                  // - this happens before we get to twiddle the occupancy
                  // slider.  So here we have to find out the index of the
                  // replaced atom and set it's fvar to whatever the slider
                  // value had been set to.

                  int fvar_number = coot::ShelxIns::shelx_occ_to_fvar(atom_occ);
                  if (fvar_number > 1) {
                     //                std::cout << "DEBUG:: replace_coords: setting fvar number "
                     //                          <<  fvar_number << " (generated from occ " << atom_occ << ") to "
                     //                          << graphics_info_t::add_alt_conf_new_atoms_occupancy << std::endl;
                     shelxins.set_fvar(fvar_number, graphics_info_t::add_alt_conf_new_atoms_occupancy);
                  }

                  if (movable_atom(mol_atom, replace_coords_with_zero_occ_flag))
                     mol_atom->SetCoordinates(atom->x,
                                              atom->y,
                                              atom->z,
                                              atom_occ,
                                              mol_atom->tempFactor);
               } else {
                  if (movable_atom(mol_atom, replace_coords_with_zero_occ_flag))
                     mol_atom->SetCoordinates(atom->x,
                                              atom->y,
                                              atom->z,
                                              atom_occ,
                                              mol_atom->tempFactor);
               }

               // similarly we adjust occupancy if this is not a shelx molecule
               if (! is_from_shelx_ins_flag) {
                  adjust_occupancy_other_residue_atoms(mol_atom, mol_atom->residue, 0);
               }
               // std::cout << atom << " coords replace " << idx << " " << mol_atom << std::endl;
            } else {
               std::cout << "ERROR:: bad atom index in replace_coords replacing atom: "
                         << atom << std::endl;
            }
         } else {

            // "don't change alt confs" mode

            if (idx != -1 ) {
               mmdb::Atom *mol_atom = atom_sel.atom_selection[idx];
               if (movable_atom(mol_atom, replace_coords_with_zero_occ_flag)) {
                  if (debug) {
                     coot::Cartesian old_pos(mol_atom->x, mol_atom->y, mol_atom->z);
                     coot::Cartesian new_pos(atom->x, atom->y, atom->z);
                     double d = (new_pos - old_pos).amplitude();
                     std::cout << "    changing coords for atom with idx " << idx << " "
                               << coot::atom_spec_t(mol_atom) << std::endl;
                     std::cout << "   " << old_pos << " " << new_pos << " moved-by " << d << std::endl;
                  }
                  mol_atom->SetCoordinates(atom->x,
                                           atom->y,
                                           atom->z,
                                           mol_atom->occupancy,
                                           mol_atom->tempFactor);
                  n_atom++;
               }
            } else {
               std::cout << "WARNING:: bad atom idx -1" << std::endl;
            }
         }
      }
   }
   std::cout << "INFO:: replace_coords: " << n_atom << " atoms updated." << std::endl;
   have_unsaved_changes_flag = 1;

   if (show_symmetry) {  // internal
      update_symmetry();
   }

   make_bonds_type_checked(__FUNCTION__);

}

// helper function for above function
bool
molecule_class_info_t::movable_atom(mmdb::Atom *mol_atom, bool replace_coords_with_zero_occ_flag) const {

   bool m = 1;

   if ((mol_atom->occupancy < 0.0001) &&
       (mol_atom->occupancy > -0.0001))
      if (replace_coords_with_zero_occ_flag == 0)
         m = 0; // zero occupancy and "dont move zero occ atoms is set"
   return m;
}


// This relies on the mol of the asc being different to the mol of the
// atom_sel.
//
// If it is the same, then error and do nothing.
//
// Typically, this is called by fit terminal residue, which has its
// mol created from a pcmmdbmanager() from the molecule of the
// residue, so this is fine in this case.
//
void
molecule_class_info_t::insert_coords(const atom_selection_container_t &asc) {

   // for each residue in the asc, do a InsResidue into its chain:

   if (! (atom_sel.n_selected_atoms > 0) ) {
      std::cout << "ERROR: Can't insert_coords this asc  - no atoms in molecule!\n";
   } else {

      // pointer comparison
      if (asc.mol == atom_sel.mol) {
         std::cout << "ERROR:: matching asc.mol and atom_sel.mol in insert_coords\n";
         std::cout << "ERROR:: new algorithm required\n";
      } else {
         make_backup("insert_coords()"); // checks backup_this_molecule
         insert_coords_internal(asc);
      }
   }
}

// This does and InsResidue or an AddResidue() as appropriate.
//
void
molecule_class_info_t::insert_coords_internal(const atom_selection_container_t &asc) {

   // run over each chain, residue of the asc (if terminal residue
   // fit only one chain, one residue, of course).

   short int inserted = 0; // not inserted yet
   mmdb::Chain *asc_chain;
   int imod = 1;
   mmdb::Model *asc_model_p = asc.mol->GetModel(imod);
   int asc_n_chains = asc_model_p->GetNumberOfChains();
   for (int i_asc_chain=0; i_asc_chain<asc_n_chains; i_asc_chain++) {
      asc_chain = asc.mol->GetChain(1,i_asc_chain);
      int nres_asc = asc_chain->GetNumberOfResidues();
//       std::cout << "DEBUG:: There are " << nres_asc << " residues in "
//                 << "asc_chain (chain id: " << asc_chain->GetChainID()
//                 << ")." << std::endl;

      int udd_atom_index = asc.UDDAtomIndexHandle;

      for (int ires_asc=0; ires_asc<nres_asc; ires_asc++) {
         mmdb::Residue *asc_residue = asc_chain->GetResidue(ires_asc);

         // Now find the corresponding chain in our atom_sel.mol:

         mmdb::Chain *chain;
         int imodel = 1;
         int n_chains = atom_sel.mol->GetNumberOfChains(imodel);
         for (int i_chain=0; i_chain<n_chains; i_chain++) {

            chain = atom_sel.mol->GetChain(1,i_chain);

            // test chains
            std::string asc_chain_str(asc_chain->GetChainID());
            std::string mol_chain_str(    chain->GetChainID());
//             std::cout << "comparig chain ids :" << asc_chain_str << ": :"
//                       << mol_chain_str << ":" << std::endl;
            if (asc_chain_str == mol_chain_str) {

               // insert that residue!
               bool embed_in_chain_flag = false;
               // use the coot-utils function, not this one
               mmdb::PResidue res = coot::deep_copy_this_residue_old_style(asc_residue, "", 1, udd_atom_index, embed_in_chain_flag);
//                std::cout  << "DEBUG:: inserting residue in chain "
//                           << mol_chain << " residue number "
//                           << asc_residue->GetSeqNum()
//                           << " i_chain = " << i_chain
//                           << " ires_asc = " << ires_asc
//                           << std::endl;

               std::pair<int, mmdb::Residue *> serial_number =
                  find_serial_number_for_insert(asc_residue->GetSeqNum(),
                                                asc_residue->GetInsCode(),
                                                mol_chain_str);

//                std::cout << "DEBUG:: returned serial_number: " << serial_number
//                          << std::endl;

               if (res) {

                  if (serial_number.first != -1) {
                     // insert at this position (other residues are
                     // shifted up).
                     chain->InsResidue(res, serial_number.first);
                     coot::copy_segid(serial_number.second, res);
                     inserted = 1;
                  } else {

                     // std::cout << "DEBUG:: insert_coords_internal() add residue\n";
                     mmdb::Residue *last_residue = last_residue_in_chain(chain);
                     if (last_residue) {

                        if (0) {   // debug::
                           int nat = last_residue->GetNumberOfAtoms();
                           for (int iat=0; iat<nat; iat++) {
                              std::cout << iat << " of " << nat << " "
                                        << last_residue->GetAtom(iat) << std::endl;
                           }
                        }

                        chain->AddResidue(res);
                        coot::copy_segid(last_residue, res);
                        inserted = 1;
                     }
                  }
               }
            }
            //if (inserted) break;
         }


         if (! inserted) {
            // OK, there was no chain in the current mol that matches
            // the chain of the asc.
            // Let's copy the asc chain and add it to atom_sel.mol
            mmdb::Chain *new_chain = new mmdb::Chain;
            int imodel = 1;
            mmdb::Model *this_model = atom_sel.mol->GetModel(imodel);
            this_model->AddChain(new_chain);
            new_chain->SetChainID(asc_chain->GetChainID());

            std::cout << "DEBUG:: Creating a new chain " << asc_chain->GetChainID()
                      << std::endl;
            bool embed_in_chain_flag = false;
            // use the coot-utils function, not this one
            mmdb::Residue *res = coot::deep_copy_this_residue_old_style(asc_residue, "", 1, udd_atom_index, embed_in_chain_flag);
            if (res) {
               new_chain->AddResidue(res);
               atom_sel.mol->FinishStructEdit(); // so that we don't keep adding a
                                                 // new Chain to atom_sel.mol
            }

         }
         //if (inserted) break;
      }
      //if (inserted) break;
   }
   atom_sel.mol->FinishStructEdit();
   update_molecule_after_additions();
}



void
molecule_class_info_t::insert_coords_change_altconf(const atom_selection_container_t &asc) {

   // There are 2 things we want to do here.
   //
   // For matching atoms, if there is only one atom that matches the
   // spec (appart from the altconf), we move the original atoms
   // altconf to "A".  If they are not "" we leave don't change the
   // altconf.
   // (see table in // molecule_class_info_t::make_new_alt_conf()
   // [molecule-class-info-other.cc]
   //

   // The second thing is the change the occ of the existing atoms:
   // 1 atom:  new_occ_existing_atom = 1 - occ_new_atom;
   // general: new_occ_existing_atom = current_occ_existing_atom - occ_new_atom/n_similar_atoms
   // where n_similar_atoms is the number of alt confs we have for that atom.
//    std::cout << "DEBUG:: ----------------- in insert_coords_change_altconf ------ " << std::endl;
//    std::cout << "DEBUG:: IN insert_coords_change_altconf" << std::endl;
   make_backup("insert_coords_change_altconf");

   // OK if we were from a shelx ins file, then we have to create a
   // new FVAR for this new alt conf thingy.
   int shelx_occ_fvar_number = -1;
   if (is_from_shelx_ins_flag) {
      // OK, what was the occupancy?
      if (asc.n_selected_atoms > 0) {
         float occ = asc.atom_selection[0]->occupancy;
//          std::cout << "DEBUG:: IN insert_coords_change_altconf adding fvar "
//                    << occ << std::endl;
         shelx_occ_fvar_number = 10 * shelxins.add_fvar(occ); // FVAR 1 is not written
         // to SHELX file so if shelx ins file has 1 FVAR value, then we've just
         // created shelx FVAR 3.
         shelx_occ_fvar_number += 1;  // so thats (e.g.) 1 x the 20th FVAR
      }
   }

   char *chain_id;
   char *atom_name;
   int  resno;
   float occ;
   mmdb::Atom *at;
   for(int i=0; i<asc.n_selected_atoms; i++) {
      at = asc.atom_selection[i];
      chain_id  = at->GetChainID();
      atom_name = at->GetAtomName();
      resno     = at->GetSeqNum();
      occ       = at->occupancy;
      char *inscode = at->GetInsCode();

      // Now find that atom the corresponding atom (with altconf "" in
      // the original atoms).  We skip over atoms that don't have
      // altconf "").

      int selHnd = atom_sel.mol->NewSelection();
      atom_sel.mol->SelectAtoms(selHnd, 0, chain_id,
                                resno, inscode,
                                resno, inscode,
                                "*", // residue name
                                atom_name,
                                "*",
                                "*"); // alt-loc
      int nSelAtoms;
      mmdb::PPAtom local_SelAtom;
      atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);

      if (nSelAtoms == 0) {
         // debugging
//                std::cout << "add alt conf: skipping " << atom_name << "/"
//                 << resno << "/" << chain_id << " in this molecule: ("
//                 <<  imol_no << ") " << name_ << std::endl;
      } else {

         // Let's deal with every atom in this residue that that has
         // the same atom name.

         // first check for atoms with alt conf "" and move them over if needed
         for (int iat=0; iat<nSelAtoms; iat++) {
            std::string current_alt_conf = local_SelAtom[iat]->altLoc;
            if (current_alt_conf == "") {
               std::string new_alt_conf("A");
               // force it down the atom's throat :)
               strcpy(local_SelAtom[0]->altLoc, new_alt_conf.c_str()); // 20220620-PE changed from strncpy (new_alt_conf length is 15
                                                                       //             says the compiler)
            }
         }

         if (shelx_occ_fvar_number == -1) {
            // i.e. was not from a shelx ins file (normal case):
            //
            // now stuff occupancies in...
            for (int iat=0; iat<nSelAtoms; iat++) {
               // local_SelAtom[0]->occupancy = 1.0 - occ; // complemetary (1+1 atom case)
               local_SelAtom[iat]->occupancy -= occ/float(nSelAtoms); // complemetary (general case)

               // 20091227 But don't add atoms with negative
               // occupancy, e.g. the residue before the split was at
               // zero occupancy.
               if (local_SelAtom[iat]->occupancy < 0.0)
                  local_SelAtom[iat]->occupancy = 0.0;
            }

         } else {

            // This is the SHELX case:
            if (nSelAtoms > 1) {
               // so we added an alt conf to a residue that already
               // has an alt conf.  This involves messing with SUMP.
               // Let's not handle that for now.
               std::cout << "WARNING:: SHELX occupancy handler under-resourced on handling "
                         << at << std::endl;
            }  else {
               local_SelAtom[0]->occupancy = -shelx_occ_fvar_number;
            }
         }
      }
      atom_sel.mol->DeleteSelection(selHnd);
   }
   insert_coords_atoms_into_residue_internal(asc, shelx_occ_fvar_number);

}

// In this instance, we don't want to install a whole residue, we want
// to install atoms in this residue (alt conf B) into a a atom_sel mol
// residue that contains (say) "" and "A".
//
// -1 is passed as shelx_occ_fvar_number if this atom was not from a
// SHELX ins file.  If shelx_occ_fvar_number > 1, then use this as the
// new atom's occupancy.
//
void
molecule_class_info_t::insert_coords_atoms_into_residue_internal(const atom_selection_container_t &asc,
                                                                 int shelx_occ_fvar_number) {

   char *chain_id;
   char *atom_name;
   int  resno;
   mmdb::Atom *at;
   int afix_handle_this_mol = -1;
   int afix_handle_intermediate_mol = -1;

   afix_handle_this_mol    = atom_sel.mol->GetUDDHandle(mmdb::UDR_ATOM, "shelx afix");
   afix_handle_intermediate_mol = asc.mol->GetUDDHandle(mmdb::UDR_ATOM, "shelx afix");

//    std::cout << "DEBUG in insert_coords_atoms_into_residue_internal afix handles:"
//              << afix_handle_this_mol << " " << afix_handle_intermediate_mol << std::endl;

   for(int i=0; i<asc.n_selected_atoms; i++) {
      at = asc.atom_selection[i];
      chain_id  = at->GetChainID();
      atom_name = at->GetAtomName();
      resno     = at->GetSeqNum();

      // Now find the corresponding residue in atom_sel.mol;

      int selHnd = atom_sel.mol->NewSelection();
      int nSelResidues;
      mmdb::PPResidue SelResidues;
      atom_sel.mol->Select(selHnd, mmdb::STYPE_RESIDUE, 1,
                           chain_id,
                           resno, "*",
                           resno, "*",
                           "*",  // residue name
                           "*",  // Residue must contain this atom name?
                           "*",  // Residue must contain this Element?
                           "*",  // altLocs
                           mmdb::SKEY_NEW // selection key
                           );
      atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

      if (nSelResidues != 1) {
         std::cout << "ERROR:: something broken in residue selection in ";
         std::cout << "insert_coords_atoms_into_residue_internal: got " << nSelResidues
                   << " residues." << std::endl;
      } else {
         mmdb::Atom *t = new mmdb::Atom;
         t->Copy(at);
         SelResidues[0]->AddAtom(t);
         // if these coords were from a shelx ins file, then we are
         // passed the free varible number (FVAR) for this atom's
         // occupancy.
         if (shelx_occ_fvar_number > 1)
            t->occupancy = shelx_occ_fvar_number;

         // If these coords were from a shelx ins file, then we need
         // to copy the AFIX numbers too:
         int afix_number; // set by getUDD
         if (afix_handle_intermediate_mol > -1) {
            int ierr = at->GetUDData(afix_handle_intermediate_mol, afix_number);
            if (ierr == mmdb::UDDATA_Ok) {
               if (afix_handle_this_mol > -1) {
                  t->PutUDData(afix_handle_this_mol, afix_number);
               } else {
                  std::cout << "ERROR:: bad afix handle for this molecule in "
                            << "insert_coords_atoms_into_residue_internal"
                            << afix_handle_this_mol << " " << at << std::endl;
               }

             } else {
               if (is_from_shelx_ins_flag)
                  std::cout << "ERROR:: attempt to get UDD afix number from "
                            << "intermediate molecule failed " << at << std::endl;
             }
         } else {
            std::cout << "ERROR:: bad afix handle for intermediate molecule in "
                      << "insert_coords_atoms_into_residue_internal"
                      << afix_handle_intermediate_mol << " " << at << std::endl;
         }
      }
      atom_sel.mol->DeleteSelection(selHnd);
   }
   atom_sel.mol->FinishStructEdit();
   update_molecule_after_additions();
}

// We need to find the serial number of the residue after the residue
// we want to insert (i.e. the new residue will be inserted just
// before the residue whose serial number we return).
//
// return -1 on error.
std::pair<int, mmdb::Residue *>
molecule_class_info_t::find_serial_number_for_insert(int seqnum_for_new,
                                                     const std::string &ins_code_for_new,
                                                     const std::string &chain_id) const {

   int iserial_no = -1;
   std::pair<int, std::string> current_diff(999999, "");
   int n_chains = atom_sel.mol->GetNumberOfChains(1);
   mmdb::Residue *res = NULL;

   for (int i_chain=0; i_chain<n_chains; i_chain++) {

      mmdb::Chain *chain_p = atom_sel.mol->GetChain(1,i_chain);

      if (chain_p) {

         std::string mol_chain(chain_p->GetChainID());

         if (chain_id == mol_chain) {

            // Find the first residue that has either the residue number or insertion code
            // greater than the passed parameters

            int nres = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<nres; ires++) { // ires is a serial number
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);

               // we need to consider insertion codes here

               int diff = residue_p->GetSeqNum() - seqnum_for_new;

               if (diff > 0) {
                  res = residue_p;
                  iserial_no = ires;
                  break;
               } else {
                  if (diff == 0) {
                     std::string ins_code_this = residue_p->GetInsCode();
                     if (ins_code_this > ins_code_for_new) {
                        res = residue_p;
                        iserial_no = ires;
                        break;
                     }
                  }
               }
            }
         }
      }
   }
   return std::pair<int, mmdb::Residue *> (iserial_no, res);
}


void
molecule_class_info_t::add_terminal_residue_wrapper(const coot::residue_spec_t &res_spec,
                                                     const std::string &residue_type) {

   mmdb::Residue *res_p = get_residue(res_spec);
   if (! res_p) return;

   if (!coot::util::is_nucleotide_by_dict_dynamic_add(res_p, graphics_info_t::Geom_p())) {
/*             g.execute_add_terminal_residue(naii.imol,
                                           term_type,
                                           res_p,
                                           chain_id,
                                           g.add_terminal_residue_type,  // eg. "ALA" or "UNK"
                                           add_it_now_flag); */
   } else {
      // g.execute_simple_nucleotide_addition(naii.imol, term_type, res_p, chain_id);
   }
}

int
molecule_class_info_t::add_terminal_residue_using_phi_psi(const std::string &chain_id,
                                                          int res_no,
                                                          const std::string &residue_type,
                                                          float phi, float psi) {

   int status = 0;
   mmdb::Residue *res = get_residue(chain_id, res_no, "");
   if (! res) {
      std::cout << "WARNING:: add_terminal_residue_using_phi_psi() residue not found for \""
                << chain_id << "\" " << res_no << std::endl;
   } else {
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      res->GetAtomTable(residue_atoms, n_residue_atoms);
      if (n_residue_atoms) {
         mmdb::Atom *at = residue_atoms[0];
         int atom_indx = get_atom_index(at);
         if (atom_indx < 0) {
            // This should not happen.
            std::cout << "WARNING:: add_terminal_residue_using_phi_psi() "
                      << "Failed to get atom index for 0th atom in \""
                      << chain_id << "\" " << res_no << std::endl;
         } else {
            std::string term_type = get_term_type(atom_indx);
            int found_atoms_count = 0;
            clipper::Coord_orth previous_ca, previous_c, previous_n;
            mmdb::Atom *C_terminal_this_residue_O = NULL;
            for (int iat=0; iat<n_residue_atoms; iat++) {
               std::string atom_name = residue_atoms[iat]->name;

               // PDBv3 FIXME.
               if (atom_name == " CA ") {
                  found_atoms_count += 1;
                  previous_ca = clipper::Coord_orth(residue_atoms[iat]->x,
                                                    residue_atoms[iat]->y,
                                                    residue_atoms[iat]->z);
               }
               if (atom_name == " C  ") {
                  found_atoms_count += 2;
                  previous_c = clipper::Coord_orth(residue_atoms[iat]->x,
                                                   residue_atoms[iat]->y,
                                                   residue_atoms[iat]->z);
               }
               if (atom_name == " N  ") {
                  found_atoms_count += 4;
                  previous_n = clipper::Coord_orth(residue_atoms[iat]->x,
                                                   residue_atoms[iat]->y,
                                                   residue_atoms[iat]->z);
               }
               if (atom_name == " O  ") {
                  // we store this because, if a residue is added to
                  // the C-terminus, we want to move the O of the
                  // current residue to be in the peptide plane that
                  // is made when we have the new residue atoms.
                  C_terminal_this_residue_O = residue_atoms[iat];
               }
            }
            coot::minimol::residue r;
            if (term_type == "N") {
               if (! (found_atoms_count&7)) {
                  std::cout << "Bad for N current atom selection " << std::endl;
               } else {
                  // happy path
                  r = coot::build_N_terminal_ALA(phi, psi, res_no-1,
                                                 previous_n,
                                                 previous_ca,
                                                 previous_c, 30);
                  std::pair<bool, clipper::Coord_orth> cb = coot::cbeta_position(r);
                  if (cb.first) {
                     coot::minimol::atom at(" CB ", " C", cb.second, "", 1.0, 30);
                     r.addatom(at);
                  }
               }
            }

            // treat singletons like a C-terminal addition.
            if (term_type == "C" || term_type == "singleton" ) {
               if (! (found_atoms_count&7)) {
                  std::cout << "Bad for N current atom selection " << std::endl;
               } else {
                  r = coot::build_C_terminal_ALA(phi, psi, res_no+1,
                                                 previous_n,
                                                 previous_ca,
                                                 previous_c,
                                                 30);
                  std::pair<bool, clipper::Coord_orth> cb = coot::cbeta_position(r);
                  if (cb.first) {
                     coot::minimol::atom at(" CB ", " C", cb.second, "", 1.0, 30);
                     r.addatom(at);

                     // Also, we need to move the O of the current
                     // atom into the plane of the CA
                     //
                     if (C_terminal_this_residue_O) {
                        std::pair<bool, clipper::Coord_orth> ca_new(false, clipper::Coord_orth(0,0,0));
                        std::pair<bool, clipper::Coord_orth>  n_new(false, clipper::Coord_orth(0,0,0));

                        // PDBv3 FIXME.
                        // now get ca_new and n_new from r:
                        for (unsigned int jat=0; jat<r.n_atoms(); jat++) {
                           if (r[jat].name == " CA ") {
                              ca_new.first = true;
                              ca_new.second = r[jat].pos;
                           }
                           if (r[jat].name == " N  ") {
                              n_new.first = true;
                              n_new.second = r[jat].pos;
                           }
                        }
                        // If we found them (we should have), generate
                        // a new pos for this O and apply it.
                        if (ca_new.first && n_new.first) {
                           double torsion_ca_n_c_o = clipper::Util::d2rad(0.0);
                           double angle_n_c_o      = clipper::Util::d2rad(132.02);
                           // bond/attach to the previous_c.
                           clipper::Coord_orth new_pos(ca_new.second,
                                                       n_new.second,
                                                       previous_c,
                                                       1.23, angle_n_c_o, torsion_ca_n_c_o);
                           C_terminal_this_residue_O->x = new_pos.x();
                           C_terminal_this_residue_O->y = new_pos.y();
                           C_terminal_this_residue_O->z = new_pos.z();
                        }
                     }
                  }
               }
            }
            if (r.atoms.size()) {
               coot::minimol::fragment f(chain_id);
               f.addresidue(r,0);
               coot::minimol::molecule m(f);
               mmdb::Manager *mol_new = m.pcmmdbmanager();

               atom_selection_container_t asc = make_asc(mol_new);
               insert_coords_internal(asc);

               update_molecule_after_additions(); // create UDDAtomIndexHandle for new atoms.
               status = 1;
            } else {
               std::cout << "No residue added for term type "
                         << term_type << std::endl;
            }
         }
      }
   }
   return status;
}


// When a new residue is added to the C-terminus of a chain/fragment, we will need to move
// the O of this one to make a proper peptide plane (the position of the next residue
// was not dependent on the position of the O of this one).
//
// find the residue in the chain that is after this one (res_p).  That will provide the N, CA
// positions that will allow us to position the O of res_p.
void
molecule_class_info_t::move_O_atom_of_added_to_residue(mmdb::Residue *res_p, const std::string &chain_id) {

   bool moved = false;
   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         std::string chain_id_this(chain_p->GetChainID());
         if (chain_id_this == chain_id) {
            for (int ires=0; ires<nres; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p == res_p) {
                  int ires_next = ires + 1;
                  if (ires_next < nres) {
                     mmdb::Residue *res_next_p = chain_p->GetResidue(ires_next);
                     if (res_next_p) {
                        mmdb::Atom *ca_this = residue_p->GetAtom(" CA ");
                        mmdb::Atom *c_this  = residue_p->GetAtom(" C  ");
                        mmdb::Atom *o_this  = residue_p->GetAtom(" O  ");
                        mmdb::Atom *ca_next = res_next_p->GetAtom(" CA ");
                        mmdb::Atom *n_next  = res_next_p->GetAtom(" N  ");

                        if (ca_this && c_this && o_this && ca_next && n_next) {

                           clipper::Coord_orth ca_this_pos = coot::co(ca_this);
                           clipper::Coord_orth c_this_pos  = coot::co(c_this);
                           clipper::Coord_orth ca_next_pos = coot::co(ca_next);
                           clipper::Coord_orth  n_next_pos = coot::co(n_next);
                           double angle   = clipper::Util::d2rad(123.0); // N-C-O
                           double tors_deg = 0.0; // O is trans to the CA of the next residue
                           // unless peptide is cis.
                           double tors_peptide = clipper::Coord_orth::torsion(ca_this_pos, c_this_pos,
                                                                              n_next_pos, ca_next_pos);
                           // cis or trans
                           if (std::abs(tors_peptide) < M_PI_2) tors_deg = 180.0;
                           double torsion = clipper::Util::d2rad(tors_deg);
                           clipper::Coord_orth new_o_pos_for_current_res_new(ca_next_pos, n_next_pos,
                                                                             c_this_pos, 1.231, angle, torsion);

                           o_this->x = new_o_pos_for_current_res_new.x();
                           o_this->y = new_o_pos_for_current_res_new.y();
                           o_this->z = new_o_pos_for_current_res_new.z();

                           moved = true;
                           make_backup("move_O_atom_of_added_to_residue()");
                        } else {
                           std::cout << "WARNING:: missing atoms in move_O_atom_of_added_to_residue " << std::endl;
                        }
                     }
                  }
                  break;
               }
            }
         }
      }
   }
   if (moved) {
      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked(__FUNCTION__);
   }

}


// Put the regularization results back into the molecule:
//
// Recall that regularized_asc contains an atom_selection_container_t
// with the new coordinates in.  the mol contains all the molecule and
// the atom_selection contains just the moving parts.
//
void
molecule_class_info_t::add_coords(const atom_selection_container_t &asc) {

   mmdb::Atom *atom;
   mmdb::Atom *mol_atom; // an atom already existing in mol
   int n_atom = 0;
   mmdb::Chain *chain;

   // std::cout << "DEBUG:: ----------------- in add_coords ----------- " << std::endl;

   make_backup("add_coords()");

   for (int i=0; i<asc.n_selected_atoms; i++) {
      int idone = 0;
      atom = asc.atom_selection[i];
      // chain = atom->GetChain();

      // run over chains of the existing mol
      int nchains = atom_sel.mol->GetNumberOfChains(1);
      for (int ichain=0; ichain<nchains; ichain++) {

         chain = atom_sel.mol->GetChain(1,ichain);
         std::string atom_chain_id(atom->GetChainID());
         std::string mol_chain_id(chain->GetChainID());

         if (atom_chain_id == mol_chain_id) {

            // int iseqno_at = atom->GetSeqNum();
            int nres = chain->GetNumberOfResidues();
            for (int ires=0; ires<nres; ires++) {
               mmdb::PResidue res = chain->GetResidue(ires);
               if (res) { // fixes bug 030813, but best way?
                          // 030814, ah, we discover FinishStructEdit().
                  if (res->GetSeqNum() == atom->GetSeqNum()) {
                     int natom = res->GetNumberOfAtoms();
                     for (int iat=0; iat<natom; iat++) {

                        mol_atom = res->GetAtom(atom->GetAtomName());
                        if (mol_atom) { // we got a match
                           // replace the coordinates then

                           // This should not happen very often
                           std::cout << "add_coords: replacing " << mol_atom
                                     << " with new atom " << atom << std::endl;
                           mol_atom->SetCoordinates(atom->x,
                                                    atom->y,
                                                    atom->z,
                                                    mol_atom->occupancy,
                                                    mol_atom->tempFactor);
                           idone = 1;
                           break;

                        } else {

                           std::cout << "adding atom to existing residue "
                                     << atom << " (already has "
                                     << res->GetNumberOfAtoms() << " atoms)"
                                     << std::endl;
                           mmdb::Atom *new_atom = new mmdb::Atom;
                           new_atom->Copy(atom);
                           res->AddAtom(new_atom);
                           new_atom->occupancy = 1.0;
                           new_atom->tempFactor = 10.0;
                           // chain id:
                           if (0)
                              std::cout << "setting chainid of this new atom from :"
                                        << new_atom->GetChainID() << ": to :"
                                        << atom->GetChainID() << ":" << std::endl;
                           new_atom->residue->chain->SetChainID(atom->GetChainID());
                           idone = 1;
                           n_atom++;
                           break;
                        }
                     }
                  }
               }
               if (idone == 1) break;
            } // residue loop
         }
      }

      if (idone == 0) {

         std::cout << "adding whole residue triggered by atom "
                   << atom << " ";
         std::cout << " with element " << atom->element << std::endl;

         // in this bit of code, atom is an atom from the asc and
         // atom_p is a new atom that we are adding to a new residue
         // (that we are adding to an existing chain (that we do a
         // lookup to find)).
         //
         mmdb::Residue *res_p = new mmdb::Residue;
         mmdb::Atom *atom_p = new mmdb::Atom;
         // mmdb::Chain *chain_p = atom_sel.mol->GetChain(1,0);
         mmdb::PChain chain_p = atom_sel.mol->GetChain(1,atom->GetChainID());
         chain_p->AddResidue(res_p);
         atom_p->SetAtomName(atom->name);
         atom_p->SetCoordinates(atom->x, atom->y, atom->z,
                                atom->occupancy, atom->tempFactor);
         atom_p->SetElementName(atom->element);
         res_p->AddAtom(atom_p);
         res_p->seqNum = atom->GetSeqNum();
         res_p->SetResID(atom->residue->name,
                         atom->GetSeqNum(),
                         atom->GetInsCode());

         // add to end:
         atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
         atom_sel.mol->FinishStructEdit();
      }
   }

   atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
   atom_sel.mol->FinishStructEdit();
   std::cout << "INFO:: " << n_atom << " atoms added to molecule." << std::endl;

   // now regenerate the atom_selection
   //

   // Uncomment when we have fixed the bug, this seems not to be it.
   //
   // clear out the old
   // if (atom_sel.atom_selection != NULL)
   // delete [] atom_sel.atom_selection;

   // and in with the new:
   int selHnd = atom_sel.mol->NewSelection();
   atom_sel.mol->SelectAtoms(selHnd, 0,"*",mmdb::ANY_RES,"*",mmdb::ANY_RES,
                            "*","*",  // EndInsertionCode, RNames
                            "*","*",  // ANames, Elements
                            "*" );    // Alternate locations.

   int old_n_atoms = atom_sel.n_selected_atoms;
   atom_sel.mol->GetSelIndex(selHnd,
                             atom_sel.atom_selection,
                             atom_sel.n_selected_atoms);

   std::cout << "INFO:: old n_atoms: " << old_n_atoms << " new: "
             << atom_sel.n_selected_atoms << std::endl;

   have_unsaved_changes_flag = 1;

   make_bonds_type_checked(__FUNCTION__);
   // std::cout << "DEBUG:: ---------------- done add_coords ----------- " << std::endl;
}

// save to default file name if has unsaved changes.  Return non-zero on problem.
int
molecule_class_info_t::quick_save() {

   if (Have_unsaved_changes_p()) {
      std::string s = stripped_save_name_suggestion();
      save_coordinates(s);
   }
   return 0;
}


void
molecule_class_info_t::close_yourself() {

   // Deletion causing problems on application closure

   bool was_map    = false;
   bool was_xmap   = false;
   bool was_nxmap  = false;
   bool was_coords = false;

   name_ = ""; // not "Baton Atoms" or anything now.

   if (atom_sel.n_selected_atoms > 0)
      was_coords = 1;

   if (has_xmap()) {
      was_xmap = true;
      was_map = true;
   }

   if (has_nxmap()) {
      was_map = true;
      was_nxmap = true;
   }

   bool delete_stored_data = true; // does this stop the crash on charybdis?
   if (delete_stored_data) {
      if (original_fphis_filled) {
	 // do I have a race condition?
	 clipper::HKL_data< clipper::datatypes::F_phi<float> >  *tmp_p = original_fphis_p;
	 original_fphis_p = 0;
         delete tmp_p;
      }

      if (original_fobs_sigfobs_filled) {
         delete original_fobs_sigfobs_p;
	 original_fobs_sigfobs_p = 0;
      }

      if (original_r_free_flags_p) { // no flag for filled?
         delete original_r_free_flags_p;
	 original_r_free_flags_p = 0;
      }
   }

   // delete from display manager combo box
   //
   graphics_info_t g;
   if (g.use_graphics_interface_flag) {
      GtkWidget *display_control_window = widget_from_builder("display_control_window_glade");
      if (display_control_window) {
         if (was_map) {
            GtkWidget *map_vbox = widget_from_builder("display_map_vbox");
            if (GTK_IS_BOX(map_vbox)) {
               GtkWidget *item_widget = gtk_widget_get_first_child(map_vbox);
               while (item_widget) {
                  GtkWidget *next_item = gtk_widget_get_next_sibling(item_widget);
                  int imol_this = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(item_widget), "imol"));
                  if (imol_this == imol_no) {
                     gtk_box_remove(GTK_BOX(map_vbox), item_widget);
                  }
                  item_widget = next_item;
               };
            }
         }

         if (was_coords) {
            GtkWidget *coords_vbox = widget_from_builder("display_molecule_vbox");
            if (GTK_IS_BOX(coords_vbox)) {
               std::cout << "in close_yourself() fix container B foreach" << std::endl;

               GtkWidget *item_widget = gtk_widget_get_first_child(coords_vbox);
               while (item_widget) {
                  GtkWidget *next_item = gtk_widget_get_next_sibling(item_widget);
                  int imol_this = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(item_widget), "imol"));
                  if (imol_this == imol_no) {
                     gtk_box_remove(GTK_BOX(coords_vbox), item_widget);
                  }
                  item_widget = next_item;
               };
            }
         }
      }
      graphics_info_t::refresh_validation_graph_model_list();
   }

   if (was_coords) {
      atom_sel.mol->DeleteSelection(atom_sel.SelectionHandle);
      delete atom_sel.mol;
   }

   if (was_xmap) {
      fc_skeleton_draw_on = 0; // turn off the skeleton

      // delete [] draw_vectors;
      // draw_vectors = NULL;
      clear_draw_vecs(); // not difference map
      clear_diff_map_draw_vecs();


      // delete [] diff_map_draw_vectors;
      // diff_map_draw_vectors = NULL;

      clipper::Xmap<float> empty;
      xmap = empty; // clear xmap
   }

   if (was_nxmap) {
      fc_skeleton_draw_on = 0; // turn off the skeleton

      // delete [] draw_vectors;
      // draw_vectors = NULL;
      clear_draw_vecs();
      clear_diff_map_draw_vecs();

      // delete [] diff_map_draw_vectors;
      // diff_map_draw_vectors = NULL;
      clipper::NXmap<float> empty;
      nxmap = empty; // clear nxmap
   }

   bonds_box.clear_up();
   // symmetry_bonds_box?  (It is a vector of pairs)
   draw_it = 0;
   draw_it_for_map = 0;
   draw_it_for_map_standard_lines = 0;

   // Do these whatever the molecule type:
   atom_sel.n_selected_atoms = 0;
   atom_sel.atom_selection = NULL;
   atom_sel.mol = NULL;

   //
   // gl widget redraw is done in close_molecule

}


// Return the atom index of the "next" atom
// -1 on failure.
int
molecule_class_info_t::intelligent_next_atom(const std::string &chain_id,
                                             int resno,
                                             const std::string &atom_name,
                                             const std::string &ins_code,
                                             const coot::Cartesian &rc) {

      // This is really a problem of "what is the next residue?", the
      // actual atom is a superficial problem that is handled by
      // intelligent_this_residue_atom().  We simply have to find this
      // residue in the chain, and return the residue after that.
      //
      // If there is no next residue, use the residue at the beginning of
      // the next chain.
      //
      // If there is no next chain, use the residue at the start of the
      // first chain.
      //
      // If this residue can't be found, then go through chain and look
      // for the first residue that has higher residue number than resno.

      int i_atom_index = -1; // failure initially.
      if (atom_sel.n_selected_atoms <= 0 || atom_sel.mol == NULL) {
         std::cout << "ERROR:: trying to move to (next) atom of a closed molecule!\n";
      } else {

         mmdb::Residue *next_residue = NULL;

         // Note: we may not be at this residue.
         //
         coot::residue_spec_t this_residue_spec(chain_id, resno, ins_code);
         mmdb::Residue *this_residue = get_residue(this_residue_spec);

         if (this_residue) {

         if (close_to_residue(this_residue, rc)) {

            // === move on to next one ===

            // Can we do that by residue index?
            //
            int ser_num = this_residue->index;
            if (ser_num != -1) {
               // yes...
               int ser_num_next = ser_num + 1;
               mmdb::Residue *rr = this_residue->chain->GetResidue(ser_num);
               if (this_residue == rr) {
                  // self reference works as it should
                  next_residue = this_residue->chain->GetResidue(ser_num_next);
               } else {
                  coot::residue_spec_t next_residue_spec(chain_id, resno+1, "");
                  mmdb::Residue *residue_p = get_residue(next_residue_spec);
                  if (residue_p)
                     next_residue = residue_p;
               }
            } else {
               // no...
               // this_residue was not properly inserted into the
               // chain for some reason.
               //
               coot::residue_spec_t next_residue_spec(chain_id, resno+1, "");
               mmdb::Residue *residue_p = get_residue(next_residue_spec);
               if (residue_p)
                  next_residue = residue_p;
            }

            if (next_residue) {

               i_atom_index = intelligent_this_residue_atom(next_residue);

            } else {

               // OK, we need to move onto the next chain.... or find
               // the residue after the gap.
               //

               next_residue = next_residue_missing_residue(coot::residue_spec_t(this_residue));
               if (next_residue) {
                  i_atom_index = intelligent_this_residue_atom(next_residue);
               } else {
                  // we are on the last atom then.
                  i_atom_index = 0;

               }
            }

         } else {
            // Go (back) to this one - because we had moved away from it.
            i_atom_index = intelligent_this_residue_atom(this_residue);
         }
      } else {
         // OK, the residue could not be found, it was a deleted atom
         // perhaps.  So find the next residue in the chain that has a
         // higher residue number than resno.

         // the residue is not found for this_residue_spec
         //
         next_residue = next_residue_missing_residue(this_residue_spec);

         if (next_residue)
            i_atom_index = intelligent_this_residue_atom(next_residue);

         }
      }
      return i_atom_index;
   }


// Return the atom index of the "previous" atom
// -1 on failure.
int
molecule_class_info_t::intelligent_previous_atom(const std::string &chain_id,
                                                 int resno,
                                                 const std::string &atom_name,
                                                 const std::string &ins_code,
                                                 const coot::Cartesian &rc) {

   // This is quite similar to intelligent_next_atom() (see comments
   // there).  However, this is a bit more complex, because we keep a
   // hold on the previous residue and only if "this" residue is found
   // do we make the previous residue previous_residue.
   //
   // We don't handle the "last_residue" in molecule case, the analog
   // of first_residue in the intelligent_next_atom() function.
   // (Maybe we should?)

   int i_atom_index = -1;
   if (atom_sel.n_selected_atoms <= 0 || atom_sel.mol == NULL) {
      std::cout << "ERROR:: trying to move to (next) atom of a closed molecule!\n";
   } else {

      mmdb::Residue *prev_residue_candidate = NULL;
      mmdb::Residue *prev_residue = NULL;
      bool found_this_residue = false;
      int imod = 1;
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      mmdb::Chain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
         chain_p = model_p->GetChain(ichain);
         std::string this_chain_id = chain_p->GetChainID();
         if ((chain_id == this_chain_id) || (found_this_residue && !prev_residue)) {
            int nres = chain_p->GetNumberOfResidues();
            mmdb::PResidue residue_p;
            for (int ires=0; ires<nres; ires++) {
               residue_p = chain_p->GetResidue(ires);
               if (residue_p->GetSeqNum() == resno) {
                  if (ins_code == residue_p->GetInsCode()) {
                     if (chain_id == this_chain_id) {
                        found_this_residue = true;
                        if (prev_residue_candidate) {
                           prev_residue = prev_residue_candidate;
                           break;
                        }
                     }
                  }
               }
               prev_residue_candidate = residue_p;
               if (prev_residue)
                  break;
            }
         }
         if (prev_residue)
            break;
      }

      // OK, we can get here by going backward through the chain
      // (Shift Spacing) and have just deleted the atom we are centred
      // on.
      //
      // So, in the chain of the atom we just deleted, find the first
      // residue that is more than the seqnum of the atom we just
      // deleted, and the residue before that (prev_residue_candidate)
      //
      if (! found_this_residue) {
         for (int ichain=0; ichain<nchains; ichain++) {
            chain_p = model_p->GetChain(ichain);
            if (chain_id == chain_p->GetChainID()) {
               int nres = chain_p->GetNumberOfResidues();
               mmdb::PResidue residue_p;
               for (int ires=0; ires<nres; ires++) {
                  residue_p = chain_p->GetResidue(ires);
                  if (residue_p->GetSeqNum() > resno) {
                     if (prev_residue_candidate)
                        prev_residue = prev_residue_candidate;
                     break;
                  }
                  prev_residue_candidate = residue_p;
               }
            }
         }
      }

      // Handle the case where we are on first atom of water chain and
      // want to go back to protein (in that case
      // prev_residue_candidate would not have been set because it
      // never passes the chain id test)
      //
      if (! prev_residue) {
         mmdb::Chain *prev_chain = NULL;
         mmdb::Chain *prev_chain_candidate = NULL;
         int nchains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<nchains; ichain++) {
            chain_p = model_p->GetChain(ichain);
            if (chain_id == chain_p->GetChainID()) {
               if (prev_chain_candidate) {
                  prev_chain = prev_chain_candidate;
                  break;
               }
            }
            prev_chain_candidate = chain_p; // for next loop
         }
         if (prev_chain) {
            // the last residue in prev_chain then:
            int nres = prev_chain->GetNumberOfResidues();
            if (nres > 0)
               prev_residue = prev_chain->GetResidue(nres-1);
         } else {
            // OK the passed (current) residue was the first in the
            // coordinates.  Pick the last residue of the last chain.
            chain_p = model_p->GetChain(nchains-1);
            if (chain_p) {
               int nres = chain_p->GetNumberOfResidues();
               if (nres > 0)
                  prev_residue = chain_p->GetResidue(nres-1);
            }
         }
      }

      if (prev_residue)
         i_atom_index = intelligent_this_residue_atom(prev_residue);
   }

   return i_atom_index;
}


// If there is a CA or C1' in this residue then return the index of that
// atom, if not, then return the index of the first atom in the
// residue.
//
// Return -1 on no atoms in residue.
//
int
molecule_class_info_t::intelligent_this_residue_atom(mmdb::Residue *res_p) const {

   mmdb::PPAtom residue_atoms;
   int nResidueAtoms;
   int ir = -1;

   if (res_p) {
      res_p->GetAtomTable(residue_atoms, nResidueAtoms);
      for (int i=0; i<nResidueAtoms; i++) {
         std::string atom_name(residue_atoms[i]->name);
         if (atom_name == " CA ") {
            ir = atom_to_atom_index(residue_atoms[i]);
            if (ir == -1)
               ir = full_atom_spec_to_atom_index(residue_atoms[i]->GetChainID(),
                                                 residue_atoms[i]->GetSeqNum(),
                                                 residue_atoms[i]->GetInsCode(),
                                                 residue_atoms[i]->name,
                                                 residue_atoms[i]->altLoc);
         }
         // likewise C1'
         if (atom_name == " C1'") {
            ir = atom_to_atom_index(residue_atoms[i]);
            if (ir == -1)
               ir = full_atom_spec_to_atom_index(residue_atoms[i]->GetChainID(),
                                                 residue_atoms[i]->GetSeqNum(),
                                                 residue_atoms[i]->GetInsCode(),
                                                 residue_atoms[i]->name,
                                                 residue_atoms[i]->altLoc);
         }
      }

      if (ir == -1) {
         if (nResidueAtoms > 0) {
            ir = atom_to_atom_index(residue_atoms[0]);
            if (ir == -1)
               ir =  full_atom_spec_to_atom_index(residue_atoms[0]->GetChainID(),
                                                  residue_atoms[0]->GetSeqNum(),
                                                  residue_atoms[0]->GetInsCode(),
                                                  residue_atoms[0]->name,
                                                  residue_atoms[0]->altLoc);
         }
      }
   }
   return ir;
}

coot::atom_spec_t
molecule_class_info_t::intelligent_this_residue_atom(const coot::residue_spec_t &rs) const {

   coot::atom_spec_t atom_spec;
   mmdb::Residue *res_p = get_residue(rs);
   if (res_p) {
      mmdb::Atom *at = intelligent_this_residue_mmdb_atom(res_p);
      if (at) {
         atom_spec = coot::atom_spec_t(at);
      }
   }
   return atom_spec;
}


// If there is a CA or a C1' in this residue then return that atom (pointer)
// atom, if not, then return the index of the first atom in the
// residue.
//
// Return NULL on no atoms in residue.
//
mmdb::Atom *
molecule_class_info_t::intelligent_this_residue_mmdb_atom(mmdb::Residue *res_p) const {

   mmdb::Atom *null_at = NULL;

   if (res_p) {
      mmdb::PPAtom residue_atoms = 0;
      int nResidueAtoms;
      res_p->GetAtomTable(residue_atoms, nResidueAtoms);
      if (nResidueAtoms > 0) {
         for (int i=0; i<nResidueAtoms; i++) {
            std::string atom_name(residue_atoms[i]->name);
            if (atom_name == " CA ") {
               return residue_atoms[i];
            }
            if (atom_name == " C1'") {
               return residue_atoms[i];
            }
         }
         return residue_atoms[0]; // ok, so any atom will do
      }
   }
   return null_at;
}

// Return pointer to atom " CA ", or the first atom in the residue, or
// null (no residue or atoms error):
//
mmdb::Atom *
molecule_class_info_t::atom_intelligent(const std::string &chain_id, int resno,
                                        const std::string &ins_code) const {

   mmdb::Atom *at = NULL;

   if (atom_sel.n_selected_atoms > 0) {
      int selHnd = atom_sel.mol->NewSelection(); // d
      mmdb::PPResidue SelResidue;
      int nSelResidues;

      atom_sel.mol->Select (selHnd, mmdb::STYPE_RESIDUE, 0,
                            chain_id.c_str(),
                            resno, ins_code.c_str(),
                            resno, ins_code.c_str(),
                            "*",  // residue name
                            "*",  // Residue must contain this atom name?
                            "*",  // Residue must contain this Element?
                            "*",  // altLocs
                            mmdb::SKEY_NEW // selection key
                            );

      atom_sel.mol->GetSelIndex(selHnd, SelResidue, nSelResidues);

      if (nSelResidues == 0) {
         std::cout << "INFO:: No selected residues" << std::endl;
      } else {

         mmdb::PPAtom residue_atoms;
         int nResidueAtoms;
         SelResidue[0]->GetAtomTable(residue_atoms, nResidueAtoms);
         if (nResidueAtoms == 0) {
            std::cout << "INFO:: No atoms in residue" << std::endl;
         } else {
            bool found_it = false;
            std::string CA       = " CA "; // PDBv3 FIXME
            std::string C1_prime = " C1'";
            for (int i=0; i<nResidueAtoms; i++) {
               if (std::string(residue_atoms[i]->name) == CA) {
                  at = residue_atoms[i];
                  found_it = true;
                  break;
               }
               if (std::string(residue_atoms[i]->name) == C1_prime) {
                  at = residue_atoms[i];
                  found_it = true;
                  break;
               }
            }
            if (! found_it)
               at = residue_atoms[0];
         }
      }
      atom_sel.mol->DeleteSelection(selHnd);
   }
   return at;
}

// is point close (< 1A) to any atom in the given residue?
bool
molecule_class_info_t::close_to_residue(mmdb::Residue *residue_p, coot::Cartesian point) const {

   bool status = false;
   if (residue_p) {
      if (atom_sel.mol) {
         //
         mmdb::PPAtom residue_atoms = 0;
         int n_residue_atoms;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            coot::Cartesian atom_pt(residue_atoms[iat]->x,
                                    residue_atoms[iat]->y,
                                    residue_atoms[iat]->z);
            double d = (atom_pt - point).amplitude();
            if (d < 1.0) {
               status = true;
               break;
            }
         }
      }
   }
   return status;
}

// residue for spec is missing, return the next residue.
//
mmdb::Residue *
molecule_class_info_t::next_residue_missing_residue(const coot::residue_spec_t &spec) const {

   mmdb::Residue *r = NULL;
   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   mmdb::Chain *chain_p;
   int n_chains = model_p->GetNumberOfChains();
   bool found_chain = false;
   for (int ichain=0; ichain<n_chains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      std::string chain_id = chain_p->GetChainID();
      if (chain_id == spec.chain_id) {
         found_chain = true;
         mmdb::Residue *residue_p;
         for (int ires=0; ires<nres; ires++) {
            residue_p = chain_p->GetResidue(ires);
            if (residue_p->GetSeqNum() > spec.res_no) {
               r = residue_p;
               break;
            }
         }
         if (r)
            break;
      } else {
         // OK, so spec was at the end of the chain, the previous
         // chain, just return the first residue of this chain.
         if (found_chain) {
            for (int ires=0; ires<nres; ires++) {
               r = chain_p->GetResidue(ires);
               break;
            }
         }
      }
      if (r)
         break;
   }
   return r;
}




// ----------------------------------------------------------------------
//               Pointer Atoms
// ----------------------------------------------------------------------

bool
molecule_class_info_t::have_atom_close_to_position(const coot::Cartesian &pos) const {

   bool r = false;

   float close_d = 0.5;
   float close_d_squared = close_d * close_d;
   if (atom_sel.mol) {
      for(int imod = 1; imod <= atom_sel.mol->GetNumberOfModels(); imod++) {
         mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
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
                     if (! at->isTer()) {
                        float dd =
                           (pos.x() - at->x) * (pos.x() - at->x) +
                           (pos.y() - at->y) * (pos.y() - at->y) +
                           (pos.z() - at->z) * (pos.z() - at->z);
                        if (dd < close_d_squared) {
                           r = true;
                           break;
                        }
                     }
                  }
                  if (r) break;
               }
               if (r) break;
            }
            if (r) break;
         }
      }
   }
   return r;
}

void
molecule_class_info_t::add_pointer_atom(coot::Cartesian pos) {


   if (atom_sel.mol) {
      mmdb::Chain *chain_p = water_chain();

      if (! chain_p) {
         // we have to make one then
         chain_p = new mmdb::Chain;
         std::pair<short int, std::string> p = unused_chain_id();
         if (p.first)
            chain_p->SetChainID(p.second.c_str());
         mmdb::Model *model_p = atom_sel.mol->GetModel(1);
         model_p->AddChain(chain_p);
      }

      make_backup("add_pointer_atom()");
      std::string mol_chain_id(chain_p->GetChainID());
      // int ires_prev = chain_p->GetNumberOfResidues();
      int ires_prev = coot::util::max_resno_in_chain(chain_p).second;

      mmdb::Residue *res_p = new mmdb::Residue;
      mmdb::Atom *atom_p = new mmdb::Atom;
      chain_p->AddResidue(res_p);
      atom_p->SetAtomName(" O  ");
      float bf = graphics_info_t::default_new_atoms_b_factor;
      atom_p->SetCoordinates(pos.x(), pos.y(), pos.z(), 1.0, bf);

      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);
      res_p->seqNum = ires_prev + 1;
      res_p->SetResName("HOH");
      coot::hetify_residue_atoms(res_p);

      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      atom_sel = make_asc(atom_sel.mol);
      std::cout << atom_p << " added to molecule" << std::endl;

      have_unsaved_changes_flag = 1;
      make_bonds_type_checked(__FUNCTION__);
   }
}

// This is a bit messy, I'm afraid - we test single atom twice. If you use this for
// a multiatom other than SO4 and P04, you will need to add it to the type test
//
std::pair<bool,std::string>
molecule_class_info_t::add_typed_pointer_atom(coot::Cartesian pos, const std::string &type) {

   bool status = false; // return this
   std::string message;

   bool single_atom = true;

   // std::cout << "INFO:: adding atom of type " << type << " at " << pos << std::endl;
   make_backup("add_typed_pointer_atom");

   if (have_atom_close_to_position(pos)) {
      return std::make_pair(false, std::string("Too close to an existing atom"));
   }

   // we get a chain pointer or NULL, if there is not yet a chain only
   // of the given type:
   coot::atom_name_bits_t bits(type);
   mmdb::Chain *single_type = coot::util::chain_only_of_type(atom_sel.mol, bits.res_name);
   // std::cout << "DEBUG:: single_type returned " << single_type << std::endl;

   // We do different things (e.g adding the chain) if this is a new
   // chain or a pre-existing one, let's set a flag.
   bool pre_existing_chain_flag;
   mmdb::Chain *chain_p;
   if (single_type) {
      chain_p = single_type;
      pre_existing_chain_flag = 1;
   } else {
      chain_p = new mmdb::Chain;
      pre_existing_chain_flag = 0;
   }

   std::pair<short int, std::string> mol_chain_id = unused_chain_id();
   mmdb::Residue *res_p = new mmdb::Residue;

   // type test
   if (type == "PO4") single_atom = 0;
   if (type == "SO4") single_atom = 0;

   if (single_atom) {
      mmdb::Atom *atom_p = new mmdb::Atom;
      float occ;
      if (is_from_shelx_ins_flag)
         occ = 11.0;
      else
         occ = 1.0;
      atom_p->SetCoordinates(pos.x(), pos.y(), pos.z(), occ,
                             graphics_info_t::default_new_atoms_b_factor);
      atom_p->Het = 1; // it's a HETATM.

      if (type == "Water") {

         // special rule for water: we add a water to a water chain if
         // possible

         atom_p->SetAtomName(" O  ");
         atom_p->SetElementName(" O");
         res_p->SetResName("HOH");

         mmdb::Chain *w = water_chain();
         int wresno = 1;

         bool ok_to_add = true;
         if (have_atom_close_to_position(pos))
            ok_to_add = false;

         if (! ok_to_add) {
            graphics_info_t g;
            std::string s = "WARNING:: new atom addition blocked by nearby atom";
            std::cout << s << std::endl;
            g.add_status_bar_text(s);
            message = s;
         } else {

            if (w) {
               // remove a TER atom if it exists on the last residue
               // prior to insertion of a new residue.
               remove_TER_on_last_residue(w);

               // Now add atom to chain w.
               std::pair<short int, int> wresno_pair = next_residue_number_in_chain(w);
               if (wresno_pair.first) {
                  wresno = wresno_pair.second;
               } else {
                  wresno = 1;
               }
               res_p->seqNum = wresno;
               res_p->AddAtom(atom_p);
               w->AddResidue(res_p);
               std::cout << "DEBUG:: " << atom_p << " added to molecule" << std::endl;
               atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
               atom_sel.mol->FinishStructEdit();
               coot::util::pdbcleanup_serial_residue_numbers(atom_sel.mol);
               atom_sel = make_asc(atom_sel.mol);
               have_unsaved_changes_flag = 1;
               make_bonds_type_checked(__FUNCTION__);
               status = true;

            } else {
               // There was no water chain
               res_p->AddAtom(atom_p);
               std::cout << "DEBUG:: " << atom_p << " in new chain added to molecule (and new chain)" << std::endl;
               if (!pre_existing_chain_flag) {
                  chain_p->SetChainID(mol_chain_id.second.c_str());
                  mmdb::Model *model_p = atom_sel.mol->GetModel(1);
                  if (model_p)
                     model_p->AddChain(chain_p);
               }
               res_p->seqNum = 1; // start of a new chain.
               chain_p->AddResidue(res_p);
               atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
               atom_sel.mol->FinishStructEdit();
               // removed replacement of the atom_sel
               coot::util::pdbcleanup_serial_residue_numbers(atom_sel.mol);
               atom_sel = make_asc(atom_sel.mol);
               have_unsaved_changes_flag = 1;
               make_bonds_type_checked(__FUNCTION__);
               status = true;
            }
         }
      } else {

         // Not water
         std::string element = "";
         if (mol_chain_id.first || pre_existing_chain_flag) {
            if (bits.filled) {

               bits.SetAtom(atom_p, res_p);
               if (true)
                  std::cout << "debug:: bits.SetAtom() called with atom " << coot::atom_spec_t(atom_p)
                            << " and residue " << coot::residue_spec_t(res_p)
                            << " with residue name \"" << res_p->GetResName() << "\"" << std::endl;
               res_p->AddAtom(atom_p);
               std::cout << atom_p << " added to molecule" << std::endl;
               if (! pre_existing_chain_flag) {
                  chain_p->SetChainID(mol_chain_id.second.c_str());
                  atom_sel.mol->GetModel(1)->AddChain(chain_p);
               }
               std::pair<short int, int> ires_prev_pair = coot::util::max_resno_in_chain(chain_p);
               int previous_max = 0;
               if (ires_prev_pair.first) { // was not an empty chain
                  previous_max =  ires_prev_pair.second;
                  res_p->seqNum = previous_max + 1;
               } else {

                  // was an empty chain.  Handle the shelx case:

                  if (! is_from_shelx_ins_flag) {
                     res_p->seqNum = 1 ; // start of a new chain.
                  } else {
                     // in a shelx molecule, we can't make the residue
                     // number 1 because there are no chains.  We need to
                     // make the residue number bigger than the biggest
                     // residue number so far.
                     ires_prev_pair = coot::util::max_resno_in_molecule(atom_sel.mol);
                     if (ires_prev_pair.first) {
                        res_p->seqNum = ires_prev_pair.second + 1;
                     } else {
                        res_p->seqNum = 1;
                     }
                  }
               }
            }

            // Add this element to the sfac (redundancy check in the addition function
            if (is_from_shelx_ins_flag) {
               shelxins.add_sfac(element);
            }
            chain_p->AddResidue(res_p);
            atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
            atom_sel.mol->FinishStructEdit();
            atom_sel = make_asc(atom_sel.mol);
            have_unsaved_changes_flag = 1;
            make_bonds_type_checked(__FUNCTION__);
            status = true;
         } else {
            std::cout << "WARNING:: Can't find new chain for new atom\n";
         }
      } // type was water, or not

   } else {
      // multi atom:

      if (mol_chain_id.first || pre_existing_chain_flag) {
         add_pointer_multiatom(res_p, pos, type);
         coot::hetify_residue_atoms(res_p);
         if (! pre_existing_chain_flag) {
            chain_p->SetChainID(mol_chain_id.second.c_str());
            atom_sel.mol->GetModel(1)->AddChain(chain_p);
         }
         std::pair<short int, int> ires_prev_pair = coot::util::max_resno_in_chain(chain_p);
         int previous_max = 0;
         if (ires_prev_pair.first) { // was not an empty chain
            previous_max =  ires_prev_pair.second;
         }
         res_p->seqNum = previous_max + 1;

         chain_p->AddResidue(res_p);
         atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
         atom_sel.mol->FinishStructEdit();
         atom_sel = make_asc(atom_sel.mol);
         have_unsaved_changes_flag = 1;
         make_bonds_type_checked(__FUNCTION__);
         status = true;
      } else {
         std::cout << "WARNING:: Can't find new chain for new atom\n";
      }
   }

   // or we could just use update_molecule_after_additions() there.
   return std::make_pair(status, message);
}


// return status [1 means "usable"] and a chain id [status = 0 when
// there are 2*26 chains...]
//
std::pair<bool, std::string>
molecule_class_info_t::unused_chain_id() const {

   std::string r("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
   std::pair<bool, std::string> s(false,"");
   mmdb::Chain *chain_p;
   if (atom_sel.n_selected_atoms > 0) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(1);
      int nchains = model_p->GetNumberOfChains();

      for (int ich=0; ich<nchains; ich++) {
         chain_p = model_p->GetChain(ich);
         std::string this_chain_id = chain_p->GetChainID();
         std::string::size_type idx = r.find(this_chain_id);
         if (idx != std::string::npos) {
            r.erase(idx,1);
         }
      }
      if (r.length()) {
         s.first = true;
         s.second = r.substr(0,1);
      }
   } else {
      s.first = true;
      s.second = "A";
   }
   return s;
}

void
molecule_class_info_t::add_pointer_multiatom(mmdb::Residue *res_p,
                                             const coot::Cartesian &pos, const std::string &type) {

   coot::Cartesian p;
   float bf = graphics_info_t::default_new_atoms_b_factor;
   res_p->SetResName(type.c_str());
   if (type == "SO4") {
      mmdb::Atom *atom_p;

      atom_p = new mmdb::Atom;
      p = pos + coot::Cartesian(0.000, 0.000, 0.088);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" S  ");
      atom_p->SetElementName(" S");
      res_p->AddAtom(atom_p);

      atom_p = new mmdb::Atom;
      p = pos + coot::Cartesian( 1.227, 0.000, -0.813);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O1 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);

      atom_p = new mmdb::Atom;
      p = pos + coot::Cartesian(-1.227, 0.000, -0.813);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O2 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);

      atom_p = new mmdb::Atom;
      p = pos + coot::Cartesian(  0.000, -1.263, 0.740);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O3 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);

      atom_p = new mmdb::Atom;
      p = pos + coot::Cartesian( 0.000, 1.263, 0.740);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O4 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);
   }

   if (type == "PO4") {
      mmdb::Atom *atom_p;

      atom_p = new mmdb::Atom;
      p = pos + coot::Cartesian(0.000, 0.021, 0.036);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" P  ");
      atom_p->SetElementName(" P");
      res_p->AddAtom(atom_p);

      atom_p = new mmdb::Atom;
      p = pos + coot::Cartesian( 1.315,   0.599,  -0.691 );
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O1 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);

      atom_p = new mmdb::Atom;
      p = pos + coot::Cartesian( -1.315,   0.599,  -0.691);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O2 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);

      atom_p = new mmdb::Atom;
      p = pos + coot::Cartesian( 0.000,  -1.587,  -0.055);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O3 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);

      atom_p = new mmdb::Atom;
      p = pos + coot::Cartesian( 0.000,   0.434,   1.457);
      atom_p->SetCoordinates(p.x(), p.y(), p.z(), 1.0, bf);
      atom_p->SetAtomName(" O4 ");
      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);
   }

}


// ----------------------------------------------------------------------
//                 Save yourself
// ----------------------------------------------------------------------
//
// return 0 on success.
//
// optional args save_hydrogens and save_aniso_records.
int
molecule_class_info_t::save_coordinates(const std::string &filename,
                                        bool save_hydrogens,
                                        bool save_aniso_records,
                                        bool save_conect_records) {

   int ierr = 0;
   std::string ext = coot::util::file_name_extension(filename);
   if (coot::util::extension_is_for_shelx_coords(ext)) {
      std::pair<int, std::string> status_pair = write_shelx_ins_file(filename);
      // we need to reverse the logic of the status, 1 is good for write_shelx_ins_file()
      if (status_pair.first != 1)
         ierr = 1;
   } else {
      mmdb::byte bz = mmdb::io::GZM_NONE;

      bool write_as_cif = false;
      if (coot::is_mmcif_filename(filename))
         write_as_cif = true;

      ierr = write_atom_selection_file(atom_sel, filename, write_as_cif, bz,
                                       save_hydrogens, save_aniso_records,
                                       save_conect_records);
   }

   if (ierr) {
      std::cout << "WARNING:: Coordinates write to " << filename
                << " failed!" << std::endl;
      std::string ws = "WARNING:: export coords: There was an error ";
      ws += "in writing ";
      ws += filename;
      graphics_info_t g;
      g.info_dialog(ws);
   } else {
      std::cout << "INFO:: saved coordinates " << filename << std::endl;
      have_unsaved_changes_flag = 0;

      // Now we have updated the molecule name, how shall we restore
      // this from the state file?
      std::vector<std::string> strings;
      strings.push_back("coot.handle-read-draw-molecule");
      strings.push_back(single_quote(coot::util::intelligent_debackslash(filename)));
      save_state_command_strings_ = strings;

      name_ = filename;  // hmm... // update go to atom widget now? FIXME.
      std::string::size_type icoot = filename.rfind("-coot-");
      if (icoot != std::string::npos) {
         coot_save_index++;
      }
      update_mol_in_display_control_widget();  // FIXME.
   }
   return ierr;
}



// Return 1 on yes, unsaved changes present,
//        0 on no
int
molecule_class_info_t::Have_unsaved_changes_p() const {
   if (has_model())
      return have_unsaved_changes_flag;
   else
      return 0;
}


// ----------------------------------------------------------------------
//               Baton Atoms
// ----------------------------------------------------------------------


// Recall that the chain is set by the creation of the empty molecule in
// graphics_info_t::baton_build_atoms_molecule();
// direction_flag is +1 for forward building, -1 for backwards direction.
//
// In the case of direction_flag being negative, I think that residues
// (atoms) will go into the chain with their seqNums in decreasing
// order.  The may need to be sorted later, I think (maybe not).
//
mmdb::Atom *
molecule_class_info_t::add_baton_atom(coot::Cartesian pos,
                                      int istart_resno,
                                      const std::string &chain_id,
                                      short int iresno_active,
                                      short int direction_flag) {

   mmdb::Model *model_p = atom_sel.mol->GetModel(1);
   int nchains = model_p->GetNumberOfChains();
   if (nchains == 0) {
      std::cout << "failed to add baton atom" << std::endl;
      return NULL;
   }

   make_backup("add_baton_atom()");
   mmdb::Chain *chain_p = NULL;

   // now look at the chains of this model, and find the chain that
   // has the same chain-id as is passed. If there is no such chain,
   // make one and set the chain_id
   //
   for (int ich=0; ich<nchains; ich++) {
      mmdb::Chain *chain_p_local = model_p->GetChain(ich);
      std::string chain_id_local = chain_p_local->GetChainID();
      if (chain_id_local == chain_id) {
         chain_p = chain_p_local;
         break;
      }
   }

   if (! chain_p) {
      // 20100125 can't do this today, the second arg to the
      // constructor is a char * (not const)
      // chain_p = new mmdb::Chain(model_p, chain_id.c_str()); // does the add for us
      chain_p = new mmdb::Chain;
      chain_p->SetChainID(chain_id.c_str());
      model_p->AddChain(chain_p);
   }


   std::string mol_chain_id(chain_p->GetChainID());
   int n_res = chain_p->GetNumberOfResidues();


   // if this is the first atom to be added in the chain, we get the
   // seqnum from the passed istart_resno otherwise it is the seqnum
   // of the previous residue, plus or minus one.
   //
   int this_res_seqnum;
   if (n_res == 0) {
      this_res_seqnum = istart_resno;
   } else {

      if (iresno_active == 0) {
         int ires_prev = chain_p->GetResidue(n_res-1)->seqNum; // seqnum of the last
                                                               // residue in chain.
         this_res_seqnum = ires_prev + 1*direction_flag;
      } else {
         this_res_seqnum = istart_resno;
      }
   }

   mmdb::Residue *res_p = new mmdb::Residue;
   mmdb::Atom *atom_p = new mmdb::Atom;
   chain_p->AddResidue(res_p);
   atom_p->SetAtomName(" CA ");
   atom_p->SetCoordinates(pos.x(), pos.y(), pos.z(), 1.0,
                          graphics_info_t::default_new_atoms_b_factor);

   atom_p->SetElementName(" C");
   res_p->AddAtom(atom_p);
   res_p->seqNum = this_res_seqnum;
   res_p->SetResName("ALA");

   atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
   atom_sel.mol->FinishStructEdit();
   atom_sel = make_asc(atom_sel.mol);
   std::cout << atom_p << " added to molecule" << std::endl;

   have_unsaved_changes_flag = 1;
   make_ca_bonds(2.4, 4.7);

   return atom_p;
}

// Return a vector of upto 3 positions of the most latestly added
// atoms with the most lastest atom addition (that is the passed atom)
// in the back() slot of the vector.
//
std::vector<clipper::Coord_orth>
molecule_class_info_t::previous_baton_atom(mmdb::Atom* latest_atom_addition,
                                           short int direction) const {
   std::vector<clipper::Coord_orth> positions;
   int direction_sign = +1;

   if (direction == 1) { // building forward, look in negative
                         // direction for previously build atoms.
      direction_sign = +1;
   } else {
      direction_sign = -1; // building backward, look in positive
                           // direction for previously build atoms.
   }
   int ires_last_atom = latest_atom_addition->GetSeqNum();

   char *chain = latest_atom_addition->GetChainID();
   // does the CA for the (ires_last_atom-2) exist?
   int selHnd = atom_sel.mol->NewSelection();

   atom_sel.mol->SelectAtoms(selHnd, 0, chain,
                             ires_last_atom-2*direction_sign, "*",
                             ires_last_atom-2*direction_sign, "*",
                             "*", " CA ", "*", "*");

   int nSelAtoms;
   mmdb::PPAtom local_SelAtom;
   atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);

   if (nSelAtoms == 0) {
      std::cout << "residue with sequence number " << ires_last_atom - 2*direction_sign
                << " not found for ires_last_atom = " << ires_last_atom
                << " with direction_sign = " << direction_sign << "\n";
   } else {
      positions.push_back(clipper::Coord_orth(local_SelAtom[0]->x,
                                              local_SelAtom[0]->y,
                                              local_SelAtom[0]->z));
   }
   atom_sel.mol->DeleteSelection(selHnd);

   // rinse, lather, repeat...

   // does the CA for the (ires_last_atom-1) exist?
   selHnd = atom_sel.mol->NewSelection();

   atom_sel.mol->SelectAtoms(selHnd, 0, chain,
                             ires_last_atom-direction_sign, "*",
                             ires_last_atom-direction_sign, "*",
                             "*", " CA ", "*", "*");

   atom_sel.mol->GetSelIndex(selHnd, local_SelAtom, nSelAtoms);

   if (nSelAtoms == 0) {
      std::cout << "residue with sequence number " << ires_last_atom - direction_sign
                << " not found\n";
   } else {
      positions.push_back(clipper::Coord_orth(local_SelAtom[0]->x,
                                              local_SelAtom[0]->y,
                                              local_SelAtom[0]->z));
   }
   atom_sel.mol->DeleteSelection(selHnd);

   // And finally this one is guaranteed to exist:
   //
   positions.push_back(clipper::Coord_orth(latest_atom_addition->x,
                                           latest_atom_addition->y,
                                           latest_atom_addition->z));

   return positions;

}

#include "build/CalphaBuild.hh"

std::vector<coot::scored_skel_coord>
molecule_class_info_t::next_ca_by_skel(const std::vector<clipper::Coord_orth> &previous_ca_positions,
                                       const clipper::Coord_grid &coord_grid_start,
                                       short int use_coord_grid_start_flag,
                                       float ca_ca_bond_length,
                                       float map_cut_off,
                                       int max_skeleton_search_depth) const {

   std::vector<coot::scored_skel_coord> t;
   coot::CalphaBuild buildca(max_skeleton_search_depth);

//     std::cout << "DEBUG:: ------ "
//               << "in molecule_class_info_t::next_ca_by_skel skeleton_treenodemap_is_filled is "
//               << skeleton_treenodemap_is_filled << " for molecule " << imol_no << std::endl;


   if (skeleton_treenodemap_is_filled) {
      t = buildca.next_ca_by_skel(previous_ca_positions,
                                  coord_grid_start,
                                  use_coord_grid_start_flag,
                                  ca_ca_bond_length,
                                  xskel_cowtan, xmap,
                                  map_cut_off,
                                  skeleton_treenodemap);
   } else {
      std::cout << "treenodemap is not filled" << std::endl;
   }
   return t;
}

#include <time.h>

// ----------------------------------------------------------------------
//               Dummy Atoms (not bonded)
// ----------------------------------------------------------------------
void
molecule_class_info_t::add_dummy_atom(coot::Cartesian pos) {

   // 20150803-PE FIXME - pass Geom_p().
   graphics_info_t g;
   coot::protein_geometry *geom_p = g.Geom_p();

   int nchains = atom_sel.mol->GetNumberOfChains(1);

   if (nchains != 1) {
      std::cout << "failed to add dummy atom" << std::endl;
      return;
   }

   make_backup("add_dummy_atom");

   mmdb::Chain *chain_p = atom_sel.mol->GetChain(1,0);

   std::string mol_chain_id(chain_p->GetChainID());
   int ires_prev = chain_p->GetNumberOfResidues();

   mmdb::Residue *res_p = new mmdb::Residue;
   mmdb::Atom *atom_p = new mmdb::Atom;
   chain_p->AddResidue(res_p);
   atom_p->SetAtomName(" DUM");
   atom_p->SetCoordinates(pos.x(), pos.y(), pos.z(), 1.0,
                          graphics_info_t::default_new_atoms_b_factor);

   atom_p->SetElementName(" O");
   res_p->AddAtom(atom_p);
   res_p->seqNum = ires_prev + 1;
   res_p->SetResName("DUM");

   atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
   atom_sel.mol->FinishStructEdit();
   atom_sel = make_asc(atom_sel.mol);
   // std::cout << atom_p << " added to molecule" << std::endl;

   have_unsaved_changes_flag = 1;
   makebonds(0.0, 0.0, geom_p);

}



// backups:

// Backup filename: return a stub.
//
// move to molecule-class-info-backups.cc
//
std::string
molecule_class_info_t::get_save_molecule_filename(const std::string &dir) {

   auto replace_char = [] (const std::string &s, char a) {
                          std::string r = s;
                          int slen = s.length();
                          for (int i=0; i<slen; i++) {
                             if (r[i] == a)
                                r[i] = '_';
                          }
                          return r;
                       };

   graphics_info_t g;
   bool decolonify = g.decoloned_backup_file_names_flag;
   std::string t_name_1 = name_;
   if (g.unpathed_backup_file_names_flag)
      t_name_1 = name_for_display_manager();
   std::string t_name_2 = replace_char(t_name_1, '/');
   std::string t_name_3 = replace_char(t_name_2, ' ');
#ifdef WINDOWS_MINGW
   std::string t_name_x = replace_char(t_name_3, '\\');
   std::string t_name_y = replace_char(t_name_x, ':');
   t_name_3 = t_name_y;
#endif

   if (save_time_string.empty()) {
      time_t t;
      time(&t);
      char *chars_time = ctime(&t);
      int l = strlen(chars_time);
      save_time_string = chars_time;
      if (! save_time_string.empty()) {
         std::string::size_type l = save_time_string.length();
         save_time_string = save_time_string.substr(0, l-1);
      }
      save_time_string = replace_char(save_time_string, ' ');
      save_time_string = replace_char(save_time_string, '/');
      if (decolonify)
         save_time_string = replace_char(save_time_string, ':');
   }
   std::string time_string = save_time_string;
   std::string t_name_4 = t_name_3 + "_" + time_string;

   std::string index_string = coot::util::int_to_string(history_index);
   std::string t_name_5 = t_name_4 + "_modification_" + index_string;

   std::string extension = ".pdb";
   if (coot::is_mmcif_filename(name_))
      extension = ".cif";
   if (is_from_shelx_ins_flag)
      extension = ".res";
   if (g.backup_compress_files_flag)
      extension += ".gz";

   std::string t_name_6 = t_name_5 + extension;

   std::string save_file_name = coot::util::append_dir_file(dir, t_name_6);
   return save_file_name;

}

// Return like mkdir: mkdir returns zero on success, or -1 if an  error  occurred
//
// if it already exists as a dir, return 0 of course.
//
int
molecule_class_info_t::make_maybe_backup_dir(const std::string &backup_dir) const {

   return coot::util::create_directory(backup_dir);
}

#include "utils/xdg-base.hh"

// Ignore return value.
//
// If successful, increase history_index and if not in a backup
// increase max_history_index too.
//
// move this function to the mci-backups file.
int
molecule_class_info_t::make_backup(const std::string &descr) { // changes history details

   graphics_info_t g;
   if (backup_this_molecule) {
      xdg_t xdg;
      std::string coot_backup_dir = xdg.get_cache_home().append("coot-backup").string();
      std::string backup_dir(coot_backup_dir);

      //shall we use the environment variable instead?
      char *env_var = getenv("COOT_BACKUP_DIR");
      if (env_var) {
         struct stat buf;

         // we better debackslash the directory (for windows)
         std::string tmp_dir = env_var;
         tmp_dir = coot::util::intelligent_debackslash(tmp_dir);
         int err = stat(tmp_dir.c_str(), &buf);

         if (!err) {
            if (! S_ISDIR(buf.st_mode)) {
               env_var = NULL;
            }
         } else {
            env_var = NULL;
         }
      }

      if (env_var)
         backup_dir = env_var;

      if (atom_sel.mol) {
         int dirstat = make_maybe_backup_dir(backup_dir);

         if (dirstat != 0) {
            // fallback to making a directory in $HOME
            std::string home_dir = coot::get_home_dir();
            if (! home_dir.empty()) {
               backup_dir = coot::util::append_dir_dir(home_dir, "coot-backup");
               dirstat = make_maybe_backup_dir(backup_dir);
               if (dirstat != 0) {
                  std::cout << "WARNING:: backup directory "<< backup_dir
                            << " failure to exist or create" << std::endl;
               } else {
                  std::cout << "INFO using backup directory " << backup_dir << std::endl;
               }
            } else {
               std::cout << "WARNING:: backup directory "<< backup_dir
                         << " failure to exist or create" << std::endl;
            }
         }

         if (dirstat == 0) {
            // all is hunkey-dorey.  Directory exists.

            std::string backup_file_name = get_save_molecule_filename(backup_dir);
            logger.log(log_t::INFO, "backup file-name:", backup_file_name);

            mmdb::byte gz;
            if (g.backup_compress_files_flag) {
               gz = mmdb::io::GZM_ENFORCE;
            } else {
               gz = mmdb::io::GZM_NONE;
            }

            // Writing out a modified binary mmdb like this results in the
            // file being unreadable (crash in mmdb read).
            //
            int istat;
            if (! is_from_shelx_ins_flag) {
               bool write_as_cif = false;
               if (coot::is_mmcif_filename(name_))
                  write_as_cif = true;

               istat = write_atom_selection_file(atom_sel, backup_file_name, write_as_cif, gz);

               // WriteMMDBF returns 0 on success, else mmdb:Error_CantOpenFile (15)
               if (istat) {
                  std::string warn;
                  warn = "WARNING:: WritePDBASCII failed! Return status ";
                  warn += istat;
                  // 2025-12-22-PE - the function should return a value
                  // with this warning message. This is not the place
                  // for GUI code.
                  g.info_dialog_and_text(warn);
               }
            } else {
               std::pair<int, std::string> p = write_shelx_ins_file(backup_file_name);
               istat = p.first;
            }

            save_history_file_name(backup_file_name, descr);
            if (history_index == max_history_index)
               max_history_index++;
            history_index++;
         }
      } else {
         std::cout << "WARNING:: BACKUP:: Ooops - no atoms to backup for this empty molecule"
                   << std::endl;
      }
   } else {
      // Occasionally useful but mostly tedious...
      // std::cout << "INFO:: backups turned off on this molecule"
      // << std::endl;
   }
   return 0;
}


void
molecule_class_info_t::save_history_file_name(const std::string &file_name,
                                              const std::string &description) {

   // this is called only from make_backup() (above).

#if 0 // 2025-12-23-PE strange logic that I now don't understand

   // First, history_index is zero and the vec is zero,
   // normal service, then another backup: history_index is 1 and vec is 1.
   //
   if (history_index == int(history_filename_vec.size())) {
      coot::backup_file_info_t(file);
      // history_filename_vec.push_back(file);
   } else {
      // we have gone back in history.
      //
      if (history_index < int(history_filename_vec.size())) {
         history_filename_vec[history_index].backup_file_name = file;
      }
   }
#endif

   coot::backup_file_info_t bfi(file_name, description);
   bfi.imol = imol_no;
   bfi.valid_status = true;
   clock_gettime(CLOCK_REALTIME, &bfi.ctime);

   history_filename_vec.push_back(bfi);
}

// restore from (previous) backup
bool molecule_class_info_t::restore_from_backup(int history_offset,
                                                const std::string &cwd) {

   // return success status
   bool status = false;

   // consider passing this:
   bool v2_convert_flag = graphics_info_t::convert_to_v2_atom_names_flag;
   bool allow_duplseqnum = graphics_info_t::allow_duplseqnum;

   int hist_vec_index = history_index + history_offset;
   if (int(history_filename_vec.size()) > hist_vec_index) {
      // std::cout << "INFO:: restoring from backup " << history_filename_vec.size()
      //           << " " << history_index << std::endl;
      logger.log(log_t::INFO, "restoring from backup", history_filename_vec.size(),
            "history index: ", history_index);
      std::string save_name = name_;
      if (hist_vec_index < int(history_filename_vec.size()) &&
          hist_vec_index >= 0) {
         std::string filename = history_filename_vec[hist_vec_index].backup_file_name;
         //      history_index = hist_index;
         short int reset_rotation_centre = 0;
         // handle_read_draw_molecule uses graphics_info_t::n_molecules
         // to determine its molecule number.  We don't want it to
         // change.
         int save_imol = imol_no;
         // similarly, it messes with the save_state_command_strings_, we
         // don't want that either:
         std::vector<std::string> save_save_state = save_state_command_strings_;
         short int is_undo_or_redo = 1;

         handle_read_draw_molecule(imol_no, filename, cwd,
                                   graphics_info_t::Geom_p(),
                                   reset_rotation_centre,
                                   is_undo_or_redo,
                                   allow_duplseqnum,
                                   v2_convert_flag,
                                   bond_width,
                                   Bonds_box_type(),
                                   false);
         save_state_command_strings_ = save_save_state;
         imol_no = save_imol;
         name_ = save_name;
         status = true;
      }
   } else {
      std::cout << "not restoring from backup because "
                << history_filename_vec.size()
                << " " << history_index << std::endl;
   }
   return status;
}


// I need to write an essay on how that backup system works.
//
// Insight: if we are at hist_index = max_hist_index (i.e. not in a
// backup) then when an undo is requested, we should make a backup.
// This makes the indexing a bit tricky,
//
// So imagine this situation:
//
//             hist_index max_hist_index                filenames filled
// pepflip         1         1                              [0]
// rotate          2         2                              [0,1]
// undo            3         3 [first step is a backup]     [0,1,2]
//
// filenames get routinely pushed back onto history_filename_vec
// (i.e. not in an undo situation)
//
// So imagine we make one mod, then the history_filename_vec size() is 1.
// on undo:
//    we restore from backup using history_filename_vec index 0.
//    we have added to history_filename_vec [now has size 2] in this proceedure
//    history_index was 1 on starting the undo
//    at end of undo it is 0.
//
// how do we redo that?
//    (obviously) we restore from backup using history_filename_vec index 1.
//    we have not added to history_filename_vec [now has size 2]
//    history_index was 0 on starting the redo.
//    at end of redo, history_index is 1
//
// So having done 2 mods:
//
// on undo:
//    restore from backup using history_filename_vec index 1
//    we have added to history_filename_vec [now has size 3] in this proceedure
//    history_index was 2 on starting undo
//    at end of undo it is 1.
//
// on redo:
//    restore from backup using history_filename_vec index 2
//    we have not added to history_filename_vec [now has size 3]
//    history_index was 1 on starting the redo
//    at end of redo, history_index was 2


//
// [It would be cool to have the Redo button greyed out when there are
// no redos availabile (set its state when either it or undo is
// pressed)]
//
// initially it should be greyed out (insensitive).

// restore from (next) backup
int
molecule_class_info_t::apply_undo(const std::string &cwd) {

   int state = 0;
//    std::cout << std::endl << "DEBUG:: in apply undo start hist_index: "
//              << history_index
//              << " max_history_index: " << max_history_index << std::endl;

   if (history_index > 0) {
      int offset = -1;
      if (history_index == max_history_index) {
         make_backup("apply_undo"); // increments history_index
         offset--;
      }
      state = 1;
      restore_from_backup(offset, cwd);
      history_index += offset;

      // So that we don't get asked to save the molecule on exist when
      // we have reverted all our modifications:
      //
      if (history_index == 0) {
         have_unsaved_changes_flag = 0;
      }
   }

   std::cout << "DEBUG:: apply_undo: (end) history_index: " <<
      history_index << " max_history_index: " << max_history_index << std::endl;

   return state;

}

int
molecule_class_info_t::apply_redo(const std::string &cwd) {

   int state = 0;
   if (history_index < max_history_index) {
      // std::cout << "DEBUG:: molecule applying redo " << history_index << std::endl;

      // When there are 3 backups made and we are viewing molecule 2,
      // we don't want to restore from history_filename_vec[3]:
      //
      if (int(history_filename_vec.size()) > (history_index + 1)) {
         restore_from_backup(+1, cwd);
         history_index++;
         state = 1;
         have_unsaved_changes_flag = 1;
      } else {
         std::cout << "Not redoing history file vec: " << history_filename_vec.size()
                   << " " << history_index << std::endl;
      }
   } else {
      std::cout << "Not redoing history: " << max_history_index
                << " " << history_index << std::endl;
   }
   return state;
}



// For model view (go to atom)
//
std::vector<coot::model_view_residue_button_info_t>
molecule_class_info_t::model_view_residue_button_labels() const {

   std::vector<coot::model_view_residue_button_info_t> v;

   if (atom_sel.n_selected_atoms > 0) {

      int nchains = atom_sel.mol->GetNumberOfChains(1);

      if (nchains < 1) {
         std::cout << "failed to find chains for atom in "
                   << " model_view_residue_button_info_t" << std::endl;
      } else {

         graphics_info_t g;
         mmdb::Model *model_p = atom_sel.mol->GetModel(1);

         mmdb::Chain *chain;
         // run over chains of the existing mol
         nchains = model_p->GetNumberOfChains();
         if (nchains <= 0) {
            std::cout << "bad nchains in model_view_residue_button_info_t: "
                      << nchains << std::endl;
         } else {
            for (int ichain=0; ichain<nchains; ichain++) {
               chain = model_p->GetChain(ichain);
               if (chain == NULL) {
                  // This should not be necessary. It seem to be a
                  // result of mmdb corruption elsewhere - possibly
                  // DeleteChain in update_molecule_to().
                  std::cout << "NULL chain in model_view_residue_button_info_t: "
                            << std::endl;
               } else {
                  int nres = chain->GetNumberOfResidues();
                  for (int ires=0; ires<nres; ires++) {
                     mmdb::PResidue residue_p = chain->GetResidue(ires);
                     std::string button_label =
                        g.int_to_string(residue_p->GetSeqNum());
                     button_label += " ";
                     button_label += residue_p->GetChainID();
                     button_label += " ";
                     button_label += residue_p->name;

                     v.push_back(coot::model_view_residue_button_info_t(button_label,
                                                                        residue_p));
                  }
               }
            }
         }
      }
   }
   return v;
}

// return vector of atom list (aka button) info for this residue
//
std::vector<coot::model_view_atom_button_info_t>
molecule_class_info_t::model_view_atom_button_labels(const std::string &chain_id,
                                                     int seqno,
                                                     const std::string &ins_code) const {

   graphics_info_t g;
   std::vector<coot::model_view_atom_button_info_t> v;

   // protection against the molecule having been deleted after the
   // gtklist widget was created:

   if (atom_sel.n_selected_atoms > 0) {

      mmdb::Chain *chain;

      // first we have to find the residue res_p (from which we will get the atoms)
      //
      int nchains = atom_sel.mol->GetNumberOfChains(1);
      for (int ichain=0; ichain<nchains; ichain++) {
         chain = atom_sel.mol->GetChain(1, ichain);
         if (chain == NULL) {
            // This should not be necessary. It seem to be a result of
            // mmdb corruption elsewhere - possibly DeleteChain in
            // update_molecule_to().
            std::cout << "ERROR getting chain in model_view_atom_button_info_t\n";
         } else {
            std::string residue_chain_id(chain_id); // passed from residue list
            std::string mol_chain_id(chain->GetChainID());
            if (residue_chain_id == mol_chain_id) {
               int nres = chain->GetNumberOfResidues();
               for (int ires=0; ires<nres; ires++) {
                  mmdb::PResidue res_p = chain->GetResidue(ires);
                  if (res_p->GetSeqNum() == seqno) {
                     std::string ins_code_res(res_p->GetInsCode());
                     if (ins_code_res == ins_code) {

                        mmdb::PPAtom residue_atoms;
                        int nResidueAtoms;

                        res_p->GetAtomTable(residue_atoms, nResidueAtoms);
                        for (int i=0; i<nResidueAtoms; i++) {
                           if (! residue_atoms[i]->isTer()) {
                              std::string button_label = residue_atoms[i]->name;
                              std::string altConf = residue_atoms[i]->altLoc;
                              if (altConf != "") {
                                 button_label += ",";
                                 button_label += altConf;
                              }
                              button_label += " occ=";
                              button_label += g.float_to_string(residue_atoms[i]->occupancy);
                              button_label += " bf=";
                              button_label += g.float_to_string(residue_atoms[i]->tempFactor);

                              v.push_back(coot::model_view_atom_button_info_t(button_label, residue_atoms[i]));
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return v;
}


std::vector<coot::model_view_atom_tree_chain_t>
molecule_class_info_t::model_view_residue_tree_labels(bool include_water_residue_flag,
                                                      bool show_ligands_only_flag) const {

   std::vector<coot::model_view_atom_tree_chain_t> v;

   if (atom_sel.n_selected_atoms > 0) {

      for(int imodel = 1; imodel<=atom_sel.mol->GetNumberOfModels(); imodel++) {
         int nchains = atom_sel.mol->GetNumberOfChains(imodel);
         for (int ichain=0; ichain<nchains; ichain++) {

            mmdb::Chain *chain_p = atom_sel.mol->GetChain(imodel, ichain);
            if (chain_p) {
               std::string chain_label("Chain ");
               chain_label += chain_p->GetChainID();
               v.push_back(coot::model_view_atom_tree_chain_t(chain_label));

               int nres = chain_p->GetNumberOfResidues();
               mmdb::PResidue residue_p;
               for (int ires=0; ires<nres; ires++) {
                  residue_p = chain_p->GetResidue(ires);
                  std::string label = residue_p->GetChainID();
                  std::string res_name(residue_p->GetResName());
                  label += " ";
                  label += coot::util::int_to_string(residue_p->GetSeqNum());
                  label += residue_p->GetInsCode();
                  label += " ";
                  label += residue_p->name;
                  bool is_water_flag = false; // gets reset in this loop
                  bool is_standard_residue_type = false;
                  if (coot::util::is_standard_residue_name(res_name)) { // amino acid and rna and dna polymer
                     is_standard_residue_type = true;
                  } else {
                     if (res_name == "HOH" || res_name == "DOD") {
                        label = std::string("<i>") + label + std::string("</i>");
                        is_water_flag = true;
                     } else {
                        label = std::string("<b>") + label + std::string("</b>");
                     }
                  }
                  if (! is_water_flag || include_water_residue_flag) {
                     if (show_ligands_only_flag) {
                        if (is_standard_residue_type) {
                           // do nothing, Clemens mode
                        } else {
                           coot::model_view_atom_tree_item_info_t res(label, residue_p);
                           v.back().add_residue(res);
                        }
                     } else {
                        coot::model_view_atom_tree_item_info_t res(label, residue_p);
                        v.back().add_residue(res);
                     }
                  }
               }
            }
         }
      }
   }
   return v;
}




// Return 0 on failure.
short int
molecule_class_info_t::move_std_residue(mmdb::Residue *moving_residue,
                                        const mmdb::Residue *reference_residue) const {

   std::map<std::string, clipper::RTop_orth> rtops =
      coot::util::get_ori_to_this_res((mmdb::Residue *)reference_residue);

   short int istat = 0; // success

   if (!reference_residue) {
      std::cout << "This should not happen!" << std::endl;
      std::cout << "null reference residue in move_std_residue" << std::endl;
   } else {

      if (rtops.size()) { // successful attempt to get the matrix
         mmdb::PPAtom residue_atoms = NULL;
         int nResidueAtoms;
         moving_residue->GetAtomTable(residue_atoms, nResidueAtoms);
         if (nResidueAtoms == 0) {
            std::cout << " something broken in atom residue selection in ";
            std::cout << "mutate, got 0 atoms" << std::endl;
            istat = 0;
         } else {
            istat = 1;

            if (false)
               std::cout << "DEBUG:: move_std_residue: " << nResidueAtoms
                         << " atoms in residue "
                         << moving_residue << " " << moving_residue->seqNum << " "
                         << moving_residue->GetChainID() << std::endl;

            for (int iat=0; iat<nResidueAtoms; iat++) {
               if (residue_atoms[iat]) {
                  clipper::Coord_orth co(residue_atoms[iat]->x,
                                         residue_atoms[iat]->y,
                                         residue_atoms[iat]->z);
                  std::string alt_conf = residue_atoms[iat]->altLoc;

                  std::map<std::string, clipper::RTop_orth>::const_iterator it = rtops.find(alt_conf);

                  if (it != rtops.end()) {
                     clipper::Coord_orth rotted = co.transform(it->second); // an rtop
                     residue_atoms[iat]->x = rotted.x();
                     residue_atoms[iat]->y = rotted.y();
                     residue_atoms[iat]->z = rotted.z();
                  }
               } else {
                  istat = 0;
                  std::cout << "ERROR:: null residue atom in moving residue in move_std_residue: iat: "
                            << iat << std::endl;
               }
            }
         }
      } else {
         istat = 0; // failure
         std::cout << "DISASTER - failed to generate RTop for move_std_residue\n";
         if (reference_residue) {
            //          molecule-class-info.cc:4184: passing `const mmdb::Residue' as `this'
            // argument of `int mmdb::Residue::GetSeqNum ()' discards qualifiers
            mmdb::Residue *tmp = (mmdb::Residue *) reference_residue;
            std::cout << "mainchain atoms missing from residue "
                      << tmp->GetSeqNum()
                      << tmp->GetChainID() << std::endl;
         } else {
            std::cout << "This should not happen!" << std::endl;
            std::cout << "null residue in move_std_residue" << std::endl;
         }
      }
   }
   return istat;
}


void
molecule_class_info_t::make_backup_from_outside() {  // when we have a multi mutate, we
                                    // want the wrapper to make a
                                    // backup when we start and set
                                    // changes when when finish.
                                    // Rather crap that this needs to
                                    // be done externally, I think.

   make_backup("from the outside");
}


void
molecule_class_info_t::set_have_unsaved_changes_from_outside() {

   have_unsaved_changes_flag = 1;

}


//
// Get a deep copy:
// return NULL on failure
//
mmdb::Residue *
molecule_class_info_t::get_standard_residue_instance(const std::string &residue_type) {

   graphics_info_t g;
   mmdb::Residue *std_residue = NULL;

   if (g.standard_residues_asc.read_success) {
      //      std::cout << "DEBUG:: There are " << g.standard_residues_asc.n_selected_atoms
      //                << " atoms in standard_residues_asc" << std::endl;
     int selHnd = g.standard_residues_asc.mol->NewSelection();
     g.standard_residues_asc.mol->Select ( selHnd,mmdb::STYPE_RESIDUE, 1, // .. TYPE, iModel
                                           "*", // Chain(s) it's "A" in this case.
                                           mmdb::ANY_RES,"*",  // starting res
                                           mmdb::ANY_RES,"*",  // ending res
                                           residue_type.c_str(),  // residue name
                                           "*",  // Residue must contain this atom name?
                                           "*",  // Residue must contain this Element?
                                           "*",  // altLocs
                                           mmdb::SKEY_NEW // selection key
                                           );
     // get the standard orientation residue for this residue type
     mmdb::PPResidue SelResidue;
     int nSelResidues;

     g.standard_residues_asc.mol->GetSelIndex(selHnd, SelResidue, nSelResidues);

     if (nSelResidues != 1) {
       std::cout << "This should never happen - ";
       std::cout << "badness in mci::get_standard_residue_instance(), we selected " << nSelResidues
                 << " residues looking for residues of type :" << residue_type << ":\n";
     } else {
       bool embed_in_chain_flag = true; // I think. Is this the one time where we *do* want embedding?
       std_residue = coot::deep_copy_this_residue_old_style(SelResidue[0], "", 1,
                                                            g.standard_residues_asc.UDDAtomIndexHandle, embed_in_chain_flag);
     }
     g.standard_residues_asc.mol->DeleteSelection(selHnd);
   }
   return std_residue;
}


// return the number of residues in chain with chain_id, return -1 on error
//
int
molecule_class_info_t::chain_n_residues(const char *chain_id) const {

   int r = -1;

   if (atom_sel.n_selected_atoms > 0) {
      int nchains = atom_sel.mol->GetNumberOfChains(1);
      for (int ichain=0; ichain<nchains; ichain++) {
         const mmdb::Chain *chain_p = (const mmdb::Chain *)atom_sel.mol->GetChain(1,ichain);
         std::string mol_chain_id(((mmdb::Chain*)chain_p)->GetChainID());
         if (mol_chain_id == std::string(chain_id)) {
            r = ((mmdb::Chain *)chain_p)->GetNumberOfResidues();
         }
      }
   }
   return r;
}


int
molecule_class_info_t::n_residues() const {

   int r = -1;
   if (atom_sel.n_selected_atoms > 0) {
      r = 0;
      for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
         mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
         mmdb::Chain *chain_p;
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            chain_p = model_p->GetChain(ichain);
            int nres = chain_p->GetNumberOfResidues();
            r += nres;
         }
      }
   }
   return r;
}

int
molecule_class_info_t::n_atoms() const {

   int r = -1;
   if (atom_sel.n_selected_atoms > 0) {
      r = 0;
      for(int imod = 1; imod<=atom_sel.mol->GetNumberOfModels(); imod++) {
         mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
         mmdb::Chain *chain_p;
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            chain_p = model_p->GetChain(ichain);
            int nres = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<nres; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               int n_atoms = residue_p->GetNumberOfAtoms();
               for (int iat=0; iat<n_atoms; iat++) {
                  mmdb::Atom *at = residue_p->GetAtom(iat);
                  if (! at->isTer())
                     if (! at->Het)
                        r++;
               }
            }
         }
      }
   }
   return r;
}



void
molecule_class_info_t::store_refmac_params(const std::string &mtz_filename,
                                           const std::string &fobs_col,
                                           const std::string &sigfobs_col,
                                           const std::string &r_free_col,
                                           int r_free_flag) {

   have_sensible_refmac_params = 1; // true
   refmac_mtz_filename = mtz_filename;
   refmac_fobs_col = fobs_col;
   refmac_sigfobs_col = sigfobs_col;
   refmac_r_free_col = r_free_col;
   refmac_r_free_flag_sensible = r_free_flag;

   if (r_free_flag) {
      // std::cout << "INFO:: Stored refmac parameters: " << refmac_fobs_col << " "  << refmac_sigfobs_col
      // << " " << refmac_r_free_col << " is sensible." << std::endl;
      logger.log(log_t::INFO, "Stored refmac parameters", refmac_fobs_col, refmac_sigfobs_col,
                 refmac_r_free_col, std::string("is sensible"));
   } else {
      // std::cout << "INFO:: Stored refmac parameters: " << refmac_fobs_col << " "  << refmac_sigfobs_col
      // << " the r-free-flag is not sensible" << std::endl;
      logger.log(log_t::INFO, "Stored refmac parameters", refmac_fobs_col, refmac_sigfobs_col,
                 refmac_r_free_col, std::string("is not sensible"));
   }
}

void
molecule_class_info_t::store_refmac_mtz_filename(const std::string &mtz_filename) {

   refmac_mtz_filename = mtz_filename;
}

void
molecule_class_info_t::store_refmac_phase_params(const std::string &phi,
                                                 const std::string &fom,
                                                 const std::string &hla,
                                                 const std::string &hlb,
                                                 const std::string &hlc,
                                                 const std::string &hld) {

   // std::cout << "BL DEBUG:: in store refmac pahse params" <<std::endl;
  have_refmac_phase_params = 1; // true
  refmac_phi_col = phi;
  refmac_fom_col = fom;
  refmac_hla_col = hla;
  refmac_hlb_col = hlb;
  refmac_hlc_col = hlc;
  refmac_hld_col = hld;
}


// return 0 on success
//
int
molecule_class_info_t::write_pdb_file(const std::string &filename) {

   int err = 1; // fail
   if (atom_sel.n_selected_atoms > 0) {
      std::string ext = coot::util::file_name_extension(filename);
      if (coot::util::extension_is_for_shelx_coords(ext)) {
         write_shelx_ins_file(filename);
      } else {
         mmdb::byte bz = mmdb::io::GZM_NONE;
         // err = write_atom_selection_file(atom_sel, filename, bz);
         err = coot::write_coords_pdb(atom_sel.mol, filename);
      }
   }
   return err;
}


// return 0 on success
//
int
molecule_class_info_t::write_cif_file(const std::string &filename) {

   int err = 1; // fail
   if (atom_sel.n_selected_atoms > 0) {
      mmdb::byte bz = mmdb::io::GZM_NONE;
      // err = write_atom_selection_file(atom_sel, filename, bz);
      err = coot::write_coords_cif(atom_sel.mol, filename);
   }
   return err;
}





// Add this molecule (typically of waters to this molecule by trying
// to put them into an already-existing solvent chain).  If a solvent
// chain does not already exist, put create a new chain id for the
// water_mol atoms.
//
// All the atoms of water_mol need to be in a
// chain that has a different chain id to all the chains in this
// molecule.  Else fail (return status 0).
//
int
molecule_class_info_t::insert_waters_into_molecule(const coot::minimol::molecule &water_mol, const std::string &res_name) {

   int istat = 0;  // set to failure initially

   // So run over the the chains of the existing molecule looking for
   // a solvent chain.  If there isn't one we simply use
   // append_to_molecule()
   //
   int nchains = atom_sel.mol->GetNumberOfChains(1);
   mmdb::Chain *chain_p = NULL;
   mmdb::Chain *solvent_chain_p = NULL;
   short int i_have_solvent_chain_flag = 0;
   for (int ichain=0; ichain<nchains; ichain++) {

      chain_p = atom_sel.mol->GetChain(1,ichain);
      if (chain_p->isSolventChain()) {
         solvent_chain_p = chain_p;
         std::string mol_chain_id(chain_p->GetChainID());
         i_have_solvent_chain_flag = 1;
      }
   }


   // For every atom in water_mol, create a new atom and a new residue
   // for it. Add the residue to our model's solvent chain and the
   // atom the the residue (of course).
   //
   if (i_have_solvent_chain_flag == 0) {

      // We didn't manage to find a solvent chain.
      // We need to create a new chain.
      chain_p = new mmdb::Chain;
      atom_sel.mol->GetModel(1)->AddChain(chain_p);
      std::pair<bool, std::string> u = unused_chain_id();
      if (u.first)
         chain_p->SetChainID(u.second.c_str());
      else
         chain_p->SetChainID("Z");
   } else {
      chain_p = solvent_chain_p; // put it back, (kludgey, should use
                                 // solvent_chain_p from here, not chain_p).
      // OK, we also need to remove any TER cards that are in that chain_p
      remove_TER_on_last_residue(solvent_chain_p);
   }

//    std::cout << "Debug:: choose chain " << chain_p->GetChainID()
//              << " with have_solvent flag: " << i_have_solvent_chain_flag
//              << std::endl;
//    std::cout << "Debug:: isSolvent for each residue of chain: " << std::endl;
//    for (int tmp_r=0; tmp_r<chain_p->GetNumberOfResidues(); tmp_r++) {
//       mmdb::Residue *rtmp = chain_p->GetResidue(tmp_r);
//       short int flag = isSolvent(rtmp->name);
//          std::cout << rtmp->name << " is solvent? " << flag << std::endl;
//    }

   std::pair<short int, int> p = coot::util::max_resno_in_chain(chain_p);
   float bf = graphics_info_t::default_new_atoms_b_factor; // 20.0 by default
   int max_resno;
   if (p.first) {
      max_resno = p.second;
   } else {
      max_resno = 0;
   }
   if (p.first || (i_have_solvent_chain_flag == 0)) {
      make_backup("insert_waters_into_molecule");
      std::cout << "INFO:: Adding to solvent chain: " << chain_p->GetChainID()
                << std::endl;
      int prev_max_resno = max_resno;
      mmdb::Residue *new_residue_p = NULL;
      mmdb::Atom    *new_atom_p = NULL;
      int water_count = 0;
      float occ = 1.0;
      if (is_from_shelx_ins_flag)
         occ = 11.0;
      for (unsigned int ifrag=0; ifrag<water_mol.fragments.size(); ifrag++) {
         for (int ires=water_mol[ifrag].min_res_no();
              ires<=water_mol[ifrag].max_residue_number();
              ires++) {
            for (unsigned int iatom=0; iatom<water_mol[ifrag][ires].atoms.size(); iatom++) {
               const coot::minimol::atom &atom = water_mol[ifrag][ires][iatom];
               new_residue_p = new mmdb::Residue;
               new_residue_p->SetResName(res_name.c_str());
               new_residue_p->seqNum = prev_max_resno + 1 + water_count;
               water_count++;
               bf = water_mol[ifrag][ires][iatom].temperature_factor;
               new_atom_p = new mmdb::Atom;
               new_atom_p->SetCoordinates(water_mol[ifrag][ires][iatom].pos.x(),
                                          water_mol[ifrag][ires][iatom].pos.y(),
                                          water_mol[ifrag][ires][iatom].pos.z(), occ, bf);
               if (false)
                  std::cout << "debug:: add water " << ires << " " << water_mol[ifrag][ires][iatom].pos.format()
                            << " with b " << bf << std::endl;
               new_atom_p->SetAtomName(water_mol[ifrag][ires][iatom].name.c_str());
               new_atom_p->Het = 1; // waters are now HETATMs
               strncpy(new_atom_p->element, water_mol[ifrag][ires][iatom].element.c_str(), 3);
               strncpy(new_atom_p->altLoc, water_mol[ifrag][ires][iatom].altLoc.c_str(), 2);

               // residue number, atom name, occ, coords, b factor

               // add the atom to the residue and the residue to the chain
               new_residue_p->AddAtom(new_atom_p);
               chain_p->AddResidue(new_residue_p);
            }
         }
      }
      atom_sel.mol->FinishStructEdit();
      update_molecule_after_additions(); // sets unsaved changes flag
      update_symmetry();
   }

   return istat;
}

// Add this molecule (typically of waters to this
// molecule... somehow).  All the atoms of water_mol need to be in a
// chain that has a different chain id to all the chains in this
// molecule.  Else fail (return status 0).
//
int
molecule_class_info_t::append_to_molecule(const coot::minimol::molecule &water_mol) {

   int istat = 0; // fail status initially.
   int n_atom = 0;  // 0 new atoms added initially.

   if (atom_sel.n_selected_atoms > 0) {

      make_backup("append_to_molecule()");

      // run over the chains in water_mol (there is only one for waters)
      //
      for (unsigned int ifrag=0; ifrag<water_mol.fragments.size(); ifrag++) {

//          std::cout << "DEBUG:: append_to_molecule: fragment id id for frag " << ifrag
//                    << " is " << water_mol[ifrag].fragment_id << std::endl;

         short int imatch = 0;

         // Run over chains of the existing mol, to see if there
         // already exists a chain with the same chain id as the
         // waters we want to add.  Only if imatch is 0 does this
         // function do anything.
         //
         int nchains = atom_sel.mol->GetNumberOfChains(1);
         mmdb::Chain *chain;
         for (int ichain=0; ichain<nchains; ichain++) {

            chain = atom_sel.mol->GetChain(1,ichain);
            std::string mol_chain_id(chain->GetChainID());

            if (water_mol.fragments[ifrag].fragment_id == mol_chain_id) {
               //
               imatch = 1;
               istat = 1;
               std::cout << "INFO:: Can't add waters from additional molecule "
                         << "chain id = " << mol_chain_id << std::endl
                         << "INFO:: That chain id already exists in this molecule"
                         << std::endl;
               break;
            }
         }

         mmdb::Model *model_p = atom_sel.mol->GetModel(1);
         if (imatch == 0) {
            // There was not already a chain in this molecule of that name.

            mmdb::Chain *new_chain_p;
            mmdb::Atom *new_atom_p;
            mmdb::Residue *new_residue_p;

            new_chain_p = new mmdb::Chain;
            std::cout << "DEBUG:: chain id of new chain :"
                      << water_mol[ifrag].fragment_id << ":" << std::endl;
            new_chain_p->SetChainID(water_mol[ifrag].fragment_id.c_str());
            model_p->AddChain(new_chain_p);

            for (int ires=water_mol[ifrag].min_res_no();
                 ires<=water_mol[ifrag].max_residue_number();
                 ires++) {

               if (water_mol[ifrag][ires].atoms.size() > 0) {
                  new_residue_p = new mmdb::Residue;
                  new_residue_p->seqNum = ires;
                  strcpy(new_residue_p->name, water_mol[ifrag][ires].name.c_str());
                  new_chain_p->AddResidue(new_residue_p);
                  for (unsigned int iatom=0; iatom<water_mol[ifrag][ires].atoms.size(); iatom++) {

                     new_atom_p = new mmdb::Atom;
                     new_atom_p->SetAtomName(water_mol[ifrag][ires][iatom].name.c_str());
                     new_atom_p->SetElementName(water_mol[ifrag][ires][iatom].element.c_str());
                     new_atom_p->SetCoordinates(water_mol[ifrag][ires][iatom].pos.x(),
                                                water_mol[ifrag][ires][iatom].pos.y(),
                                                water_mol[ifrag][ires][iatom].pos.z(),
                                                1.0, graphics_info_t::default_new_atoms_b_factor);
                     new_residue_p->AddAtom(new_atom_p);
                     n_atom++;
                  }
               }
            }
         }
      }

      std::cout << "INFO:: " << n_atom << " atoms added to molecule." << std::endl;
      if (n_atom > 0) {
         atom_sel.mol->FinishStructEdit();
         update_molecule_after_additions(); // sets unsaved changes flag
      }
   }

   return istat;
}


void
molecule_class_info_t::update_molecule_after_additions() {

   atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
   atom_sel = make_asc(atom_sel.mol); // does the udd stuff too.
   have_unsaved_changes_flag = 1;
   make_bonds_type_checked(__FUNCTION__);
}

std::string
molecule_class_info_t::Refmac_in_name() const {

   return Refmac_name_stub() + "-pre.pdb";

}

std::string
molecule_class_info_t::Refmac_out_name() const {

   return Refmac_name_stub() + ".pdb";

}

std::string
molecule_class_info_t::Refmac_mtz_out_name() const {

   return Refmac_name_stub() + ".mtz";

}

// combine and strip molecule and refmac count to come up with a pdb filename
// for refmac
std::string
molecule_class_info_t::Refmac_name_stub() const {

   // shall we try to take into account ccp4i refmac naming?
   // OK:
   // Here is an example:
   // demo.pdb         -> demo_refmac1.pdb
   // demo_refmac1.pdb -> demo_refmac2.pdb

   std::string refmac_name = "pre-refmac.pdb"; // default

   // First strip off the path of name_:
   std::string stripped_name;
   // /a/b.mtz -> b.mtz
#ifdef WINDOWS_MINGW
   std::string::size_type islash = coot::util::intelligent_debackslash(name_).find_last_of("/");
#else
   std::string::size_type islash = name_.find_last_of("/");
#endif // MINGW
   if (islash == std::string::npos) {
      // std::cout << "DEBUG:: slash not found in " << name_ << std::endl;
      stripped_name = name_;
   } else {
      // std::cout << "DEBUG:: slash found at " << islash << std::endl;
      // stripped_name = name_.substr(islash+1, name_.length());
      stripped_name = name_.substr(islash+1);
   }
   // std::cout << "DEBUG:: stripped_name: " << stripped_name << std::endl;


   std::string::size_type irefmac = stripped_name.rfind("-refmac");
   std::string::size_type irefmac_ccp4i = stripped_name.rfind("_refmac");

   if (irefmac == std::string::npos) { // not found

      // so was it a ccp4i refmac pdb file?

      if ( ! (irefmac_ccp4i == std::string::npos) ) {
         // it *was* a ccp4i pdb file:
         //
         refmac_name = stripped_name.substr(0,irefmac_ccp4i) + "_refmac";
         refmac_name += graphics_info_t::int_to_string(refmac_count);
      }
      // std::cout << "DEBUG:: irefmac not found in " << stripped_name  << std::endl;
      // lets strip off ".pdb", ".pdb.gz"
      std::string::size_type ipdb = stripped_name.rfind(".pdb");

      if (ipdb == std::string::npos) { // not a pdb

         // std::cout << "DEBUG:: ipdb not found" << std::endl;
         // just tack "refmac-2.pdb" on to the name then
         refmac_name = stripped_name + "_refmac";
         refmac_name += graphics_info_t::int_to_string(refmac_count);

      } else {
         // is a pdb:

         // std::cout << "DEBUG:: ipdb *was* found" << std::endl;
         refmac_name = stripped_name.substr(0,ipdb) + "_refmac";
         refmac_name += graphics_info_t::int_to_string(refmac_count);
      }
   } else {

      // refmac *was* found as part of the name
      // std::cout << "DEBUG:: irefmac *was* found in " << stripped_name << std::endl;
      refmac_name = stripped_name.substr(0,irefmac) + "_refmac";
      refmac_name += graphics_info_t::int_to_string(refmac_count);
   }

   // std::cout << "DEBUG:: returning refmac_name: " << refmac_name << std::endl;
   return refmac_name;

}


std::string
molecule_class_info_t::name_sans_extension(short int include_path_flag) const {

   std::string outstring = name_;

   std::string::size_type ipdb = name_.rfind(".pdb");
   if (ipdb != std::string::npos)
      outstring = name_.substr(0, ipdb);
   if (include_path_flag != 1) {
#ifdef WINDOWS_MINGW
     outstring = coot::util::intelligent_debackslash(outstring);
#endif // MINGW

     std::string::size_type islash = outstring.rfind("/");
     if (islash != std::string::npos)
       outstring = outstring.substr(islash+1);
   }

   return outstring;
}

void
molecule_class_info_t::update_molecule_to(std::vector<coot::scored_skel_coord> &pos_position) {

   std::cout << "DEBUG:: molecule_class_info_t update_molecule_to() with " << pos_position.size()
             << " skeleton positions" << std::endl;

   if (has_model()) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(1);

      if (! model_p) {
         std::cout << "ERROR:: Disaster in finding model_p in update_molecule_to"
                   << std::endl;
      } else {
         mmdb::Chain *chain_p;
         int n_chains = atom_sel.mol->GetNumberOfChains(1);
         for (int i_chain=0; i_chain<n_chains; i_chain++) {
            model_p->DeleteChain(i_chain);
         }

         // Now add new chain to model:
         chain_p = new mmdb::Chain;
         if (chain_p) {
            model_p->AddChain(chain_p);
            add_multiple_dummies(chain_p, pos_position);
         } else {
            std::cout << "ERROR:: creating chain in mol::update_molecule_to" << std::endl;
         }
      }
   } else {
      std::cout << "WARNING:: strange! This is not a valid model molecule. " << std::endl;
   }
}

// function callable from graphics_info_t::create_molecule_and_display()
//
void
molecule_class_info_t::add_multiple_dummies(const std::vector<coot::scored_skel_coord> &pos_position) {

   if (has_model()) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(1);
      if (model_p) {
         int n_chains = atom_sel.mol->GetNumberOfChains(1);
         if (n_chains > 0) {
            mmdb::Chain *chain_p = model_p->GetChain(0);
            add_multiple_dummies(chain_p, pos_position);
         }
      }
   }
}


// we presume that the chain exists.  This exists so that we dont do a
// backup every time we add a dummy atom (as is done using
// add_dummy_atom().
//
void
molecule_class_info_t::add_multiple_dummies(mmdb::Chain *chain_p,
                                            const std::vector<coot::scored_skel_coord> &pos_position) {


   // 20150803-PE FIXME - pass Geom_p().
   graphics_info_t g;
   coot::protein_geometry *geom_p = g.Geom_p();

   if (pos_position.size() > 0) {
      make_backup("add_multiple_dummies"); // maybe
   }

   for (unsigned int i=0; i<pos_position.size(); i++) {
      mmdb::Residue *res_p = new mmdb::Residue;
      mmdb::Atom *atom_p = new mmdb::Atom;
      chain_p->AddResidue(res_p);
      atom_p->SetAtomName(" DUM");
      atom_p->SetCoordinates(pos_position[i].position.x(),
                             pos_position[i].position.y(),
                             pos_position[i].position.z(), 1.0,
                             graphics_info_t::default_new_atoms_b_factor);

      atom_p->SetElementName(" O");
      res_p->AddAtom(atom_p);
      res_p->seqNum = i + 1;
      res_p->SetResName("DUM");

      // std::cout << atom_p << " added to molecule" << std::endl;
   }


   if (false)
      std::cout << "DEBUG:: add_multiple_dummies finishing.. "
                << pos_position.size() << std::endl;

   // Actually, we want to run this code when there are no new guide
   // points too.  This sets atom_sel.SelectionHandle properly, which
   // is needed in close_yourself, where a DeleteSelection() is done
   // to give back the memory.

   atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
   atom_sel.mol->FinishStructEdit();
   atom_sel = make_asc(atom_sel.mol);
   have_unsaved_changes_flag = 1;
   makebonds(0.0, 0.0, geom_p);
}

void
molecule_class_info_t::add_multiple_dummies(const std::vector<coot::Cartesian> &pos_position) {

   // 20150803-PE FIXME - pass Geom_p().
   graphics_info_t g;
   coot::protein_geometry *geom_p = g.Geom_p();

   if (atom_sel.mol) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(1);
      int n_chains = atom_sel.mol->GetNumberOfChains(1);
      if (n_chains > 0) {
         mmdb::Chain *chain_p = model_p->GetChain(0);
         if (pos_position.size() > 0) {
            make_backup("add_multiple_dummies"); // maybe

            for (unsigned int i=0; i< pos_position.size(); i++) {
               mmdb::Residue *res_p = new mmdb::Residue;
               mmdb::Atom *atom_p = new mmdb::Atom;
               chain_p->AddResidue(res_p);
               atom_p->SetAtomName(" DUM");
               atom_p->SetCoordinates(pos_position[i].x(),
                                      pos_position[i].y(),
                                      pos_position[i].z(), 1.0,
                                      graphics_info_t::default_new_atoms_b_factor);

               atom_p->SetElementName(" O");
               res_p->AddAtom(atom_p);
               res_p->seqNum = i + 1;
               res_p->SetResName("DUM");
            }
            atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
            atom_sel.mol->FinishStructEdit();
            atom_sel = make_asc(atom_sel.mol);
            have_unsaved_changes_flag = 1;
            makebonds(0.0, 0.0, geom_p);
         }
      }
   }
}


// return an empty vector on failure, a vector of size 6 on success:
//
std::pair<std::vector<float>, std::string>
molecule_class_info_t::get_cell_and_symm() const {

   std::pair<std::vector<float>, std::string> cell_spgr;

   mmdb::mat44 my_matt;
   if (atom_sel.mol) {
      int err = atom_sel.mol->GetTMatrix(my_matt, 0, 0, 0, 0);
      if (err != 0) {
         std::cout << "!! Warning:: No symmetry available for this template molecule"
                   << std::endl;
      } else {
         mmdb::realtype a[6];
         mmdb::realtype vol;
         int orthcode;
         atom_sel.mol->GetCell(a[0], a[1], a[2], a[3], a[4], a[5], vol, orthcode);
         for (int i=0; i<6; i++) cell_spgr.first.push_back(a[i]);
         cell_spgr.second = std::string(atom_sel.mol->GetSpaceGroup());
      }
   }
   return cell_spgr;
}

void
molecule_class_info_t::set_mmdb_cell_and_symm(std::pair<std::vector<float>, std::string> cell_spgr) {

   if (cell_spgr.first.size() == 6) {
      std::vector<float> a = cell_spgr.first; // short name
      atom_sel.mol->SetCell(a[0], a[1], a[2], a[3], a[4], a[5]);
      atom_sel.mol->SetSpaceGroup(cell_spgr.second.c_str());
      std::cout << "successfully set cell and symmetry" << std::endl;
   } else {
      std::cout << "WARNING:: failure to set cell on this molecule" << std::endl;
   }
}

bool
molecule_class_info_t::set_mmdb_symm(const std::string &spg) {

   atom_sel.mol->SetSpaceGroup(spg.c_str());
   std::string new_sg;
   const char *new_sg_chars = atom_sel.mol->GetSpaceGroup();
   if (new_sg_chars)
      new_sg = new_sg_chars;
   return (new_sg == spg);
}

// Return atom_index of -1 when no nearest atom.
//
std::pair<float, int>
molecule_class_info_t::nearest_atom(const coot::Cartesian &pos) const {

   float min_dist = 999999999.9;
   float d;
   int atom_index = -1;

   for(int i=0; i<atom_sel.n_selected_atoms; i++) {
      coot::Cartesian a(atom_sel.atom_selection[i]->x, atom_sel.atom_selection[i]->y, atom_sel.atom_selection[i]->z);
      d =  fabs((pos-a).length());
      if (d < min_dist) {
         min_dist = d;
         atom_index = i;
      }
   }

   std::pair<float, int> r;
   r.first = min_dist;
   r.second = atom_index;
   return r;
}

// return an empty
// vector
// if closed or is a coords
// mol, 3 elements and 6
// elements for a
// difference map
#ifndef EMSCRIPTEN
std::pair<GdkRGBA, GdkRGBA>
molecule_class_info_t::get_map_colours() const {

   return std::pair<GdkRGBA, GdkRGBA> (map_colour, map_colour_negative_level);
}
#else
std::pair<coot::colour_holder, coot::colour_holder>
molecule_class_info_t::get_map_colours() const {

   return std::pair<coot::colour_holder, coot::colour_holder> (map_colour, map_colour_negative_level);
}
#endif

// perhaps there is a better place for this?
//
std::vector<std::string>
molecule_class_info_t::set_map_colour_strings() const {

   // return something like
   // (list "set_last_map_colour" "0.2" "0.3" "0.4")


   std::vector<std::string> r;

   r.push_back("coot.set-last-map-colour");
   r.push_back(graphics_info_t::float_to_string(map_colour.red));
   r.push_back(graphics_info_t::float_to_string(map_colour.green));
   r.push_back(graphics_info_t::float_to_string(map_colour.blue));

   return r;
}


// and symm labels.
void
molecule_class_info_t::remove_atom_labels() {

   labelled_atom_index_list.clear();
   labelled_symm_atom_index_list.clear();
}

coot::minimol::molecule
molecule_class_info_t::eigen_flip_residue(const std::string &chain_id, int resno) {


   coot::minimol::molecule m;

   mmdb::Residue *res = get_residue(chain_id, resno, "");
   if (!res) {
      std::cout << "DEBUG:: residue not found " << chain_id << " " << resno
                << " in molecule number " << MoleculeNumber()
                << std::endl;
   } else {
      // make_backup();

      coot::ligand lig;
      coot::minimol::residue r(res);
      coot::minimol::fragment f(res->GetChainID());
      f.residues.push_back(coot::minimol::residue(res));
      coot::minimol::molecule ligand;
      ligand.fragments.push_back(f);

      ligand_flip_number++;
      if (ligand_flip_number == 4)
         ligand_flip_number = 0;

      lig.install_ligand(ligand);
      m = lig.flip_ligand(ligand_flip_number);

      make_bonds_type_checked(__FUNCTION__);
      have_unsaved_changes_flag = 1;

      replace_coords(make_asc(m.pcmmdbmanager()), 0, 1);
   }
   return m;
}

// return non blank string on a problem (so that it can be written to the status bar)
//
std::string
molecule_class_info_t::jed_flip(coot::residue_spec_t &spec,
                                const std::string &atom_name,
                                const std::string &alt_conf,
                                bool invert_selection,
                                coot::protein_geometry *geom) {

   // This function was copied to coot-utils - don't edit this, edit the coot-utils
   // version and call it from here - possibly delete this.
   //
   // But this does have_unsaved_changes_flag and make_backup.

   std::string problem_string;

   mmdb::Residue *residue = get_residue(spec);
   if (! residue) {
      std::cout << "WARNING:: no residue " << spec << " found in molecule" << std::endl;
   } else {

      // Does atom atom_name with given alt_conf exist in this residue?
      mmdb::Atom *clicked_atom = 0;
      int clicked_atom_idx = -1;
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      residue->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         std::string an(residue_atoms[iat]->name);
         if (an == atom_name) {
            std::string ac(residue_atoms[iat]->altLoc);
            if (ac == alt_conf) {
               clicked_atom = residue_atoms[iat];
               clicked_atom_idx = iat;
               break;
            }
         }
      }

      if (! clicked_atom) {
         std::cout << "WARNING:: atom \"" << atom_name << "\" not found in residue " << std::endl;
      } else {
         std::string monomer_type = residue->GetResName();

         std::pair<bool, coot::dictionary_residue_restraints_t> p =
            geom->get_monomer_restraints(monomer_type, imol_no);

         if (! p.first) {
            std::cout << "WARNING residue type " << monomer_type << " not found in dictionary" << std::endl;
         } else {
            bool iht = false;  // include_hydrogen_torsions_flag
            std::vector<coot::dict_torsion_restraint_t> all_torsions = p.second.get_non_const_torsions(iht);

            if (all_torsions.size() == 0) {
               problem_string = "There are no non-CONST torsions for this residue type";
            } else {
               std::vector<std::vector<std::string> > ring_atoms_sets = p.second.get_ligand_ring_list();
               std::vector<coot::dict_torsion_restraint_t> interesting_torsions;
               for (unsigned int it=0; it<all_torsions.size(); it++) {

                  bool is_ring_torsion_flag = all_torsions[it].is_ring_torsion(ring_atoms_sets);
                  if (! all_torsions[it].is_ring_torsion(ring_atoms_sets)) {
                     if (all_torsions[it].atom_id_2_4c() == atom_name)
                        interesting_torsions.push_back(all_torsions[it]);
                     if (all_torsions[it].atom_id_3_4c() == atom_name)
                        interesting_torsions.push_back(all_torsions[it]);
                  }
               }

               if (interesting_torsions.size() == 0) {
                  problem_string = "There are no non-CONST non-ring torsions for this atom";
               } else {

                  // make a constructor?
                  atom_selection_container_t residue_asc;
                  residue_asc.n_selected_atoms = n_residue_atoms;
                  residue_asc.atom_selection = residue_atoms;
                  residue_asc.mol = 0;

                  coot::contact_info contact = coot::getcontacts(residue_asc, monomer_type, imol_no, geom);
                  std::vector<std::vector<int> > contact_indices =
                     contact.get_contact_indices_with_reverse_contacts();

                  try {
                     coot::atom_tree_t tree(contact_indices, clicked_atom_idx, residue, alt_conf);
                     problem_string = jed_flip_internal(tree, interesting_torsions,
                                                        atom_name, clicked_atom_idx, invert_selection);
                     atom_sel.mol->FinishStructEdit();
                  }
                  catch (const std::runtime_error &rte) {
                     std::cout << "RUNTIME ERROR:: " << rte.what() << " - giving up" << std::endl;
                  }
                  make_bonds_type_checked(__FUNCTION__);
               }
            }
         }
      }
   }
   return problem_string;
}


// return a non-empty string on a problem
//
std::string
molecule_class_info_t::jed_flip_internal(coot::atom_tree_t &tree,
                                         const std::vector<coot::dict_torsion_restraint_t> &interesting_torsions,
                                         const std::string &atom_name,
                                         int atom_idx,
                                         bool invert_selection) {

   // This function was copied to coot-utils - don't edit this, edit the coot-utils
   // version and call it from here - possibly delete this.
   //
   // But this does have_unsaved_changes_flag and make_backup.

   std::string problem_string;
   unsigned int selected_idx = 0;

   if (interesting_torsions.size() > 0) {

      unsigned int best_fragment_size = 9999;
      if (interesting_torsions.size() > 1) {
         // select the best torsion based on fragment size.
         for (unsigned int i=0; i<interesting_torsions.size(); i++) {
            std::string atn_1 = interesting_torsions[i].atom_id_2_4c();
            std::string atn_2 = interesting_torsions[i].atom_id_3_4c();
            bool reverse = false; // dummy value

            std::pair<unsigned int, unsigned int> p = tree.fragment_sizes(atn_1, atn_2, reverse);
            if (p.first < best_fragment_size) {
               best_fragment_size = p.first;
               selected_idx = i;
            }
            if (p.second < best_fragment_size) {
               best_fragment_size = p.second;
               selected_idx = i;
            }
         }
      }

      problem_string = jed_flip_internal(tree, interesting_torsions[selected_idx],
                                         atom_name, atom_idx, invert_selection);
   }
   return problem_string;
}


// return a non-null string on a problem
//
std::string
molecule_class_info_t::jed_flip_internal(coot::atom_tree_t &tree,
                                         const coot::dict_torsion_restraint_t &torsion,
                                         const std::string &atom_name,
                                         int clicked_atom_idx,
                                         bool invert_selection) {

   // This function was copied to coot-utils - don't edit this, edit the coot-utils
   // version and call it from here - possibly delete this.
   //
   // But this does have_unsaved_changes_flag and make_backup.

   std::string problem_string;

   make_backup("jed_flip_internal");

   bool reverse = false; // reverse the moving dog<->tail fragment?

   if (invert_selection)
      reverse = true;

   std::string atn_1 = torsion.atom_id_2_4c();
   std::string atn_2 = torsion.atom_id_3_4c();

   if (torsion.atom_id_3_4c() == atom_name) {
      atn_1 = torsion.atom_id_3_4c();
      atn_2 = torsion.atom_id_2_4c();
   }

   int period = torsion.periodicity();

   if (period > 1) {

      double angle = 360/double(period);
      std::pair<unsigned int, unsigned int> p = tree.fragment_sizes(atn_1, atn_2, false);

      if (false) {  // debug
         std::cout << "flip this torsion: " << torsion << std::endl;
         std::cout << "DEBUG:: jed_flip_internal() fragment sizes: " << p.first << " " << p.second
                   << std::endl;
      }

      if (p.first > p.second)
         reverse = !reverse;

      tree.rotate_about(atn_1, atn_2, angle, reverse);
      have_unsaved_changes_flag = 1;
   } else {
      problem_string = "Selected torsion had a periodicity of ";
      problem_string += clipper::String(period);
   }
   return problem_string;
}



// So that we can move around all the atoms of a ligand (typically)
//
void
molecule_class_info_t::translate_by(float x, float y, float z) {

   if (has_model()) {
      make_backup("translate_by");
      for (int i=0; i<atom_sel.n_selected_atoms; i++) {
         atom_sel.atom_selection[i]->x += x;
         atom_sel.atom_selection[i]->y += y;
         atom_sel.atom_selection[i]->z += z;
      }
      make_bonds_type_checked(__FUNCTION__);
      have_unsaved_changes_flag = 1;
   }
}

void
molecule_class_info_t::translate_by_internal(const clipper::Coord_orth &co, mmdb::Residue *residue_p) {

   if (residue_p) {
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         at->x += co.x();
         at->y += co.y();
         at->z += co.z();
      }
   }
}


// Sets coot_save_index maybe (using set_coot_save_index()).
//
std::string
molecule_class_info_t::stripped_save_name_suggestion() {

   std::string s;

   std::string stripped_name1;
#ifdef WINDOWS_MINGW
   std::string::size_type islash = coot::util::intelligent_debackslash(name_).find_last_of("/");
#else
   std::string::size_type islash = name_.find_last_of("/");
#endif // MINGW
   if (islash == std::string::npos) {
      stripped_name1 = name_;
   } else {
      stripped_name1 = name_.substr(islash+1, name_.length());
   }
   // so we have got rid of the pathname.
   // now lets get rid of the extension
   //
   std::string::size_type ibrk   = stripped_name1.rfind(".brk");
   std::string::size_type ibrkgz = stripped_name1.rfind(".brk.gz");
   std::string::size_type ipdb   = stripped_name1.rfind(".pdb");
   std::string::size_type ires   = stripped_name1.rfind(".res");
   std::string::size_type ipdbgz = stripped_name1.rfind(".pdb.gz");
   std::string::size_type icoot  = stripped_name1.rfind("-coot-");

   std::string stripped_name2;
   if (icoot == std::string::npos) {
      if (ibrk == std::string::npos) {
         if (ibrkgz == std::string::npos) {
            if (ipdb == std::string::npos) {
               if (ires == std::string::npos) {
                  if (ipdbgz == std::string::npos) {
                     stripped_name2 = stripped_name1;
                  } else {
                     stripped_name2 = stripped_name1.substr(0, ipdbgz);
                  }
               } else {
                  stripped_name2 = stripped_name1.substr(0, ires);
               }
            } else {
               stripped_name2 = stripped_name1.substr(0, ipdb);
            }
         } else {
            stripped_name2 = stripped_name1.substr(0, ibrkgz);
         }
      } else {
         stripped_name2 = stripped_name1.substr(0, ibrk);
      }
   } else {
      set_coot_save_index(stripped_name1.substr(icoot));
      stripped_name2 = stripped_name1.substr(0,icoot);
   }

   stripped_name2 += "-coot-";
   stripped_name2 += graphics_info_t::int_to_string(coot_save_index);
   // As per George Sheldrick's suggestion, if this was from shelx,
   // suggest a .ins extension, not .pdb
   if (!is_from_shelx_ins_flag) {
      if (input_molecule_was_in_mmcif)
         stripped_name2 += ".cif";
      else
         stripped_name2 += ".pdb";
   } else {
      stripped_name2 += ".ins";
   }

//    std::cout << "DEBUG:: stripped_save_name_suggestion: "
//              << stripped_name2 << std::endl;

   return stripped_name2;
}

int
molecule_class_info_t::set_coot_save_index(const std::string &filename) {

   // std::cout << "extracting from :" << filename << std::endl;

   // filename is something like: "-coot-12.pdb".
   //
   // We want to find 12 and set coot_save_index to 12 + 1, which gets
   // used to suggest the next saved filename.
   //

   std::string twelve_pdb = filename.substr(6);
   // std::cout << "twelve_pdb:"<< twelve_pdb << std::endl;

   std::string::size_type ipdb   = twelve_pdb.rfind(".pdb");
   if (ipdb != std::string::npos) {
      // .pdb was found
      std::string twelve = twelve_pdb.substr(0,ipdb);
      int i = atoi(twelve.c_str());
      // std::cout << "found i: " << i << std::endl;
      if (i >= 0 && i<100000)
         coot_save_index = i+1;
   }
   return coot_save_index;
}


void
molecule_class_info_t::transform_by(mmdb::mat44 mat) {

   if (has_model()) {
      clipper::Coord_orth co;
      clipper::Coord_orth trans_pos;
      make_backup("transform-by");
      clipper::Mat33<double> clipper_mat(mat[0][0], mat[0][1], mat[0][2],
                                         mat[1][0], mat[1][1], mat[1][2],
                                         mat[2][0], mat[2][1], mat[2][2]);
      clipper::Coord_orth cco(mat[0][3], mat[1][3], mat[2][3]);
      clipper::RTop_orth rtop(clipper_mat, cco);
      std::cout << "INFO:: coordinates transformed by orthogonal matrix: \n"
                << rtop.format() << std::endl;
      clipper::Rotation rtn( clipper_mat );
      clipper::Polar_ccp4 polar = rtn.polar_ccp4();
      clipper::Euler_ccp4 euler = rtn.euler_ccp4();
      std::cout << "  Rotation - polar (omega,phi,kappa)  " << clipper::Util::rad2d(polar.omega()) << " " << clipper::Util::rad2d(polar.phi()) << " " << clipper::Util::rad2d(polar.kappa()) << std::endl;
      std::cout << "  Rotation - euler (alpha,beta,gamma) " << clipper::Util::rad2d(euler.alpha()) << " " << clipper::Util::rad2d(euler.beta()) << " " << clipper::Util::rad2d(euler.gamma()) << std::endl;
      std::cout << "  Translation - Angstroms             " << cco.x() << " " << cco.y() << " " << cco.z() << " " << std::endl;
      for (int i=0; i<atom_sel.n_selected_atoms; i++) {
         // atom_sel.atom_selection[i]->Transform(mat); // doesn't compile!
         // Argh.  sigh.  Use clipper. c.f. graphics_info_t::fill_hybrid_atoms()
         co = clipper::Coord_orth(atom_sel.atom_selection[i]->x,
                                  atom_sel.atom_selection[i]->y,
                                  atom_sel.atom_selection[i]->z);
         trans_pos = co.transform(rtop);
         atom_sel.atom_selection[i]->x = trans_pos.x();
         atom_sel.atom_selection[i]->y = trans_pos.y();
         atom_sel.atom_selection[i]->z = trans_pos.z();
      }
      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked(__FUNCTION__);
   }

}


void
molecule_class_info_t::transform_by(const clipper::RTop_orth &rtop) {

   make_backup("transform-by-clipper-rtop");
   std::cout << "INFO:: coordinates transformed by orthogonal matrix: \n"
             << rtop.format() << std::endl;
   if (have_unit_cell) {

      mmdb::realtype cell_params[6];
      mmdb::realtype vol;
      int orthcode;
      atom_sel.mol->GetCell(cell_params[0], cell_params[1], cell_params[2],
                            cell_params[3], cell_params[4], cell_params[5],
                            vol, orthcode);

      clipper::Cell cell(clipper::Cell_descr(cell_params[0],
                                             cell_params[1],
                                             cell_params[2],
                                             clipper::Util::d2rad(cell_params[3]),
                                             clipper::Util::d2rad(cell_params[4]),
                                             clipper::Util::d2rad(cell_params[5])));
      std::cout << "INFO:: fractional coordinates matrix:" << std::endl;
      std::cout << rtop.rtop_frac(cell).format() << std::endl;
   } else {
      std::cout << "No unit cell for this molecule, hence no fractional matrix." << std::endl;
   }
   clipper::Coord_orth co;
   clipper::Coord_orth trans_pos;
   if (has_model()) {
      for (int i=0; i<atom_sel.n_selected_atoms; i++) {
         co = clipper::Coord_orth(atom_sel.atom_selection[i]->x,
                                  atom_sel.atom_selection[i]->y,
                                  atom_sel.atom_selection[i]->z);
         trans_pos = co.transform(rtop);
         atom_sel.atom_selection[i]->x = trans_pos.x();
         atom_sel.atom_selection[i]->y = trans_pos.y();
         atom_sel.atom_selection[i]->z = trans_pos.z();
      }
      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked(__FUNCTION__);
   }
}

void
molecule_class_info_t::transform_by(const clipper::RTop_orth &rtop, mmdb::Residue *residue_moving) {

   make_backup("transform-by-with-residue-moving");
   std::cout << "INFO:: coordinates transformed_by: \n"
             << rtop.format() << std::endl;
   if (has_model()) {
      transform_by_internal(rtop, residue_moving);
      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked(__FUNCTION__);
   }
}



void
molecule_class_info_t::transform_by_internal(const clipper::RTop_orth &rtop, mmdb::Residue *residue_moving) {

   if (has_model()) {
      clipper::Coord_orth co;
      clipper::Coord_orth trans_pos;
      mmdb::PPAtom residue_atoms;
      int n_residue_atoms;
      residue_moving->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iatom=0; iatom<n_residue_atoms; iatom++) {
         clipper::Coord_orth p(residue_atoms[iatom]->x,
                               residue_atoms[iatom]->y,
                               residue_atoms[iatom]->z);
         clipper::Coord_orth p2 = p.transform(rtop);
         residue_atoms[iatom]->x = p2.x();
         residue_atoms[iatom]->y = p2.y();
         residue_atoms[iatom]->z = p2.z();
      }
   }
}


void
molecule_class_info_t::transform_zone_by(const std::string &chain_id, int resno_start, int resno_end,
                                         const std::string &ins_code,
                                         const clipper::RTop_orth &rtop,
                                         bool make_backup_flag) {

   if (make_backup_flag)
      make_backup("transform_zone_by");

   bool transformed_something = 0;
   if (resno_end < resno_start)
      std::swap(resno_end, resno_start);

   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   mmdb::Chain *chain_p;
   // run over chains of the existing mol
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      if (chain_id == chain_p->GetChainID()) {
         int nres = chain_p->GetNumberOfResidues();
         mmdb::Residue *residue_p;
         mmdb::Atom *at;
         for (int ires=0; ires<nres; ires++) {
            residue_p = chain_p->GetResidue(ires);
            int this_resno = residue_p->GetSeqNum();
            std::string this_ins_code = residue_p->GetInsCode();
            if ((this_resno >= resno_start) &&
                (this_resno <= resno_end)) {
               if (this_ins_code == ins_code) {
                  transform_by_internal(rtop, residue_p);
                  transformed_something = 1;
               }
            }
         }
      }
   }

   if (transformed_something) {
      atom_sel.mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      atom_sel.mol->FinishStructEdit();
      have_unsaved_changes_flag = 1;
      make_bonds_type_checked();
   }
}


float
molecule_class_info_t::get_contour_level_by_sigma() const {
   float cl = get_contour_level();
   return cl/map_sigma_;
}


void
molecule_class_info_t::set_contour_level(float f) {

   if (has_xmap()  || has_nxmap()) {
      contour_level = f;
      update_map(true);
   }
}

void
molecule_class_info_t::set_contour_level_by_sigma(float f) {
   if (has_xmap() || has_nxmap()) {
      contour_level = f * map_sigma_;
      update_map(true);
   }
}

std::vector <std::string>
molecule_class_info_t::get_map_contour_strings() const {

   std::vector <std::string> s;
   s.push_back("coot.set-last-map-contour-level");
   char cs[100];
   snprintf(cs, 99, "%e", contour_level);
   s.push_back(cs);

   return s;
}

std::vector <std::string>
molecule_class_info_t::get_map_contour_sigma_step_strings() const {

   std::vector <std::string> s;
   s.push_back("coot.set-last-map-sigma-step");
   s.push_back(graphics_info_t::float_to_string(contour_sigma_step));

//    s.push_back("set_contour_by_sigma_step_by_mol");
//    s.push_back(coot::util::float_to_string(contour_sigma_step));
//    s.push_back(coot::util::int_to_string(1));
//    s.push_back(coot::util::(imol_no));

   return s;
}

short int
molecule_class_info_t::contoured_by_sigma_p() const {
   return contour_by_sigma_flag;
}



mmdb::Chain *
molecule_class_info_t::water_chain() const {

   mmdb::Chain *water_chain = 0;

   if (has_model()) {

      mmdb::Model *model_p = atom_sel.mol->GetModel(1);
      mmdb::Residue *residue_p;
      mmdb::Chain *chain_p;

      if (model_p) {

         if (is_from_shelx_ins_flag) {
            water_chain = water_chain_from_shelx_ins();
         } else {
            int nchains = model_p->GetNumberOfChains();
            for (int ich=0; ich<nchains; ich++) {
               chain_p = model_p->GetChain(ich);
               int nres = chain_p->GetNumberOfResidues();
               short int all_water_flag = 1;
               for (int ires=0; ires<nres; ires++) {
                  residue_p = chain_p->GetResidue(ires);
                  std::string resname(residue_p->name);
                  if (! ( (resname == "WAT") || (resname == "HOH"))) {
                     all_water_flag = 0;
                     break;
                  }
               }
               if (all_water_flag) {
                  water_chain = chain_p;
                  break;
               }
            }
         }
      }
   }
   return water_chain;
}


// there is only one chain from a shelxl ins file.
mmdb::Chain *
molecule_class_info_t::water_chain_from_shelx_ins() const {

   mmdb::Chain *water_chain = 0;
   mmdb::Model *model_p = atom_sel.mol->GetModel(1);

   if (has_model()) {
      int nchains = model_p->GetNumberOfChains();
      for (int ich=0; ich<nchains; ich++) {
         water_chain = model_p->GetChain(ich);
      }
   }
   return water_chain;
}

bool
molecule_class_info_t::is_het_residue(mmdb::Residue *residue_p) const {

   bool status = false;

   if (residue_p) {
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for(int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            if (at->Het) {
               status =  true;
               break;
            }
         }
      }
   }
   return status;
}



// return state, max_resno + 1, or 0, 1 of no residues in chain.
//
// new_res_no_by_hundreds is default false
std::pair<short int, int>
molecule_class_info_t::next_residue_number_in_chain(mmdb::Chain *w,
                                                    bool new_res_no_by_hundreds) const {

   std::pair<short int, int> p(0,1);
   int max_res_no = -9999;

   if (w) {
      int nres = w->GetNumberOfResidues();
      mmdb::Residue *residue_p;
      if (nres > 0) {
         for (int ires=nres-1; ires>=0; ires--) {
            residue_p = w->GetResidue(ires);
            if (residue_p->seqNum > max_res_no) {
               max_res_no = residue_p->seqNum;
               bool is_het_residue_flag = is_het_residue(residue_p);
               if (is_het_residue_flag) {
                  p = std::pair<short int, int>(1, residue_p->seqNum+1);
               } else {
                  if (new_res_no_by_hundreds) {
                     if (max_res_no < 9999) {
                        int res_no = coot::util::round_up_by_hundreds(max_res_no+1);
                        p = std::pair<short int, int>(1, res_no+1);
                     }
                  } else {
                     if (max_res_no < 9999) {
                        p = std::pair<short int, int>(1, max_res_no+1);
                     }
                  }
               }
            }
         }
         if (! p.first) {
            //  first the first space starting from the front
            int test_resno_start = 1001;
            bool is_clear = false;
            while (! is_clear) {
               is_clear = true;
               for (int iser=0; iser<nres; iser++) {
                  int resno_res = w->GetResidue(iser)->seqNum;
                  if (resno_res >= test_resno_start) {
                     if (resno_res <= (test_resno_start+10)) {
                        is_clear = false;
                     }
                  }
                  if (! is_clear)
                     break;
               }
               test_resno_start += 100;
            }
            p = std::pair<short int, int> (1, test_resno_start);
         }
      }
   }
   return p;
}



// add a factor to scale the colours in b factor representation:.
// It goes into the atom_sel.mol
void
molecule_class_info_t::set_b_factor_bonds_scale_factor(float f) {

   if (atom_sel.mol) {
      int udd_handle =
         atom_sel.mol->RegisterUDReal(mmdb::UDR_HIERARCHY,
                                      coot::b_factor_bonds_scale_handle_name.c_str());
      if (udd_handle > 0) {
         atom_sel.mol->PutUDData(udd_handle, f);

         // test getting the uddata:
         int udd_b_factor_handle =
            atom_sel.mol->GetUDDHandle(mmdb::UDR_HIERARCHY,
                                       coot::b_factor_bonds_scale_handle_name.c_str());
         if (udd_b_factor_handle > 0) {
            mmdb::realtype scale;
            if (atom_sel.mol->GetUDData(udd_b_factor_handle, scale) == mmdb::UDDATA_Ok) {
               //                std::cout << " test got b factor scale: " << scale << std::endl;
            } else {
                std::cout << "ERROR:: bad get b factor scale " << std::endl;
            }
         }
      }
   }
   make_bonds_type_checked(__FUNCTION__);
}

std::pair<bool, std::string>
molecule_class_info_t::chain_id_for_shelxl_residue_number(int shelxl_resno) const {

   int imod = 1;
   bool found_it = 0;
   std::string chain_id_unshelxed = "not-found";

   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   mmdb::Chain *chain_p;
   // run over chains of the existing mol
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      mmdb::PResidue residue_p;
      for (int ires=0; ires<nres; ires++) {
         residue_p = chain_p->GetResidue(ires);
         int resno = residue_p->GetSeqNum();
         if (resno == shelxl_resno) {
            chain_id_unshelxed = chain_p->GetChainID();
            found_it = 1;
         }

         if (found_it)
            break;
      }
      if (found_it)
         break;
   }

   return std::pair<bool, std::string> (found_it, chain_id_unshelxed);
}


// default arg debug_atoms_also_flag=false
void
molecule_class_info_t::debug(bool debug_atoms_also_flag) const {

   int imod = 1;

   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   mmdb::Chain *chain_p;
   // run over chains of the existing mol
   int nchains = model_p->GetNumberOfChains();
   std::cout << "debug:: debug(): model 1 has " << nchains << " chains" << std::endl;
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      mmdb::PResidue residue_p;
      for (int ires=0; ires<nres; ires++) {
         residue_p = chain_p->GetResidue(ires);
         if (residue_p) {
            std::cout << "debug():  " << residue_p->GetResName() << " "
                      << chain_p->GetChainID() << " " << residue_p->GetSeqNum()
                      << " \"" << residue_p->GetInsCode() << "\" index: "
                      << residue_p->index << std::endl;

            if (debug_atoms_also_flag) {
               mmdb::Atom **residue_atoms = 0;
               int n_residue_atoms;
               residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
               for (int iat=0; iat<n_residue_atoms; iat++) {
                  mmdb::Atom *at = residue_atoms[iat];
                  std::cout << "     " << std::setw(2) << iat << " " << coot::atom_spec_t(at)
                            << " " << at->x << " " << at->y << " " << at->z
                            << std::endl;
               }
            }
         }
      }
   }
}

void
molecule_class_info_t::clear_all_fixed_atoms() {

   std::cout << "m::clear_all_fixed_atoms() " << fixed_atom_specs.size() << std::endl;

   for (unsigned int i=0; i<fixed_atom_specs.size(); i++) {
      mark_atom_as_fixed(fixed_atom_specs[i], false);
   }
   fixed_atom_specs.clear();
   fixed_atom_positions.clear();
}

void
molecule_class_info_t::mark_atom_as_fixed(const coot::atom_spec_t &atom_spec, bool state) {

   std::cout << "--------------------- mci: mark_atom_as_fixed() " << atom_spec << " " <<  state << std::endl;

   if (has_model()) {
      int imod = 1;
      bool found = 0;
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         std::string chain_id_model = chain_p->GetChainID();
         if (atom_spec.chain_id == chain_id_model) {
            int nres = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<nres; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               int resno_model = residue_p->GetSeqNum();
               if (resno_model == atom_spec.res_no) {
                  int n_atoms = residue_p->GetNumberOfAtoms();

                  for (int iat=0; iat<n_atoms; iat++) {
                     mmdb::Atom *at = residue_p->GetAtom(iat);
                     if (atom_spec.matches_spec(at)) {
                        if (state) {
                           // try to get the atom index of this atom
                           // (at) and make it part of the spec.  So
                           // that we can use it again in
                           // update_fixed_atom_positions().  If the
                           // atom index is correct, we don't need to
                           // search the molecule for the spec.
                           coot::atom_spec_t atom_spec_local = atom_spec;
                           int idx = get_atom_index(at);
                           atom_spec_local.int_user_data = idx;
                           fixed_atom_specs.push_back(atom_spec_local);
                           std::cout << "INFO:: " << atom_spec << " marked as fixed"
                                     << std::endl;
                           found = 1;
                        } else {
                           //  try to remove at from marked list
                           if (fixed_atom_specs.size() > 0) {
                              std::vector<coot::atom_spec_t>::iterator it;
                              for (it=fixed_atom_specs.begin();
                                   it != fixed_atom_specs.end();
                                   ++it) {
                                 if (atom_spec == *it) {
                                    std::cout << "INFO:: removed " << atom_spec
                                              << " from fixed atom." << std::endl;
                                    fixed_atom_specs.erase(it);
                                    found = 1;
                                    break;
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
      if (found) {
         update_fixed_atom_positions();
      }
   }
}


int
molecule_class_info_t::move_waters_to_around_protein() {

   make_backup("move_waters_to_around_protein");
   int r = coot::util::move_waters_around_protein(atom_sel.mol);
   have_unsaved_changes_flag = 1;
   make_bonds_type_checked(__FUNCTION__);
   return r;
}

void
molecule_class_info_t::move_hetgroups_to_around_protein() {

   make_backup(__FUNCTION__);
   coot::util::move_hetgroups_around_protein(atom_sel.mol);
   have_unsaved_changes_flag = 1;
   make_bonds_type_checked(__FUNCTION__);
}




// Return the maximum minimum distance of waters to protein atoms.
// return something negative when we can't do above (no protein
// atoms or no water atoms).
float
molecule_class_info_t::max_water_distance() {

   // Do not account for alt confs of the water positions. That
   // sophistication can ome later if we use this function for
   // something other than testing a greg-test.

   float f = -1.0;

   std::vector<clipper::Coord_orth> protein_positions;
   std::vector<clipper::Coord_orth> water_positions;

   for (int iat=0; iat<atom_sel.n_selected_atoms; iat++) {
      clipper::Coord_orth pt(atom_sel.atom_selection[iat]->x,
                             atom_sel.atom_selection[iat]->y,
                             atom_sel.atom_selection[iat]->z);
      std::string res_name(atom_sel.atom_selection[iat]->GetResName());
      if (res_name == "HOH" || res_name == "WAT") {
         water_positions.push_back(pt);
      } else {
         protein_positions.push_back(pt);
      }
   }

   // now protein_positions and water_positions are filled.
   if (protein_positions.size() > 0) {
      if (water_positions.size() > 0) {
         double max_water_dist_2 = -1.0;
         for (unsigned int iw=0; iw<water_positions.size(); iw++) {
            double best_dist_2 = 999999999.9;
            for (unsigned int ip=0; ip<protein_positions.size(); ip++) {
               double d2 = (water_positions[iw]-protein_positions[ip]).lengthsq();
               if (d2 < best_dist_2) {
                  best_dist_2 = d2;
               }
            }
            if (best_dist_2 > max_water_dist_2) {
               max_water_dist_2 = best_dist_2;
            }
         }
         if (max_water_dist_2 > 0.0)
            f = sqrt(max_water_dist_2);
      }
   }
   return f;
}


// ---- utility function --- (so that we know to delete hydrogens
// from HETATM molecule before merging with this one)
//
bool
molecule_class_info_t::molecule_has_hydrogens() const {

   bool r = 0;
   for (int i=0; i<atom_sel.n_selected_atoms; i++) {
      std::string ele(atom_sel.atom_selection[i]->element);
      if (ele == " H") {
         r = 1;
         break;
      }
      if (ele == " D") {
         r = 1;
         break;
      }
   }
   return r;
}

// -------- simply print it (at the moment) --------------
void
molecule_class_info_t::print_secondary_structure_info() {

   int n_models = atom_sel.mol->GetNumberOfModels();
   for (int imod=1; imod<=n_models; imod++) {
      mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
      coot::util::print_secondary_structure_info(model_p);
   }
}


// static
int
molecule_class_info_t::watch_coordinates_file(gpointer data) {

   int status = 1;

   updating_coordinates_molecule_parameters_t *ucp = static_cast<updating_coordinates_molecule_parameters_t *>(data);
   const updating_coordinates_molecule_parameters_t &rucp = *ucp;

   std::cout << "DEBUG:: watching " << rucp.imol << " " << rucp.pdb_file_name << std::endl;
   status = graphics_info_t::molecules[rucp.imol].update_coordinates_molecule_if_changed(rucp);
   return status;
}

// bool continue_watching_coordinates_file;
// updating_coordinates_molecule_parameters_t updating_coordinates_molecule_previous;
int
molecule_class_info_t::update_coordinates_molecule_if_changed(const updating_coordinates_molecule_parameters_t &ucp_in) {

   int status = 1;
   if (continue_watching_coordinates_file) {
      bool update_it = false;

      updating_coordinates_molecule_parameters_t ucp = ucp_in;
      struct stat s;
      int status = stat(ucp.pdb_file_name.c_str(), &s);
      if (status != 0) {
         std::cout << "WARNING:: update_map_from_mtz_if_changed() Error reading "
                   << ucp.pdb_file_name << std::endl;
      } else {
         if (!S_ISREG (s.st_mode)) {
            std::cout << "WARNING:: update_map_from_mtz_if_changed() not a reguular file: "
                      << ucp.pdb_file_name << std::endl;
            continue_watching_coordinates_file = false;
         } else {
            // happy path
            ucp.update_from_stat_info(s);
         }
      }

      if (false)
         std::cout << "#### ctime comparision: imol << " << imol_no << " was "
            << updating_coordinates_molecule_previous.ctime.tv_sec << " " << updating_coordinates_molecule_previous.ctime.tv_nsec
            << " now " << ucp.ctime.tv_sec << " " << ucp.ctime.tv_nsec
            << std::endl;

      if (ucp.ctime.tv_sec > updating_coordinates_molecule_previous.ctime.tv_sec) {
         update_it = true;
      } else {
         if (ucp.ctime.tv_sec == updating_coordinates_molecule_previous.ctime.tv_sec) {
            if (ucp.ctime.tv_nsec > updating_coordinates_molecule_previous.ctime.tv_nsec) {
               update_it = true;
            }
         }
      }
      if (update_it) {

         std::string cwd = coot::util::current_working_dir();
         short int reset_rotation_centre = 0;
         short int is_undo_or_redo = 0;
         bool allow_duplseqnum = true;
         bool v2_convert_flag = false;

         handle_read_draw_molecule(imol_no, ucp.pdb_file_name, cwd,
                                   graphics_info_t::Geom_p(),
                                   reset_rotation_centre,
                                   is_undo_or_redo,
                                   allow_duplseqnum,
                                   v2_convert_flag,
                                   bond_width,
                                   Bonds_box_type(),
                                   false);
         updating_coordinates_molecule_previous = ucp;
         graphics_info_t::graphics_draw();
      }
   } else {
      status = 0;
   }
   return status;
}

// no redraw
void
molecule_class_info_t::update_self_from_file(const std::string &pdb_file_name) {

   std::string cwd = coot::util::current_working_dir();
   short int reset_rotation_centre = 0;
   short int is_undo_or_redo = 0;
   bool allow_duplseqnum = true;
   bool v2_convert_flag = false;

   handle_read_draw_molecule(imol_no, pdb_file_name, cwd,
                             graphics_info_t::Geom_p(),
                             reset_rotation_centre,
                             is_undo_or_redo,
                             allow_duplseqnum,
                             v2_convert_flag,
                             bond_width,
                             Bonds_box_type(),
                             false);

}

void
molecule_class_info_t::update_self(const coot::mtz_to_map_info_t &mmi) {

   // 20240702-PE make this part of mmi?
   float n_sd = 12.0f; // number of suggested standard deviations for the contour level of a newly-created map

   bool previous_map_was_sane = true;
   if (xmap.is_null()) previous_map_was_sane = false;

   std::cout << "############### --- start --- update_self() xmap is sane: " << previous_map_was_sane << std::endl;

   float sr = graphics_info_t::map_sampling_rate;
   std::string cwd = coot::util::current_working_dir(); // why is this needed?
   map_fill_from_mtz(mmi.mtz_file_name, cwd, mmi.f_col, mmi.phi_col, mmi.w_col, mmi.use_weights, mmi.is_difference_map, sr, true);

   bool ipz = true; // ignore_pseudo zeros
   mean_and_variance<float> mv = map_density_distribution(xmap, 20, false, ipz);

   if (! previous_map_was_sane) {
      contour_level = mv.mean + n_sd * std::sqrt(mv.variance);
      std::cout << "-------- new map contour level " << contour_level << std::endl;
      update_map_in_display_control_widget();
   } else {
      std::cout << "--------- using old map contour level " << contour_level << std::endl;
   }
   update_map_internal();

}

#include "coot-utils/diff-diff-map-peaks.hh"

// update this difference map if the coordinates of the other molecule have changed.
//
// static (!) - this is a callback function
int
molecule_class_info_t::watch_coordinates_updates(gpointer data) {

   // The bulk of this function should be in graphics_info_t.
   // This function should merely call that function

   int status = 1; // continue

   if (data) {
      updating_model_molecule_parameters_t *ummp_p = static_cast<updating_model_molecule_parameters_t *> (data);
      int imol_coords    = ummp_p->imol_coords;
      int imol_2fofc_map = ummp_p->imol_2fofc_map;
      int imol_diff_map  = ummp_p->imol_fofc_map;
      int imol_data = ummp_p->imol_map_with_data_attached;

      // std::cout << "DEBUG:: in watch_coordinates_updates() ummp: " << ummp_p->format() << std::endl;

      graphics_info_t g;
      if (g.is_valid_map_molecule(imol_data)) {
         if (g.is_valid_map_molecule(imol_diff_map)) {
            if (g.is_valid_model_molecule(imol_coords)) {
               int backup_index_current = g.molecules[imol_diff_map].get_other_molecule_backup_index();
               // Do we need to update the map?
               // We don't want to update the map the first time around
               int backup_index_for_molecule = g.molecules[imol_coords].get_history_index();
               if (false)
                  std::cout << "DEBUG:: watch_coordinates_updates() backup_index_current " << backup_index_current
                            << " backup_index_for_molecule " << backup_index_for_molecule << std::endl;
               if (backup_index_current != backup_index_for_molecule) {

                  // 20230430-PE We *do* want the maps to be updated at the first time.
                  // if (backup_index_current == -1) {
                  //    std::cout << "DEBUG:: watch_coordinates_updates() First time, do nothing " << std::endl;
                  // } else {
                  //    std::cout << "DEBUG:: watch_coordinates_updates() Update the map " << imol_map << std::endl;
                  //    g.sfcalc_genmap(imol_coords, imol_data, imol_map);
                  // }
                  // g.sfcalc_genmap(imol_coords, imol_data, imol_map);

                  clipper::Xmap<float> *xmap_2fofc_p = &g.molecules[imol_data].xmap;
                  clipper::Xmap<float> *xmap_fofc_p  = &g.molecules[imol_diff_map].xmap;

                  // save the old (current) difference map
                  g.molecules[imol_diff_map].updating_map_previous_difference_map = g.molecules[imol_diff_map].xmap;

                  coot::util::sfcalc_genmap_stats_t stats =
                     g.sfcalc_genmaps_using_bulk_solvent(imol_coords, imol_data, xmap_2fofc_p, xmap_fofc_p);

                  g.latest_sfcalc_stats = stats;

                  // ------------ gone diff map peaks ---------------
                  float base_level = 0.2;
                  clipper::Coord_orth screen_centre(g.X(), g.Y(), g.Z());
                  auto diff_diff_map_peaks = coot::diff_diff_map_peaks(g.molecules[imol_diff_map].updating_map_previous_difference_map,
                                                                       g.molecules[imol_diff_map].xmap, base_level);
                  clipper::Cell       cell       = g.molecules[imol_diff_map].xmap.cell();
                  clipper::Spacegroup spacegroup = g.molecules[imol_diff_map].xmap.spacegroup();
                  auto moved_peaks = coot::move_peaks_to_around_position(screen_centre, spacegroup, cell, diff_diff_map_peaks);
                  std::cout << "INFO:: moved peaks " << moved_peaks.size() << std::endl;
                  // the first one, where we shift from Refmac map to Clipper map has many thousands
                  // of peaks. So ignore that one.
                  if (moved_peaks.size() < 1000) {
                     std::vector<std::pair<glm::vec3, float> > positions(moved_peaks.size());
                     for (unsigned int i=0; i<moved_peaks.size(); i++) {
                        const auto &p = moved_peaks[i].first;
                        float f = moved_peaks[i].second;
                        positions[i] = std::make_pair(glm::vec3(p.x(), p.y(), p.z()), f);
                     }
                     g.setup_draw_for_particles_for_gone_diff_map_peaks(positions);
                  }
                  // ------------ done gone diff map peaks ---------------

                  // 20230501-PE tweak the contour levels so that the levels (number of lines in the mesh)
                  // are about the same
                  //
                  if (backup_index_current == -1) {
                     g.molecules[imol_data].set_contour_level(    1.3  *     g.molecules[imol_data].get_contour_level());
                     g.molecules[imol_diff_map].set_contour_level(0.55 * g.molecules[imol_diff_map].get_contour_level());
                  }

                  g.molecules[imol_diff_map].update_map(true);
                  g.molecules[imol_data].update_map(true);

                  g.calculate_new_rail_points(*ummp_p);
                  g.updating_maps_update_the_coot_points_overlay();

                  g.molecules[imol_diff_map].other_molecule_backup_index = backup_index_for_molecule;

                  g.graphics_draw(); // so that we can see the results of the recontouring.

               } else {
                  if (false)
                     std::cout << "DEBUG:: watch_coordinates_updates() No need for an update "
                               << backup_index_current << std::endl;
               }
            } else {
               std::cout << "ERROR:: bad model index in watch_coordinates_updates() " << imol_coords << std::endl;
            }
         } else {
            std::cout << "ERROR:: bad diff map index in watch_coordinates_updates() " << imol_diff_map << std::endl;
         }
      } else {
         std::cout << "ERROR:: bad 2fofc map index in watch_coordinates_updates() " << imol_2fofc_map << std::endl;
      }
   }
   return status;
}

// update this map if the coordinates of the other molecule have changed. Both maps use this function callback
// looking for changes in the model molecule
//
// static (!) - this is a callback function
int
molecule_class_info_t::updating_coordinates_updates_genmaps(gpointer data) {

   int status = 1; // continue

   bool debug = false;

   if (debug)
      std::cout << "------- updating_coordinates_updates_genmaps() called " << std::endl;

   if (data) {
      updating_model_molecule_parameters_t *ummp_p = static_cast<updating_model_molecule_parameters_t *> (data);
      int imol_coords = ummp_p->imol_coords;
      int imol_data = ummp_p->imol_map_with_data_attached;
      graphics_info_t g;
      if (is_valid_model_molecule(imol_coords)) {
         if (g.is_valid_map_molecule(ummp_p->imol_fofc_map)) {
            if (g.is_valid_model_molecule(imol_coords)) {
               int backup_index_current  = g.molecules[imol_coords].get_history_index();
               int backup_index_previous = g.molecules[imol_coords].previous_backup_index;
               if (debug)
                  std::cout << "DEBUG:: updating_coordinates_updates_genmaps() backup_index_current "
                            << backup_index_current << " backup_index_previous " << backup_index_previous
                            << std::endl;
               if (backup_index_current != backup_index_previous) {
                  if (is_valid_map_molecule(ummp_p->imol_2fofc_map)) {
                     if (is_valid_map_molecule(ummp_p->imol_fofc_map)) {
                        clipper::Xmap<float> *xmap_1_p = &g.molecules[ummp_p->imol_2fofc_map].xmap;
                        clipper::Xmap<float> *xmap_2_p = &g.molecules[ummp_p->imol_fofc_map ].xmap;
                        if (debug)
                           std::cout << ":::::::::::::::::: genmaps! " << imol_coords << " " << imol_data << " "
                                     << xmap_1_p << " " << xmap_2_p << std::endl;
                        g.sfcalc_genmaps_using_bulk_solvent(imol_coords, imol_data, xmap_1_p, xmap_2_p);
                        g.update_maps();
                        g.molecules[imol_coords].previous_backup_index = backup_index_current; // for next time
                        g.graphics_draw();
                     }
                  }
               } else {
                  if (debug)
                     std::cout << "DEBUG:: updating_coordinates_updates_genmaps() No need for an update " << backup_index_current << std::endl;
               }
            } else {
               status = 0;
            }
         } else {
            status = 0;
         }
      } else {
         status = 0;
      }
   }
   return status;
}



// Don't forget to call graphics_info_t::attach_buffers() before calling this function
void
molecule_class_info_t::add_ribbon_representation_with_user_defined_residue_colours(const std::vector<std::pair<unsigned int, coot::colour_holder> > &user_defined_colours,
                                                                                   const std::string &mesh_name) {

   int secondary_structure_usage_flag = CALC_SECONDARY_STRUCTURE; // pass this.

   if (true) {
      std::cout << "DEBUG:: in add_ribbon_represenation_with_user_defined_residue_colours....................." << std::endl;
      std::cout << "DEBUG:: in add_ribbon_represenation_with_user_defined_residue_colours user_defined_colours size "
                << user_defined_colours.size() << std::endl;
      for (size_t i = 0; i < user_defined_colours.size(); i++) {
         unsigned int idx = user_defined_colours[i].first;
         const auto &col  = user_defined_colours[i].second;
         std::cout << "DEBUG:: add_ribbon_represenation_with_user_defined_residue_colours() col-index: " << idx
                   << " col: " << col << std::endl;
      }
   }

   molecular_mesh_generator_t mmg;
   Material material;

   material.do_specularity = true;
   material.shininess = 256;
   material.specular_strength = 0.55;

   int imod = 1;
   mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int n_res = chain_p->GetNumberOfResidues();
         if (n_res > 1) {
            // the indexing into the user_defined_colours vector is in the UDD data of the residue
            std::pair<std::vector<s_generic_vertex>, std::vector<g_triangle> > verts_and_tris =
               mmg.get_molecular_triangles_mesh_for_ribbon_with_user_defined_residue_colours(atom_sel.mol, chain_p,
                                                                                             user_defined_colours,
                                                                                             secondary_structure_usage_flag,
                                                                                             M2T_float_params, M2T_int_params);
            Mesh mesh(verts_and_tris);
            mesh.set_name(mesh_name);
            meshes.push_back(mesh);
            meshes.back().setup(material);
         }
      }
   }

}
