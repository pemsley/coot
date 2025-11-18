/* src/pick.cc
 *
 * Copyright 2002, 2003, 2004 by The University of York
 * Author: Paul Emsley
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
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"

#include <string>
#include <vector> // for mmdb-crystal

#include <math.h>
#include <string.h>  // strncpy

#include <mmdb2/mmdb_manager.h>

#include "coords/cos-sin.h"
#include "coords/mmdb-extras.hh"
#include "coords/mmdb.hh"
#include "coords/mmdb-crystal.hh" //need for Bond_lines now
#include "coords/Cartesian.hh"
#include "coords/Bond_lines.hh"

#include "graphics-info.h"

// #if __APPLE__
// #   include <OpenGL/glu.h>
// #else
// #   include <GL/glu.h>
// #endif


#include "molecule-class-info.h"
#include "globjects.h"
#include "cc-interface.hh" // for status bar text
#include "c-interface-generic-objects.h" // for tomo_pick

// this should be in graphics_info_t also.
pick_info
pick_atom_from_atom_selection(const atom_selection_container_t &SelAtom, int imol,
                              const coot::Cartesian &front,
                              const coot::Cartesian &back,
                              short int pick_mode, bool verbose_mode) {

   float min_dist = 0.6; // should depend on zoom so that we can pick atom in faraway molecules
   int nearest_atom_index = -1;
   float dist = -999.9;
   pick_info p_i;

   if (false)
      std::cout << "pick_atom_from_atom_selection() imol " << imol
                << " n_selected_atoms " << SelAtom.n_selected_atoms << " "
                << front << " " << back << " pick_mode: " << pick_mode << std::endl;

   for (int i=0; i< SelAtom.n_selected_atoms; i++) {

      if (! SelAtom.atom_selection[i]->isTer()) {
         mmdb::Atom *at = SelAtom.atom_selection[i];
         coot::Cartesian atom_pos(at->x, at->y, at->z);
         if (atom_pos.within_box(front,back)) {
            dist = atom_pos.distance_to_line(front, back);
            // std::cout << "debug " << coot::atom_spec_t(at) << " "
            //           << atom_pos << "dist to line: " << dist << std::endl;

            if (dist < min_dist) {

               if ((pick_mode != PICK_ATOM_CA_ONLY) ||
                   (std::string(SelAtom.atom_selection[i]->name) == " CA ") ||
                   (std::string(SelAtom.atom_selection[i]->name) == " P  ")) {

                  std::string ele(SelAtom.atom_selection[i]->element);

                  if (((pick_mode == PICK_ATOM_NON_HYDROGEN) && (ele != " H")) ||
                      (pick_mode != PICK_ATOM_NON_HYDROGEN)) {

                     bool allow_pick = true;

                     // std::cout << "pick_mode: " << pick_mode << std::endl;

                     // 20101211 stop picking on regular residue atoms
                     // in CA+ligand mode
                     //
                     if (pick_mode == PICK_ATOM_CA_OR_LIGAND) {
                        std::string res_name = SelAtom.atom_selection[i]->GetResName();
                        std::string atom_name(SelAtom.atom_selection[i]->name);
                        // std::cout << "res_name: " << res_name << std::endl;
                        if (coot::util::is_standard_residue_name(res_name))
                        // no CAs in RNA/DNA and no Ps in protein.
                        if ((atom_name != " CA ") && (atom_name != " P  "))
                        allow_pick = 0;
                     }

                     if (pick_mode == PICK_ATOM_CA_OR_SIDECHAIN_OR_LIGAND) {
                        std::string res_name = SelAtom.atom_selection[i]->GetResName();
                        std::string atom_name(SelAtom.atom_selection[i]->name);
                        // std::cout << "res_name: " << res_name << std::endl;
                        if (coot::util::is_standard_residue_name(res_name))
                        // no CAs in RNA/DNA and no Ps in protein.
                        // Ignoring NA at the moment
                        if ((atom_name == " C  ") && (atom_name == " O  ") &&
                        (atom_name == " N  "))
                        allow_pick = 0;
                     }

                     if (allow_pick) {
                        min_dist = dist;
                        nearest_atom_index = i;
                        p_i.success = GL_TRUE;
                        p_i.atom_index = nearest_atom_index;
                        p_i.model_number = SelAtom.atom_selection[i]->GetModelNum();
                        p_i.imol = imol;
                        p_i.min_dist = dist;

                        if (verbose_mode) {
                           std::cout << "   DEBUG:: imol " << imol << " "
                                     << " atom index " << nearest_atom_index << std::endl;
                           std::cout << "   DEBUG:: imol " << imol << " "
                                     << SelAtom.atom_selection[i] << " " << min_dist
                                     << std::endl;
                        }
                     }
                  }
               } else {
                  if (verbose_mode) {
                     std::cout << "CA pick mode:" << std::endl;
                  }
               }
            }
         }
      }
   }
   return p_i;
}

// use the correct include file in the correct place
// glm::mat4 get_molecule_mvp();
#include <glm/gtx/string_cast.hpp>  // to_string()


std::pair<coot::Cartesian, coot::Cartesian>
graphics_info_t::get_front_and_back_for_pick() const {

   // modern version of getting front and back (the position in 3D space of the mouse on
   // the front clipping plane and the back clipping plane)
   GtkAllocation allocation = get_glarea_allocation();
   int w = allocation.width;
   int h = allocation.height;
   float mouseX = GetMouseBeginX() / (w * 0.5f) - 1.0f;
   float mouseY = GetMouseBeginY() / (h * 0.5f) - 1.0f;
   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 vp_inv = glm::inverse(mvp);
   // std::cout << "get_front_and_back_for_pick() mvp " << glm::to_string(mvp) << " back " << glm::to_string(vp_inv) << std::endl;
   float real_y = - mouseY; // in range -1 -> 1
   glm::vec4 screenPos_f = glm::vec4(mouseX, real_y, -1.0f, 1.0f);
   glm::vec4 screenPos_b = glm::vec4(mouseX, real_y,  1.0f, 1.0f); // or other way round?
   glm::vec4 worldPos_f = vp_inv * screenPos_f;
   glm::vec4 worldPos_b = vp_inv * screenPos_b;
   // std::cout << "get_front_and_back_for_pick() worldPos_f " << glm::to_string(worldPos_f) << " back " << glm::to_string(worldPos_b) << std::endl;
   float w_scale_f = 1.0/worldPos_f.w;
   float w_scale_b = 1.0/worldPos_b.w;
   coot::Cartesian front(worldPos_f.x * w_scale_f, worldPos_f.y * w_scale_f, worldPos_f.z * w_scale_f);
   coot::Cartesian  back(worldPos_b.x * w_scale_b, worldPos_b.y * w_scale_b, worldPos_b.z * w_scale_b);

   return std::pair<coot::Cartesian, coot::Cartesian>(front, back);
}

bool
graphics_info_t::tomo_pick(double x, double y, gint n_press, bool shift_is_pressed) {

   bool state = true; // unless we miss the box (currently not tested)

   if (shift_is_pressed) { // remove the previous pick

      std::string object_name =  "TomoPick " + std::to_string(tomo_view_info.section_index);
      int object_number = generic_object_index(object_name);
      from_generic_object_remove_last_item(object_number);

   } else {

      // normal path

      // let front F, back B:
      // delta = (B - F).uv()
      // equation of pick vector: F + t * delta because Ax + By + Cz = D and we know A and B are 0 and C is 1
      // For a pick on a z-section: P_z = F_z + t * delta_z
      // t = (P_z - F_z)/delta_z
      // x,y,z = F + t * delta

      std::pair<coot::Cartesian, coot::Cartesian> front_and_back = get_front_and_back_for_pick();
      float P_z = tomo_view_info.get_P_z();
      coot::Cartesian &front = front_and_back.first;
      coot::Cartesian &back  = front_and_back.second;
      coot::Cartesian delta = back - front;
      float t = (P_z - front.z())/delta.z();
      coot::Cartesian pick_point = front + delta * t;
      clipper::Coord_orth pt(pick_point.x(), pick_point.y(), pick_point.z());
      // std::cout << "pt: " << pt.format() << std::endl;
      coot::colour_holder ch(0.3, 0.3, 0.9);
      std::string object_name =  "TomoPick " + std::to_string(tomo_view_info.section_index);
      int object_number = generic_object_index(object_name);
      std::cout << "in tomo_pick A with object_number " << object_number << std::endl;
      if (object_number == -1)
         object_number = new_generic_object_number(object_name);
      std::cout << "in tomo_pick B with object_number " << object_number << std::endl;

      to_generic_object_add_point_internal(object_number, "dummy", ch, 3000.0, pt);
      set_display_generic_object(object_number, 1);
   }

   return state;
}

// Put these in graphics_info_t. Move this function into graphics-info-pick.cc
pick_info
graphics_info_t::atom_pick_gtk3(bool intermediate_atoms_only_flag) const {

   pick_info p_i;

   //GLenum err = glGetError(); std::cout << "atom_pick_gtk3() A err " << err << std::endl;
   // coot::Cartesian front = unproject(0.0);
   // coot::Cartesian back  = unproject(1.0);

   // modern version of getting front and back (the position in 3D space of the mouse on
   // the front clipping plane and the back clipping plane)
   GtkAllocation allocation = get_glarea_allocation();
   int w = allocation.width;
   int h = allocation.height;
   float screen_ratio = static_cast<float>(w)/static_cast<float>(h);
   float mouseX = GetMouseBeginX() / (w * 0.5f) - 1.0f;
   float mouseY = GetMouseBeginY() / (h * 0.5f) - 1.0f;
   glm::mat4 mvp = get_molecule_mvp();
   glm::mat4 vp_inv = glm::inverse(mvp);
   float real_y = - mouseY; // in range -1 -> 1
   glm::vec4 screenPos_f = glm::vec4(mouseX, real_y, -1.0f, 1.0f);
   glm::vec4 screenPos_b = glm::vec4(mouseX, real_y,  1.0f, 1.0f); // or other way round?
   glm::vec4 worldPos_f = vp_inv * screenPos_f;
   glm::vec4 worldPos_b = vp_inv * screenPos_b;
   float w_scale_f = 1.0/worldPos_f.w;
   float w_scale_b = 1.0/worldPos_b.w;
   coot::Cartesian front(worldPos_f.x * w_scale_f, worldPos_f.y * w_scale_f, worldPos_f.z * w_scale_f);
   coot::Cartesian  back(worldPos_b.x * w_scale_b, worldPos_b.y * w_scale_b, worldPos_b.z * w_scale_b);
   if (false)
      std::cout << mouseX << " " << real_y << " screen pos "
                << glm::to_string(screenPos_f) << " " << glm::to_string(screenPos_b)
                << " " << std::endl;

   // atom_pick() allows event to be null, in that case we don't check pick.
   // I don't follow what that is about at the moment.
   float dist_closest = 999999999999999999.9;
   int max_mol_no = n_molecules() - 1;

   auto l = [front, back](const molecule_class_info_t &m, int imol) {
               short int pick_mode = PICK_ATOM_ALL_ATOM;
               if (m.Bonds_box_type() == coot::CA_BONDS)	                     pick_mode = PICK_ATOM_CA_ONLY;
               if (m.Bonds_box_type() == coot::BONDS_NO_HYDROGENS)		     pick_mode = PICK_ATOM_NON_HYDROGEN;
               if (m.Bonds_box_type() == coot::CA_BONDS_PLUS_LIGANDS)                pick_mode = PICK_ATOM_CA_OR_LIGAND;
               if (m.Bonds_box_type() == coot::COLOUR_BY_RAINBOW_BONDS)		     pick_mode = PICK_ATOM_CA_OR_LIGAND; // yes, this mode shows ligands
               if (m.Bonds_box_type() == coot::CA_BONDS_PLUS_LIGANDS_AND_SIDECHAINS) pick_mode = PICK_ATOM_CA_OR_SIDECHAIN_OR_LIGAND;
               // this should be combined with the above, not override it.
               if (m.draw_hydrogens() == 0) pick_mode = PICK_ATOM_NON_HYDROGEN;
               bool verbose_mode = graphics_info_t::debug_atom_picking;
               pick_info mpi = pick_atom_from_atom_selection(m.atom_sel, imol, front, back, pick_mode, verbose_mode);
               return mpi;
            };

   if (intermediate_atoms_only_flag) {
      if (moving_atoms_asc) {
         if (moving_atoms_asc->n_selected_atoms > 0) {
            moving_atoms_molecule.atom_sel = *moving_atoms_asc; // why is this needed here? Actually, is it needed here?
                                                                // it shold be set when we start refinement.
            p_i = l(moving_atoms_molecule, -1);
         }
      }
   } else {
      // was a normal molecule
      for (int ii=max_mol_no; ii>=0; ii--) {
         if (molecules[ii].has_model()) {
            if (molecules[ii].atom_selection_is_pickable()) {
               const molecule_class_info_t &m = graphics_info_t::molecules[ii];
               pick_info mpi = l(m, ii);
               if (mpi.success) {
                  if (m.no_bonds_to_these_atom_indices.find(mpi.atom_index) ==  m.no_bonds_to_these_atom_indices.end()) {
                     if (mpi.min_dist < dist_closest) {
                        p_i = mpi;
                        dist_closest = mpi.min_dist;
                     }
                  }
               }
            }
         }
      }
      if (p_i.success) {
         // mmdb::Atom *at = molecules[p_i.imol].atom_sel.atom_selection[p_i.atom_index];
         // std::cout << "INFO:: picked atom: " << coot::atom_spec_t(at) << std::endl;
      }
   }
   return p_i;
}


// event can be null. if so Crtl key press check is not made.
pick_info
atom_pick() {

   graphics_info_t g;
   return g.atom_pick_gtk3(false);
}



pick_info
pick_intermediate_atom(const atom_selection_container_t &SelAtom) {

   graphics_info_t g;
   glm::vec4 f = g.unproject(0.0);
   glm::vec4 b = g.unproject(1.0);

   coot::Cartesian front(f.x, f.y, f.z);
   coot::Cartesian back(b.x, b.y, b.z);

   short int pick_mode = PICK_ATOM_ALL_ATOM;
   return pick_atom_from_atom_selection(SelAtom, -1, front, back, pick_mode, 0);
}



// Bang a point in the centre of the screen.  Bang another
// with x shifted by lets say 20 pixels.
//
// Find the sum of the differences of the fronts and the backs.
// That then is the vector that we move the rotation centre by.
// Plus some sort of scaling factor.
//
// We will use unproject, like in atom picking, but the code is
// different, because we are passing the position of the "mouse"
// x, y, not reading where the mouse was pressed.
//
coot::CartesianPair
screen_x_to_real_space_vector(GtkWidget *widget) {

   coot::ScreenVectors sv;
   return coot::CartesianPair(sv.screen_x, sv.screen_y);
}

coot::Cartesian
screen_z_to_real_space_vector(GtkWidget *widget) {

   coot::ScreenVectors sv;
   return sv.screen_z;

}
