/* src/pick.cc
 * 
 * Copyright 2002, 2003, 2004 by The University of York
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
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"

#include <string>
#include <vector> // for mmdb-crystal

#include <math.h>
#include <string.h>  // strncpy

#include <mmdb2/mmdb_manager.h>

#include "coords/cos-sin.h"
#include "coords/mmdb-extras.h"
#include "coords/mmdb.h"
#include "coords/mmdb-crystal.h" //need for Bond_lines now
#include "coords/Cartesian.h"
#include "coords/Bond_lines.h"

#include "graphics-info.h"

#if __APPLE__
#   include <OpenGL/glu.h>
#else
#   include <GL/glu.h>
#endif


#include "molecule-class-info.h"
#include "globjects.h"
#include "pick.h"
#include "cc-interface.hh" // for status bar text

pick_info
pick_atom(const atom_selection_container_t &SelAtom, int imol,
	  const coot::Cartesian &front, const coot::Cartesian &back, short int pick_mode,
	  bool verbose_mode) {

   float min_dist = 0.6;
   int nearest_atom_index = -1;
   float dist = -999.9;
   pick_info p_i;
   p_i.min_dist = 0; // keep compiler happy
   p_i.atom_index = -1; // ditto
   p_i.imol = -1; // ditto
   p_i.model_number = mmdb::MinInt4; // unset
   p_i.success = GL_FALSE; 
   for (int i=0; i< SelAtom.n_selected_atoms; i++) {

      if (! SelAtom.atom_selection[i]->isTer()) {

	 coot::Cartesian atom(SelAtom.atom_selection[i]->x,
			      SelAtom.atom_selection[i]->y,
			      SelAtom.atom_selection[i]->z);
	 //
	 if (atom.within_box(front,back)) { 

	    dist = atom.distance_to_line(front, back);

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
			      allow_pick = false;
		     }

                     if (pick_mode == PICK_ATOM_CA_OR_SIDECHAIN_OR_LIGAND) {
			std::string res_name = SelAtom.atom_selection[i]->GetResName();
			std::string atom_name(SelAtom.atom_selection[i]->name);
			// std::cout << "res_name: " << res_name << std::endl;
			if (coot::util::is_standard_residue_name(res_name))
			   // no CAs in RNA/DNA and no Ps in protein.
                           // Ignoring NA at the moment
			   if ((atom_name == " C  ") || (atom_name == " O  ") ||
                               (atom_name == " N  ") || (atom_name == " H  ") ||
                               (atom_name == " D  ") || (atom_name == " HA "))
			      allow_pick = false;
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

// event can be null. if so Crtl key press check is not made.
pick_info
atom_pick(GdkEventButton *event) { 

   coot::Cartesian front = unproject(0.0);

   coot::Cartesian back  = unproject(1.0);

   float dist_closest = 999999999999999999.9;
   int nearest_atom_index = -1;

   //cout << "front: " << front << endl;
   //cout << "back:  " << back  << endl;

   pick_info p_i;
   p_i.min_dist = 0; // keep compiler happy
   p_i.atom_index = -1; // ditto
   p_i.imol = -1; // ditto
   p_i.model_number = mmdb::MinInt4; // unset
   p_i.success = GL_FALSE;

   // 20210412 new style - let's check the intermediate atoms first

   graphics_info_t g;
   p_i = g.pick_moving_atoms(front, back); // will often/usually fail

   if (! p_i.success) {

      // Consider this senario: 2 molecules that have (at least in one
      // region) atom in the same place (a result of merging molecules,
      // for example):
      // 
      // Note that we want to select the atom that is "on top" (that is
      // the one the user is looking at) i.e. later on in the molecule
      // list - so we run through the list backwards.  This is
      // particularly important for deleting, because otherwise we often
      // don't see anything happen when something is deleted.
      //

      short int check_pick = 0;
      if (graphics_info_t::control_key_for_rotate_flag == 0) {

         if (event) {
            // i.e. control_key is for picking
            if (event->state & GDK_CONTROL_MASK) {
               check_pick = 1;
            }
         }
      } else {

         if (event) {
            // control_key is for rotation
            if (! (event->state & GDK_CONTROL_MASK)) {
               check_pick = 1;
            }
         }
      }

      if (check_pick) { 

         if (graphics_info_t::debug_atom_picking) {
            std::cout << "   == Level 2 atom picking diagnostic (send to Paul) ==\n";
         }
      
         int n_pickable = 0;
         int max_mol_no = graphics_info_t::n_molecules() - 1;
         for (int ii=max_mol_no; ii>=0; ii--) {

            if (graphics_info_t::molecules[ii].has_model()) { 
               if (graphics_info_t::molecules[ii].atom_selection_is_pickable()) {

                  n_pickable++;

                  atom_selection_container_t SelAtom = graphics_info_t::molecules[ii].atom_sel;
                  short int pick_mode = PICK_ATOM_ALL_ATOM;
                  if (graphics_info_t::molecules[ii].Bonds_box_type() == coot::CA_BONDS)
                     pick_mode = PICK_ATOM_CA_ONLY;
                  if (graphics_info_t::molecules[ii].Bonds_box_type() == coot::BONDS_NO_HYDROGENS)
                     pick_mode = PICK_ATOM_NON_HYDROGEN;
                  if (graphics_info_t::molecules[ii].Bonds_box_type() == coot::CA_BONDS_PLUS_LIGANDS)
                     pick_mode = PICK_ATOM_CA_OR_LIGAND;
                  if (graphics_info_t::molecules[ii].Bonds_box_type() == coot::CA_BONDS_PLUS_LIGANDS_AND_SIDECHAINS)
                     pick_mode = PICK_ATOM_CA_OR_SIDECHAIN_OR_LIGAND;
                  if (graphics_info_t::molecules[ii].Bonds_box_type() == coot::COLOUR_BY_RAINBOW_BONDS)
                     pick_mode = PICK_ATOM_CA_OR_LIGAND; // yes, this mode shows ligands

                  bool verbose_mode = graphics_info_t::debug_atom_picking;
                  pick_info mpi = pick_atom(SelAtom, ii, front, back, pick_mode, verbose_mode);
                  if (mpi.success) {
                     if (mpi.min_dist < dist_closest) {
                        p_i = mpi;
                        dist_closest = mpi.min_dist;
                     }
                  }
               }
            }
         }

         if (graphics_info_t::debug_atom_picking) {
            for (int ii=max_mol_no; ii>=0; ii--) {
               std::cout << "   MolNo " << ii << " of "
                         << graphics_info_t::n_molecules() << ":  " 
                         << graphics_info_t::molecules[ii].has_model() << " " 
                         << graphics_info_t::molecules[ii].is_displayed_p() << " " 
                         << graphics_info_t::molecules[ii].atom_selection_is_pickable() << "  "
                         << graphics_info_t::molecules[ii].atom_sel.n_selected_atoms << "  "
                         << graphics_info_t::molecules[ii].name_ << " "
                         << "\n";
            }
         }

         // we don't want to do this now that we have middle mouse panning
         //
         if (false) {
            if (n_pickable == 0) {
               std::string s = "There were no pickable (\"Active\") molecules!";
               GtkWidget *w = graphics_info_t::wrapped_nothing_bad_dialog(s);
               gtk_widget_show(w);
            }
         }

         //cout << "There were " << i_outside_count << " atoms outside "
         //	<< "the limits" << endl;

         // still something strange going on here, because p_i.success = -1073747613
         // on failure to find a (direct) hit (which is interpretted in glarea_button_press()
         // as failure fortuneately(?)).
         //
         //    cout << "p_i.imol    = " << p_i.imol << endl;
         //    cout << "p_i.success = " << p_i.success << endl;
         //    cout << "GL_FALSE    = " << GL_FALSE    << endl;

         if (p_i.success) {
            std::string ai;
            mmdb::Atom *at =
               graphics_info_t::molecules[p_i.imol].atom_sel.atom_selection[p_i.atom_index];
            std::string alt_conf_bit("");
            std::string segid = at->segID;
            if (strncmp(at->altLoc, "", 1))
               alt_conf_bit=std::string(",") + std::string(at->altLoc);
            atom_selection_container_t SelAtom = graphics_info_t::molecules[p_i.imol].atom_sel;
            nearest_atom_index = p_i.atom_index;

            std::cout << "(" << p_i.imol << ") \""
                      << (SelAtom.atom_selection)[nearest_atom_index]->name
                      << alt_conf_bit << "\"/"
                      << (SelAtom.atom_selection)[nearest_atom_index]->GetModelNum()
                      << "/chainid=\""
                      << (SelAtom.atom_selection)[nearest_atom_index]->GetChainID()
                      << "\"/"
                      << (SelAtom.atom_selection)[nearest_atom_index]->GetSeqNum()
                      << (SelAtom.atom_selection)[nearest_atom_index]->GetInsCode()
                      << "/"
                      << (SelAtom.atom_selection)[nearest_atom_index]->GetResName()
                      << ", "
                      << (SelAtom.atom_selection)[nearest_atom_index]->segID
                      << " occ: "
                      << (SelAtom.atom_selection)[nearest_atom_index]->occupancy
                      << " with B-factor: "
                      << (SelAtom.atom_selection)[nearest_atom_index]->tempFactor
                      << " element: \""
                      << (SelAtom.atom_selection)[nearest_atom_index]->element
                      << "\" at " << "("
                      << (SelAtom.atom_selection)[nearest_atom_index]->x << ","
                      << (SelAtom.atom_selection)[nearest_atom_index]->y << ","
                      << (SelAtom.atom_selection)[nearest_atom_index]->z << ")"
                      << " : " << dist_closest << std::endl;

            ai = atom_info_as_text_for_statusbar(nearest_atom_index, p_i.imol);

            gtk_statusbar_push(GTK_STATUSBAR(graphics_info_t::statusbar),
                               graphics_info_t::statusbar_context_id,
                               ai.c_str());
         }
      }
   }

   return p_i;
}



pick_info
pick_intermediate_atom(const atom_selection_container_t &SelAtom) {
   coot::Cartesian front = unproject(0.0);
   coot::Cartesian back  = unproject(1.0);

   std::cout << "---- here in pick_intermediate_atom() " << front << " " << back << std::endl;
   short int pick_mode = PICK_ATOM_ALL_ATOM;
   return pick_atom(SelAtom, -1, front, back, pick_mode, 0);
}


// Return the real world coordinates corresponding to the mouse
// postion, for various values of screen z (actually only 1.0 and 0.0);
// 
coot::Cartesian unproject(float screen_z) {

   graphics_info_t info;

   int x = int (rint(info.GetMouseBeginX()));
   int y = int (rint(info.GetMouseBeginY()));

   return unproject_xyz(x,y,screen_z);

}

coot::Cartesian unproject_xyz(int x, int y, float screen_z) {
   
   //cout << "using mouse coords: "
   //	<< x << " (" << info.GetMouseBeginX() << "), "
   //	<< y << " (" << info.GetMouseBeginY() << ") "
   //	<< endl;

   // This bit copied out the book
   // (but we change it so that it uses the z that we pass)
   //
   GLint viewport[4];
   GLdouble mvmatrix[16], projmatrix[16];
   GLint realy;  /*  OpenGL y coordinate position  */
   GLdouble wx, wy, wz;  /*  returned world x, y, z coords  */

      
   glGetIntegerv (GL_VIEWPORT, viewport);
   glGetDoublev (GL_MODELVIEW_MATRIX, mvmatrix);
   glGetDoublev (GL_PROJECTION_MATRIX, projmatrix);
   realy = viewport[3] - y - 1;
   double x_as_double = x;
   double realy_as_double = realy;
   /*  note viewport[3] is height of window in pixels  */
   //printf ("Coordinates at cursor are (%4d, %4d)\n", x, realy);
   gluUnProject (x_as_double, realy_as_double, screen_z, 
		 mvmatrix, projmatrix, viewport, &wx, &wy, &wz); 

   // you called glBegin(GL_LINES) before unprojecting, didn't you? :-)
   if (false) {
      std::cout << "unproject results in "
                << x_as_double << " " << realy_as_double << " " << screen_z << std::endl;
      std::cout << "unproject results out " << wx << " " << wy << " " << wz << std::endl;
   }

   return coot::Cartesian(wx, wy, wz);
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

   int x0 = widget->allocation.width/2;
   int y0 = widget->allocation.height/2;

   int x1 = x0 + 20;
   int y1 = y0 + 20;

   coot::Cartesian front_0 = unproject_xyz(x0, y0, 0.0);
   coot::Cartesian front_1 = unproject_xyz(x1, y0, 0.0);

   coot::Cartesian back_0 =  unproject_xyz(x0, y0, 1.0);
   coot::Cartesian back_1 =  unproject_xyz(x1, y0, 1.0);

   coot::Cartesian  back_diff_x =  back_1 -  back_0;
   coot::Cartesian front_diff_x = front_1 - front_0;

   coot::Cartesian sum_diff_x = back_diff_x + front_diff_x;

   // and same for y:
   //
   coot::Cartesian front_2 = unproject_xyz(x0, y1, 0.0);
   coot::Cartesian  back_2 = unproject_xyz(x0, y1, 1.0);

   coot::Cartesian  back_diff_y =  back_2 -  back_0;
   coot::Cartesian front_diff_y = front_2 - front_0;

   coot::Cartesian sum_diff_y = back_diff_y + front_diff_y;
   
   // This is a pair of coot::Cartesians.
   // 
   return coot::CartesianPair(sum_diff_x, sum_diff_y);

}

coot::Cartesian
screen_z_to_real_space_vector(GtkWidget *widget) { 

   int x0 = widget->allocation.width/2;
   int y0 = widget->allocation.height/2;

   coot::Cartesian front_0 = unproject_xyz(x0, y0, 0.0);
   coot::Cartesian back_0 =  unproject_xyz(x0, y0, 1.0);

   return (front_0 - back_0);
} 

