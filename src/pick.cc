 
#include <string>
#include <vector> // for mmdb-crystal

#include <math.h>
#include <gtkgl/gtkglarea.h>
#include <GL/gl.h>
#include <GL/glu.h>

#include "cos-sin.h"

#include "mmdb_manager.h"
#include "mmdb-extras.h"
#include "mmdb.h"
#include "mmdb-crystal.h" //need for Bond_lines now

#include "Cartesian.h"
#include "Bond_lines.h"

#include "graphics-info.h"

#include "molecule-class-info.h"

#include "cos-sin.h"

#include "globjects.h"

#include "pick.h"

//
pick_info
atom_pick(GdkEventButton *event) { 

   int i_outside_count = 0;
   
   
   coot::Cartesian front = unproject(0.0);

   coot::Cartesian back  = unproject(1.0);

   float dist, min_dist = 0.4;
   int nearest_atom_index = 0;

   //cout << "front: " << front << endl;
   //cout << "back:  " << back  << endl;

   pick_info p_i;
   p_i.success = GL_FALSE; 

//     std::cout << "There are " << graphics_info_t::n_molecules << " molecules"
// 	      << " to check in picking " << std::endl;


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
      // i.e. control_key is for picking
      if (event->state & GDK_CONTROL_MASK) {
	 check_pick = 1;
      }
   } else {
      // control_key is for rotation
      if (! (event->state & GDK_CONTROL_MASK)) {
	 check_pick = 1;
      }
   }

   if (check_pick) { 
      int n_pickable = 0;
      for (int ii=(graphics_info_t::n_molecules-1); ii>=0; ii--) {

	 if (graphics_info_t::molecules[ii].has_model()) { 
	    if (graphics_info_t::molecules[ii].atom_selection_is_pickable()) {

	       n_pickable++;

	       atom_selection_container_t SelAtom = graphics_info_t::molecules[ii].atom_sel; 
	 
	       for (int i=0; i< SelAtom.n_selected_atoms; i++) {
	 
		  coot::Cartesian atom( (SelAtom.atom_selection)[i]->x,
					(SelAtom.atom_selection)[i]->y,
					(SelAtom.atom_selection)[i]->z);
		  //
		  if (atom.within_box(front,back)) { 
		     dist = atom.distance_to_line(front, back);

		     if (dist < min_dist) {

			min_dist = dist;
			nearest_atom_index = i;
			p_i.success = GL_TRUE;
			p_i.atom_index = nearest_atom_index;
			p_i.imol = ii;

			// make sure that we are looking at the molecule that had
			// the nearest contact:
			// 
			SelAtom = graphics_info_t::molecules[p_i.imol].atom_sel; 

			std::string alt_conf_bit("");
			CAtom *at = SelAtom.atom_selection[nearest_atom_index];
			if (strncmp(at->altLoc, "", 1))
			   alt_conf_bit=std::string(",") + std::string(at->altLoc);
		  
			cout << "(" << p_i.imol << ") " 
			     << (SelAtom.atom_selection)[nearest_atom_index]->name 
			     << alt_conf_bit << "/"
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
			     << " element: "
			     << (SelAtom.atom_selection)[nearest_atom_index]->element
			     << " at " << "("
			     << (SelAtom.atom_selection)[nearest_atom_index]->x << ","
			     << (SelAtom.atom_selection)[nearest_atom_index]->y << ","
			     << (SelAtom.atom_selection)[nearest_atom_index]->z << ")"
			     << " : " << min_dist << endl;
			// std::cout << "DEBUG:: atom index " << nearest_atom_index << std::endl;

		     }
	    
		  } else {
		     i_outside_count ++;
		     // 		  std::cout << "atom " << i << " outside box " 
		     // 			    << front << " " << back << "\n";

		  }
	       }
	    }
	 }
      }

      if (n_pickable == 0) {
	 std::string s = "There were no pickable (\"Active\") molecules!";
	 GtkWidget *w = graphics_info_t::wrapped_nothing_bad_dialog(s);
	 gtk_widget_show(w);
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
	 CAtom *at =
	    graphics_info_t::molecules[p_i.imol].atom_sel.atom_selection[p_i.atom_index];
	 std::string alt_conf_bit("");
	 std::string segid = at->segID;
	 if (strncmp(at->altLoc, "", 1))
	    alt_conf_bit=std::string(",") + std::string(at->altLoc);
	 ai += "(mol. no: ";
	 ai += graphics_info_t::int_to_string(p_i.imol);
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
	 if (segid != "")
	    ai += segid;
	 ai += " occ: ";
	 ai += graphics_info_t::float_to_string(at->occupancy);
	 ai += " bf: ";
	 ai += graphics_info_t::float_to_string(at->tempFactor);
	 ai += " ele: ";
	 ai += at->element;
	 ai += " pos: (";
	 ai += graphics_info_t::float_to_string(at->x);
	 ai += ",";
	 ai += graphics_info_t::float_to_string(at->y);
	 ai += ",";
	 ai += graphics_info_t::float_to_string(at->z);
	 ai += ")";
	 gtk_statusbar_push(GTK_STATUSBAR(graphics_info_t::statusbar),
			    graphics_info_t::statusbar_context_id,
			    ai.c_str());
      }
   }
   
   return p_i;
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
   // printf ("World coords at z=%f are (%f, %f, %f)\n", 
   // screen_z, wx, wy, wz);

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

