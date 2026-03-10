/* src/graphics-info.cc
 *
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by the University of York
 * Copyright 2007, 2008, 2009 by the University of Oxford
 * Copyright 2014, 2015, 2016 by Medical Research Council
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



#include "coords/Cartesian.hh"
#ifdef USE_PYTHON
#include "python-3-interface.hh"
#endif

#include "compat/coot-sysdep.h"

#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif

#include <gtk/gtk.h>  // must come after mmdb_manager on MacOS X Darwin

#include <iostream>
#ifdef _MSC_VER
#include "compat/dirent.h"
#else
#include <dirent.h>   // for refmac dictionary files
#endif

#include <sys/types.h> // for stating
#include <sys/stat.h>

#if !defined _MSC_VER && !defined WINDOWS_MINGW
#include <unistd.h>
#else
//#include "coot-sysdep.h"
#endif

#include <mmdb2/mmdb_manager.h>
#include "coords/mmdb.hh"
#include "coords/Bond_lines.hh"

#include "clipper/core/map_utils.h" // Map_stats
#include "skeleton/graphical_skel.h"

#include "graphics-info.h"

#include "molecule-class-info.h"
#include "skeleton/BuildCas.h"
#include "utils/coot-utils.hh"
#include "manipulation-modes.hh"

#ifdef USE_GUILE

#include <cstdio> /* for std::FILE in gmp.h for libguile.h */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wvolatile"
#include <libguile.h>
#pragma GCC diagnostic pop
#endif

#include "cc-interface.hh" // needed for display_all_model_molecules()

#ifdef USE_PYTHON
#include "cc-interface-scripting.hh"
#endif

#include "geometry/dict-utils.hh"

#include "interface.h"
#include "widget-from-builder.hh"
#include "draw-2.hh"
#include "pick.hh"
#include "utils/logging.hh"
extern logging logger;

// static
GtkWidget *
graphics_info_t::get_widget_from_builder(const std::string &w_name) { // use gtkbuilder to do new-style lookup_widget();

   GtkWidget *w = GTK_WIDGET(gtk_builder_get_object(gtkbuilder, w_name.c_str()));
   return w;
}

// static
GObject *
graphics_info_t::get_gobject_from_builder(const std::string &w_name) { // use gtkbuilder but return a gobject (for menus)

   GObject *o = G_OBJECT(gtk_builder_get_object(gtkbuilder, w_name.c_str()));
   return o;
}


// static
GtkWidget *
graphics_info_t::get_widget_from_preferences_builder(const std::string &w_name) { // use gtkbuilder to do new-style lookup_widget();

   if (false)
      std::cout << "debug:: in get_widget_from_preferences_builder() using builder " << preferences_gtkbuilder
                << " to lookup " << w_name << std::endl;
   GtkWidget *w = GTK_WIDGET(gtk_builder_get_object(preferences_gtkbuilder, w_name.c_str()));
   return w;
}

// return a vector of the current valid map molecules
std::vector<int>
graphics_info_t::valid_map_molecules() const {

   std::vector<int> v;
   for (unsigned int i=0; i<molecules.size(); i++)
      if (is_valid_map_molecule(i))
    v.push_back(i);
   return v;
}

// static
GtkAllocation graphics_info_t::get_glarea_allocation() {
   GtkAllocation allocation;
   if (!glareas.empty())
      gtk_widget_get_allocation(glareas[0], &allocation);
   return allocation;
}


// return the new molecule number
// static
int graphics_info_t::create_molecule() {
   int imol = molecules.size();
   try {
      molecules.push_back(molecule_class_info_t(imol));
   }
   catch (const std::bad_alloc &ba) {
      std::cout << "ERROR:: bad_alloc: " << ba.what() << std::endl;
      imol = -1;
   }
   return imol;
}

void
graphics_info_t::post_recentre_update_and_redraw() {

   //
   int t0 = 0; // glutGet(GLUT_ELAPSED_TIME);
   std::cout << "Fix timer in post_recentre_update_and_redraw()\n";
   for (int ii=0; ii<n_molecules(); ii++) {
      molecules[ii].update_clipper_skeleton();
      molecules[ii].update_map(auto_recontour_map_flag);  // uses statics in graphics_info_t
                                                          // and redraw the screen using the new map
   }

   // int t1 = 0; // glutGet(GLUT_ELAPSED_TIME);
   // std::cout << "Elapsed time for map contouring: " << t1-t0 << "ms" << std::endl;

   for (int ii=0; ii<n_molecules(); ii++) {
      // std::cout << "update symmetry  for ii " << ii << std::endl;
      molecules[ii].update_symmetry();
   }
   make_pointer_distance_objects();
   graphics_draw();
}


#ifndef EMSCRIPTEN
GdkRGBA colour_by_distortion(float dist) {

   GdkRGBA col;

   col.alpha = 1;
   col.blue  = 0;

   if (dist < 0.0) {
      // black for negative numbers
      col.red   = 0;
      col.green = 0;
   } else {
      if (dist < 1.4 /* was 2.0 before Tickle-fix */) {
         col.red   = 0;
         col.green = 55535;
      } else {
         if (dist < 2.2 /* was 5.0 */ ) {
            col.red   = 55000;
            col.green = 55000;
            // col.blue  = 22000;
         } else {
            if (dist < 3.0 /* was 8.0 */ ) {
               col.red   = 64000;
               col.green = 32000;
            } else {
               col.red   = 65535;
               col.green = 0;
            }
         }
      }
   }
   return col;
}
#endif

GdkRGBA colour_by_rama_plot_distortion(float plot_value, int rama_type) {

   if (true)
      std::cout << "in colour_by_rama_plot_distortion plot_value "
                << plot_value << " rama_type " << rama_type
                << " c.f. coot::RAMA_TYPE_LOGRAMA " << coot::RAMA_TYPE_LOGRAMA
                << " coot::RAMA_TYPE_ZO " << coot::RAMA_TYPE_ZO
                << std::endl;

   // ZO type data need to scaled to match
   // 20*zo_type_data-80 = log_rama_type_data
   //
   // if (rama_type == coot::RAMA_TYPE_ZO)
   //   plot_value = 20*plot_value -80;

   GdkRGBA col;

   col.alpha = 1;
   col.blue  = 0;

   auto rotation_size_raw_to_gdkcol = [] (float rotation_size_raw) {
                                         float rotation_size = -0.33f * rotation_size_raw; // cooked
                                         // std::vector<float> orig_colours = { 0.0f,  0.8f, 0.0f };
                                         // std::vector<float> rgb_new = rotate_rgb(orig_colours, rotation_size);
                                         coot::colour_holder ch(0.0f, 0.8f, 0.0f);
                                         ch.rotate_by(rotation_size);
                                         GdkRGBA col;
                                         col.alpha = 1;
                                         col.red   = ch.red * 255.0 * 255.0;
                                         col.green = ch.green * 255.0 * 255.0;
                                         col.blue  = ch.blue * 255.0 * 255.0;
                                         return col;
                                      };

   if (rama_type == coot::RAMA_TYPE_LOGRAMA) {
      // This used to be true:
      // now that RAMA_TYPE_ZO are on the same scale this colour
      // scheme will do for both
      // But then I changed the weight on ZO rama
      // So colours need to be different

      // the range of good to bad rama plot score is -18 to -8. That should be mapped to
      // rotation_size_raw of 0.0 to 1.0.
      float rotation_size_raw = 0.0;
      if (plot_value > -18.0) {
         rotation_size_raw = (plot_value + 18.0f) / (-8.0f - -18.0f);
         if (rotation_size_raw > 1.0f)
            rotation_size_raw = 1.0f;
      }
      col = rotation_size_raw_to_gdkcol(rotation_size_raw);

   } else {
      // RAMA_TYPE_ZO.  -2.5 is bad. -5 is good
      //
      // That should be mapped to rotation_size_raw of 0.0 to 1.0.
      float rotation_size_raw = 0.0;
      if (plot_value > -5.0f) {
         rotation_size_raw = (plot_value + 5.0f) / (-2.5f - -5.0f);
         if (rotation_size_raw > 1.0f)
            rotation_size_raw = 1.0f;
      }
      col = rotation_size_raw_to_gdkcol(rotation_size_raw);

   }
   return col;
}


// static
int
graphics_info_t::get_latest_model_molecule() {

   int imol = -1;
   int n = n_molecules();
   for(int ii=0; ii<n; ii++) {
      if (is_valid_model_molecule(ii)) {
         if (ii > imol) {
            imol = ii;
         }
      }
   }
   return imol;
}

//static
int
graphics_info_t::get_biggest_model_molecule() {

   int imol = -1;
   int n_atoms_max = -1;
   int n = n_molecules();
   for(int ii=0; ii<n; ii++) {
      if (is_valid_model_molecule(ii)) {
         int n_atoms_mol = molecules[ii].atom_sel.n_selected_atoms;
         if (n_atoms_mol > n_atoms_max) {
            imol = ii;
            n_atoms_max = n_atoms_mol;
         }
      }
   }
   return imol;
}




double graphics_info_t::GetMouseBeginX() const { return mouse_begin.first; };

double graphics_info_t::GetMouseBeginY() const { return mouse_begin.second; };

void graphics_info_t::SetMouseBegin(double x, double y) {
   mouse_begin.first  = x;
   mouse_begin.second = y;
}

void graphics_info_t::SetMouseClicked(double x, double y) {
   mouse_clicked_begin.first  = x;
   mouse_clicked_begin.second = y;
}


// static
GtkWidget *graphics_info_t::wrapped_nothing_bad_dialog(const std::string &label, bool use_markup) {

   auto add_image_widgets_if_needed = [] (GtkWidget *box) {
      GtkWidget *ch = gtk_widget_get_first_child(box);
      if (ch) {
         // already added - so do nothing
      } else {
         std::string prefix = coot::prefix_dir();
         std::string dir = coot::util::append_dir_dir(prefix, "share/icons/hicolor/scalable/actions");
         std::string fn_1 = coot::util::append_dir_file(dir, "Stock-dialog-information.svg");
         std::string fn_2 = coot::util::append_dir_file(dir, "Stock-dialog-warning.svg");
         GtkWidget *image_1 = gtk_image_new_from_file(fn_1.c_str());
         GtkWidget *image_2 = gtk_image_new_from_file(fn_2.c_str());
         gtk_box_append(GTK_BOX(box), image_1);
         gtk_box_append(GTK_BOX(box), image_2);
         g_object_set_data(G_OBJECT(box), "information", image_1);
         g_object_set_data(G_OBJECT(box), "warning",     image_2);
         gtk_widget_set_size_request(image_1, 80, 80);
         gtk_widget_set_size_request(image_2, 80, 80);
      }
   };

   GtkWidget *dialog = NULL;
   if (use_graphics_interface_flag) {

      dialog = widget_from_builder("nothing_bad_dialog");
      GtkWidget *box = widget_from_builder("nothing_bad_image_box");
      add_image_widgets_if_needed(box);

      GtkWidget *label_widget = widget_from_builder("nothing_bad_label");

      gtk_widget_set_visible(label_widget, TRUE);
      gtk_label_set_text(GTK_LABEL(label_widget), label.c_str());

      // are these correct?
      gtk_label_set_xalign(GTK_LABEL(label_widget), 0.0);
      gtk_label_set_use_markup(GTK_LABEL(label_widget), TRUE);

      if (use_markup) {
	 gtk_label_set_justify(GTK_LABEL(label_widget), GTK_JUSTIFY_LEFT);
	 gtk_label_set_markup(GTK_LABEL(label_widget), label.c_str());
      }

      GtkWidget *main_window = graphics_info_t::get_main_window();
      gtk_window_set_transient_for(GTK_WINDOW(dialog), GTK_WINDOW(main_window));
      gtk_widget_set_visible(dialog, TRUE);


      // Handle the info and warning icon
      //
      bool warning = false;
      if (label.find(std::string("WARNING")) != std::string::npos) warning = true;
      if (label.find(std::string("warning")) != std::string::npos) warning = true;
      if (label.find(std::string("Warning")) != std::string::npos) warning = true;
      if (label.find(std::string("Oops!"))   != std::string::npos) warning = true;
      GtkWidget *info_image = GTK_WIDGET(g_object_get_data(G_OBJECT(box), "information"));
      GtkWidget *warn_image = GTK_WIDGET(g_object_get_data(G_OBJECT(box), "warning"));
      if (warning) {
         gtk_widget_set_visible(GTK_WIDGET(info_image), FALSE);
         gtk_widget_set_visible(GTK_WIDGET(warn_image), TRUE);
      } else {
         gtk_widget_set_visible(GTK_WIDGET(info_image), TRUE);
         gtk_widget_set_visible(GTK_WIDGET(warn_image), FALSE);
      }
   }
   return dialog;
}

void
graphics_info_t::set_do_anti_aliasing(int state) {

  short int old_flag = graphics_info_t::do_anti_aliasing_flag;
  graphics_info_t::do_anti_aliasing_flag = state;
  if (do_anti_aliasing_flag != old_flag) {
    draw_anti_aliasing();
  }
}

// static
bool
graphics_info_t::background_is_black_p() {

   bool v = false;
   if (background_colour[0] < 0.3)
      if (background_colour[1] < 0.3)
    if (background_colour[2] < 0.3)
       v = true;

   return v;
}

void
graphics_info_t::draw_anti_aliasing() {

}

// This addresses the "everything is an INH" problem.
//
// imol_enc can be a specific model molecule number or
// IMOL_ENC_AUTO, IMOL_ENC_ANY are the interesting values
// otherwise mol number.
//
// if imol_enc_in is IMOL_ENC_AUTO, then try to find to which
// molecule this dictionary refers.
// If the residue type is on the non-auto load list, simply go through
// the molecule list backwards, starting from the hightest molecule number looking for
// a molecule that is a valid model molecule - that's the one.
// If the residue type is not in the non-auto list, then
// it is a dictionary for all molecules, i.e. IMOL_ENC_ANY.
//
// return the index of the monomer in the geometry store. Return -1 on failure
//
coot::read_refmac_mon_lib_info_t
graphics_info_t::add_cif_dictionary(std::string cif_dictionary_filename,
                                    int imol_enc_in,
                                    short int show_no_bonds_dialog_maybe_flag) {

   if (false)
      std::cout << "::: add_cif_dictionary() called with "
                << cif_dictionary_filename << " " << imol_enc_in << " "
                << show_no_bonds_dialog_maybe_flag << std::endl;

   int imol_enc = imol_enc_in;

   if (imol_enc_in == coot::protein_geometry::IMOL_ENC_AUTO) {
      std::vector<std::string> comp_ids = coot::comp_ids_in_dictionary_cif(cif_dictionary_filename);
      bool is_non_auto_load_comp_id = false;  // because it is ATP, not LIG
      for (unsigned int i=0; i<comp_ids.size(); i++) {
         if (geom_p->is_non_auto_load_ligand(comp_ids[i])) {
            // imol_enc is the latest model added that contains this comp_id
            //
            is_non_auto_load_comp_id = true;

            for (int ii=(n_molecules()-1); ii>=0; ii--){
               if (is_valid_model_molecule(ii)) {
                  imol_enc = ii;
                  break;
               }
            }
            break;
         }
      }
      if (! is_non_auto_load_comp_id)
         imol_enc = coot::protein_geometry::IMOL_ENC_ANY;
   }

   coot::read_refmac_mon_lib_info_t rmit =
      geom_p->init_refmac_mon_lib(cif_dictionary_filename,
                                  cif_dictionary_read_number,
                                  imol_enc);

   cif_dictionary_read_number++;
   if (rmit.success > 0) {
      cif_dictionary_filename_vec->push_back(cif_dictionary_filename);
      if (show_no_bonds_dialog_maybe_flag) {
         display_density_level_this_image = 1;
         std::string s;
         s = "Read ";
         s += int_to_string(rmit.n_atoms + rmit.n_links);
         s += " atoms/links in restraints from ";
         s += cif_dictionary_filename;
         display_density_level_screen_string = s;
         add_status_bar_text(s);
         graphics_draw();
      }
      std::cout << display_density_level_screen_string << std::endl;
   } else {
      std::cout << "init_refmac_mon_lib "  << cif_dictionary_filename
                << " had no bond restraints\n";
      if (use_graphics_interface_flag) {
         if (show_no_bonds_dialog_maybe_flag) {
            // GtkWidget *widget = create_no_cif_dictionary_bonds_dialog();
            GtkWidget *widget = widget_from_builder("no_cif_dictionary_bonds_dialog");
            gtk_widget_set_visible(widget, TRUE);
         }
      }

      std::string s;
      for (unsigned int i=0; i<rmit.error_messages.size(); i++) {
         s += rmit.error_messages[i];
         s += "\n";
      }
      info_dialog(s);
   }

   // Redraw all molecules! Yikes!
   for (unsigned int i=0; i<molecules.size(); i++) {
      if (is_valid_model_molecule(i)) {
         molecules[i].make_bonds_type_checked(__FUNCTION__);
      }
   }
   // return rmit.n_atoms;
   return rmit;
}


// Ideally don't redraw everything, just those that have a residue with name res_name
//
void
graphics_info_t::redraw_molecules_with_residue(const std::string &res_name) {

   for (unsigned int i=0; i<molecules.size(); i++) {
      if (is_valid_model_molecule(i)) {
    if (molecules[i].has_residue_with_name(res_name)) {
       molecules[i].make_bonds_type_checked(__FUNCTION__);
    }
      }
   }
   graphics_draw();
}


void
graphics_info_t::import_all_refmac_cifs() {

   char *env = getenv("COOT_REFMAC_LIB_DIR");
   if (! env) {
      std::cout << "Can't import dictionary because COOT_REFMAC_LIB_DIR is not defined\n";
   } else {

      std::string coot_refmac_lib_dir(env);

      // is coot_refmac_lib_dir a directory?

      struct stat buf;
      int status = stat(coot_refmac_lib_dir.c_str(), &buf);
      if (status != 0) {
         std::cout << "Error finding directory " << coot_refmac_lib_dir << std::endl;
      } else {
         if (S_ISDIR(buf.st_mode)) {
            std::cout << coot_refmac_lib_dir << " is a directory (good). " << std::endl;

            std::string data_dir = add_dir_file(coot_refmac_lib_dir, "data");
            std::string monomer_dir = add_dir_file(data_dir, "monomers");

            // good

            DIR *lib_dir = opendir(monomer_dir.c_str());
            if (lib_dir == NULL) {
               std::cout << "An ERROR occured on opening the directory "
                         << monomer_dir << std::endl;
            } else {

               struct dirent *dir_ent;

               // loop until the end of the filelist (readdir returns NULL)
               //
               while (1) {
                  dir_ent = readdir(lib_dir);
                  if (dir_ent == NULL) {
                     break;
                  } else {

                     std::string sub_dir_part(std::string(dir_ent->d_name));

                     if ( ! (sub_dir_part == ".") ) {
                        std::string subdirname = add_dir_file(monomer_dir, sub_dir_part);

                        // we need to test that sub_dir_part is a directory:
                        // (if not, silently skip over it)
                        //
                        status = stat(subdirname.c_str() , &buf);
                        if (S_ISDIR(buf.st_mode)) {

                           DIR *sub_dir = opendir(subdirname.c_str());

                           if (sub_dir == NULL) {
                              std::cout << "An ERROR occured on opening the subdirectory "
                                        << subdirname << std::endl;
                           } else {

                              struct dirent *sub_dir_ent;

                              while (1) {
                                 sub_dir_ent = readdir(sub_dir);
                                 if (sub_dir_ent == NULL) {
                                    break;
                                 } else {
                                    std::string cif_filename =
                                       add_dir_file(subdirname, std::string(sub_dir_ent->d_name));
                                    status = stat(cif_filename.c_str(), &buf);
                                    if (status == 0) {
                                       if (S_ISREG(buf.st_mode)) {
                                          add_cif_dictionary(cif_filename,
                                                             coot::protein_geometry::IMOL_ENC_ANY, 0);
                                       }
                                    }
                                 }
                              }
                           }
                           closedir(sub_dir);
                        }
                     } // not "."
                  }
               }
               closedir(lib_dir);
            }
         } else {
            std::cout << "Failure to import - " << coot_refmac_lib_dir
                      << " is not a directory\n";
         }
      }
   }
}

// a static
std::string
graphics_info_t::add_dir_file(const std::string &dirname, const std::string &filename) {

   std::string r = dirname;
   r += "/";
   r += filename;
   return r;
}

// if dir is true, we want to go forward
void
graphics_info_t::reorienting_next_residue(bool dir) {

   std::pair<int, mmdb::Atom *> atom_pair = get_active_atom();

   bool done_it = false;

   if (atom_pair.second) {

      int imol = atom_pair.first;
      mmdb::Residue *residue_current = atom_pair.second->residue;

      // check direction dir here
      mmdb::Residue *residue_next = 0;
      if (dir)
         residue_next = coot::util::next_residue(residue_current);
      else
         residue_next = coot::util::previous_residue(residue_current);

      if (residue_next) {
         std::pair<bool, clipper::RTop_orth> ro =
            coot::util::get_reorientation_matrix(residue_current, residue_next);
         std::pair<bool, clipper::Coord_orth> residue_centre =
            molecules[imol].residue_centre(residue_next);

         if (ro.first) {

            // Happy case

            // ----------- target rotation -------------

            clipper::Mat33<double> cr = ro.second.rot();
            glm::mat3 m(cr(0,0), cr(0,1), cr(0,2),
                        cr(1,0), cr(1,1), cr(1,2),
                        cr(2,0), cr(2,1), cr(2,2));
            glm::quat qq(m);
            glm::quat target_quat = view_quaternion * qq;

            //  --------- target position - try the residue and change it if we started on a CA

            const clipper::Coord_orth &rc = residue_centre.second;
            coot::Cartesian res_centre(rc.x(), rc.y(), rc.z());
            coot::Cartesian rot_centre = RotationCentre();
            coot::Cartesian target_pos = res_centre;

            // however, if the next residue has a CA atom
            // then we should centre on the CA of the next
            // residue (rather than the mean position)

            mmdb::Atom *at = atom_pair.second;
            std::string atom_name(at->GetAtomName());
            if (atom_name == " CA ") { // PDBv3 FIXME
               coot::Cartesian ca_pos(at->x, at->y, at->z);
               std::pair<bool, clipper::Coord_orth> ca_next_pos =
                  coot::util::get_CA_position_in_residue(residue_next);
               if (ca_next_pos.first) {
                  target_pos = coot::Cartesian(ca_next_pos.second.x(),
                                               ca_next_pos.second.y(),
                                               ca_next_pos.second.z());
               }
            }

            smooth_scroll_delta = target_pos - rot_centre;
            set_old_rotation_centre(RotationCentre());

            if (smooth_scroll == 1) {

               coot::view_info_t view1(view_quaternion,    rot_centre, zoom, "current");
               coot::view_info_t view2(target_quat, target_pos, zoom, "next");
               int nsteps = smooth_scroll_n_steps * 2;

               coot::view_info_t::interpolate(view1, view2, nsteps); // sets up a gtk_widget_add_tick_callback with
                                                                     // a (no-capture) lambda function, so we don't
                                                                     // want to change thew view here.

            } else {

               // "snap" the view
               rotation_centre_x = target_pos.x();
               rotation_centre_y = target_pos.y();
               rotation_centre_z = target_pos.z();

               // Replace this with updating the glm_quat, then redraw

               // bleugh again
               // for(int i=0; i<4; i++) quat[i] = vqf[i];

               try_centre_from_new_go_to_atom();
               update_things_on_move_and_redraw(); // (symmetry, environment, rama, map) and draw it

            }

            done_it = true;
            go_to_atom_chain_       = residue_next->GetChainID();
            go_to_atom_residue_     = residue_next->GetSeqNum();
            go_to_atom_inscode_     = residue_next->GetInsCode();


            // if (go_to_atom_window) {
               // what is next_atom here? Hmm
               // update_widget_go_to_atom_values(go_to_atom_window, next_atom);
            // }
         }
      }
   }
   if (! done_it) {
      // Oops! Next residue was not found, back to normal/standard/old mode
      if (dir)
         intelligent_next_atom_centring();
      else
         intelligent_previous_atom_centring();
   }
   // graphics_draw();
}

// static
void
graphics_info_t::set_rotation_centre(const clipper::Coord_orth &pt) {

   graphics_info_t g;
   coot::Cartesian centre(pt.x(), pt.y(), pt.z());
   bool done_centre_jump = g.setRotationCentre(centre);
   if (done_centre_jump)
      g.update_things_on_move_and_redraw();
}


void
graphics_info_t::setRotationCentre(int index, int imol) {

   mmdb::Atom *atom = molecules[imol].atom_sel.atom_selection[index];
   float x = atom->x;
   float y = atom->y;
   float z = atom->z;
   clipper::Coord_orth pt(x,y,z);
   set_rotation_centre(pt);

   if (environment_distance_label_atom) {
      molecules[imol].unlabel_last_atom();
      molecules[imol].add_to_labelled_atom_list(index);
   }

}

// update the green square, where we are.
#if DO_RAMA_PLOT
void
graphics_info_t::update_ramachandran_plot_point_maybe(int imol, mmdb::Atom *atom) {

   coot::residue_spec_t r(atom->residue);
   update_ramachandran_plot_point_maybe(imol, r);
}
#endif

#if DO_RAMA_PLOT
void
graphics_info_t::update_ramachandran_plot_point_maybe(int imol, const coot::residue_spec_t &res_spec) {

}
#endif // DO_RAMA_PLOT


#ifdef DO_RAMA_PLOT
// called from accept_moving_atoms()
void
graphics_info_t::update_ramachandran_plot_point_maybe(int imol, atom_selection_container_t moving_atoms) {

   // get the centre residue from moving atoms when 1 or 3 residue are
   // being refined and call update_ramachandran_plot_point_maybe()
   // with the residue spec.
   std::pair<bool, std::pair<int, coot::atom_spec_t> > aa_spec_pair = active_atom_spec();
   if (aa_spec_pair.first) {
      coot::residue_spec_t r(aa_spec_pair.second.second);
      update_ramachandran_plot_point_maybe(imol, r);
   }
}
#endif



// We need to know the atom index and imol of the last centred atom
// before this gets activated. The atom index and imol must be right
// because they are used without protection i) to get to
// the right molecule_class_info_t and then the right atom.
//
void
graphics_info_t::update_environment_distances_maybe(int index, int imol) {

   if (environment_show_distances) {
      if (go_to_atom_molecule() < n_molecules()) {
 	 if (is_valid_model_molecule(imol)) {
 	    update_environment_graphics_object(index, imol);
 	    if (show_symmetry)
  	       update_symmetry_environment_graphics_object(index, imol);
 	 }
      }
   }
}

// Return imol = -1 if no (close) atoms found.
//
// Ignore molecules that are not displayed.
//
// index, imol
std::pair<int, int>
graphics_info_t::get_closest_atom() const {
   // int index, int imol

   std::pair <float, int> dist_info;
   float dist_min = 999999999.0;
   coot::Cartesian rc = RotationCentre();
   int imol_close = -1;
   int index_close = -1;

   for (int imol=0; imol<n_molecules(); imol++) {

      if (molecules[imol].has_model()) {
    if (molecules[imol].is_displayed_p()) {
       dist_info = molecules[imol].nearest_atom(rc);
       if (dist_info.first < dist_min) {
          imol_close = imol;
          index_close = dist_info.second;
          dist_min = dist_info.first;
       }
    }
      }
   }
   return std::pair<int, int>(index_close, imol_close);
}

void
graphics_info_t::setRotationCentre(const symm_atom_info_t &symm_atom_info) {

   std::cout << "setRotationCentre by symmetry atom" << std::endl;

   // Invalid read according to valgrind
   mmdb::PAtom atom = symm_atom_info.trans_sel[symm_atom_info.atom_index];

   if (atom) {
      float x = atom->x; // invalid read according to valgrind.
      float y = atom->y; // ditto
      float z = atom->z; // ditto

      rotation_centre_x = x;
      rotation_centre_y = y;
      rotation_centre_z = z;
   } else {
      std::cout << "ERROR:: NULL atom in setRotationCentre(symm_atom_info_t)\n";
   }
}

void
graphics_info_t::setRotationCentre(const coot::clip_hybrid_atom &hybrid_atom) {

   if (false)
      // std::cout << "INFO:: setRotationCentre by symmetry hybrid atom "
      //           << hybrid_atom.atom << " at "
      //           << hybrid_atom.pos << std::endl;
      logger.log(log_t::INFO, "setRotationCentre by symmetry hybrid atom at " +
                 std::to_string(hybrid_atom.pos.x()) + " " +
                 std::to_string(hybrid_atom.pos.y()) + " " +
                 std::to_string(hybrid_atom.pos.z()));

   rotation_centre_x = hybrid_atom.pos.x();
   rotation_centre_y = hybrid_atom.pos.y();
   rotation_centre_z = hybrid_atom.pos.z();

}


bool
graphics_info_t::smooth_scroll_maybe(float x, float y, float z,
                                     bool do_zoom_and_move_flag,
                                     float target_zoom) {

   bool done = false;
   if ( (x - rotation_centre_x) != 0.0 ||
        (y - rotation_centre_y) != 0.0 ||
        (z - rotation_centre_z) != 0.0) {
      done = smooth_scroll_maybe_sinusoidal_acceleration(x,y,z,do_zoom_and_move_flag, target_zoom);
   }

   return done;
}

#include <glm/gtx/string_cast.hpp>

#ifndef EMSCRIPTEN
// static
gboolean
graphics_info_t::smooth_scroll_animation_func(GtkWidget *widget,
                                              GdkFrameClock *frame_clock,
                                              gpointer data) {

   // this is not sinusoidal. The first step is 1.

   float frac = 1.0;
   if (graphics_info_t::smooth_scroll_n_steps > 0)
      frac = 1.0/static_cast<float>(graphics_info_t::smooth_scroll_n_steps);
   smooth_scroll_current_step += 1;

   if (smooth_scroll_current_step <= smooth_scroll_n_steps) {

      double theta = 2.0 * M_PI * frac * smooth_scroll_current_step; // not used!
      coot::Cartesian this_step_delta = smooth_scroll_delta * frac;
      add_vector_to_rotation_centre(this_step_delta);
      if (false)
         std::cout << "animation this_step " << smooth_scroll_current_step
                   << " this_step_delta: " << this_step_delta
                   << " for frac " << frac
                   << " Rotation centre now " << glm::to_string(get_rotation_centre()) << std::endl;

      graphics_draw(); // adds to the queue

      return G_SOURCE_CONTINUE;
   } else {
      graphics_info_t g;
      g.update_things_on_move_and_redraw();
      return G_SOURCE_REMOVE;
   }
}
#endif

#ifndef EMSCRIPTEN
// static
gboolean
graphics_info_t::smooth_sinusoidal_scroll_animation_func(GtkWidget *widget,
                                                         GdkFrameClock *frame_clock,
                                                         gpointer data) {

   // smooth_scroll_n_steps = 20; // make this user-defined, or use frame_clock

   smooth_scroll_current_step++;
   if (smooth_scroll_current_step <= smooth_scroll_n_steps) {
      double frac_now  = static_cast<double>(smooth_scroll_current_step  )/static_cast<double>(smooth_scroll_n_steps);
      double frac_next = static_cast<double>(smooth_scroll_current_step+1)/static_cast<double>(smooth_scroll_n_steps);
      double theta_now  = M_PI * frac_now;
      double theta_next = M_PI * frac_next;
      double cos_theta_now  = cos(theta_now);
      double cos_theta_next = cos(theta_next);
      double fp_1 = 0.5 * (1.0 - cos_theta_now);
      double fp_2 = 0.5 * (1.0 - cos_theta_next);
      coot::Cartesian full_delta = smooth_scroll_target_point - smooth_scroll_start_point;
      coot::Cartesian delta_path = full_delta * fp_2 - full_delta * fp_1;
      add_vector_to_rotation_centre(delta_path);
      graphics_draw();
      return G_SOURCE_CONTINUE;
   } else {
      // finished moving
      graphics_info_t g;
      g.update_things_on_move_and_redraw();
      g.update_environment_distances_by_rotation_centre_maybe(g.go_to_atom_molecule());
      return G_SOURCE_REMOVE;
   }
}
#endif

#ifndef EMSCRIPTEN
bool
graphics_info_t::smooth_scroll_maybe_sinusoidal_acceleration(float x, float y, float z,
                                                             short int do_zoom_and_move_flag,
                                                             float target_zoom) {

   bool done_the_move = false; // well, "set it up to go" to be more accurate

   // This is more like how PyMOL does it (and is better than stepped
   // acceleration).

   // acceleration between istep 0 and smooth_scroll_steps (n_steps) is:
   //
   // acc = sin(istep/nsteps * 2 pi)
   //
   // v   = -cos(istep/nsteps * pi)

   // for theta (0->2pi) for frac (0,1)
   // how about acc = sin(theta) + pi * 0.032 * sin(3 * theta)
   // so v  = ? -cos(theta) + pi * 0.032 * 3 * -cos(3*theta)

   float xd = x - rotation_centre_x;
   float yd = y - rotation_centre_y;
   float zd = z - rotation_centre_z;

   smooth_scroll_start_point  = get_rotation_centre_cart();
   smooth_scroll_target_point = coot::Cartesian(x,y,z);

   if (false)
      std::cout << "debug:: in smooth_scroll_maybe_sinusoidal_acceleration "
                << "current centre " << X() << " " << Y() << " " << Z()
                << " delta from target " << xd << " " << yd << " " << zd << std::endl;
   if ( (xd*xd + yd*yd + zd*zd) < smooth_scroll_limit*smooth_scroll_limit ) {

      float pre_zoom = zoom;

      float frac = 1;
      if (smooth_scroll_n_steps > 0)
         frac = 1/float (smooth_scroll_n_steps);
      float stepping_x = frac*xd;
      float stepping_y = frac*yd;
      float stepping_z = frac*zd;

      float rc_x_start = rotation_centre_x;
      float rc_y_start = rotation_centre_y;
      float rc_z_start = rotation_centre_z;

      smooth_scroll_on = 1; // flag to stop wirecube being drawn.
      double v_acc = 0; // accumulated distance
      gpointer user_data = 0;
      smooth_scroll_current_step = -1; // first thing the function does is add 1 to smooth_scroll_current_step.
      smooth_scroll_delta = coot::Cartesian(xd, yd, zd);
      if (false) {
         std::cout << "in smooth_scroll_maybe_sinusoidal_acceleration() with set smooth_scroll_delta "
                   << smooth_scroll_delta << " length " << smooth_scroll_delta.amplitude()
                   << std::endl;

         std::cout << "debug:: in smooth_scroll_maybe_sinusoidal_acceleration "
                   << "current centre " << X() << " " << Y() << " " << Z() << " about to add tick"
                   << std::endl;
      }

      // put this in glarea_tick_func() ?
      //
      gtk_widget_add_tick_callback(glareas[0], smooth_sinusoidal_scroll_animation_func, user_data, NULL);
      done_the_move = true;

      // restore state
      smooth_scroll_on = 0;
   }
   return done_the_move;
}
#endif

#ifndef EMSCRIPTEN
void
graphics_info_t::smooth_scroll_maybe_stepped_acceleration(float x, float y, float z,
     short int do_zoom_and_move_flag,
     float target_zoom) {
          // defunct
}
#endif

std::vector<int>
graphics_info_t::displayed_map_imols() const {

   std::vector<int> is;
   for (int i=0; i<n_molecules(); i++) {
      if (molecules[i].has_xmap()) {
    if (molecules[i].is_displayed_p()) {
       is.push_back(i);
    }
      }
   }
   return is;
}

// for CFC
void
graphics_info_t::display_all_model_molecules() {

   int n = n_molecules();

   for (int i=0; i<n; i++) {
      int state = 1;
      if (is_valid_model_molecule(i)) {
         molecules[i].set_mol_is_displayed(state);
         set_display_control_button_state(i, "Displayed", state);
      }
   }
}

// and unactivate
void
graphics_info_t::undisplay_all_model_molecules_except(int imol) {

   int n = n_molecules();

   for (int i=0; i<n; i++) {
      int state = 0;
      if (i == imol)
         state = 1;
      if (is_valid_model_molecule(i)) {
         molecules[i].set_mol_is_displayed(state); // raw, no callbacks
         molecules[i].set_mol_is_active(state);    //
         set_display_control_button_state(imol, "Displayed", state);
         set_display_control_button_state(imol, "Active",   state);
      }
   }
}

void
graphics_info_t::undisplay_all_model_molecules_except(const std::vector<int> &keep_these) {

   int n = n_molecules();

   for (int i=0; i<n; i++) {
      int state = 0;
      bool found_in_keep_these = false;
      for (unsigned int j=0; j<keep_these.size(); j++) {
         if (keep_these[j] == i) {
            found_in_keep_these = true;
            break;
         }
      }
      if (found_in_keep_these)
         state = 1;
      if (is_valid_model_molecule(i)) {
         molecules[i].set_mol_is_displayed(state);
         molecules[i].set_mol_is_active(state);
         set_display_control_button_state(i, "Displayed", state);
         set_display_control_button_state(i, "Active", state);
      }
   }
}


void
graphics_info_t::setRotationCentreSimple(const coot::Cartesian &c) {

   rotation_centre_x = c.get_x();
   rotation_centre_y = c.get_y();
   rotation_centre_z = c.get_z();

}

// return true if this function did the (simple and immediate) jump to the centre
// force_jump is default false
bool
graphics_info_t::setRotationCentre(coot::Cartesian new_centre, bool force_jump) {

   class pulse_data_t {
   public:
      int n_pulse_steps;
      int n_pulse_steps_max;
      pulse_data_t(int n1, int n2) {
         n_pulse_steps = n1;
         n_pulse_steps_max = n2;
      }
   };

   auto centre_identification_pulse = [] (coot::Cartesian current_centre) {

      auto identification_pulse_func = [] (GtkWidget *widget,
                                           GdkFrameClock *frame_clock,
                                           gpointer data) {

                           gboolean continue_status = 1;
                           pulse_data_t *pulse_data = reinterpret_cast<pulse_data_t *>(data);
                           pulse_data->n_pulse_steps += 1;
                           if (pulse_data->n_pulse_steps > pulse_data->n_pulse_steps_max) {
                              continue_status = 0;
                              lines_mesh_for_identification_pulse.clear();
                           } else {
                              float ns = pulse_data->n_pulse_steps;
                              lines_mesh_for_identification_pulse.update_buffers_for_pulse(ns);
                           }
                           graphics_draw();
                           return gboolean(continue_status);
      };

      // Here I need to check that there isn't already a pulse running! (e.g. from atom selection pulse)
      // (happens when mouse double clicked)
      if (lines_mesh_for_identification_pulse.empty()) {
         if (generic_pulse_centres.empty()) {
            pulse_data_t *pulse_data = new pulse_data_t(0, 30);
            gpointer user_data = reinterpret_cast<void *>(pulse_data);
            identification_pulse_centre = cartesian_to_glm(current_centre);
            gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0]));
            bool broken_line_mode = true;
            lines_mesh_for_identification_pulse.setup_green_pulse(broken_line_mode);
            gtk_widget_add_tick_callback(glareas[0], identification_pulse_func, user_data, NULL);
         }
      }
   };

   // std::cout << "----------- in setRotationCentre() " << std::endl;
   bool needs_centre_jump = true;

   coot::Cartesian current_centre = RotationCentre();
   set_old_rotation_centre(current_centre);

   if (! use_graphics_interface_flag) {
      setRotationCentreSimple(new_centre);
      return true;
   }

   // smooth_scroll_maybe

   // smooth_scroll_maybe now sets up a timeout to move
   // to the given location, so we can't just jump
   // to it here as we used to at the end of this function,
   // because the timeout function  will use "centre"
   // as the place from which to start moving :-)

   // If we are already here then don't animate a move.
   //
   bool already_here = false;
   coot::Cartesian position_delta = new_centre - current_centre;
   if (position_delta.amplitude() < 0.3) {

      centre_identification_pulse(current_centre);

      already_here = true;
      needs_centre_jump = false;
   }

   if (!already_here) {
      if (force_jump) {
         setRotationCentreSimple(new_centre);
         run_post_set_rotation_centre_hook();
      } else {
         if (graphics_info_t::smooth_scroll == 1) {
            // don't zoom and dummy value
            bool status = smooth_scroll_maybe(new_centre.x(), new_centre.y(), new_centre.z(), 0, 100.0);
            if (status) needs_centre_jump = false; // all in hand
         }

         if (needs_centre_jump)  {
            setRotationCentreSimple(new_centre);
            run_post_set_rotation_centre_hook();
         }
      }
   }

   return needs_centre_jump;
}

void
graphics_info_t::setRotationCentreAndZoom(coot::Cartesian centre,
                                          float target_zoom) {

   set_old_rotation_centre(RotationCentre());

   if (graphics_info_t::smooth_scroll == 1)
      smooth_scroll_maybe(centre.x(), centre.y(), centre.z(), 1, target_zoom);

   rotation_centre_x = centre.get_x();
   rotation_centre_y = centre.get_y();
   rotation_centre_z = centre.get_z();
   zoom = target_zoom;
   run_post_set_rotation_centre_hook();
}




void
graphics_info_t::ShowFPS(){

   std::cout << "............. in ShowFPS()" << std::endl;

   long t = 0;

   t = 0; // glutGet(GLUT_ELAPSED_TIME);
   std::cout << "Fix timer in ShowFPS()\n";
   if (t - graphics_info_t::T0 >= 5000) {
      GLfloat seconds = (t-T0)/1000.0;
      GLfloat fps = GLfloat (Frames)/seconds;

      std::string s = "INFO:: ";
      s += int_to_string(Frames);
      s += " frames in ";
      s += float_to_string(seconds);
      s += " seconds = ";
      s += float_to_string(fps);
      s += " frames/sec";

      graphics_info_t g;
      g.add_status_bar_text(s);
      std::cout << s << std::endl;
      graphics_info_t::T0=t;
      graphics_info_t::Frames=0;
   }
}

// We need to reset the Frames so that the first time we get a FPS
// response we are not including all those frames that were made
// without the timer being on.
//
void
graphics_info_t::SetShowFPS(int t) {

#ifndef EMSCRIPTEN
   show_fps_flag = t;
   Frames = 0;
   if (t == 0) {
      // turn it off
      do_tick_constant_draw = false;
   } else {
      // turn it on
      if (! tick_function_is_active()) {
         int new_tick_id = gtk_widget_add_tick_callback(glareas[0], glarea_tick_func, 0, 0);
         idle_function_spin_rock_token = new_tick_id;  // fix this name!
      }
      do_tick_constant_draw = true;
   }
#endif
}

//
void
graphics_info_t::SetActiveMapDrag(int t) {

   active_map_drag_flag = t;

}


void
graphics_info_t::set_font_size(int size) {

//    cout << "graphics_info_t::setting atom_label_font_size to "
// 	<< size << endl;

   atom_label_font_size = size;

#if 0
   if (size == 2) {
      atom_label_font = GLUT_BITMAP_HELVETICA_12;
   } else {
      if (size < 2) {
      atom_label_font = GLUT_BITMAP_HELVETICA_10;
      } else {
   if (size == 3) {
     atom_label_font = GLUT_BITMAP_HELVETICA_18;
   } else {
     // no all other fonts
     if (size == 4) {
       atom_label_font = GLUT_BITMAP_TIMES_ROMAN_10;
     } else if (size == 5) {
       atom_label_font = GLUT_BITMAP_TIMES_ROMAN_24;
     } else if (size == 6) {
       atom_label_font = GLUT_BITMAP_8_BY_13;
     } else if (size == 7) {
       atom_label_font = GLUT_BITMAP_9_BY_15;
     } else {
       // somethign above 7 -> reset default
       atom_label_font = GLUT_BITMAP_HELVETICA_12;
     }
   }
      }
   }
#endif
   // make the labels (if there are any) change now

   graphics_draw();
}


// For the bin?
void
graphics_info_t::update_map_colour_menu()
{
   for (int ii=0; ii<n_molecules(); ii++)
      molecules[ii].update_map_colour_menu_maybe(ii);
}


// virtual trackball
void
graphics_info_t::set_vt_surface(int v){

   if (v == 1) { // VT_FLAT
      //
      trackball_size = 8.8;

   } else {  // VT_SPHERICAL
      trackball_size = 0.8; // as it was in the original code
   }
}

int
graphics_info_t::vt_surface_status() const {
   int status = 2;
   if (trackball_size > 5)
      status = 1;
   return status;
}

// phs reading
//
std::string
graphics_info_t::get_phs_filename() const {

   return phs_filename;

}

void
graphics_info_t::set_phs_filename(std::string filename) {

   phs_filename = filename;

}




// static
void
graphics_info_t::skeletonize_map(int imol, short int prune_it) {

   graphics_info_t g;

   if (is_valid_map_molecule(imol)) {
      // so that we don't do this when the skeleton is on already:
      //
      if (g.molecules[imol].fc_skeleton_draw_on == 0) {
         g.molecules[imol].fc_skeleton_draw_on = 1;

         //       mean_and_variance<float> mv =
         // 	 map_density_distribution(g.molecules[imol].xmap,0);

         clipper::Map_stats stats(g.molecules[imol].xmap);

         // std::cout << "INFO:: Mean and sigma of map: " << stats.mean() << " and "
         //           << stats.std_dev() << std::endl;
         logger.log(log_t::INFO, "Mean and sigma of map:", stats.mean(), "and", stats.std_dev());

         float map_cutoff = stats.mean() + 1.5*stats.std_dev();
         g.skeleton_level = map_cutoff;

         // derived from sktest:
         //
         g.molecules[imol].xskel_cowtan.init(g.molecules[imol].xmap.spacegroup(),
                                             g.molecules[imol].xmap.cell(),
                                             g.molecules[imol].xmap.grid_sampling());

         // std::cout << "INFO:: making skeleton cowtan..." << std::endl;
         logger.log(log_t::INFO, "making skeleton cowtan...");
         GraphicalSkel cowtan(g.molecules[imol].xmap,
                              g.molecules[imol].xskel_cowtan); //fill xskel_cowtan

         g.molecules[imol].xskel_is_filled = 1; // TRUE

         // various experiments....

         // cowtan.tip_filter(xmap, &xskl); // tinker with xskel_cowtan

         //cowtan.prune(g.molecules[imol].xmap_list[imap],
         //	 &g.molecules[imol].xskel_cowtan);

         //
         cowtan.Pprune(g.molecules[imol].xmap,
                       &g.molecules[imol].xskel_cowtan,
                       map_cutoff);

         if (prune_it) {
            BuildCas bc(g.molecules[imol].xmap, map_cutoff);

            // mark segments by connectivity
            //
            int nsegments = bc.count_and_mark_segments(g.molecules[imol].xskel_cowtan,
                                                       g.molecules[imol].xmap, map_cutoff);

            // std::cout << "INFO:: There were " << nsegments << " different segments" << std::endl;
            logger.log(log_t::INFO, "There were", nsegments, "different segments");

            bc.transfer_segment_map(&g.molecules[imol].xskel_cowtan);
            g.molecules[imol].set_colour_skeleton_by_segment(); // use random colours

         } else {
            g.molecules[imol].set_colour_skeleton_by_level(); // use conventional
            // colouring, (just
            // sets a flag)
         }


         // now display the skeleton

         g.molecules[imol].update_clipper_skeleton();
         graphics_draw();

      } else {
         std::cout << "This map has a skeleton already" << std::endl;
      }
   }
}

// static
void
graphics_info_t::set_initial_map_for_skeletonize() {

   // Initially map_for_skeletonize is -1;

   if (graphics_info_t::map_for_skeletonize == -1) {
      for (int imol=0; imol<n_molecules();imol++) {
         if (graphics_info_t::molecules[imol].has_xmap()) {
            graphics_info_t::map_for_skeletonize = imol;
            break;
         }
      }
   }
}

// static
void
graphics_info_t::unskeletonize_map(int imol) {

   graphics_info_t g;

   if (imol >= 0) {
      g.molecules[imol].unskeletonize_map();
      graphics_draw();
   } else {
      std::cout << "Map skeleton not selected from optionmenu." << std::endl;
      std::cout << "Please try again and this time, select "
                << "a map from the optionmenu" << std::endl;
   }
}




// Do we need to delete the old regularize_object_bonds_box?
// Yes, well, clear_up() it.
//
void
graphics_info_t::clear_moving_atoms_object() {

   in_edit_chi_mode_flag = 0;
   in_edit_torsion_general_flag = 0;
   in_moving_atoms_drag_atom_mode_flag = 0; // no more dragging atoms
   have_fixed_points_sheared_drag_flag = 0;
   // and set the rotation translation atom index to unknown again:
   // rot_trans_atom_index_rotation_origin_atom = -1;
   rot_trans_rotation_origin_atom = NULL;

   // std::cout << "clearing intermediate object..." << std::endl;
   graphical_bonds_container empty_box;
   regularize_object_bonds_box.clear_up();
   regularize_object_bonds_box = empty_box;

   dynamic_distances.clear();

   graphics_draw();
}

void
graphics_info_t::set_refinement_map(int i) {
   imol_refinement_map = i;
}



// Recall that moving_atoms_asc contains an atom_selection_container_t
// with the new coordinates in.  the mol contains all the molecule and
// the atom_selection contains just the moving parts.
//
// We want to put the moving parts back into the object that we were
// regularizing (imol_moving_atoms).
//
coot::refinement_results_t
graphics_info_t::accept_moving_atoms() {

   auto debug_moving_atoms = [] () {
      std::cout << "::::::::: debug_moving_atoms() moving_atoms_asc:" << moving_atoms_asc << std::endl;
      if (! moving_atoms_asc) {
         std::cout << "ERROR:: null moving_atoms_asc in accept_moving_atoms() debug_moving_atoms()" << std::endl;
         return;
      }
      std::cout << "::::::::: debug_moving_atoms() moving_atoms_asc mol: " << moving_atoms_asc->mol << std::endl;
      mmdb::Manager *mol = moving_atoms_asc->mol;
      if (! mol) {
         std::cout << "ERROR:: null moving_atoms_asc mol in accept_moving_atoms() " << std::endl;
         return;
      }
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
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
                        std::cout << "moving atom: " << at << " " << coot::atom_spec_t(at) << std::endl;
                     }
                  }
               }
            }
         }
      }
   };

   coot::refinement_results_t rr;

   while (continue_threaded_refinement_loop) {
      // std::cout << "waiting for continue_threaded_refinement_loop to be false..." << std::endl;
      std::this_thread::sleep_for(std::chrono::milliseconds(200));
   }

   bool debug = false;

   if (debug) {
      std::cout << ":::: INFO:: accept_moving_atoms() imol moving atoms is " << imol_moving_atoms
                << std::endl;
      std::cout << ":::: INFO:: accept_moving_atoms() imol moving atoms type is "
                << moving_atoms_asc_type << " vs " << coot::NEW_COORDS_REPLACE << std::endl;
      debug_moving_atoms();
   }

   if (! moving_atoms_asc) {
      std::cout << "ERROR:: null moving_atoms_asc in accept_moving_atoms() " << std::endl;
      return rr;
   }
   if (! moving_atoms_asc->mol) {
      std::cout << "ERROR:: null moving_atoms_asc mol in accept_moving_atoms() " << std::endl;
      return rr;
   }

   rr = get_refinement_results();

   if (moving_atoms_asc_type == coot::NEW_COORDS_ADD) { // not used!
      molecules[imol_moving_atoms].add_coords(*moving_atoms_asc);
   } else {
      bool mzo = refinement_move_atoms_with_zero_occupancy_flag;
      if (moving_atoms_asc_type == coot::NEW_COORDS_REPLACE_CHANGE_ALTCONF) {
         molecules[imol_moving_atoms].replace_coords(*moving_atoms_asc, 1, mzo); // doesn't dealloc moving_atoms_asc
         update_validation(imol_moving_atoms);
      } else {
         if (moving_atoms_asc_type == coot::NEW_COORDS_REPLACE) {
            molecules[imol_moving_atoms].replace_coords(*moving_atoms_asc, 0, mzo);
            // update_validation(imol_moving_atoms);
         } else {
            if (moving_atoms_asc_type == coot::NEW_COORDS_INSERT) {
               molecules[imol_moving_atoms].insert_coords(*moving_atoms_asc);
            } else {
               if  (moving_atoms_asc_type == coot::NEW_COORDS_INSERT_CHANGE_ALTCONF) {
                  molecules[imol_moving_atoms].insert_coords_change_altconf(*moving_atoms_asc);
               } else {
                  std::cout << "------------ ERROR! -------------------" << std::endl;
                  std::cout << "       moving_atoms_asc_type not known: ";
                  std::cout << moving_atoms_asc_type << std::endl;
                  std::cout << "------------ ERROR! -------------------" << std::endl;
               }
            }
         }
      }
   }


   // reset the b-factor?
   if (graphics_info_t::reset_b_factor_moved_atoms_flag) {
     molecules[imol_moving_atoms].set_b_factor_atom_selection(*moving_atoms_asc, graphics_info_t::default_new_atoms_b_factor, 1);
   }

   if (do_probe_dots_post_refine_flag) {
      setup_for_probe_dots_on_chis_molprobity(imol_moving_atoms);
   }

   std::cout << "DEBUG:: accept_moving_atoms() GTK4 update: update rama plot for " << imol_moving_atoms << std::endl;

   clear_all_atom_pull_restraints(false); // no re-refine
   clear_up_moving_atoms();
   update_environment_distances_by_rotation_centre_maybe(imol_moving_atoms);

   hide_atom_pull_toolbar_buttons();

   normal_cursor(); // we may have had fleur cursor.
   // and set the rotation translation atom index to unknown again:
   // rot_trans_atom_index_rotation_origin_atom = -1;
   rot_trans_rotation_origin_atom = NULL;
   in_moving_atoms_drag_atom_mode_flag = 0; // no more dragging atoms
   have_fixed_points_sheared_drag_flag = 0;
   in_edit_chi_mode_view_rotate_mode = 0;

   // Why do I want to do this when the intermediate atoms are undisplayed?
   // To *clear up* the interactive probe dots.
   if (do_coot_probe_dots_during_refine_flag)
      do_interactive_coot_probe(); // in accept_moving_atoms() to turn off the dots

   // Hmm... this won't work as expected because the difference map is not updated yet!
   // I need to hook into the end of a difference map update.
   // fill_difference_map_peaks_button_box(); // update the difference map peaks if the dialog is open

   rama_plot_boxes_handle_molecule_update(imol_moving_atoms);
   //    draw_rama_plots(); // 20230526-PE should this be here or elsewhere? Don't rama graphs now
                            //  get drawn in graphics_draw()?

   // 20230527-PE does this belong here? - lets see....
   // update_active_validation_graph_model(imol_moving_atoms);
   update_validation(imol_moving_atoms);

   int mode = MOVINGATOMS;

   run_post_manipulation_hook(imol_moving_atoms, mode);

   if (debug) {
      std::cout << ":::: INFO:: accept_moving_atoms() imol moving atoms finishes " << std::endl;
   }

   return rr;
}

void
graphics_info_t::run_post_read_model_hook(int imol) {

#if 0 // 20220606-PE Python doesn't work yet

   std::string s;

#ifdef USE_GUILE

   s = "post-read-model-hook";
   SCM v = safe_scheme_command(s.c_str());
   // std::cout << "scm v " << v << std::endl;
   if (scm_is_true(scm_procedure_p(v))) {
      s += "(" + s + " " + int_to_string(imol) + ")";
      SCM result = safe_scheme_command(s);
   } else {
   }
#endif

#ifdef USE_PYTHON
   s = "post_read_model_hook";
   PyObject *pName_coot = myPyString_FromString("__main__");  // not "coot" at the moment
   PyObject *pModule = PyImport_Import(pName_coot);
   PyObject *pDict = PyModule_GetDict(pModule);
   PyObject *pFunc = PyDict_GetItemString(pDict, s.c_str());

   if (false) {
      PyObject *keys = PyDict_Keys(pDict);
      unsigned int l = PyObject_Length(keys);
      for (std::size_t i=0; i<l; i++) {
	 PyObject *item = PyList_GetItem(keys, i);
	 std::cout << "key " << i << " " << myPyString_AsString(item) << std::endl;
      }
   }

   if (PyCallable_Check(pFunc)) {
      PyObject *arg_list = PyTuple_New(1);
      PyObject *imol_py = PyLong_FromLong(imol);
      PyTuple_SetItem(arg_list, 0, imol_py);
      PyObject *result_py = PyEval_CallObject(pFunc, arg_list);
      std::cout << "DEBUG:: post_read_model_hook() got result " << result_py << std::endl;
   } else {
      // std::cout << "INFO:: in run_post_read_model_hook() pFunc " << pFunc << " is not callable" << std::endl;
      logger.log(log_t::INFO, "in run_post_read_model_hook() pFunc", pFunc, "is not callable");
      // std::cout << "INFO:: in run_post_read_model_hook() pDict " << pDict << " " << std::endl;
      logger.log(log_t::INFO, "in run_post_read_model_hook() pDict", pDict);
      // std::cout << "INFO:: in run_post_read_model_hook() pModule " << pModule << " " << std::endl;
      logger.log(log_t::INFO, "in run_post_read_model_hook() pModule", pModule);
   }
#endif

#endif

}

void
graphics_info_t::run_post_manipulation_hook(int imol, int mode) {


#ifdef USE_GUILE
   run_post_manipulation_hook_scm(imol, mode);
#endif // GUILE
#ifdef USE_PYTHON
   // turn this off for the moment.
   run_post_manipulation_hook_py(imol, mode);
#endif
}

#ifdef USE_GUILE
void
graphics_info_t::run_post_manipulation_hook_scm(int imol, int mode) {

   std::string pms = "post-manipulation-hook";
   SCM v = safe_scheme_command(pms);

   if (scm_is_true(scm_procedure_p(v))) {
      std::string ss = "(";
      ss += pms;
      ss += " ";
      ss += int_to_string(imol);
      ss += " ";
      ss += int_to_string(mode);
      ss += ")";
      SCM res = safe_scheme_command(ss);
      SCM dest = SCM_BOOL_F;
      SCM mess =  scm_from_locale_string("result: ~s\n");
      SCM p = scm_simple_format(dest, mess, scm_list_1(res));
      std::cout << scm_to_locale_string(p);
   }
}
#endif


#ifdef USE_PYTHON
void
graphics_info_t::run_post_manipulation_hook_py(int imol, int mode) {

   // 20230527-PE exiciting dangerous times - turning this on again:
   // std::cout << "FIXME:: ----- due to python setup problems not running run_post_manipulation_hook_py()"
   //           << std::endl;
   // return;

   if (false) // debugging
      std::cout << "debug run_post_manipulation_hook_py() --- start --- " << std::endl;

   PyObject *error_thing = PyErr_Occurred();
   if (error_thing) {
      std::cout << "ERROR:: while executing run_post_manipulation_hook_py() a python error occured before start " << std::endl;
      PyObject *type, *value, *traceback;
      PyErr_Fetch(&type, &value, &traceback);
      PyErr_NormalizeException(&type, &value, &traceback);
      PyObject *exception_string = PyObject_Repr(value);
      const char *em = myPyString_AsString(exception_string);
      std::cout << "ERROR:: " << em << std::endl;
      Py_XDECREF(value);
      Py_XDECREF(traceback);
      Py_XDECREF(type);
   }

   PyObject *result = nullptr;
   PyObject *pModule = PyImport_ImportModule("coot_utils");
   if (!pModule) {
      std::cout << "ERROR:: run_post_manipulation_hook_py() no pModule " << std::endl;
      PyErr_Print();
      return;
   } else {
      PyObject *attr = PyObject_GetAttrString(pModule, "post_manipulation_script");
      if (attr && PyCallable_Check(attr)) {

         if (false) { // debugging calling this function
            // attr is a new/borrowed reference you already obtained
            PyObject *repr = PyObject_Repr(attr);  // new reference
            if (repr) {
               const char *s = PyUnicode_AsUTF8(repr); // NULL if not a unicode
               if (s)
                  std::cout << "callable repr: " << s << std::endl;
               else
                  std::cout << "callable repr: <non-unicode repr>" << std::endl;
               Py_XDECREF(repr);
            } else {
               // If PyObject_Repr fails, print the Python error so you can debug
               PyErr_Print();
            }
         }

         PyObject *arg_1 = PyLong_FromLong(imol);
         PyObject *arg_2 = PyLong_FromLong(mode);
         // result = PyObject_CallObject(attr, t);
         result = PyObject_CallFunctionObjArgs(attr, arg_1, arg_2, NULL);
         Py_XDECREF(result);
      } else {
         if (false) // useful when debugging.
            std::cout << "DEBUG:: coot_utils.post_manipulation_script is not callable" << std::endl;
         return;
      }
      Py_XDECREF(attr);
      Py_XDECREF(pModule);
   }

   // the above function can set an error  - that's bad news for the python wrapping
   // of accept_moving_atoms(). So instead of properly handling the error, or investigating
   // why it is happening, let's just clear it.
   //
   error_thing = PyErr_Occurred();
   if (! error_thing) {
      if (false)
         std::cout << "DEBUG:: run_post_manipulation_hook_py() No Python error on callable check" << std::endl;
   } else {
      std::cout << "ERROR:: while executing run_post_manipulation_hook_py() a python error occured " << std::endl;
      PyObject *type, *value, *traceback;
      PyErr_Fetch(&type, &value, &traceback);
      PyErr_NormalizeException(&type, &value, &traceback);
      PyObject *exception_string = PyObject_Repr(value);
      const char *em = myPyString_AsString(exception_string);
      std::cout << "ERROR:: " << em << std::endl;
      Py_XDECREF(value);
      Py_XDECREF(traceback);
      Py_XDECREF(type);
   }

   PyErr_Clear();

   // Py_XDECREF(v);
}
#endif


void
graphics_info_t::run_post_set_rotation_centre_hook() {

   // Don't run the post update hook at the moment.
  // Python is not wired up (correctly)

#if defined USE_GUILE
   // run_post_set_rotation_centre_hook_scm();
#endif // GUILE

#ifdef USE_PYTHON
   // run_post_set_rotation_centre_hook_py();
#endif

}

#ifdef USE_GUILE
void
graphics_info_t::run_post_set_rotation_centre_hook_scm() {
   std::string s = "post-set-rotation-centre-hook";
   SCM v = safe_scheme_command(s);

   if (scm_is_true(scm_procedure_p(v))) {
      std::string ss = "(";
      ss += s;
      ss += ")";
      SCM res = safe_scheme_command(ss);
      if (false) {  // too noisy
    SCM dest = SCM_BOOL_F;
    SCM mess = scm_from_locale_string("result: ~s\n");
    SCM p = scm_simple_format(dest, mess, scm_list_1(res));
    std::cout << scm_to_locale_string(p);
      }
   }
}
#endif

#ifdef USE_PYTHON
void
graphics_info_t::run_post_set_rotation_centre_hook_py() {

      // Same as for the post manipulation "hook", i.e. script (maybe
      // could become an extra function then...
      // BL says:: we can do it all in python API or use the 'lazy' method
      // and check in the python layer (which we will do...)
      PyObject *v;
      int ret;
      std::string ps = "post_set_rotation_centre_script";
      std::string check_ps = "callable(";
      check_ps += ps;
      check_ps += ")";
      v = safe_python_command_with_return(check_ps);
      ret = PyLong_AsLong(v);
      if (ret == 1) {
        std::string ss = ps;
        ss += "()";
        PyObject *res = safe_python_command_with_return(ss);
        PyObject *fmt =  myPyString_FromString("result: \%s");
        PyObject *tuple = PyTuple_New(1);
        PyTuple_SetItem(tuple, 0, res);
        PyObject *msg = PyUnicode_Format(fmt, tuple);

        std::cout << PyUnicode_AsUTF8String(msg)<<std::endl;;
        Py_DECREF(msg);
      }
      Py_XDECREF(v);
}
#endif


void
graphics_info_t::pull_restraint_neighbour_displacement_change_max_radius(bool up_or_down) {

   if (last_restraints) {
      if (up_or_down)
         pull_restraint_neighbour_displacement_max_radius -= 1.0;
      else
         pull_restraint_neighbour_displacement_max_radius += 1.0;

      if (pull_restraint_neighbour_displacement_max_radius < 0.0)
         pull_restraint_neighbour_displacement_max_radius = 0.0;

      // std::cout << "debug:: pull_restraint_neighbour_displacement_max_radius "
      // << pull_restraint_neighbour_displacement_max_radius << std::endl;

      float r = pull_restraint_neighbour_displacement_max_radius;
      attach_buffers(); // because we touch some GL buffers
      lines_mesh_for_pull_restraint_neighbour_displacement_max_radius_ring.update_radius_ring_vertices(r);

      if (pull_restraint_neighbour_displacement_max_radius > 1.99) {
         last_restraints->set_use_proportional_editing(true);
         last_restraints->pull_restraint_neighbour_displacement_max_radius =
            pull_restraint_neighbour_displacement_max_radius;
      } else {
         last_restraints->set_use_proportional_editing(false);
      }
   }
}



void
graphics_info_t::update_environment_distances_by_rotation_centre_maybe(int imol_in) {

   // Oh this is grimly "long hand".
   graphics_info_t g;
   if (g.environment_show_distances) {
      coot::at_dist_info_t at_d_i = g.molecules[imol_in].closest_atom(RotationCentre());
      if (at_d_i.atom) {
    int atom_index;
    if (at_d_i.atom->GetUDData(g.molecules[imol_in].atom_sel.UDDAtomIndexHandle,
       atom_index) == mmdb::UDDATA_Ok) {
       g.mol_no_for_environment_distances = imol_in;
       g.update_environment_distances_maybe(atom_index, imol_in);
    }
      }
   }
}


void
graphics_info_t::clear_up_moving_atoms() {

   // this function is not just graphics, so it needs to check use_graphics_interface_flag
   // internally.

   // Note to self: why don't I do a delete moving_atoms_asc somewhere here?
   // Where does the moving_atoms_asc->mol go?

   moving_atoms_asc_type = coot::NEW_COORDS_UNSET; // unset
   in_moving_atoms_drag_atom_mode_flag = 0; // no more dragging atoms
   have_fixed_points_sheared_drag_flag = 0;
   // and take out any drag refine idle function:

   // We can't clear up until we have the lock and the lock here means that refinement
   // can't continue (lock is restraints_lock).

   // bool compare_exchange_weak (T& expected, T desired);
   // when "expected" equals the value of restraint_lock, compare_exchange_weak() returns
   // true and replaces "desired" as the value of restraints lock - atomically of course.
   // If the values are not equal, then we want to sleep/wait for a bit - hence
   // enter the while loop.
   // Dealing with spurious failure by using '&& !unlocked' in the while test:
   // it seems that the test is always true and we never enter the while loop
   // and wait - even if restraints_lock is true when we start.

   // is this useful?
   // std::this_thread::sleep_for(std::chrono::milliseconds(100));

   get_restraints_lock(__FUNCTION__);

   // We must not delete the moving atoms if they are being used to manipulate pull restraints
   //
   bool unlocked_atoms = false;
   while (! moving_atoms_lock.compare_exchange_weak(unlocked_atoms, true) && !unlocked_atoms) {
      std::this_thread::sleep_for(std::chrono::milliseconds(1));
      unlocked_atoms = false;
   }

   moving_atoms_visited_residues.clear();
   continue_update_refinement_atoms_flag = false;
   continue_threaded_refinement_loop = false;
   if (moving_atoms_asc) {

      if (moving_atoms_asc->atom_selection != NULL) {
         if (moving_atoms_asc->n_selected_atoms > 0) {
            moving_atoms_asc->mol->DeleteSelection(moving_atoms_asc->SelectionHandle);
            moving_atoms_asc->atom_selection = NULL;
         } else {
            std::cout << "WARNING:: attempting to delete non-NULL ";
            std::cout << "moving_atoms_asc.atom_selection" << std::endl;
            std::cout << "but moving_atoms_asc.n_selected_atoms == 0" << std::endl;
            std::cout << "ignoring " << std::endl;
         }
      } else {
         // std::cout << "WARNING:: Ignoring attempt to delete NULL moving_atoms_asc.atom_selection"
         // << std::endl;
      }
      if (moving_atoms_asc->mol != NULL) {
         if (moving_atoms_asc->n_selected_atoms > 0) {
            moving_atoms_asc->mol = NULL;
         } else {
            std::cout << "WARNING:: attempting to delete non-NULL moving_atoms_asc.mol" << std::endl;
            std::cout << "but moving_atoms_asc.n_selected_atoms == 0" << std::endl;
            std::cout << "ignoring " << std::endl;
         }
      } else {
         if (false) {
            std::cout << "attempting to delete NULL moving_atoms_asc.mol" << std::endl;
            std::cout << "ignoring " << std::endl;
         }
      }
      moving_atoms_asc->n_selected_atoms = 0;
   }

   dynamic_distances.clear();


   // std::cout << "------------------------ clear_up_moving_atoms(): setting moving_atoms_asc to null" << std::endl;
   // and now the signal that moving_atoms_asc has been cleared:
   //
   moving_atoms_asc = NULL; // 20200412-PE. Why was this not done years ago?
                            // Suspicious.
                            // OK, so let the test be on moving_atoms_asc->mol
                            // moving_atoms_asc is set in init()

   if (last_restraints) {
      last_restraints->clear();
      delete last_restraints;
      last_restraints = 0;
      unset_moving_atoms_currently_dragged_atom_index();
   }

   release_restraints_lock(__FUNCTION__); // refinement ended and cleared up.

   moving_atoms_lock  = false;

   // std::cout << "calling rebond_molecule_corresponding_to_moving_atoms() " << std::endl;
   // 20220220-PE I will comment this out (because I think the answer to the below question is "yes"
   // graphics_info_t::rebond_molecule_corresponding_to_moving_atoms(); // haven't we done this?

   if (use_graphics_interface_flag) {

      draw_gl_ramachandran_plot_flag = false;

      hydrogen_bonds_atom_position_pairs.clear();
      update_hydrogen_bond_mesh("");

      // now the diegos
      bad_nbc_atom_pair_marker_positions.clear(); // this should be in the update function, surely?
      update_bad_nbc_atom_pair_marker_positions();
      update_bad_nbc_atom_pair_dashed_lines();
      update_chiral_volume_outlier_marker_positions();

   }
}


// if the imol for moving atoms is imol, delete the moving atoms (called from close_molecule)
void
graphics_info_t::clear_up_moving_atoms_maybe(int imol) {

   // clear up moving atoms for this molecule if they exist for this given molecule

   if (imol_moving_atoms == imol) {
      if (moving_atoms_asc) {
         if (moving_atoms_asc->n_selected_atoms > 0) {
            clear_up_moving_atoms();
            clear_moving_atoms_object();
         }
      }
   }
}



void
graphics_info_t::set_dynarama_is_displayed(GtkWidget *dyna_toplev, int imol) {

#ifdef HAVE_GOOCANVAS

   // first delete the old plot for this molecule (if it exists)
   //
   if (is_valid_model_molecule(imol)) {

      // Clear out the old one if it was there.
      GtkWidget *w = coot::get_validation_graph(imol, coot::RAMACHANDRAN_PLOT);
      if (w) {
         coot::rama_plot *plot = (coot::rama_plot *) g_object_get_data(G_OBJECT(w), "rama_plot");
         delete plot;
      }
      coot::set_validation_graph(imol, coot::RAMACHANDRAN_PLOT, dyna_toplev);
   } else {
      std::cout << "DEBUG:: in graphics_info_t::set_dynarama_is_displayed() imol " << imol
                << " is not valid" << std::endl;
   }
#endif
}

// Not used.
void
graphics_info_t::delete_molecule_from_display_manager(int imol, bool was_map) {

   if (! use_graphics_interface_flag) return;

   // 20230427-PE these should be "display_control_model_vbox" and "display_control_map_vbox"
   //             to be more clear.
   GtkWidget *vbox_for_molecules = widget_from_builder("display_molecule_vbox");
   if (was_map)
      vbox_for_molecules = widget_from_builder("display_map_vbox");

   for (GtkWidget* child = gtk_widget_get_first_child(vbox_for_molecules);
       child != nullptr;
       child = gtk_widget_get_next_sibling(child)) {

      // child is the molecule hbox
      int imol_for_child = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(child), "imol"));

      if (imol_for_child == imol) {
         gtk_box_remove(GTK_BOX(vbox_for_molecules), child);
         break; // (child is no longer valid - so calling gtk_widget_get_next_sibling() is an error
      }
   }

}


// As these are moving atoms, they cannot point to atoms of a real
// molecule in the molecules array.  [Otherwise as we move moving
// atoms, the static molecule's atoms move (but note that their bonds
// are not updated)].
//
// 20230212-PE This function changes moving_atoms_asc!
//
void
graphics_info_t::make_moving_atoms_graphics_object(int imol,
                                                   const atom_selection_container_t &asc,
                                                   unsigned int do_rama_markup_in, // default MOVING_ATOMS_DO_RAMA_MARKUP_USE_INTERNAL_SETTING,
                                                   unsigned int do_rota_markup_in  // default MOVING_ATOMS_DO_ROTA_MARKUP_USE_INTERNAL_SETTING
                                                   ) {

   if (! moving_atoms_asc) {
      if (false)
         std::cout << "info:: make_moving_atoms_graphics_object() makes a new moving_atoms_asc" << std::endl;
      moving_atoms_asc = new atom_selection_container_t;
   } else {
      // moving_atoms_asc->clear_up(); // crash.  Much complexity to fix the crash, I think.
      // i.e. the clear_up() should be here - I think the problem lies elsewhere.
      // Not clearing up here produces a memory leak, I think (not a bad one (for some reason!)).
   }

   // 20230212-PE this needs to happen in --no-graphics mode
   *moving_atoms_asc = asc;

   if (! use_graphics_interface_flag) return;

   // --------------- also do the restraints -----------------------
   //

   make_moving_atoms_restraints_graphics_object();

   int do_disulphide_flag = 0;

   //    Needed that to debug doubly clear_up on a bonds box.  Now points
   //    reset in clear_up().

   if (false)
      std::cout << "DEBUG:: make_moving_atoms_graphics_object() bonds box type of molecule "
                << imol_moving_atoms << " is " << molecules[imol_moving_atoms].Bonds_box_type()
                << std::endl;

   bool do_ca_mode = false;

   if (molecules[imol_moving_atoms].Bonds_box_type() == coot::CA_BONDS ||
       molecules[imol_moving_atoms].Bonds_box_type() == coot::CA_BONDS_PLUS_LIGANDS ||
       molecules[imol_moving_atoms].Bonds_box_type() == coot::CA_BONDS_PLUS_LIGANDS_AND_SIDECHAINS ||
       molecules[imol_moving_atoms].Bonds_box_type() == coot::CA_BONDS_PLUS_LIGANDS_SEC_STRUCT_COLOUR ||
       molecules[imol_moving_atoms].Bonds_box_type() == coot::COLOUR_BY_RAINBOW_BONDS)
      do_ca_mode = true;

   if (residue_type_selection_was_user_picked_residue_range)
      do_ca_mode = false;

   if (do_ca_mode) {

      if (molecules[imol_moving_atoms].Bonds_box_type() == coot::CA_BONDS_PLUS_LIGANDS) {

         Bond_lines_container bonds;
         bool draw_hydrogens_flag = false;
         if (molecules[imol_moving_atoms].draw_hydrogens())
         draw_hydrogens_flag = true;
         bonds.do_Ca_plus_ligands_bonds(*moving_atoms_asc, imol, Geom_p(), 1.0, 4.7,
                                        draw_missing_loops_flag, draw_hydrogens_flag);

         // 20210725-PE I got a lock-up with one of these locks (not sure which) when playing with
         //             the fun start GM demo model and dragging a CA model around. Something
         //             somewhere else was not unlocking? I think that it may be a pull atom
         //             restraint.

         unsigned int unlocked = 0;
         // Neither of these seems to make a difference re: the intermediate atoms python representation
         // while (! moving_atoms_bonds_lock.compare_exchange_weak(unlocked, 1)) {
         // For now, we don't print the python variable - now convert it to json - that seems
         // to work OK ~20MB/s for a chain.
         //
         while (! moving_atoms_bonds_lock.compare_exchange_weak(unlocked, 1) && !unlocked) {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
            unlocked = 0;
         }

         // we shouldn't draw bonds if they have been deleted in clear_up_moving_atoms():
         if (moving_atoms_asc->atom_selection) {

            // moving_atoms_lock is a bool
            bool unlocked_ma = false;
            while (! moving_atoms_lock.compare_exchange_weak(unlocked_ma, 1) && !unlocked_ma) {
                 std::this_thread::sleep_for(std::chrono::milliseconds(1));
                 unlocked_ma = false;
            }

            regularize_object_bonds_box.clear_up();
            regularize_object_bonds_box = bonds.make_graphical_bonds();


            moving_atoms_lock = 0; // unlock them
         }
         moving_atoms_bonds_lock = 0;

      } else {

         Bond_lines_container bonds;
         bonds.do_Ca_bonds(*moving_atoms_asc, 1.0, 4.7, draw_missing_loops_flag);
         unsigned int unlocked = false;
         while (! moving_atoms_bonds_lock.compare_exchange_weak(unlocked, 1) && !unlocked) {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
            unlocked = 0;
         }
         // we shouldn't draw bonds if they have been deleted in clear_up_moving_atoms():
         if (moving_atoms_asc->atom_selection) {

            // moving_atoms_lock is a bool
            bool unlocked_ma = false;
            while (! moving_atoms_lock.compare_exchange_weak(unlocked_ma, 1) && !unlocked) {
                 std::this_thread::sleep_for(std::chrono::milliseconds(1));
                 unlocked_ma = false;
            }

            regularize_object_bonds_box.clear_up();
            regularize_object_bonds_box = bonds.make_graphical_bonds();
            moving_atoms_lock = 0; // unlock them
         }
         moving_atoms_bonds_lock = 0;
      }

   } else {

      // normal bond representation

      bool do_rama_markup = graphics_info_t::do_intermediate_atoms_rama_markup;
      bool do_rota_markup = graphics_info_t::do_intermediate_atoms_rota_markup;
      // check the args:
      if (do_rama_markup_in == MOVING_ATOMS_DO_RAMA_MARKUP_FALSE) do_rama_markup = false;
      if (do_rama_markup_in == MOVING_ATOMS_DO_RAMA_MARKUP_TRUE ) do_rama_markup = true;
      if (do_rota_markup_in == MOVING_ATOMS_DO_ROTA_MARKUP_FALSE) do_rota_markup = false;
      if (do_rota_markup_in == MOVING_ATOMS_DO_ROTA_MARKUP_TRUE ) do_rota_markup = true;

      // wrap the filling of the rotamer probability tables
      //
      coot::rotamer_probability_tables *tables_pointer = NULL;
      if (do_rota_markup) {
         if (! rot_prob_tables.tried_and_failed()) {
            if (rot_prob_tables.is_well_formatted()) {
               tables_pointer = &rot_prob_tables;
            } else {
               rot_prob_tables.fill_tables();
               if (rot_prob_tables.is_well_formatted()) {
                  tables_pointer = &rot_prob_tables;
               }
            }
         } else {
            do_rota_markup = false;
         }
      }

      int draw_hydrogens_flag = 0;
      if (molecules[imol_moving_atoms].draw_hydrogens())
         draw_hydrogens_flag = 1;
      std::set<int> dummy;
      bool do_sticks_for_waters = false; // 20220508-PE Now we have modern waters.
      bool draw_missing_loops_flag_local = false;
      Bond_lines_container bonds(*moving_atoms_asc, imol_moving_atoms, dummy, Geom_p(),
                                 do_disulphide_flag, draw_hydrogens_flag,
                                 draw_missing_loops_flag_local, 0, "dummy",
                                 do_rama_markup, do_rota_markup, do_sticks_for_waters, tables_pointer);
      unsigned int unlocked = false;
      while (! moving_atoms_bonds_lock.compare_exchange_weak(unlocked, 1) && !unlocked) {
         std::this_thread::sleep_for(std::chrono::milliseconds(1));
         unlocked = 0;
      }
      regularize_object_bonds_box.clear_up();
      regularize_object_bonds_box = bonds.make_graphical_bonds(ramachandrans_container,
                                                               do_rama_markup, do_rama_markup);
      moving_atoms_bonds_lock = 0; // unlocked
   }

   // 20220209-PE adding an atom selection to moving_atoms_molecule. Is this safe?
   // It's needed because I use
   //       int udd_handle_bonded_type = atom_sel.mol->GetUDDHandle(mmdb::UDR_ATOM, "found bond");
   // for the sphere vs hemisphere test
   //
   moving_atoms_molecule.atom_sel = asc;

   moving_atoms_molecule.bonds_box = regularize_object_bonds_box; // needed? or does
                                                                  // make_glsl_bonds_type_checked()
                                                                  // update this?
   float bond_width = molecules[imol].get_bond_thickness();
   moving_atoms_molecule.set_bond_thickness(0.7f * bond_width); // why scaled like this is needed!?. At some
                                                                 // stage remove this and intermediate atom
                                                                 // radius adjustment in make_glsl_bonds_type_checked()
   moving_atoms_molecule.is_intermediate_atoms_molecule = true;

   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0])); // needed?
   shader_for_models.Use();
   moving_atoms_molecule.make_glsl_bonds_type_checked(__FUNCTION__);

   setup_atom_pull_restraints_glsl();

   { // put this somewhere
      std::vector<Instanced_Markup_Mesh_attrib_t> balls;
      update_rama_balls(&balls);
      rama_balls_mesh.update_instancing_buffers(balls);
   }

}

#include "utils/dodec.hh"

// delete this function?
std::vector<coot::old_generic_display_object_t::dodec_t>
graphics_info_t::get_rotamer_dodecs() {

   std::vector<coot::old_generic_display_object_t::dodec_t> dodecs;
   if (regularize_object_bonds_box.num_colours > 0) {
      if (regularize_object_bonds_box.n_rotamer_markups > 0) {
         dodec d;
         glm::vec4 top    = unproject(100, 100, 0.5);
         glm::vec4 bottom = unproject(100,   0, 0.5);
         clipper::Coord_orth screen_y(top.x-bottom.x,
                                      top.y-bottom.y,
                                      top.z-bottom.z);
         for (int i=0; i<regularize_object_bonds_box.n_rotamer_markups; i++) {

            clipper::Coord_orth pos = regularize_object_bonds_box.rotamer_markups[i].pos;
            double size = 0.52;
            pos -= screen_y * double(1.5 * size * 22.0/double(graphics_info_t::zoom));
            coot::old_generic_display_object_t::dodec_t dodec(d, size, pos);
            dodec.col = regularize_object_bonds_box.rotamer_markups[i].col;
            dodecs.push_back(dodec);
         }
      } else {
         std::cout << "WARNING:: in " << __FUNCTION__ << " regularize_object_bonds_box.n_rotamer_markups was 0" << std::endl;
      }
   } else {
      std::cout << "WARNING:: in " << __FUNCTION__ << " regularize_object_bonds_box.num_colours was 0" << std::endl;
   }
   return dodecs;
}


// Merge weirdness
//
// This does (draws) symmetry too.
//
// static
// void
// graphics_info_t::draw_environment_graphics_object() {

//    graphics_info_t g;
//    if (is_valid_model_molecule(mol_no_for_environment_distances)) {
//       if (g.molecules[mol_no_for_environment_distances].is_displayed_p()) {
// 	 g.environment_graphics_object_internal(environment_object_bonds_box);
// 	 if (g.show_symmetry)
// 	    g.environment_graphics_object_internal(symmetry_environment_object_bonds_box);
//       }
//    }
// }


mmdb::Atom *
graphics_info_t::get_moving_atom(const pick_info &pi) const {
   mmdb::Atom *at = 0;
   if (moving_atoms_asc) {
      if (moving_atoms_asc->mol) {
         at = moving_atoms_asc->atom_selection[pi.atom_index];
      }
   }
   return at;
}


// static
void
graphics_info_t::picked_intermediate_atom_graphics_object() {

   if (flash_intermediate_atom_pick_flag) {
      glPointSize(12.0);
      glColor3f(0.99, 0.99, 0.2);
      // for some reason, I have to add the point twice and use
      // GL_POINTS, GL_POINT with one vertex does not work.
      glBegin(GL_POINTS);
      glVertex3f(intermediate_flash_point.x(),
    intermediate_flash_point.y(),
    intermediate_flash_point.z());
      glVertex3f(intermediate_flash_point.x(),
    intermediate_flash_point.y(),
    intermediate_flash_point.z());
      glEnd();
   }
}

// This is the GL rendering of the environment bonds box
//
void
graphics_info_t::environment_graphics_object_internal(const graphical_bonds_container &env_bonds_box) const {

   if (! display_environment_graphics_object_as_solid_flag) {
      environment_graphics_object_internal_lines(env_bonds_box); // GL lines
   } else {
      glEnable(GL_LIGHTING);
      // glEnable(GL_LIGHT0);
      // glEnable(GL_LIGHT1);
      environment_graphics_object_internal_tubes(env_bonds_box); // GL cylinders and disks
      glDisable(GL_LIGHTING);
   }
}

// This is the GL rendering of the environment bonds box
//
void
graphics_info_t::environment_graphics_object_internal_lines(const graphical_bonds_container &env_bonds_box) const {

#if 0 // old - delete this one day.

   if (environment_show_distances == 1) {

      if (env_bonds_box.num_colours > 0) {

         graphical_bonds_lines_list<graphics_line_t> ll;
         coot::Cartesian text_pos;
         float dist;

         float dark_bg_cor = 0.0;
         if (! background_is_black_p())
            dark_bg_cor = 0.29;

         glLineStipple (1, 0x00FF);
         for (int i=0; i< env_bonds_box.num_colours; i++) {

            bool display_these_distances_flag = 1;
            if (i==0)
               if (!environment_distances_show_bumps)
                  display_these_distances_flag = 0;
            if (i==1)
               if (!environment_distances_show_h_bonds)
                  display_these_distances_flag = 0;

            if (display_these_distances_flag) {
               ll = env_bonds_box.bonds_[i]; // lightweight
               float it = float(i);
               if (it > 1.0)
                  it = 1.0;

               // now we want to draw out our bonds in various colour,
               // according to if they have a carbon or not.
               //
               glColor3f (0.8-dark_bg_cor, 0.8-0.4*it-dark_bg_cor, 0.4+0.5*it-dark_bg_cor);
               // These get turned off and set to 1 when writing stroke characters
               glEnable(GL_LINE_STIPPLE);
               glLineWidth(2.0);

               glBegin(GL_LINES);
               for (int j=0; j< env_bonds_box.bonds_[i].num_lines; j++) {
                  const coot::CartesianPair &pair = ll.pair_list[j].positions;
                  glVertex3f(pair.getStart().get_x(),
                             pair.getStart().get_y(),
                             pair.getStart().get_z());
                  glVertex3f(pair.getFinish().get_x(),
                             pair.getFinish().get_y(),
                             pair.getFinish().get_z());
               }
               glEnd();

               for (int j=0; j< env_bonds_box.bonds_[i].num_lines; j++) {
                  const coot::CartesianPair &pair = ll.pair_list[j].positions;
                  text_pos = pair.getFinish().mid_point(pair.getStart()) +
                     coot::Cartesian(0.0, 0.1, 0.1);
                  dist = (pair.getStart() - pair.getFinish()).amplitude();
                  printString(float_to_string(dist), text_pos.x(), text_pos.y(), text_pos.z());
               }
            }
         }
         glDisable(GL_LINE_STIPPLE);
      }
   }

#endif
}

// This is the GL rendering of the environment bonds box
//
void
graphics_info_t::environment_graphics_object_internal_tubes(const graphical_bonds_container &env_bonds_box) const {

   if (environment_show_distances == 1) {

      if (env_bonds_box.num_colours > 0) {

    coot::Cartesian text_pos;
    float dist;

    float dark_bg_cor = 0.0;
    if (! background_is_black_p())
       dark_bg_cor = 0.29;

    glEnable(GL_COLOR_MATERIAL);

    for (int i=0; i< env_bonds_box.num_colours; i++) {

       bool display_these_distances_flag = 1;
       if (i==0)
          if (!environment_distances_show_bumps)
     display_these_distances_flag = 0;
       if (i==1)
          if (!environment_distances_show_h_bonds)
     display_these_distances_flag = 0;

       if (display_these_distances_flag) {
          graphical_bonds_lines_list<graphics_line_t> ll = env_bonds_box.bonds_[i]; // lightweight
          float it = float(i);
          if (it > 1.0)
     it = 1.0;

          // now we want to draw out our bonds in various colour,
          // according to if they have a carbon or not.
          //
          glColor3f (0.8-dark_bg_cor, 0.8-0.4*it-dark_bg_cor, 0.4+0.5*it-dark_bg_cor);

          for (int j=0; j< env_bonds_box.bonds_[i].num_lines; j++) {
     const coot::CartesianPair &pair = ll.pair_list[j].positions;

     unsigned int n_parts = 15;
     for (unsigned int ipart=0; ipart<n_parts; ipart++) {
        if (coot::util::even_p(ipart)) {
   environment_graphics_object_internal_tube(pair, ipart, n_parts);
        }
     }

     // the distance text
     glDisable(GL_LIGHTING);
     text_pos = pair.getFinish().mid_point(pair.getStart()) +
        coot::Cartesian(0.0, 0.2, 0.2);
     // glRasterPos3f();
     dist = (pair.getStart() - pair.getFinish()).amplitude();
     printString(float_to_string(dist), text_pos.x(), text_pos.y(), text_pos.z());
     glEnable(GL_LIGHTING);
          }
       }
    }
      }
   }
}



// ----------------------------------------------------------
//
// Remember, being in GL_LINES mode will cause this to fail silently.
//

// #define HACK_OUT_GLUTBITMAPSCHARS

// static
void
graphics_info_t::printString(const std::string &s,
        const double &x, const double &y, const double &z) {

   double sf = 0.00008 * graphics_info_t::zoom;
   graphics_info_t::printString_internal(s, x, y, z, true, false, sf);
}


// static
void
graphics_info_t::printString_for_axes(const std::string &s,
         const double &x, const double &y, const double &z) {

   double sf = 0.00078;
   graphics_info_t::printString_internal(s, x, y, z, true, true, sf);
}

// static
void
graphics_info_t::printString_for_density_level(const std::string &s,
          const double &x, const double &y, const double &z) {

   double sf = 0.00028;
   sf = 0.00030;
   graphics_info_t::printString_internal(s, x, y, z, false, false, sf);
}

// static
void
graphics_info_t::printString_internal(const std::string &s,
         const double &x, const double &y, const double &z,
         bool do_unproject, bool mono_font, double sf) {

   // Don't do this. Delete this function.
}


void
graphics_info_t::environment_graphics_object_internal_tube(const coot::CartesianPair &pair,
                                                           int ipart, int n_parts) const {

   coot::Cartesian bond_vec = pair.getFinish() - pair.getStart();
   coot::Cartesian bond_frag = bond_vec * (1.0/double(n_parts));
   coot::Cartesian base_point = pair.getStart() + (bond_frag * float(ipart));
   double radius = 0.04;

   graphics_object_internal_single_tube(base_point, base_point + bond_frag,
   radius, coot::FLAT_ENDS);
}

void
graphics_info_t::graphics_object_internal_single_tube(const coot::Cartesian &base_point,
         const coot::Cartesian &end_point,
         const double &radius,
         const coot::tube_end_t &end_type) const {

#if 0 // Old OpenGL
   double top =  radius;
   double base = radius;
   int slices  = 12;
   int stacks  = 2;

   glPushMatrix();

   coot::Cartesian bond_frag = end_point - base_point;
   double height = bond_frag.amplitude();

   glTranslatef(base_point.x(), base_point.y(), base_point.z());


   // 	    This code from ccp4mg's cprimitive.cc (but modified)
   //  	    -----
   double ax;
   double rx = 0;
   double ry = 0;
   double length = height;
   double vz = bond_frag.z();

   bool rot_x = false;
   if(fabs(vz)>1e-7){
      ax = 180.0/M_PI*acos(vz/length);
      if(vz<0.0) ax = -ax;
      rx = -bond_frag.y()*vz;
      ry = bond_frag.x()*vz;
   }else{
      double vx = bond_frag.x();
      double vy = bond_frag.y();
      ax = 180.0/M_PI*acos(vx/length);
      if(vy<0) ax = -ax;
      rot_x = true;
   }

   if (rot_x) {
      glRotated(90.0, 0.0, 1.0, 0.0);
      glRotated(ax,  -1.0, 0.0, 0.0);
   } else {
      glRotated(ax, rx, ry, 0.0);
   }
   // 	    --------

   GLUquadric* quad = gluNewQuadric();
   gluCylinder(quad, base, top, height, slices, stacks);

   if (end_type == coot::FLAT_ENDS) {
      glScalef(1.0, 1.0, -1.0);
      gluDisk(quad, 0, base, slices, 2);
      glScalef(1.0, 1.0, -1.0);
      glTranslated(0,0,height);
      gluDisk(quad, 0, base, slices, 2);
   }

   if (end_type == coot::ROUND_ENDS) {
      GLUquadric* sphere_quad = gluNewQuadric();
      int sphere_slices = 10;
      int sphere_stacks = 10;
      gluSphere(sphere_quad, top, sphere_slices, sphere_stacks);
      glTranslated(0,0,height);
      gluSphere(sphere_quad, top, sphere_slices, sphere_stacks);
      gluDeleteQuadric(sphere_quad);
   }

   gluDeleteQuadric(quad);
   glPopMatrix();
#endif
}

void
graphics_info_t::graphics_object_internal_arrow(const coot::Cartesian &base_point,
   const coot::Cartesian &end_point,
   float fraction_head_size,
   const double &radius) const {

#if 0
   double top =  radius;
   double base = radius;
   int slices  = 12;
   int stacks  = 2;

   glPushMatrix();

   coot::Cartesian bond_frag = end_point - base_point;
   double height = bond_frag.amplitude();

   glTranslatef(base_point.x(), base_point.y(), base_point.z());


   // 	    This code from ccp4mg's cprimitive.cc (but modified)
   //  	    -----
   double ax;
   double rx = 0;
   double ry = 0;
   double length = height;
   double vz = bond_frag.z();

   bool rot_x = false;
   if(fabs(vz)>1e-7){
      ax = 180.0/M_PI*acos(vz/length);
      if(vz<0.0) ax = -ax;
      rx = -bond_frag.y()*vz;
      ry = bond_frag.x()*vz;
   }else{
      double vx = bond_frag.x();
      double vy = bond_frag.y();
      ax = 180.0/M_PI*acos(vx/length);
      if(vy<0) ax = -ax;
      rot_x = true;
   }

   if (rot_x) {
      glRotated(90.0, 0.0, 1.0, 0.0);
      glRotated(ax,  -1.0, 0.0, 0.0);
   } else {
      glRotated(ax, rx, ry, 0.0);
   }
   // 	    --------

   GLUquadric* quad_1 = gluNewQuadric();
   GLUquadric* quad_2 = gluNewQuadric();
   GLUquadric* quad_3 = gluNewQuadric();

   gluCylinder(quad_1, base, top, height, slices, stacks);
   glTranslated(0, 0, height);
   gluCylinder(quad_2, 2*base, 0, fraction_head_size * height, slices, stacks);

   glScalef(1.0, 1.0, -1.0);
   gluDisk(quad_3, 0, base, slices, 2);


   gluDeleteQuadric(quad_1);
   gluDeleteQuadric(quad_2);
   gluDeleteQuadric(quad_3);
   glPopMatrix();
#endif
}


void
graphics_info_t::graphics_object_internal_torus(const coot::Cartesian &base_point,
                                                const coot::Cartesian &end_point,
                                                const double &radius_1,
                                                const double &radius_2,
                                                int n_ring_atoms) const {
#if 0
   double top =  0.2;
   double base = 0.2;
   double fraction_head_size = 0.3;
   int slices  = 12;
   int stacks  = 2;

   coot::Cartesian bond_frag = end_point - base_point;
   double height = bond_frag.amplitude();

   if (height > 0.0) {
      glPushMatrix();
      glTranslatef(base_point.x(), base_point.y(), base_point.z());

      // we need to rotate the world so that the normal (bond_frag) is
      // along the z axis.  (maybe it already is - say, we are looking at
      // the feats from a raw prodrg result - which is aligned in the XY
      // plane)
      //
      coot::Cartesian normal = bond_frag * (1.0/height);
      if (normal.z() > 0.9999) {
         // std::cout << "      no rotation needed" << normal << std::endl;
      } else {

         double cos_theta_y = normal.z();
         double theta_y_rad = acos(cos_theta_y);
         double theta_z_rad = atan2(normal.y(), normal.x());
         double theta_z = clipper::Util::rad2d(theta_z_rad);
         double theta_y = clipper::Util::rad2d(theta_y_rad);
         glRotated(theta_z, 0, 0, 1); // not negative.  I don't know why.
         glRotated(theta_y, 0, 1, 0); //   ditto.
      }

      glTranslated(0, 0, 1.3 * height);
      if (n_ring_atoms == 5)
         // this makes the ring brighter.  I don't know why.
         glScalef(0.95, 0.95, 0.48);
      else
         glScalef(1.1, 1.1, 0.55);
      glutSolidTorus(radius_1, radius_2, 20, 32);
      glPopMatrix();
   }
#endif
}


void
graphics_info_t::graphics_object_internal_arc(float start_angle,
                                              float end_angle,
                                              const coot::Cartesian &start_point,
                                              const coot::Cartesian &start_dir,
                                              const coot::Cartesian &normal,
                                              float radius, float radius_inner) {

#if 0
   glPushMatrix();

   double cos_theta_y = normal.z();
   double theta_y_rad = acos(cos_theta_y);
   double theta_z_rad = atan2(normal.y(), normal.x());
   double theta_z = clipper::Util::rad2d(theta_z_rad);
   double theta_y = clipper::Util::rad2d(theta_y_rad);

   glTranslatef(start_point.x(), start_point.y(), start_point.z());

   glRotated(theta_z, 0, 0, 1); // not negative.  I don't know why.
   glRotated(theta_y, 0, 1, 0); //   ditto.

   // where is the normal after these rotations?

   // sin_t = sin(-theta_x_rad);
   // cos_t = cos(-theta_x_rad);
   // clipper::Mat33<double> x_mat(1,0,0, 0,cos_t,-sin_t, 0,sin_t,cos_t);

   double sin_y_t = sin(-theta_y_rad);
   double cos_y_t = cos(-theta_y_rad);

   double sin_z_t = sin(-theta_z_rad);
   double cos_z_t = cos(-theta_z_rad);

   clipper::Mat33<double> y_mat(cos_y_t,0,sin_y_t, 0,1,0, -sin_y_t,0,cos_y_t);
   clipper::RTop_orth y_rtop(y_mat, clipper::Coord_orth(0,0,0));

   clipper::Mat33<double> z_mat(cos_z_t,-sin_z_t,0, sin_z_t, cos_z_t,0, 0,0,1);
   clipper::RTop_orth z_rtop(z_mat, clipper::Coord_orth(0,0,0));

   clipper::Coord_orth norm_rot_0(normal.x(), normal.y(), normal.z());
   clipper::Coord_orth norm_rot_1 = norm_rot_0.transform(z_rtop);
   clipper::Coord_orth norm_rot_2 = norm_rot_1.transform(y_rtop);

   clipper::Coord_orth start_dir_rot_0(start_dir.x(), start_dir.y(), start_dir.z());
   clipper::Coord_orth start_dir_rot_1 = start_dir_rot_0.transform(z_rtop);
   clipper::Coord_orth start_dir_rot_2 = start_dir_rot_1.transform(y_rtop);

   // OK now we know the coordinates of the normal after these rotations.
   //

   double theta_z_prime_rad = atan2(start_dir_rot_2.y(), start_dir_rot_2.x());
   double theta_z_prime = clipper::Util::rad2d(theta_z_prime_rad);

   start_angle += theta_z_prime;
   end_angle   += theta_z_prime;

   // draw arc here (it's at the origin now) in the XY plane
   //

   // end_angle is less than start angle.
   double angle_step = -6.0; // degrees

   glBegin(GL_QUADS);

   for (float a=start_angle; a>=end_angle; a+=angle_step) {
      float b = a + angle_step;
      if (b < end_angle) b = end_angle;

      // these are the points on the centre of the axis of the curve
      clipper::Coord_orth pt_this_unit(cos(clipper::Util::d2rad(a)), sin(clipper::Util::d2rad(a)), 0);
      clipper::Coord_orth pt_next_unit(cos(clipper::Util::d2rad(b)), sin(clipper::Util::d2rad(b)), 0);

      // perpendicular to pt_this_unit is the direction of the circle tangent.
      // To get that, we rotate about the z-axis:
      clipper::Mat33<double> z_90_mat(0, -1, 0, 1, 0, 0, 0, 0, 1);
      clipper::RTop_orth z_90_rtop(z_mat, clipper::Coord_orth(0,0,0));

      clipper::Coord_orth pt_this(r*cos(clipper::Util::d2rad(a)), r*sin(clipper::Util::d2rad(a)), 0);
      clipper::Coord_orth pt_next(r*cos(clipper::Util::d2rad(b)), r*sin(clipper::Util::d2rad(b)), 0);

      clipper::Coord_orth this_circle_tangent(clipper::Coord_orth::cross(pt_this, norm_rot_2));
      clipper::Coord_orth next_circle_tangent(clipper::Coord_orth::cross(pt_next, norm_rot_2));

      // and a point on the ring that we shall rotate about the tangent
      clipper::Coord_orth pt_ring_this = pt_this + radius_inner * pt_this_unit;
      clipper::Coord_orth pt_ring_next = pt_next + radius_inner * pt_next_unit;

      float inner_angle_step = 6;
      for (float iangle=0; iangle<=360.1; iangle+=inner_angle_step) {

    // rotate_around_vector() args: direction position origin_shift angle
    clipper::Coord_orth pt_multi_this =
       coot::util::rotate_around_vector(this_circle_tangent,
        pt_ring_this,
        pt_this,
        clipper::Util::d2rad(iangle));
    clipper::Coord_orth pt_multi_next =
       coot::util::rotate_around_vector(next_circle_tangent,
        pt_ring_next,
        pt_next,
        clipper::Util::d2rad(iangle));

    // and the inner next of those points, i.e. neighbours
    // on the same inner ring.
    //
    clipper::Coord_orth pt_multi_this_plus =
       coot::util::rotate_around_vector(this_circle_tangent,
        pt_ring_this,
        pt_this,
        clipper::Util::d2rad(iangle+inner_angle_step));

    clipper::Coord_orth pt_multi_next_plus =
       coot::util::rotate_around_vector(next_circle_tangent,
        pt_ring_next,
        pt_next,
        clipper::Util::d2rad(iangle+inner_angle_step));

    // 4 GL_QUAD vertices
    // small radius direction
    clipper::Coord_orth smd_this(pt_multi_this - pt_this);
    clipper::Coord_orth smd_next(pt_multi_next - pt_next);
    clipper::Coord_orth smd_both((smd_this + smd_next).unit());

    // get glNormal() wrong (as this is) and things turn dark brown or just matt
    // the normal needs to be normalized!
    //
    // std::cout << "calling glNormal() with " << smd_both.format() << std::endl;
    glNormal3f(smd_both.x(), smd_both.y(), smd_both.z());

    glVertex3f(pt_multi_this.x(),      pt_multi_this.y(),      pt_multi_this.z());
    glVertex3f(pt_multi_this_plus.x(), pt_multi_this_plus.y(), pt_multi_this_plus.z());
    glVertex3f(pt_multi_next_plus.x(), pt_multi_next_plus.y(), pt_multi_next_plus.z());
    glVertex3f(pt_multi_next.x(),      pt_multi_next.y(),      pt_multi_next.z());

    // clipper::Coord_orth diff = pt_multi_this - pt_ring_this;
    // glNormal3f(diff.x(), diff.y(), diff.z());
    // glVertex3f(pt_multi_this.x(), pt_multi_this.y(), pt_multi_this.z());
      }
   }
   glEnd();

   glPopMatrix();
#endif
}

void
graphics_info_t::graphics_object_internal_dodec(const coot::old_generic_display_object_t::dodec_t &dodec) {

#if 0
   glPushMatrix();

   glTranslated(dodec.position.x(), dodec.position.y(), dodec.position.z());
   glScaled(dodec.size, dodec.size, dodec.size);

   std::vector<clipper::Coord_orth> v = dodec.d.coords();

   if (false) {
      glBegin(GL_POINTS);
      for (unsigned int i=0; i<v.size(); i++) {
         glVertex3d(v[i].x(), v[i].y(), v[i].z());
      }
      glEnd();
   }

   for (unsigned int i=0; i<12; i++) {
      glBegin(GL_TRIANGLE_FAN);
      const std::vector<unsigned int> &face = dodec.d.face(i);
      clipper::Coord_orth sum_vertex(0,0,0);
      for (unsigned int j=0; j<5; j++)
         sum_vertex += v[face[j]];
      clipper::Coord_orth face_normal(sum_vertex.unit());
      for (unsigned int j=0; j<5; j++) {
         glNormal3d(face_normal.x(), face_normal.y(), face_normal.z());
         glVertex3d(v[face[j]].x(),  v[face[j]].y(),  v[face[j]].z());
      }
      glEnd();
   }
   glPopMatrix();
#endif
}

#if 0
void
graphics_info_t::graphics_object_internal_pentakis_dodec(const coot::generic_display_object_t::pentakis_dodec_t &penta_dodec) {

   glPushMatrix();

   glTranslated(penta_dodec.position.x(), penta_dodec.position.y(), penta_dodec.position.z());
   glScaled(penta_dodec.size, penta_dodec.size, penta_dodec.size);

   std::vector<clipper::Coord_orth> v = penta_dodec.pkdd.d.coords();
   const std::vector<clipper::Coord_orth> &pv = penta_dodec.pkdd.pyrimid_vertices;

   // bool smooth_shading = true;
   bool smooth_shading = false;

   if (smooth_shading) {

      for (unsigned int i=0; i<12; i++) {

        std::vector<unsigned int> face = penta_dodec.pkdd.d.face(i);

        glBegin(GL_TRIANGLE_FAN);

        // first the base point (tip of the triangles/pyrimid)
        clipper::Coord_orth pvu(pv[i].unit());
        glNormal3d(pvu.x(), pvu.y(), pvu.z());
        glVertex3d(pv[i].x(), pv[i].y(), pv[i].z());

        for (unsigned int j=0; j<=4; j++) {
           const clipper::Coord_orth &pt = v[face[j]];
           clipper::Coord_orth ptu(pt.unit());
           glNormal3d(ptu.x(), ptu.y(), ptu.z());
           glVertex3d(pt.x(),  pt.y(),  pt.z());
        }
        const clipper::Coord_orth &pt = v[face[0]];
        clipper::Coord_orth ptu(pt.unit());
        glNormal3d(ptu.x(), ptu.y(), ptu.z());
        glVertex3d(pt.x(),  pt.y(),  pt.z());
        glEnd();
          }

   } else {

      // the surfaces of the triangles are flat/sharp and don't blend into each other.

      for (unsigned int i=0; i<12; i++) {

    std::vector<unsigned int> face = penta_dodec.pkdd.d.face(i);
    for (unsigned int j=0; j<=4; j++) {
       unsigned int idx_1 = j;
       unsigned int idx_2 = j+1;
       if (j == 4)
          idx_2 = 0;
       clipper::Coord_orth v1(v[face[idx_1]] - pv[i]);
       clipper::Coord_orth v2(v[face[idx_2]] - pv[i]);
       clipper::Coord_orth cp(clipper::Coord_orth::cross(v2, v1));
       clipper::Coord_orth cpu(cp.unit());

       glNormal3d(cpu.x(), cpu.y(), cpu.z());

       glBegin(GL_TRIANGLES);
          glVertex3d(pv[i].x(), pv[i].y(), pv[i].z());
          glVertex3d(v[face[idx_1]].x(), v[face[idx_1]].y(), v[face[idx_1]].z());
          glVertex3d(v[face[idx_2]].x(), v[face[idx_2]].y(), v[face[idx_2]].z());
       glEnd();
    }
      }
   }
   glPopMatrix();
}
#endif

// static
void
graphics_info_t::from_generic_object_remove_last_item(int object_number) {

   if (! use_graphics_interface_flag) return;

   int ngos = generic_display_objects.size();
   if (object_number >= 0) {
      if (object_number < ngos) {
         generic_display_objects[object_number].remove_last_object();
      }
   }
   graphics_draw();
}

// static
bool
graphics_info_t::is_valid_generic_display_object_number(int obj_no) {

   bool status = false;
   if (obj_no >= 0) {
      int ss = generic_display_objects.size();
      if (obj_no < ss)
         status = true;
   }
   return status;
}


void
graphics_info_t::add_label(const std::string &l, const glm::vec3 &p, const glm::vec4 &c) {
   labels.push_back(atom_label_info_t(l, p, c));
}



void
graphics_info_t::update_environment_graphics_object(int atom_index, int imol) {

   environment_object_bonds_box = molecules[imol].make_environment_bonds_box(atom_index, geom_p);

#ifndef EMSCRIPTEN
   gtk_gl_area_attach_buffers(GTK_GL_AREA(graphics_info_t::glareas[0]));
   mesh_for_environment_distances.init(environment_object_bonds_box, background_is_black_p());

   Material material;
   // mesh_for_environment_distances.mesh.setup(&shader_for_moleculestotriangles, material); 20210910-PE
   mesh_for_environment_distances.mesh.setup(material);

   labels.clear(); // remove everything

   add_distance_labels_for_environment_distances();

#endif

}

void
graphics_info_t::update_symmetry_environment_graphics_object(int atom_index, int imol) {

   symmetry_environment_object_bonds_box =
      molecules[imol].make_symmetry_environment_bonds_box(atom_index, geom_p);
}

// Just a stub.  Not used currently.
//
std::string graphics_info_t::make_mmdb_atom_string_from_go_to_atom() {

   std::string a;
   return a;

}




// ------------------------------------------------------------------
//                   density level
// ------------------------------------------------------------------
//
void
graphics_info_t::set_density_level_string(int imol, float dlevel) {

   graphics_info_t g;
   float map_sigma = g.molecules[imol].map_sigma();

   display_density_level_screen_string = "map " + int_to_string(imol);
   display_density_level_screen_string += " level = ";
   display_density_level_screen_string += float_to_string_using_dec_pl(dlevel, 3);
   // std::string units = "e/A^3";
   std::string units = molecules[imol].map_units();
   display_density_level_screen_string += units;
   display_density_level_screen_string += " (";
   display_density_level_screen_string += float_to_string(dlevel/map_sigma);
   display_density_level_screen_string += "rmsd)";

}

// ------------------------------------------------------------------
//                   geometry
// ------------------------------------------------------------------
//
float
graphics_info_t::add_measure_distance(const coot::Cartesian &p1,
                                      const coot::Cartesian &p2) {

   auto coord_orth_to_glm = [] (const clipper::Coord_orth &co) {
                               return glm::vec3(co.x(), co.y(), co.z());
                            };

   auto add_measure_distance_label = [coord_orth_to_glm] (const coot::simple_distance_object_t &sdo,
                                                          const double &dist,
                                                          const glm::vec4 &col) {
                                        clipper::Coord_orth mid_point(0.5 * (sdo.start_pos + sdo.end_pos));
                                        glm::vec3 offset_mid_point = coord_orth_to_glm(mid_point) + glm::vec3(0.15, 0.05, 0.05);
                                        unsigned int n_dec_pl = 3; // 20220503-PE was 2
                                        std::string label_str = float_to_string_using_dec_pl(static_cast<float>(dist), 3);
                                        // label_str += "ÅAÅ"; A-ring How to do this!?
                                        // degree symbol: &#176; A-ring symbol: &#197
                                        unsigned char c = 197;
                                        label_str += c;
                                        atom_label_info_t ali(label_str, offset_mid_point, col);
                                        labels_for_measure_distances_and_angles.push_back(ali);
                                     };

#ifndef EMSCRIPTEN
   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0]));
#endif

   clipper::Coord_orth cp1(p1.x(), p1.y(), p1.z());
   clipper::Coord_orth cp2(p2.x(), p2.y(), p2.z());
   double dist = sqrt((cp2-cp1).lengthsq());
   coot::simple_distance_object_t p(geometry_atom_index_1_mol_no, cp1, geometry_atom_index_2_mol_no, cp2);
   measure_distance_object_vec.push_back(p);
   Material mat;
   glm::vec4 col(0.72, 0.79, 0.72, 1.0);
#ifndef EMSCRIPTEN
   mesh_for_measure_distance_object_vec.add_dashed_line(p, mat, col);
#endif
   add_measure_distance_label(p, dist, col);
   
   graphics_draw();

   // std::cout << "INFO:: distance: " << dist << " Angstroems" << std::endl;
   logger.log(log_t::INFO, "distance:", dist, "Angstroems");
   std::string s = "Distance: ";
   s += float_to_string(dist);
   s += " A";
   add_status_bar_text(s);
   return dist;
}

// double
// graphics_info_t::display_geometry_distance_symm(int imol1, const coot::Cartesian &p1,
// 						int imol2, const coot::Cartesian &p2) {


//    coot::simple_distance_object_t p(imol1, clipper::Coord_orth(p1.x(), p1.y(), p1.z()),
// 				    imol2, clipper::Coord_orth(p2.x(), p2.y(), p2.z()));
//    distance_object_vec->push_back(p);
//    graphics_draw();
//    double d =  (p1-p2).length();
//    return d;
// }

void
graphics_info_t::add_measure_angle() const {

#ifndef EMSCRIPTEN
   gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0]));
#endif

   clipper::Coord_orth p1(angle_tor_pos_1.x(), angle_tor_pos_1.y(), angle_tor_pos_1.z());
   clipper::Coord_orth p2(angle_tor_pos_2.x(), angle_tor_pos_2.y(), angle_tor_pos_2.z());
   clipper::Coord_orth p3(angle_tor_pos_3.x(), angle_tor_pos_3.y(), angle_tor_pos_3.z());

   clipper::Coord_orth v1 = p2 - p1;
   clipper::Coord_orth v2 = p2 - p3;

   double dp = clipper::Coord_orth::dot(v1,v2);
   double len_v1 = sqrt(v1.lengthsq());
   double len_v2 = sqrt(v2.lengthsq());
   len_v1 = len_v1 < 0.0001 ? 0.0001 : len_v1;
   len_v2 = len_v2 < 0.0001 ? 0.0001 : len_v2;
   double cos_theta = dp/(len_v1 * len_v2);
   double theta = acos(cos_theta);

   auto coord_orth_to_glm = [] (const clipper::Coord_orth &co) {
                               return glm::vec3(co.x(), co.y(), co.z());
                            };

   auto add_measure_angle_label = [] (const glm::vec3 &p, double theta) {
                                     double theta_deg = (180.0/M_PI) * theta;
                                     std::string label_str = float_to_string_using_dec_pl(static_cast<float>(theta_deg), 2);
                                     // degree symbol: &#176; A-ring symbol: &#197
                                     unsigned char c = 176;
                                     label_str += c;
                                     glm::vec4 col(0.72, 0.72, 0.72, 1.0);
                                     atom_label_info_t ali(label_str, p, col); // not an atom label of course
                                     labels_for_measure_distances_and_angles.push_back(ali);
                                  };

   // p2 is the middle atom
   Material mat;
   glm::vec4 colour(0.6, 0.7, 0.5, 1.0); // 20211007-PE same as in add_dashed_line();
#ifndef EMSCRIPTEN
   mesh_for_measure_angle_object_vec.add_dashed_angle_markup(coord_orth_to_glm(p1),
                                                             coord_orth_to_glm(p2),
                                                             coord_orth_to_glm(p3), colour, mat);
#endif

   clipper::Coord_orth mid_point(0.3333 * (p1+p2+p3));
   clipper::Coord_orth centre_atom_to_mid_point_uv((mid_point-p2).unit());
   clipper::Coord_orth adjusted_mid_point(mid_point + 0.2 * centre_atom_to_mid_point_uv);
   add_measure_angle_label(coord_orth_to_glm(adjusted_mid_point), theta);

   // std::cout << "INFO:: angle: " << theta << " radians " << theta*57.29578 << " degrees " << std::endl;
   logger.log(log_t::INFO, "angle:", theta, "radians", theta*57.29578, "degrees");

   display_density_level_this_image = 1;
   display_density_level_screen_string = "  Angle:  ";
   display_density_level_screen_string += float_to_string(theta*57.29578);
   display_density_level_screen_string += " degrees";
   add_status_bar_text(display_density_level_screen_string);
   // redraw is in calling function for angles.
}


// ---- simple torsion ------
bool
graphics_info_t::set_angle_tors(int imol,
   const coot::atom_spec_t &as1,
   const coot::atom_spec_t &as2,
   const coot::atom_spec_t &as3,
   const coot::atom_spec_t &as4) {

   bool r = 0;
   if (is_valid_model_molecule(imol)) {
      mmdb::Atom *at1 = molecules[imol].get_atom(as1);
      mmdb::Atom *at2 = molecules[imol].get_atom(as2);
      mmdb::Atom *at3 = molecules[imol].get_atom(as3);
      mmdb::Atom *at4 = molecules[imol].get_atom(as4);
      if (! at1)
    std::cout << "   WARNING:: atom not found in molecule #"
      << imol << " " << as1 << std::endl;
      if (! at2)
    std::cout << "   WARNING:: atom not found in molecule #"
      << imol << " " << as2 << std::endl;
      if (! at3)
    std::cout << "   WARNING:: atom not found in molecule #"
      << imol << " " << as3 << std::endl;
      if (! at4)
    std::cout << "   WARNING:: atom not found in molecule #"
      << imol << " " << as4 << std::endl;

      if (at1 && at2 && at3 && at4) {
    angle_tor_pos_1 = coot::Cartesian(at1->x, at1->y, at1->z);
    angle_tor_pos_2 = coot::Cartesian(at2->x, at2->y, at2->z);
    angle_tor_pos_3 = coot::Cartesian(at3->x, at3->y, at3->z);
    angle_tor_pos_4 = coot::Cartesian(at4->x, at4->y, at4->z);
    r = 1;
      }
   }
   return r;
}

void
graphics_info_t::display_geometry_torsion() const {

   double torsion = get_geometry_torsion();

   display_density_level_this_image = 1;
   display_density_level_screen_string = "  Torsion:  ";
   display_density_level_screen_string += float_to_string(torsion);
   display_density_level_screen_string += " degrees";
   add_status_bar_text(display_density_level_screen_string);
   graphics_draw();

}

double
graphics_info_t::get_geometry_torsion() const {

   clipper::Coord_orth p1(angle_tor_pos_1.x(), angle_tor_pos_1.y(), angle_tor_pos_1.z());
   clipper::Coord_orth p2(angle_tor_pos_2.x(), angle_tor_pos_2.y(), angle_tor_pos_2.z());
   clipper::Coord_orth p3(angle_tor_pos_3.x(), angle_tor_pos_3.y(), angle_tor_pos_3.z());
   clipper::Coord_orth p4(angle_tor_pos_4.x(), angle_tor_pos_4.y(), angle_tor_pos_4.z());

   clipper::Coord_orth v1  = p2 - p1;
   clipper::Coord_orth v2  = p2 - p3;
   clipper::Coord_orth v3  = p3 - p4;
   clipper::Coord_orth v3r = p4 - p3;

//    std::cout << ":p1: " << p1.format() << std::endl;
//    std::cout << ":p2: " << p2.format() << std::endl;
//    std::cout << ":p3: " << p3.format() << std::endl;
//    std::cout << ":p4: " << p4.format() << std::endl;

//     std::cout << "display_geometry_angle: " << std::endl
//  	     << "      " << v1.format() << std::endl
//  	     << "      " << v2.format() << std::endl;

   // double dp = clipper::Coord_orth::dot(v1,v2);
   double len_v1 = sqrt(v1.lengthsq());
   double len_v2 = sqrt(v2.lengthsq());
   double len_v3 = sqrt(v3.lengthsq());
   len_v1 = len_v1 < 0.0001 ? 0.0001 : len_v1;
   len_v2 = len_v2 < 0.0001 ? 0.0001 : len_v2;
   len_v3 = len_v3 < 0.0001 ? 0.0001 : len_v3;
   // double cos_theta = dp/(len_v1 * len_v2);
   // double theta = acos(cos_theta);

   // we could do this when we kept the atom indcies for the torsion.
   // Now we just save the positions, so we can't give this nice
   // geometry list, ho hum...
//    std::cout << "       angle atom 1: "
// 	     << "(" << geometry_atom_index_1_mol_no << ") "
// 	     << atom1->name << "/"
// 	     << atom1->GetChainID()  << "/"
// 	     << atom1->GetSeqNum()   << "/"
// 	     << atom1->GetResName() << std::endl;
//    std::cout << "       angle atom 2: "
// 	     << "(" << geometry_atom_index_2_mol_no << ") "
// 	     << atom2->name << "/"
// 	     << atom2->GetChainID()  << "/"
// 	     << atom2->GetSeqNum()   << "/"
// 	     << atom2->GetResName() << std::endl;
//    std::cout << "       angle atom 3: "
// 	     << "(" << geometry_atom_index_3_mol_no << ") "
// 	     << atom3->name << "/"
// 	     << atom3->GetChainID()  << "/"
// 	     << atom3->GetSeqNum()   << "/"
// 	     << atom3->GetResName() << std::endl;
//    std::cout << "       angle atom 4: "
// 	     << "(" << geometry_atom_index_4_mol_no << ") "
// 	     << atom4->name << "/"
// 	     << atom4->GetChainID()  << "/"
// 	     << atom4->GetSeqNum()   << "/"
// 	     << atom4->GetResName() << std::endl;

   double tors = clipper::Coord_orth::torsion(p1, p2, p3, p4);
   double torsion = clipper::Util::rad2d(tors);
   std::cout << "       torsion: " << torsion << " degrees "
        << std::endl;

   return torsion;
}

// old interface
void
graphics_info_t::pepflip() {

   if (is_valid_model_molecule(imol_pepflip)) {
      molecules[imol_pepflip].pepflip(atom_index_pepflip);
      normal_cursor();
      model_fit_refine_unactive_togglebutton("model_refine_dialog_pepflip_togglebutton");
   }
}

void
graphics_info_t::pepflip(int imol, const coot::atom_spec_t &spec) {

   if (is_valid_model_molecule(imol)) {
      molecules[imol].pepflip(spec);
   }
}


// ----------------------------------------------------------
//                     some utilities:
// ----------------------------------------------------------

std::string
graphics_info_t::int_to_string(int i) {
   char s[100];
   for (int ii=0; ii<100; ii++) s[ii]=0;
   snprintf(s, 99, "%d", i);
   return std::string(s);
}

std::string
graphics_info_t::float_to_string(float f) {
   return coot::util::float_to_string(f);
}

std::string
graphics_info_t::float_to_string_using_dec_pl(float f, unsigned short int n_dec_pl) {
   return coot::util::float_to_string_using_dec_pl(f, n_dec_pl);
}



int
graphics_info_t::Imol_Refinement_Map() const {

   // maybe sets imol_refinement_map (if it was -1 or if the map it
   // was previously is now closed).

   int only_map = -1;  // gets set, maybe (when not a difference map)

   if (imol_refinement_map != -1) { // has been set already (or was reset)
      if (imol_refinement_map < n_molecules())
         if (imol_refinement_map >= 0)
            if (molecules[imol_refinement_map].has_xmap())
               return imol_refinement_map;
   }

   // check the molecules for maps - we can assign if there is only
   // one non-difference map.
   //
   std::vector<int> direct_maps;
   for (int imol=0; imol<n_molecules(); imol++) {
      if (molecules[imol].has_xmap()) {
         if (! molecules[imol].is_difference_map_p()) {
            direct_maps.push_back(imol);
         }
      }
   }

   if (direct_maps.size() == 1) {
      // let's set it then
      imol_refinement_map = direct_maps[0];
   } else {
      imol_refinement_map = -1;
   }

   return imol_refinement_map;
}

int
graphics_info_t::set_imol_refinement_map(int imol) {

   int r = -1;
   if (molecules[imol].has_xmap()) {
      imol_refinement_map = imol;
      r = imol;
   }
   return r;
}

// for threading, static
void
graphics_info_t::update_maps_for_mols(const std::vector<int> &mol_idxs) {
   for (unsigned int i=0; i<mol_idxs.size(); i++)
      graphics_info_t::molecules[mol_idxs[i]].update_map(auto_recontour_map_flag);
}

#include <thread>


void
graphics_info_t::update_maps() {

   if (GetActiveMapDrag() == 1) {

      // now map updates are internally threaded - we don't need
      // this extra threaded mechanism (makes things crash?)

      bool do_threaded_map_updates = false;

      if (! do_threaded_map_updates) {
         for (int ii=0; ii<n_molecules(); ii++) {
            if (molecules[ii].has_xmap()) {
               molecules[ii].update_map(auto_recontour_map_flag); // to take account
                                                           // of new rotation centre.
            }
         }

      } else {

         // unsigned int n_threads = 4;
         unsigned int n_threads = coot::get_max_number_of_threads();
         // std::cout << "got n_threads: " << n_threads << std::endl;

         if (n_threads == 0) {
            for (int ii=0; ii<n_molecules(); ii++) {
               if (molecules[ii].has_xmap()) {
                  molecules[ii].update_map(graphics_info_t::auto_recontour_map_flag); // to take account
                                                                              // of new rotation centre.
               }
            }
         } else {
            std::vector<std::thread> threads;
            std::vector<int> molecules_with_maps;
            for (int ii=0; ii<n_molecules(); ii++) {
               if (molecules[ii].has_xmap()) { // or nxmap
                  molecules_with_maps.push_back(ii);
               }
            }

            // we must make sure that the threads don't update the same map
            //

            std::vector<std::vector<int> > maps_vec_vec(n_threads);
            unsigned int thread_idx = 0;
            // put the maps in maps_vec_vec
            for (unsigned int ii=0; ii<molecules_with_maps.size(); ii++) {
               maps_vec_vec[thread_idx].push_back(molecules_with_maps[ii]);
               thread_idx++;
               if (thread_idx == n_threads) thread_idx = 0;
            }

            for (unsigned int i_thread=0; i_thread<n_threads; i_thread++) {
               const std::vector<int> &mv = maps_vec_vec[i_thread];
               threads.push_back(std::thread(update_maps_for_mols, mv));
            }
            for (unsigned int i_thread=0; i_thread<n_threads; i_thread++)
               threads.at(i_thread).join();

         }
      }
   } // active map drag test
}

// simple
// static
void
graphics_info_t::add_vector_to_rotation_centre(const coot::Cartesian &vec) {

   rotation_centre_x += vec.x();
   rotation_centre_y += vec.y();
   rotation_centre_z += vec.z();
}

void
graphics_info_t::add_vector_to_RotationCentre(const coot::Cartesian &vec) {

   rotation_centre_x += vec.x();
   rotation_centre_y += vec.y();
   rotation_centre_z += vec.z();

   update_maps();
   for (int ii=0; ii<n_molecules(); ii++) {
      molecules[ii].update_symmetry();
   }
   graphics_draw();
}



// return NULL on not found:
mmdb::Atom *
graphics_info_t::find_atom_in_moving_atoms(const coot::atom_spec_t &at) const {

   mmdb::Atom *cat = NULL;
   if (moving_atoms_asc->mol != NULL) {

      int SelHnd = coot::get_selection_handle(moving_atoms_asc->mol, at);
      int nSelAtoms;
      mmdb::PPAtom local_SelAtom = NULL;
      moving_atoms_asc->mol->GetSelIndex(SelHnd, local_SelAtom, nSelAtoms);
      if (nSelAtoms > 0)
    cat = local_SelAtom[0];
       std::cout << "DEBUG:: in find_atom_in_moving_atoms: here are the "
 		<< nSelAtoms << " qualifying atoms..." << std::endl;
       for(int i=0; i<nSelAtoms; i++)
 	 std::cout << "      " << i << "  " << local_SelAtom[i] << std::endl;
      moving_atoms_asc->mol->DeleteSelection(SelHnd);
   } else {
     std::cout << "WARNING:: OOps: moving_atoms_asc->mol is NULL" << std::endl;
   }
   return cat;
}



int
graphics_info_t::load_db_main() {
   return 1;
}



//
int
graphics_info_t::create_pointer_atom_molecule_maybe() const {

   int i = -1; // must be changed by this function.

   if (user_pointer_atom_molecule >= 0) {
      if (user_pointer_atom_molecule < n_molecules()) {
    if (molecules[user_pointer_atom_molecule].open_molecule_p()) {
       i = user_pointer_atom_molecule;
    }
      }
   }

   if (i == -1) {
      // user did not explicictly set the molecule, the usual case:

      for (int imol=0; imol<n_molecules(); imol++) {
    if (molecules[imol].open_molecule_p()) { // not closed
       if (molecules[imol].name_ == "Pointer Atoms") {
          return imol;
       }
    }
      }

      // If we get here, it was not found, let's create one:

      std::cout << "Creating a molecule for Pointer Atoms" << std::endl;
      mmdb::Manager *MMDBManager = new mmdb::Manager();

      // do we attach a model, chain etc here?
      // Yes.
      mmdb::Model *model_p = new mmdb::Model;
      mmdb::Chain *chain_p = new mmdb::Chain;

      model_p->AddChain(chain_p);
      MMDBManager->AddModel(model_p);

      atom_selection_container_t asc = make_asc(MMDBManager);
      int imol = create_molecule();
      graphics_info_t g;
      molecules[imol].install_model(imol, asc, g.Geom_p(), "Pointer Atoms", 1);
      return imol;
   }
   return i;
}

void
graphics_info_t::update_things_on_move_and_redraw() {

   update_things_on_move();
   graphics_draw();
}


void
graphics_info_t::update_things_on_move() {

   for (int ii=0; ii<n_molecules(); ii++) {
      if (GetActiveMapDrag())
         molecules[ii].update_map(auto_recontour_map_flag);
      molecules[ii].update_clipper_skeleton();
      molecules[ii].update_symmetry();
   }
   make_pointer_distance_objects();
   setup_graphics_ligand_view_using_active_atom();
}

// return the state whether to really show the baton.
bool
graphics_info_t::start_baton_here() {

   baton_root = RotationCentre();

   int imol_for_skel = imol_for_skeleton(); // if unset, sets and
   // returns if only one
   // map, else return -1
   if (imol_for_skel < 0) {

      std::cout << "WARNING: no skeleton found " << std::endl;

      std::vector<int> map_molecules = valid_map_molecules();

      if (map_molecules.size() > 0) {
         GtkWidget *w = wrapped_create_skeleton_dialog(1);
         gtk_widget_set_visible(w, TRUE);
         return 0;

      } else {

         // 20091218 It is as it was - No map.
         //
         // GtkWidget *w = create_baton_mode_make_skeleton_dialog();

         GtkWidget *w = widget_from_builder("baton_mode_make_skeleton_dialog");
         g_object_set_data(G_OBJECT(w), "imol", GINT_TO_POINTER(imol_for_skel));
         gtk_widget_set_visible(w, TRUE);

         return 0;
      }

   } else {

      molecules[imol_for_skel].fill_skeleton_treenodemap(); // filled only if not filled.
      clipper::Coord_grid dummy_cg;
      short int use_dummy_cg = 0;
      baton_next_directions(imol_for_skel, NULL, baton_root, dummy_cg, use_dummy_cg);
      baton_next_ca_options_index = 0;
      baton_tip = baton_tip_by_ca_option(baton_next_ca_options_index);
      return 1;
   }
}

// Fill baton_next_ca_options, used in accept_baton_position() and
// several other functions.
//
// Generate the "Atom Guide Points" molecule (based on
// baton_next_ca_options) and display it.
//
void
graphics_info_t::baton_next_directions(int imol_for_skel, mmdb::Atom *latest_atom,
          const coot::Cartesian &baton_root_in,
          const clipper::Coord_grid &cg_start,
          short int use_cg_start) {


//    std::cout << "DEBUG in baton_next_directions imol_for_skel is "
// 	     << imol_for_skel << std::endl;

   std::vector<clipper::Coord_orth> previous_ca_positions;
   // store the position of the just accepted atom as a previous atom
   //
   // previous_ca_positions.push_back(to_coord_orth(baton_root_in));
   int imol_baton_atoms = baton_build_atoms_molecule();


   // std::cout << "DEBUG INFO:::::: latest_atom is " << latest_atom << std::endl;

   if (latest_atom == NULL) {
      previous_ca_positions.push_back(to_coord_orth(baton_root_in));
   } else {
      previous_ca_positions = molecules[imol_baton_atoms].previous_baton_atom(latest_atom,
         baton_build_direction_flag);
   }
//    std::cout << "DEBUG: in graphics: cg_start: " << cg_start.format() << "  "
//          << use_cg_start << std::endl;
   baton_next_ca_options = molecules[imol_for_skel].next_ca_by_skel(previous_ca_positions,
                                                                    cg_start,
                                                                    use_cg_start,
                                                                    3.8,
                                                                    skeleton_level,
                                                                    max_skeleton_search_depth);

   // Print out the baton_next_ca_options
   //
   std::cout << "-- baton_next_ca_options" << std::endl;
   for(unsigned int i=0; i<baton_next_ca_options.size(); i++) {
      std::cout << "   " << baton_next_ca_options[i].score  << "  "
                << baton_next_ca_options[i].position.format() << std::endl;
   }
   std::cout << "--" << std::endl;

   // Graphics the baton_next_ca_options
   //
   std::string molname("Baton Atom Guide Points");
   if (baton_tmp_atoms_to_new_molecule) {
      create_molecule_and_display(baton_next_ca_options, molname);
   } else {
      update_molecule_to(baton_next_ca_options, molname);
   }
}



void
graphics_info_t::draw_baton_object() {

   if (graphics_info_t::draw_baton_flag) {

      if (true)
         std::cout << "baton from " << baton_root << " to " << baton_tip
                   << " draw_baton_flag: " << draw_baton_flag << std::endl;

#if 0
      // needs replacing

      coot::Cartesian centre = unproject_xyz(0, 0, 0.5);
      coot::Cartesian front  = unproject_xyz(0, 0, 0.0);
      // coot::Cartesian right  = unproject_xyz(1, 0, 0.5);
      // coot::Cartesian screen_x = (right - centre);
      coot::Cartesian screen_z = (front - centre);
      screen_z.unit_vector_yourself();

      coot::Cartesian baton_vec = baton_tip - baton_root;
      coot::Cartesian buv = baton_vec;
      baton_vec.unit_vector_yourself(); // bleugh

      coot::Cartesian arrow_head_pt = baton_tip - baton_vec * 0.5;
      // rotate arrow_head_pt_1 30 degrees about screen z
      coot::Cartesian rp_1 = arrow_head_pt.rotate_about_vector(screen_z, baton_tip,  30.0*M_PI/180.0);
      coot::Cartesian rp_2 = arrow_head_pt.rotate_about_vector(screen_z, baton_tip, -30.0*M_PI/180.0);

      glLineWidth(5);
      glColor3f (0.8, 0.8, 0.9);
      glBegin(GL_LINES);
      glVertex3f(baton_root.x(), baton_root.y(), baton_root.z());
      glVertex3f(baton_tip.x(),  baton_tip.y(),  baton_tip.z());
      glVertex3f(baton_tip.x(),  baton_tip.y(),  baton_tip.z());
      glVertex3f(rp_1.x(), rp_1.y(), rp_1.z());
      glVertex3f(baton_tip.x(), baton_tip.y(), baton_tip.z());
      glVertex3f(rp_2.x(), rp_2.y(), rp_2.z());
      glEnd();
#endif

   }

}

// aka imol_for_baton_atoms()
int
graphics_info_t::baton_build_atoms_molecule() const {

   int imol = -1;

   /////////////////////////////////////////////////// FIXME ////////////////////////
   for (int im=0; im<n_molecules(); im++) {
      if (molecules[im].name_ == "Baton Atoms") {
    return im;
      }
   }

   if (0) {
      std::cout << "------------- existing molecule names: " << std::endl;
      for (int im=0; im<n_molecules(); im++)
    std::cout << im << " " << molecules[im].name_ << std::endl;
   }


   // std::cout << "INFO:: Creating a molecule for Baton Atoms" << std::endl;
   logger.log(log_t::INFO, "Creating a molecule for Baton Atoms");
   // not found, let's create one:
   mmdb::Manager *MMDBManager = new mmdb::Manager();

   // do we attach a model, chain etc here?
   // Yes.
   mmdb::Model *model_p = new mmdb::Model;
   mmdb::Chain *chain_p = new mmdb::Chain;
   chain_p->SetChainID(baton_build_chain_id.c_str());

   model_p->AddChain(chain_p);
   MMDBManager->AddModel(model_p);

   // now lets add an mmdbcryst to the mmdbmanger, which is generated
   // from the skeleton map.
   //
   int  imol_for_skel = imol_for_skeleton();
   if (imol_for_skel >= 0) {
      // mmdb::CMMDBCryst *cryst = new mmdb::CMMDBCryst;
      MMDBManager->SetCell(molecules[imol_for_skel].xskel_cowtan.cell().descr().a(),
      molecules[imol_for_skel].xskel_cowtan.cell().descr().b(),
      molecules[imol_for_skel].xskel_cowtan.cell().descr().c(),
      clipper::Util::rad2d(molecules[imol_for_skel].xskel_cowtan.cell().descr().alpha()),
      clipper::Util::rad2d(molecules[imol_for_skel].xskel_cowtan.cell().descr().beta()),
      clipper::Util::rad2d(molecules[imol_for_skel].xskel_cowtan.cell().descr().gamma()), 1);

      std::string spacegroup = molecules[imol_for_skel].xskel_cowtan.spacegroup().symbol_hm();

      // std::cout << "INFO:: setting spacegroup of Baton Atoms to be: " << spacegroup << std::endl;
      logger.log(log_t::INFO, "setting spacegroup of Baton Atoms to be:", spacegroup);
      // std::cout << "INFO:: setting cell of Baton Atoms to be: "
      //    << molecules[imol_for_skel].xskel_cowtan.cell().format() << std::endl;
      logger.log(log_t::INFO, "setting cell of Baton Atoms to be:", molecules[imol_for_skel].xskel_cowtan.cell().format());

      int istat_spgr = MMDBManager->SetSpaceGroup(spacegroup.c_str());
      if (istat_spgr != 0) {
    std::cout << "WARNING:: Problem:: mmdb does not understand space group: " << spacegroup << std::endl;
      }

   } else {
      std::cout << "WARNING: skeleton not found - no symmetry for Baton Atoms " << std::endl;
   }

   atom_selection_container_t asc = make_asc(MMDBManager);
   asc.SelectionHandle = -1;
   imol = create_molecule();
   molecules[imol].install_model(imol, asc, graphics_info_t::Geom_p(), "Baton Atoms", 1);

   // std::cout << "INFO:: returning baton atom molecule " << imol << std::endl;
   logger.log(log_t::INFO, "returning baton atom molecule", imol);
   return imol;
}

void
graphics_info_t::accept_baton_position() {

   int imol_for_skel = imol_for_skeleton();

   // First add the atom to the baton build molecule
   //
   //
   mmdb::Atom *baton_atom = NULL; // for baton_next_directions() usage?
   int imol = baton_build_atoms_molecule();

   std::cout << "--------------------- in accept_baton_position() imol is " << imol << std::endl;
   if (imol >= 0) {
      baton_atom = molecules[imol].add_baton_atom(baton_tip,
                                                  baton_build_start_resno,
                                                  baton_build_chain_id,
                                                  baton_build_params_active,
                                                  baton_build_direction_flag);

      if (baton_atom == 0) {
         // we didn't add one (because there were no chains?)
         mmdb::Model *model_p = molecules[imol].atom_sel.mol->GetModel(1);
         if (! model_p) {
            std::cout << "in accept_baton_position fallback: no model " << std::endl;
         } else {
            mmdb::Chain *chain_p = new mmdb::Chain;
            chain_p->SetChainID("A");
            model_p->AddChain(chain_p);
            baton_atom = molecules[imol].add_baton_atom(baton_tip,
                                                        baton_build_start_resno,
                                                        baton_build_chain_id,
                                                        baton_build_params_active,
                                                        baton_build_direction_flag);
         }
      }
      baton_build_params_active = 0; // This flag was set after
      // set_baton_build_params.  We
      // clear it now so that we don't
      // any more force the start resno -
      // molecule_class_info_t::add_baton_atom
      // can work it out.
   }
   std::cout << "setting screen rotation centre to " << baton_tip << std::endl;
   setRotationCentre(baton_tip);
   for(int ii=0; ii<n_molecules(); ii++) {
      // but not skeleton, lets do skeleton only on a middle-mouse recentre
      molecules[ii].update_map(true);
      molecules[ii].update_symmetry();
   }

   // debug
   //
   // std::cout << "DEBUG:: in accept_baton_position baton_atom is "
   // << baton_atom << std::endl;

   // But first show us the options for the next point as dummy atoms
   //

   // int imol_for_skel = imol_for_skeleton();
   // internally set the baton_next_ca_options
   //
   if (imol_for_skel < 0) {
      std::cout << "Ooops:: must have a skeleton first" << std::endl;
   } else {
      short int use_cg = 1;
      std::cout << "DEBUG:: accept_baton_position: " << baton_next_ca_options.size() << " "
                << baton_next_ca_options_index << std::endl;
      if (baton_next_ca_options.size() > 0) {
         clipper::Coord_grid cg = baton_next_ca_options[baton_next_ca_options_index].near_grid_pos;
         baton_next_directions(imol_for_skel, baton_atom, baton_tip, cg, use_cg); // old tip
      } else {
         clipper::Coord_grid cg;
         use_cg = 0;
         baton_next_directions(imol_for_skel, baton_atom, baton_tip, cg, use_cg);
      }
   }

   // Now set the baton tip to the next point:

   // set these for redraw
   baton_root = baton_tip;
   baton_next_ca_options_index = 0; // for next
   baton_length = 3.8; // reset to to optimal

   baton_tip = baton_tip_by_ca_option(baton_next_ca_options_index);

   graphics_draw();
}



coot::Cartesian
graphics_info_t::baton_tip_by_ca_option(int index) const {

   coot::Cartesian tip_pos(0.0, 0.0, 0.0);
   unsigned int uindex = index;

   {
      if (uindex >= baton_next_ca_options.size()) {
         if ((uindex == 0) && (baton_next_ca_options.size() == 0)) {
            // std::cout << "INFO:: no baton next positions from here\n";
            logger.log(log_t::INFO, "no baton next positions from here");
            tip_pos = non_skeleton_tip_pos();
         } else {
            std::cout << "ERROR: bad baton_next_ca_options index: "
                      << index << " size " << baton_next_ca_options.size()
                      << std::endl;
         }
      } else {
         // now we want a vector baton_length in the direction starting
         // at baton_root to baton_next_ca_options[index]
         //
         coot::Cartesian target_point = to_cartesian(baton_next_ca_options[index].position);
         std::cout << "Ca option " << index << " score: "
                   << baton_next_ca_options[index].score << std::endl;
         coot::Cartesian target_dir = target_point - baton_root;
         target_dir.unit_vector_yourself();
         target_dir *= baton_length;
         tip_pos = target_dir + baton_root;
      }
   }
   return tip_pos;
}

coot::Cartesian
graphics_info_t::non_skeleton_tip_pos() const {

   double l=1.56;
   coot::Cartesian new_tip_pos = baton_root + coot::Cartesian(l,l,l);
   // consider using baton_previous_ca_positions, if/when it get set properly.
   return new_tip_pos;
}



// Gawd, Kevin's having an effect on me.. (my programming at least...)
// "const double &" indeed?  Tsk! Whatever next?
//
void
graphics_info_t::rotate_baton(const double &x, const double &y) {

   mouse_current_x = x;
   mouse_current_y = y;
   double diff;

   diff  = mouse_current_x - GetMouseBeginX();
   diff += mouse_current_y - GetMouseBeginY();

#if 0

   // needs replacing

   coot::Cartesian centre = unproject_xyz(0, 0, 0.5);
   coot::Cartesian front  = unproject_xyz(0, 0, 0.0);
   coot::Cartesian screen_z = (front - centre);

   clipper::Coord_orth new_pos = coot::util::rotate_around_vector(to_coord_orth(screen_z),
     to_coord_orth(baton_tip),
     to_coord_orth(baton_root),
     0.01*diff);

   baton_tip = to_cartesian(new_pos);
   graphics_draw();
#endif

}

void
graphics_info_t::toggle_baton_mode() {

   if (baton_mode == 0) {
      baton_mode = 1;
      // std::cout << "INFO::baton rotation mode on." << std::endl;
      logger.log(log_t::INFO, "baton rotation mode on.");
   } else {
      baton_mode = 0;
      // std::cout << "INFO::baton rotation mode off." << std::endl;
      logger.log(log_t::INFO, "baton rotation mode off.");
   }
}

void
graphics_info_t::baton_tip_try_another() {

   baton_next_ca_options_index++;

   // make baton_next_ca_options_index an unsigned int
   if (baton_next_ca_options_index >= int(baton_next_ca_options.size())) {
      std::cout << "info: cycling back to start of ca options" << std::endl;
      baton_next_ca_options_index = 0;
   }
   baton_tip = baton_tip_by_ca_option(baton_next_ca_options_index);
   graphics_draw();
}



void
graphics_info_t::baton_tip_previous() {

   if (baton_next_ca_options_index == 0) {
      baton_next_ca_options_index = int(baton_next_ca_options.size()-1);
   } else {
      baton_next_ca_options_index--;
   }
   baton_tip = baton_tip_by_ca_option(baton_next_ca_options_index);
   graphics_draw();
}



void
graphics_info_t::shorten_baton() {

   double short_factor = 0.952;
   baton_length *= short_factor;
//    baton_tip = baton_tip_by_ca_option(baton_next_ca_options_index);
   coot::Cartesian baton_vec = baton_tip - baton_root;
   baton_vec *= short_factor;
   baton_tip = baton_root + baton_vec;
   graphics_draw();
}


void
graphics_info_t::lengthen_baton() {

   double lengthen_factor = 1.05;
   baton_length *= lengthen_factor;
   coot::Cartesian baton_vec = baton_tip - baton_root;
   baton_vec *= lengthen_factor;
   baton_tip = baton_root + baton_vec;
   graphics_draw();
}

void
graphics_info_t::baton_build_delete_last_residue() {

   int imol = baton_build_atoms_molecule();
   if (is_valid_model_molecule(imol)) {
      std::pair<short int, mmdb::Atom *> new_centre
    = molecules[imol].baton_build_delete_last_residue();
      if (new_centre.first) {
    coot::Cartesian new_centre_cart(new_centre.second->x,
    new_centre.second->y,
    new_centre.second->z);
    setRotationCentre(new_centre_cart);
    baton_root = new_centre_cart;
    int imol_for_skel = imol_for_skeleton();
    if (imol_for_skel >= 0) {
       short int use_cg = 1;
       // recall std::vector<coot::scored_skel_coord> *baton_next_ca_options;
       //
       std::pair <short int,clipper::Coord_grid> cg =
          molecules[imol_for_skel].search_for_skeleton_near(new_centre_cart);

       if (cg.first)
          baton_next_directions(imol_for_skel, new_centre.second, new_centre_cart, cg.second, use_cg);
    }
    baton_next_ca_options_index = 0;
    baton_length = 3.8;
    baton_tip = baton_tip_by_ca_option(baton_next_ca_options_index);
      }
   }
   graphics_draw();
}


// should be similar to Imol for refinement?
int
graphics_info_t::imol_for_skeleton() const {

//    for (int imol=0; imol<n_molecules; imol++) {
//       if (molecules[imol].xskel_is_filled) {
// 	 return imol;
//       }
//    }

//    int nmaps = 0;
//    int only_map;

//     if (map_for_skeletonize == -1) {
//        for (int imol=0; imol<n_molecules; imol++) {
//  	 if (molecules[imol].has_map()) {
//  	    nmaps++;
//  	 only_map = imol;
//  	 }
//        }
//        if (nmaps == 1) {
//  	 map_for_skeletonize = only_map;
//        }
//     }

   return map_for_skeletonize;
}

;
void
graphics_info_t::create_molecule_and_display(std::vector<coot::scored_skel_coord> &pos_position,
        const std::string &molname) {

   int imol = create_empty_molecule(molname);
   std::vector<coot::Cartesian> cv;
   // now add atoms:
   for (unsigned int i=0; i<pos_position.size(); i++) {
      coot::Cartesian c(pos_position[i].position.x(),
   pos_position[i].position.y(),
   pos_position[i].position.z());
      cv.push_back(c);
   }
   molecules[imol].add_multiple_dummies(cv);

}

// as above, except we update molecule with name molname to the
// pos_positions (and delete everything else).
//
void
graphics_info_t::update_molecule_to(std::vector<coot::scored_skel_coord> &pos_position,
     const std::string &molname) {

   int imol = lookup_molecule_name(molname);

   if (pos_position.size() > 0) {
      if (is_valid_model_molecule(imol)) {
    graphics_info_t::molecules[imol].update_molecule_to(pos_position);
      } else {
    create_molecule_and_display(pos_position, molname);
      }
   } else {
      std::cout << "WARNING:: No atoms guide points in update_molecule_to."
   << "  Not updating guide points molecule" << std::endl;
   }
}

// return -1 on no such map.
int
graphics_info_t::lookup_molecule_name(const std::string &molname) const {

   for (int imol=0; imol<n_molecules(); imol++) {
      if (is_valid_map_molecule(imol) || (is_valid_model_molecule(imol))) {
    if (0)
       std::cout << "comparing map names:\n     :"
         << graphics_info_t::molecules[imol].name_ << ":\n   "
         << "  :" << molname << ":" << std::endl;
    if (graphics_info_t::molecules[imol].name_ == molname) {
       return imol;
    }
      }
   }
   return -1;
}



int
graphics_info_t::create_empty_molecule(const std::string &molname) {

   std::cout << "Creating a molecule for " << molname << std::endl;

   mmdb::Manager *MMDBManager = new mmdb::Manager();

   mmdb::Model *model_p = new mmdb::Model;
   mmdb::Chain *chain_p = new mmdb::Chain;

   model_p->AddChain(chain_p);
   MMDBManager->AddModel(model_p);

   atom_selection_container_t asc = make_asc(MMDBManager);
   int imol = create_molecule();
   molecules[imol].install_model(imol, asc, graphics_info_t::Geom_p(), molname, 1);
   asc.read_error_message = "No error";
   asc.read_success = 1;
   return imol;
}


// ------------------------------------------------------------------
//                        undo functions
// ------------------------------------------------------------------

// This is the callback when "Undo" is pressed:
//
// There is a problem now that we are including redo code: the
// question is what to do when we reach the end of the modifications
// in the current undo_molecule.  This is hard... let's go with the
// current functionality for the moment, which is to unset it when we
// get to the end of the modifications list.
//
// But we need to modify the response of Undo_molecule() depending on
// whether it was asked by apply_undo or apply_redo, we do that by
// passing an enumerated type.
//
//
int
graphics_info_t::apply_undo() {

   int state = 0;
   // use (class) enum for the return value?
   int umol = Undo_molecule(coot::UNDO); // return -2 on uncertainty. Return -1 on "unset"/nothing
   if (umol == -2) {
      if (use_graphics_interface_flag) {

         // GtkWidget *dialog = create_undo_molecule_chooser_dialog();
         GtkWidget *dialog = widget_from_builder("undo_molecule_chooser_dialog");
         GtkWidget *combobox = widget_from_builder("undo_molecule_chooser_comboboxtext");
         if (false) {
            std::cout << "DEBUG:: apply_undo(): dialog: " << dialog << std::endl;
            std::cout << "DEBUG:: apply_undo(): combobox: " << combobox << std::endl;
         }
         fill_combobox_with_undo_options(combobox);
         gtk_widget_set_visible(dialog, TRUE);
      }
   } else {
      if (umol == -1) {
         std::string mess = "There are no molecules with modifications that can be undone";
         logger.log(log_t::INFO, mess);
      } else {

         std::string cwd = coot::util::current_working_dir();
         if (molecules[umol].Have_modifications_p()) {
            if (molecules[umol].is_displayed_p()) {
               state = molecules[umol].apply_undo(cwd);
               if (use_graphics_interface_flag) {
                  graphics_draw();

                  // need to update the atom and residue list in Go To Atom widget
                  // (maybe)
                  update_go_to_atom_window_on_changed_mol(umol);

                  // update the ramachandran, if there was one
                  rama_plot_boxes_handle_molecule_update(umol);
                  draw_rama_plots();

                  // now update the geometry graphs, so get the asc
                  atom_selection_container_t u_asc = molecules[umol].atom_sel;

                  update_validation(umol);

                  run_post_manipulation_hook(umol, 0);
               }
            } else {
               if (use_graphics_interface_flag) {
                  std::string s = "WARNING:: Coot will not undo modifications on a \n";
                  s += "molecule that is not displayed";
                  info_dialog(s);
               }
            }
         } else {
            undo_molecule = -1; // reset it
            if (use_graphics_interface_flag) {
               std::cout << "WARNING:: !!!  Changing the molecule to which "
                         << "\"Undo\"s are done." << std::endl;
               std::string s = "WARNING:: Changing to Undo molecule";
               add_status_bar_text(s);
            }
            apply_undo();       // find another molecule to undo
         }
      }
   }

   // and now tinker with the Redo button to make it active
   //
   activate_redo_button();  // has protection for --no-graphics
   return state;
}

int
graphics_info_t::apply_redo() {

   int state = 0;

   int umol = Undo_molecule(coot::REDO);
   if (umol == -2) { // ambiguity
      // GtkWidget *dialog = create_undo_molecule_chooser_dialog();
      GtkWidget *dialog = widget_from_builder("undo_molecule_chooser_dialog");
      GtkWidget *combobox = widget_from_builder("undo_molecule_chooser_combobox");
      fill_combobox_with_undo_options(combobox);
      gtk_widget_set_visible(dialog, TRUE);
   } else {
      if (umol == -1) { // unset
         std::cout << "There are no molecules with modifications "
                   << "that can be re-done" << std::endl;
      } else {

         if (molecules[umol].Have_redoable_modifications_p()) {
            // std::cout << "DEBUG:: applying redo" << std::endl;
            std::string cwd = coot::util::current_working_dir();
            state = molecules[umol].apply_redo(cwd);
            graphics_draw();

            // need to update the atom and residue list in Go To Atom widget
            // (maybe)
            update_go_to_atom_window_on_changed_mol(umol);
            // BL says:: from undo, maybe more should be updated!?!
            // update the ramachandran, if there was one

            // update the ramachandran, if there was one
            rama_plot_boxes_handle_molecule_update(umol);
            draw_rama_plots();

            // now update the geometry graphs, so get the asc
            atom_selection_container_t u_asc = molecules[umol].atom_sel;

            update_validation(umol);

            run_post_manipulation_hook(umol, 0);

         } else {
            // std::cout << "DEBUG:: not applying redo" << std::endl;
         }
      }
   }
   return state;
}

// This is a noddy - better is needed.
//
// No account is taken to deactivate the button again when we have run
// out of redos.
//
// When the widget gets created, the redo button is always insensitve
// - it should depend on if there are redoable molecules.
//
void
graphics_info_t::activate_redo_button() {

#if 0
   GtkWidget *dialog = model_fit_refine_dialog;
   if (dialog) {
      // which it should be!
      GtkWidget *button = widget_from_builder("model_refine_dialog_redo_button");
      gtk_widget_set_sensitive(button, TRUE);
   }
#endif

}



// Return -2 on ambiguity, -1 on unset and a molecule number >=0 for
// no ambiguity (or undo_molecule has been set already).
//
int
graphics_info_t::Undo_molecule(coot::undo_type undo_type) const {

   int r = -1;
   if (is_valid_model_molecule(undo_molecule)) {
      r = undo_molecule;
   } else {
      int n_mol = 0;
      for (int imol=0; imol<n_molecules(); imol++) {

         if (molecules[imol].open_molecule_p()) {

            if (undo_type == coot::UNDO) {
               if (molecules[imol].Have_modifications_p()) {
                  n_mol++;
                  r = imol;
               }
            }

            if (undo_type == coot::REDO) {
               if (molecules[imol].Have_redoable_modifications_p()) {
                  n_mol++;
                  r = imol;
               }
            }
         }
      }
      if (n_mol > 1) {
         r = -2;
      }
   }
   return r;
}


void
graphics_info_t::set_bond_thickness(int imol, float t) {

   auto close_float_p = [] (float f1, float f2) {
                           return (fabsf(f1-f2) < 0.001);
                        };

   std::cout << "debug:: graphics_info_t::set_bond_thickness() called with imol " << imol << " thickness " << t << std::endl;
   if (is_valid_model_molecule(imol)) {
      if (molecules[imol].has_model()) {
         if (! close_float_p(molecules[imol].get_bond_thickness(), t)) {
            molecules[imol].set_bond_thickness(t);
            molecules[imol].make_bonds_type_checked(__FUNCTION__);
            graphics_draw();
         }
      }
   }
}

void
graphics_info_t::set_bond_thickness_intermediate_atoms(float t) {

   bond_thickness_intermediate_atoms = t;

}

void
graphics_info_t::crosshairs_text() const {

   if (draw_crosshairs_flag > 0) {
      std::cout << "Crosshair ticks: 1.54A (C-C bond), 2.7A (H-bond), 3.8A (Ca-Ca)\n";
   }
}

// Thank you for this idea Stuart Makay
//
void
graphics_info_t::pick_cursor_maybe() {

   if (control_key_for_rotate_flag) {
      pick_cursor_real();
   }
}

void
graphics_info_t::pick_cursor_real() {

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
      // 20220528-PE FIXME cursor
#else
   if (use_graphics_interface_flag) {
      //    GdkCursorType c = GDK_CROSSHAIR;
      GdkCursorType c = pick_cursor_index;
      GdkCursor *cursor;
      // cursor = gdk_cursor_new (c);
      cursor = gdk_cursor_new_for_display (gdk_display_get_default(), c);
      GdkWindow *window = gtk_widget_get_window(glareas[0]);
      gdk_window_set_cursor(window, cursor);
      // gdk_cursor_destroy(cursor);
   }
#endif
}

// static
void
graphics_info_t::normal_cursor() {

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
      // 20220528-PE FIXME cursor
#else
   if (use_graphics_interface_flag) {
      if (control_key_for_rotate_flag) {
         GdkCursorType c = GDK_LEFT_PTR;
         GdkCursor *cursor;
         //cursor = gdk_cursor_new (c);
         cursor = gdk_cursor_new_for_display (gdk_display_get_default(), c);
         GdkWindow *window = gtk_widget_get_window(glareas[0]);
         gdk_window_set_cursor(window, cursor);
         // gdk_cursor_destroy (cursor);
      }
   }
#endif
}

// static
void
graphics_info_t::watch_cursor() {

}

// static
void
graphics_info_t::fleur_cursor() {

}



// static
short int
graphics_info_t::alt_conf_split_type_number() {

   return graphics_info_t::alt_conf_split_type;

}



// The start of the edit one residue with the pullable point in the
// Ramachandran.  Currently disabled, but will be a nice feature in future.
//
void
graphics_info_t::execute_edit_phi_psi(int atom_index, int imol) {

}

void
graphics_info_t::rama_plot_for_single_phi_psi(int imol, int atom_index) {

}

void
graphics_info_t::rama_plot_for_2_phi_psis(int imol, int atom_index) {

}

// activated from the edit torsion angles cancel button (and OK button, I
// suppose).
//
void
graphics_info_t::destroy_edit_backbone_rama_plot() {  // only one of these.

}



// This is called as part of the callback of moving a point in the ramachandran
// edit window
//
// We (Kevin and I) currently don't like the way the atoms move, so the button
// is invisible now.  However, it will be reinstated in the future, when we
// change residue_edit_phi_psi to use torsion restraints.  residue_edit_phi_psi
// should return a double which is the distortion at the end of the refinement.
// The box should be coloured according to distortion.  We should add a help
// button underneath the canvas (perhaps not visible unless made so, when we are
// in a edit-one-residue situation) [25Feb2004].
//
void
graphics_info_t::set_edit_phi_psi_to(double phi, double psi) {

#ifdef HAVE_GOOCANVAS
   // tinker with the coordinates of the moving_atoms_asc

   short int istat = molecules[imol_moving_atoms].residue_edit_phi_psi(*moving_atoms_asc, edit_phi_psi_atom_index, phi, psi);

   if (istat) {
      // You would do this much less heay-weightedly if you have your
      // bonding already sorted out here.
      //
      int do_disulphide_flag = 0;
      Bond_lines_container bonds(*moving_atoms_asc, do_disulphide_flag);
      regularize_object_bonds_box.clear_up();
      regularize_object_bonds_box = bonds.make_graphical_bonds();
      graphics_draw();
   }
#endif
}


// There are 2 places that have to be united in how they get their torsions.
//
// Here is 1) - (or more precisely in
// wrapped_create_edit_chi_angles_dialog) i.e. where we make create
// the buttons for the torsions.
//
// These torsion labels/indices must match:
// 2) update_residue_by_chi_change()/ chi_angles::change_by()
//
// Particularly we should exclude hydrogen rotations consistently.
// Note that setup_flash_bond_internal() needs to be consistent too.
//
void
graphics_info_t::execute_edit_chi_angles(int atom_index, int imol) {

   // check that we have chis for this residue:
   int n_chis = molecules[imol].N_chis(atom_index);

   // set the static variable for the alt conf of this atom: used when
   // we actually move the atoms.
   chi_angle_alt_conf = molecules[imol].atom_sel.atom_selection[atom_index]->altLoc;

   if (n_chis) {

      std::string res_type(molecules[imol].atom_sel.atom_selection[atom_index]->residue->GetResName());
      chi_angles_clicked_atom_spec =
         coot::atom_spec_t(molecules[imol].atom_sel.atom_selection[atom_index]);
      chi_angles_clicked_atom_spec.int_user_data = 1; // not magic "don't use" value

      // Make Phil Evans happy (well, slightly happier.. :-)
      if (res_type == "MSE")
         chi_angles_clicked_atom_spec.atom_name = " C  ";
      if (res_type == "ARG")
         chi_angles_clicked_atom_spec.atom_name = " C  ";
      if (res_type == "PHE")
         chi_angles_clicked_atom_spec.atom_name = " C  ";
      if (res_type == "TYR")
         chi_angles_clicked_atom_spec.atom_name = " C  ";

      // belt and braces:
      if ( (res_type == "GLY") || (res_type == "ALA") ) {
         std::cout << "This residue does not have chi angles (GLY/ALA)." << std::endl;
      } else {

         // copy the residue, just like we do in execute_edit_phi_psi:
         //
         moving_atoms_asc_type = coot::NEW_COORDS_REPLACE;
         imol_moving_atoms = imol;
         short int whole_res_flag = 0; // We only want to pull the
         // atoms of *this* alternative
         // conformation (and this
         // includes atoms with altconf
         // "").
         atom_selection_container_t residue_asc =
            graphics_info_t::molecules[imol].edit_residue_pull_residue(atom_index,
                                                                       whole_res_flag);

         regularize_object_bonds_box.clear_up();

         edit_chi_edit_type mode = EDIT_CHI;
         int ires = wrapped_create_edit_chi_angles_dialog(res_type, mode);
         if (ires > 0) {
            // std::cout << "Use the 1,2,3,4 keys to select rotamers, 0 for "
            //           << "normal rotation mode" << std::endl;
            make_moving_atoms_graphics_object(imol, residue_asc);

            if (do_probe_dots_on_rotamers_and_chis_flag) {
               setup_for_probe_dots_on_chis_molprobity(imol);
            }
         } else {
            std::cout << "WARNING:: couldn't find torsions in the dictionary "
                      << "for this residue: " << res_type << std::endl;
         }
         graphics_draw();
      }
   } else {
      std::cout << "WARNING:: This residue does not have chi angles." << std::endl;
      std::cout << "Missing dictionary, perhaps? " << std::endl;
      std::string s = "WARNING:: This residue does not have assigned torsions/chi angles.\n";
      s += "Missing dictionary, perhaps?\n";
      info_dialog(s); // checks use_graphics_interface_flag
   }
}

void
graphics_info_t::residue_partial_alt_locs_split_residue(int i_bond, bool wag_the_dog) {

   if (is_valid_model_molecule(imol_residue_partial_alt_locs)) {
      molecules[imol_residue_partial_alt_locs].residue_partial_alt_locs_split_residue(residue_partial_alt_locs_spec, i_bond, residue_partial_alt_locs_rotate_fragment_angle, wag_the_dog, geom_p);
   }
}


void
graphics_info_t::setup_for_probe_dots_on_chis_molprobity(int imol) {

   if (moving_atoms_asc->n_selected_atoms) {

      // make a directory where the probe dots will go:
      int dir_status = coot::util::create_directory("coot-molprobity");

      int n_atoms = moving_atoms_asc->n_selected_atoms;
      molecules[imol].atom_sel.mol->WritePDBASCII("molprobity-tmp-reference-file.pdb");

      // so find the centre and radius of the set of moving atoms:
      coot::Cartesian acc(0,0,0);
      for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
    coot::Cartesian pt(moving_atoms_asc->atom_selection[i]->x,
       moving_atoms_asc->atom_selection[i]->y,
       moving_atoms_asc->atom_selection[i]->z);
    acc += pt;
      }
      coot::Cartesian av(acc.x()/float(n_atoms),
    acc.y()/float(n_atoms),
    acc.z()/float(n_atoms));
      probe_dots_on_chis_molprobity_centre = av;
      float max_d = 0;
      for (int i=0; i<moving_atoms_asc->n_selected_atoms; i++) {
    coot::Cartesian pt(moving_atoms_asc->atom_selection[i]->x,
       moving_atoms_asc->atom_selection[i]->y,
       moving_atoms_asc->atom_selection[i]->z);
    float this_d = (pt - av).amplitude();
    if (this_d > max_d)
       max_d = this_d;
      }

      // so we have the maximum radius of the atoms from the centre
      // point, but we should enlargen the because we are swinging chi
      // and we need to make contact to the reference protein atoms:
      probe_dots_on_chis_molprobity_radius = (max_d + 2) * 1.7;
      if (dir_status == 0) { // success
    do_probe_dots_on_rotamers_and_chis();
      }
   }
}


void
graphics_info_t::setup_flash_bond_using_moving_atom_internal(int i_torsion_index) {

   // turn it off first, only enable it if we find a pair:
   draw_chi_angle_flash_bond_flag = 0; // member data item

   // std::cout << "flash bond i_torsion_index: " << i_torsion_index << std::endl;

   // get the residue type and from that the atom name pairs:
   //
   if (! moving_atoms_asc) {
      // std::cout << "ERROR:: moving_atoms_asc is NULL" << std::endl;
      logger.log(log_t::ERROR, logging::function_name_t(__FUNCTION__), "moving_atoms_asc is NULL");
   } else {
      if (moving_atoms_asc->n_selected_atoms == 0) {
         std::cout << "ERROR: no atoms in moving_atoms_asc" << std::endl;
      } else {
         mmdb::Model *model_p = moving_atoms_asc->mol->GetModel(1);
         if (model_p) {
            mmdb::Chain *chain_p = model_p->GetChain(0);
            if (chain_p) {
               mmdb::Residue *residue_p = chain_p->GetResidue(0);
               if (residue_p) {

                  std::string residue_type(residue_p->GetResName());

                  std::pair<std::string, std::string> atom_names;

                  std::pair<short int, coot::dictionary_residue_restraints_t> r =
                     geom_p->get_monomer_restraints(residue_type, imol_moving_atoms);

                  if (r.first) {
                     std::vector <coot::dict_torsion_restraint_t> torsion_restraints =
                        r.second.get_non_const_torsions(find_hydrogen_torsions_flag);

                     if (i_torsion_index >= 0 && i_torsion_index < int(torsion_restraints.size())) {

                        atom_names.first  = torsion_restraints[i_torsion_index].atom_id_2_4c();
                        atom_names.second = torsion_restraints[i_torsion_index].atom_id_3_4c();

                        if ((atom_names.first != "") &&
                            (atom_names.second != "")) {

                           mmdb::PPAtom residue_atoms;
                           int nResidueAtoms;
                           residue_p->GetAtomTable(residue_atoms, nResidueAtoms);

                           if (nResidueAtoms > 0) { // of course it is!
                              for (int iat1=0; iat1<nResidueAtoms; iat1++) {
                                 std::string ra1=residue_atoms[iat1]->name;
                                 if (ra1 == atom_names.first) {
                                    for (int iat2=0; iat2<nResidueAtoms; iat2++) {
                                       std::string ra2=residue_atoms[iat2]->name;
                                       if (ra2 == atom_names.second) {

                                          draw_chi_angle_flash_bond_flag = 1;
                                          clipper::Coord_orth p1(residue_atoms[iat1]->x,
                                                                 residue_atoms[iat1]->y,
                                                                 residue_atoms[iat1]->z);
                                          clipper::Coord_orth p2(residue_atoms[iat2]->x,
                                                                 residue_atoms[iat2]->y,
                                                                 residue_atoms[iat2]->z);

                                          std::pair<clipper::Coord_orth, clipper::Coord_orth> cp(p1, p2);
                                          graphics_info_t g;
                                          g.add_flash_bond(cp);
                                          graphics_draw();
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
         }
      }
   }
}

void
graphics_info_t::add_flash_bond(const std::pair<clipper::Coord_orth, clipper::Coord_orth> &bond) {

   draw_chi_angle_flash_bond_flag = 1;
   flash_bond = bond;
}





// and the pull restraint neighbour displacement radius (mabye)
void graphics_info_t::draw_pull_restraint_neighbour_displacement_max_radius_circle() {

   if (control_is_pressed) {
      if (pull_restraint_neighbour_displacement_max_radius > 1.0) {
         if (moving_atoms_displayed_p()) {

#if 0 // OpenGL v1
            // there should be a function for this?
            clipper::Coord_orth rc(RotationCentre_x(), RotationCentre_y(), RotationCentre_z());

            coot::Cartesian centre = unproject_xyz(0, 0, 0.5);
            coot::Cartesian front  = unproject_xyz(0, 0, 0.0);
            coot::Cartesian right  = unproject_xyz(1, 0, 0.5);
            coot::Cartesian top    = unproject_xyz(0, 1, 0.5);

            coot::Cartesian screen_x = (right - centre);
            coot::Cartesian screen_y = (top   - centre);
            coot::Cartesian screen_z = (front - centre);

            screen_x.unit_vector_yourself();
            screen_y.unit_vector_yourself();
            screen_z.unit_vector_yourself();

            clipper::Coord_orth screen_x_co(screen_x.x(), screen_x.y(), screen_x.z());
            clipper::Coord_orth screen_y_co(screen_y.x(), screen_y.y(), screen_y.z());

            glColor3f(0.6, 0.6, 0.6);
            glLineWidth(3.0);
            glBegin(GL_LINES); // call this after unproject_xyz().
            for (unsigned int i=0; i<50; i++) {
               float theta_this = 0.02 * static_cast<float>(i)   * M_PI * 2.0;
               float theta_next = 0.02 * static_cast<float>(i+1) * M_PI * 2.0;
               float stt = sinf(theta_this);
               float ctt = cosf(theta_this);
               float stn = sinf(theta_next);
               float ctn = cosf(theta_next);
               float r = pull_restraint_neighbour_displacement_max_radius;

               clipper::Coord_orth p1 = r * ctt * screen_x_co + r * stt * screen_y_co;
               clipper::Coord_orth p2 = r * ctn * screen_x_co + r * stn * screen_y_co;

               p1 += rc;
               p2 += rc;

               glVertex3f(p1.x(), p1.y(), p1.z());
               glVertex3f(p2.x(), p2.y(), p2.z());
            }
            glEnd();
#endif
         }
      }
   }
}

// Question to self? Have I locked the restraints before I call this?
//
void
graphics_info_t::clear_all_atom_pull_restraints(bool refine_again_flag) {

   //std::cout << "debug:: in clear_all_atom_pull_restraints() " << refine_again_flag << std::endl;

   all_atom_pulls_off();
   if (last_restraints) {
      last_restraints->clear_all_atom_pull_restraints();
      if (refine_again_flag)
         drag_refine_refine_intermediate_atoms();
   }
}

// this is not static
void
graphics_info_t::clear_atom_pull_restraint(const coot::atom_spec_t &spec, bool refine_again_flag) {

   // clear_atom_pull_restraints() is a simple wrapper around this (currently in the header)

   if (last_restraints) {
      last_restraints->clear_atom_pull_restraint(spec);
      atom_pull_off(spec);
      if (refine_again_flag)
         drag_refine_refine_intermediate_atoms();
   }
}

// setup and draw
//
void
graphics_info_t::setup_flash_bond(int imol,
     coot::residue_spec_t spec,
     int i_bond) {

   if (is_valid_model_molecule(imol)) {
      mmdb::Residue *residue_p = molecules[imol].get_residue(spec);
      if (residue_p) {
    std::string residue_type = residue_p->GetResName();
    std::pair<short int, coot::dictionary_residue_restraints_t> r =
       geom_p->get_monomer_restraints(residue_type, imol);

    if (r.first) {
       std::vector <coot::dict_torsion_restraint_t> torsion_restraints =
          r.second.get_non_const_torsions(find_hydrogen_torsions_flag);

       if (i_bond >= 0 && i_bond < int(torsion_restraints.size())) {

          std::pair<std::string, std::string> atom_names;
          atom_names.first  = torsion_restraints[i_bond].atom_id_2_4c();
          atom_names.second = torsion_restraints[i_bond].atom_id_3_4c();

          if ((atom_names.first != "") && (atom_names.second != "")) {

     mmdb::PPAtom residue_atoms;
     int nResidueAtoms;
     residue_p->GetAtomTable(residue_atoms, nResidueAtoms);

     if (nResidueAtoms > 0) {
        for (int iat1=0; iat1<nResidueAtoms; iat1++) {
   std::string ra1=residue_atoms[iat1]->name;
   std::string alt_conf_1 = residue_atoms[iat1]->altLoc;
   if (ra1 == atom_names.first) {
      for (int iat2=0; iat2<nResidueAtoms; iat2++) {
         std::string ra2=residue_atoms[iat2]->name;
         std::string alt_conf_2 = residue_atoms[iat2]->altLoc;
         if (ra2 == atom_names.second) {

    if (alt_conf_1 == alt_conf_2) {

       draw_chi_angle_flash_bond_flag = 1;
       clipper::Coord_orth p1(residue_atoms[iat1]->x,
      residue_atoms[iat1]->y,
      residue_atoms[iat1]->z);
       clipper::Coord_orth p2(residue_atoms[iat2]->x,
      residue_atoms[iat2]->y,
      residue_atoms[iat2]->z);

       if (false)
          std::cout << "flash bond: "
    << coot::atom_spec_t(residue_atoms[iat1])
    << " - "
    << coot::atom_spec_t(residue_atoms[iat2])
    << std::endl;

       std::pair<clipper::Coord_orth, clipper::Coord_orth> cp(p1, p2);
       add_flash_bond(cp);
       graphics_draw();
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
   }

}

// static
void graphics_info_t::draw_chi_angles_flash_bond() {

#if 0
   if (draw_chi_angle_flash_bond_flag) {
      glLineWidth(10);
      glColor3f(0.3,1.0,0.3);
      glBegin(GL_LINES);
      glVertex3f(graphics_info_t::flash_bond.first.x(),
    graphics_info_t::flash_bond.first.y(),
    graphics_info_t::flash_bond.first.z());
      glVertex3f(graphics_info_t::flash_bond.second.x(),
    graphics_info_t::flash_bond.second.y(),
    graphics_info_t::flash_bond.second.z());
      glEnd();
   }
#endif
}


// ----------------------------------------------------------------------------
//                          map colour stuff
// ----------------------------------------------------------------------------

void
graphics_info_t::set_last_map_colour(double f1, double f2, double f3) const {

#ifndef EMSCRIPTEN // 20220724-PE for now. This should be restored

   // first find the last map:
   int imap = -1;
   for (int i=0; i<n_molecules(); i++) {
      if (molecules[i].has_xmap()) {
         imap = i;
      }
   }

   if (imap == -1) {
      std::cout << "No maps available for the setting of colour" << std::endl;
   } else {
      GdkRGBA colour;
      if (f1 > 1.0) f1 = 1.0;
      if (f2 > 1.0) f2 = 1.0;
      if (f3 > 1.0) f3 = 1.0;
      if (f1 < 0.0) f1 = 0.0;
      if (f2 < 0.0) f2 = 0.0;
      if (f3 < 0.0) f3 = 0.0;
      colour.red   = f1;
      colour.green = f2;
      colour.blue  = f3;

      if (use_graphics_interface_flag) {
         molecules[imap].handle_map_colour_change(colour, swap_difference_map_colours,
                                                  GL_CONTEXT_MAIN,
                                                  get_rotation_centre_co(),
                                                  box_radius_xray);
         if (display_mode_use_secondary_p()) {
            make_gl_context_current(GL_CONTEXT_SECONDARY);
            molecules[imap].handle_map_colour_change(colour, swap_difference_map_colours,
                                                     GL_CONTEXT_SECONDARY,
                                                     get_rotation_centre_co(),
                                                     box_radius_xray);
            make_gl_context_current(GL_CONTEXT_MAIN);
         }
      }
   }
#endif
}

void
graphics_info_t::set_last_map_contour_level(float level) {

   int imap = -1;
   for (int i=0; i<n_molecules(); i++) {
      if (molecules[i].has_xmap()) {
    imap = i;
      }
   }

   if (imap == -1) {
      std::cout << "No maps available for the setting of contour" << std::endl;
   } else {
      molecules[imap].set_contour_level(level);
   }
}

void
graphics_info_t::set_last_map_contour_level_by_sigma(float f) {

   int imap = -1;
   for (int i=0; i<n_molecules(); i++) {
      if (molecules[i].has_xmap()) {
    imap = i;
      }
   }

   if (imap == -1) {
      std::cout << "No maps available for the setting of contour" << std::endl;
   } else {
      molecules[imap].set_contour_level_by_sigma(f);
   }
}


// And turn it on.
void
graphics_info_t::set_last_map_sigma_step(float f) {

   int imap = -1;
   for (int i=0; i<n_molecules(); i++) {
      if (molecules[i].has_xmap()) {
    imap = i;
      }
   }

   if (imap == -1) {
      std::cout << "No maps available for the setting of contour step"
                << std::endl;
   } else {
//       molecules[imap].contour_by_sigma_flag = 1;
//       molecules[imap].contour_sigma_step = f;
      molecules[imap].set_contour_by_sigma_step(f, 1);
   }
}




// ---------------------- geometry objects -----------------------------
void
graphics_info_t::draw_geometry_objects() {

#if 0 // old
   // 20090715 We change the type of distance_object_vec, and attach a
   // molecule from which the distance was made.  Don't display the
   // distance if the molecule corresponding to the start or end point
   // is not displayed.

   int ndist = distance_object_vec->size();
   double dist;
   clipper::Coord_orth text_pos;

   if (ndist > 0) {
      glEnable(GL_LINE_STIPPLE);
      glLineStipple (1, 0x00FF);
      glLineWidth(2.0);
      glColor3f(0.5, 0.8, 0.6);
      glBegin(GL_LINES);

      const std::vector<coot::simple_distance_object_t> &d = *distance_object_vec;
      for (int i=0; i<ndist; i++) {
         if (is_valid_model_molecule(d[i].imol_start)) {
            if (is_valid_model_molecule(d[i].imol_end)) {
               if (molecules[d[i].imol_start].is_displayed_p()) {
                  if (molecules[d[i].imol_end].is_displayed_p()) {
                     glVertex3d( d[i].start_pos.x(),
                                 d[i].start_pos.y(),
                                 d[i].start_pos.z());
                     glVertex3d( d[i].end_pos.x(),
                                 d[i].end_pos.y(),
                                 d[i].end_pos.z());
                  }
               }
            }
         }
      }
      glEnd();
      glDisable(GL_LINE_STIPPLE);

      for (int i=0; i<ndist; i++) {
         if (is_valid_model_molecule(d[i].imol_start)) {
            if (is_valid_model_molecule(d[i].imol_end)) {
               if (molecules[d[i].imol_start].is_displayed_p()) {
                  if (molecules[d[i].imol_end].is_displayed_p()) {
                     text_pos = d[i].start_pos +
                        0.5 * ( d[i].end_pos - d[i].start_pos +
                                clipper::Coord_orth(0.0, 0.1, 0.1));
                     dist = clipper::Coord_orth::length( d[i].start_pos, d[i].end_pos);
                     printString(float_to_string(dist), text_pos.x(), text_pos.y(), text_pos.z());
                  }
               }
            }
         }
      }


   }

   if (dynamic_distances.size() > 0) {
      draw_dynamic_distances();
   }
#endif
}

// static
void
graphics_info_t::draw_dynamic_distances() {

   std::cout << "graphics_info_t:: draw_dynamic_distances() needs to be replaced " << std::endl;

#if 0
   if (dynamic_distances.size() > 0) {
      glLineWidth(2.0);
      glColor3f(0.5, 0.8, 0.6);

      glEnable(GL_LINE_STIPPLE);
      for (unsigned int i=0; i<dynamic_distances.size(); i++) {
    dynamic_distances[i].draw_dynamic_distance();
      }
      glDisable(GL_LINE_STIPPLE);
   }
#endif
}

void
coot::intermediate_atom_distance_t::draw_dynamic_distance() const {

   std::cout << "graphics_info_t:: draw_dynamic_distance() needs to be replaced " << std::endl;

#if 0
   //glEnable(GL_LINE_STIPPLE);
   glBegin(GL_LINES);
   glVertex3d(dynamic_atom->x,
         dynamic_atom->y,
         dynamic_atom->z);

   glVertex3d(static_position.x(),
         static_position.y(),
         static_position.z());
   glEnd();
   // glDisable(GL_LINE_STIPPLE);

   coot::Cartesian at_pt(dynamic_atom->x,
    dynamic_atom->y,
    dynamic_atom->z);

   // The length of the intermediate distance:
   coot::Cartesian text_pos = at_pt + static_position;
   text_pos *= 0.5;
   text_pos += coot::Cartesian(0.0, 0.1, 0.1);
   coot::Cartesian vec_diff = at_pt - static_position;
   float dist = vec_diff.length();
   // glRasterPos3d();
   std::string t = coot::util::float_to_string(dist);
   graphics_info_t::printString(t, text_pos.x(), text_pos.y(), text_pos.z());
#endif
}

// update_pointer_distances() you might say
void
graphics_info_t::make_pointer_distance_objects() {

   auto coord_orth_to_glm = [] (const clipper::Coord_orth &co) {
                               return glm::vec3(co.x(), co.y(), co.z());
                            };

   auto add_pointer_distance_label = [coord_orth_to_glm] (const std::pair<clipper::Coord_orth, clipper::Coord_orth> &coord_pair,
                                                          const double &dist,
                                                          const glm::vec4 &col) {
                                        clipper::Coord_orth mid_point(0.5 * (coord_pair.first + coord_pair.second));
                                        glm::vec3 offset_mid_point =
                                           coord_orth_to_glm(mid_point) + glm::vec3(0.15, 0.05, 0.05);
                                        std::string label_str = float_to_string_using_dec_pl(static_cast<float>(dist), 2);
                                        unsigned char c = 197;
                                        label_str += c;
                                        atom_label_info_t ali(label_str, offset_mid_point, col);
                                        labels_for_pointer_distances.push_back(ali);
                                     };
   

   clipper::Coord_orth cen(rotation_centre_x,
                           rotation_centre_y,
                           rotation_centre_z);

   std::vector<clipper::Coord_orth> distances;

   if (show_pointer_distances_flag) {
      for (int imol=0; imol<n_molecules(); imol++) {
         if (molecules[imol].has_model()) {
            if (molecules[imol].is_displayed_p()) {
               if (molecules[imol].atom_selection_is_pickable()) {
                  std::vector<clipper::Coord_orth> mol_distances =
                     molecules[imol].distances_to_point(cen, pointer_min_dist, pointer_max_dist);
                  if (mol_distances.size() > 0) {
                     // use insert
                     for (unsigned int id=0; id<mol_distances.size(); id++)
                        distances.push_back(mol_distances[id]);
                  }
               }
            }
         }
      }

      pointer_distances_object_vec.clear();
      for (unsigned int id=0; id<distances.size(); id++) {
         pointer_distances_object_vec.push_back(std::pair<clipper::Coord_orth, clipper::Coord_orth> (distances[id], cen));
      }

      // Now add to the meshed-generic-display-object

      gtk_gl_area_attach_buffers(GTK_GL_AREA(glareas[0]));
      mesh_for_pointer_distances.clear();
      labels_for_pointer_distances.clear();
      unsigned int n_segments = 10;
      coot::colour_holder colour_holder(0.5, 0.6, 0.75);
      std::string colour_name = "LightBlue";
      glm::vec4 glm_colour = colour_holder_to_glm(colour_holder);
      float line_width_scale = 1.0;
      Material material;
      for (unsigned int id=0; id<distances.size(); id++) {
         auto coord_pair = std::make_pair(cen, distances[id]);
         mesh_for_pointer_distances.add_dashed_line(colour_holder, colour_name, coord_pair, material,
                                                    line_width_scale, n_segments);

         double bl = std::sqrt((cen - distances[id]).lengthsq());
         add_pointer_distance_label(coord_pair, bl, glm_colour);
      }
   }
}

void
graphics_info_t::clear_pointer_distances() {

   pointer_distances_object_vec.clear();
   graphics_draw();

}

std::ostream&
coot::operator<<(std::ostream &s, simple_distance_object_t o) {

   s << "simple-distance: start-mol " << o.imol_start << " end-mol " << o.imol_end << " "
     << o.start_pos.format() << " " << o.end_pos.format();
   return s;
}


void
graphics_info_t::clear_measure_distances() {

   measure_distance_object_vec.clear();
   mesh_for_measure_distance_object_vec.clear();
   mesh_for_measure_distance_object_vec.setup_buffers(); // above function should do this I think
   labels_for_measure_distances_and_angles.clear();
   graphics_draw();
}

void
graphics_info_t::clear_last_measure_distance() {

   unsigned int n = measure_distance_object_vec.size();
   std::cout << "debug:: graphics_info_t::clear_last_measure_distance() " << n << std::endl;

   if (n > 0) {
      measure_distance_object_vec.pop_back();

      // a hack that will often work.
      if (labels_for_measure_distances_and_angles.size() > 0)
         labels_for_measure_distances_and_angles.pop_back();

      // rebuild the mesh for measure_distance_object_vec

      mesh_for_measure_distance_object_vec.clear();

      Material material;
      glm::vec4 col(0.72, 0.79, 0.72, 1.0); // same as add_measure_distance()

      for (unsigned int i=0; i<measure_distance_object_vec.size(); i++) {
         const auto &sdo = measure_distance_object_vec[i];
         mesh_for_measure_distance_object_vec.add_dashed_line(sdo, material, col);
      }
      graphics_draw();
   }
}


void
graphics_info_t::draw_generic_text() {

   // should be const, I think.

#if 0
   if (generic_texts_p->size() > 0 ) {
      // GLfloat pink[3] =  { 1.0, 0.8, 0.8 };
      GLfloat pink[3] =  { font_colour.red, font_colour.green, font_colour.blue };
      glColor3fv(pink);
      glPushAttrib (GL_LIST_BIT);
      void *font = graphics_info_t::atom_label_font;
      font = GLUT_BITMAP_TIMES_ROMAN_24;
      glDisable(GL_FOG);
      for (unsigned int i=0; i<generic_texts_p->size(); i++) {
    const coot::generic_text_object_t &gto = (*generic_texts_p)[i];
    glRasterPos3f(gto.x, gto.y, gto.z);
    for (unsigned int is = 0; is < gto.s.length(); is++)
       glutBitmapCharacter (font, gto.s[is]);
      }
      glEnable(GL_FOG);
      glPopAttrib ();
   }
#endif
}

// static
coot::colour_holder
coot::old_generic_display_object_t::colour_values_from_colour_name(const std::string &c) {

   coot::colour_holder colour;
   colour.red = 0.4;
   colour.green = 0.4;
   colour.blue = 0.4;

   if (c.length() == 7) {
      if (c[0] == '#') {
         return coot::colour_holder(c); // hex colour string
      }
   }

   if (c == "blue") {
      colour.red = 0.1; colour.green = 0.1;
      colour.blue = 0.8;
   } else {
      if (c == "sky") {
         colour.red = 0.53 * 0.6;
         colour.green = 0.81 * 0.6;
         colour.blue = 0.92 * 0.6;
      } else {
         if (c == "green") {
            colour.red   = 0.05;
            colour.green = 0.8;
            colour.blue  = 0.05;
         } else {
            if (c == "greentint") {
               colour.red = 0.45;
               colour.green = 0.63;
               colour.blue = 0.45;
            } else {
               if (c == "sea") {
                  colour.red = 0.1;
                  colour.green = 0.6;
                  colour.blue = 0.6;
               } else {
                  if (c == "yellow") {
                     colour.red = 0.8;
                     colour.green = 0.8;
                     colour.blue = 0.0;
                  } else {
                     if (c == "yellowtint") {
                        colour.red = 0.65;
                        colour.green = 0.65;
                        colour.blue = 0.4;
                     } else {
                        if (c == "orange") {
                           colour.red = 0.9;
                           colour.green = 0.6;
                           colour.blue = 0.1;
                        } else {
                           if (c == "red") {
                              colour.red = 0.9;
                              colour.green = 0.1;
                              colour.blue = 0.1;
                           } else {
                              if (c == "hotpink") {
                                 colour.red = 0.9;
                                 colour.green = 0.2;
                                 colour.blue = 0.6;
                              } else {
                                 if (c == "pink") {
                                    colour.red = 0.9;
                                    colour.green = 0.3;
                                    colour.blue = 0.3;
                                 } else {
                                    if (c == "cyan") {
                                       colour.red = 0.1;
                                       colour.green = 0.7;
                                       colour.blue = 0.7;
                                    } else {
                                       if (c == "aquamarine") {
                                          colour.red = 0.1;
                                          colour.green = 0.8;
                                          colour.blue = 0.6;
                                       } else {
                                          if (c == "forestgreen") {
                                             colour.red   = 0.6;
                                             colour.green = 0.8;
                                             colour.blue  = 0.1;
                                          } else {
                                             if (c == "yellowgreen") {
                                                colour.red   = 0.6;
                                                colour.green = 0.8;
                                                colour.blue  = 0.2;
                                             } else {
                                                if (c == "goldenrod") {
                                                   colour.red   = 0.85;
                                                   colour.green = 0.65;
                                                   colour.blue  = 0.12;
                                                } else {
                                                   if (c == "orangered") {
                                                      colour.red   = 0.9;
                                                      colour.green = 0.27;
                                                      colour.blue  = 0.0;
                                                   } else {
                                                      if (c == "magenta") {
                                                         colour.red   = 0.7;
                                                         colour.green = 0.2;
                                                         colour.blue  = 0.7;
                                                      } else {
                                                         if (c == "cornflower") {
                                                            colour.red   = 0.38;
                                                            colour.green = 0.58;
                                                            colour.blue  = 0.93;
                                                         } else {
                                                            if (c == "royalblue") {
                                                               colour.red   = 0.25;
                                                               colour.green = 0.41;
                                                               colour.blue  = 0.88;
                                                            } else {
                                                               if (c == "darkpurple") {
                                                                  colour.red   = 0.5;
                                                                  colour.green = 0.0;
                                                                  colour.blue  = 0.5;
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
   }

//    std::cout << "debug:: in colour_values_from_colour_name from colour " << c
// 	     << " we assign colour values "
// 	     << colour[0] << " "
// 	     << colour[1] << " "
// 	     << colour[2] << "\n";
   return colour;
}


// ----------------- done generic objects -----------------------------


void
graphics_info_t::remove_all_atom_labels() {

   for (int i=0; i<n_molecules(); i++) {
      if (molecules[i].has_model()) {
    molecules[i].remove_atom_labels();
      }
   }
   graphics_draw();

}


// static
std::string
graphics_info_t::ccp4_defs_file_name() {

#if defined WIN32
#ifdef WINDOWS_MINGW
// BL says:: in my windows it's found in $USERPROFILE
// would guess that's true for other win32 too....
    char *home = getenv("USERPROFILE");
#else
    char *home = getenv("HOMEPATH");
#endif // WINDOWS_MINGW
#else
    char *home = getenv("HOME");
#endif // WIN32

#if defined(WINDOWS_MINGW)|| defined(_MSC_VER)
    std::string path = "/CCP4/windows/directories.def";
#else
    std::string path = "/.CCP4/unix/directories.def";
#endif
   std::string filename = home + path;

   return filename;
}


void
graphics_info_t::add_coordinates_glob_extension(const std::string &extension) {

   coordinates_glob_extensions->push_back(extension);
}


void
graphics_info_t::add_data_glob_extension(const std::string &extension) {
   data_glob_extensions->push_back(extension);

}

void
graphics_info_t::add_map_glob_extension(const std::string &extension) {
   map_glob_extensions->push_back(extension);

}

void
graphics_info_t::add_dictionary_glob_extension(const std::string &extension) {
   dictionary_glob_extensions->push_back(extension);

}

void
graphics_info_t::remove_coordinates_glob_extension(const std::string &extension) {

  std::vector<std::string>::iterator it;
  for (it = coordinates_glob_extensions->begin(); it<coordinates_glob_extensions->end(); ++it) {
    if (*it == extension) {
      coordinates_glob_extensions->erase(it);
      // could put in break here!?
      // avoid since it could happen that you have multiples of same entry
    }
  }
}


void
graphics_info_t::remove_data_glob_extension(const std::string &extension) {

  std::vector<std::string>::iterator it;
  for (it = data_glob_extensions->begin(); it<data_glob_extensions->end(); it++) {
    if (*it == extension) {
      data_glob_extensions->erase(it);
    }
  }
}

void
graphics_info_t::remove_map_glob_extension(const std::string &extension) {

  std::vector<std::string>::iterator it;
  for (it = map_glob_extensions->begin(); it<map_glob_extensions->end(); it++) {
    if (*it == extension) {
      map_glob_extensions->erase(it);
    }
  }
}

void
graphics_info_t::remove_dictionary_glob_extension(const std::string &extension) {

  std::vector<std::string>::iterator it;
  for (it = dictionary_glob_extensions->begin(); it<dictionary_glob_extensions->end(); it++) {
    if (*it == extension) {
      dictionary_glob_extensions->erase(it);
    }
  }
}


void
graphics_info_t::check_chiral_volumes(int imol) {

#ifndef EMSCRIPTEN
   if (imol < n_molecules()) {
      if (molecules[imol].has_model()) {
         // return a pair: first is the residues for which no
         // restraints were found second is a vector of atom specs
         // that violate chiral volume constraint.
         std::pair<std::vector<std::string>, std::vector <coot::atom_spec_t> > v =
         molecules[imol].bad_chiral_volumes();
         GtkWidget *w = wrapped_check_chiral_volumes_dialog(v.second, imol);
         if (w)
         gtk_widget_set_visible(w, TRUE);
         if (v.first.size() != 0) { // bad, there was at least one residue not found in dic.
            GtkWidget *wcc = wrapped_create_chiral_restraints_problem_dialog(v.first);
            gtk_widget_set_visible(wcc, TRUE);
         }
      }
   }
#endif
}


void
graphics_info_t::set_moving_atoms(atom_selection_container_t asc,
     int imol, int new_coords_type) {

   imol_moving_atoms = imol;
   make_moving_atoms_graphics_object(imol, asc);
   moving_atoms_asc_type = new_coords_type;
}

// //   static
// void graphics_info_t::bond_parameters_molecule_menu_item_select(GtkWidget *item, GtkPositionType pos) {

//    graphics_info_t g;
//    g.bond_parameters_molecule = pos;
//    GtkWidget *w = lookup_widget(GTK_WIDGET(item), "bond_parameters_dialog");
//    fill_bond_parameters_internals(w, pos); // pos is imol
// }

#ifndef EMSCRIPTEN
// static
void graphics_info_t::bond_parameters_molecule_combobox_changed(GtkWidget *combobox_molecule, gpointer data) {

   std::cout << "-------------------- bond_parameters_molecule_combobox_changed() "
             << combobox_molecule << std::endl;

   graphics_info_t g;
   int imol = g.combobox_get_imol(GTK_COMBO_BOX(combobox_molecule)); // not static
   bond_parameters_molecule = imol;
   // GtkWidget *w = widget_from_builder("bond_parameters_dialog");
   fill_bond_parameters_internals(combobox_molecule, imol);

}
#endif


void
graphics_info_t::clear_diff_map_peaks() {

   diff_map_peaks->resize(0);
   max_diff_map_peaks = 0;

   // Also need to clear the user data on the buttons of the widget -
   // so I should pass the widget pointer then, shouldn't I?

}



// -------- keyboard rotamer control: ---------
// static
void
graphics_info_t::rotamer_dialog_neighbour_rotamer(int istep) {

   graphics_info_t g;

   GtkWidget *rotamer_dialog = widget_from_builder("rotamer_selection_dialog");
   if (rotamer_dialog) {
      // void *t  = (void *) (gtk_object_get_user_data(GTK_OBJECT(g.rotamer_dialog)));
      // std::cout << "user data: " << t << std::endl;
      int n_rotamers = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(rotamer_dialog), "n_rotamers"));
      // std::cout << "We find " << n_rotamers << " rotamers in the widget\n";
      GtkWidget *button;
      short int ifound_active_button = 0;
      int active_button_number = 0;
      int new_active_button_number;
      for (int i=0; i<n_rotamers; i++) {
         std::string button_name = "rotamer_selection_button_rot_";
         button_name += int_to_string(i);
         // button = lookup_widget(g.rotamer_dialog, button_name.c_str());
         button = widget_from_builder(button_name.c_str());
         if (button) {
            if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button))) {
               ifound_active_button = 1;
               active_button_number = i;
               break;
            }
         } else {
            std::cout << "ERROR:: rotamer button not found " << button_name << std::endl;
         }
      }
      if (ifound_active_button) {
         if (istep == 1) {
            new_active_button_number = active_button_number + 1;
            if (new_active_button_number == n_rotamers) {
               new_active_button_number = 0;
            }
         } else {
            new_active_button_number = active_button_number - 1;
            if (new_active_button_number < 0) {
               new_active_button_number = n_rotamers -1;
            }
         }
         std::string button_name = "rotamer_selection_button_rot_";
         button_name += int_to_string(new_active_button_number);
         // GtkWidget *new_button = lookup_widget(g.rotamer_dialog, button_name.c_str());
         GtkWidget *new_button = widget_from_builder(button_name.c_str());

         std::cout << "GTK-FIXME rotamer_dialog_neighbour_rotamer() gtk_signal_emit_by_name()" << std::endl;
         //gtk_signal_emit_by_name(GTK_OBJECT(new_button), "clicked");

      } else {
         std::cout << "ERROR:: not active rotamer button found " << std::endl;
      }
   }

}

void
graphics_info_t::rotamer_dialog_next_rotamer() {

   graphics_info_t::rotamer_dialog_neighbour_rotamer(+1);
}



// static
void
graphics_info_t::rotamer_dialog_previous_rotamer() {

   graphics_info_t::rotamer_dialog_neighbour_rotamer(-1);
}


// -------- keyboard difference map peak control: ------------
// static
void graphics_info_t::difference_map_peaks_next_peak() {
   graphics_info_t::difference_map_peaks_neighbour_peak(1);
}

// static
void graphics_info_t::difference_map_peaks_previous_peak() {
   graphics_info_t::difference_map_peaks_neighbour_peak(-1);
}

// static
void graphics_info_t::difference_map_peaks_neighbour_peak(int istep) { // could be private

#ifndef EMSCRIPTEN
   graphics_info_t g;
   if (g.difference_map_peaks_dialog) {
      int n_peaks = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(g.difference_map_peaks_dialog), "n_peaks"));
      short int ifound_active_button = 0;
      int active_button_number = -99;     // set later
      int new_active_button_number = -99; // set later
      for (int i=0; i<n_peaks; i++) {
         std::string button_name = "difference_map_peaks_button_";
         button_name +=  int_to_string(i);
         // GtkWidget *button = lookup_widget(g.difference_map_peaks_dialog, button_name.c_str());
         GtkWidget *button = nullptr;
         std::cout << "FIXME in difference_map_peaks_neighbour_peak() set the button correctly" << std::endl;
         if (button) {
            if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button))) {
               ifound_active_button = 1;
               active_button_number = i;
            }
         } else {
            std::cout << "DEBUG:: Failed to find button " << button_name << "\n";
         }
      }
      if (ifound_active_button) {
         if (istep == 1) {
            new_active_button_number = active_button_number +1;
            if (new_active_button_number == n_peaks)
            new_active_button_number = 0;
         } else {
            new_active_button_number = active_button_number - 1;
            if (new_active_button_number < 0)
            new_active_button_number = n_peaks -1;
         }
      }
      std::string button_name = "difference_map_peaks_button_";
      button_name += int_to_string(new_active_button_number);
      // GtkWidget *new_button = lookup_widget(g.difference_map_peaks_dialog, button_name.c_str());
      GtkWidget *new_button = 0;
      std::cout << "FIXME in difference_map_peaks_neighbour_peak() set the button 2 correctly" << std::endl;
      std::cout << "GTK-FIXME difference_map_peaks_neighbour_peak() gtk_signal_emit_by_name() " << std::endl;
         // gtk_signal_emit_by_name(GTK_OBJECT(new_button), "clicked");

   } else {
         std::cout << "ERROR:: difference_map_peaks_neighbour_peak called in error\n";
   }
#endif
}

// static
void
graphics_info_t::checked_waters_next_baddie(int dir) {

#ifndef EMSCRIPTEN
   graphics_info_t g;
   GtkWidget *dialog = g.checked_waters_baddies_dialog;
   if (dialog) {
      int n_baddies = GPOINTER_TO_INT(g_object_get_data(G_OBJECT(dialog), "n_baddies"));
      GtkWidget *button;
      bool ifound_active_button = 0;
      int active_button_number = -99; // set later
      int new_active_button_number = -99; // set later

      for (int i=0; i<n_baddies; i++) {
         std::string button_name = "checked_waters_baddie_button_";
         button_name += int_to_string(i);
         // button = lookup_widget(dialog, button_name.c_str());
         button = nullptr;
         std::cout << "FIXME in checked_waters_next_baddie() set the button correctly " << std::endl;
         if (button) {
            if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(button))) {
               ifound_active_button = 1;
               active_button_number = i;
            }
         } else {
            std::cout << "failed to find button " << button_name
                      << std::endl;
         }
      }
      if (ifound_active_button) {
         if (dir == 1) {
            new_active_button_number = active_button_number + 1;
            if (new_active_button_number == n_baddies) {
               new_active_button_number = 0;
            }
         } else {
            new_active_button_number = active_button_number - 1;
            if (new_active_button_number < 0)
               new_active_button_number = n_baddies - 1;
         }
         std::string active_button_name = "checked_waters_baddie_button_";
         active_button_name += int_to_string(new_active_button_number);
         // GtkWidget *new_active_button = lookup_widget(dialog, active_button_name.c_str());
         GtkWidget *new_active_button = 0;
         std::cout << "----- GTK-FIXME checked_waters_next_baddie() gtk_signal_emit_by_name()" << std::endl;
         std::cout << "FIXME in checked_waters_next_baddie() set the button correctly 2 " << std::endl;
         // gtk_signal_emit_by_name(GTK_OBJECT(new_active_button), "clicked");
      } else {
         std::cout << "active button not found" << std::endl;
      }
   }
#endif
}



#ifdef USE_GUILE
// static
SCM
graphics_info_t::safe_scheme_command(const std::string &scheme_command) {

   // if this happens, it's because reading the PDB file (from the command line) is happening too early
   // i.e. before I've called my_wrap_scm_boot_guile()
   //
   if (! scm_boot_guile_booted) return SCM_BOOL_F;

   if (scheme_command.empty()) {
      std::cout << "ERROR:: in safe_scheme_command() empty scheme_command" << std::endl;
      return SCM_BOOL_F;
   }
   
   // std::cout << "starting safe_scheme_command() with scheme_command " << scheme_command << std::endl;

   // return SCM_BOOL_F;

   // FIXME!
   SCM handler = scm_c_eval_string ("(lambda (key . args) (display (list \"(safe_scheme_command) Error in proc: key: \" key \" args: \" args)) (newline))");

   // I am undecided if I want this or not:
   // std::cout << "debug:: safe_scheme_command(): :" << scheme_command << ":" << std::endl;
   std::string thunk("(lambda() ");
   thunk += scheme_command;
   thunk += " )";
   SCM scm_thunk = SCM_BOOL_F;

   // try/catch does not make flow control come back here when bad thunk.
   //
   scm_thunk = scm_c_eval_string(thunk.c_str());
   SCM v = scm_catch(SCM_BOOL_T, scm_thunk, handler);

   SCM dest = SCM_BOOL_F;
   SCM mess = scm_from_locale_string("scm_catch returns: ~s\n");
   SCM sf = scm_simple_format(dest, mess, scm_list_1(v));
   std::string bad_str = scm_to_locale_string(sf);

   return v;
}
#endif // USE_GUILE

void graphics_info_t::run_user_defined_click_func() {

#if defined USE_GUILE && ! defined WINDOWS_MINGW

   if (scm_is_true(scm_procedure_p(user_defined_click_scm_func))) {
      SCM arg_list = SCM_EOL;
      for (unsigned int i=0; i<user_defined_atom_pick_specs.size(); i++) {
         SCM spec_scm = atom_spec_to_scm(user_defined_atom_pick_specs[i]);
         SCM spec_with_model_num = scm_cons(scm_from_int(user_defined_atom_pick_specs[i].model_number), spec_scm);
         arg_list = scm_cons(spec_with_model_num, arg_list);
      }
      arg_list = scm_reverse(arg_list);

      // what are we running? Print it out.
      SCM dest = SCM_BOOL_F;
      SCM mess = scm_from_locale_string("~s");
      SCM ds = scm_simple_format(dest, mess, scm_list_1(user_defined_click_scm_func));
      SCM da = scm_simple_format(dest, mess, scm_list_1(arg_list));
      // std::cout << "INFO:: run_user_defined_click_func() applying " << scm_to_locale_string(ds) << " on "
      //           << scm_to_locale_string(da) << std::endl;
      logger.log(log_t::INFO, "run_user_defined_click_func() applying", scm_to_locale_string(ds), "on", scm_to_locale_string(da));

      SCM rest = SCM_EOL;
      SCM v = scm_apply_1(user_defined_click_scm_func, arg_list, rest);
   }

#endif // USE_GUILE

#ifdef USE_PYTHON

   if (user_defined_click_py_func) {

      if (!PyCallable_Check(user_defined_click_py_func)) {
         std::cout<<"(PYTHON) ERROR:: user_defined_click function must be callable, is "
                  << user_defined_click_py_func->ob_type->tp_name<<std::endl;
      } else {
         // what are we running? Print it out.
         // std::cout << "INFO:: (py) run_user_defined_click_func() applying > "
         //           << PyEval_GetFuncName(user_defined_click_py_func) << " < on:\n";
         logger.log(log_t::INFO, "(py) run_user_defined_click_func() applying >", PyEval_GetFuncName(user_defined_click_py_func), "< on:");

         PyObject *arg_list_py = PyTuple_New(user_defined_atom_pick_specs.size());
         for (unsigned int i=0; i<user_defined_atom_pick_specs.size(); i++) {
            PyObject *spec_py = atom_spec_to_py(user_defined_atom_pick_specs[i]);
            // we need to add the model number too
            PyObject *model_number_py = PyLong_FromLong(user_defined_atom_pick_specs[i].model_number);
            PyList_Insert(spec_py, 0, model_number_py);

            // continue output from above
            PyObject *fmt = myPyString_FromString("[%i,%i,'%s',%i,'%s','%s','%s']");
            PyObject *msg = PyUnicode_Format(fmt, PyList_AsTuple(spec_py));
            std::cout << "   " << myPyString_AsString(msg) << "\n";
            PyTuple_SetItem(arg_list_py, i, spec_py);
            Py_DECREF(fmt);
            Py_DECREF(msg);
         }

         if (PyTuple_Check(arg_list_py)) {
            if (!PyCallable_Check(user_defined_click_py_func)) {
               std::cout << "WARNING:: python user click function should have been callable." << std::endl;
               std::cout << "WARNING:: Ignoring it." << std::endl;
               return;
            }
            // PyObject *result = PyEval_CallObject(user_defined_click_py_func, arg_list_py);
            PyObject *kwargs = nullptr;
            PyObject *result = PyObject_Call(user_defined_click_py_func, arg_list_py, kwargs);
            PyObject *error_thing = PyErr_Occurred();
            if (! error_thing) {
               std::cout << "No Python error" << std::endl;
            } else {
               std::cout << "ERROR:: while executing py run_user_defined_click_func() a python error occured "
                         << error_thing << std::endl;
               PyObject *type, *value, *traceback;
               PyErr_Fetch(&type, &value, &traceback);
               PyErr_NormalizeException(&type, &value, &traceback);
               PyObject *exception_string = PyObject_Repr(value);
               const char *em = myPyString_AsString(exception_string);
               std::cout << "ERROR:: " << em << std::endl;

               Py_XDECREF(value);
               Py_XDECREF(traceback);
               Py_XDECREF(type);
            }
            Py_DECREF(arg_list_py);
            if (result) {
               Py_DECREF(result);
            }
         } else {
            Py_DECREF(arg_list_py);
            std::cout<<"ERROR:: executing user_defined_click" <<std::endl;
         }
      }
   }

   std::cout << "DEBUG:: --------------- run_user_defined_click_func() --- finished " << std::endl;

#endif // USE_PYTHON
}


#ifdef USE_GUILE
SCM
graphics_info_t::atom_spec_to_scm(const coot::atom_spec_t &spec) const {

   SCM r = SCM_EOL;
   r = scm_cons(scm_from_locale_string(spec.alt_conf.c_str()), r);
   r = scm_cons(scm_from_locale_string(spec.atom_name.c_str()), r);
   r = scm_cons(scm_from_locale_string(spec.ins_code.c_str()), r);
   r = scm_cons(scm_from_int(spec.res_no), r);
   r = scm_cons(scm_from_locale_string(spec.chain_id.c_str()), r);
   r = scm_cons(scm_from_int(spec.int_user_data), r); // not the model number? Urgh (unexpected).
                                                     // Where is this user_data used?

   return r;
}
#endif

#ifdef USE_PYTHON
// lets have it as a tuple not a list
PyObject *
graphics_info_t::atom_spec_to_py(const coot::atom_spec_t &spec) const {

  //  PyObject *r = PyTuple_New(6);
  PyObject *r = PyList_New(6);
  PyList_SetItem(r, 0, PyLong_FromLong(spec.int_user_data));
  PyList_SetItem(r, 1, myPyString_FromString(spec.chain_id.c_str()));
  PyList_SetItem(r, 2, PyLong_FromLong(spec.res_no));
  PyList_SetItem(r, 3, myPyString_FromString(spec.ins_code.c_str()));
  PyList_SetItem(r, 4, myPyString_FromString(spec.atom_name.c_str()));
  PyList_SetItem(r, 5, myPyString_FromString(spec.alt_conf.c_str()));

  return r;
}
#endif

#ifdef USE_GUILE
// static
SCM
graphics_info_t::process_socket_string_waiting() {

   SCM r = SCM_BOOL_F;   // was unitiailized
   if (graphics_info_t::have_socket_string_waiting_flag) {

      graphics_info_t::have_socket_string_waiting_flag = 0; // draw() looks here
      std::string ss = graphics_info_t::socket_string_waiting;

      // really the right way?  Perhaps we should just stick to scheme
      // internals?
      std::vector<std::string> v;
      v.push_back("eval-socket-string");
      v.push_back(coot::util::single_quote(ss));

      graphics_info_t g;
      std::string s = g.state_command(v, coot::STATE_SCM);
      r = safe_scheme_command(s);

   }
   return r;
}
#endif

#ifndef EMSCRIPTEN
// static
gboolean
graphics_info_t::process_socket_string_waiting_bool(gpointer user_data) {

#ifdef USE_GUILE
   if (graphics_info_t::have_socket_string_waiting_flag) {
      graphics_info_t::have_socket_string_waiting_flag = 0; // draw() looks here
      std::string ss = graphics_info_t::socket_string_waiting;

      // try internal evaluation:
      if (1) {
    SCM ss_scm = scm_from_locale_string(ss.c_str());

    std::cout << "DEBUG: evaluting :" << ss << ":" << std::endl;
    // SCM r = safe_scheme_command(ss);
    // SCM r = scm_eval_string(ss_scm);
    scm_eval_string(ss_scm);
    // should store r.
    std::cout << "DEBUG: done evaluating" << std::endl;
      }

      // really the right way?  Perhaps we should just stick to scheme
      // internals?
      if (0) {
    std::vector<std::string> v;
    v.push_back("eval-socket-string");
    v.push_back(coot::util::single_quote(ss));

    graphics_info_t g;
    std::string s = g.state_command(v, coot::STATE_SCM);
    safe_scheme_command(s);
      }
   }
   std::cout << " =============== unsetting mutex lock (scheme version) =========" << std::endl;
   graphics_info_t::socket_string_waiting_mutex_lock = 0; // we're done.  release lock.

#endif // USE_GUILE

   return FALSE;
}
#endif

#ifndef EMSCRIPTEN
// static
gboolean
graphics_info_t::process_socket_python_string_waiting_bool(gpointer user_data) {

#ifdef USE_PYTHON

   if (graphics_info_t::have_socket_python_string_waiting_flag) {
      graphics_info_t::have_socket_python_string_waiting_flag = false; // draw() looks here
      std::string ss = graphics_info_t::socket_python_string_waiting;
      safe_python_command(ss);
   }
#endif // USE_PYTHON
   return FALSE;
}
#endif


// static
std::string
graphics_info_t::backslash_filename(const std::string &s) { // needed for windows?

   std::string r = s;

   for (unsigned int i=0; i<s.size(); i++) {
      if (s[i] == '/')
    r[i] = '\\';
   }
   return r;
}




void
graphics_info_t::render_lsq_plane_atoms() {  // put a blob at atoms in lsq_plane_atom_positions

   if (lsq_plane_atom_positions->size() > 0) {
      glColor3f(0.6, 0.6, 0.9);
      glPointSize(8.0);
      glBegin(GL_POINTS);
      for (unsigned int i=0; i<lsq_plane_atom_positions->size(); i++) {
    glVertex3f((*lsq_plane_atom_positions)[i].x(),
       (*lsq_plane_atom_positions)[i].y(),
       (*lsq_plane_atom_positions)[i].z());
      }
      glEnd();
   }
}

int
graphics_info_t::measure_lsq_plane_deviant_atom(int imol, int atom_index) {

   int r = 0;
   if (molecules[imol].has_model()) {
      mmdb::Atom *at = molecules[imol].atom_sel.atom_selection[atom_index];
      clipper::Coord_orth p(at->x, at->y, at->z);

      if (lsq_plane_atom_positions->size() > 2) {

    graphics_draw();
    std::pair<float,float> d_pair =
       coot::lsq_plane_deviation(*lsq_plane_atom_positions, p);
    float d = d_pair.first;

    std::string s("Atom ");
    s += at->name;
    std::string a(at->altLoc);
    if (a != "") {
       s += ",";
       s += a;
    }
    s += " ";
    s += int_to_string(at->GetSeqNum());
    s += at->GetChainID();
    s += " is ";
    s += float_to_string_using_dec_pl(d, 3);
    s += "A from the least squares plane";
    add_status_bar_text(s);
      } else {
    std::string s("Not enough atoms to find plane");
    std::cout << s << "\n";
    add_status_bar_text(s);
      }
   }
   return r;
}


int
graphics_info_t::add_lsq_plane_atom(int imol, int atom_index) {

   if (molecules[imol].has_model()) {
      mmdb::Atom *at = molecules[imol].atom_sel.atom_selection[atom_index];
      clipper::Coord_orth p(at->x, at->y, at->z);
      std::string s("Added plane atom ");
      s += at->name;
      s += " ";
      s += int_to_string(at->GetSeqNum());
      s += at->GetChainID();
      std::cout << s << std::endl;
      add_status_bar_text(s);
      lsq_plane_atom_positions->push_back(p);
      graphics_draw();
   }
   return 0;
}

int
graphics_info_t::remove_last_lsq_plane_atom() {

   if (lsq_plane_atom_positions->size() > 1) {
      lsq_plane_atom_positions->resize(lsq_plane_atom_positions->size()-1);
      graphics_draw();
   }
   return 0;
}



// molecule_info_class_t draw_bonds() function use this to see if the point is within the
// distance from the screen centre.
// Maybe this is not the best place for this function?
// static
bool
graphics_info_t::is_within_display_radius(const coot::CartesianPair &p) {

   coot::Cartesian c(graphics_info_t::RotationCentre_x(),
        graphics_info_t::RotationCentre_y(),
        graphics_info_t::RotationCentre_z());
   float d_sqrd = graphics_info_t::model_display_radius.second * graphics_info_t::model_display_radius.second;

   coot::Cartesian delta_1 = p.getStart() - c;
   if (delta_1.amplitude_squared() > d_sqrd) {
      return false;
   } else {
      coot::Cartesian delta_2 = p.getFinish() - c;
      return (delta_2.amplitude_squared() <= d_sqrd);
   }

}

// molecule_info_class_t draw_bonds() function use this to see if the point is within the
// distance from the screen centre.
// Maybe this is not the best place for this function?
// static
bool
graphics_info_t::is_within_display_radius(const coot::Cartesian &p) {

   coot::Cartesian c(graphics_info_t::RotationCentre_x(),
        graphics_info_t::RotationCentre_y(),
        graphics_info_t::RotationCentre_z());
   float d_sqrd = graphics_info_t::model_display_radius.second * graphics_info_t::model_display_radius.second;

   coot::Cartesian delta = p - c;
   return (delta.amplitude_squared() <= d_sqrd);

}


void
graphics_info_t::set_merge_molecules_ligand_spec(const coot::residue_spec_t &spec_in) {

   merge_molecules_ligand_spec = spec_in;
}



/*! \brief Calculate structure factors from the model and update the given difference
           map accordingly */
void
graphics_info_t::sfcalc_genmap(int imol_model,
                               int imol_map_with_data_attached,
                               int imol_updating_difference_map) {

   // I am keen for this function to be fast - so that it can be used with cryo-EM structures
   //
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map_with_data_attached)) {
         if (true) {
            if (is_valid_map_molecule(imol_updating_difference_map)) {
               if (molecules[imol_updating_difference_map].is_difference_map_p()) {
                  clipper::Xmap<float> *xmap_p = &molecules[imol_updating_difference_map].xmap;
                  try {
                     if (! on_going_updating_map_lock) {
                        on_going_updating_map_lock = true;
                        float cls = molecules[imol_updating_difference_map].get_contour_level_by_sigma();
                        molecules[imol_map_with_data_attached].fill_fobs_sigfobs();
                        const clipper::HKL_data<clipper::data32::F_sigF> *fobs_data =
                           molecules[imol_map_with_data_attached].get_original_fobs_sigfobs();
                        const clipper::HKL_data<clipper::data32::Flag> *free_flag =
                           molecules[imol_map_with_data_attached].get_original_rfree_flags();
                        if (fobs_data && free_flag) {
                           molecules[imol_model].sfcalc_genmap(*fobs_data, *free_flag, xmap_p);
                           molecules[imol_updating_difference_map].set_mean_and_sigma(false, ignore_pseudo_zeros_for_map_stats);
                           molecules[imol_updating_difference_map].set_contour_level_by_sigma(cls); // does an update
                           fill_difference_map_peaks_button_box();
                        }
                        on_going_updating_map_lock = false;
                     } else {
                        std::cout << "DEBUG:: on_going_updating_map_lock was set! - aborting map update." << std::endl;
                     }
                     graphics_draw();
                  }
                  catch (const std::runtime_error &rte) {
                     std::cout << rte.what() << std::endl;
                  }
               }
            }
         }
      }
   }
}


// household function - maybe make a separate file for this sort of function?
void
graphics_info_t::delete_pointers_to_map_in_other_molecules(int imol_map) {

   // std::cout << "---------------------------------------- delete_pointers_to_map_in_other_molecules " << imol_map << std::endl;

   if (is_valid_map_molecule(imol_map)) { // it is at the moment, not for long though!
      clipper::Xmap<float> *xmap_p = &molecules[imol_map].xmap;
      for (int i=0; i<n_molecules(); i++) {
         if (is_valid_map_molecule(i)) {
            if (molecules[i].other_map_for_colouring_p) {
               if (molecules[i].other_map_for_colouring_p == xmap_p) {
                  molecules[i].turn_off_other_map_for_colouring();
               }
            }
         }
      }
   }
}

coot::util::sfcalc_genmap_stats_t
graphics_info_t::sfcalc_genmaps_using_bulk_solvent(int imol_model,
                                                   int imol_map_with_data_attached,
                                                   clipper::Xmap<float> *xmap_2fofc_p, // 2mFo-DFc I mean, of course
                                                   clipper::Xmap<float> *xmap_fofc_p) {

   coot::util::sfcalc_genmap_stats_t stats;
   if (is_valid_model_molecule(imol_model)) {
      if (is_valid_map_molecule(imol_map_with_data_attached)) {
         try {
            if (! on_going_updating_map_lock) {
               on_going_updating_map_lock = true;
               molecules[imol_map_with_data_attached].fill_fobs_sigfobs();

               // 20210815-PE used to be const reference (get_original_fobs_sigfobs() function changed too)
               // const clipper::HKL_data<clipper::data32::F_sigF> &fobs_data = molecules[imol_map_with_data_attached].get_original_fobs_sigfobs();
               // const clipper::HKL_data<clipper::data32::Flag> &free_flag = molecules[imol_map_with_data_attached].get_original_rfree_flags();
               // now the full object (40us for RNAse test).
               // 20210815-PE OK, the const reference was not the problem. But we will leave it as it is now, for now.
               //
               clipper::HKL_data<clipper::data32::F_sigF> *fobs_data_p = molecules[imol_map_with_data_attached].get_original_fobs_sigfobs();
               clipper::HKL_data<clipper::data32::Flag>   *free_flag_p = molecules[imol_map_with_data_attached].get_original_rfree_flags();

               if (fobs_data_p && free_flag_p) {

                  if (true) {

                     // sanity check data
                     const clipper::HKL_info &hkls_check = fobs_data_p->base_hkl_info();
                     const clipper::Spacegroup &spgr_check = hkls_check.spacegroup();
                     const clipper::Cell &cell_check = fobs_data_p->base_cell();
                     const clipper::HKL_sampling &sampling_check = fobs_data_p->hkl_sampling();

                     std::cout << "DEBUG:: in sfcalc_genmaps_using_bulk_solvent() imol_map_with_data_attached "
                               << imol_map_with_data_attached << std::endl;

                     std::cout << "DEBUG:: Sanity check in graphics_info_t:sfcalc_genmaps_using_bulk_solvent(): HKL_info: "
                               << "base_cell: " << cell_check.format() << " "
                               << "spacegroup: " << spgr_check.symbol_xhm() << " "
                               << "sampling-is-null?: " << sampling_check.is_null() << " "
                               << "resolution: " << hkls_check.resolution().limit() << " "
                               << "invsqreslim: " << hkls_check.resolution().invresolsq_limit() << " "
                               << "num_reflections: " << hkls_check.num_reflections()
                               << std::endl;
                  }

                  stats = molecules[imol_model].sfcalc_genmaps_using_bulk_solvent(*fobs_data_p, *free_flag_p, xmap_2fofc_p, xmap_fofc_p);

                  // call fill_difference_map_peaks_button_box() from here like above.

               } else {
                  std::cout << "ERROR:: null data pointer in graphics_info_t::sfcalc_genmaps_using_bulk_solvent() " << std::endl;
               }
               on_going_updating_map_lock = false;
            }
         }
         catch (const std::runtime_error &rte) {
            std::cout << rte.what() << std::endl;
         }
      }
   }
   return stats;
}

#include "utils/xdg-base.hh"

// static
void graphics_info_t::ephemeral_overlay_label_from_id(const std::string &overlay_label) {

   GtkWidget *w = widget_from_builder(overlay_label.c_str());
   if (w) {
      gtk_widget_set_visible(w, TRUE);

      auto label_callback = +[] (gpointer user_data) {
         GtkWidget *w = GTK_WIDGET(user_data);
         gtk_widget_set_visible(w, FALSE);
         return 0;
      };
      g_timeout_add(2000, G_SOURCE_FUNC(label_callback), w);
   }
}

// and the generalization of that! Just pass the text of the ephemeral overlay label
// static
void graphics_info_t::ephemeral_overlay_label(const std::string &overlay_label_text) {

   GtkWidget *w = widget_from_builder("general_use_overlay_label");
   if (w) {
      gtk_widget_set_visible(w, TRUE);
      gtk_label_set_text(GTK_LABEL(w), overlay_label_text.c_str());

      auto label_callback = +[] (gpointer user_data) {
         GtkWidget *w = GTK_WIDGET(user_data);
         gtk_widget_set_visible(w, FALSE);
         return 0;
      };
      g_timeout_add(2000, G_SOURCE_FUNC(label_callback), w);
   }
}


void
graphics_info_t::quick_save() {

   std::cout << "Quick Save!" << std::endl;

   for (int imol=0; imol<n_molecules(); imol++)
      molecules[imol].quick_save();

   short int il = coot::SCRIPT_UNSET;

   xdg_t xdg;
   std::filesystem::path path;

#ifdef USE_GUILE
   il = coot::SCHEME_SCRIPT;
   path = xdg.get_state_home().append(save_state_file_name);
   save_state_file(path.string(), il);
#endif

#ifdef USE_PYTHON
   il = coot::PYTHON_SCRIPT;
   path = xdg.get_state_home().append("0-coot.state.py");
   save_state_file(path.string(), il);
#endif

   add_status_bar_text("Quick Saved");

   GtkWidget *w = widget_from_builder("session_saved_label");
   if (w) {
      gtk_widget_set_visible(w, TRUE);

      auto label_callback = +[] (gpointer user_data) {
         GtkWidget *w = GTK_WIDGET(user_data);
         gtk_widget_set_visible(w, FALSE);
         return 0;
      };
      g_timeout_add(2000, G_SOURCE_FUNC(label_callback), w);
   }
}


// run glColor3f())
// static
void
graphics_info_t::set_bond_colour_from_user_defined_colours(int icol) {

   std::cout << "Don't call this function " <<  __FUNCTION__ << std::endl;
}

// static
void
graphics_info_t::set_user_defined_colours(const std::vector<std::pair<unsigned int, coot::colour_holder> > &user_defined_colours_in) {

   user_defined_colours = user_defined_colours_in;

#if 0
   // (2026-02-04-PE I don't understand what this does)
   // texture colours:
   if (! user_defined_colours.empty()) {
      std::vector<glm::vec4> t_cols(user_defined_colours.size());
      for (unsigned int i=0; i<user_defined_colours.size(); i++) {
         unsigned int idx = user_defined_colours[i].first;
         const auto &col  = user_defined_colours[i].second;
         float alpha = 1.0; // put alpha into coot::colour_holder
         if (idx < t_cols.size())
            t_cols[idx] = glm::vec4(col.red, col.green, col.blue, alpha);
         else
            std::cout << "ERROR:: idx out of range in set_user_defined_colours() " << idx << std::endl;
      }
      texture_for_hud_colour_bar = Texture(400, 200, t_cols, 5);
   }
#endif
}


// static
void graphics_info_t::print_user_defined_colour_table() {

   if (user_defined_colours.empty()) {
      // std::cout << "INFO:: no user-defined colours" << std::endl;
      logger.log(log_t::INFO, "no user-defined colours");
      return;
   }

   for (unsigned int i=0; i<user_defined_colours.size(); i++) {
      unsigned int icol = user_defined_colours[i].first;
      const auto &col   = user_defined_colours[i].second;
      std::cout << "   user-defined colour-table: " << icol << " col " << col << std::endl;
   }
}



// static
void
graphics_info_t::check_keyboard_history_for_easter_egg_codes() {

   std::vector<std::pair<unsigned int, int> > idkfa_pairs = {std::make_pair(5, GDK_KEY_A),
                                                             std::make_pair(4, GDK_KEY_F),
                                                             std::make_pair(3, GDK_KEY_K),
                                                             std::make_pair(2, GDK_KEY_D),
                                                             std::make_pair(1, GDK_KEY_I), };

   size_t l = keyboard_key_history.size();
   if (l >= idkfa_pairs.size()) {

      bool all_matched = true;
      for (const auto &item : idkfa_pairs) {
         if (keyboard_key_history[l-item.first].gdk_key != item.second) {
            all_matched = false;
            break;
         }
      }
      if (all_matched) {
      }
   }
}

GtkWidget *
graphics_info_t::wrapped_create_display_control_window() {

   GtkWidget *widget = widget_from_builder("display_control_window_glade");
   // 20220808-PE unhide the dialog here maybe.
   return widget;
}

//static
void
graphics_info_t::update_symmetry() { // of models

   for (int i=0; i<n_molecules(); i++) {
      if (is_valid_model_molecule(i)) {
         molecules[i].update_symmetry();
      }
   }
}


//static
GdkRGBA
graphics_info_t::symmetry_colour_to_rgba() {

   GdkRGBA rgba;
   rgba.red   = symmetry_colour.r;
   rgba.green = symmetry_colour.g;
   rgba.blue  = symmetry_colour.b;
   rgba.alpha = symmetry_colour.a;

   if (rgba.red   < 0.0) rgba.red   = 0.0;
   if (rgba.green < 0.0) rgba.green = 0.0;
   if (rgba.blue  < 0.0) rgba.blue  = 0.0;
   if (rgba.alpha < 0.0) rgba.alpha = 0.0;

   if (rgba.red   > 1.0) rgba.red   = 1.0;
   if (rgba.green > 1.0) rgba.green = 1.0;
   if (rgba.blue  > 1.0) rgba.blue  = 1.0;
   if (rgba.alpha > 1.0) rgba.alpha = 1.0;

   return rgba;
}

//static
void
graphics_info_t::rgba_to_symmetry_colour(GdkRGBA rgba) {

   symmetry_colour.r = rgba.red;
   symmetry_colour.g = rgba.green;
   symmetry_colour.b = rgba.blue;
   symmetry_colour.a = rgba.alpha;

}

void graphics_info_t::hide_vertical_validation_frame_if_appropriate() {

   // Paul style:
   auto get_n_children = [] (GtkWidget *box) {
      int n_children = 0;
      GtkWidget *item_widget = gtk_widget_get_first_child(box);
      while (item_widget) {
         n_children++;
         item_widget = gtk_widget_get_next_sibling(item_widget);
      }
      return n_children;
   };

   GtkWidget *vbox = widget_from_builder("validation_boxes_vbox");
   // Jakub style:  :-)
   bool should_show_vbox = false;
   for (GtkWidget *w = gtk_widget_get_first_child(vbox); w != nullptr; w = gtk_widget_get_next_sibling(w)) {
      if (gtk_widget_get_visible(w)) {
         should_show_vbox = true;
      }
   }

   GtkWidget *scrolled        = widget_from_builder("ramachandran_plots_scrolled_window");
   GtkWidget *rama_plots_vbox = widget_from_builder("ramachandran_plots_vbox");
   int n_children = get_n_children(rama_plots_vbox);

   // 20230910-PE I don't think that this is right
   // bool rama_plot_shown = gtk_widget_get_visible(scrolled);
   bool rama_plot_shown = false;
   if (n_children > 0) rama_plot_shown = true;

   bool should_hide = !rama_plot_shown && !should_show_vbox;

   std::cout << "here in hide_vertical_validation_frame_if_appropriate rama_plot_shown : " << rama_plot_shown << std::endl;
   std::cout << "here in hide_vertical_validation_frame_if_appropriate should_show_vbox : " << should_show_vbox << std::endl;
   std::cout << "here in hide_vertical_validation_frame_if_appropriate should_hide: " << should_hide << std::endl;

   if(should_hide) {
      GtkWidget* pane = widget_from_builder("main_window_ramchandran_and_validation_pane");
      gtk_widget_set_visible(pane, FALSE);
   }
}


void
graphics_info_t::update_scroll_wheel_map_on_molecule_close() {

   int imol_start = scroll_wheel_map;

   if (is_valid_map_molecule(imol_start)) {
      // nothing need be done
   } else {
      bool changed = false; // to a higher number
      int m = molecules.size() - 1;
      for(int imol=m; imol>=0; imol--) {
         if (imol > imol_start) {
            if (is_valid_map_molecule(imol)) {
               scroll_wheel_map = imol;
               changed = true;
            }
         } else {
            if (! changed) {
               if (is_valid_map_molecule(imol))
                  scroll_wheel_map = imol;
            }
         }
      }
      // nothing was satisfactory then.
      scroll_wheel_map = -1;
   }
}

int
graphics_info_t::get_n_pressed_for_leftquote_tap(std::chrono::time_point<std::chrono::high_resolution_clock> tp) {

   unsigned int s = leftquote_press_times.size();
   unsigned int r = s % 5 + 1;
   if (s != 0) {
      auto tpl = leftquote_press_times.back();
      auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp - tpl).count();
      if (d10 > 2000) {
         leftquote_press_times.clear();
         r = 1;
      }
   }
   leftquote_press_times.push_back(tp);
   return r;
}

void
graphics_info_t::display_next_map() { // one at a time, all, none.

   auto find_index = [] (const std::vector<int> &v, int item) {
      int idx = -1;
      for (size_t i = 0; i < v.size(); i++) {
         if (v[i] == item) {
            idx = i;
            break;
         }
      }
      return idx;
   };

   // one, all or none.
   std::vector<int> maps_to_be_displayed;

   std::vector<int> mm;
   std::vector<int> dm;
   int n_mol = molecules.size();
   for (int imol=0; imol<n_mol; imol++) {
      if (is_valid_map_molecule(imol)) {
         mm.push_back(imol);
         if (molecules[imol].is_displayed_p()) {
            dm.push_back(imol);
         }
      }
   }
   if (mm.empty()) return;
   int mm_size = mm.size();

   if (mm.size() > 1) {
      int idx = -1;
      if (! dm.empty()) {
         if (dm.back() == mm.back() && dm.size() == 1) {
            // the last map was displayed, so now we turn off all maps
         } else {
            int imol_first = dm[0];
            if (dm.size() > 1) {
               maps_to_be_displayed.push_back(imol_first);
            } else {
               int idx_first = find_index(mm, imol_first);
               if (idx_first != -1) {
                  int next_index = idx_first + 1;
                  if (next_index >= mm_size) {
                     next_index = 0;
                  }
                  maps_to_be_displayed.push_back(mm[next_index]);
               } else {
                  // this cannot happen
                  maps_to_be_displayed.push_back(mm[0]);
               }
            }
         }
      } else {
         // none of them were displayed, now display all.
         maps_to_be_displayed = mm;
      }

      for (int imol=0; imol<n_mol; imol++) {
         if (std::find(maps_to_be_displayed.begin(), maps_to_be_displayed.end(), imol) == maps_to_be_displayed.end()) {
            molecules[imol].set_map_is_displayed(0);
         } else {
            molecules[imol].set_map_is_displayed(1);
         }
      }
   } else {
      // toggle the display of the one displayed map
      int imol_map = mm[0];
      if (dm.empty())
         molecules[imol_map].set_map_is_displayed(1);
      else
         molecules[imol_map].set_map_is_displayed(0);
   }
}

// static
void
graphics_info_t::toggle_display_of_last_model() {

   int imol = -1;
   int n_mols = n_molecules();
   for (int i=0; i<n_mols; i++) {
      if (is_valid_model_molecule(i))
         imol = i;
   }

   if (imol > -1) {
      if (molecules[imol].is_displayed_p())
         molecules[imol].set_mol_is_displayed(0);
      else
         molecules[imol].set_mol_is_displayed(1);
   }
}
