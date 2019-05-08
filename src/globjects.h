/* src/globjects.h
 * 
 * Copyright 2002, 2003 by The University of York
 * Copyright 2009 by The University of Oxford
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
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifndef GLOBJECTS_H
#define GLOBJECTS_H
// contains the c interface for the graphics, i.e. those 
// needed in main.c, used by 

#include <string>
#include <vector>

#include <gtk/gtk.h>

#include "gl-matrix.h"

#define GRAPHICS_WINDOW_X_START_SIZE 500
#define GRAPHICS_WINDOW_Y_START_SIZE 500


// Gtk2
GtkWidget* gl_extras(GtkWidget* window, short int try_hardware_stero);
// Gtk3
GtkWidget* gl_gtk3_widget(GtkWidget* window, short int try_hardware_stero);

enum { IN_STEREO_MONO = 0, 
       IN_STEREO_HARDWARE_STEREO=1, 
       IN_STEREO_ZALMAN_RIGHT=5, 
       IN_STEREO_ZALMAN_LEFT=6, 
       IN_STEREO_SIDE_BY_SIDE_LEFT=10,
       IN_STEREO_SIDE_BY_SIDE_RIGHT=11};

gint init(GtkWidget *widget); 
gint init_gl_widget(GtkWidget *widget);
void gdkglext_finish_frame(GtkWidget *widget);

void init_surface_wrapper(); 


void init_surface(std::string file);


void read_triangles (std::string tri_file);


void draw_triangles();

void init_molecule();

gint draw(GtkWidget *widget, GdkEventExpose *event);
gint expose(GtkWidget *widget, GdkEventExpose *event);
gint draw_mono(GtkWidget *widget, GdkEventExpose *event, short int in_stereo_flag);
void debug_eye_position(GtkWidget *widget);
gint draw_hardware_stereo(GtkWidget *widget, GdkEventExpose *event);
gint draw_zalman_stereo(GtkWidget *widget, GdkEventExpose *event);
void stereo_projection_setup_maybe(GtkWidget *widget, short int in_stereo_flag);
coot::Cartesian eye_position();

void do_drag_pan(gdouble x, gdouble y, GtkWidget *widget);
void do_button_zoom(gdouble x, gdouble y);
void do_screen_z_rotate(gdouble x, gdouble y);
void do_ztrans_and_clip(gdouble x, gdouble y);
void adjust_clipping(double d);

gint reshape(GtkWidget *widget, GdkEventConfigure *event); 
gint glarea_motion_notify (GtkWidget *widget, GdkEventMotion *event);
gint glarea_button_press(GtkWidget *widget, GdkEventButton *event);
#if (GTK_MAJOR_VERSION == 2) || defined(WINDOWS_MINGW) || defined(_MSC_VER)
gint glarea_scroll_event(GtkWidget *widget, GdkEventScroll *event);
#endif
gint glarea_button_release(GtkWidget *widget, GdkEventButton *event);
gint key_press_event(GtkWidget *widget, GdkEventKey *event);
gint key_release_event(GtkWidget *widget, GdkEventKey *event);
void handle_scroll_density_level_event(int scroll_up_down_flag);

void debug_draw_rotation_axes(float y1, float y2, float x1, float x2); 
void rotate_baton(); 

void set_bond_colour(int i);
// void set_symm_bond_colour(int i);
// void set_skeleton_bond_colour(float f);

void display_density_level_maybe(); 

// colour helper function
// double combine_colour(double v, int col_part_index); 

std::vector<float> rotate_rgb(std::vector<float> &in_vals, 
			      float amount); 
std::vector<float> convert_rgb_to_hsv(const std::vector<float> &rgb);
std::vector<float> convert_hsv_to_rgb(const std::vector<float> &hsv);

//

void setup_for_mol_triangles();

void setup_lighting(short int do_lighting_flag); 
void show_lighting();

void draw_surface_as_display_list();

void draw_surface_not_display_list();


void aniso_atom(int i, float x, float y, float z); 


float r_50(std::string ele);

float rad_50_and_prob_to_radius(float rad_50, float prob);

void draw_axes(GL_matrix &m);
void graphics_ligand_view();

void draw_crosshairs_maybe();

gint animate_idle_spin(gpointer user_data);
gboolean animate_idle_rock(gpointer user_data);
// gint animate_idle_ligand_interactions(GtkWidget *widget);
gboolean animate_idle_ligand_interactions(gpointer data);


void update_things_on_move_and_redraw(); 

// delete me when fixed
void set_skeleton_bond_colour_random(int i, const std::vector< std::vector<float> > &colour_table); 
void myWireCube(float size); 

void keypad_translate_xyz(short int axis, short int direction);

void test_object();

#endif // GLOBJECTS_H
