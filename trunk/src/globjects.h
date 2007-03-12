
#ifndef GLOBJECTS_H
#define GLOBJECTS_H
// contains the c interface for the graphics, i.e. those 
// needed in main.c, used by 

#include <string>
#include <vector>

#include <gtk/gtk.h>

#include "gl-matrix.h"

GtkWidget* gl_extras(GtkWidget* window, short int try_hardware_stero); 


gint init(GtkWidget *widget); 
gint init_gl_widget(GtkWidget *widget);
int make_current_gl_context(GtkWidget *widget);
void gdkglext_finish_frame(GtkWidget *widget);

void init_surface_wrapper(); 


void init_surface(std::string file);


void read_triangles (std::string tri_file);


void draw_triangles();

void init_molecule();

gint draw(GtkWidget *widget, GdkEventExpose *event);
gint draw_mono(GtkWidget *widget, GdkEventExpose *event, short int in_stereo_flag);
gint draw_hardware_stereo(GtkWidget *widget, GdkEventExpose *event);


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
void handle_scroll_event(int scroll_up_down_flag);

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


void setup_lighting(short int do_lighting_flag); 


void draw_surface_as_display_list();

void draw_surface_not_display_list();


void aniso_atom(int i, float x, float y, float z); 


float r_50(std::string ele);

float rad_50_and_prob_to_radius(float rad_50, float prob);

void draw_axes(GL_matrix &m);

void draw_crosshairs_maybe();

void printString(std::string s); 

gint animate_idle(GtkWidget *widget);

void update_things_on_move_and_redraw(); 

// delete me when fixed
void set_skeleton_bond_colour_random(int i, const std::vector< std::vector<float> > &colour_table); 
void myWireCube(float size); 

void keypad_translate_xyz(short int axis, short int direction);

#endif // GLOBJECTS_H
