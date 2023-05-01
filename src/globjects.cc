/* src/globjects.cc
 *
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2006 by Bernhard Lohkamp
 * Copyright 2007 by Paul Emsley
 * Copyright 2008, 2009 by The University of Oxford
 * Copyright 2014, 2016 by Medical Research Council
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

// we need the functions:
//
// draw, reshape, init, mouse_move and mouse_button press
// (and animate(for idle)).

// 20220528-PE question to self - is there actually anything in this file that we now need?
//             .. let's presume not for now
#if 0 

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"


#ifndef NULL
#define NULL 0
#endif

#include <gtk/gtk.h>
#include <gdk/gdkkeysyms.h> // for keyboarding.

#include <epoxy/gl.h>

#include "compat/sleep-fixups.h"

#include <string.h> // strncmp

#include <math.h>
#ifndef HAVE_VECTOR
#include <vector>
#endif // HAVE_VECTOR

#ifndef HAVE_STRING
#include <string>
#endif // HAVE_STRING

#include "drag-and-drop.hh"
#include "interface.h"

#include <mmdb2/mmdb_manager.h>
#include "coords/mmdb-extras.h"
#include "coords/mmdb.hh"
#include "coords/mmdb-crystal.h"

#include "coords/cos-sin.h"
#include "coords/Cartesian.h"
#include "coords/Bond_lines.h"

#include "graphics-info.h"

#include "gl-matrix.h"

// #include "xmap-interface.h" // is this necessary? nope


#include "trackball.h"

#include "c-interface.h" // for toggle idle function
#include "c-interface-gtk-widgets.h"
#include "cc-interface.hh"
#include "cc-interface-scripting.hh"
#include "globjects.h"  // string

#include "coot-database.hh"

#include "rotate-translate-modes.hh"
#include "rotamer-search-modes.hh"

#include "draw.hh"


// Initialize the graphics_info_t mouse positions
// and rotation centre.


GtkWidget *
gl_gtk3_widget(GtkWidget *vbox, short int try_stereo_flag) {

   GtkWidget *drawing_area = gtk_gl_area_new();
   // GtkWidget *main_window_graphics_hbox = lookup_widget(vbox, "main_window_graphics_hbox");
   GtkWidget *main_window_graphics_hbox = 0; // 20220309-PE removing lookup-widget()s
   gtk_container_add(GTK_CONTAINER(main_window_graphics_hbox), drawing_area);
   gtk_widget_set_size_request(drawing_area, 500, 500);
   gtk_widget_show(drawing_area);
   gtk_gl_area_make_current(GTK_GL_AREA(drawing_area));

   return drawing_area;
}

// GTK2 code
//
// if try_stereo_flag is 1, then hardware stereo.
// if try_stereo_flag is 2, the side-by-side stereo
//
GtkWidget *
gl_extras(GtkWidget* vbox1, short int try_stereo_flag) { // rename gl_extras_gtk2

#if 0
   graphics_info_t g;
   GdkGLConfig *glconfig = 0;
   bool got_hardware_stereo_flag = 0;

//    GdkGLConfigMode mode = static_cast<GdkGLConfigMode>
//       (GDK_GL_MODE_RGB    |
//        GDK_GL_MODE_DEPTH  |
//        GDK_GL_MODE_DOUBLE |
//        GDK_GL_MODE_MULTISAMPLE |
//        /* 2x FSAA */
//        (2 << GDK_GL_MODE_SAMPLES_SHIFT)       );

   GdkGLConfigMode mode;

#ifdef GDKGLEXT_HAVE_MODE_SAMPLES_SHIFT

   mode = static_cast<GdkGLConfigMode>
      (GDK_GL_MODE_RGB    |
       GDK_GL_MODE_DEPTH  |
       GDK_GL_MODE_MULTISAMPLE |
       (4 << GDK_GL_MODE_SAMPLES_SHIFT) |
       GDK_GL_MODE_DOUBLE);

#else

   mode = static_cast<GdkGLConfigMode>
      (GDK_GL_MODE_RGB    |
       GDK_GL_MODE_DEPTH  |
       GDK_GL_MODE_MULTISAMPLE |
       GDK_GL_MODE_DOUBLE);

#endif

   if (try_stereo_flag == coot::HARDWARE_STEREO_MODE) {
      mode = static_cast<GdkGLConfigMode>
	 (GDK_GL_MODE_RGB    |
	  GDK_GL_MODE_DEPTH  |
	  GDK_GL_MODE_DOUBLE |
	  GDK_GL_STEREO);
      /* Try double-buffered visual */
      glconfig = gdk_gl_config_new_by_mode(mode);
      if (glconfig) {
	 got_hardware_stereo_flag = 1;// for message later
      } else {
	 std::cout << "WARNING:: Can't enable stereo visual - falling back"
		   << std::endl;
	 mode = static_cast<GdkGLConfigMode>
	    (GDK_GL_MODE_RGB    |
	     GDK_GL_MODE_DEPTH  |
	     GDK_GL_MODE_DOUBLE);
      }
   }

   if (try_stereo_flag == coot::ZALMAN_STEREO) {
     mode = static_cast<GdkGLConfigMode>
       (GDK_GL_MODE_RGB    |
	GDK_GL_MODE_DEPTH  |
	// GDK_GL_MODE_MULTISAMPLE |
	GDK_GL_MODE_STENCIL |
	GDK_GL_MODE_DOUBLE
	);
      /* Try stencil buffer */
      glconfig = gdk_gl_config_new_by_mode(mode);
      if (glconfig == NULL) {
	g_print ("\n*** Cannot use stencil buffer.\n");
	g_print ("\n*** Sorry no Zalman setting possible.\n");
	mode = static_cast<GdkGLConfigMode>
	  (GDK_GL_MODE_RGB    |
	   GDK_GL_MODE_DEPTH  |
	   GDK_GL_MODE_DOUBLE);
      } else {
        graphics_info_t::display_mode = coot::ZALMAN_STEREO;
      }
   }



   /* Try double-buffered visual */
   glconfig = gdk_gl_config_new_by_mode(mode);
   if (glconfig == NULL) {
      g_print ("\n*** Cannot find the double-buffered visual.");
      g_print ("\n*** Trying single-buffered visual.\n\n");

      mode =static_cast<GdkGLConfigMode>
	 (GDK_GL_MODE_RGB   |
	  GDK_GL_MODE_DEPTH);

      /* Try single-buffered visual */
      glconfig = gdk_gl_config_new_by_mode (mode);
      if (glconfig == NULL) {
	 g_print ("*** No appropriate OpenGL-capable visual found.\n");
	 exit (1);
      }
   }

   // BL debugging
   //if (gdk_gl_config_has_stencil_buffer(glconfig)) {
   //  g_print("BL DEBUG:: have stencil\n");
   //} else {
   //  g_print("BL DEBUG:: dont have stencil\n");
   //}

  /*
   * Drawing area to draw OpenGL scene.
   */
   GtkWidget *drawing_area = NULL; // the returned thing

  int n_contexts = 1;
  if ((try_stereo_flag == coot::SIDE_BY_SIDE_STEREO) ||
      (try_stereo_flag == coot::SIDE_BY_SIDE_STEREO_WALL_EYE) ||
      (try_stereo_flag == coot::DTI_SIDE_BY_SIDE_STEREO)) {
     n_contexts = 2;
     graphics_info_t::display_mode = try_stereo_flag;
  }

  int context_count = 0;

  while(context_count < n_contexts) {

     context_count++;

     GtkWidget *drawing_area_tmp = gtk_drawing_area_new();
     if (context_count == 1) {
	drawing_area = drawing_area_tmp;
     } else {
	graphics_info_t::glarea[1] = drawing_area_tmp;
     }

     {
//	int gl_context_x_size = GRAPHICS_WINDOW_X_START_SIZE;
//	int gl_context_y_size = GRAPHICS_WINDOW_Y_START_SIZE;
	int gl_context_x_size = graphics_info_t::graphics_x_size;
	int gl_context_y_size = graphics_info_t::graphics_y_size;
	// int gl_context_y_size = g.graphics_y_size; old

 	if (context_count > 1) { // more than the first context
	   // std::cout << " =============== " << context_count << std::endl;
	  if (graphics_info_t::glareas[0]) {
 	   gl_context_x_size = graphics_info_t::glareas[0]->allocation.width;
 	   gl_context_y_size = graphics_info_t::glareas[0]->allocation.height;
	  }
	   // std::cout << " ===============" << gl_context_x_size
	   // << " "<< gl_context_y_size << std::endl;
 	}

// 	gtk_widget_set_size_request (drawing_area_tmp,
// 				     gl_context_x_size,
// 				     gl_context_y_size);

	// GtkWindow *window1 = GTK_WINDOW(lookup_widget(vbox1,"window1"));
	GtkWindow *window1 = 0; // 20220309-PE so that it compiles at least.
	gtk_window_set_default_size(window1,
				    n_contexts * gl_context_x_size,
				    n_contexts * gl_context_y_size);

     }

     /* Set OpenGL-capability to the widget */
     gtk_widget_set_gl_capability (drawing_area_tmp,
				   glconfig,
				   NULL,
				   TRUE,
				   GDK_GL_RGBA_TYPE); // render_type (only choice)
     if (drawing_area_tmp) {

	if (try_stereo_flag == coot::HARDWARE_STEREO_MODE) {
	   if (got_hardware_stereo_flag) {
	      std::cout << "INFO:: Hardware stereo widget opened successfully"
			<< std::endl;
	      graphics_info_t::display_mode = coot::HARDWARE_STEREO_MODE;
	   }
	}

	/* Events for widget must be set before X Window is created */
	gtk_widget_set_events(GTK_WIDGET(drawing_area_tmp),
			      GDK_EXPOSURE_MASK      |
			      GDK_BUTTON_PRESS_MASK  |
			      GDK_BUTTON_RELEASE_MASK|
			      GDK_SCROLL_MASK        |
			      GDK_POINTER_MOTION_MASK|
			      GDK_KEY_PRESS_MASK     |
			      GDK_KEY_RELEASE_MASK   |
// 			      GDK_DRAG_ENTER         |
// 			      GDK_DRAG_LEAVE         |
// 			      GDK_DRAG_STATUS        |
// 			      GDK_DROP_START         |
// 			      GDK_DROP_FINISHED      |
			      GDK_POINTER_MOTION_HINT_MASK);

	/* Connect signal handlers */
	/* Redraw image when exposed. */
	g_signal_connect(G_OBJECT(drawing_area_tmp), "expose_event",
			 G_CALLBACK(expose), NULL);
	/* When window is resized viewport needs to be resized also. */
	g_signal_connect(G_OBJECT(drawing_area_tmp), "configure_event",
			 G_CALLBACK(reshape), NULL);
	/* Do initialization when widget has been realized. */
	g_signal_connect(G_OBJECT(drawing_area_tmp), "realize",
			 G_CALLBACK(init), NULL);

	/* pressed a button? */
	g_signal_connect (G_OBJECT(drawing_area_tmp), "button_press_event",
			    G_CALLBACK(glarea_button_press), NULL);
	g_signal_connect (G_OBJECT(drawing_area_tmp), "button_release_event",
			  G_CALLBACK(glarea_button_release), NULL);
	/* mouse in motion! */
	g_signal_connect (G_OBJECT(drawing_area_tmp), "motion_notify_event",
			  G_CALLBACK(glarea_motion_notify), NULL);
	// mouse wheel scrolled:
	g_signal_connect (G_OBJECT(drawing_area_tmp), "scroll_event",
			  G_CALLBACK(glarea_scroll_event), NULL);

	/* put glarea into vbox */
	// GtkWidget *main_window_graphics_hbox = lookup_widget(vbox1, "main_window_graphics_hbox");
	GtkWidget *main_window_graphics_hbox = 0; // removing lookup_widget()s.
	gtk_container_add(GTK_CONTAINER(main_window_graphics_hbox), GTK_WIDGET(drawing_area_tmp));

// #endif here? confliced merge

	/* Capture keypress events */
	g_signal_connect(G_OBJECT(drawing_area_tmp), "key_press_event",
			 G_CALLBACK(key_press_event), NULL);
	g_signal_connect(G_OBJECT(drawing_area_tmp), "key_release_event",
			 G_CALLBACK(key_release_event), NULL);

	// setup drag and drop
	int n_dnd_targets = 4;
	GtkTargetEntry target_list[] = {
	   { (gchar *) "text/uri-list", 0, TEXT_URI },
	   { (gchar *) "_NETSCAPE_URL", 0, TEXT_URL },
	   { (gchar *) "STRING",        0, TARGET_STRING },
	   { (gchar *) "text/plain",    0, TARGET_STRING }}; // deprecated conversions from string constant
                                                // to char *.  GTK+ problem AFAICS.
	GtkDestDefaults dest_defaults = GTK_DEST_DEFAULT_ALL; // we don't need ALL, but this
	                                                      // removes int->GtkDestDefaults
	                                                      // conversion problems with |.
	gtk_drag_dest_set(GTK_WIDGET(drawing_area_tmp), /* widget that will accept a drop */
			  dest_defaults,
			  target_list,            /* lists of target to support */
			  n_dnd_targets,
			  GDK_ACTION_COPY);       /* what to do with data after dropped */


	// 2.6? - but... what does it give us?  Something useful?
	// More documenation-reading required...
	//
	gtk_drag_dest_add_uri_targets(GTK_WIDGET(drawing_area_tmp));


	// if something was dropped
        g_signal_connect (GTK_WIDGET(drawing_area_tmp), "drag-drop",
			  G_CALLBACK (on_gl_canvas_drag_drop), NULL);
	// what to we do with dropped data...
	g_signal_connect (GTK_WIDGET(drawing_area_tmp), "drag-data-received",
			  G_CALLBACK(on_drag_data_received), NULL);

	/* set focus to glarea widget - we need this to get key presses. */
	GTK_WIDGET_SET_FLAGS(drawing_area_tmp, GTK_CAN_FOCUS);
	gtk_widget_grab_focus(GTK_WIDGET(drawing_area_tmp));

     } else {
	std::cout << "CATASTROPHIC ERROR:: in gl_extras  no GtkGL widget!"
		  << std::endl;
     }
  } // end while
  //   std::cout << "DEBUG:: gl_extras returns " << drawing_area << std::endl;
  return drawing_area;
#endif

  return 0; // for now
}


gint
init(GtkWidget *widget)
{

   graphics_info_t g;

   // The cosine->sine lookup table, used in picking.
   //
   // The data in it are static, so we can get to them anywhere
   // now that we have run this
   cos_sin cos_sin_table(1000);

   //
   // what about glXMakeCurrent(dpy, Drawable, Context) here?
   // Bool glXMakeCurrent(Display * dpy,
   //                     GLXDrawable  Drawable,
   //                     GLXContext  Context)

   init_gl_widget(widget);

   return GL_TRUE;
}


gint
init_gl_widget(GtkWidget *widget) {

   GtkAllocation allocation;
   gtk_widget_get_allocation(widget, &allocation);
   glViewport(0,0, allocation.width, allocation.height);
   glMatrixMode(GL_PROJECTION);
   glLoadIdentity();
   glOrtho(-10,30, 10,-20, -20,20); // change clipping
   glMatrixMode(GL_MODELVIEW);
   glLoadIdentity();
   graphics_info_t g;

   // set the background colour
   // glClearColor(0,0,0,1);
   // these statics are set at the begining.
   glClearColor(g.background_colour[0],
		g.background_colour[1],
		g.background_colour[2],
		1);

   // Enable fog
   //
   glEnable(GL_FOG);

   // Linear fog
   glFogi(GL_FOG_MODE, GL_LINEAR);
   glFogf(GL_FOG_START, -20.0);
   // glFogf(GL_FOG_START, 0.0);
   glFogf(GL_FOG_END, 20.0);

   // glFogfv(GL_FOG_COLOR, g.background_colour);

   // Exponential fog
   //
   // glFogi(GL_FOG_MODE, GL_EXP2);
   // glFogf(GL_FOG_DENSITY, 0.10);

   //
   glEnable(GL_DEPTH_TEST);

   // glDepthFunc(GL_LESS); what does this do?

   glEnable(GL_POINT_SMOOTH); // circles not squares

   if (g.do_anti_aliasing_flag)
      glEnable(GL_LINE_SMOOTH);

   // glEnable(GL_POLYGON_SMOOTH);

   glEnable(GL_MULTISAMPLE_ARB);

   // Mac play
   // glEnable(GL_LINE_SMOOTH);
   // glEnable(GL_BLEND);
   // glBlendFunc(GL_SRC_ALPHA, GL_ZERO);

   // Solid model lighting
   //
   setup_lighting(graphics_info_t::do_lighting_flag);


   // set the skeleton to initially be yellow
   graphics_info_t::skeleton_colour[0] = 0.7;
   graphics_info_t::skeleton_colour[1] = 0.7;
   graphics_info_t::skeleton_colour[2] = 0.2;


   // last value does not matter because zero rotation
   // trackball(graphics_info_t::quat , 0.0, 0.0, 0.0, 0.0, 0.0);

   // initialize also the baton quaternion
   // trackball(graphics_info_t::baton_quat , 0.0, 0.0, 0.0, 0.0, 0.0);


   // gtk_idle_add((GtkFunction)animate, widget);

   // should be in graphics_info_t?
   //
   setup_molecular_triangles();

   // setup_for_single_triangle();

   return TRUE;
}

void
setup_molecular_triangles() {
   //g goodby innards
}

void
setup_lighting(short int do_lighting_flag) {

   std::cout << "------------------------------------------------------------------------------------------------------------------"
             << " setup_lighting() -------------------------------------------------------------------" << std::endl;

   // I'm not sure that this does anything other than enable the lights.
   // The light positions are properly set in draw_mono().
   // Needs rationalization.
   // FIXME-lighting

   if (do_lighting_flag) { // set this to 1 to light a surface currently.

      // w = 0.0 means directional light
      //
      // GL_LIGHT2 is for cut-glass mode
      //
      glClearColor(0.0, 0.0, 0.0, 0.0);
      glShadeModel(GL_SMOOTH);

      GLfloat light_ambient[]  = { 0.2, 0.2, 0.2, 1.0 };
      GLfloat light_diffuse[]  = { 0.2, 0.2, 0.2, 1.0 };
      GLfloat light_specular[] = { 0.2, 0.2, 0.2, 1.0 };
      GLfloat light_position[] = { 0.2, 0.2, 0.2, 1.0 };

      GLfloat  light_0_position[] = { 1.0,  1.0, 1.0, 0.0};
      GLfloat  light_1_position[] = {-1.0,  0.0, 1.0, 0.0};
      GLfloat  light_2_position[] = { 0.0,  0.0, 0.0, 0.0};

      // Light 1 position is set elsewhere.

      glLightfv(GL_LIGHT0,   GL_POSITION, light_0_position);
//      glLightfv(GL_LIGHT1,   GL_POSITION, light_1_position);
      glLightfv(GL_LIGHT2,   GL_POSITION, light_2_position);

      glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
      glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
      glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
      glLightfv(GL_LIGHT0, GL_POSITION, light_position);

      glLightfv(GL_LIGHT1, GL_AMBIENT,  light_ambient);
      glLightfv(GL_LIGHT1, GL_DIFFUSE,  light_diffuse);
      glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
      glLightfv(GL_LIGHT1, GL_POSITION, light_position);

      glEnable(GL_LIGHT0);
      glEnable(GL_LIGHT1);
      // GLfloat  light_0_position[] = { 1.0,  1.0,  1.0, 1.0};
      // GLfloat  light_1_position[] = {-1.0,  0.0,  1.0, 1.0};
      // GLfloat  light_2_position[] = { 0.0,  0.0, -1.0, 1.0};

      // glLightfv(GL_LIGHT0,   GL_POSITION, light_0_position);
      // glLightfv(GL_LIGHT1,   GL_POSITION, light_1_position);
      // glLightfv(GL_LIGHT2,   GL_POSITION, light_2_position);

      // glEnable(GL_LIGHT0);
      // glEnable(GL_LIGHT1);

      glEnable(GL_LIGHTING);
      glEnable(GL_DEPTH_TEST);

      glPopMatrix();
   } else {
      glDisable(GL_LIGHTING);
      glDisable(GL_LIGHT0);
      glDisable(GL_LIGHT1);
      glDisable(GL_LIGHT2);
   }
}

void show_lighting() {

   // Is this a useful function?

   std::cout << "------------------------------------------------------------------------------------------------------------------"
             << " show_lighting() -------------------------------------------------------------------" << std::endl;

   if (true) {

      glEnable (GL_BLEND); // these 2 lines are needed to make the transparency work.
      glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

      glPushMatrix();
      glLoadIdentity();

      GLfloat  light_0_position[] = { 1.0,  1.0, 1.0, 0.0};
      GLfloat  light_1_position[] = { 1.0, -1.0, 1.0, 0.0};

      glLightfv(GL_LIGHT0, GL_POSITION, light_0_position);
      glLightfv(GL_LIGHT1, GL_POSITION, light_1_position);
      glPopMatrix();
   }
}


/* When glarea widget size changes, viewport size is set to match the
   new size */

gint reshape(GtkWidget *widget, GdkEventConfigure *event) { return 0; }

/* When widget is exposed it's contents are redrawn. */
#define DEG_TO_RAD .01745327
gint expose(GtkWidget *widget, GdkEventExpose *event) { return 0; }


void myWireCube(float size) {

	 float corners[8][3] = {
	    {-0.5,-0.5,-0.5}, //0
	    {-0.5,-0.5,0.5}, //1
	    {-0.5,0.5,-0.5}, //2
	    {-0.5,0.5,0.5}, //3
	    {0.5,-0.5,-0.5}, //4
	    {0.5,-0.5,0.5}, //5
	    {0.5,0.5,-0.5}, //6
	    {0.5,0.5,0.5}};//7

	 glBegin(GL_LINES);

	 // bottom left connections
	   glVertex3f(corners[0][0], corners[0][1], corners[0][2]);
	   glVertex3f(corners[1][0], corners[1][1], corners[1][2]);

	   glVertex3f(corners[0][0], corners[0][1], corners[0][2]);
	   glVertex3f(corners[2][0], corners[2][1], corners[2][2]);

	   glVertex3f(corners[0][0], corners[0][1], corners[0][2]);
	   glVertex3f(corners[4][0], corners[4][1], corners[4][2]);

	 // top right front connections
	   glVertex3f(corners[6][0], corners[6][1], corners[6][2]);
	   glVertex3f(corners[4][0], corners[4][1], corners[4][2]);

	   glVertex3f(corners[6][0], corners[6][1], corners[6][2]);
	   glVertex3f(corners[2][0], corners[2][1], corners[2][2]);

	   glVertex3f(corners[6][0], corners[6][1], corners[6][2]);
	   glVertex3f(corners[7][0], corners[7][1], corners[7][2]);

	 // from 5
	   glVertex3f(corners[5][0], corners[5][1], corners[5][2]);
	   glVertex3f(corners[7][0], corners[7][1], corners[7][2]);

	   glVertex3f(corners[5][0], corners[5][1], corners[5][2]);
	   glVertex3f(corners[4][0], corners[4][1], corners[4][2]);

	   glVertex3f(corners[5][0], corners[5][1], corners[5][2]);
	   glVertex3f(corners[1][0], corners[1][1], corners[1][2]);

	 // from 3
	   glVertex3f(corners[3][0], corners[3][1], corners[3][2]);
	   glVertex3f(corners[1][0], corners[1][1], corners[1][2]);

	   glVertex3f(corners[3][0], corners[3][1], corners[3][2]);
	   glVertex3f(corners[7][0], corners[7][1], corners[7][2]);

	   glVertex3f(corners[3][0], corners[3][1], corners[3][2]);
	   glVertex3f(corners[2][0], corners[2][1], corners[2][2]);

	   glEnd();

}

void
draw_crosshairs_maybe() {

   if (graphics_info_t::draw_crosshairs_flag) {
      glPushMatrix();
      float s = 0.033 * 100/graphics_info_t::zoom;
      //       std::cout << " zoom: " << graphics_info_t::zoom << std::endl;

      glLoadIdentity();
      GLfloat grey[3] = {0.7, 0.7, 0.7};
      glColor3fv(grey);
      glMatrixMode(GL_PROJECTION);
      glPushMatrix();

      glLoadIdentity();

      glBegin(GL_LINES);

      // screen y axis

      float val = 3.8;
      glVertex3f(0.0, -val*s, 0.0);
      glVertex3f(0.0,  val*s, 0.0);

      val = 3.8;
      glVertex3f( 0.0,   -val*s, 0.0);
      glVertex3f(-0.3*s, -val*s, 0.0);
      glVertex3f( 0.0,    val*s, 0.0);
      glVertex3f(-0.3*s,  val*s, 0.0);

      val = 1.54;
      glVertex3f( 0.0,   -val*s, 0.0);
      glVertex3f(-0.3*s, -val*s, 0.0);
      glVertex3f( 0.0,    val*s, 0.0);
      glVertex3f(-0.3*s,  val*s, 0.0);

      val = 2.7;
      glVertex3f( 0.0,   -val*s, 0.0);
      glVertex3f(-0.3*s, -val*s, 0.0);
      glVertex3f( 0.0,    val*s, 0.0);
      glVertex3f(-0.3*s,  val*s, 0.0);


      // screen x axis
      //
      // adjust for the width being strange
      GtkAllocation allocation;
      gtk_widget_get_allocation(graphics_info_t::glareas[0], &allocation);
      float adjustment = float(allocation.height) / float(allocation.width);
      s *= adjustment;

      val = 3.8;
      glVertex3f(-val*s, 0.0, 0.0);
      glVertex3f( val*s, 0.0, 0.0);

      val = 3.8;
      glVertex3f(-val*s,  0.0,   0.0);
      glVertex3f(-val*s, -0.3*s, 0.0);
      glVertex3f( val*s,  0.0,   0.0);
      glVertex3f( val*s, -0.3*s, 0.0);

      val = 2.7;
      glVertex3f(-val*s,  0.0,   0.0);
      glVertex3f(-val*s, -0.3*s, 0.0);
      glVertex3f( val*s,  0.0,   0.0);
      glVertex3f( val*s, -0.3*s, 0.0);

      val = 1.54;
      glVertex3f(-val*s,  0.0,   0.0);
      glVertex3f(-val*s, -0.3*s, 0.0);
      glVertex3f( val*s,  0.0,   0.0);
      glVertex3f( val*s, -0.3*s, 0.0);

      glEnd();
      glPopMatrix();
      glMatrixMode(GL_MODELVIEW);
      glPopMatrix();
   }
}

void
debug_draw_rotation_axes(float y_x, float y_z, float x_y, float x_z) {

   glBegin(GL_LINES);
   // x
   glVertex3f(0.0, 0.0, 0.0);
   glVertex3f(y_x*10, 0.0, y_z*10);
   // y
   glVertex3f(0.0, 0.0, 0.0);
   glVertex3f(x_y*10, 0.0, x_z*10);
   // z
   glVertex3f(0.0, 0.0, 0.0);
   glVertex3f(0.0, 1.0*10, 0.0);
   //
   glEnd();

   // glRasterPos3f();
   graphics_info_t::printString("x", 12, 0.0, 0.0);
   // glRasterPos3f();
   graphics_info_t::printString("y", 0, 12.0, 0.0);
   // glRasterPos3f();
   graphics_info_t::printString("z", 0, 0.0, 12.0);
}

gint glarea_motion_notify (GtkWidget *widget, GdkEventMotion *event) {

   graphics_info_t info;
   gdouble x;
   gdouble y;
   int x_as_int=0, y_as_int=0;

   gdouble x_diff, y_diff;

   //GdkRectangle area; // needed for redrawing a fraction of the area
//    area.x = 0;
//    area.y = 0;
//    area.width  = widget->allocation.width;
//    area.height = widget->allocation.height;

   GtkAllocation allocation;
   gint width;
   gint height;
   gtk_widget_get_allocation (widget, &allocation);
   width = allocation.width;
   height = allocation.height;

   GdkModifierType my_button1_mask = info.gdk_button1_mask();
   GdkModifierType my_button2_mask = info.gdk_button2_mask();
   GdkModifierType my_button3_mask = info.gdk_button3_mask();

   GdkModifierType state = GDK_SHIFT_MASK; // 20211124-PE give it an initial value, stops the compiler complaining
   short int button_was_pressed = 0;

   if (event->is_hint) {
      // gdk_window_get_pointer(event->window, &x_as_int, &y_as_int, &state);
      GdkModifierType mask;
      GdkSeat *seat = gdk_display_get_default_seat(gdk_display_get_default());
      GdkDevice *mouse = gdk_seat_get_pointer(seat);
      gdk_window_get_device_position(event->window, mouse, &x_as_int, &y_as_int, &mask);
      x = x_as_int;
      y = y_as_int;
   } else {
      // when does this happen?  Never as far as I can see.
      // Actually Bernie says that it happens for him on WindowsXX
      //
      // std::cout << "!!!!!!!!! changing state!" << std::endl;
      x_as_int = int(event->x);
      y_as_int = int(event->y);
      x = event->x;
      y = event->y;
      state = (GdkModifierType) event->state;
   }

   // Try to correct cntrl and shift anomalies:
   //
   if (event->state & GDK_CONTROL_MASK)
      info.control_is_pressed = 1;
   else
      info.control_is_pressed = 0;
   if (event->state & GDK_SHIFT_MASK)
      info.shift_is_pressed = 1;
   else
      info.shift_is_pressed = 0;

   if (state & my_button1_mask) {

      short int local_rotate_view_mode = 0;

      button_was_pressed = 1;
      if (info.baton_mode) {

	 info.rotate_baton(x,y);

      } else {

	 if (info.in_edit_chi_mode_flag ||
	     info.in_edit_torsion_general_flag ||
	     info.in_multi_residue_torsion_mode) {

	    if (info.in_edit_chi_mode_flag) {

	       if (info.in_edit_chi_mode_view_rotate_mode) {
		  local_rotate_view_mode = 1;
	       } else {
		  if (graphics_info_t::do_probe_dots_on_rotamers_and_chis_flag) {
		     graphics_info_t g;
		     g.do_probe_dots_on_rotamers_and_chis();
		  }
		  info.rotate_chi(x, y);
	       }
	    }
	    if (info.in_edit_torsion_general_flag) {
	       if (info.in_edit_chi_mode_view_rotate_mode) {
		  local_rotate_view_mode = 1;
	       } else {
		  info.rotate_chi_torsion_general(x,y);
	       }
	    }
	    if (info.in_multi_residue_torsion_mode) {
	       if (info.in_edit_chi_mode_view_rotate_mode) {
		  local_rotate_view_mode = 1;
	       } else {
		  info.rotate_multi_residue_torsion(x,y);
	       }
	    }

	 } else {

	    // not rotating chi or torsion_general
	    if (info.in_moving_atoms_drag_atom_mode_flag) {

	       if (info.control_is_pressed) {

		  // single atom move, except if we are in
		  // rotate/translate zone, where a Ctrl-click does a
		  // intermediate-atom-screen-z-rotate not a single
		  // atom move.

		  info.move_single_atom_of_moving_atoms(x_as_int, y_as_int);
		  // info.move_atom_pull_target_position(x_as_int ,y_as_int);

	       } else {

		  // multi atom move

#ifdef HAVE_GSL
		  if (info.last_restraints_size() > 0) {

                    info.move_atom_pull_target_position(x_as_int, y_as_int);

		  } else {
		     // don't allow translation drag of the
		     // intermediate atoms when they are a rotamer:
		     //
		     if (! info.rotamer_dialog) {
			// e.g. translate an added peptide fragment.
			info.move_moving_atoms_by_simple_translation(x_as_int,
								     y_as_int);
		     }
		  }
#endif // HAVE_GSL
	       }

	    } else {
	       // info.in_moving_atoms_drag_atom_mode_flag test

	       bool handled_non_atom_drag_event = false;
           handled_non_atom_drag_event = info.rotate_intermediate_atoms_maybe(width, height);

	       if (! handled_non_atom_drag_event) {

		  bool do_drag_pan_flag = false;
		  if (info.static_graphics_pick_pending()) {
		     if (info.control_key_for_rotate_flag) {
			if (event->state & GDK_CONTROL_MASK) {
			   local_rotate_view_mode = 1;
			} else {
			   //
			}
		     } else {
			// ctrl is needed for picking, not mouse motion
			local_rotate_view_mode = 1;
		     }
		  } else {
		     if (event->state & GDK_CONTROL_MASK) {
			do_drag_pan_flag = true;
		     } else {
			local_rotate_view_mode = 1;
		     }
		  }

		  if (do_drag_pan_flag) {
		     do_drag_pan(x, y, widget);
		  } // control_is_pressed
	       }
	    }
	 }

      } // baton_mode

      // This is factored out and put here because there are 2 ways of
      // getting here now: either by conventional left mouse drag or
      // by cntl left mouse drag when in edit_chi_mode.
      //
      // Perhaps it should have its own function.

      if (local_rotate_view_mode) {

	 /* Mouse button 1 is engaged (and we are in motion) */

	 // save it for the quaternion construction
	 //
	 info.mouse_current_x = x;
	 info.mouse_current_y = y;

	 // 	 x_diff = x - info.GetMouseBeginX();
	 // 	 y_diff = y - info.GetMouseBeginY();

	 // Stop the nasty orientation glitch that we get when
	 // dragging the mouse from one gl context to the other in
	 // side-by-side stereo.
	 //
	 if (fabs(info.GetMouseBeginX()-info.mouse_current_x) < 300) {

	    float spin_quat[4];

	    // map the coordinates of the mouse to the range -1:1, (so
	    // top right is at (1,1) and the centre of the window is at
	    // (0,0).
	    //
	    // modify spin_quat:
	    GtkAllocation allocation;
	    gtk_widget_get_allocation(widget, &allocation);
	    trackball(spin_quat,
		      (2.0*info.GetMouseBeginX() - allocation.width)/allocation.width,
		      (allocation.height - 2.0*info.GetMouseBeginY())/allocation.height,
		      (2.0*info.mouse_current_x - allocation.width)  /allocation.width,
		      (allocation.height -  2.0*info.mouse_current_y)/allocation.height,
		      info.get_trackball_size() );

	    // 	 cout << (2.0*info.GetMouseBeginX() - widget->allocation.width) /widget->allocation.width
	    // 	      << " "
	    // 	      << (2.0*info.mouse_current_x - widget->allocation.width)  /widget->allocation.width
	    // 	      <<  "      "
	    // 	      << (widget->allocation.height - 2.0*info.GetMouseBeginY())/widget->allocation.height
	    // 	      << " "
	    // 	      << (widget->allocation.height -  2.0*info.mouse_current_y)/widget->allocation.height
	    // 	      << endl;

	    // args: q1, q2, destination
	    // add_quats(spin_quat, info.quat, info.quat);

	    // 	 cout << "spin_quat " << spin_quat[0] << " " << spin_quat[1] << " "
	    // 	      << spin_quat[2] << " " << spin_quat[3] << endl;
	    // 	 cout << "info.quat " << info.quat[0] << " " << info.quat[1] << " "
	    // 	      << info.quat[2] << " " << info.quat[3] << endl;


	    // 	 mol_rot.x_axis_angle += y_diff*0.14;  // we don't like x
	    //                                                // rotation much.
	    // 	 mol_rot.y_axis_angle += x_diff*0.3;

	    // info.SetMouseBegin(x,y);  // for next round

	    // redraw because the scene rotation centre has changed.
	    //
	    info.graphics_draw();
	 }

      }

   } // my_button1_mask

   if (state & my_button3_mask) {

      button_was_pressed = 1;
      /* Mouse button 3 is engaged */
//       cout << "Button 3 motion: Zoom related "
//        << event->x << " " << event->y << endl;

      if (fabs(info.GetMouseBeginX()-x) < 300) {
	 if (info.control_is_pressed) {
	    if (info.shift_is_pressed) {
	       do_screen_z_rotate(x, y); // z-rotation
	    } else {
	       // Pymol-like z-translation or slab adjustment
	       // requested by FvD 20050304
	       //
	       // Diagonal up-right/down-left changes the clipping
	       //
	       // Diagonal up-left/down-righ changes the z-translation
	       // (c.f keypad 3 and .: keypad_translate_xyz(3, 1)).
	       //
	       do_ztrans_and_clip(x, y);
	    }
	 } else {
	    do_button_zoom(x, y);
	 }
      }
   }

   if (state & my_button2_mask) {
      button_was_pressed = 1;
      /* Mouse button 2 is engaged */
      // cout << "Button 2 motion " << event->x << " " << event->y;

      do_drag_pan(x, y, widget);
   }

   if (! button_was_pressed) {
      if (info.quanta_like_zoom_flag) {
	 if (info.shift_is_pressed) {
	    do_button_zoom(x,y);
	 }
      }
   }

   // because we are only really interested in how the mouse has
   // moved from the previous position, not from when the mouse was
   // clicked.
   //
   info.SetMouseBegin(x,y);

   return TRUE;

}

void
do_drag_pan(gdouble x, gdouble y, GtkWidget *widget) {

   graphics_info_t info;

   // ----- DRAG -----

   // We are in (density) drag mode:
   //
   // We need to do two things:
   // i) change the rotation axis centre.
   //    Question: which direction are we dragging in...?
   //    That is, what is the vector in coordinate space
   //    that corresponds to the screen x axis?
   //
   // ii) update the map about the new rotation centre
   //

   gdouble x_diff = x - info.GetMouseBeginX();
   gdouble y_diff = y - info.GetMouseBeginY();

   coot::CartesianPair vec_x_y = screen_x_to_real_space_vector(widget);

   // x_diff and y_diff are the scale factors to the x and y
   // drag vectors.
   //
   info.add_to_RotationCentre(vec_x_y, -x_diff*0.02, -y_diff*0.02);

//    if (info.GetActiveMapDrag() == 1) {
//       for (int ii=0; ii<info.n_molecules(); ii++) {
// 	 info.molecules[ii].update_map(); // to take account
// 	 // of new rotation centre.
//       }
//    }

   info.update_maps();

   for (int ii=0; ii<info.n_molecules(); ii++) {
      info.molecules[ii].update_symmetry();
   }
   info.graphics_draw();

}

void
do_button_zoom(gdouble x, gdouble y) {
      /* Mouse button 3 is engaged */
      //cout << "Button 3 motion: Zoom related "
      // << event->x << " " << event->y << endl;

   graphics_info_t info;

   gdouble x_diff = x - info.GetMouseBeginX();
   scale_zoom_internal(1 -  x_diff/200.0);

   gdouble y_diff = y - info.GetMouseBeginY();
   scale_zoom_internal(1 -  y_diff/200.0);

   //cout << "difference on zoom: " << x << " - "
   //	   << info.GetMouseBeginX() << " = " << x_diff << endl;
   //cout << "Zoom : " << info.zoom << endl;

   if (info.dynamic_map_resampling == 1) {

      // is a new new resampling triggered by this new zoom?
      //
      // int iv = 1 + int (0.009*info.zoom) + info.dynamic_map_zoom_offset;
      int iv = 1 + int (0.009*(info.zoom + info.dynamic_map_zoom_offset));
      if (iv != info.graphics_sample_step) {

         for (int imap=0; imap<info.n_molecules(); imap++) {
            info.molecules[imap].update_map(true); // uses g.zoom
         }
         info.graphics_sample_step = iv;
      }
   }
   // redraw it
   info.graphics_draw();
}

void
do_screen_z_rotate(gdouble x, gdouble y) {

#if 0
   // c. f. animate_idle()

   float spin_quat[4];
   graphics_info_t g;


   // int x0 = graphics_info_t::glarea->allocation.width;
   // int y0 = graphics_info_t::glarea->allocation.height;

   // gdouble x_diff = x - g.GetMouseBeginX();
   // gdouble y_diff = y - g.GetMouseBeginY();
   // double y2 = 1.0 + (x_diff + y_diff)*0.01;

//    double a_x = x - x0;
//    double a_y = y - y0;
//    double b_x = g.GetMouseBeginX() - x0;
//    double b_y = g.GetMouseBeginY() - y0;
   // double alen = sqrt(a_x*a_x + a_y*a_y);
   // double blen = sqrt(b_x*b_x + b_y*b_y);
   // double y2 = 0.001 * (a_x*b_x + a_y*b_y)/(alen*blen);

   gdouble x_diff = 0.01 * (x - g.GetMouseBeginX());
   gdouble y_diff = 0.01 * (y - g.GetMouseBeginY());
   x_diff += y_diff;

//    trackball(spin_quat,
// 	     1, 1.0,
// 	     y2, ,
// 	     9.0);
//    trackball(spin_quat,
// 	     (2.0*g.GetMouseBeginX() - x0)                    /x0,
// 	     (y0                     - 2.0*g.GetMouseBeginY())/y0,
// 	     (2.0*x                  - x0)                    /x0,
// 	     (y0                     - 2.0*y)                 /y0,
// 	     1.0);
   trackball(spin_quat,  0.0, 1.0, x_diff, 1.0, 0.5);
   add_quats(spin_quat, g.quat, g.quat);
   g.graphics_draw();
#endif
}

void
do_ztrans_and_clip(gdouble x, gdouble y) {

   graphics_info_t g;
   gdouble x_diff = x - g.GetMouseBeginX();
   gdouble y_diff = y - g.GetMouseBeginY();

   coot::Cartesian v = screen_z_to_real_space_vector(graphics_info_t::glareas[0]);

//    // Like Frank showed me in Pymol?
//    double slab_change = 0.05  * (x_diff + y_diff); // about 0.2?
//    double ztr_change  = 0.005 * (x_diff - y_diff);

   // I prefer straight left/right and up and down:
   // (Up and down move the individual clipping planes in PyMol.
   // I'm not keen on allowing that in Coot).
   double slab_change = -0.02 * x_diff; // about 0.2?
   double ztr_change  = 0.001 * y_diff;

   // Do *either* z-trans or clipping, not both.  Doing both confuses
   // me.

   v *= ztr_change;
   if (fabs(x_diff) > fabs(y_diff))
      adjust_clipping(slab_change);
   else
      g.add_vector_to_RotationCentre(v);

   g.graphics_draw();

}


// move this to graphics_info_t and its call from key-bindings.
void
adjust_clipping(double d) {

   graphics_info_t g;
   g.adjust_clipping(d);
}

#include "c-interface-ligands-swig.hh"

gint key_press_event(GtkWidget *widget, GdkEventKey *event)
{

   // Try to correct cntrl and shift anomalies:
   // (I don't think this code does anything useful...)
   //
   if (event->state & GDK_CONTROL_MASK) {
      graphics_info_t::control_is_pressed = 1;
   } else {
      graphics_info_t::control_is_pressed = 0;
   }
   if (event->state & GDK_SHIFT_MASK) {
      graphics_info_t::shift_is_pressed = 1;
   } else {
      graphics_info_t::shift_is_pressed = 0;
   }

   gint handled = 0; // initially unhandled

   // std::cout << "keyval: " << event->keyval << " " << std::hex << event->keyval << std::endl;

   switch (event->keyval) {
   case GDK_KEY_Control_L:
   case GDK_KEY_Control_R:

//       std::cout << "DEBUG ctrl key press  : graphics_info_t::control_is_pressed "
// 		<< graphics_info_t::control_is_pressed
// 		<< " graphics_info_t::pick_pending_flag "
// 		<< graphics_info_t::pick_pending_flag << std::endl;

      graphics_info_t::control_is_pressed = 1; // TRUE.
      if (graphics_info_t::in_edit_chi_mode_flag)
	 graphics_info_t::in_edit_chi_mode_view_rotate_mode = 1;
      if (graphics_info_t::in_edit_torsion_general_flag)
	 graphics_info_t::in_edit_chi_mode_view_rotate_mode = 1;
      if (graphics_info_t::in_multi_residue_torsion_mode)
	 graphics_info_t::in_edit_chi_mode_view_rotate_mode = 1;

      if (graphics_info_t::control_key_for_rotate_flag) {
	 normal_cursor();
      } else {
	 // control key is for pick
	 if (graphics_info_t::pick_pending_flag) {
	    pick_cursor_maybe();
	 }
      }
      handled = TRUE;
      break;

   case GDK_KEY_Alt_L:
   case GDK_KEY_Alt_R:
   case GDK_KEY_Meta_L:
   case GDK_KEY_Meta_R:

      handled = TRUE; // stops ALT key getting through to key-press-hook
      break;

   case GDK_KEY_Shift_L: // stops Shift key getting through to key-press-hook
   case GDK_KEY_Shift_R:
      graphics_info_t::shift_is_pressed = 1;
      handled = TRUE;
      break;

   case GDK_KEY_Return:

      if (graphics_info_t::accept_reject_dialog) {

	 if (graphics_info_t::continue_threaded_refinement_loop) {
	    graphics_info_t::continue_threaded_refinement_loop = false;
	    std::cout << ".... Return key tells refinement to accept and clean up" << std::endl;
	    graphics_info_t::threaded_refinement_needs_to_clear_up = true;
	    graphics_info_t::threaded_refinement_needs_to_accept_moving_atoms = true;
	 } else {
	    accept_regularizement();

	    if (graphics_info_t::accept_reject_dialog_docked_flag == coot::DIALOG) {
	       save_accept_reject_dialog_window_position(graphics_info_t::accept_reject_dialog);
	       gtk_widget_destroy(graphics_info_t::accept_reject_dialog);
	    } else {
	       // have docked dialog
	       if (graphics_info_t::accept_reject_dialog_docked_show_flag == coot::DIALOG_DOCKED_HIDE) {
	          gtk_widget_hide(graphics_info_t::accept_reject_dialog);
	       } else {
	          gtk_widget_set_sensitive(graphics_info_t::accept_reject_dialog, FALSE);
	       }
	    }
	    graphics_info_t::accept_reject_dialog = 0;
	 }
      }

      if (graphics_info_t::rotamer_dialog) {
	 accept_regularizement();
	 gtk_widget_destroy(graphics_info_t::rotamer_dialog);
	 set_graphics_rotamer_dialog(NULL);
      }

      handled = TRUE;
      break;

   case GDK_KEY_Escape:

      graphics_info_t::rebond_molecule_corresponding_to_moving_atoms();

      // poke a value into the threaded refinement loop, to stop
      if (graphics_info_t::continue_threaded_refinement_loop) {
	 // and tell it to clear up the moving atoms
	 graphics_info_t::threaded_refinement_needs_to_clear_up = true;
	 std::cout << ".... Esc key tells refinement to clean up" << std::endl;
	 graphics_info_t::continue_threaded_refinement_loop = false;
      } else {

	 // refinement was not running. we can clear up the atoms ourselves
	 graphics_info_t g;
	 g.clear_up_moving_atoms();
	 g.clear_moving_atoms_object();

	 if (graphics_info_t::accept_reject_dialog) {
	    if (graphics_info_t::accept_reject_dialog_docked_flag == coot::DIALOG) {
	       save_accept_reject_dialog_window_position(graphics_info_t::accept_reject_dialog);
	       // this calls clear_up_moving_atoms() and clears atom pull restraint.
	       gtk_widget_destroy(graphics_info_t::accept_reject_dialog);
	       graphics_info_t::accept_reject_dialog = 0;
	    } else {
	       gtk_widget_set_sensitive(graphics_info_t::accept_reject_dialog, FALSE);
	    }
	 }
      }

      handled = TRUE;
      break;

   case GDK_KEY_space:
      handled = TRUE;
      break; // stops Space key getting through to key-press-hook

   case GDK_KEY_g:
      // say I want to go to residue 1G: first time Ctrl-G (second if)
      // and then the first if.
      if (graphics_info_t::control_is_pressed) {
	 show_go_to_residue_keyboarding_mode_window();
	 handled = TRUE;
      }
      break;

   case GDK_KEY_i:
      // throw away i key pressed (we act on i key released).
      handled = TRUE;
      break;

   case GDK_KEY_a:

      if (graphics_info_t::in_range_define_for_refine == 2) {

	 graphics_info_t::a_is_pressed = 1;

	 int rot_trans_rotation_origin_atom = 0; // flag for Ctrl left
	 // mouse behaviour (we don't want to rotate
	 // the atoms)
	 graphics_info_t::watch_cursor();
	 int auto_range_flag = 1;
	 graphics_info_t g;
	 g.refine(graphics_info_t::residue_range_mol_no,
		  auto_range_flag,
		  graphics_info_t::residue_range_atom_index_1,
		  graphics_info_t::residue_range_atom_index_1);
	 g.in_range_define_for_refine = 0;
	 normal_cursor();
	 g.pick_pending_flag = 0;
	 g.model_fit_refine_unactive_togglebutton("model_refine_dialog_refine_togglebutton");
	 handled = TRUE;
      }
      break;

   case GDK_KEY_b:
      break;

   case GDK_KEY_c:
      break; // stops C key getting through to key-press-hook

   case GDK_KEY_u:
      undo_last_move();
      handled = TRUE;
      break;

   case GDK_KEY_r:
      if (graphics_info_t::control_is_pressed) {
	 toggle_idle_rock_function();
	 handled = TRUE;
      }
      break;


   case GDK_KEY_s:
      if (graphics_info_t::control_is_pressed) {
	 quick_save();
	 handled = TRUE;
      }
      break;

   case GDK_KEY_d:

      if (graphics_info_t::control_is_pressed) {
	 graphics_info_t g;
	 g.delete_active_residue();
      } else {
	 if (graphics_info_t::clipping_back < 15.0) {
	    set_clipping_front(graphics_info_t::clipping_front + 0.4);
	    set_clipping_back (graphics_info_t::clipping_front + 0.4);
	    // std::cout << "INFO:: clipping " << graphics_info_t::clipping_front << " "
	    // << graphics_info_t::clipping_back << std::endl;
	 }
      }
      handled = TRUE;
      break;

   case GDK_KEY_e:
      if (graphics_info_t::control_is_pressed) {

      }
      break;

   case GDK_KEY_E:
      break;

   case GDK_KEY_f:

      if (graphics_info_t::clipping_back > -15.2) {
	 set_clipping_front(graphics_info_t::clipping_front - 0.4);
	 set_clipping_back (graphics_info_t::clipping_front - 0.4);
	 // std::cout << "INFO:: clipping " << graphics_info_t::clipping_front << " "
	 // << graphics_info_t::clipping_back << std::endl;
      }
      handled = TRUE;
      break;

   case GDK_KEY_n:

      for (int i=0; i<5; i++) {
	 graphics_info_t::zoom *= 1.01;
	 graphics_info_t::graphics_draw();
      }
      handled = TRUE;
      break;

   case GDK_KEY_m:

      for (int i=0; i<5; i++) {
	 graphics_info_t::zoom *= 0.99;
	 graphics_info_t::graphics_draw();
      }
      handled = TRUE;
      break;

   case GDK_KEY_1:
   case GDK_KEY_KP_1:
      if (graphics_info_t::control_is_pressed) {
	 std::pair<bool, std::pair<int, coot::atom_spec_t> > aa = active_atom_spec();
	 if (aa.first) {
            coot_all_atom_contact_dots(aa.second.first);
         }
         handled = true;
         break;
      } else {
         if (graphics_info_t::moving_atoms_move_chis_flag) {
            graphics_info_t g;
            g.setup_flash_bond_using_moving_atom_internal(0);
            graphics_info_t::edit_chi_current_chi = 1;
            graphics_info_t::in_edit_chi_mode_flag = 1; // on
         }
         handled = TRUE;
         break;
      }
   case GDK_KEY_2:
   case GDK_KEY_KP_2:
      if (graphics_info_t::moving_atoms_move_chis_flag) {
         graphics_info_t g;
         g.setup_flash_bond_using_moving_atom_internal(1);
	 graphics_info_t::edit_chi_current_chi = 2;
	 graphics_info_t::in_edit_chi_mode_flag = 1; // on
      }
      handled = TRUE;
      break;
   case GDK_KEY_3:
   case GDK_KEY_KP_3:
      if (graphics_info_t::moving_atoms_move_chis_flag) {
         graphics_info_t g;
         g.setup_flash_bond_using_moving_atom_internal(2);
	 graphics_info_t::edit_chi_current_chi = 3;
	 graphics_info_t::in_edit_chi_mode_flag = 1; // on
      } else {
	 keypad_translate_xyz(3, 1);
      }
      handled = TRUE;
      break;
   case GDK_KEY_4:
   case GDK_KEY_KP_4:
      if (graphics_info_t::moving_atoms_move_chis_flag) {
         graphics_info_t g;
         g.setup_flash_bond_using_moving_atom_internal(3);
	 graphics_info_t::edit_chi_current_chi = 4;
	 graphics_info_t::in_edit_chi_mode_flag = 1; // on
      }
      handled = TRUE;
      break;
   case GDK_KEY_5:
   case GDK_KEY_KP_5:
      if (graphics_info_t::moving_atoms_move_chis_flag) {
         graphics_info_t g;
         g.setup_flash_bond_using_moving_atom_internal(4);
	 graphics_info_t::edit_chi_current_chi = 5;
	 graphics_info_t::in_edit_chi_mode_flag = 1; // on
      }
      handled = TRUE;
      break;
   case GDK_KEY_6:
   case GDK_KEY_KP_6:
      if (graphics_info_t::moving_atoms_move_chis_flag) {
         graphics_info_t g;
         g.setup_flash_bond_using_moving_atom_internal(5);
	 graphics_info_t::edit_chi_current_chi = 6;
	 graphics_info_t::in_edit_chi_mode_flag = 1; // on
      }
      handled = TRUE;
      break;

   case GDK_KEY_7:
   case GDK_KEY_KP_7:
      handled = TRUE;
      break;

   case GDK_KEY_8:
   case GDK_KEY_KP_8:
      handled = TRUE;
      break;

   case GDK_KEY_9:
   case GDK_KEY_KP_9:
      handled = TRUE;
      break;

   case GDK_KEY_0:
   case GDK_KEY_KP_0:
      graphics_info_t::edit_chi_current_chi = 0;
      graphics_info_t::in_edit_chi_mode_flag = 0; // off
      handled = TRUE;
      break;

   case GDK_KEY_l:
      {
	 graphics_info_t g;

	 if (graphics_info_t::control_is_pressed) {
	    go_to_ligand();
	 } else {

	    std::pair<int, int> cl_at = g.get_closest_atom();
	    if (is_valid_model_molecule(cl_at.second)) {
	       g.molecules[cl_at.second].add_to_labelled_atom_list(cl_at.first);
	       // shall add to status_bar ? maybe this should be a function?
	       // it is now
	       std::string at_info = atom_info_as_text_for_statusbar(cl_at.first, cl_at.second);
	       g.add_status_bar_text(at_info);
	       mmdb::Atom *at = g.molecules[cl_at.second].atom_sel.atom_selection[cl_at.first];
	       mmdb::Residue *residue_p  = at->residue;
	       std::string alt_conf = at->altLoc;
	       g.setup_graphics_ligand_view(cl_at.second, residue_p, alt_conf);
	    }
	    graphics_info_t::graphics_draw();
	 }
      }
      handled = TRUE;
      break;

   case GDK_KEY_x: // delete active residue

      if (graphics_info_t::control_is_pressed) {
	 std::pair<bool, std::pair<int, coot::atom_spec_t> > aa = active_atom_spec();
	 if (aa.first) {
	    delete_residue(aa.second.first, aa.second.second.chain_id.c_str(),
			   aa.second.second.res_no,
			   aa.second.second.ins_code.c_str());
	 }
	 handled = TRUE;
      }
      break;

   case GDK_KEY_z:
      graphics_info_t::z_is_pressed = 1;
      if (graphics_info_t::control_is_pressed) {
	 if (graphics_info_t::draw_baton_flag)
	    baton_build_delete_last_residue();
	 else
	    apply_undo();
	 handled = TRUE;
      }
      break;

   case GDK_KEY_y:
      graphics_info_t::y_is_pressed = 1;
      if (graphics_info_t::control_is_pressed) {
	 apply_redo();
	 handled = TRUE;
      }
      break;

   case GDK_KEY_F5:
      post_model_fit_refine_dialog();
      handled = TRUE;
      break;
   case GDK_KEY_F6:
      post_go_to_atom_window();
      handled = TRUE;
      break;
   case GDK_KEY_F7:
      post_display_control_window();
      handled = TRUE;
      break;
   case GDK_KEY_F8:
      // case GDK_3270_PrintScreen: // that the gnome screenshot
      raster_screen_shot();
      handled = TRUE;
      break;
   case GDK_KEY_KP_Down:
   case GDK_KEY_KP_Page_Down:
      keypad_translate_xyz(3, -1);
      handled = TRUE;
      break;
   case GDK_KEY_KP_Up:
   case GDK_KEY_KP_Page_Up:
      keypad_translate_xyz(3, 1);
      handled = TRUE;
      break;
   case GDK_KEY_KP_Decimal:
   case GDK_KEY_KP_Delete:
      keypad_translate_xyz(3, -1);
      handled = TRUE;
      break;

   case GDK_KEY_period:
      // std::cout << "Got a .\n";
      if (graphics_info_t::rotamer_dialog) {
	 // We have to make a synthetic keypress on the "next" rotamer.
	 // So, what is the next rotamer (bear in mind wrapping).
	 graphics_info_t::rotamer_dialog_next_rotamer();
      } else {
	 if (graphics_info_t::difference_map_peaks_dialog) {
	    graphics_info_t::difference_map_peaks_next_peak();
	 } else {
	    if (graphics_info_t::checked_waters_baddies_dialog) {
	       graphics_info_t::checked_waters_next_baddie(+1);
	    } else {
#ifdef USE_GUILE
	       std::string scheme_command("(graphics-dot-key-pressed-hook)");
	       safe_scheme_command(scheme_command);
#endif // USE_GUILE
#ifdef USE_PYTHON
	       std::string python_command("graphics_dot_key_pressed_hook()");
	       safe_python_command(python_command);
#endif // PYTHON
	    }
	 }
      }
      handled = TRUE;
      break;

   case GDK_KEY_comma:
      //       std::cout << "Got an ,\n";
      if (graphics_info_t::rotamer_dialog) {
	 graphics_info_t::rotamer_dialog_previous_rotamer();
      } else {
	 if (graphics_info_t::difference_map_peaks_dialog) {
	    graphics_info_t::difference_map_peaks_previous_peak();
	 } else {
	    if (graphics_info_t::checked_waters_baddies_dialog) {
	       graphics_info_t::checked_waters_next_baddie(-1); // prev
	    } else {
#ifdef USE_GUILE
	       std::string scheme_command("(graphics-comma-key-pressed-hook)");
	       safe_scheme_command(scheme_command);
#endif // USE_GUILE
#ifdef USE_PYTHON
	       std::string python_command("graphics_comma_key_pressed_hook()");
	       safe_python_command(python_command);
#endif // USE_PYTHON
	    }
	 }
      }
      handled = TRUE;
      break;

   case GDK_KEY_minus:
      handled = TRUE;
      break;
   case GDK_KEY_plus:
      handled = TRUE;
      break;
   case GDK_KEY_equal:
      handled = TRUE;
      break;

   case GDK_KEY_Left:
      if (graphics_info_t::control_is_pressed) {
         if (graphics_info_t::shift_is_pressed)
            graphics_info_t::nudge_active_residue_by_rotate(GDK_KEY_Left);
         else
            graphics_info_t::nudge_active_residue(GDK_KEY_Left);
      } else {
         keypad_translate_xyz(1, 1);
      }
      handled = TRUE;
      break;
   case GDK_KEY_Right:
      if (graphics_info_t::control_is_pressed) {
         if (graphics_info_t::shift_is_pressed)
            graphics_info_t::nudge_active_residue_by_rotate(GDK_KEY_Right);
         else
            graphics_info_t::nudge_active_residue(GDK_KEY_Right);
      } else {
         keypad_translate_xyz(1, -1);
      }
      handled = TRUE;
      break;
   case GDK_KEY_Up:
      if (graphics_info_t::control_is_pressed) {
	 if (graphics_info_t::shift_is_pressed)
	    graphics_info_t::nudge_active_residue_by_rotate(GDK_KEY_Up);
	 else
	    graphics_info_t::nudge_active_residue(GDK_KEY_Up);
      } else {
         keypad_translate_xyz(2, 1);
      }
      handled = TRUE;
      break;
   case GDK_KEY_Down:
      if (graphics_info_t::control_is_pressed) {
	 if (graphics_info_t::shift_is_pressed)
	    graphics_info_t::nudge_active_residue_by_rotate(GDK_KEY_Down);
	 else
	    graphics_info_t::nudge_active_residue(GDK_KEY_Down);
      } else {
         keypad_translate_xyz(2, -1);
      }

      handled = TRUE;
      break;

   }

   // Now to test event->keyval against dynamic "reprogramable" key
   // bindings.

   if (handled == 0) { // initial value
      if (int(event->keyval) == graphics_info_t::ncs_next_chain_skip_key) {
	 if (graphics_info_t::prefer_python) {
#if defined USE_PYTHON
	    std::string python_command("skip_to_next_ncs_chain('forward')");
	    safe_python_command(python_command);
#endif // USE_PYTHON
	 } else {
#if defined USE_GUILE
	    std::string scheme_command("(skip-to-next-ncs-chain 'forward)");
	    safe_scheme_command(scheme_command);
#endif // USE_GUILE
	 }
	 handled = TRUE;
      }

      if (int(event->keyval) == graphics_info_t::ncs_prev_chain_skip_key) {
	 if (graphics_info_t::prefer_python) {
#if defined USE_PYTHON
	    std::string python_command("skip_to_next_ncs_chain('backward')");
	    safe_python_command(python_command);
#endif // USE_PYTHON
	 } else {
#if defined USE_GUILE
	    std::string scheme_command("(skip-to-next-ncs-chain 'backward)");
	    safe_scheme_command(scheme_command);
#endif // USE_GUILE
	 }
	 handled = TRUE;
      }


      if (int(event->keyval) == graphics_info_t::update_go_to_atom_from_current_residue_key) {
	 update_go_to_atom_from_current_position();
	 handled = TRUE;
      }

      if (handled == 0) {
         if (event->keyval != GDK_KEY_backslash) {
	         int ikey = event->keyval;

	    if (graphics_info_t::prefer_python) {
#ifdef USE_PYTHON
#endif
	    } else {
#ifdef USE_GUILE
#endif
	    }

#if defined USE_GUILE && !defined WINDOWS_MINGW
	    std::string scheme_command("(graphics-general-key-press-hook ");
	    // scheme_command += "\"";
	    scheme_command += graphics_info_t::int_to_string(ikey);
	    // scheme_command += "\"";
	    scheme_command += ")";
	    if (false)
	       std::cout << "::::::: key press event not handled~! running scheme command: "
			 << scheme_command << std::endl;
	    safe_scheme_command(scheme_command);
#else // not GUILE
#ifdef USE_PYTHON
	    std::string python_command("graphics_general_key_press_hook(");
	    python_command += graphics_info_t::int_to_string(ikey);
	    python_command += ",";
	    python_command += graphics_info_t::int_to_string(graphics_info_t::control_is_pressed);
	    python_command += ")";
	    safe_python_command(python_command);
#endif // USE_PYTHON
#endif // USE_GUILE

	 } else {
	    std::cout << "INFO:: Ignoring GDK_backslash key press event" << std::endl;
	 }
      }
   }

   return TRUE;
}

// Direction is either +1 or -1 (in or out)
//
void keypad_translate_xyz(short int axis, short int direction) {

  graphics_info_t g;
  if (axis == 3) {
    coot::Cartesian v = screen_z_to_real_space_vector(graphics_info_t::glareas[0]);
    v *= 0.05 * float(direction);
    g.add_vector_to_RotationCentre(v);
  } else {
    gdouble x_diff, y_diff;
    x_diff = y_diff = 0;
    coot::CartesianPair vec_x_y = screen_x_to_real_space_vector(graphics_info_t::glareas[0]);
    if (axis == 1) x_diff = 1;
    if (axis == 2) y_diff = 1;
    g.add_to_RotationCentre(vec_x_y, x_diff * 0.1 * float(direction),
                            y_diff * 0.1 * float(direction));
    if (g.GetActiveMapDrag() == 1) {
      for (int ii=0; ii<g.n_molecules(); ii++) {
        g.molecules[ii].update_map(true); // to take account
        // of new rotation centre.
      }
    }
    for (int ii=0; ii<g.n_molecules(); ii++) {
      g.molecules[ii].update_symmetry();
    }
    g.graphics_draw();

  }
}

#include "idles.hh"

gint key_release_event(GtkWidget *widget, GdkEventKey *event)
{

//    // try to correct cntrl and shift anomalies:
//    //
//    if (event->state & GDK_CONTROL_MASK)
//       graphics_info_t::control_is_pressed = 1;
//    else
//       graphics_info_t::control_is_pressed = 0;
//    if (event->state & GDK_SHIFT_MASK)
//       graphics_info_t::shift_is_pressed = 1;
//    else
//       graphics_info_t::shift_is_pressed = 0;

//    // we are getting key release events as the key is pressed (but not
//    // released) :-).
//    //
//    // std::cout << "key release event" << std::endl;
//    graphics_info_t g;
//    int s = graphics_info_t::scroll_wheel_map;
//    if (s < graphics_info_t::n_molecules()) {
//       if (s >= 0) {
// 	 if (! graphics_info_t::molecules[s].has_xmap()) {  // NXMAP-FIXME
// 	    s = -1; // NO MAP
// 	 }
//       }
//    } else {
//       s = -1; // NO MAP
//    }


//    std::vector<int> num_displayed_maps = g.displayed_map_imols();
//    if (num_displayed_maps.size() == 1)
//       s = num_displayed_maps[0];
//    short int istate = 0;

//    switch (event->keyval) {
//    case GDK_KEY_Control_L:
//    case GDK_KEY_Control_R:

// //       std::cout << "DEBUG ctrl key release: graphics_info_t::control_is_pressed "
// // 		<< graphics_info_t::control_is_pressed
// // 		<< " graphics_info_t::pick_pending_flag "
// // 		<< graphics_info_t::pick_pending_flag << std::endl;

//       graphics_info_t::control_is_pressed = 0; // FALSE.
//       if (graphics_info_t::in_edit_chi_mode_flag)
// 	 graphics_info_t::in_edit_chi_mode_view_rotate_mode = 0;
//       if (graphics_info_t::in_edit_torsion_general_flag)
// 	 graphics_info_t::in_edit_chi_mode_view_rotate_mode = 0;
//       if (graphics_info_t::in_multi_residue_torsion_mode)
// 	 graphics_info_t::in_edit_chi_mode_view_rotate_mode = 0;

//       if (graphics_info_t::control_key_for_rotate_flag) {
// 	 if (graphics_info_t::pick_pending_flag) {
// 	    pick_cursor_maybe();
// 	 } else {
// 	    normal_cursor();
// 	 }
//       } else {
// 	 // control key is for pick
// 	 if (graphics_info_t::pick_pending_flag) {
// 	    normal_cursor();
// 	 } else {
// 	    pick_cursor_maybe();
// 	 }
//       }


//       // When we don't have active map dragging, we want to renew the map at
//       // control release, not anywhere else.
//       //
//       for (int ii=0; ii<graphics_info_t::n_molecules(); ii++) {
// 	 graphics_info_t::molecules[ii].update_map();
// 	 graphics_info_t::molecules[ii].update_clipper_skeleton();
//       }
//       g.make_pointer_distance_objects();
//       g.graphics_draw();
//       break;
//    case GDK_KEY_minus:
//       //
//       // let the object decide which level change it needs (and if it needs it)
//       // using graphics_info_t static members
//       if (s >= 0) {

// 	 std::cout << "here in key_release_event for -" << std::endl;
// 	 graphics_info_t::molecules[s].pending_contour_level_change_count--;
//          std::cout << "graphics_info_t::molecules[s].pending_contour_level_change_count now "
//                    << graphics_info_t::molecules[s].pending_contour_level_change_count << std::endl;
// 	 int contour_idle_token = g_idle_add(idle_contour_function, g.glareas[0]);
// 	 g.set_density_level_string(s, g.molecules[s].contour_level);
// 	 g.display_density_level_this_image = 1;

// 	 // g.graphics_draw();
//       } else {
// 	 std::cout << "WARNING: No map - Can't change contour level.\n";
//       }
//       break;
//    case GDK_KEY_plus:
//    case GDK_KEY_equal:  // unshifted plus, usually.
//       //

//       // let the object decide which level change it needs:
//       //
//       if (s >= 0) {

// 	 graphics_info_t::molecules[s].pending_contour_level_change_count++;
// 	 int contour_idle_token = g_idle_add(idle_contour_function, g.glareas[0]);

// 	 // graphics_info_t::molecules[s].change_contour(1); // positive change
// 	 // graphics_info_t::molecules[s].update_map();

// 	 g.set_density_level_string(s, g.molecules[s].contour_level);
// 	 g.display_density_level_this_image = 1;

// 	 // g.graphics_draw();
//       } else {
// 	 std::cout << "WARNING: No map - Can't change contour level.\n";
//       }
//       break;

//    case GDK_KEY_A:
//    case GDK_KEY_a:
//       graphics_info_t::a_is_pressed = 0;
//       break;

//    case GDK_KEY_B:
//    case GDK_KEY_b:
//       // Only toggle baton mode if we are showing a baton!
//       // (Otherwise confusion reigns!)
//       if (g.draw_baton_flag)
// 	 g.toggle_baton_mode();
//       break;

//    case GDK_KEY_c:
//       if (graphics_info_t::control_is_pressed) {
//          g.copy_active_atom_molecule();
//       } else {
// 	 if (! graphics_info_t::shift_is_pressed) {
// 	    g.draw_crosshairs_flag = 1 - g.draw_crosshairs_flag;
// 	    g.crosshairs_text();
// 	 }
//       }
//       g.graphics_draw();
//       break;

//    case GDK_KEY_e:
//       if (graphics_info_t::control_is_pressed) {
// 	 std::pair<bool, std::pair<int, coot::atom_spec_t> > active_atom = active_atom_spec();
// 	 if (active_atom.first) {
// 	    if (is_valid_model_molecule(active_atom.second.first)) {
// 	       // environment graphics object doesn't work this way.  There is
// 	       // only one of it and its display is controlled by having its bonds box
// 	       // filled or not.
// 	       // graphics_info_t::molecules[imol].toggle_display_environment_graphics_object();
// 	       graphics_draw();
// 	    }
// 	 }
//       }
//       break;

//    case GDK_KEY_s:
//    case GDK_KEY_S:
//       if (graphics_info_t::control_is_pressed) {
// 	 // quick_save() is on button-press
//       } else {
// 	 // as it used to be
// 	 for (int ii = 0; ii< graphics_info_t::n_molecules(); ii++)
// 	    graphics_info_t::molecules[ii].update_clipper_skeleton();
// 	 g.graphics_draw();
//       }
//       break;

//    case GDK_KEY_i:
//    case GDK_KEY_I:
//       if (! graphics_info_t::control_is_pressed) {
// 	 toggle_idle_spin_function();
//       } else {
// 	 std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = active_atom_spec();
// 	 if (pp.first) {
// 	    const coot::atom_spec_t &spec = pp.second.second;
// 	    residue_info_dialog(pp.second.first, spec.chain_id.c_str(), spec.res_no,
// 				spec.ins_code.c_str());
// 	 }
//       }
//       break;

//    case GDK_KEY_l:
//    case GDK_KEY_L:
//       // something here, L is released.
//       break;

//    case GDK_KEY_Shift_L:
//    case GDK_KEY_Shift_R:
//       graphics_info_t::shift_is_pressed = 0;
//       break;

//    case GDK_KEY_z:
//    case GDK_KEY_Z:
//       graphics_info_t::z_is_pressed = 0;
//       break;

//    case GDK_KEY_y:
//    case GDK_KEY_Y:
//       graphics_info_t::y_is_pressed = 0;
//       break;

//    case GDK_KEY_space:
//       // go to next residue
// //       int next = 1;
// //       if (graphics_info_t::shift_is_pressed == 1)
// // 	 next = -1;  // i.e. previous

//       // std::cout << "got a space acting on " << g.go_to_atom_chain()
// // 		<< " " << g.go_to_atom_residue()+ next
// // 		<<  " " << g.go_to_atom_atom_name() << std::endl;

//       // commented 030805:
// //       set_go_to_atom_chain_residue_atom_name(g.go_to_atom_chain(),
// // 					     g.go_to_atom_residue()+next,
// // 					     g.go_to_atom_atom_name());

//       bool reorienting = graphics_info_t::reorienting_next_residue_mode;
//       if (reorienting) {
// 	 if (graphics_info_t::shift_is_pressed) {
// 	    g.reorienting_next_residue(false); // backwards
// 	 } else {
// 	    g.reorienting_next_residue(true); // forwards
// 	 }
//       } else {
// 	 // old/standard simple translation
// 	 if (graphics_info_t::shift_is_pressed) {
// 	    g.intelligent_previous_atom_centring(g.go_to_atom_window);
// 	 } else {
// 	    g.intelligent_next_atom_centring(g.go_to_atom_window);
// 	 }
//       }
//       break;
//    }

//    /* prevent the default handler from being run */
//    // gtk_signal_emit_stop_by_name(GTK_OBJECT(widget), "key_release_event");

//    std::cout << "---------- GTK-FIXME gtk_signal_emit_stop_by_name() " << std::endl;
//    g_signal_emit_by_name(G_OBJECT(widget), "stop", "key_release_event");
//    return TRUE;

  return TRUE;
}

// widget is the glarea.
//
gint idle_contour_function(gpointer data) {

   gint continue_status = 0;
   bool something_changed = false;

   bool is_from_contour_level_change(GPOINTER_TO_INT(data));

   // when there's nothing else to do, update the contour levels
   //
   // then update maps

   for (int imol=0; imol<graphics_info_t::n_molecules(); imol++) {
      if (graphics_info_t::molecules[imol].has_xmap()) { // FIXME or nxmap : needs test for being a map molecule
         int &cc = graphics_info_t::molecules[imol].pending_contour_level_change_count;

         if (cc != 0) {

	          if (cc < 0) {
	             while (cc != 0) {
	                cc++;
	                graphics_info_t::molecules[imol].change_contour(-1);
	             }
	          }

	          if (cc > 0) {
	              while (cc != 0) {
	                 cc--;
	                 graphics_info_t::molecules[imol].change_contour(1);
	              }
	          }

           graphics_info_t g;
           bool really_change_the_map_contours = true;
           if (! is_from_contour_level_change) really_change_the_map_contours = false;
	   g.molecules[imol].update_map(really_change_the_map_contours);
           float map_rmsd = g.molecules[imol].map_sigma();
	   continue_status = 0;
           float cl = g.molecules[imol].contour_level;
           float r = cl/map_rmsd;
           std::cout << "DEBUG:: idle_contour_function() imol: " << imol << " contour level: "
                     << g.molecules[imol].contour_level << " n-rmsd: " << r << std::endl;
           g.set_density_level_string(imol, g.molecules[imol].contour_level);
           std::string s = "Map " + std::to_string(imol) + "  contour_level " +
              coot::util::float_to_string_using_dec_pl(cl, 3) + "  n-rmsd: " +
              coot::util::float_to_string_using_dec_pl(r, 3);
           add_status_bar_text(s.c_str());
           g.display_density_level_this_image = 1;
           something_changed = true;
         }
      }
   }
   // std::cout << "Here with something_changed: " << something_changed << std::endl;

   // is this needed?
   // if (something_changed)
   //    graphics_draw();
   // std::cout << "--- debug:: idle_contour_function() done " << continue_status << std::endl;
   return continue_status;
}



// widget is the glarea.
//
gboolean
animate_idle_ligand_interactions(gpointer data) {

   graphics_info_t g;
   for (int imol=0; imol<graphics_n_molecules(); imol++) {
      if (is_valid_model_molecule(imol)) {
	 if (g.molecules[imol].is_displayed_p()) {
	    g.molecules[imol].draw_animated_ligand_interactions_flag = 1;
	 }
      }
   }
   g.graphics_draw();
   return 1; // don't stop calling this idle function
}


void
draw_axes(GL_matrix &m) {

   graphics_info_t g;


   float mat_p[16];

   if (g.draw_axes_flag) {

      float f = (float) graphics_info_t::graphics_y_size/(float) graphics_info_t::graphics_x_size;

      glPushMatrix();
      glMultMatrixf(m.transpose().get());

      // Guess and fiddle debugging... Now that the scale is before
      // the load (now multiplication) of mat_p we have to adjust
      // where the axes are.  At least now the axes along screen x
      // don't expand as the window is widened.
      float mf = f * 0.121;

      // top left corner
      glTranslatef(-0.26/(mf), +0.85/0.4, 0);

      // This orients the axes to match world orientation
      //
      glMultMatrixf(m.get());
      if (graphics_info_t::use_axes_orientation_matrix_flag)
	 glMultMatrixf(graphics_info_t::axes_orientation.get());

      glMatrixMode(GL_MODELVIEW);

      // If we don't have this, we can't see the axes

      glGetFloatv(GL_MODELVIEW_MATRIX, mat_p);

      glPushMatrix();
      glLoadIdentity();

      glMatrixMode(GL_PROJECTION);
      glPushMatrix();

      glLoadIdentity();

      glScalef(0.4*f, 0.4, 0.4);

      // This doesn't do much/anything.
      // glTranslatef(0.1, 0.1, 0.1);

      // without this, axes are at the centre of the screen and don't rotate
      // glLoadMatrixf(mat_p);
      glMultMatrixf(mat_p);


      glLineWidth(1.0);

//       std::cout << endl;
//       std::cout << "( " << mat_p[0] << " " << mat_p[1]
// 		<< " " <<  mat_p[2] << " " << mat_p[3] << " )\n";
//       std::cout << "( " << mat_p[4] << " " << mat_p[5]
// 		<< " " <<  mat_p[6] << " " << mat_p[7] << " )\n";
//       std::cout << "( " << mat_p[8] << " " << mat_p[9]
// 		<< " " <<  mat_p[10] << " " << mat_p[11] << " )\n";
//       std::cout << "( " << mat_p[12] << " " << mat_p[13]
// 		<< " " <<  mat_p[14] << " " << mat_p[15] << " )\n";
//       std::cout << endl;

      GLfloat col[3] = { 0.55, 0.7, 0.55 };
      glColor3fv(col);

      glBegin(GL_LINES);

      // axes
      glVertex3f(0.0f, 0.0f, 0.0f);
      glVertex3f(0.2f, 0.0f, 0.0f);

      glVertex3f(0.0f, 0.0f, 0.0f);
      glVertex3f(0.0f, 0.2f, 0.0f);

      glVertex3f(0.0f, 0.0f, 0.0f);
      glVertex3f(0.0f, 0.0f, 0.2f);

      // arrowheads:
      glVertex3f(0.2f,   0.0f,  0.0f);
      glVertex3f(0.18f,  0.02f, 0.0f);
      glVertex3f(0.2f,   0.0f,  0.0f);
      glVertex3f(0.18f, -0.02f, 0.0f);

      glVertex3f(0.0f,  0.2f,   0.0f);
      glVertex3f(0.0f,  0.18f,  0.02f);
      glVertex3f(0.0f,  0.2f,   0.0f);
      glVertex3f(0.0f,  0.18f, -0.02f);

      glVertex3f(0.0f,   0.0f,   0.2f);
      glVertex3f(0.02f,  0.0f,   0.18f);
      glVertex3f(0.0f,   0.0f,   0.2f);
      glVertex3f(-0.02f, 0.0f,   0.18f);

      glEnd();

      // glRasterPos3f();
      graphics_info_t::printString_for_axes("x", 0.22, 0.0, 0.0);
      // glRasterPos3f();
      graphics_info_t::printString_for_axes("y", 0.0, 0.22, 0.0);
      // glRasterPos3f();
      graphics_info_t::printString_for_axes("z", 0.0, 0.0, 0.22);

      glMatrixMode(GL_PROJECTION);
      glPopMatrix();
      glMatrixMode(GL_MODELVIEW);
      glPopMatrix();
      glPopMatrix();
   }
}



float r_50(std::string ele) {

   // fix later
   //
   return 1.2;
}

float rad_50_and_prob_to_radius(float rad_50, float prob) {

   // Note the prob*(prob-100.0) part is scaled up by a factor of
   // 100*100... and the (50-prob) is scaled by a factor of 100, so
   // multiply the top by 100.0 to compensate.
   //
   // The leading factor is to stop the radius moving too fast.
   //
   return rad_50*exp(0.3*100.0*(50.0-prob)/(prob*(prob-100.0)));

}

// This function is way too long.  It needs to be split it into
// internal bits.
//
gint glarea_button_press(GtkWidget *widget, GdkEventButton *event) {

   graphics_info_t info;

   bool was_a_double_click = false;
   if (event->type==GDK_2BUTTON_PRESS)
      was_a_double_click = true;

   int x_as_int, y_as_int;
   GdkModifierType state = GDK_MODIFIER_RESERVED_19_MASK; // something

   GdkWindow *window = 0; // was: gtk_widget_get_window(GTK_WIDGET(widget));

   // gdk_window_get_pointer(window, &x_as_int, &y_as_int, &state);

   // GdkWindow *window = gtk_widget_get_window(widget);
   // GdkDeviceManager *device_manager = gdk_display_get_default_seat(gtk_widget_get_display(widget));
   // GdkDevice *device = gdk_device_manager_get_client_pointer(device_manager);
   // gdk_window_get_device_position(window, device, &x_as_int, &y_as_int, &state);

   std::cout << "FIX button press (frustrated) " << std::endl;

   info.SetMouseBegin(event->x, event->y);
   info.SetMouseClicked(event->x, event->y);

   GdkModifierType my_button1_mask = info.gdk_button1_mask();
   GdkModifierType my_button2_mask = info.gdk_button2_mask();
   // GdkModifierType my_button3_mask = info.gdk_button3_mask();

//    std::cout << "debug:: mouse button " << event->button << " was pressed at ("
// 	     << event->x << "," << event->y << ")" << std::endl;

   // cout << gtk_events_pending() << " events pending..." << endl;

   // On refinement atom selection it seems double clicks generate an
   // extra event (i.e. we get 3 click events of which the last one
   // seems to be a synthetic double click event).  Let's bin the last
   // event.
   //
   if (event->type==GDK_2BUTTON_PRESS ||
       event->type==GDK_3BUTTON_PRESS) {

      if (false)
	 printf("INFO:: ---- Considering this %s click on button %d\n",
		event->type==GDK_2BUTTON_PRESS ? "double" : "triple",
		event->button);

      info.check_if_moving_atom_pull(true); // doesn't return a value, but can remove pull restraints

      if (state & my_button1_mask) {
	 // instead do a label atom:
	 pick_info nearest_atom_index_info = atom_pick(event);
	 if (nearest_atom_index_info.success == GL_TRUE) {
	    int im = nearest_atom_index_info.imol;
            if (is_valid_model_molecule(im)) {
               info.molecules[im].add_to_labelled_atom_list(nearest_atom_index_info.atom_index);
               mmdb::Residue *r = info.molecules[im].atom_sel.atom_selection[nearest_atom_index_info.atom_index]->residue;
               std::string alt_conf = info.molecules[im].atom_sel.atom_selection[nearest_atom_index_info.atom_index]->altLoc;
               info.setup_graphics_ligand_view(im, r, alt_conf);
            } else {
               if (im == -1) {
                  // this is intermediate atoms double click, but let's add protection any way
                  mmdb::Atom *at = info.get_moving_atom(nearest_atom_index_info);
                  if (at) {
                     std::string s= "Atom picked: ";
                     s += at->GetChainID();
                     s += " ";
                     s += std::to_string(at->residue->GetSeqNum());
                     s += " ";
                     s += at->residue->GetInsCode();
                     s += " (";
                     s += at->residue->GetResName();
                     s += ") ";
                     s += at->GetAtomName();
                     add_status_bar_text(s.c_str());
                  }
               }
            }
	 } else {
	    if (graphics_info_t::show_symmetry) {
	       coot::Symm_Atom_Pick_Info_t symm_atom_info = info.symmetry_atom_pick();
	       if (symm_atom_info.success == GL_TRUE) {
		  info.molecules[symm_atom_info.imol].add_atom_to_labelled_symm_atom_list(symm_atom_info.atom_index,
											  symm_atom_info.symm_trans,
											  symm_atom_info.pre_shift_to_origin);
	       }
	    }
	 }
	 graphics_draw();
      }
      return TRUE;
   }


   if (state & my_button1_mask) {

      if (info.shift_is_pressed == 1) {

	 // label atom
	 //
	 pick_info nearest_atom_index_info;
	 nearest_atom_index_info = atom_pick(event);

	 if ( nearest_atom_index_info.success == GL_TRUE ) {

	    int im = nearest_atom_index_info.imol;
            if (is_valid_model_molecule(im)) {
               info.molecules[im].add_to_labelled_atom_list(nearest_atom_index_info.atom_index);
               mmdb::Residue          *r = info.molecules[im].atom_sel.atom_selection[nearest_atom_index_info.atom_index]->residue;
               std::string alt_conf = info.molecules[im].atom_sel.atom_selection[nearest_atom_index_info.atom_index]->altLoc;
               info.setup_graphics_ligand_view(im, r, alt_conf);
               info.graphics_draw();
            }

	 } else {

	    // we should only try to pick symmetry atoms if they are being
	    // displayed.
	    // (note we are still in shift left mouse mode)

	    if ( info.show_symmetry == 1 ) {

	       std::cout << "Trying symmetry label\n";
	       coot::Symm_Atom_Pick_Info_t symm_atom_info = info.symmetry_atom_pick();

	       if (symm_atom_info.success == GL_TRUE) {

		  info.molecules[symm_atom_info.imol].add_atom_to_labelled_symm_atom_list(symm_atom_info.atom_index,
											  symm_atom_info.symm_trans,
											  symm_atom_info.pre_shift_to_origin);
		  info.graphics_draw();

		  // need equivalent for this?
		  // clear_symm_atom_info(symm_atom_info);

	       } else {
		  std::cout << "Symmetry label failed\n";
	       }
	    }
	 }

      } else {

	 // Left mouse, and not shift-left-mouse:

	 // (this is the conventional case)
	 //
	 // and if so, do the regularization or refinement
	 // or angle and distance geometries.
	 //
	 int iv = info.check_if_in_range_defines(event, state);
	 if (! iv)
	    info.check_if_moving_atom_pull(false); // and if so, set it up (it
	                                           // executes on *motion* not a button press event).  Also,
                                                   // remove any on-going drag-refine-idle-function.
                                                   // Note that this doesn't only apply while pulling on
                                                   // an atom during refinement, it is for Rotation/Translation too

	 // check_if_moving_atom_pull sets in_moving_atoms_drag_atom_mode_flag.

      }  // shift is pressed
   }     // button 1


   if (state & my_button2_mask) {
      // std::cout << "Nothing doin' at the moment" << std::endl;
   }

   return TRUE;
}

gint glarea_button_release(GtkWidget *widget, GdkEventButton *event) {

   graphics_info_t g;
   if (graphics_info_t::in_moving_atoms_drag_atom_mode_flag) {
      g.unset_moving_atoms_currently_dragged_atom_index();
      g.do_post_drag_refinement_maybe();
   } else {

      int x_as_int = 0, y_as_int = 0;
      GdkModifierType state;
      GdkWindow *window = 0; // was gtk_widget_get_window(widget);

      std::cout << "FIX button release (frustrated) " << std::endl;
      // GdkDeviceManager *device_manager = gdk_display_get_default_seat(gtk_widget_get_display(widget));
      // GdkDevice *device = gdk_device_manager_get_client_pointer(device_manager);
      // gdk_window_get_device_position(window, device, &x_as_int, &y_as_int, &state);

      GdkModifierType my_button2_mask = g.gdk_button2_mask();

      if (event->button == 2) {

	 // move this block into glarea_button_release (now done of course)

	 // Atom picking (recentre)

	 // check that the mouse hasn't moved much
	 // c.f. event->x vs mouse_clicked_begin.first
	 //      event->y vs mouse_clicked_begin.second

	 pick_info nearest_atom_index_info = atom_pick(event);

	 double delta_x = g.GetMouseClickedX() - x_as_int;
	 double delta_y = g.GetMouseClickedY() - y_as_int;

	 if (false) {
	    std::cout << "MouseBegin " << g.GetMouseBeginX() << " " << g.GetMouseBeginY()
		      << " now " << x_as_int << " " << y_as_int << std::endl;
	    std::cout << "mouse deltas " << delta_x << " " << delta_y << std::endl;
	 }

	 if (std::abs(delta_x) < 10.0) {
	    if (std::abs(delta_y) < 10.0) {

	       if (nearest_atom_index_info.success == GL_TRUE) {

		  int im = nearest_atom_index_info.imol;
                  if (is_valid_model_molecule(nearest_atom_index_info.imol)) {
                     std::cout << "INFO:: recentre: clicked on imol: " << im << std::endl;
                     mmdb::Atom *at = g.molecules[im].get_atom(nearest_atom_index_info);
                     g.setRotationCentre(nearest_atom_index_info.atom_index, nearest_atom_index_info.imol);
                     g.sequence_view_highlight_residue_maybe(at, g.get_sequence_view_is_displayed(im));
                  } else {
                     mmdb::Atom *at = g.get_moving_atom(nearest_atom_index_info);
                     if (at) {
                        coot::Cartesian c(at->x, at->y, at->z);
                        g.setRotationCentre(c);
                     }
                  }

		  // Lets display the coordinate centre change
		  // *then* update the map, so we can see how fast
		  // the map updating is.
		  //
		  g.graphics_draw();

		  g.post_recentre_update_and_redraw();

	       } else {

		  std::cout << "Model atom pick failed. " << std::endl;

		  // Only try to pick symmetry atoms if the symmetry atoms
		  // are being displayed.

		  if (g.show_symmetry == 1) {

		     std::cout << "Trying symmetry pick" << std::endl;
		     coot::Symm_Atom_Pick_Info_t symm_atom_info = g.symmetry_atom_pick();
		     if (symm_atom_info.success == GL_TRUE) {

			std::cout << "Found Symmetry atom pick" << std::endl;

			// info.setRotationCentre(symm_atom_info.Hyb_atom());
			std::pair<symm_trans_t, Cell_Translation> symtransshiftinfo(symm_atom_info.symm_trans,
										    symm_atom_info.pre_shift_to_origin);
			g.setRotationCentre(translate_atom_with_pre_shift(g.molecules[symm_atom_info.imol].atom_sel,
									  symm_atom_info.atom_index,
									  symtransshiftinfo));

			// clear_symm_atom_info(symm_atom_info);

			for (int ii=0; ii<g.n_molecules(); ii++) {
			   g.molecules[ii].update_symmetry();
			}
			for (int ii=0; ii<g.n_molecules(); ii++) {
			   g.molecules[ii].update_clipper_skeleton();
			   g.molecules[ii].update_map(graphics_info_t::auto_recontour_map_flag);
			}
			g.graphics_draw();

		     } else {

			std::cout << "Symmetry atom pick failed." << std::endl;

		     }
		  }
	       }
	    }
	 }
      }
   }
   graphics_info_t::in_moving_atoms_drag_atom_mode_flag = 0;
   return TRUE;
}

// Francois says to remove this:
//
// #if defined(WINDOWS_MINGW) || defined(_MSC_VER)

// gint glarea_scroll_event(GtkWidget *widget, GdkEventScroll *event) {

//    if (event->direction == SCROLL_UP)
//       handle_scroll_event(1);
//    if (event->direction == SCROLL_DOWN)
//       handle_scroll_event(0);
//    return TRUE;
// }
// #endif


gint glarea_scroll_event(GtkWidget *widget, GdkEventScroll *event) {

   graphics_info_t info;
   bool handled = false;
   if (info.control_is_pressed) {
      if (info.shift_is_pressed) {
	 if (event->direction == GDK_SCROLL_UP)
	    change_model_molecule_representation_mode(1);
	 if (event->direction == GDK_SCROLL_DOWN)
	    change_model_molecule_representation_mode(0);
	 handled = true;
      } else {
         // maybe we change the proportional editing (pull atom neighbour displacement) radius
         bool dir = true;
         if (event->direction == GDK_SCROLL_DOWN)
            dir = false;
         info.pull_restraint_neighbour_displacement_change_max_radius(dir);
         info.graphics_draw();
         handled = true;
      }
   } else {
      // std::cout << "control is not pressed " << std::endl;
   }

   if (! handled) {
      if (event->direction == GDK_SCROLL_UP)
	 handle_scroll_density_level_event(1);
      if (event->direction == GDK_SCROLL_DOWN)
	 handle_scroll_density_level_event(0);
   }
   return TRUE;
}

void handle_scroll_density_level_event(int scroll_up_down_flag) {

   std::cout << "debug:: handle_scroll_density_level_event() " << scroll_up_down_flag << std::endl;

   graphics_info_t info;

   GdkEvent *peek_event = gdk_event_peek();
   if (peek_event) {
      std::cout << "peaking found an event!" << std::endl;
   }

   int imol_map_for_scroll = info.scroll_wheel_map;
   std::vector<int> dm = info.displayed_map_imols();
   if (std::find(dm.begin(), dm.end(), imol_map_for_scroll) == dm.end()) {
      // imol_map_for_scroll is not visible, choose another one
      if (dm.size() > 0)
      imol_map_for_scroll = dm[0];
   }

   if (scroll_up_down_flag == 1) {
      if (graphics_info_t::do_scroll_by_wheel_mouse_flag) {
         if (imol_map_for_scroll>=0) {
            // short int istate = info.molecules[s].change_contour(1);
            info.molecules[imol_map_for_scroll].pending_contour_level_change_count++;
            GSourceFunc f = idle_contour_function;
            // int contour_idle_token = g_idle_add(f, info.glarea);
            float cl = info.molecules[imol_map_for_scroll].contour_level;
            info.set_density_level_string(imol_map_for_scroll, cl);
            info.display_density_level_this_image = 1;
         } else {
            std::cout << "WARNING: No map - Can't change contour level.\n";
         }
      }
   }


   if (scroll_up_down_flag == 0) {
      if (graphics_info_t::do_scroll_by_wheel_mouse_flag) {
         if (imol_map_for_scroll>=0) {
            // short int istate = info.molecules[s].change_contour(-1);
            info.molecules[imol_map_for_scroll].pending_contour_level_change_count--;
            int contour_idle_token = g_idle_add(idle_contour_function, info.glareas[0]);
            float cl = info.molecules[imol_map_for_scroll].contour_level;
            info.set_density_level_string(imol_map_for_scroll, cl);
            info.display_density_level_this_image = 1;
         } else {
            std::cout << "WARNING: No map - Can't change contour level.\n";
         }
      }
   }

}


void test_object() {

   // a solid triangle
   //
   if (1) {


      glEnable (GL_BLEND);
      glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

      glBegin(GL_TRIANGLE_STRIP);

      glColor4f(0.6, 0.6, 0.7, 0.8);
      glNormal3f(1.0, 0.0, 0.0);
      glVertex3f(10.0,  0.0,  0.0);
      glVertex3f( 0.0, 10.0,  0.0);
      glVertex3f(10.0, 10.0,  0.0);

      glNormal3f(0.0, 1.0, 0.0);
      glVertex3f(-5.0, 0.0,  -10.0);
      glVertex3f( 0.0, 10.0,  0.0);
      glVertex3f(10.0, 10.0,  0.0);
      glEnd();

   }

   if (0) {
      int n = 40;
      float f = 1.0;

      std::vector<int> colours;
      colours.push_back(GREEN_BOND);
      colours.push_back(BLUE_BOND);
      colours.push_back(GREY_BOND);

      for (unsigned int ic=0; ic<colours.size(); ic++) {
	 set_bond_colour(colours[ic]);
	 if (ic == 2)
	    glLineWidth(2.0);
	 else
	    glLineWidth(5.0);

	 glBegin(GL_LINES);
	 for (int i=0; i<n; i++) {
	    for (int j=0; j<n; j++) {
	       glVertex3f(i*f, j*f, ic);
	       glVertex3f(i*f, j*f, ic+1);
	    }
	 }
	 glEnd();
      }
   }
}


void
set_bond_colour(int i) {

   if (false)
      std::cout << "globjects.cc set_bond_colour() idx: " << i << " vs "
		<< " green "   << GREEN_BOND << " "
		<< " blue "    << BLUE_BOND << " "
		<< " red "     << RED_BOND << " "
		<< " yellow "  << YELLOW_BOND << " "
		<< " grey "    << GREY_BOND << " "
		<< " H-grey "  << HYDROGEN_GREY_BOND << " "
		<< " magenta " << MAGENTA_BOND << " "
		<< std::endl;

   if (background_is_black_p()) {
      switch (i) {
      case CARBON_BOND:
	 glColor3f (0.2, 0.7, 0.1);
	 break;
      case GREEN_BOND:
	 glColor3f (0.0, 0.7, 0.0);
	 break;
      case BLUE_BOND:
	 glColor3f (0.2, 0.2, 0.8);
	 break;
      case RED_BOND:
	 glColor3f (0.8, 0.1, 0.1);
	 break;
      case YELLOW_BOND:
	 glColor3f (0.7, 0.7, 0.0);
	 break;
      case GREY_BOND:
	 glColor3f (0.7, 0.7, 0.7);
	 break;
      case HYDROGEN_GREY_BOND:
	 glColor3f (0.6, 0.6, 0.6);
	 break;
      case MAGENTA_BOND:
	 glColor3f (0.8, 0.1, 0.8);
	 break;
      case DARK_GREEN_BOND:
	 glColor3f (0.05, 0.69, 0.05);
	 break;
      case DARK_ORANGE_BOND:
	 glColor3f (0.7, 0.7, 0.05);
	 break;
      case DARK_BROWN_BOND:
	 glColor3f (0.5, 0.5, 0.1);
	 break;
      default:
	 glColor3f (0.7, 0.8, 0.8);
      }
   } else {
      // Are you sure that this is begin executed (and the colour-state not overwritten?)
      // How about set_bond_colour_by_mol_no()?
      switch (i) {
      case CARBON_BOND:
	 glColor3f (0.2, 0.6, 0.0);
	 break;
      case GREEN_BOND:
	 glColor3f (0.05, 0.6, 0.05);
	 break;
      case BLUE_BOND:
	 glColor3f (0.1, 0.1, 0.8);
	 break;
      case RED_BOND:
	 glColor3f (0.6, 0.05, 0.5);
	 break;
      case YELLOW_BOND:
	 glColor3f (0.35, 0.35, 0.0);
	 break;
      case GREY_BOND:
	 glColor3f (0.3, 0.3, 0.3);
	 break;
      case HYDROGEN_GREY_BOND:
	 glColor3f (0.4, 0.4, 0.4);
	 break;
      case MAGENTA_BOND:
	 glColor3f (0.7, 0.1, 0.7);
	 break;
      case DARK_GREEN_BOND:
	 glColor3f (0.05, 0.69, 0.05);
	 break;
      case DARK_ORANGE_BOND:
	 glColor3f (0.7, 0.7, 0.05);
	 break;
      case DARK_BROWN_BOND:
	 glColor3f (0.5, 0.5, 0.1);
	 break;
      default:
	 glColor3f (0.3, 0.4, 0.4);
      }
   }
}



// // we get passed a number that is 0.0 - 1.0 inclusive (the higher the
// // number, the more "core" is the skeleton).
// //
// void set_skeleton_bond_colour(float f)
// {
//    // glColor3f (0.2+0.8*f, 0.4+0.5*f, 0.8-0.8*f);
//    // glColor3f (0.4+0.6*f, f, 0.3);

//    graphics_info_t g;
//    glColor3f(0.1+0.6*f*g.skeleton_colour[0],
// 	     0.1+0.9*f*g.skeleton_colour[1],
// 	     0.1+0.2*f*g.skeleton_colour[2]);

// //    glColor3f(0.2+0.8*g.skeleton_colour[0],
// // 	     0.3+0.5*g.skeleton_colour[1],
// // 	     0.0+0.2*g.skeleton_colour[2]);
// }

// i goes upto bond_box.num_colours
//
void set_skeleton_bond_colour_random(int i, const std::vector<std::vector<float> > &colour_table) {

   glColor3f(0.2+0.8*colour_table[i][0],
	     0.2+0.8*colour_table[i][1],
	     0.2+0.8*colour_table[i][2]);

}

// double
// combine_colour(double v, int col_part_index) {
//    // col_part_index is 0,1,2 for red gree blue components of the colour

//    graphics_info_t g;
//    double w = g.symm_colour_merge_weight[0];

//    return w*g.symm_colour[0][col_part_index] + v*(1.0-w);
// }


// // redundant now, I hope.  It's internal to molecule_class_info_t.
// void
// set_symm_bond_colour(int i) {

//    switch (i) {
//    case green:
//       glColor3f (combine_colour(0.1,0),
// 		 combine_colour(0.8,1),
// 		 combine_colour(0.1,2));
//       break;
//    case blue:
//       glColor3f (combine_colour(0.2,0),
// 		 combine_colour(0.2,1),
// 		 combine_colour(0.8,2));
//       break;
//    case red:
//       glColor3f (combine_colour(0.8,0),
// 		 combine_colour(0.1,1),
// 		 combine_colour(0.1,2));
//       break;
//    case yellow:
//       glColor3f (combine_colour(0.7,0),
// 		 combine_colour(0.7,1),
// 		 combine_colour(0.0,2));
//       break;

//    default:
//       glColor3f (combine_colour(0.7, 0),
// 		 combine_colour(0.8, 1),
// 		 combine_colour(0.8, 2));
//    }
// }

#endif // 0

