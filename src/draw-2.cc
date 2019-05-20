
#include <iostream>
#include <string>
#include <fstream>
#include <gtk/gtk.h>
#include <epoxy/gl.h>

#include "draw.hh"

// static int programID_global = -1;
// static int location_global = -1;
// GLuint VertexArrayID = -1;

void init_shaders() {

   std::cout << "----------- parse and create shader " << std::endl;
   shader_program_source sps = parse_shader("Basic.shader");
   unsigned int programID = CreateShader(sps.VertexSource, sps.FragmentSource);
   programID_global = programID;
   std::cout << "----------- created shader program " << programID << std::endl;

   glBindAttribLocation(programID, 0, "position");

}

void init_buffers() {

   {
      float positions[12] = {
	 -0.5,  -0.5, -0.0,
   	 -0.5,   0.5, -0.0,
  	  0.5,   0.5, -0.0,
	  0.5,  -0.5, -0.0
      };

      unsigned int indices[8] { 0,1,1,2,2,3,3,0 };

      // GLuint VertexArrayID;
      glGenVertexArrays(1, &VertexArrayID);
      glBindVertexArray(VertexArrayID);

      GLuint vertexbuffer;
      glGenBuffers(1, &vertexbuffer);
      glBindBuffer(GL_ARRAY_BUFFER, vertexbuffer);
      glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 12, &positions[0], GL_STATIC_DRAW);
      glEnableVertexAttribArray(0);
      glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

      unsigned int ibo;
      glGenBuffers(1, &ibo);
      glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned int) * 8, &indices[0], GL_STATIC_DRAW);

      // int ul = glGetUniformLocation(programID, "u_Color");
      // std::cout << "Got glGetUniformLocation for u_Color " << ul << std::endl;
      // location_global = ul;

   }
}

void draw_triangle(GtkGLArea *glarea) {

   glLineWidth(3.0);  // GLv4 antialiasing - OpenGL implementations are not required to support this

   // To see the possible values of the line width in aliased mode:
   // GLfloat line_width_max_min[2] = {0.0f, 0.0f};
   // glGetFloatv(GL_ALIASED_LINE_WIDTH_RANGE, lineWidthRange);
   // This may not be possible in GL_LINE_SMOOTH mode.

   glBindVertexArray(VertexArrayID);
   glUseProgram(programID_global);
   glDrawElements(GL_LINES, 8, GL_UNSIGNED_INT, nullptr);
   glBindVertexArray(0); // unbind
   glUseProgram(0);

}

GtkWidget *my_gtkglarea(GtkWidget *vbox) {

   GtkWidget *w = gtk_gl_area_new();
   gtk_widget_set_size_request(w, 400, 400);
   gtk_box_pack_start(GTK_BOX(vbox), w, TRUE, TRUE, 2);
   return w;
}

void
on_glarea_realize(GtkGLArea *glarea) {

   std::cout << "realize!" << std::endl;
   gtk_gl_area_make_current(glarea);
   init_shaders();
   init_buffers();
   glEnable(GL_LINE_SMOOTH);

}


gboolean
on_glarea_render(GtkGLArea *glarea) {

   std::cout << "render!" << std::endl;
   glClearColor (0.5, 0.5, 0.3, 1.0);
   glClear (GL_COLOR_BUFFER_BIT);
   draw_triangle(glarea);
   glFlush ();

  return FALSE;
}

void
on_glarea_resize(GtkGLArea *glarea, gint width, gint height) {

   std::cout << "resize!" << std::endl;
}

gboolean
on_glarea_scroll(GtkWidget *widget, GdkEventScroll *event) {

   std::cout << "scroll!" << std::endl;
   return TRUE;
}

gboolean
on_glarea_button_press(GtkWidget *widget, GdkEventButton *event) {

   std::cout << "press!" << std::endl;
   // Use gtk_widget_queue_draw(GtkGLArea *gl_drawing_area);
   return TRUE;
}

gboolean
on_glarea_button_release(GtkWidget *widget, GdkEventButton *event) {

   std::cout << "release!" << std::endl;
   return TRUE;
}

gboolean
on_glarea_motion_notify(GtkWidget *wdiget, GdkEventMotion *event) {

   std::cout << "motion!" << std::endl;

   return TRUE;
}

void my_glarea_add_signals_and_events(GtkWidget *glarea) {

   gtk_widget_add_events(glarea, GDK_SCROLL_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON_PRESS_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON_RELEASE_MASK);
   gtk_widget_add_events(glarea, GDK_BUTTON1_MOTION_MASK);

   g_signal_connect(glarea, "realize", G_CALLBACK(on_glarea_realize), NULL);
   g_signal_connect(glarea, "render",  G_CALLBACK(on_glarea_render),  NULL);
   g_signal_connect(glarea, "resize",  G_CALLBACK(on_glarea_resize),  NULL);
   g_signal_connect(glarea, "scroll-event",          G_CALLBACK(on_glarea_scroll),         NULL);
   g_signal_connect(glarea, "button-press-event",    G_CALLBACK(on_glarea_button_press),   NULL);
   g_signal_connect(glarea, "button-release-event",  G_CALLBACK(on_glarea_button_release), NULL);
   g_signal_connect(glarea, "motion-notify-event",   G_CALLBACK(on_glarea_motion_notify),  NULL);

}
