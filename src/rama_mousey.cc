/* src/rama_mousey.cc
 * 
 * Copyright 2001, 2002, 2003, 2004, 2005, 2006 by The University of York
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#ifdef HAVE_GOOCANVAS

#include <iostream>

#include "rama_mousey.hh" // has gtk/gtk.h
#include "interface.h"    // needs gtk/gtk.h

#include "c-interface-gtk-widgets.h"
#include "c-interface.h"

#include "rama_plot.hh"


extern "C" G_MODULE_EXPORT void
on_dynarama2_window_destroy(GObject *caller, gpointer user_data) {

   // 20220315-PE this is no longer connected.

   // std::cout << "on_dynarama2_window_destroy() caller-object: " << caller << std::endl;

   // maybe no callback from builder for mainwindow!?
   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "rama_plot"));
   if (!plot) {
      std::cout<<"debug:: in on_dynarama2_window_destroy() failed to get the plot from " << canvas <<std::endl;
   } else {
      if (plot->is_stand_alone()) {
         // gtk_exit(0);
         std::cout << "exit here " << std::endl;
      }  else {
         int imol = plot->molecule_number();
         if (imol >= 0) {
            set_dynarama_is_displayed(0, imol); // which frees/deletes the
            // memory of the user data.
         }
      }
   }

}

extern "C" G_MODULE_EXPORT gboolean
on_dynarama2_window_configure_event(GtkWidget       *widget,
                                   GdkEventConfigure *event,
                                   gpointer         user_data) {

   // std::cout << "DEBUG:: on_dynarama2_window_configure_event() --- start --- " << std::endl;

   // maybe no callback from builder for mainwindow!?
   // or use the "new one"
   // do we need this then?

   return FALSE;
}


extern "C" G_MODULE_EXPORT void
on_dynarama2_ok_button_clicked(GtkButton *button, gpointer user_data) {

   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "rama_plot"));
   if (!plot) {
      std::cout<<"debug:: on-dynarama_ok_button() failed to get the plot from " << canvas <<std::endl;
   } else {
      if (plot->is_stand_alone()) {
         // gtk_exit(0);
         std::cout << "exit here" << std::endl;
      } else {
         int imol = plot->molecule_number();
         if (imol == -9999)
            accept_phi_psi_moving_atoms();
         gtk_widget_hide(plot->dynawin); // 20220315-PE hide, not destroy.
      }
   }

}

// extern "C" G_MODULE_EXPORT void
// on_dynarama2_cancel_button_clicked(GtkButton *button, gpointer user_data) {

//    GtkWidget *canvas = GTK_WIDGET(user_data);
//    coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "rama_plot"));
//    if (!plot) {
//       std::cout<<"debug:: on_dynarama2_cancel_button_clicked() failed to get the plot from " << canvas <<std::endl;
//    } else {
//       if (plot->is_stand_alone()) {
//          // gtk_exit(0);
//          std::cout << "exit here in on_dynarama2_cancel_button_clicked() " << std::endl;
//       } else {
//          int imol = plot->molecule_number();
//          if (imol == -9999)
//             clear_moving_atoms_object();
//          gtk_widget_destroy(plot->dynawin);
//       }
//    }
// }

extern "C" G_MODULE_EXPORT void
on_dynarama2_close_button_clicked_gtkbuilder_callback(GtkButton *button, gpointer user_data) {

   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(button), "rama_plot"));
   if (!plot) {
      std::cout << "debug:: on_dynarama2_close_button_clicked() failed to get the plot from button "
                << button << std::endl;
   } else {
      if (plot->is_stand_alone()) {
         // gtk_exit(0);
         std::cout << "exit here in on_dynarama2_cancel_button_clicked() " << std::endl;
      } else {
         int imol = plot->molecule_number();
         if (imol == -9999)
            clear_moving_atoms_object();
         gtk_widget_hide(plot->dynawin);

         std::cout << "Hiding Rama plot dynawin for imol " << imol << " " << plot->dynawin << std::endl;
      }
   }
}


extern "C" G_MODULE_EXPORT void
on_kleywegt_apply_chain_button_clicked_gtkbuilder_callback(GtkButton *button, gpointer user_data) {

   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "rama_plot"));
   if (!plot) {
      std::cout << "error:: on_kleywegt_apply_chain_button_clicked()))  failed to get the plot from "
                << canvas <<std::endl;
   } else {
      plot->update_kleywegt_plot();
   }
}

extern "C" G_MODULE_EXPORT void
on_dynarama2_outliers_only_togglebutton_toggled_gtkbuilder_callback(GtkToggleButton *button, gpointer user_data) {

   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "rama_plot"));

   std::cout << "debug:: on_dynarama2_outliers_only_togglebutton_toggled_gtkbuilder_callback canvas " << canvas
             << " plot " << plot << std::endl;
   if (plot) {
      int state = gtk_toggle_button_get_active(button);
      plot->show_outliers_only(state);
   } else {
      std::cout << "debug:: in on_dynarama2_outliers_only_togglebutton_toggled_gtkbuilder_callback() null plot"
                << " from canvas " << canvas << std::endl;
   }
}

extern "C" G_MODULE_EXPORT void
on_dynarama_selection_checkbutton_toggled_gtkbuilder_callback(GtkToggleButton *button, gpointer user_data){

   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "rama_plot"));
   if (plot) {
      int state = gtk_toggle_button_get_active(button);
      plot->show_selection_widget(state);
   }
}

extern "C" G_MODULE_EXPORT void
on_dynarama_selection_entry_activate_gtkbuilder_callback(GtkEntry *entry, gpointer  user_data) {

   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "rama_plot"));
   if (plot) {
      plot->apply_selection_from_widget();
   }
}

extern "C" G_MODULE_EXPORT void
on_dynarama_selection_apply_button_clicked_gtkbuilder_callback(GtkButton *button, gpointer user_data){

   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "rama_plot"));
   if (plot) {
      plot->apply_selection_from_widget();
   }
}

extern "C" G_MODULE_EXPORT void
on_psi_axis_classic_radioitem_toggled(GtkToggleButton *button, gpointer user_data) {

   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "rama_plot"));
   if (plot) {
      plot->psi_axis_changed();
   }
}

extern "C" G_MODULE_EXPORT void
on_psi_axis_paule_radioitem_toggled(GtkToggleButton *button, gpointer user_data) {

   // shouldnt be needed.
//   GtkWidget *canvas = GTK_WIDGET(user_data);
//   coot::rama_plot *plot = static_cast<coot::rama_plot *> (gtk_object_get_user_data(GTK_OBJECT(canvas)));
//   if (plot) {
//      plot->axis_type_change();
//   }
}

extern "C" G_MODULE_EXPORT void
on_dynarama2_zoom_resize_togglebutton_toggled(GtkToggleButton *button, gpointer user_data) {

   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "rama_plot"));
   if (plot) {
      int state = gtk_toggle_button_get_active(button);
      plot->resize_mode_changed(state);
   }
}

// Menu callbacks
extern "C" G_MODULE_EXPORT void
on_rama_open_menuitem_activate(GtkMenuItem *item, gpointer user_data) {


   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "rama_plot"));
   if (plot) {
      gtk_widget_show(plot->rama_open_filechooserdialog);
   } else {
      std::cout<< "BL ERROR:: failed to get a plot" <<std::endl;
   }
}

extern "C" G_MODULE_EXPORT void
on_rama_print_menuitem_activate(GtkMenuItem *item, gpointer user_data) {

   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "rama_plot"));
   if (plot) {
      const gchar *file_name = "dynarama.pdf";
      gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(plot->rama_export_as_pdf_filechooserdialog),
                                        file_name);
      gtk_widget_show(plot->rama_export_as_pdf_filechooserdialog);
   } else {
      std::cout<< "BL ERROR:: failed to get a plot" <<std::endl;
   }
}

extern "C" G_MODULE_EXPORT void
on_rama_save_as_png_menuitem_activate(GtkMenuItem *item, gpointer user_data) {

   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "rama_plot"));
   if (plot) {
      const gchar *file_name = "dynarama.png";
      gtk_file_chooser_set_current_name(GTK_FILE_CHOOSER(plot->rama_export_as_png_filechooserdialog),
                                        file_name);
      gtk_widget_show(plot->rama_export_as_png_filechooserdialog);
   } else {
      std::cout<< "BL ERROR:: failed to get a plot" <<std::endl;
   }
}


extern "C" G_MODULE_EXPORT void
on_rama_close_menuitem_activate(GtkMenuItem *item, gpointer user_data) {

   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "rama_plot"));
   if (plot)
      gtk_widget_destroy(plot->dynawin);
   else
      std::cout<< "BL ERROR:: failed to get a plot" <<std::endl;
}

extern "C" G_MODULE_EXPORT void
on_rama_radiomenuitem_toggled(GtkCheckMenuItem *checkmenuitem, gpointer user_data) {

   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "rama_plot"));
   if (plot) {
      plot->plot_type_changed();
   }
}

extern "C" G_MODULE_EXPORT void
on_kleywegt_radiomenuitem_toggled(GtkCheckMenuItem *checkmenuitem, gpointer user_data) {

   // do we need this then?

}

extern "C" G_MODULE_EXPORT void
on_outliers_only_menuitem_toggled(GtkCheckMenuItem *checkmenuitem, gpointer user_data) {

   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "rama_plot"));
   if (plot) {
      int state;
      state = gtk_check_menu_item_get_active (checkmenuitem);
      plot->show_outliers_only(state);
   } else {
      std::cout<< "BL ERROR:: failed to get a plot" <<std::endl;
   }
}

extern "C" G_MODULE_EXPORT void
on_zoom_in_activate(GtkMenuItem *item, gpointer user_data) {

   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "rama_plot"));
   if (plot) {
      plot->zoom_in();
   } else {
      std::cout<< "BL ERROR:: failed to get a plot" <<std::endl;
   }
}

extern "C" G_MODULE_EXPORT void
on_zoom_out_activate(GtkMenuItem *item, gpointer user_data) {

   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "rama_plot"));
   if (plot) {
      plot->zoom_out();
   } else {
      std::cout<< "BL ERROR:: failed to get a plot" <<std::endl;
   }
}


extern "C" G_MODULE_EXPORT void
on_zoom_resize_menuitem_toggled(GtkCheckMenuItem *checkmenuitem, gpointer user_data) {

   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "rama_plot"));
   if (plot) {
      int state;
      state = gtk_check_menu_item_get_active (checkmenuitem);
      plot->resize_mode_changed(state);
   } else {
      std::cout<< "BL ERROR:: failed to get a plot" <<std::endl;
   }
}

extern "C" G_MODULE_EXPORT void
on_rama_about_menuitem_activate(GtkMenuItem *item, gpointer user_data) {

   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "rama_plot"));
   if (plot) {
      GtkWidget *about = plot->about_dialog;
      gtk_widget_show(about);
   } else {
      std::cout<< "BL ERROR:: failed to get a plot" <<std::endl;
   }
}

// About dialog
extern "C" G_MODULE_EXPORT void
on_rama_aboutdialog1_close(GtkDialog *dialog, gpointer user_data){

   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "rama_plot"));
   if (plot) {
      GtkWidget *about = plot->about_dialog;
      gtk_widget_hide(about);
   } else {
      std::cout<< "BL ERROR:: failed to get a plot" <<std::endl;
   }
}

extern "C" G_MODULE_EXPORT void
on_rama_aboutdialog1_response(GtkDialog *dialog, gint response_id, gpointer user_data) {

   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "rama_plot"));
   if (plot) {
      GtkWidget *about = plot->about_dialog;
      gtk_widget_hide(about);
   } else {
      std::cout<< "BL ERROR:: failed to get a plot" <<std::endl;
   }
}


// file chooser responses
extern "C" G_MODULE_EXPORT void
on_rama_export_as_pdf_filechooserdialog_close(GtkDialog *dialog, gpointer user_data){

   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "rama_plot"));
   if (plot) {
      GtkWidget *w = plot->rama_export_as_pdf_filechooserdialog;
      gtk_widget_hide(w);
   } else {
      std::cout<< "BL ERROR:: failed to get a plot" <<std::endl;
   }
}

extern "C" G_MODULE_EXPORT void
on_rama_export_as_pdf_filechooserdialog_response(GtkDialog *dialog, gint response_id, gpointer user_data) {

   if (response_id == GTK_RESPONSE_OK) {
      GtkWidget *canvas = GTK_WIDGET(user_data);
      coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas), "user_data"));
      if (plot) {
         std::string file_name =
               gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(plot->rama_export_as_pdf_filechooserdialog));
         plot->write_pdf(file_name);
      } else {
         std::cout<< "BL ERROR:: failed to get a plot" <<std::endl;
      }
   }
   gtk_widget_hide(GTK_WIDGET(dialog));
}

extern "C" G_MODULE_EXPORT void
on_rama_export_as_png_filechooserdialog_close(GtkDialog *dialog, gpointer user_data){

   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas),
                                                                             "user_data"));
   if (plot) {
      GtkWidget *w = plot->rama_export_as_png_filechooserdialog;
      gtk_widget_hide(w);
   } else {
      std::cout<< "BL ERROR:: failed to get a plot" <<std::endl;
   }
}

extern "C" G_MODULE_EXPORT void
on_rama_export_as_png_filechooserdialog_response(GtkDialog *dialog, gint response_id, gpointer user_data) {

   if (response_id == GTK_RESPONSE_OK) {
      GtkWidget *canvas = GTK_WIDGET(user_data);
      coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas),
                                                                                "user_data"));
      if (plot) {
         std::string file_name =
               gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(plot->rama_export_as_png_filechooserdialog));
         plot->write_png(file_name);
      } else {
         std::cout<< "BL ERROR:: failed to get a plot" <<std::endl;
      }
   }
   gtk_widget_hide(GTK_WIDGET(dialog));
}

extern "C" G_MODULE_EXPORT void
on_rama_open_filechooserdialog_close(GtkDialog *dialog, gpointer user_data){

   GtkWidget *canvas = GTK_WIDGET(user_data);
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas),
                                                                             "user_data"));
   if (plot) {
      GtkWidget *w = plot->rama_open_filechooserdialog;
      gtk_widget_hide(w);
   } else {
      std::cout<< "BL ERROR:: failed to get a plot" <<std::endl;
   }
}

extern "C" G_MODULE_EXPORT void
on_rama_open_filechooserdialog_response(GtkDialog *dialog, gint response_id, gpointer user_data) {

   if (response_id == GTK_RESPONSE_OK) {
      GtkWidget *canvas = GTK_WIDGET(user_data);
      coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(canvas),
                                                                                "user_data"));
      if (plot) {
         std::string file_name =
               gtk_file_chooser_get_filename(GTK_FILE_CHOOSER(plot->rama_open_filechooserdialog));
         plot->open_pdb_file(file_name);
      } else {
         std::cout<< "BL ERROR:: failed to get a plot" <<std::endl;
      }
   }
   gtk_widget_hide(GTK_WIDGET(dialog));
}

#ifdef HAVE_GOOCANVAS

// Canvas and item callbacks
gboolean rama_item_button_press (GooCanvasItem *item,
                                 GooCanvasItem *target,
                                 GdkEventButton *event,
                                 gpointer data) {

   gchar *id;
   id = (gchar*)g_object_get_data (G_OBJECT (item), "id");
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(item), "rama_plot"));

   //g_print ("BL DEBUG:: %s received button-press event\n", id ? id : "unknown");

   plot->button_item_press(item, event);
   return TRUE;
}
#endif

#ifdef HAVE_GOOCANVAS

gboolean rama_item_button_release (GooCanvasItem *item,
                      GooCanvasItem *target,
                      GdkEventButton *event,
                      gpointer data) {

   gchar *id;
   id = (gchar*)g_object_get_data (G_OBJECT (item), "id");
   coot::rama_plot *plot =
      static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(item), "rama_plot"));

   //g_print ("BL DEBUG:: %s received button-release event\n", id ? id : "unknown");

   plot->button_item_release(item, event);
   return TRUE;
}
#endif

#ifdef HAVE_GOOCANVAS
gboolean rama_item_enter_event (GooCanvasItem *item,
                                GooCanvasItem *target,
                                GdkEventCrossing *event,
                                gpointer data) {

   coot::rama_plot *plot =
      static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(item), "rama_plot"));

   plot->item_enter_event(item, event);

   // do something, maybe pass some data for usefullness.

   return TRUE;
}
#endif

#ifdef HAVE_GOOCANVAS
gboolean rama_item_motion_event (GooCanvasItem *item,
                                GooCanvasItem *target,
                                GdkEventMotion *event,
                                gpointer data) {


   void *p = g_object_get_data(G_OBJECT(item), "rama_plot");
   coot::rama_plot *plot = static_cast<coot::rama_plot *> (p);

   plot->item_motion_event(item, event);

   // do something, maybe pass some data for usefullness.

   return TRUE;

}
#endif

// The motion callback was attached at the canvas, so widget is a
// canvas here.
// 
gint rama_motion_notify(GtkWidget *widget, GdkEventMotion *event) {

   double x,y;
   int x_as_int, y_as_int;
   GdkModifierType state;

   if (event->is_hint) {
      // gdk_window_get_pointer(event->window, &x_as_int, &y_as_int, &state);
      GdkModifierType mask;
      GdkSeat *seat = gdk_display_get_default_seat(gdk_display_get_default());
      GdkDevice *mouse = gdk_seat_get_pointer(seat);
      gdk_window_get_device_position(event->window, mouse, &x_as_int, &y_as_int, &mask);
      x = x_as_int;
      y = y_as_int;
   } else {
      x = event->x;
      y = event->y;
      state = (GdkModifierType) event->state;
   }

   return 0;
}


gint rama_button_press (GtkWidget *widget, GdkEventButton *event) {

   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(widget), "rama_plot"));

   plot->button_press(widget, event);
   //g_print("BL DEBUG:: button press notify\n");

   return 0;
   
}

gint rama_key_release_event(GtkWidget *widget, GdkEventKey *event) {

   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(widget), "rama_plot"));

   gint i = plot->key_release_event(widget, event);
   //g_print("BL DEBUG:: key release event\n");
   return i;
}

gint rama_key_press_event(GtkWidget *widget, GdkEventKey *event) {

   // needed?
   //g_print("BL DEBUG:: key press event\n");
   return 0;

}

//FIXME - maybe should have at some point....
//void rama_show_preferences() {
 
//   //    GtkWidget *widget = create_propertybox1();
//   GtkWidget *widget = create_dynarama_properties_window();
//   gtk_widget_show(widget);

//}

void rama_zoom_out(GtkWidget *widget) {

   coot::rama_plot *plot = static_cast<coot::rama_plot *> (g_object_get_data(G_OBJECT(widget), "rama_plot"));

   plot->zoom_out(); 

}

void rama_zoom_in(GtkWidget *widget) {

   coot::rama_plot *plot = (coot::rama_plot *) g_object_get_data(G_OBJECT(widget), "rama_plot");

   plot->zoom_in(); 

}

gboolean rama_configure_event(GtkWidget *widget, GdkEventConfigure *event, gpointer user_data){

   coot::rama_plot *plot = static_cast<coot::rama_plot *>(user_data);
   // std::cout << "rama_resize() " << user_data << " " << plot << std::endl;
   plot->resize_rama_canvas_internal(widget, event);
   return FALSE;
}

#endif // HAVE_GTK_CANVAS or HAVE_GNOME_CANVAS
