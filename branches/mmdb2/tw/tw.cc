/* src/tw.c
 * 
 * Copyright 2004 by The University of York
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
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */

#include <iostream>
#include "coot-tw.hh"

#if defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)

#ifdef HAVE_GNOME_CANVAS
  #define gtk_canvas_init gnome_canvas_init
  #define gtk_canvas_new  gnome_canvas_new
  #define gtk_canvas_root gnome_canvas_root
  #define gtk_canvas_item_new gnome_canvas_item_new
  #define gtk_canvas_item_w2i gnome_canvas_item_w2i
  #define gtk_canvas_item_grab gnome_canvas_item_grab
  #define gtk_canvas_item_lower_to_bottom gnome_canvas_item_lower_to_bottom
  #define gtk_canvas_item_lower gnome_canvas_item_lower
  #define gtk_canvas_item_raise_to_top gnome_canvas_item_raise_to_top
  #define gtk_canvas_item_raise gnome_canvas_item_raise
  #define gtk_canvas_item_move gnome_canvas_item_move
  #define gtk_canvas_item_ungrab gnome_canvas_item_ungrab
  #define gtk_canvas_rect_get_type gnome_canvas_rect_get_type
  #define GTK_CANVAS      GNOME_CANVAS
#endif

GtkCanvas *
coot::tw::create_canvas() {

   int usize_x = 650;
   int usize_y = 350;

   // From rama_plot:
   // gtk_widget_set_usize(GTK_WIDGET(canvas), 400, 400);
   // gtk_widget_set_usize(app1, 400, 500);

#ifndef HAVE_GNOME_CANVAS
   gdk_imlib_init();
   gtk_canvas_init(); 
#endif

   /* Get gdk to use imlib's visual and colormap */
   gtk_widget_push_visual(gdk_imlib_get_visual());
#ifndef HAVE_GNOME_CANVAS
   gtk_widget_push_colormap(gdk_imlib_get_colormap());
#endif
   canvas = GTK_CANVAS(gtk_canvas_new());
   gtk_widget_pop_colormap();
   gtk_widget_pop_visual();
   
   gtk_widget_set_usize(GTK_WIDGET(canvas), usize_x, usize_y);
   gtk_widget_show(GTK_WIDGET(canvas));
   gtk_object_set_user_data(GTK_OBJECT(canvas), (void *) this);

   gtk_widget_set_events(GTK_WIDGET(canvas),
			 GDK_EXPOSURE_MASK      |
			 GDK_BUTTON_PRESS_MASK  |
			 GDK_BUTTON_RELEASE_MASK|
			 GDK_POINTER_MOTION_MASK|
			 GDK_KEY_RELEASE_MASK   |
			 GDK_POINTER_MOTION_HINT_MASK);


   return canvas;

};


GtkWidget *
coot::tw::create_filled_toplevel() {

   GtkWidget *toplevel = create_tw_dialog();
   GtkWidget *scrolled_window = tw_lookup_widget(toplevel, "tw_scrolledwindow");
   canvas = create_canvas();
   gtk_container_add(GTK_CONTAINER(scrolled_window), GTK_WIDGET(canvas));
   gtk_widget_show(GTK_WIDGET(canvas));
   add_elements_to_canvas();

   /* mouse in motion! */
//    gtk_signal_connect (GTK_OBJECT(canvas), "event",
// 		       GTK_SIGNAL_FUNC(tw_motion_notify), NULL);

   gtk_widget_show(toplevel);
   return toplevel;

};

// static 
gint
coot::tw::tw_motion_notify (GtkWidget *widget, GdkEventMotion *event) {

   std::cout << "event on canvas! " << std::endl;
   return TRUE;
}



GtkWidget*
coot::tw::create_tw_dialog (void)
{
  GtkWidget *tw_dialog;
  GtkWidget *dialog_vbox1;
  GtkWidget *tw_scrolledwindow;
  GtkWidget *dialog_action_area1;
  GtkWidget *hbox1;
  GtkWidget *tw_close_button;

  tw_dialog = gtk_dialog_new ();
  gtk_object_set_data (GTK_OBJECT (tw_dialog), "tw_dialog", tw_dialog);
  gtk_window_set_title (GTK_WINDOW (tw_dialog), "Cooters");
  gtk_window_set_policy (GTK_WINDOW (tw_dialog), TRUE, TRUE, FALSE);

  dialog_vbox1 = GTK_DIALOG (tw_dialog)->vbox;
  gtk_object_set_data (GTK_OBJECT (tw_dialog), "dialog_vbox1", dialog_vbox1);
  gtk_widget_show (dialog_vbox1);

  tw_scrolledwindow = gtk_scrolled_window_new (0, 0);
  gtk_widget_ref (tw_scrolledwindow);
  gtk_object_set_data_full (GTK_OBJECT (tw_dialog), "tw_scrolledwindow", tw_scrolledwindow,
                            (GtkDestroyNotify) gtk_widget_unref);
  gtk_widget_show (tw_scrolledwindow);
  gtk_box_pack_start (GTK_BOX (dialog_vbox1), tw_scrolledwindow, TRUE, TRUE, 0);

  dialog_action_area1 = GTK_DIALOG (tw_dialog)->action_area;
  gtk_object_set_data (GTK_OBJECT (tw_dialog), "dialog_action_area1", dialog_action_area1);
  gtk_widget_show (dialog_action_area1);
  gtk_container_set_border_width (GTK_CONTAINER (dialog_action_area1), 10);

  hbox1 = gtk_hbox_new (FALSE, 0);
  gtk_widget_ref (hbox1);
  gtk_object_set_data_full (GTK_OBJECT (tw_dialog), "hbox1", hbox1,
                            (GtkDestroyNotify) gtk_widget_unref);
  gtk_widget_show (hbox1);
  gtk_box_pack_start (GTK_BOX (dialog_action_area1), hbox1, TRUE, TRUE, 0);

  tw_close_button = gtk_button_new_with_label ("  Close ");
  gtk_widget_ref (tw_close_button);
  gtk_object_set_data_full (GTK_OBJECT (tw_dialog), "tw_close_button", tw_close_button,
                            (GtkDestroyNotify) gtk_widget_unref);

  gtk_signal_connect (GTK_OBJECT(tw_close_button), "clicked",
		      GTK_SIGNAL_FUNC(on_tw_close_button_clicked),
		      NULL);

  gtk_widget_show (tw_close_button);
  gtk_box_pack_start (GTK_BOX (hbox1), tw_close_button, TRUE, TRUE, 0);


//   float v = 0.4;
//   for (int i=0; i<2; i++) {
//      if (i==1) v = 0.7;
//      GtkWidget *hscale = gtk_hscale_new (GTK_ADJUSTMENT (gtk_adjustment_new (0, 0, 0, 0, 0, 0)));
//      // 0.0 to 0.2 steps 0.1 start 0.4
//      // GtkAdjustment *adj = GTK_ADJUSTMENT(gtk_adjustment_new(v, 0.0, 1.6, 0.01, 0.2, 1.4));
//      // 0.0 to 1.0 steps 0.1 start 0.4 [hooray]
//      // GtkAdjustment *adj = GTK_ADJUSTMENT(gtk_adjustment_new(v, 0.0, 1.6, 0.01, 0.2, 0.6));
//      // steps in size 0.1
//      // GtkAdjustment *adj = GTK_ADJUSTMENT(gtk_adjustment_new(v, 0.0, 1.0, 0.001, 0.01, 0.0));

//      // The max value is 3rd arg - 6th arg
//      GtkAdjustment *adj = GTK_ADJUSTMENT(gtk_adjustment_new(v, 0.0, 2.01, 0.01, 0.1, 1.01));
//      gtk_range_set_adjustment(GTK_RANGE(hscale), GTK_ADJUSTMENT(adj));
//      gtk_box_pack_start(GTK_BOX(hbox1), hscale, TRUE, TRUE, 0);
//      gtk_widget_show(hscale);
//   }


  return tw_dialog;
}


void
on_tw_close_button_clicked     (GtkButton *button,
				gpointer         user_data)
{
   GtkWidget *w = tw_lookup_widget(GTK_WIDGET(button),
				   "tw_dialog");
   // gtk_widget_destroy(w);
   gtk_exit(0);
}

GtkWidget*
tw_lookup_widget                          (GtkWidget       *widget,
					   const gchar     *widget_name)
{
  GtkWidget *parent, *found_widget;

  for (;;)
    {
      if (GTK_IS_MENU (widget))
        parent = gtk_menu_get_attach_widget (GTK_MENU (widget));
      else
        parent = widget->parent;
      if (parent == NULL)
        break;
      widget = parent;
    }

  found_widget = (GtkWidget*) gtk_object_get_data (GTK_OBJECT (widget),
                                                   widget_name);
  if (!found_widget)
    g_warning ("Widget not found: %s", widget_name);
  return found_widget;
}


void
coot::tw::add_elements_to_canvas() {

   // put some thingies in canvas
   //
   // First let's draw a box
   GtkCanvasItem *item;
   double d = 150;

   item = gtk_canvas_item_new(gtk_canvas_root(canvas),
			      // GTK_CANVAS_TYPE_CANVAS_RECT,
			      gtk_canvas_rect_get_type(),
			      "x1", -250.0,
			      "y1", -70.0,
			      "x2", 300.0,
			      "y2", 220.0,
			      "fill_color", "grey60",
			      "outline_color", "black",
			      NULL);
   
   canvas_item_vec.push_back(item);

   std::string coot_pixmap_dir = PKGDATADIR;  
   coot_pixmap_dir += "/";
   coot_pixmap_dir += "pixmaps";
   
   // an image
   std::string image_path;
   image_path = coot_pixmap_dir + "/eugene.jpg";
   add_image_to_canvas(image_path, 0);
   image_path = coot_pixmap_dir + "/kevin.jpg";
   add_image_to_canvas(image_path, 1);
   image_path = coot_pixmap_dir + "/emsley.jpg";
   add_image_to_canvas(image_path, 2);

}

int
coot::tw::add_image_to_canvas(const std::string &filename,
			      int position_index) {

   int istat = 0;
   // ArtPixBuf *aa_image; // for aa-mode
#ifndef HAVE_GNOME_CANVAS
   GdkImlibImage *im = gdk_imlib_load_image((char *)filename.c_str());
#else
   GdkPixbuf *im = gdk_pixbuf_new_from_file((char *)filename.c_str(), NULL);
#endif
   GtkAnchorType anchor = GTK_ANCHOR_NW;
   GtkCanvasItem *item; 

   if (im) { 


#ifdef HAVE_GTK_CANVAS
      item = gtk_canvas_item_new(gtk_canvas_root(canvas),
				 gtk_canvas_image_get_type(),
				 "image", im, 
				 "x", -220.0 + 200.0 * position_index,
				 "y", 0.0,
				 "width", (double) im->rgb_width/2.0,
				 "height", (double) im->rgb_height/2.0,
				 "anchor", anchor,
				 NULL);
#else
      item = gtk_canvas_item_new(gnome_canvas_root(canvas),
				 gnome_canvas_pixbuf_get_type(),
				 "pixbuf", im, 
				 "x", -220.0 + 200.0 * position_index,
				 "y", 0.0,
				 "width", (double) gdk_pixbuf_get_width(im)/2.0,
				 "height", (double) gdk_pixbuf_get_height(im)/2.0,
				 "anchor", anchor,
				 NULL);
#endif

      gtk_signal_connect (GTK_OBJECT (item), "destroy",
			  (GtkSignalFunc) free_imlib_image,
			  im); 

      gtk_signal_connect(GTK_OBJECT(item), "event", 
			 (GtkSignalFunc) on_tw_canvas_item_clicked, 
			 gpointer(item));
   
      canvas_item_vec.push_back(item);
      istat = 1;
   } else {
      std::cout << "No image " << filename << std::endl;
   }
   return istat;
}

// static
void
coot::tw::free_imlib_image (GtkObject *object, gpointer data)
{
#ifdef HAVE_GTK_CANVAS
   gdk_imlib_destroy_image ((GdkImlibImage *)data);
#else
   gdk_pixbuf_unref ((GdkPixbuf *)data);
#endif
}


// static
gint
coot::tw::on_tw_canvas_item_clicked(GtkCanvasItem *item, GdkEvent *event, gpointer data) {

   if (event->motion.state & GDK_BUTTON1_MASK) {
      // std::cout << "click!" << std::endl;
   } else {
      // std::cout << "something else!" << std::endl;
   }

   static double x, y;
   double new_x, new_y;
   GdkCursor *fleur;
   static int dragging;
   double item_x, item_y;

   /* set item_[xy] to the event x,y position in the parent's item-relative coordinates */
   item_x = event->button.x;
   item_y = event->button.y;
   gtk_canvas_item_w2i (item->parent, &item_x, &item_y);

   switch (event->type) {
   case GDK_BUTTON_PRESS:
      switch (event->button.button) {
      case 1:
	 if (event->button.state & GDK_SHIFT_MASK)
	    gtk_object_destroy (GTK_OBJECT (item));
	 else {
	    x = item_x;
	    y = item_y;

	    fleur = gdk_cursor_new (GDK_FLEUR);
	    gtk_canvas_item_grab (item,
				  GDK_POINTER_MOTION_MASK | GDK_BUTTON_RELEASE_MASK,
				  fleur,
				  event->button.time);
	    gdk_cursor_destroy (fleur);
	    dragging = TRUE;
	 }
	 break;

      case 2:
	 if (event->button.state & GDK_SHIFT_MASK)
	    gtk_canvas_item_lower_to_bottom (item);
	 else
	    gtk_canvas_item_lower (item, 1);
	 break;

      case 3:
	 if (event->button.state & GDK_SHIFT_MASK)
	    gtk_canvas_item_raise_to_top (item);
	 else
	    gtk_canvas_item_raise (item, 1);
	 break;

      default:
	 break;
      }

      break;

   case GDK_MOTION_NOTIFY:
      if (dragging && (event->motion.state & GDK_BUTTON1_MASK)) {
	 new_x = item_x;
	 new_y = item_y;

	 gtk_canvas_item_move (item, new_x - x, new_y - y);
	 x = new_x;
	 y = new_y;
      }
      break;

   case GDK_BUTTON_RELEASE:
      gtk_canvas_item_ungrab (item, event->button.time);
      dragging = FALSE;
      break;

   default:
      break;
   }

   return FALSE;

}

#endif //  defined(HAVE_GNOME_CANVAS) || defined(HAVE_GTK_CANVAS)
