/* src/c-interface-build.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 The University of York
 * Author: Paul Emsley
 * Copyright 2007 by Paul Emsley
 * Copyright 2007 by the University of Oxford
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

#include <string>
#include <vector>

#include <gtk/gtk.h>
#include "graphics-info.h"
// Including python needs to come after graphics-info.h, because
// something in Python.h (2.4 - chihiro) is redefining FF1 (in
// ssm_superpose.h) to be 0x00004000 (Grrr).
// BL says:: and (2.3 - dewinter), i.e. is a Mac - Python issue
// since the follwing two include python graphics-info.h is moved up
#include "c-interface.h"

void do_rot_trans_adjustments(GtkWidget *dialog) { 
   graphics_info_t g;
   g.do_rot_trans_adjustments(dialog);
}

short int delete_item_widget_is_being_shown() {
   short int r = 0; 
   if (graphics_info_t::delete_item_widget != NULL) {
      r = 1;
   }
   return r;
}

short int delete_item_widget_keep_active_on() {
   short int r = 0;
   if (delete_item_widget_is_being_shown()) { 
      GtkWidget *checkbutton = lookup_widget(graphics_info_t::delete_item_widget,
					     "delete_item_keep_active_checkbutton");
      if (GTK_TOGGLE_BUTTON(checkbutton)->active) {
	 r = 1;
      }
   }
   return r;
}

void store_delete_item_widget_position() {

   gint upositionx, upositiony;
   gdk_window_get_root_origin (graphics_info_t::delete_item_widget->window,
			       &upositionx, &upositiony);
   graphics_info_t::delete_item_widget_x_position = upositionx;
   graphics_info_t::delete_item_widget_y_position = upositiony;
   gtk_widget_destroy(graphics_info_t::delete_item_widget);
   clear_delete_item_widget();
}

void clear_delete_item_widget() {

   graphics_info_t::delete_item_widget = NULL;
}


void store_delete_item_widget(GtkWidget *widget) {
   graphics_info_t::delete_item_widget = widget;
}



/*  find the molecule that the single map dialog applies to and set
    the contour level and redraw */
void single_map_properties_apply_contour_level_to_map(GtkWidget *w) {

   int imol = GPOINTER_TO_INT(gtk_object_get_user_data(GTK_OBJECT(w)));

   if (is_valid_map_molecule(imol)) { 
      GtkToggleButton *toggle_button =
	 GTK_TOGGLE_BUTTON(lookup_widget(w, "single_map_properties_sigma_radiobutton"));

      GtkWidget *entry = lookup_widget(w, "single_map_properties_contour_level_entry");
      const char *txt = gtk_entry_get_text(GTK_ENTRY(entry));
      float level = atof(txt);
      if (toggle_button->active) {
	 set_contour_level_in_sigma(imol, level);
      } else {
	 set_contour_level_absolute(imol, level);
      }
   }
}

							 

/*! \brief a gui dialog showing remarks header info (for a model molecule). */
void remarks_dialog(int imol) { 

   if (graphics_info_t::use_graphics_interface_flag) { 
      if (is_valid_model_molecule(imol)) {
	 CMMDBManager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
	 if (mol) {
	    CTitleContainer *tc_p = mol->GetRemarks();
	    int l = tc_p->Length();
	    std::map<int, std::vector<std::string> > remarks;
	    for (unsigned int i=0; i<l; i++) { 
	       CRemark *cr = static_cast<CRemark *> (tc_p->GetContainerClass(i));
	       int rn = cr->remarkNum;
	       std::string s = cr->Remark;
	       remarks[rn].push_back(s);
	    }
	    if (! remarks.size()) {
	       info_dialog("No REMARKS");
	    } else { 
	       GtkWidget *d = gtk_dialog_new();
	       gtk_object_set_data(GTK_OBJECT(d), "remarks_dialog", d);
	       GtkWidget *vbox = GTK_DIALOG(d)->vbox;
	       GtkWidget *vbox_inner = gtk_vbox_new(FALSE, 2);
	       GtkWidget *scrolled_window = gtk_scrolled_window_new (NULL, NULL);
	       gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrolled_window),
						     GTK_WIDGET(vbox_inner));
	       gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(scrolled_window), TRUE, TRUE, 2);
	       gtk_widget_show(scrolled_window);
	       gtk_widget_show(vbox_inner);
	       
	       std::map<int, std::vector<std::string> >::const_iterator it;
	       for (it=remarks.begin(); it != remarks.end(); it++) {
		  std::string remark_name = "REMARK ";
		  remark_name += coot::util::int_to_string(it->first);
		  GtkWidget *frame = gtk_frame_new(remark_name.c_str());
		  gtk_box_pack_start(GTK_BOX(vbox_inner), frame, FALSE, FALSE, 1);
		  gtk_widget_show(frame);
		  // std::cout << "REMARK number " << it->first << std::endl;
		  GtkTextBuffer *text_buffer = gtk_text_buffer_new(NULL);
		  GtkWidget *text_view = gtk_text_view_new();
		  gtk_text_view_set_border_window_size(GTK_TEXT_VIEW(text_view),
						       GTK_TEXT_WINDOW_RIGHT, 10);
		  gtk_widget_set_usize(GTK_WIDGET(text_view), 400, -1);
		  gtk_container_add(GTK_CONTAINER(frame), GTK_WIDGET(text_view));
		  gtk_widget_show(GTK_WIDGET(text_view));
		  gtk_text_view_set_buffer(GTK_TEXT_VIEW(text_view), text_buffer);
		  gtk_text_view_set_wrap_mode(GTK_TEXT_VIEW(text_view), GTK_WRAP_WORD);

		  GdkColor colour = remark_number_to_colour(it->first); 
		  gtk_widget_modify_base(GTK_WIDGET(text_view), GTK_STATE_NORMAL, &colour);
		  
		  GtkTextIter end_iter;
		  for (unsigned int itext=0; itext<it->second.size(); itext++) { 
		     gtk_text_buffer_get_end_iter(text_buffer, &end_iter);
		     std::string s = it->second[itext];
		     s += "\n";
		     gtk_text_buffer_insert(text_buffer, &end_iter, s.c_str(), -1);
		  }
	       }


	       GtkWidget *close_button = gtk_button_new_with_label("  Close   ");
	       GtkWidget *aa = GTK_DIALOG(d)->action_area;
	       gtk_box_pack_start(GTK_BOX(aa), close_button, FALSE, FALSE, 2);
	       
	       gtk_signal_connect(GTK_OBJECT(close_button), "clicked",
				  GTK_SIGNAL_FUNC(on_remarks_dialog_close_button_clicked), NULL);
	       gtk_widget_show(close_button);
	       gtk_widget_set_usize(d, 500, 400);
	       gtk_widget_show(d);
	    } 
	 } 
      }
   }
} 


void
on_remarks_dialog_close_button_clicked     (GtkButton *button,
					    gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button), "remarks_dialog");
   gtk_widget_destroy(window);
}


GdkColor remark_number_to_colour(int remark_number) {

   GdkColor colour;
   colour.red   = 65535;
   colour.green = 65535;
   colour.blue  = 65535; 
   colour.pixel = 65535; 
   if (remark_number == 2) { 
      colour.blue   = 60000;
   }
   if (remark_number == 3) { 
      colour.red   = 60000;
   }
   if (remark_number == 4) { 
      colour.green   = 60000;
   }
   if (remark_number == 5) { 
      colour.green   = 62000;
      colour.blue    = 62000;
   }
   if (remark_number == 280) { 
      colour.green   = 61000;
      colour.red    = 62500;
   }
   if (remark_number == 465) { 
      colour.blue    = 60000;
      colour.green   = 60000;
   }
   return colour;
} 


void on_simple_text_dialog_close_button_pressed( GtkWidget *button,
						 GtkWidget *dialog) {
   gtk_widget_destroy(dialog);
}


void simple_text_dialog(const std::string &dialog_title, const std::string &text,
			std::pair<int, int> geom) {

   if (graphics_info_t::use_graphics_interface_flag) {

      GtkWidget *d = gtk_dialog_new();
      gtk_object_set_data(GTK_OBJECT(d), "simple_text_dialog", d);
      gtk_window_set_title (GTK_WINDOW (d), _(dialog_title.c_str()));
      GtkWidget *vbox = GTK_DIALOG(d)->vbox;
      GtkWidget *vbox_inner = gtk_vbox_new(FALSE, 2);
      GtkWidget *scrolled_window = gtk_scrolled_window_new (NULL, NULL);
      gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrolled_window),
					    GTK_WIDGET(vbox_inner));
      gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(scrolled_window), TRUE, TRUE, 2);
      gtk_widget_show(scrolled_window);
      gtk_widget_show(vbox_inner);
      
      GtkWidget *text_widget = gtk_text_view_new ();
      gtk_widget_show (text_widget);
      gtk_container_add (GTK_CONTAINER (vbox_inner), text_widget);
      gtk_text_view_set_editable (GTK_TEXT_VIEW (text_widget), FALSE);
      gtk_text_view_set_wrap_mode (GTK_TEXT_VIEW (text_widget), GTK_WRAP_WORD);
      gtk_text_buffer_set_text (gtk_text_view_get_buffer(GTK_TEXT_VIEW (text_widget)),
				text.c_str(), -1);
      gtk_window_set_default_size(GTK_WINDOW(d), geom.first, geom.second);

      GtkWidget *close_button = gtk_dialog_add_button(GTK_DIALOG(d), "Close", 2);
      gtk_widget_show(close_button);

       g_signal_connect(G_OBJECT(close_button), "clicked",
 		       G_CALLBACK(on_simple_text_dialog_close_button_pressed),
 		       (gpointer) d);
      
      gtk_widget_show(d);

   }
}
