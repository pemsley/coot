/* src/ncs.cc
 * 
 * Copyright 2009 by The University of Oxford
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


#if defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)

// Don't forget to enable
// g.set_sequence_view_is_displayed(seq_view->Canvas(), imol) in nsv()
// when this is working.

#include <vector>
#include <iostream>
#include <algorithm>
#include <map>

#include "nsv.hh"
#include "support.h"
#include "utils/coot-utils.hh"

#include "c-interface.h" // for sequence_view_is_displayed()
#include "cc-interface.hh" // for set_go_to_atom_from_spec()

#include "get-residue.hh"

#include "canvas-fixes.hh"

exptl::nsv::nsv(CMMDBManager *mol,
		const std::string &molecule_name,
		int molecule_number_in,
		bool use_graphics_interface_in) {

   molecule_number = molecule_number_in;
   use_graphics_interface_flag = use_graphics_interface_in;

   GtkWidget *top_lev = gtk_dialog_new();
   gtk_object_set_data(GTK_OBJECT(top_lev), "nsv_dialog", top_lev);
   GtkWidget *vbox = GTK_DIALOG(top_lev)->vbox;
   canvas = GTK_CANVAS(gtk_canvas_new());
   GtkWidget *name_label = gtk_label_new(molecule_name.c_str());
   gtk_box_pack_start(GTK_BOX(vbox), name_label, FALSE, FALSE, 1);
   GtkWidget *scrolled_window = gtk_scrolled_window_new(NULL, NULL);
   gtk_box_pack_start(GTK_BOX(vbox), GTK_WIDGET(scrolled_window), TRUE, TRUE, 1);
   gtk_scrolled_window_add_with_viewport(GTK_SCROLLED_WINDOW(scrolled_window),
					 GTK_WIDGET(canvas));
   
   GtkWidget *close_button = gtk_button_new_with_label("  Close   ");
   GtkWidget *aa = GTK_DIALOG(top_lev)->action_area;
   gtk_box_pack_start(GTK_BOX(aa), close_button, FALSE, FALSE, 2);

   gtk_signal_connect(GTK_OBJECT(close_button), "clicked",
		      GTK_SIGNAL_FUNC(on_nsv_close_button_clicked), NULL);

   gtk_signal_connect(GTK_OBJECT(top_lev), "destroy",
		      GTK_SIGNAL_FUNC(on_nsv_dialog_destroy), top_lev);

   // used on destroy
   gtk_object_set_user_data(GTK_OBJECT(top_lev),
			    GINT_TO_POINTER(molecule_number));

   setup_canvas(mol, scrolled_window);
   gtk_widget_show(scrolled_window);
		      
   if (use_graphics_interface_flag) { 
      gtk_widget_show(aa);
      gtk_widget_show(close_button);
      gtk_widget_show(GTK_WIDGET(name_label));
      gtk_widget_show(GTK_WIDGET(canvas));
      gtk_widget_show(top_lev);
   }

   // set_sequence_view_is_displayed(canvas, molecule_number) is run
   // by nsv() (currently that's its name)
}

// static 
void
exptl::nsv::on_nsv_close_button_clicked     (GtkButton *button,
					     gpointer         user_data)
{
   GtkWidget *window = lookup_widget(GTK_WIDGET(button),
				     "nsv_dialog");
   gtk_widget_destroy(window);
}

// static
void
exptl::nsv::on_nsv_dialog_destroy (GtkObject *obj,
				   gpointer user_data) {

   GtkWidget *dialog = (GtkWidget *) user_data;
   int imol = GPOINTER_TO_INT(gtk_object_get_user_data(GTK_OBJECT(dialog)));
   set_sequence_view_is_displayed(0, imol);
} 



void
exptl::nsv::setup_canvas(CMMDBManager *mol, GtkWidget *scrolled_window) {

#if defined(WINDOWS_MINGW) || defined(_MSC_VER)
   fixed_font_str = "monospace";
#else
   fixed_font_str = "fixed";
   fixed_font_str = "Sans 9";
#endif
   pixels_per_letter = 10; // 10 for my F10 box
   pixels_per_chain  = 12;

   bool debug = true;

#ifdef HAVE_GTK_CANVAS
   gtk_canvas_init();
#endif   

   std::vector<chain_length_residue_units_t> rcv = get_residue_counts(mol);
   if (rcv.size() > 0) {
      // if the chains start within 20 residues of the chain with the
      // lowest first residue number, then make all the chains start
      // at the same residue number.  They should be in the same graph
      // group.  If not, then a different graph group.
      //
      //
      int lowest_resno = rcv[0].first_res_no;
      int biggest_res_count = 0;
      int biggest_res_number = 0;
      for (unsigned int i=0; i<rcv.size(); i++) {
	 if (rcv[i].first_res_no < lowest_resno)
	    lowest_resno = rcv[i].first_res_no;
	 if (rcv[i].n_residue_count > biggest_res_count)
	    biggest_res_count = rcv[i].n_residue_count;
	 if (rcv[i].max_resno > biggest_res_number)
	    biggest_res_number = rcv[i].max_resno;
      }
      int total_res_range = biggest_res_number - lowest_resno;

      if (debug) {
	 std::cout << "DEBUG:: resnos: biggest_res_number "
		   << biggest_res_number << " lowest_resno: "
		   << lowest_resno << " total_res_range: " << total_res_range << std::endl;
      }
      
      std::vector<CResidue*> ins_code_residues =
	 coot::util::residues_with_insertion_codes(mol);


      // If there are no insertion codes, let's go.  Work directly
      // from the residue number
      
      if (ins_code_residues.size() == 0) { 

	 int n_chains = rcv.size();
	 int n_limited_chains = n_chains;
	 if (n_limited_chains > 30)
	    n_limited_chains = 30;
	 // e.g. chain_id_graph_map["B"] -> (graph_group_no, sequence_line_no)
	 // 
	 std::map <std::string, std::pair<int, int> > chain_id_graph_map;

	 // 100 is good for 2 chains
	 // 

	 int canvas_x_size =  5 + total_res_range * pixels_per_letter + 140;
	 // 50 is good for 2 chains.
	 // int canvas_y_size =  5 + n_limited_chains * 50; 
	 int canvas_y_size = 65 + n_limited_chains * 20; 

	 // the size of the widget on the screens
	 gtk_widget_set_usize(GTK_WIDGET(scrolled_window), 700, 20*(3+n_limited_chains));
	 // the size of the canvas (e.g. long chain, we see only part
	 // of it at one time).

	 if (debug) { 
	    std::cout << "DEBUG:: in setup_canvas(), total_res_range: " << total_res_range
		      << " canvas_x_size " << canvas_x_size << " and canvas_y_size "
		      << canvas_y_size << std::endl;
	    std::cout << "DEBUG:: n_limited_chains: " << n_limited_chains << std::endl;
	 }
	 gtk_widget_set_usize(GTK_WIDGET(canvas), canvas_x_size, canvas_y_size);

	 double left_limit = 0.0;
	 double upper_limit = 0.0;
	 double scroll_width;   // right 
	 double scroll_height;  // lower
	 double scroll_width_max = 32700;  // not sure this is needed
					   // (I guess that it is
					   // though)

	 left_limit =    biggest_res_count * 10;
	 upper_limit =   0.0;

	 // scroll_height affects the position of the canvas origin
	 // relative to the widget. scroll_width likewise affects the
	 // position of the origin on the widget.  Taking off a few
	 // pixels (currently 110) moves the graph rightwards on the
	 // page, a bit more esthetically pleasing, and we will need
	 // room to write Chain "A" etc on the left there too.
	 //

	 // for scroll_height, -100 is good for 2 chains, -150 is good
	 // for 4 chains.

	 // the casting of canvas_y_size is needed!  Without it
	 // scroll_height is 4.3e9.  Why!? canvas_y_size is an int.
	 // 
	 scroll_height = double(canvas_y_size) - 50 - 25 * rcv.size();
	 scroll_width  = double(canvas_x_size) - 130 - left_limit ; // -130 for fat font on jackal


	 if (debug) { 
	    std::cout << "DEBUG:: rcv.size(): " << rcv.size() << std::endl;
	    std::cout << "DEBUG:: canvas_y_size: " << canvas_y_size << std::endl;
	    std::cout << "DEBUG:: scroll_height: " << scroll_height << std::endl;
	    std::cout << "DEBUG:: gtk_canvas_set_scroll_region "
		      << left_limit << " " << upper_limit << " "
		      << scroll_width << " " << scroll_height << std::endl;
	 }
	 
	 if (scroll_width > scroll_width_max) {
	    // std::cout << "setting scroll width max to " << scroll_width_max << std::endl;
	    scroll_width= scroll_width_max;
	 }

	 // bring the items on the canvas leftwards, no empty big
	 // space on the left before sequences are shown.
	 //
	 // The more large negative, the more the sequence items are
	 // render on the canvas the the rightwards direction.
	 //
	 // double x_offset = -7 * pixels_per_letter + (biggest_res_number - lowest_resno) * 0.65;

	 // double x_offset = - total_res_range/15.0 * pixels_per_letter + total_res_range * 0.65;

	 // if you want to change this multiplier, test vs 3GP.pdb,
	 // 3mdo, tut-modern, 2vtu, 3g5u, others?
	 // 
	 double x_offset = - 0.0666 * total_res_range * pixels_per_letter + total_res_range * 0.65;

	 if (debug)
	    std::cout << " -> x_offset:" << x_offset << std::endl;

 	 gtk_canvas_set_scroll_region(canvas, left_limit, upper_limit,
 				      scroll_width, scroll_height);
	 origin_marker();
	 draw_axes(rcv, lowest_resno, biggest_res_number, x_offset);
	 mol_to_canvas(mol, lowest_resno, x_offset);
	 

      } else {
	 std::cout << "Some residues have insertion codes, this is not coded for"
		   << " yet\n";
      }
   }
}

void
exptl::nsv::mol_to_canvas(CMMDBManager *mol, int lowest_resno, double x_offset) {

   int imod = 1;
   CModel *model_p = mol->GetModel(imod);
   int ss_status = model_p->CalcSecStructure(1);
   CChain *chain_p;
   // run over chains of the existing mol
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int position_number = nchains - ichain - 1;
      chain_to_canvas(chain_p, position_number, lowest_resno, x_offset);
   }
}

void
exptl::nsv::chain_to_canvas(CChain *chain_p, int position_number, int lowest_resno, double x_offset) {

   int nres = chain_p->GetNumberOfResidues();
   GtkCanvasItem *item;

   for (int ires=0; ires<nres; ires++) {
      CResidue *residue_p = chain_p->GetResidue(ires);
      add_text_and_rect(residue_p, position_number, lowest_resno, x_offset);  // adds items to canvas_item_vec
   }

   // now the chain label:
   std::string chain_label = "Chain ";
   chain_label += chain_p->GetChainID();
   double x = -50;
   // On jackal it needs to be more left? Why? The font is wider?
   // How do I check the font width?
   x -= 10 + x_offset;
   double y = -pixels_per_chain * position_number - 6;
   item = gtk_canvas_item_new(gtk_canvas_root(canvas),
			      GTK_CANVAS_TYPE_CANVAS_TEXT,
			      "text", chain_label.c_str(),
			      "x", x,
			      "y", y,
			      "anchor", GTK_ANCHOR_WEST,
			      "font", fixed_font_str.c_str(),
			      NULL);
   canvas_item_vec.push_back(item);
}


bool
exptl::nsv::add_text_and_rect(CResidue *residue_p,
			      int position_number,
			      int lowest_resno,
			      double x_offset) {

   bool too_wide = false;
   if (residue_p) { 
      CAtom *at = coot::util::intelligent_this_residue_mmdb_atom(residue_p);
      coot::atom_spec_t at_spec(at);
      std::string res_code =
	 coot::util::three_letter_to_one_letter_with_specials(residue_p->GetResName());
      std::string colour = "black";
      double x = (residue_p->GetSeqNum() - lowest_resno + 1) * pixels_per_letter - x_offset;
      double y = - pixels_per_chain * position_number - 6;
      
      exptl::nsv::spec_and_object *so = 
	 new exptl::nsv::spec_and_object(molecule_number, at_spec, position_number);
      
      // It seems that I don't need the following, the letter seems to
      // be being clicked even if we click off the actual black of the
      // letter.
      // 
      // Hmm... maybe I do.
      // 
      // double x1 = x - 2;
      double x1 = x -3; // otherwise rect too far to right relative to letter
                         // (on fed10 home).
      double y1 = y + 5;
      
      // double x2 = x1 + 7; // pixels_per_letter; 8 is too many for
                          // jackal, it partially covers the
                          // previous letter.

      double x2 = x + pixels_per_letter - 2 - 3; // how about that then?
      
      double y2 = y1 - 11; // pixels_per_letter;
      std::string rect_colour = "grey85";

      double x1_max = 22500; // magic number
      if (x1 < x1_max) { 
	 GtkCanvasItem *rect_item;
	 rect_item = gtk_canvas_item_new(gtk_canvas_root(canvas),
					 GTK_CANVAS_TYPE_CANVAS_RECT,
					 "x1", x1,
					 "y1", y1,
					 "x2", x2,
					 "y2", y2,
					 "fill_color", rect_colour.c_str(),
					 "width_pixels", 1,
					 NULL);
	 so->obj = rect_item;
	 canvas_item_vec.push_back(rect_item);
	 gtk_signal_connect(GTK_OBJECT(rect_item), "event",
			    GTK_SIGNAL_FUNC(rect_event), so);

	 GtkCanvasItem *text_item;
	 std::string colour = colour_by_secstr(residue_p);
	 text_item = gtk_canvas_item_new(gtk_canvas_root(canvas),
					 GTK_CANVAS_TYPE_CANVAS_TEXT,
					 "text", res_code.c_str(),
					 "x", x,
					 "y", y,
					 "anchor", GTK_ANCHOR_CENTER,
					 "fill_color", colour.c_str(),
					 "font", fixed_font_str.c_str(),
					 NULL);
	 gtk_signal_connect(GTK_OBJECT(text_item), "event",
			    GTK_SIGNAL_FUNC(letter_event), so);
	 canvas_item_vec.push_back(text_item);
      } else {
	 too_wide = true;
      } 
   }
   return too_wide;
}



// static
gint
exptl::nsv::letter_event (GtkObject *obj,
			    GdkEvent *event,
			    gpointer data) {

   exptl::nsv::spec_and_object spec_obj = *((exptl::nsv::spec_and_object *) data);

   // I had been checking on mouse motion (it seems).  But lots of
   // mouse motion causes many calls to set_go_to_atom_from_spec()
   // (which is expensive).
   // So, only go to an atom on button release.
   // 

   // std::cout << "event type: " << event->type << std::endl;
   
   // if (event->motion.state & GDK_BUTTON1_MASK) {
   if (event->type == GDK_BUTTON_RELEASE) {
      set_go_to_atom_molecule(spec_obj.mol_no);
      set_go_to_atom_from_spec(spec_obj.atom_spec);
   } else {
      if (event->type == GDK_ENTER_NOTIFY) {
	 // std::cout << "Entering a text " << obj << " " << obj << std::endl;
	 gtk_canvas_item_set(GTK_CANVAS_ITEM(spec_obj.obj), "fill_color", "moccasin", NULL);
      }

      if (event->type == GDK_LEAVE_NOTIFY) {
	 // std::cout << " Leaving a text " << obj << " " << obj << std::endl;
	 gtk_canvas_item_set(GTK_CANVAS_ITEM(spec_obj.obj), "fill_color", "grey85", NULL);
      }
   }
   return 1;
}
   
// static
gint
exptl::nsv::rect_event (GtkObject *obj,
			GdkEvent *event,
			gpointer data) {
   
   exptl::nsv::spec_and_object spec_obj = *((exptl::nsv::spec_and_object *) data);
   // if (event->motion.state & GDK_BUTTON1_MASK) {
   if (event->type == GDK_BUTTON_RELEASE) {
      // std::cout << "box clicked" << std::endl;
      set_go_to_atom_molecule(spec_obj.mol_no);
      set_go_to_atom_from_spec(spec_obj.atom_spec);
      
   } else {

      if (event->type == GDK_ENTER_NOTIFY) {
	 // std::cout << "Entering a box " << obj << std::endl;
	 gtk_canvas_item_set(GTK_CANVAS_ITEM(obj), "fill_color", "moccasin", NULL);
      }

      if (event->type == GDK_LEAVE_NOTIFY) {
	 // std::cout << " Leaving a box " << obj << std::endl;
	 gtk_canvas_item_set(GTK_CANVAS_ITEM(obj), "fill_color", "grey85", NULL);
      } 
   } 
   return 1;
}


std::vector<exptl::nsv::chain_length_residue_units_t> 
exptl::nsv::get_residue_counts(CMMDBManager *mol) const {

   std::vector<exptl::nsv::chain_length_residue_units_t> chain_length_residue_units;
   int imod = 1;
   CModel *model_p = mol->GetModel(imod);
   CChain *chain_p;
   // run over chains of the existing mol
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      int n_residue_count = 0; // counting residues, resiudes with
			       // insertion codes and gaps.
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      PCResidue residue_p = 0;
      CAtom *at;
      int lowest_resno = 9999999;
      int highest_resno = -9999999;
      std::vector<coot::residue_spec_t> specs;
      for (int ires=0; ires<nres; ires++) {
	 residue_p = chain_p->GetResidue(ires);
	 int res_no = residue_p->GetSeqNum();
	 if (res_no < lowest_resno)
	    lowest_resno = res_no;
	 if (res_no > highest_resno)
	    highest_resno = res_no;
	 specs.push_back(residue_p);
      }
      std::sort(specs.begin(), specs.end());
      int prev_resno = lowest_resno - 1;
      for (unsigned int i=0; i<specs.size(); i++) {
	 if (specs[i].resno > (prev_resno + 1)) {
	    n_residue_count++; // add a spacer for a gap
	 }
	 prev_resno = specs[i].resno; // for next round
	 n_residue_count++;
      }
      std::cout << "   Chain " << chain_p->GetChainID() << " "
		<< "min resno " << lowest_resno << " n_residue_count "
		<< n_residue_count << std::endl;
      exptl::nsv::chain_length_residue_units_t p(chain_p->GetChainID(),
						 n_residue_count,
						 lowest_resno,
						 highest_resno);
      chain_length_residue_units.push_back(p);
   }

   return chain_length_residue_units;
}

// I feel this is the wrong thing to pass here, because we haven't
// checked the graph group.
// 
void
exptl::nsv::draw_axes(std::vector<chain_length_residue_units_t> clru,
		      int lrn, int brn, double x_offset) {

   GtkCanvasItem *item;
   GtkCanvasPoints *points = gtk_canvas_points_new(2);
   float font_scaler = pixels_per_letter;


   // ticks and text
   int irn_start = tick_start_number(lrn);
   // 2 sets of x-axis and ticks and tick numbers, one above and one
   // below the sequences
   for (int i_ax_pos=0; i_ax_pos<2; i_ax_pos++) {
      // x axis

      double y_value = 0;
      if (i_ax_pos == 1)
	 y_value = -1.0 * double(clru.size() * pixels_per_chain) - 2.0;

      // old values, bad with offset resnos.
      points->coords[0] = lrn*font_scaler;
      points->coords[1] = y_value;
      points->coords[2] = brn*font_scaler;
      points->coords[3] = y_value;

      points->coords[0] = 5 - x_offset; // don't extend too far to the left
      points->coords[1] = y_value;
      points->coords[2] = (brn-lrn)*font_scaler - x_offset;
      points->coords[3] = y_value;

      double points_max = 22500; // 23000 too many (strangely)
      if (points->coords[2] > points_max)
	 points->coords[2] = points_max;

      double tick_length = 3.0;
      if (i_ax_pos == 1)
	 tick_length = -3.0;
      
      double y_for_text = 10.0;
      if (i_ax_pos == 1)
	 y_for_text = y_value - 10.0;
      
      item = gtk_canvas_item_new(gtk_canvas_root(canvas),
				 GTK_CANVAS_TYPE_CANVAS_LINE,
				 "width_pixels", 1,
				 "points", points,
				 "fill_color", "black",
				 NULL);
      canvas_item_vec.push_back(item);

      // tick marks and tick labels
      for (int irn=irn_start; irn<=brn; irn+=5) {

	 points->coords[0] = (irn-lrn+1)*font_scaler - x_offset;
	 points->coords[1] = y_value;
	 points->coords[2] = (irn-lrn+1)*font_scaler - x_offset;
	 points->coords[3] = double(y_value + tick_length);

	 // Don't draw things that are too wide (i.e. to much x) for X to handle.
	 if (points->coords[2] < points_max) { 
	       
	    double x = (irn-lrn+1)*font_scaler - x_offset; // x for resno label
	    std::string lab = coot::util::int_to_string(irn);
	    item = gtk_canvas_item_new(gtk_canvas_root(canvas),
				       GTK_CANVAS_TYPE_CANVAS_LINE,
				       "width_pixels", 1,
				       "points", points,
				       "fill_color", "black",
				       NULL);
	    canvas_item_vec.push_back(item);
	    item = gtk_canvas_item_new(gtk_canvas_root(canvas),
				       GTK_CANVAS_TYPE_CANVAS_TEXT,
				       "text", lab.c_str(),
				       "x", x,
				       "y", y_for_text,
				       "anchor", GTK_ANCHOR_CENTER,
				       "font", fixed_font_str.c_str(),
				       "fill_color", "black",
				       NULL);
	    canvas_item_vec.push_back(item);
	 } 
      }
   }

}

int 
exptl::nsv::tick_start_number(int low_res_no) const {
   // 1 -> 5
   // 12 -> 15
   int j = (low_res_no-1)/5;
   return (j+1)*5;
}

void
exptl::nsv::origin_marker() {

   if (0) { 
      GtkCanvasItem *item;
      GtkCanvasPoints *points = gtk_canvas_points_new(5);
      
      points->coords[0] = 5;
      points->coords[1] = 5;
      points->coords[2] = 5;
      points->coords[3] = -5;
      points->coords[4] = -5;
      points->coords[5] = -5;
      points->coords[6] = -5;
      points->coords[7] = 5;
      points->coords[8] = 5;
      points->coords[9] = 5;
      item = gtk_canvas_item_new(gtk_canvas_root(canvas),
				 GTK_CANVAS_TYPE_CANVAS_LINE,
				 "width_pixels", 1,
				 "points", points,
				 "fill_color", "blue",
				 NULL);
   }
}


std::string
exptl::nsv::colour_by_secstr(CResidue *residue_p) const {

   std::string s("black");

   switch (residue_p->SSE)  {

   case SSE_Strand : s = "firebrick3";  break;
   case SSE_Bulge  : s = "firebrick1";  break;
   case SSE_3Turn  : s = "MediumBlue";  break;
   case SSE_4Turn  : s = "SteelBlue4";  break;
   case SSE_5Turn  : s = "DodgerBlue4"; break;
   case SSE_Helix  : s = "navy";        break;
   case SSE_None   : s = "black";       break;
   } 
   return s;
}


#endif // defined(HAVE_GTK_CANVAS) || defined(HAVE_GNOME_CANVAS)
