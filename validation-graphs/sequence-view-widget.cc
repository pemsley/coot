
#include <gtk/gtk.h>

#include "cairo.h"
#include "sequence-view-widget.hh"
#include "coot-utils/coot-coord-utils.hh"

const float RESIDUE_BOX_HEIGHT   = 12.0;
const float RESIDUE_BOX_WIDTH    = 12.0;
const float Y_OFFSET_PER_CHAIN   = 16.0;
const float X_OFFSET_PER_RESIDUE = 12.0;
const float TICK_LINE_WIDTH  = 2.0;
const float TICK_LINE_LENGTH = 8.0;

struct _CootSequenceView {
   GtkWidget parent;
   int imol;
   mmdb::Manager *mol;
   std::vector<box_info_t> box_info_store;
};

static guint sequence_view_residue_clicked_signal;

G_BEGIN_DECLS

G_DEFINE_TYPE(CootSequenceView, coot_sequence_view, GTK_TYPE_WIDGET)

void coot_sequence_view_snapshot(GtkWidget *widget, GtkSnapshot *snapshot) {

   CootSequenceView* self = COOT_COOT_SEQUENCE_VIEW(widget);

   const float x_offset_base = 30.0;
   const float y_offset_base = 30.0;

   self->box_info_store.clear();

   // This function is used from the loop below.
   //
   auto add_box_letter_code_label = [widget, snapshot] (mmdb::Residue *residue_p, float x_base, float y_base) {
      std::string res_name = residue_p->GetResName();
      std::string slc = coot::util::three_letter_to_one_letter_with_specials(res_name);

      float w_pixels_label = 40.0;
      float h_pixels_label = 40.0;
      
      graphene_rect_t m_graphene_rect = GRAPHENE_RECT_INIT(x_base, y_base, x_base + w_pixels_label, y_base + h_pixels_label);
      cairo_t *cairo_canvas = gtk_snapshot_append_cairo(snapshot, &m_graphene_rect);
      GdkRGBA residue_color; // maybe per-residue colouring later
      gdk_rgba_parse(&residue_color, "#222222");
      cairo_set_source_rgb(cairo_canvas, residue_color.red, residue_color.green, residue_color.blue);

      PangoLayout* pango_layout = pango_layout_new(gtk_widget_get_pango_context(widget));
      std::string label_markup = std::string("<tt>") + slc + std::string("</tt>");
      pango_layout_set_markup(pango_layout, label_markup.c_str(), -1);
      int layout_width, layout_height;
      pango_layout_get_pixel_size(pango_layout, &layout_width, &layout_height);

      cairo_move_to(cairo_canvas, x_base +2 , y_base - 2);
      pango_cairo_show_layout(cairo_canvas, pango_layout);
   };

   auto add_chain_label = [widget, snapshot] (int i_chain, const std::string &chain_id, float x_base, float y_base) {
      
      float w_pixels_label = 12.0;
      float h_pixels_label = 12.0;

      graphene_rect_t m_graphene_rect = GRAPHENE_RECT_INIT(x_base, y_base, x_base + w_pixels_label, y_base + h_pixels_label);
      cairo_t* cairo_canvas = gtk_snapshot_append_cairo(snapshot, &m_graphene_rect);
      GdkRGBA residue_color; // maybe per-residue colouring later
      gdk_rgba_parse(&residue_color, "#222222");
      cairo_set_source_rgb(cairo_canvas, residue_color.red, residue_color.green, residue_color.blue);

      PangoLayout* pango_layout = pango_layout_new(gtk_widget_get_pango_context(widget));
      std::string label_markup = std::string("<span weight=\"bold\">") + chain_id + std::string("</span>");
      pango_layout_set_markup(pango_layout, label_markup.c_str(), -1);
      int layout_width, layout_height;
      pango_layout_get_pixel_size(pango_layout, &layout_width, &layout_height);

      cairo_move_to(cairo_canvas, x_base +2 , y_base - 2);
      pango_cairo_show_layout(cairo_canvas, pango_layout);
   };

   auto get_min_max_residue_number = [] (mmdb::Model *model_p) {
      int n_chains = model_p->GetNumberOfChains();
      std::pair<int, int> min_max(10000, -10000);
      bool min_max_was_set = false;
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         std::pair<bool, std::pair<int, int> > min_max_for_chain = coot::util::min_max_residues_in_polymer_chain(chain_p);
         if (min_max_for_chain.first) {
            min_max_was_set = true;
            if (min_max_for_chain.second.first  < min_max.first)  min_max.first  = min_max_for_chain.second.first;
            if (min_max_for_chain.second.second > min_max.second) min_max.second = min_max_for_chain.second.second;
         }
      }
      return std::pair<bool, std::pair<int, int> > (min_max_was_set, min_max);
   };

   auto calculate_min_max_by_5s = +[] (const std::pair<int, int> &mm) {
      int f5 = mm.first/5;
      int f = f5 * 5;
      int s5 = mm.second/5;
      int s = s5 * 5;
      if (f > mm.first)  f -= 5;
      if (s < mm.second) f += 5;
      return std::pair<int, int>(f, s);
   };

   auto add_tick_marks = [calculate_min_max_by_5s, get_min_max_residue_number, x_offset_base, y_offset_base]
      (cairo_t *cairo_canvas, mmdb::Model *model_p) {

      std::pair<bool, std::pair<int, int> > mm = get_min_max_residue_number(model_p);
      int n_chains = model_p->GetNumberOfChains();
      bool min_max_was_set = mm.first;
      if (min_max_was_set) {
         // calculate these, rounding up to the nearest 5
         std::pair<int, int> mmr = calculate_min_max_by_5s(mm.second);
         int res_no_start =  mmr.first;
         int res_no_end   = mmr.second; // inclusive, less than or equal to the max residue number
         cairo_set_line_width(cairo_canvas, TICK_LINE_WIDTH);
         for (int ires = res_no_start;  ires<=res_no_end; ires+=5) {

            // above the top line
            double pos_x = x_offset_base + X_OFFSET_PER_RESIDUE * static_cast<double>(ires-mm.second.first+1) + X_OFFSET_PER_RESIDUE/2;
            // std::cout << "add_tick_marks(): ires: " << ires << " pos_x " << pos_x << std::endl;
            double pos_y_start = y_offset_base - 1.0;
            double pos_y_end   = pos_y_start - TICK_LINE_LENGTH;
            cairo_set_source_rgb(cairo_canvas, 0.1, 0.1, 0.1);
            cairo_move_to(cairo_canvas, pos_x, pos_y_start);
            cairo_line_to(cairo_canvas, pos_x, pos_y_end);
            cairo_stroke(cairo_canvas);

            // below the bottom line
            pos_x = x_offset_base + X_OFFSET_PER_RESIDUE * static_cast<double>(ires-mm.second.first+1) + X_OFFSET_PER_RESIDUE/2;
            pos_y_start = y_offset_base + static_cast<double>(n_chains) * Y_OFFSET_PER_CHAIN + -2.0;
            pos_y_end   = pos_y_start + TICK_LINE_LENGTH;
            cairo_set_source_rgb(cairo_canvas, 0.1, 0.1, 0.1);
            cairo_move_to(cairo_canvas, pos_x, pos_y_start);
            cairo_line_to(cairo_canvas, pos_x, pos_y_end);
            cairo_stroke(cairo_canvas);
         }
      }
   };
   
   auto add_tick_labels = [calculate_min_max_by_5s, get_min_max_residue_number, x_offset_base, y_offset_base]
      (cairo_t *cairo_canvas, mmdb::Model *model_p) {

      int n_chains = model_p->GetNumberOfChains();
      std::pair<bool, std::pair<int, int> > mm = get_min_max_residue_number(model_p);
      bool min_max_was_set = mm.first;
      if (min_max_was_set) {
         // as above
         // calculate these, rounding up to the nearest 5
         std::pair<int, int> mmr = calculate_min_max_by_5s(mm.second);
         int res_no_start = mmr.first;
         int res_no_end   = mmr.second; // inclusive, less than or equal to the max residue number
         for (int ires = res_no_start;  ires<=res_no_end; ires+=5) {

            // above the top line
            double pos_x = x_offset_base + X_OFFSET_PER_RESIDUE * static_cast<double>(ires-mm.second.first+1) + X_OFFSET_PER_RESIDUE/2;
            // float x_1 = static_cast<float>(res_no-mm.second.first+1) * X_OFFSET_PER_RESIDUE + x_offset_base;
            std::string text = std::to_string(ires);
            float l = text.length();
            pos_x -= 3.5 * l ; // so that the text is centred on the tick.
            double pos_y = y_offset_base - 2.0 - TICK_LINE_LENGTH;
            cairo_set_source_rgb(cairo_canvas, 0.1, 0.1, 0.1);
            cairo_move_to(cairo_canvas, pos_x, pos_y);
            cairo_show_text(cairo_canvas, text.c_str());
            // std::cout << "text top " << text << " at " << pos_x << " " << pos_y<< std::endl;

            // below the bottom line
            pos_y = y_offset_base + 8.0 + TICK_LINE_LENGTH + Y_OFFSET_PER_CHAIN * n_chains;
            cairo_set_source_rgb(cairo_canvas, 0.1, 0.1, 0.1);
            cairo_move_to(cairo_canvas, pos_x, pos_y);
            cairo_show_text(cairo_canvas, text.c_str());
            // std::cout << "text bot " << text << " at " << pos_x << " " << pos_y<< std::endl;
         }
      }
   };

   auto n_residues_and_n_chains = [] (mmdb::Model *model_p) {

      int n_chains = model_p->GetNumberOfChains();
      bool min_max_is_set = false;
      std::pair<int, int> min_max(10000, -10000);
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         std::pair<bool, std::pair<int, int> > mm = coot::util::min_max_residues_in_polymer_chain(chain_p);
         if (mm.first) {
            min_max_is_set = true;
            if (mm.second.first  < min_max.first)  min_max.first  = mm.second.first;
            if (mm.second.second > min_max.second) min_max.second = mm.second.second;
         }
      }
      if (min_max_is_set) {
         return std::make_pair(min_max.second - min_max.first + 1, n_chains);
      } else {
         return std::make_pair(0,0);
      }
   };

   if (self->mol) {

      int imod = 1;
      mmdb::Model *model_p = self->mol->GetModel(imod);
      if (model_p) {

         GtkAllocation allocation;
         gtk_widget_get_allocation(GTK_WIDGET(widget), &allocation);
         float w = allocation.width;
         float h = allocation.height;
         GdkRGBA attribute_color;
         gdk_rgba_parse(&attribute_color, "#223322");

         std::pair<int, int> n_res_and_n_chains = n_residues_and_n_chains(model_p);
         float w_pixels_rect = n_res_and_n_chains.first * X_OFFSET_PER_RESIDUE + x_offset_base ; // w and h of the box it sits in
         float h_pixels_rect = 100 + n_res_and_n_chains.second * Y_OFFSET_PER_CHAIN + y_offset_base;

         graphene_rect_t m_graphene_rect = GRAPHENE_RECT_INIT(0, 0, w_pixels_rect, h_pixels_rect);
         cairo_t* cairo_canvas = gtk_snapshot_append_cairo(snapshot,&m_graphene_rect);

         // make tick marks
         add_tick_marks(cairo_canvas, model_p);

         //
         add_tick_labels(cairo_canvas, model_p);

         // Make the labels for the chains
         //
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            std::string chain_id = chain_p->GetChainID();
            float y_offset = y_offset_base + static_cast<float>(ichain) * Y_OFFSET_PER_CHAIN;
            float x_offset_base = 15.0;
            add_chain_label(ichain, chain_id, x_offset_base, y_offset);
         }
         

         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            std::string chain_id = chain_p->GetChainID();
            int n_res = chain_p->GetNumberOfResidues();
            float y_offset = y_offset_base + static_cast<float>(ichain) * Y_OFFSET_PER_CHAIN;

            std::pair<bool, std::pair<int, int> > mm = get_min_max_residue_number(model_p);

            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  int res_no = residue_p->GetSeqNum();

                  float x_1 = static_cast<float>(res_no-mm.second.first+1) * X_OFFSET_PER_RESIDUE + x_offset_base;
                  float y_1 = y_offset;
                  float x_2 = RESIDUE_BOX_WIDTH;  // a delta
                  float y_2 = RESIDUE_BOX_HEIGHT; // a delta
                  m_graphene_rect = GRAPHENE_RECT_INIT(x_1, y_1, x_2, y_2);

                  GdkRGBA residue_color; // maybe per-residue colouring later
                  gdk_rgba_parse(&residue_color, "#eeeeee");
                  gtk_snapshot_append_color(snapshot, &residue_color, &m_graphene_rect);

                  if (false)
                     std::cout << "ires " << ires << " res_no " << res_no << " mm.first " << mm.second.first  << " x_1 " << x_1 << " "
                               << coot::residue_spec_t(residue_p) << " " << residue_p->GetResName() << std::endl;

                  box_info_t box_info(residue_p, x_1, y_1);
                  self->box_info_store.push_back(box_info);
                  add_box_letter_code_label(residue_p, x_1, y_1);
               }
            }
         }
         //    g_object_unref(pango_layout); // for the title
         cairo_destroy(cairo_canvas);
      }
   } else {
      std::cout << "error in coot_sequence_view_snapshot() null mol " << std::endl;
   }

}

gboolean sequence_view_query_tooltip(CootSequenceView* self,
                                     gint x,
                                     gint y,
                                     gboolean keyboard_mode,
                                     GtkTooltip* tooltip,
                                     gpointer user_data) {
   return FALSE; // for now
}

void
find_the_clicked_residue(CootSequenceView *sv, float x, float y) {

   coot::residue_spec_t best_spec;
   float best_dist_initial = 14.0;
   float best_dist = best_dist_initial;
   for (unsigned int i=0; i<sv->box_info_store.size(); i++) {
      const auto &box = sv->box_info_store[i];
      float d_x = box.x_base - x + RESIDUE_BOX_WIDTH/2;  // Offset from corner to middle of box applied here.
      float d_y = box.y_base - y + RESIDUE_BOX_HEIGHT/2; // ditto.
      float dd = d_x * d_x + d_y * d_y;
      float d = sqrtf(dd);
      if (d < best_dist) {
         best_dist = d;
         best_spec = box.residue_spec;
      }
   }
   if (best_dist < best_dist_initial) {
      mmdb::Residue *residue_p = coot::util::get_residue(best_spec, sv->mol);
      if (residue_p) {
         std::string rn = residue_p->GetResName();
         std::cout << "   closest clicked residue " << best_spec << " " << rn << std::endl;
      }
   }
}

// static
void on_sequence_view_left_click(GtkGestureClick* gesture_click,
                                 gint n_press,
                                 gdouble x,
                                 gdouble y,
                                 gpointer user_data) {

   // user_data is the CootSequenceView* self.

   CootSequenceView* self = COOT_COOT_SEQUENCE_VIEW(user_data);
   // std::cout << "--- on_sequence_view_left_click() " << self << " x " << x << " y " << y << std::endl;
   find_the_clicked_residue(self, x, y);
}


static void coot_sequence_view_init(CootSequenceView* self) {

   gtk_widget_set_has_tooltip(GTK_WIDGET(self),TRUE);
   g_signal_connect(self, "query-tooltip", G_CALLBACK(sequence_view_query_tooltip), NULL);

   GtkGesture* click_controller = gtk_gesture_click_new();
   // GtkEventController* hover_controller = gtk_event_controller_motion_new();

   // left mouse button
   gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(click_controller), GDK_BUTTON_PRIMARY);
   g_signal_connect(click_controller, "pressed", G_CALLBACK(on_sequence_view_left_click), self);

   // g_signal_connect(hover_controller,"motion",G_CALLBACK(on_hover),self);

   gtk_widget_add_controller(GTK_WIDGET(self), GTK_EVENT_CONTROLLER(click_controller));
   // gtk_widget_add_controller(GTK_WIDGET(self),GTK_EVENT_CONTROLLER(hover_controller));

}


static void coot_sequence_view_dispose(GObject* _self) {
   CootSequenceView* self = COOT_COOT_SEQUENCE_VIEW(_self);
   // clean up self here
   G_OBJECT_CLASS(coot_sequence_view_parent_class)->dispose(_self);
}

static void coot_sequence_view_class_init(CootSequenceViewClass* klass) {

    // I think that this is a GObject class constructor that sets up the GObject class at runtime.
   sequence_view_residue_clicked_signal = g_signal_new("residue-clicked",
        G_TYPE_FROM_CLASS (klass),
        (GSignalFlags) (G_SIGNAL_RUN_LAST | G_SIGNAL_NO_RECURSE | G_SIGNAL_NO_HOOKS),
        0 /* class offset.Subclass cannot override the class handler (default handler). */,
        NULL /* accumulator */,
        NULL /* accumulator data */,
        NULL /* C marshaller. g_cclosure_marshal_generic() will be used */,
        G_TYPE_NONE /* return_type */,
        1     /* n_params */,
        G_TYPE_POINTER
    );
    GTK_WIDGET_CLASS(klass)->snapshot = coot_sequence_view_snapshot;
    G_OBJECT_CLASS(klass)->dispose    = coot_sequence_view_dispose;

}


CootSequenceView *
coot_sequence_view_new() {

    return COOT_COOT_SEQUENCE_VIEW(g_object_new (COOT_SEQUENCE_VIEW_TYPE, NULL));
}

G_END_DECLS

void coot_sequence_view_set_structure(CootSequenceView* self, int imol, mmdb::Manager *mol) {


   if (false) { // debug
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            std::string chain_id = chain_p->GetChainID();
            std::cout << "------ " << chain_id << " ----" << std::endl;
            int n_res = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  std::string res_name = residue_p->GetResName();
                  std::string slc = coot::util::three_letter_to_one_letter_with_specials(res_name);
                  std::cout << slc;
               }
            }
            std::cout << "\n";
         }
      }
   }

   self->imol = imol;
   self->mol = mol;

}
