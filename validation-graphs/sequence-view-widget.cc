
#include <gtk/gtk.h>

#include "sequence-view-widget.hh"
#include "coot-utils/coot-coord-utils.hh"

/// For drawing the main title
const int TITLE_HEIGHT = -25; // 30;
const float RESIDUE_BOX_HEIGHT   = 12.0;
const float RESIDUE_BOX_WIDTH    = 12.0;
const float Y_OFFSET_PER_CHAIN   = 16.0;
const float X_OFFSET_PER_RESIDUE = 12.0;

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

   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(widget), &allocation);
   float w = allocation.width;
   float h = allocation.height;
   GdkRGBA attribute_color;
   gdk_rgba_parse(&attribute_color, "#223322");

   float w_pixels_title = 600; // w and h of the box it sits in
   float h_pixels_title = 100;

   // title
   graphene_rect_t m_graphene_rect = GRAPHENE_RECT_INIT(0, 0, w_pixels_title, h_pixels_title);
   cairo_t *cairo_canvas = gtk_snapshot_append_cairo(snapshot, &m_graphene_rect);
   cairo_set_source_rgb(cairo_canvas, attribute_color.red, attribute_color.green, attribute_color.blue);

   PangoLayout *pango_layout = pango_layout_new(gtk_widget_get_pango_context(widget));
   std::string name = "Coooooooooooooot Sequence View";
   std::string title_markup = std::string("<span size=\"large\" weight=\"bold\">") + name + std::string("</span>");
   pango_layout_set_markup(pango_layout, title_markup.c_str(), -1);
   int layout_width, layout_height;
   pango_layout_get_pixel_size(pango_layout, &layout_width, &layout_height);
   // std::cout << "layout_width " << layout_width << " layout_height " << layout_height << std::endl;
   cairo_move_to(cairo_canvas, (w - layout_width) / 2.f, (TITLE_HEIGHT + layout_height) / 2.f);
   pango_cairo_show_layout(cairo_canvas, pango_layout);


   float x_offset_base = 30.0;
   float y_offset_base = 30.0; // down because of title - this should be changed later

   self->box_info_store.clear();

   // This function is used from the loop below.
   //
   auto add_box_letter_code_label = [widget, snapshot] (mmdb::Residue *residue_p, float x_base, float y_base) {
      std::string res_name = residue_p->GetResName();
      std::string slc = coot::util::three_letter_to_one_letter_with_specials(res_name);

      float w_pixels_label = 40.0;
      float h_pixels_label = 40.0;
      
      graphene_rect_t m_graphene_rect = GRAPHENE_RECT_INIT(x_base, y_base, x_base + w_pixels_label, y_base + h_pixels_label);
      cairo_t* cairo_canvas = gtk_snapshot_append_cairo(snapshot, &m_graphene_rect);
      GdkRGBA residue_color; // maybe per-residue colouring later
      gdk_rgba_parse(&residue_color, "#222222");
      cairo_set_source_rgb(cairo_canvas, residue_color.red, residue_color.green, residue_color.blue);

      PangoLayout* pango_layout = pango_layout_new(gtk_widget_get_pango_context(widget));
      std::string label_markup = slc;
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
   


   if (self->mol) {

      int imod = 1;
      mmdb::Model *model_p = self->mol->GetModel(imod);
      if (model_p) {

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
            for (int ires=0; ires<n_res; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               if (residue_p) {
                  int res_no = residue_p->GetSeqNum();

                  float x_offset = x_offset_base;
                  float x_1 = static_cast<float>(res_no) * X_OFFSET_PER_RESIDUE + x_offset;
                  float y_1 = y_offset;
                  float x_2 = RESIDUE_BOX_WIDTH;  // a delta
                  float y_2 = RESIDUE_BOX_HEIGHT; // a delta
                  m_graphene_rect = GRAPHENE_RECT_INIT(x_1, y_1, x_2, y_2);

                  GdkRGBA residue_color; // maybe per-residue colouring later
                  gdk_rgba_parse(&residue_color, "#eeeeee");
                  gtk_snapshot_append_color(snapshot, &residue_color, &m_graphene_rect);

                  box_info_t box_info(residue_p, x_1, y_1);
                  self->box_info_store.push_back(box_info);
                  add_box_letter_code_label(residue_p, x_1, y_1);
               }
            }
         }
      }
   } else {
      std::cout << "error in coot_sequence_view_snapshot() null mol " << std::endl;
   }
   g_object_unref(pango_layout);
   cairo_destroy(cairo_canvas);
   
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
