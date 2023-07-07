
#include <gtk/gtk.h>

#include "sequence-view-widget.hh"
#include "coot-utils/coot-coord-utils.hh"

struct _CootSequenceView {
   GtkWidget parent;
   mmdb::Manager *mol;
};
static guint sequence_view_residue_clicked_signal;

G_BEGIN_DECLS

G_DEFINE_TYPE(CootSequenceView, coot_sequence_view, GTK_TYPE_WIDGET)

void coot_sequence_view_snapshot (GtkWidget *widget, GtkSnapshot *snapshot) {

   std::cout << "sequence view snapshot" << std::endl;

   GtkAllocation allocation;
   gtk_widget_get_allocation(GTK_WIDGET(widget), &allocation);
   float w = allocation.width;
   float h = allocation.height;
   GdkRGBA attribute_color;
   gdk_rgba_parse (&attribute_color, "#aaddff");
   graphene_rect_t m_graphene_rect = GRAPHENE_RECT_INIT(0, 0, w, h);
   cairo_t* cairo_canvas = gtk_snapshot_append_cairo(snapshot,&m_graphene_rect);
   cairo_set_source_rgb(cairo_canvas, attribute_color.red, attribute_color.green, attribute_color.blue);
   
}

gboolean sequence_view_query_tooltip (CootSequenceView* self,
                                      gint x,
                                      gint y,
                                      gboolean keyboard_mode,
                                      GtkTooltip* tooltip,
                                      gpointer user_data) {
   return FALSE; // for now
}

static void on_sequence_view_left_click(GtkGestureClick* gesture_click,
                                        gint n_press,
                                        gdouble x,
                                        gdouble y,
                                        gpointer user_data ) {

   CootSequenceView* self = COOT_COOT_SEQUENCE_VIEW(user_data);
   g_debug("On click at widget: %p, at x: %f, y: %f", self, x, y);
   // something here
}


static void coot_sequence_view_init(CootSequenceView* self) {

   gtk_widget_set_has_tooltip(GTK_WIDGET(self),TRUE);
   g_signal_connect(self,"query-tooltip",G_CALLBACK(sequence_view_query_tooltip),NULL);

   GtkGesture* click_controller = gtk_gesture_click_new();
   // GtkEventController* hover_controller = gtk_event_controller_motion_new();

   // left mouse button
   gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(click_controller),GDK_BUTTON_PRIMARY);
   g_signal_connect(click_controller, "pressed", G_CALLBACK(on_sequence_view_left_click), self);

   // g_signal_connect(hover_controller,"motion",G_CALLBACK(on_hover),self);

   gtk_widget_add_controller(GTK_WIDGET(self),GTK_EVENT_CONTROLLER(click_controller));
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

void coot_sequence_view_set_structure(CootSequenceView* self, mmdb::Manager *mol) {

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
