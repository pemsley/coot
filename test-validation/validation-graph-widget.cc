#include "validation-graph-widget.hh"

G_BEGIN_DECLS
struct _CootValidationGraph {
    GtkWidget parent;
};

G_DEFINE_TYPE(CootValidationGraph, coot_validation_graph, GTK_TYPE_WIDGET)

// not sure what this is for or whether it is going to be needed at all
// struct _CootValidationGraphClass {
//     GObjectClass parent_class;
// };


/// This is a "demo" function copied from here: https://blog.gtk.org/2020/04/
void
demo_snapshot (GtkWidget *widget, GtkSnapshot *snapshot)
{
  GdkRGBA red, green, yellow, blue;
  float w, h;

  gdk_rgba_parse (&red, "red");
  gdk_rgba_parse (&green, "green");
  gdk_rgba_parse (&yellow, "yellow");
  gdk_rgba_parse (&blue, "blue");

  w = gtk_widget_get_width (widget) / 2.0;
  h = gtk_widget_get_height (widget) / 2.0;

  graphene_rect_t m_graphene_rect;

  m_graphene_rect = GRAPHENE_RECT_INIT(0, 0, w, h);
  gtk_snapshot_append_color (snapshot, &red,
                             &m_graphene_rect);
  m_graphene_rect = GRAPHENE_RECT_INIT(w, 0, w, h);
  gtk_snapshot_append_color (snapshot, &green,
                             &m_graphene_rect);
  m_graphene_rect = GRAPHENE_RECT_INIT(0, h, w, h);
  gtk_snapshot_append_color (snapshot, &yellow,
                             &m_graphene_rect);
  m_graphene_rect = GRAPHENE_RECT_INIT(w, h, w, h);
  gtk_snapshot_append_color (snapshot, &blue,
                             &m_graphene_rect);
}

void
demo_measure (GtkWidget      *widget,
              GtkOrientation  orientation,
              int             for_size,
              int            *minimum_size,
              int            *natural_size,
              int            *minimum_baseline,
              int            *natural_baseline)
{
  *minimum_size = 100;
  *natural_size = 200;
}

static void on_left_click (
  GtkGestureClick* gesture_click,
  gint n_press,
  gdouble x,
  gdouble y,
  gpointer user_data
) {
    CootValidationGraph* self = COOT_COOT_VALIDATION_GRAPH(user_data);
    g_debug("On click at widget: %p, at x: %f, y: %f",self,x,y);
    gtk_gesture_set_state(GTK_GESTURE(gesture_click),GTK_EVENT_SEQUENCE_CLAIMED);
}

static void coot_validation_graph_init(CootValidationGraph* self) {
    // I think that this is the primary constructor


    GtkGesture* click_controller = gtk_gesture_click_new();
    // left mouse button
    gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(click_controller),GDK_BUTTON_PRIMARY);

    g_signal_connect(click_controller,"pressed",G_CALLBACK(on_left_click),self);

    gtk_widget_add_controller(GTK_WIDGET(self),GTK_EVENT_CONTROLLER(click_controller));
}

static void coot_validation_graph_class_init(CootValidationGraphClass* klass) {
    // I think that this is a GObject class constructor that sets up the GObject class at runtime.

    GTK_WIDGET_CLASS(klass)->snapshot = demo_snapshot;
    GTK_WIDGET_CLASS(klass)->measure = demo_measure;
}

CootValidationGraph* 
coot_validation_graph_new()
{
    return COOT_COOT_VALIDATION_GRAPH(g_object_new (COOT_VALIDATION_GRAPH_TYPE, NULL));
}



G_END_DECLS
