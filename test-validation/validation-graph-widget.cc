#include "validation-graph-widget.hh"
#include <algorithm>
struct _CootValidationGraph {
    GtkWidget parent;

    std::unique_ptr<coot::validation_information_t> _vi;
};

const int COOT_CHAIN_HEIGHT = 80;
const int COOT_RESIDUE_WIDTH = 5;

size_t max_chain_residue_count(CootValidationGraph* self) {
    return std::max_element(self->_vi->cviv.begin(),self->_vi->cviv.end(),
        [](const auto& lhs, const auto& rhs){
            return lhs.rviv.size() < rhs.rviv.size();
        })->rviv.size();
}

G_BEGIN_DECLS

G_DEFINE_TYPE(CootValidationGraph, coot_validation_graph, GTK_TYPE_WIDGET)

// not sure what this is for or whether it is going to be needed at all
// struct _CootValidationGraphClass {
//     GObjectClass parent_class;
// };


void coot_validation_graph_snapshot (GtkWidget *widget, GtkSnapshot *snapshot)
{
    GdkRGBA green, yellow;
    float w, h;

    //   gdk_rgba_parse (&red, "red");
    gdk_rgba_parse (&green, "green");
    gdk_rgba_parse (&yellow, "yellow");
    //   gdk_rgba_parse (&blue, "blue");


    w = (float) gtk_widget_get_width (widget);
    h = (float) gtk_widget_get_height (widget);


    CootValidationGraph* self = COOT_COOT_VALIDATION_GRAPH(widget);
    if(self->_vi) {
        graphene_rect_t m_graphene_rect = GRAPHENE_RECT_INIT(0, 0, w, h);
        cairo_t* cairo_canvas = gtk_snapshot_append_cairo(snapshot,&m_graphene_rect);
        // 1. Draw text
        cairo_move_to(cairo_canvas,0,0);
        cairo_set_source_rgb(cairo_canvas, 0.1, 0.1, 0.1);
        cairo_set_line_width(cairo_canvas,1);
        cairo_line_to(cairo_canvas, w, h);
        cairo_set_font_size(cairo_canvas,14);
        cairo_select_font_face (cairo_canvas, "sans-serif", CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
        const char* title = self->_vi->name.c_str();
        g_debug("title: %s",title);
        cairo_show_text(cairo_canvas, title);
        // 2. Draw axes
        g_free(cairo_canvas);

        

        float base_height = 0;
        float width_step = w / (float) max_chain_residue_count(self);
        float height_step = h / (float) self->_vi->cviv.size();
        for(const auto& cvi: self->_vi->cviv) {
            float base_width = 0;
            base_height += height_step;
            for(const auto& ri: cvi.rviv) {
                float bar_height = height_step * ri.distortion / 30.d;
                m_graphene_rect = GRAPHENE_RECT_INIT(base_width + 1, base_height - bar_height, width_step - 1, bar_height);
                float border_width = 1;
                GskRoundedRect outline;
                gsk_rounded_rect_init_from_rect(
                    &outline,
                    &m_graphene_rect,
                    0
                );
                gtk_snapshot_append_color(snapshot, &green, &m_graphene_rect);
                // border is a bit buggy
                gtk_snapshot_append_border(snapshot, &outline , &border_width, &yellow);
                base_width += width_step;
            }
        }
    } else {
        // do nothing
    }
}

void coot_validation_graph_measure
    (GtkWidget      *widget,
    GtkOrientation  orientation,
    int             for_size,
    int            *minimum_size,
    int            *natural_size,
    int            *minimum_baseline,
    int            *natural_baseline)
{
    CootValidationGraph* self = COOT_COOT_VALIDATION_GRAPH(widget);
    if (self->_vi) {
        switch (orientation)
        {
        case GTK_ORIENTATION_HORIZONTAL:{
            auto max_chain_residues = max_chain_residue_count(self);
            *minimum_size = max_chain_residues * COOT_RESIDUE_WIDTH;
            *natural_size = max_chain_residues * COOT_RESIDUE_WIDTH;
            break;
        }
        case GTK_ORIENTATION_VERTICAL:{
            auto num_of_chains = self->_vi->cviv.size();
            //g_debug("Num of chains: %u",num_of_chains);
            auto size = num_of_chains * COOT_CHAIN_HEIGHT;
            //g_debug("Vertical size: %u",size);
            *minimum_size = size;
            *natural_size = size;
            break;
        }
        default:
            break;
        }
    } else {
        // do nothing
    }
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

static void coot_validation_graph_dispose(GObject* _self) {

    CootValidationGraph* self = COOT_COOT_VALIDATION_GRAPH(_self);
    self->_vi.reset(nullptr);
    G_OBJECT_CLASS(coot_validation_graph_parent_class)->dispose(_self);
}


static void coot_validation_graph_class_init(CootValidationGraphClass* klass) {
    // I think that this is a GObject class constructor that sets up the GObject class at runtime.

    GTK_WIDGET_CLASS(klass)->snapshot = coot_validation_graph_snapshot;
    GTK_WIDGET_CLASS(klass)->measure = coot_validation_graph_measure;
    G_OBJECT_CLASS(klass)->dispose = coot_validation_graph_dispose;
    
}

CootValidationGraph* 
coot_validation_graph_new()
{
    return COOT_COOT_VALIDATION_GRAPH(g_object_new (COOT_VALIDATION_GRAPH_TYPE, NULL));
}



G_END_DECLS

void coot_validation_graph_set_validation_information(CootValidationGraph* self, std::unique_ptr<coot::validation_information_t> vi) {
    self->_vi = std::move(vi);
}