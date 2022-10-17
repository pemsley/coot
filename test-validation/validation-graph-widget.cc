#include "validation-graph-widget.hh"
#include <algorithm>
struct _CootValidationGraph {
    GtkWidget parent;

    std::unique_ptr<coot::validation_information_t> _vi;
};

/// Basis for max bar height
const int CHAIN_HEIGHT = 120;
/// Used for allocating space for axes and labels
const int CHAIN_SPACING = 40;
const int RESIDUE_WIDTH = 3;
/// Breathing space for residue rectangle's borders
const int RESIDUE_SPACING = 3;
/// For drawing the main title
const int TITLE_HEIGHT = 0;
/// Space for the axis to be drawn on the left side of the graph
const int AXIS_MARGIN = 20;
const double AXIS_LINE_WIDTH = 2;
const float RESIDUE_BORDER_WIDTH = 1;

size_t max_chain_residue_count(CootValidationGraph* self) {
    return std::max_element(self->_vi->cviv.cbegin(),self->_vi->cviv.cend(),
        [](const auto& lhs, const auto& rhs){
            return lhs.rviv.size() < rhs.rviv.size();
        })->rviv.size();
}

double max_chain_residue_distortion(const std::vector<coot::residue_validation_information_t>& rviv) {
    return std::max_element(rviv.cbegin(),rviv.cend(),
    [](const auto& lhs, const auto& rhs){
        return lhs.distortion < rhs.distortion;
    }
    )->distortion;
}

G_BEGIN_DECLS

G_DEFINE_TYPE(CootValidationGraph, coot_validation_graph, GTK_TYPE_WIDGET)

// not sure what this is for or whether it is going to be needed at all
// struct _CootValidationGraphClass {
//     GObjectClass parent_class;
// };


void coot_validation_graph_snapshot (GtkWidget *widget, GtkSnapshot *snapshot)
{
    // attribute_color is used for drawing labels and axes
    GdkRGBA residue_color, border_color, attribute_color;

    gdk_rgba_parse (&residue_color, "#008000");
    gdk_rgba_parse (&border_color, "#002000");

    // gdk_rgba_parse (&attribute_color, "#ffffff");
    // Gtk 4.10 ?
    // gtk_widget_get_style_color(widget,&attribute_color);
    GtkStyleContext* style_context = gtk_widget_get_style_context(widget);
    gtk_style_context_get_color(style_context,&attribute_color);



    float w = (float) gtk_widget_get_width (widget);
    float h = (float) gtk_widget_get_height (widget);


    CootValidationGraph* self = COOT_COOT_VALIDATION_GRAPH(widget);
    if(self->_vi) {
        // 1. Draw title
        graphene_rect_t m_graphene_rect = GRAPHENE_RECT_INIT(0, 0, w, h);
        cairo_t* cairo_canvas = gtk_snapshot_append_cairo(snapshot,&m_graphene_rect);
        cairo_set_source_rgb(cairo_canvas, attribute_color.red, attribute_color.green, attribute_color.blue);
        
        // This does not respect GTK theming
        // PangoLayout* pango_layout = pango_cairo_create_layout(cairo_canvas);
        PangoLayout* pango_layout = pango_layout_new(gtk_widget_get_pango_context(widget));
        pango_layout_set_text(pango_layout,self->_vi->name.c_str(),-1);
        int layout_width, layout_height;
        pango_layout_get_pixel_size(pango_layout,&layout_width,&layout_height);
        cairo_move_to(cairo_canvas,(w - layout_width) / 2,(TITLE_HEIGHT + layout_height) / 2);
        pango_cairo_show_layout(cairo_canvas, pango_layout);

        // I can't get this to render the text where it needs to be, so I'm using cairo directly
        // gtk_snapshot_append_layout(snapshot,pango_layout,&attribute_color);
        // A GtkLabel as a child widget could also be used, but I have no idea how to manage layout inside widgets
        g_object_unref(pango_layout);


        float base_height = TITLE_HEIGHT;
        float width_step = (w - (float) AXIS_MARGIN) / (float) max_chain_residue_count(self);
        float height_diff = (h - (float)TITLE_HEIGHT - (float) CHAIN_SPACING / 2.f) / ((float) self->_vi->cviv.size()) - (CHAIN_HEIGHT + CHAIN_SPACING);

        cairo_set_line_width(cairo_canvas,AXIS_LINE_WIDTH);

        for(const auto& chain: self->_vi->cviv) {
            m_graphene_rect = GRAPHENE_RECT_INIT(0, 0, w, h);

            // Label chain
            pango_layout = pango_layout_new(gtk_widget_get_pango_context(widget));
            std::string chain_label = "Chain " + chain.chain_id;
            pango_layout_set_text(pango_layout,chain_label.c_str(),-1);
            pango_layout_get_pixel_size(pango_layout,&layout_width,&layout_height);
            cairo_move_to(cairo_canvas,0,base_height + layout_height / 2);
            pango_cairo_show_layout(cairo_canvas, pango_layout);

            g_object_unref(pango_layout);
            // Draw axes
            float axis_y_offset = base_height + CHAIN_SPACING / 2.f;
            cairo_move_to(cairo_canvas,0, axis_y_offset);
            
            cairo_line_to(cairo_canvas, 0, axis_y_offset + CHAIN_HEIGHT + CHAIN_SPACING / 2.f);
            cairo_stroke(cairo_canvas);

            
            base_height += CHAIN_HEIGHT + CHAIN_SPACING / 2.f;
            cairo_move_to(cairo_canvas,0, base_height + CHAIN_SPACING / 2.f);
            
            cairo_line_to(cairo_canvas, w, base_height + CHAIN_SPACING / 2.f);
            cairo_stroke(cairo_canvas);

            float base_width = AXIS_MARGIN;
            const double normalization_divisor = max_chain_residue_distortion(chain.rviv);
            for(const auto& residue: chain.rviv) {
                float bar_height = CHAIN_HEIGHT * residue.distortion / normalization_divisor;
                float bar_y_offset = base_height;
                m_graphene_rect = GRAPHENE_RECT_INIT(base_width, bar_y_offset - bar_height, width_step, bar_height);
                float border_thickness[] = {RESIDUE_BORDER_WIDTH,RESIDUE_BORDER_WIDTH,RESIDUE_BORDER_WIDTH,RESIDUE_BORDER_WIDTH};
                GskRoundedRect outline;
                gsk_rounded_rect_init_from_rect(
                    &outline,
                    &m_graphene_rect,
                    0
                );
                GdkRGBA residue_color_computed = residue_color;
                GdkRGBA border_color_computed = border_color;
                border_color_computed.red = 0.4 * residue.distortion / normalization_divisor;
                residue_color_computed.red = residue.distortion / normalization_divisor;
                GdkRGBA border_colors[] = {border_color_computed,border_color_computed,border_color_computed,border_color_computed};
                gtk_snapshot_append_color(snapshot, &residue_color_computed, &m_graphene_rect);
                gtk_snapshot_append_border(snapshot, &outline , border_thickness, border_colors);
                base_width += width_step;
            }
            base_height += CHAIN_SPACING/2.f + height_diff;
        }

        cairo_destroy(cairo_canvas);
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
            *minimum_size = max_chain_residues * (RESIDUE_WIDTH + RESIDUE_SPACING) + AXIS_MARGIN;
            *natural_size = max_chain_residues * (RESIDUE_WIDTH + RESIDUE_SPACING) + AXIS_MARGIN;
            break;
        }
        case GTK_ORIENTATION_VERTICAL:{
            auto num_of_chains = self->_vi->cviv.size();
            //g_debug("Num of chains: %u",num_of_chains);
            auto size = num_of_chains * (CHAIN_HEIGHT + CHAIN_SPACING) + TITLE_HEIGHT;
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