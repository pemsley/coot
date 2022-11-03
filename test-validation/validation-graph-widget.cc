#include "validation-graph-widget.hh"
#include "residue-validation-information.hh"
#include "validation-information.hh"
#include <algorithm>
#include <memory>
#include <string>
#include <vector>

// class GrapheneRectCompare {
//     public:
//     bool operator()(const graphene_rect_t& lhs,const graphene_rect_t& rhs) const noexcept {
//         if (graphene_rect_get_x(&lhs) < graphene_rect_get_x(&rhs)) {
//             return true;
//         } else if (graphene_rect_get_y(&lhs) < graphene_rect_get_y(&rhs)) {
//             return true;
//         } else if (graphene_rect_get_width(&lhs) < graphene_rect_get_width(&rhs)) {
//             return true;
//         } else {
//             return (graphene_rect_get_height(&lhs) < graphene_rect_get_height(&rhs));
//         }
//         //return lhs < rhs;
//     }
// };

//typedef std::map<graphene_rect_t,const coot::residue_validation_information_t*, GrapheneRectCompare> coord_cache_t;
typedef std::vector<std::pair<graphene_rect_t,const coot::residue_validation_information_t*>> coord_cache_t;

static guint residue_clicked_signal;
struct _CootValidationGraph {
    GtkWidget parent;

    std::shared_ptr<const coot::validation_information_t> _vi;
    std::unique_ptr<coord_cache_t> coordinate_cache;
};

/// Basis for max bar height
const int CHAIN_HEIGHT = 120;
/// Used for allocating space for axes and labels
const int CHAIN_SPACING = 60;
const int RESIDUE_WIDTH = 9;
/// Breathing space for residue rectangle's borders
const int RESIDUE_SPACING = 1;
/// For drawing the main title
const int TITLE_HEIGHT = 30;
/// Space between the y-axis and the left-most bar in the graph
const int GRAPH_Y_AXIS_SEPARATION = 10;
/// Space reserved for the y-axis and its' labels. Axis is drawn at this X offset.
const int AXIS_MARGIN = 25;
const double AXIS_LINE_WIDTH = 2;
const float RESIDUE_BORDER_WIDTH = 1;
const int MARKER_LENGTH = 3;
const unsigned int VERTICAL_MARKER_COUNT = 4;

// COMPUTED VALUES:

const int GRAPH_HORIZ_OFFSET = AXIS_MARGIN + GRAPH_Y_AXIS_SEPARATION;
const float CHAIN_LABEL_VERT_OFFSET = CHAIN_SPACING * 2.f / 5.f;
/// Space between the x-axis and the bottom of the graph
const float GRAPH_X_AXIS_SEPARATION = CHAIN_SPACING / 5.f;
const float AXIS_HEIGHT = GRAPH_X_AXIS_SEPARATION + CHAIN_HEIGHT;
const float AXIS_VERT_OFFSET = CHAIN_SPACING * 4.f / 5.f;
const float GRAPH_VERT_OFFSET = AXIS_VERT_OFFSET - GRAPH_X_AXIS_SEPARATION;
const int MARKER_VERT_PLACEMENT = AXIS_MARGIN - MARKER_LENGTH;

size_t max_chain_residue_count(CootValidationGraph* self) {
    using it_t = coot::chain_validation_information_t;
    return std::max_element(self->_vi->cviv.cbegin(),self->_vi->cviv.cend(),
        [](const it_t& lhs, const it_t& rhs){
            return lhs.rviv.size() < rhs.rviv.size();
    })->rviv.size();
}

double max_chain_residue_distortion(const std::vector<coot::residue_validation_information_t>& rviv) {
    using it_t = coot::residue_validation_information_t;
    return std::max_element(rviv.cbegin(),rviv.cend(),
    [](const it_t& lhs, const it_t& rhs){
        return lhs.distortion < rhs.distortion;
    })->distortion;
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

    self->coordinate_cache->clear();
    if(self->_vi) {
        // 1. Draw title
        graphene_rect_t m_graphene_rect = GRAPHENE_RECT_INIT(0, 0, w, h);
        cairo_t* cairo_canvas = gtk_snapshot_append_cairo(snapshot,&m_graphene_rect);
        cairo_set_source_rgb(cairo_canvas, attribute_color.red, attribute_color.green, attribute_color.blue);
        
        // This does not respect GTK theming
        // PangoLayout* pango_layout = pango_cairo_create_layout(cairo_canvas);
        PangoLayout* pango_layout = pango_layout_new(gtk_widget_get_pango_context(widget));
        std::string title_markup = "<span size=\"large\" weight=\"bold\">" + self->_vi->name + "</span>";
        pango_layout_set_markup(pango_layout,title_markup.c_str(),-1);
        int layout_width, layout_height;
        pango_layout_get_pixel_size(pango_layout,&layout_width,&layout_height);
        cairo_move_to(cairo_canvas,(w - layout_width) / 2,(TITLE_HEIGHT + layout_height) / 2);
        pango_cairo_show_layout(cairo_canvas, pango_layout);

        // I can't get this to render the text where it needs to be, so I'm using cairo directly
        // gtk_snapshot_append_layout(snapshot,pango_layout,&attribute_color);
        // A GtkLabel as a child widget could also be used, but I have no idea how to manage layout inside widgets

        float base_height = TITLE_HEIGHT;
        float width_step = (w - (float) AXIS_MARGIN) / (float) max_chain_residue_count(self);
        const int chain_count = self->_vi->cviv.size();
        float height_diff = 0;
        if (chain_count != 1) {
            height_diff = (h - (float) TITLE_HEIGHT - (float) chain_count * (CHAIN_HEIGHT + CHAIN_SPACING)) / ((float) chain_count - 1.f);
        }
        
        cairo_set_line_width(cairo_canvas,AXIS_LINE_WIDTH);

        for(const auto& chain: self->_vi->cviv) {
            m_graphene_rect = GRAPHENE_RECT_INIT(0, 0, w, h);

            // Label chain
            std::string chain_markup = "<span size=\"medium\" weight=\"bold\">Chain " + chain.chain_id + "</span>";
            pango_layout_set_markup(pango_layout,chain_markup.c_str(),-1);
            pango_layout_get_pixel_size(pango_layout,&layout_width,&layout_height);
            cairo_move_to(cairo_canvas,0,base_height - layout_height / 2.f + CHAIN_LABEL_VERT_OFFSET);
            pango_cairo_show_layout(cairo_canvas, pango_layout);

            // Draw axes

            const float axis_y_offset = base_height + AXIS_VERT_OFFSET;

            // main vertical axis
            cairo_move_to(cairo_canvas, AXIS_MARGIN, axis_y_offset);
            cairo_line_to(cairo_canvas, AXIS_MARGIN, axis_y_offset + AXIS_HEIGHT);
            cairo_stroke(cairo_canvas);

            const double normalization_divisor = max_chain_residue_distortion(chain.rviv);
            // vertical axis markers
            for(unsigned int m = 0; m <= VERTICAL_MARKER_COUNT; m++) {
                float marker_offset = m * CHAIN_HEIGHT / (float) VERTICAL_MARKER_COUNT;
                cairo_move_to(cairo_canvas, MARKER_VERT_PLACEMENT, axis_y_offset + marker_offset);
                cairo_line_to(cairo_canvas, AXIS_MARGIN, axis_y_offset + marker_offset);
                cairo_stroke(cairo_canvas);
                
                double marker_level = (1 - m / (float) VERTICAL_MARKER_COUNT) * normalization_divisor;
                std::string marker_label = "<span size=\"x-small\" >" + std::to_string(marker_level).erase(4) + "</span>";
                pango_layout_set_markup(pango_layout,marker_label.c_str(),-1);
                pango_layout_get_pixel_size(pango_layout,&layout_width,&layout_height);
                cairo_move_to(cairo_canvas, 0,axis_y_offset + marker_offset - layout_height/2.f);
                pango_cairo_show_layout(cairo_canvas, pango_layout);
            }
            
            // horizontal axis
            cairo_move_to(cairo_canvas, AXIS_MARGIN, axis_y_offset + AXIS_HEIGHT);
            cairo_line_to(cairo_canvas, w, axis_y_offset + AXIS_HEIGHT);
            cairo_stroke(cairo_canvas);

            base_height += AXIS_HEIGHT + GRAPH_VERT_OFFSET;
            float base_width = GRAPH_HORIZ_OFFSET;
            for(const auto& residue: chain.rviv) {
                float bar_height = CHAIN_HEIGHT * residue.distortion / normalization_divisor;
                float bar_y_offset = base_height;
                m_graphene_rect = GRAPHENE_RECT_INIT(base_width, bar_y_offset - bar_height, RESIDUE_WIDTH, bar_height);
                //self->coordinate_cache->operator[](m_graphene_rect) = &residue;
                self->coordinate_cache->push_back(std::pair<graphene_rect_t,const coot::residue_validation_information_t*>{m_graphene_rect,&residue});
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
            base_height += GRAPH_X_AXIS_SEPARATION + height_diff;
        }
        g_object_unref(pango_layout);
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
            *minimum_size = max_chain_residues * (RESIDUE_WIDTH + RESIDUE_SPACING) + GRAPH_HORIZ_OFFSET;
            *natural_size = max_chain_residues * (RESIDUE_WIDTH + RESIDUE_SPACING) + GRAPH_HORIZ_OFFSET;
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

inline coord_cache_t::const_iterator residue_from_coords(CootValidationGraph* self, gdouble x, gdouble y) {

    graphene_point_t point;
    graphene_point_init(&point,x,y);
    coord_cache_t::const_iterator clicked = std::find_if(self->coordinate_cache->cbegin(),self->coordinate_cache->cend(),[point](const std::pair<graphene_rect_t,const coot::residue_validation_information_t*>& it){
        return graphene_rect_contains_point(&it.first,&point);
    });
    return clicked;
}

// static void on_hover (
//     GtkEventControllerMotion* hover_controller,
//     gdouble x,
//     gdouble y,
//     gpointer user_data
// ) {

//     CootValidationGraph* self = COOT_COOT_VALIDATION_GRAPH(user_data);
//     coord_cache_t::const_iterator hovered = residue_from_coords(self,x,y);
//     if(self->coordinate_cache->cend() != hovered) {
//         const auto* residue_ptr = hovered->second;
//         g_debug("Hover over residue: %s, at x: %f, y: %f",residue_ptr->label.c_str(),x,y);
//     }
// }

static void on_left_click (
  GtkGestureClick* gesture_click,
  gint n_press,
  gdouble x,
  gdouble y,
  gpointer user_data
) {
    CootValidationGraph* self = COOT_COOT_VALIDATION_GRAPH(user_data);
    //g_debug("On click at widget: %p, at x: %f, y: %f",self,x,y);
    coord_cache_t::const_iterator clicked = residue_from_coords(self,x,y);
    if(self->coordinate_cache->cend() != clicked) {
        const auto* residue_ptr = clicked->second;
        gtk_gesture_set_state(GTK_GESTURE(gesture_click),GTK_EVENT_SEQUENCE_CLAIMED);
        g_signal_emit(self,residue_clicked_signal,0,residue_ptr);
    } else {
        gtk_gesture_set_state(GTK_GESTURE(gesture_click),GTK_EVENT_SEQUENCE_NONE);
    }
}

gboolean query_tooltip (
    CootValidationGraph* self,
    gint x,
    gint y,
    gboolean keyboard_mode,
    GtkTooltip* tooltip,
    gpointer user_data
) {
    coord_cache_t::const_iterator hovered = residue_from_coords(self,x,y);
    if(self->coordinate_cache->cend() != hovered) {
        const auto* residue_ptr = hovered->second;
        // g_debug("Hover over residue: %s, at x: %f, y: %f",residue_ptr->label.c_str(),x,y);
        gtk_tooltip_set_text(tooltip,residue_ptr->label.c_str());
        // todo: remove magic numbers
        GdkRectangle rect = {x,y - 20,100,100};
        gtk_tooltip_set_tip_area(tooltip,&rect);
        return TRUE;
    } else
        return FALSE;
}

static void coot_validation_graph_init(CootValidationGraph* self) {
    // I think that this is the primary constructor
    gtk_widget_set_has_tooltip(GTK_WIDGET(self),TRUE);
    g_signal_connect(self,"query-tooltip",G_CALLBACK(query_tooltip),NULL);

    // I don't know how g_object_new initializes C++ stuff. Better set those up manually
    self->_vi = std::shared_ptr<const coot::validation_information_t>(nullptr);
    self->coordinate_cache = std::make_unique<coord_cache_t>();

    GtkGesture* click_controller = gtk_gesture_click_new();
    // GtkEventController* hover_controller = gtk_event_controller_motion_new();

    // left mouse button
    gtk_gesture_single_set_button(GTK_GESTURE_SINGLE(click_controller),GDK_BUTTON_PRIMARY);
    g_signal_connect(click_controller,"pressed",G_CALLBACK(on_left_click),self);

    // g_signal_connect(hover_controller,"motion",G_CALLBACK(on_hover),self);

    gtk_widget_add_controller(GTK_WIDGET(self),GTK_EVENT_CONTROLLER(click_controller));
    // gtk_widget_add_controller(GTK_WIDGET(self),GTK_EVENT_CONTROLLER(hover_controller));
}

static void coot_validation_graph_dispose(GObject* _self) {

    CootValidationGraph* self = COOT_COOT_VALIDATION_GRAPH(_self);
    self->_vi.reset();
    self->coordinate_cache.reset(nullptr);
    G_OBJECT_CLASS(coot_validation_graph_parent_class)->dispose(_self);
}

static void coot_validation_graph_class_init(CootValidationGraphClass* klass) {
    // I think that this is a GObject class constructor that sets up the GObject class at runtime.
    residue_clicked_signal = g_signal_new("residue-clicked",
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

void coot_validation_graph_set_validation_information(CootValidationGraph* self, std::shared_ptr<coot::validation_information_t> vi) {
    // The stored pointers become invalidated
    self->coordinate_cache->clear();
    self->_vi = vi;
}
