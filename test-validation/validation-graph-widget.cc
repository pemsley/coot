#include "validation-graph-widget.hh"
#include "residue-validation-information.hh"
#include "validation-information.hh"
#include <algorithm>
#include <memory>
#include <string>
#include <vector>
#include <cmath>

//typedef std::map<graphene_rect_t,const coot::residue_validation_information_t*, GrapheneRectCompare> coord_cache_t;
typedef std::vector<std::pair<graphene_rect_t,const coot::residue_validation_information_t*>> coord_cache_t;

static guint residue_clicked_signal;
struct _CootValidationGraph {
    GtkWidget parent;

    std::shared_ptr<const coot::validation_information_t> _vi;
    std::unique_ptr<coord_cache_t> coordinate_cache;
    /// Single-chain mode if not null
    std::unique_ptr<std::string> single_chain_id;
    float horizontal_scale;
};

/// Basis for max bar height
const int CHAIN_HEIGHT = 250;
/// Used for allocating space for axes and labels.
/// Fractions of this value can be used as a basis for computing other constants
/// in order to achieve desired layout proportions.
const int CHAIN_SPACING = 70;
const int RESIDUE_WIDTH = 9;
/// Breathing space for residue rectangle's borders
const int RESIDUE_SPACING = 3;
/// For drawing the main title
const int TITLE_HEIGHT = 30;
/// Space between the y-axis and the left-most bar in the graph
const int GRAPH_Y_AXIS_SEPARATION = 10;
/// Space reserved for the y-axis and its' labels. Axis is drawn at this X offset.
const int AXIS_MARGIN = 25;
/// Space left between the last residue bar and the right edge of the graph
const int RIGHT_SIDE_MARGIN = 10;
const double AXIS_LINE_WIDTH = 2;
const float RESIDUE_BORDER_WIDTH = 1;
const int MARKER_LENGTH = 3;
const unsigned int VERTICAL_MARKER_COUNT = 12;
const unsigned int HORIZONTAL_MARKER_INTERVAL = 10;

// COMPUTED VALUES:

const int GRAPH_HORIZ_OFFSET = AXIS_MARGIN + GRAPH_Y_AXIS_SEPARATION;
const float CHAIN_LABEL_VERT_OFFSET = CHAIN_SPACING * 1.f / 5.f;
/// Space between the x-axis and the bottom of the widget
const float BOTTOM_MARGIN = CHAIN_SPACING * 4.f / 15.f;
/// Space between the x-axis and the bottom of the graph (bars)
const float GRAPH_X_AXIS_SEPARATION = CHAIN_SPACING * 2.f / 15.f;
const float AXIS_HEIGHT = GRAPH_X_AXIS_SEPARATION + CHAIN_HEIGHT;
const float AXIS_VERT_OFFSET = CHAIN_LABEL_VERT_OFFSET + CHAIN_SPACING * 2.f / 5.f;
const float GRAPH_VERT_OFFSET = AXIS_VERT_OFFSET - GRAPH_X_AXIS_SEPARATION;
const int MARKER_VERT_PLACEMENT = AXIS_MARGIN - MARKER_LENGTH;

/// Returns the maximum number of residues in a chain (maximum among all the chains)
size_t max_chain_residue_count(CootValidationGraph* self) {
    using it_t = coot::chain_validation_information_t;
    auto biggest_chain = std::max_element(self->_vi->cviv.cbegin(),self->_vi->cviv.cend(),
        [](const it_t& lhs, const it_t& rhs){
            return lhs.rviv.size() < rhs.rviv.size();
    });
    if (biggest_chain != self->_vi->cviv.cend()) {
        return biggest_chain->rviv.size();
    } else {
        return 0;
    }
}

const coot::chain_validation_information_t* get_chain_with_id(CootValidationGraph* self,const std::string& chain_id) {
    auto ret = self->_vi->cviv.cend();
    ret = std::find_if(self->_vi->cviv.cbegin(), self->_vi->cviv.cend(), 
        [&](const coot::chain_validation_information_t& chain){
            //g_debug("cmp \"%s\" == \"%s\" gives %s",chain.chain_id.c_str(),chain_id.c_str(),chain.chain_id == chain_id ? "true" : "false");
            return chain.chain_id == chain_id;
        }
    );
    if(ret == self->_vi->cviv.cend()) {
        //g_debug("Chain with id \"%s\" not found!",chain_id.c_str());
        return nullptr;
    } else {
        const auto* retaddr = &*ret;
        //g_debug("Returning %p",retaddr);
        return retaddr;
    }
}

/// Returns maximum function value found for the given chain
double max_residue_function_value_for_chain(const std::vector<coot::residue_validation_information_t>& rviv) {
    using it_t = coot::residue_validation_information_t;
    auto iter = std::max_element(rviv.cbegin(),rviv.cend(),
    [](const it_t& lhs, const it_t& rhs){
        return lhs.function_value < rhs.function_value;
    });
    if(iter != rviv.end()) {
        return iter->function_value;
    } else {
        g_warning("Returning 0 as max value for an empty chain");
        return 0.0;
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

/// Computes the amplitude between the minimum and maximum (distortion) value representable on the graph
double compute_amplitude(coot::graph_data_type type, const std::vector<coot::residue_validation_information_t>& rviv) {
    using ty = coot::graph_data_type;
    switch (type) {
        case ty::Distortion: {
            return 100.f;
        }
        case ty::Probability:
        case ty::LogProbability: {
            return 1.f;
        }
        case ty::Correlation: {
            return 1.2f;
        }
        default: {
            // This assumes that by default, the lowest value is 0
            return max_residue_function_value_for_chain(rviv);
        }
    }
}

/// Compute the lowest expected value
inline double compute_floor_value(coot::graph_data_type type) {
    using ty = coot::graph_data_type;
    switch (type) {
        case ty::Correlation: {
            return -0.2f;
        }
        default: {
            return -0.f;
        }
    }
}

double map_value_to_bar_proportion(double residue_value, double amplitude, coot::graph_data_type type) {
    using ty = coot::graph_data_type;
    double floor = compute_floor_value(type);
    double base_proportion = (residue_value - floor) / amplitude;
    switch (type) {
        case ty::LogProbability: { 
            return std::log10(base_proportion * 9.f + 1.f);
        }
        default: {
            return base_proportion;
        }
    }
}

double map_bar_proportion_to_value(double bar_height_ratio, double amplitude, coot::graph_data_type type) {
    using ty = coot::graph_data_type;
    double floor = compute_floor_value(type);
    switch (type) {
        case ty::LogProbability: { 
            return (std::pow(10,bar_height_ratio) - 1.f) / 9.f * amplitude + floor;
            
        }
        default: {
            return bar_height_ratio * amplitude + floor;
        }
    }
}

G_BEGIN_DECLS

G_DEFINE_TYPE(CootValidationGraph, coot_validation_graph, GTK_TYPE_WIDGET)

// not sure what this is for or whether it is going to be needed at all
// struct _CootValidationGraphClass {
//     GObjectClass parent_class;
// };

void coot_validation_graph_snapshot (GtkWidget *widget, GtkSnapshot *snapshot)
{
    CootValidationGraph* self = COOT_COOT_VALIDATION_GRAPH(widget);
    self->coordinate_cache->clear();
    
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
    cairo_move_to(cairo_canvas,(w - layout_width) / 2.f,(TITLE_HEIGHT + layout_height) / 2.f);
    pango_cairo_show_layout(cairo_canvas, pango_layout);

    if(self->_vi) {
        if(self->single_chain_id) {
            if (! get_chain_with_id(self,*self->single_chain_id)) {
                return;
            }
        }

        // I can't get this to render the text where it needs to be, so I'm using cairo directly
        // gtk_snapshot_append_layout(snapshot,pango_layout,&attribute_color);
        // A GtkLabel as a child widget could also be used, but I have no idea how to manage layout inside widgets

        float base_height = TITLE_HEIGHT;
        const float _max_chain_residue_count = self->single_chain_id ? get_chain_with_id(self,*self->single_chain_id)->rviv.size() : max_chain_residue_count(self);
        float width_step = (w - (float) (GRAPH_HORIZ_OFFSET + RIGHT_SIDE_MARGIN)) / _max_chain_residue_count;
        const int chain_count =  self->single_chain_id ? 1 : self->_vi->cviv.size();
        float height_diff = 0;
        if (chain_count != 1) {
            height_diff = (h - (float) TITLE_HEIGHT - (float) chain_count * (CHAIN_HEIGHT + CHAIN_SPACING)) / ((float) chain_count - 1.f);
        }
        
        cairo_set_line_width(cairo_canvas,AXIS_LINE_WIDTH);

        for(const auto& chain: self->_vi->cviv) {
            if(self->single_chain_id) {
                if(*self->single_chain_id != chain.chain_id) {
                    continue;
                }
            }
            m_graphene_rect = GRAPHENE_RECT_INIT(0, 0, w, h);

            // Label chain
            std::string chain_markup = "<span size=\"medium\" weight=\"bold\">";
            // if(self->single_chain_id) {
            //     // we wanna actually show the type of the graph, not the name of the chain
            //     chain_markup += self->_vi->name;
            // } else { 
            chain_markup += "Chain " + chain.chain_id;
            //}
            chain_markup += "</span>";
            pango_layout_set_markup(pango_layout,chain_markup.c_str(),-1);
            pango_layout_get_pixel_size(pango_layout,&layout_width,&layout_height);
            cairo_move_to(cairo_canvas,0,base_height + CHAIN_LABEL_VERT_OFFSET);
            pango_cairo_show_layout(cairo_canvas, pango_layout);

            // Draw axes

            const float axis_y_offset = base_height + AXIS_VERT_OFFSET;

            // main vertical axis
            cairo_move_to(cairo_canvas, AXIS_MARGIN, axis_y_offset);
            cairo_line_to(cairo_canvas, AXIS_MARGIN, axis_y_offset + AXIS_HEIGHT);
            cairo_stroke(cairo_canvas);

            const double amplitude = compute_amplitude(self->_vi->type,chain.rviv);
            // vertical axis markers
            for(unsigned int m = 0; m <= VERTICAL_MARKER_COUNT; m++) {
                float marker_offset = m * CHAIN_HEIGHT / (float) VERTICAL_MARKER_COUNT;
                cairo_move_to(cairo_canvas, MARKER_VERT_PLACEMENT, axis_y_offset + marker_offset);
                cairo_line_to(cairo_canvas, AXIS_MARGIN, axis_y_offset + marker_offset);
                cairo_stroke(cairo_canvas);
                
                double marker_level;
                if (coot::should_hang_down(self->_vi->type)) {
                    marker_level = map_bar_proportion_to_value(m / (float) VERTICAL_MARKER_COUNT, amplitude, self->_vi->type);
                } else {
                    marker_level = map_bar_proportion_to_value(1 - m / (float) VERTICAL_MARKER_COUNT, amplitude, self->_vi->type);
                }
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
            unsigned int idx = 0;
            for(const auto& residue: chain.rviv) {
                /// draw bar

                float bar_height = CHAIN_HEIGHT * std::max(std::min(map_value_to_bar_proportion(residue.function_value, amplitude, self->_vi->type),1.0),0.0);
                float bar_y_offset = base_height;
                if(coot::should_hang_down(self->_vi->type)) {
                    m_graphene_rect = GRAPHENE_RECT_INIT(base_width, bar_y_offset - CHAIN_HEIGHT, RESIDUE_WIDTH * self->horizontal_scale, bar_height);
                } else {
                    m_graphene_rect = GRAPHENE_RECT_INIT(base_width, bar_y_offset - bar_height, RESIDUE_WIDTH * self->horizontal_scale, bar_height);
                }
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
                auto green_to_red = [&](double bar_proportion){
                    border_color_computed.red = 0.6 * bar_proportion;
                    border_color_computed.green = (1.f - bar_proportion) * residue_color.green;
                    //border_color_computed.blue = std::pow(0.9 * bar_proportion,5);
                    residue_color_computed.red = bar_proportion;
                    residue_color_computed.green = (1.f - std::pow(bar_proportion,3)) * residue_color.green;
                    //residue_color_computed.blue = std::pow(bar_proportion,5);

                };
                auto red_to_green = [&](double bar_proportion){
                    // dirty trick
                    green_to_red(1-bar_proportion);
                };
                switch (self->_vi->type) {
                    case coot::graph_data_type::LogProbability:
                    case coot::graph_data_type::Probability: {
                        red_to_green(map_value_to_bar_proportion(residue.function_value, amplitude, self->_vi->type));
                        break;
                    }
                    default: {
                        green_to_red(map_value_to_bar_proportion(residue.function_value, amplitude, self->_vi->type));
                    }
                }
                GdkRGBA border_colors[] = {border_color_computed,border_color_computed,border_color_computed,border_color_computed};
                gtk_snapshot_append_color(snapshot, &residue_color_computed, &m_graphene_rect);
                gtk_snapshot_append_border(snapshot, &outline , border_thickness, border_colors);

                // draw horizontal markers
                if(++idx % HORIZONTAL_MARKER_INTERVAL == 0) {
                    cairo_move_to(cairo_canvas,base_width + width_step/2.f, axis_y_offset + AXIS_HEIGHT);
                    cairo_line_to(cairo_canvas,base_width + width_step/2.f, axis_y_offset + AXIS_HEIGHT + MARKER_LENGTH * 2);
                    cairo_stroke(cairo_canvas);
                    std::string marker_label = "<span size=\"x-small\" >" + std::to_string(idx) + "</span>";
                    pango_layout_set_markup(pango_layout,marker_label.c_str(),-1);
                    pango_layout_get_pixel_size(pango_layout,&layout_width,&layout_height);
                    cairo_move_to(cairo_canvas, base_width + width_step/2.f - layout_width/2.f,axis_y_offset + AXIS_HEIGHT + MARKER_LENGTH * 2 );
                    pango_cairo_show_layout(cairo_canvas, pango_layout);
                }

                base_width += width_step;
            }
            base_height += GRAPH_X_AXIS_SEPARATION + BOTTOM_MARGIN + height_diff;
        }
    }
    g_object_unref(pango_layout);
    cairo_destroy(cairo_canvas);
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
            if (self->single_chain_id) {
                const auto* chain = get_chain_with_id(self, *self->single_chain_id);
                if(chain) {
                    *minimum_size = chain->rviv.size() * (RESIDUE_WIDTH + RESIDUE_SPACING) * self->horizontal_scale + GRAPH_HORIZ_OFFSET + RIGHT_SIDE_MARGIN;
                    *natural_size = chain->rviv.size() * (RESIDUE_WIDTH + RESIDUE_SPACING) * self->horizontal_scale + GRAPH_HORIZ_OFFSET + RIGHT_SIDE_MARGIN;
                }
            } else {
                auto max_chain_residues = max_chain_residue_count(self);
                *minimum_size = max_chain_residues * (RESIDUE_WIDTH + RESIDUE_SPACING) * self->horizontal_scale + GRAPH_HORIZ_OFFSET + RIGHT_SIDE_MARGIN;
                *natural_size = max_chain_residues * (RESIDUE_WIDTH + RESIDUE_SPACING) * self->horizontal_scale + GRAPH_HORIZ_OFFSET + RIGHT_SIDE_MARGIN;
            }
            break;
        }
        case GTK_ORIENTATION_VERTICAL:{
            auto num_of_chains = self->single_chain_id ? 1 : self->_vi->cviv.size();
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
        std::string tooltip_text = residue_ptr->label + ", value=" + std::to_string(residue_ptr->function_value).erase(5);
        gtk_tooltip_set_text(tooltip,tooltip_text.c_str());
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
    self->horizontal_scale = 1.f;
    self->single_chain_id = nullptr;

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
    self->single_chain_id.reset(nullptr);
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

void 
coot_validation_graph_set_horizontal_zoom_scale(CootValidationGraph* self, float scale) 
{
    self->horizontal_scale = scale;
    gtk_widget_queue_draw(GTK_WIDGET(self));
    gtk_widget_queue_resize(GTK_WIDGET(self));
}

void
coot_validation_graph_set_single_chain_mode(CootValidationGraph* self, const char* chain_id)
{
    if(chain_id) {
        self->single_chain_id.reset(new std::string(chain_id));
    } else {
        self->single_chain_id.reset(nullptr);
    }
    gtk_widget_queue_draw(GTK_WIDGET(self));
    gtk_widget_queue_resize(GTK_WIDGET(self));
}


G_END_DECLS

void coot_validation_graph_set_validation_information(CootValidationGraph* self, std::shared_ptr<coot::validation_information_t> vi) {
    // The stored pointers become invalidated
    self->coordinate_cache->clear();
    self->_vi = vi;
    gtk_widget_queue_draw(GTK_WIDGET(self));
    gtk_widget_queue_resize(GTK_WIDGET(self));
}



