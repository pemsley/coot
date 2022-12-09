#include "text-icon-button.hh"
#include <cstring>

static guint clicked_signal;

struct _CootTextIconButton {
    GtkWidget parent;

    GtkWidget* internal_button;
    GtkWidget* container_box;
    GtkWidget* label;
    GtkWidget* icon;
    bool show_label;
    bool show_icon;
    gchar* label_text;
    gchar* icon_name;
};

G_BEGIN_DECLS

G_DEFINE_TYPE(CootTextIconButton, coot_text_icon_button, GTK_TYPE_WIDGET)

enum {
    PROP_0,
    PROP_LABEL,
    PROP_ICON,
    PROP_SHOW_LABEL,
    PROP_SHOW_ICON,
    NUM_PROPERTIES
};
static GParamSpec *label_property = NULL;
static GParamSpec *icon_name_property = NULL;
static GParamSpec *show_label_property = NULL;
static GParamSpec *show_icon_property = NULL;


const char* widget_template = " \
<interface> \
  <template class=\"CootTextIconButton\" parent=\"GtkWidget\"> \
    <child> \
      <object class=\"GtkButton\" id=\"internal_button\"> \
        <child> \
          <object class=\"GtkBox\" id=\"container_box\"> \
            <property name=\"orientation\">horizontal</property> \
            <property name=\"spacing\">5</property> \
            <child> \
              <object class=\"GtkImage\" id=\"icon\"> \
              </object> \
            </child> \
            <child> \
              <object class=\"GtkLabel\" id=\"label\"> \
              </object> \
            </child> \
          </object> \
        </child> \
      </object> \
    </child> \
  </template> \
</interface> \
";

static void coot_text_icon_button_init(CootTextIconButton* self) {
    // I think that this is the primary constructor
    gtk_widget_init_template (GTK_WIDGET (self));

    g_signal_connect (G_OBJECT (self), "notify::label", G_CALLBACK (+[](GObject* _self, GParamSpec* pspec, gpointer user_data){
        CootTextIconButton* self = COOT_TEXT_ICON_BUTTON(_self);
        gtk_label_set_text(GTK_LABEL(self->label), self->label_text);
    }), NULL);
    g_signal_connect (G_OBJECT (self), "notify::icon-name", G_CALLBACK (+[](GObject* _self, GParamSpec* pspec, gpointer user_data){
        CootTextIconButton* self = COOT_TEXT_ICON_BUTTON(_self);
        gtk_image_set_from_icon_name(GTK_IMAGE(self->icon), self->icon_name);
    }), NULL);
    g_signal_connect (G_OBJECT (self), "notify::show-icon", G_CALLBACK (+[](GObject* _self, GParamSpec* pspec, gpointer user_data){
        CootTextIconButton* self = COOT_TEXT_ICON_BUTTON(_self);
        gtk_widget_set_visible(self->label, self->show_label);
    }), NULL);
    g_signal_connect (G_OBJECT (self), "notify::show-label", G_CALLBACK (+[](GObject* _self, GParamSpec* pspec, gpointer user_data){
        CootTextIconButton* self = COOT_TEXT_ICON_BUTTON(_self);
        gtk_widget_set_visible(self->icon, self->show_icon);
    }), NULL);
    g_signal_connect(self->internal_button, "clicked", G_CALLBACK(+[](GtkButton *button, gpointer user_data){
        CootTextIconButton* self = COOT_TEXT_ICON_BUTTON(user_data);
        g_signal_emit(self,clicked_signal,0,NULL);
    }), self);
}

static void coot_text_icon_button_dispose(GObject* _self) {
    gtk_widget_dispose_template (GTK_WIDGET (_self), COOT_TEXT_ICON_BUTTON_TYPE);
    CootTextIconButton* self = COOT_TEXT_ICON_BUTTON(_self);
    g_free(self->icon_name);
    g_free(self->label_text);
    G_OBJECT_CLASS(coot_text_icon_button_parent_class)->dispose(_self);
}

static void coot_text_icon_button_set_property(GObject        *object,
                                         guint           property_id,
                                         const GValue   *value,
                                         GParamSpec     *pspec) {
    CootTextIconButton *self = COOT_TEXT_ICON_BUTTON (object);
    switch (property_id) {
        case PROP_ICON:{
            //self->icon_name = g_value_get_string(value);
            break;
        }
        case PROP_SHOW_ICON:{
            self->show_icon = g_value_get_boolean(value);
            break;
        }
        case PROP_LABEL:{
            //self->label_text = g_value_get_string (value);
            break;
        }
        case PROP_SHOW_LABEL:{
            self->show_label = g_value_get_boolean(value);
            break;
        }
        default: {
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
        }
    }
}

static void coot_text_icon_button_get_property(GObject        *object,
                                         guint           property_id,
                                         GValue         *value,
                                         GParamSpec     *pspec) {
    CootTextIconButton *self = COOT_TEXT_ICON_BUTTON (object);
    switch (property_id) {
        case PROP_ICON:{
            g_value_set_string (value, self->icon_name);
            break;
        }
        case PROP_SHOW_ICON:{
            g_value_set_boolean(value, self->show_icon);
            break;
        }
        case PROP_LABEL:{
            g_value_set_string (value, self->label_text);
            break;
        }
        case PROP_SHOW_LABEL:{
            g_value_set_boolean(value, self->show_label);
            break;
        }
        default: {
            G_OBJECT_WARN_INVALID_PROPERTY_ID (object, property_id, pspec);
        }
    }
}

static void coot_text_icon_button_measure(
                                        GtkWidget* _self, 
                                        GtkOrientation orientation, 
                                        int for_size,
                                        int *minimum_size,
                                        int *natural_size,
                                        int *minimum_baseline,
                                        int *natural_baseline) {
    CootTextIconButton *self = COOT_TEXT_ICON_BUTTON (_self);
    gtk_widget_measure(self->container_box,orientation,for_size,minimum_size,natural_size,minimum_baseline,natural_baseline);
    // switch (orientation) {
    //     case GTK_ORIENTATION_HORIZONTAL:{
            
    //         break;
    //     }
    //     case GTK_ORIENTATION_VERTICAL:{
            
    //         break;
    //     }
    //     default:{
    //         break;
    //     }
    // }
}

// static void
// coot_text_icon_button_size_allocate (
//   GtkWidget* widget,
//   int a,
//   int b,
//   int baseline
// ) {
//      CootTextIconButton *self = COOT_TEXT_ICON_BUTTON (widget);
//      gtk_widget_size_allocate(self->internal_button, *allocation, baseline);
// }

static void coot_text_icon_button_class_init(CootTextIconButtonClass* klass) {
    GtkWidgetClass *widget_class = GTK_WIDGET_CLASS (klass);
    GBytes* template_bytes = g_bytes_new_static(widget_template, strlen(widget_template));
    gtk_widget_class_set_template(widget_class, template_bytes);

    gtk_widget_class_bind_template_child (widget_class, CootTextIconButton, internal_button);
    gtk_widget_class_bind_template_child (widget_class, CootTextIconButton, container_box);
    gtk_widget_class_bind_template_child (widget_class, CootTextIconButton, label);
    gtk_widget_class_bind_template_child (widget_class, CootTextIconButton, icon);

    widget_class->measure = coot_text_icon_button_measure;
    // widget_class->size_allocate = coot_text_icon_button_size_allocate;
    GObjectClass *gobject_class = G_OBJECT_CLASS (klass);
    gobject_class->set_property = coot_text_icon_button_set_property;
    gobject_class->get_property = coot_text_icon_button_get_property;
    gobject_class->dispose = coot_text_icon_button_dispose;

    label_property = g_param_spec_string("label", "label", "Displayed label", "e", G_PARAM_READWRITE);
    icon_name_property = g_param_spec_string("icon-name", "icon-name", "Icon name", "", G_PARAM_READWRITE);
    show_icon_property = g_param_spec_boolean("show-icon", "show-icon", "Controls visibility of icon", TRUE, G_PARAM_READWRITE);
    show_label_property = g_param_spec_boolean("show-label", "show-label", "Controls visibility of label", TRUE, G_PARAM_READWRITE);

    g_object_class_install_property (gobject_class, PROP_LABEL, label_property);
    g_object_class_install_property (gobject_class, PROP_ICON, icon_name_property);
    g_object_class_install_property (gobject_class, PROP_SHOW_LABEL, show_label_property);
    g_object_class_install_property (gobject_class, PROP_SHOW_ICON, show_icon_property);

    clicked_signal = g_signal_new("clicked",
        G_TYPE_FROM_CLASS (klass),
        (GSignalFlags) (G_SIGNAL_RUN_LAST | G_SIGNAL_NO_RECURSE | G_SIGNAL_NO_HOOKS),
        0 /* class offset.Subclass cannot override the class handler (default handler). */,
        NULL /* accumulator */,
        NULL /* accumulator data */,
        NULL /* C marshaller. g_cclosure_marshal_generic() will be used */,
        G_TYPE_NONE /* return_type */,
        0     /* n_params */
        // G PARAMETER TYPES
    );
    
}

CootTextIconButton* 
coot_text_icon_button_new()
{
    return COOT_TEXT_ICON_BUTTON(g_object_new (COOT_TEXT_ICON_BUTTON_TYPE, NULL));
}


G_END_DECLS
