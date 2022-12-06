#include "text-icon-button.hh"
#include <cstring>

struct _CootTextIconButton {
    GtkButton parent;

    GtkWidget* container_box;
    GtkWidget* label;
    GtkWidget* icon;
};

G_BEGIN_DECLS

G_DEFINE_TYPE(CootTextIconButton, coot_text_icon_button, GTK_TYPE_BUTTON)

const char* widget_template = " \
<interface> \
  <template class=\"CootTextIconButton\" parent=\"GtkButton\"> \
    <child> \
      <object class=\"GtkBox\" id=\"container_box\"> \
        <property name=\"orientation\">horizontal</property> \
        <child> \
          <object class=\"GtkLabel\" id=\"label\"> \
          </object> \
        </child> \
        <child> \
          <object class=\"GtkImage\" id=\"icon\"> \
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
   
}

static void coot_text_icon_button_dispose(GObject* _self) {
    gtk_widget_dispose_template (GTK_WIDGET (_self), COOT_TEXT_ICON_BUTTON_TYPE);
    G_OBJECT_CLASS(coot_text_icon_button_parent_class)->dispose(_self);
}

static void coot_text_icon_button_class_init(CootTextIconButtonClass* klass) {
    GtkWidgetClass *widget_class = GTK_WIDGET_CLASS (klass);
    GBytes* template_bytes = g_bytes_new_static(widget_template, strlen(widget_template));
    gtk_widget_class_set_template(widget_class, template_bytes);

    gtk_widget_class_bind_template_child (widget_class, CootTextIconButton, container_box);
    gtk_widget_class_bind_template_child (widget_class, CootTextIconButton, label);
    gtk_widget_class_bind_template_child (widget_class, CootTextIconButton, icon);
    G_OBJECT_CLASS(klass)->dispose = coot_text_icon_button_dispose;
    
}

CootTextIconButton* 
coot_text_icon_button_new()
{
    return COOT_COOT_TEXT_ICON_BUTTON(g_object_new (COOT_TEXT_ICON_BUTTON_TYPE, NULL));
}


G_END_DECLS
