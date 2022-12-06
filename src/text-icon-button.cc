#include "text-icon-button.hh"

struct _CootTextIconButton {
    GtkButton parent;

};

G_BEGIN_DECLS

G_DEFINE_TYPE(CootTextIconButton, coot_text_icon_button, GTK_TYPE_BUTTON)


static void coot_text_icon_button_init(CootTextIconButton* self) {
    // I think that this is the primary constructor
   
}

static void coot_text_icon_button_dispose(GObject* _self) {

    G_OBJECT_CLASS(coot_text_icon_button_parent_class)->dispose(_self);
}

static void coot_text_icon_button_class_init(CootTextIconButtonClass* klass) {
    G_OBJECT_CLASS(klass)->dispose = coot_text_icon_button_dispose;
    
}

CootTextIconButton* 
coot_text_icon_button_new()
{
    return COOT_COOT_TEXT_ICON_BUTTON(g_object_new (COOT_TEXT_ICON_BUTTON_TYPE, NULL));
}


G_END_DECLS
