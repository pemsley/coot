#ifndef TEXT_ICON_BUTTON_HH
#define TEXT_ICON_BUTTON_HH
#include <gtk/gtk.h>


G_BEGIN_DECLS   

#define COOT_TEXT_ICON_BUTTON_TYPE (coot_text_icon_button_get_type ())
G_DECLARE_FINAL_TYPE  (CootTextIconButton, coot_text_icon_button, COOT, TEXT_ICON_BUTTON, GtkWidget)


CootTextIconButton *coot_text_icon_button_new();

// C methods go here

G_END_DECLS

// C++ methods go here


#endif // TEXT_ICON_BUTTON_HH