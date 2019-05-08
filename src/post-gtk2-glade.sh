# ok we leave interface.h alone and put the functions in add-c-extern.sh
# we only manipulate gtk2-interface.c

SED=sed

if [ $(uname) = Darwin ] ; then
   SED=gsed
fi

$SED -e 's/#include "callbacks.h.gtk2"/#include "callbacks.h"/' \
    -e 's/interface.h.gtk2/interface.h/' \
    -e 's/support.h.gtk2/support.h/' \
    -e 's/#include <unistd.h>/#ifndef _MSC_VER\n#include <unistd.h>\n#endif/' \
    -e '/gtk_about_dialog_set_comments/s/_(/_(g_strconcat(/' \
    -e '/gtk_about_dialog_set_comments/s/));/, NULL)));/' \
    -e '
        /gtk_image_menu_item_new/,/gtk_image_menu_item_set_image/ {
        /svg/s/create_pixmap (window1, /gtk_image_new_from_stock (/
        /svg/s/);/, GTK_ICON_SIZE_MENU);/
        /png/s/create_pixmap (window1, /gtk_image_new_from_stock (/
        /png/s/);/, GTK_ICON_SIZE_MENU);/
        }' \
    -e '
        /model_toolbar = gtk_toolbar_new ();/,/model_toolbar_menubar1 = gtk_menu_bar_new ();/ {
        /svg/s/create_pixmap (window1, /gtk_image_new_from_stock (/
        /svg/s/);/, tmp_toolbar_icon_size);/
        /png/s/create_pixmap (window1, /gtk_image_new_from_stock (/
        /png/s/);/, tmp_toolbar_icon_size);/
        }' \
    -e '
        /create_model_refine_dialog (void)/,/}/ {
        /svg/s/create_pixmap (model_refine_dialog, /gtk_image_new_from_stock (/
        /svg/s/);/, GTK_ICON_SIZE_BUTTON);/
        /png/s/create_pixmap (model_refine_dialog, /gtk_image_new_from_stock (/
        /png/s/);/, GTK_ICON_SIZE_BUTTON);/
        }' \
    -e '
        /create_other_model_tools_dialog (void)/,/}/ {
        /svg/s/create_pixmap (other_model_tools_dialog, /gtk_image_new_from_stock (/
        /svg/s/);/, GTK_ICON_SIZE_BUTTON);/
        /png/s/create_pixmap (other_model_tools_dialog, /gtk_image_new_from_stock (/
        /png/s/);/, GTK_ICON_SIZE_BUTTON);/
        }' \
    -e '
        /main_toolbar = gtk_toolbar_new ();/,/model_toolbar = gtk_toolbar_new ();/ {
        /svg/s/create_pixmap (window1, /gtk_image_new_from_stock (/
        /svg/s/);/, tmp_toolbar_icon_size);/
        /png/s/create_pixmap (window1, /gtk_image_new_from_stock (/
        /png/s/);/, tmp_toolbar_icon_size);/
        }' \
    -e '
        /toolbar1 = gtk_toolbar_new ();/,/accept_reject_dialog_frame_docked = gtk_frame_new (NULL);/ {
        /svg/s/create_pixmap (window1, /gtk_image_new_from_stock (/
        /svg/s/);/, tmp_toolbar_icon_size);/
        /png/s/create_pixmap (window1, /gtk_image_new_from_stock (/
        /png/s/);/, tmp_toolbar_icon_size);/
        }' \
    -e '
        /create_preferences (void)/,/}/ {
        /svg/s/create_pixmap (preferences, /gtk_image_new_from_stock (/
        /svg/s/);/, GTK_ICON_SIZE_BUTTON);/
        /png/s/create_pixmap (preferences, /gtk_image_new_from_stock (/
        /png/s/);/, GTK_ICON_SIZE_BUTTON);/
        }' \
    -e 's/tmp_image = .*rtz.svg/#ifdef GTK_TYPE_MENU_TOOL_BUTTON\n  &/' \
    -e 's/set_tooltip .*model_toolbar_rot_trans_toolbutton.*;/&\n#endif\n/' \
    -e 's/ *GLADE_HOOKUP_OBJECT .*model_toolbar_rot_trans_toolbutton.*/#ifdef GTK_TYPE_MENU_TOOL_BUTTON\n  &\n#endif/' \
    gtk2-interface.c > gtk2-interface.post-sed

cp gtk2-interface.post-sed gtk2-interface.c
sh fixup-gtk2-interface.sh
bash fixup-interface.h.sh

cp gtk2-interface.c gtk2-interface.c-orig

