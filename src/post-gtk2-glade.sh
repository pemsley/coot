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
    -e 's/ gtk_combo_box_append_text .GTK_COMBO_BOX / gtk_combo_box_text_append_text (GTK_COMBO_BOX_TEXT /' \
    -e 's/ gtk_combo_box_new_text / gtk_combo_box_text_new /' \
    -e 's?  GtkTooltips *?  // GtkTooltips ?' \
    -e 's? tooltips = gtk_tooltips_new ? // tooltips = gtk_tooltips_new ?' \
    -e 's? gtk_tooltips_set_tip ? // gtk_tooltips_set_tip ?' \
    -e 's? gtk_tool_item_set_tooltip ? // gtk_tool_item_set_tooltip ?' \
    -e 's? gtk_dialog_set_has_separator ? // gtk_dialog_set_has_separator ?' \
    -e 's? gtk_toolbar_set_orientation (GTK_TOOLBAR ? gtk_orientable_set_orientation (GTK_ORIENTABLE ?' \
    -e 's/gtk_widget_ref (widget), (GDestroyNotify) gtk_widget_unref)/g_object_ref (widget), (GDestroyNotify) g_object_unref)/' \
    -e 's/ gtk_about_dialog_set_name / gtk_about_dialog_set_program_name /' \
    -e 's/GDK_F7/GDK_KEY_F7/' \
    -e 's/GDK_F6/GDK_KEY_F6/' \
    -e 's/GDK_D/GDK_KEY_D/'   \
    -e 's/GDK_U/GDK_KEY_U/'   \
    -e 's/ gtk_vbox_new .FALSE/ gtk_box_new (GTK_ORIENTATION_VERTICAL/'   \
    -e 's/ gtk_vbox_new .TRUE/  gtk_box_new (GTK_ORIENTATION_VERTICAL/'   \
    -e 's/ gtk_hbox_new .FALSE/ gtk_box_new (GTK_ORIENTATION_HORIZONTAL/' \
    -e 's/ gtk_hbox_new .TRUE/  gtk_box_new (GTK_ORIENTATION_HORIZONTAL/' \
    -e 's? GLADE_HOOKUP_OBJECT_NO_REF .*tooltip? // GLADE_HOOKUP_OBJECT_NO_REF tooltip thing?' \
    gtk2-interface.c \
    | awk '
/ = GTK_DIALOG/ {f=$4; gsub("[(]", "", f); gsub("[)].*", "", f); print(" ", $1, "=", "gtk_dialog_get_content_area(", f, ");")}
/ GTK_WIDGET_SET_FLAGS .* GTK_CAN_DEFAULT/ {print "  gtk_widget_set_can_default " $2, "1);" ; next }
/ GTK_WIDGET_SET_FLAGS .* GTK_CAN_FOCUS/   {print "  gtk_widget_set_can_focus " $2, "1);" ; next }
/ GTK_WIDGET_UNSET_FLAGS .* GTK_CAN_FOCUS/   {print "  gtk_widget_set_can_focus " $2, "0);" ; next }
$0 !~ "= GTK_DIALOG"
    ' \
    > gtk2-interface.post-sed

sh fixup-gtk2-interface.sh gtk2-interface.post-sed
bash fixup-interface.h.sh

#    -e 's? GTK_WIDGET_SET_FLAGS ? // GTK_WIDGET_SET_FLAGS ?'
#   -e 's? GTK_WIDGET_UNSET_FLAGS ? // GTK_WIDGET_UNSET_FLAGS ?' \
#    -e 's? GTK_WIDGET_SET_FLAGS ? gtk_widget_set_can_default ?' \
 
# / GTK_WIDGET_SET_FLAGS .* GTK_CAN_DEFAULT / {print "  gtk_widget_set_can_default " $2, "1);" }

#    -e 's? GLADE_HOOKUP_OBJECT_NO_REF ? // GLADE_HOOKUP_OBJECT_NO_REF tooltip thing?'

# We don't care about old gtks now
#    -e 's/tmp_image = .*rtz.svg/#ifdef GTK_TYPE_MENU_TOOL_BUTTON\n  &/'
#    -e 's/set_tooltip .*model_toolbar_rot_trans_toolbutton.*;/&\n#endif\n/' 
#    -e 's/ *GLADE_HOOKUP_OBJECT .*model_toolbar_rot_trans_toolbutton.*/#ifdef GTK_TYPE_MENU_TOOL_BUTTON\n  &\n#endif/' 

