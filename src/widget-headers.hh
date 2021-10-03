
#include <vector>
#include <string>
#include <utility>
#include <gtk/gtk.h>

std::string pre_directory_file_selection(GtkWidget *sort_button);
void filelist_into_fileselection_clist(GtkWidget *fileselection, const std::vector<std::string> &v);

GtkWidget *wrapped_nothing_bad_dialog(const std::string &label);

std::pair<short int, float> float_from_entry(GtkWidget *entry);
std::pair<short int, int>   int_from_entry(GtkWidget *entry);

GtkWidget *
add_validation_mol_menu_item(int imol, const std::string &name, GtkWidget *menu, GCallback callback);
void create_initial_validation_graph_submenu_generic(GtkWidget *window1,
                                                     const std::string &menu_name,
                                                     const std::string &sub_menu_name);



// To be used to (typically) get the menu item text label from chain
// option menus (rather than the ugly/broken casting of
// GtkPositionType data.  A wrapper to a static graphics_info_t
// function.
std::string menu_item_label(GtkWidget *menu_item);

void
add_map_colour_mol_menu_item(int imol, const std::string &name,
                                  GtkWidget *sub_menu, GCallback callback);

void add_map_scroll_wheel_mol_menu_item(int imol, 
                                        const std::string &name,
                                        GtkWidget *menu, 
                                        GCallback callback);

