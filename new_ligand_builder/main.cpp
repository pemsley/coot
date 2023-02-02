#include "ligand-builder.hpp"
#include "ligand_editor_canvas.hpp"
#include <gtk/gtk.h>



void build_main_window(GtkWindow* win) {
    GtkWidget* mainbox = gtk_box_new(GTK_ORIENTATION_VERTICAL,0);
    // Top toolbar
    GtkWidget* top_toolbar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL,5);
    gtk_box_append(GTK_BOX(mainbox), top_toolbar);

    GtkWidget* undo_button = gtk_button_new_with_label("Undo");
    gtk_box_append(GTK_BOX(top_toolbar), undo_button);
    GtkWidget* redo_button = gtk_button_new_with_label("Redo");
    gtk_box_append(GTK_BOX(top_toolbar), redo_button);
    GtkWidget* single_bond_button = gtk_button_new_with_label("Single Bond");
    gtk_box_append(GTK_BOX(top_toolbar), single_bond_button);
    GtkWidget* double_bond_button = gtk_button_new_with_label("Double Bond");
    gtk_box_append(GTK_BOX(top_toolbar), double_bond_button);
    GtkWidget* triple_bond_button = gtk_button_new_with_label("Triple Bond");
    gtk_box_append(GTK_BOX(top_toolbar), triple_bond_button);
    GtkWidget* stereo_out_modifier_button = gtk_button_new_with_label("Stereo Out Tool");
    gtk_box_append(GTK_BOX(top_toolbar), stereo_out_modifier_button);
    GtkWidget* charge_modifier_button = gtk_button_new_with_label("Charge Tool");
    gtk_box_append(GTK_BOX(top_toolbar), charge_modifier_button);
    GtkWidget* delete_button = gtk_button_new_with_label("Delete");
    gtk_box_append(GTK_BOX(top_toolbar), delete_button);
    GtkWidget* delete_hydrogens_button = gtk_button_new_with_label("Delete Hydrogens");
    gtk_box_append(GTK_BOX(top_toolbar), delete_hydrogens_button);
    GtkWidget* smiles_button = gtk_button_new_with_label("SMILES");
    gtk_box_append(GTK_BOX(top_toolbar), smiles_button);
    GtkWidget* format_button = gtk_button_new_with_label("Format");
    gtk_box_append(GTK_BOX(top_toolbar), format_button);
    GtkWidget* info_button = gtk_button_new_with_label("Info");
    gtk_box_append(GTK_BOX(top_toolbar), info_button);
    
    // Carbon ring picker
    GtkWidget* carbon_ring_picker = gtk_box_new(GTK_ORIENTATION_HORIZONTAL,5);
    gtk_box_append(GTK_BOX(mainbox), carbon_ring_picker);

    GtkWidget* buttom_3C = gtk_button_new_with_label("3-C");
    gtk_box_append(GTK_BOX(carbon_ring_picker), buttom_3C);
    GtkWidget* buttom_4C = gtk_button_new_with_label("4-C");
    gtk_box_append(GTK_BOX(carbon_ring_picker), buttom_4C);
    GtkWidget* buttom_5C = gtk_button_new_with_label("5-C");
    gtk_box_append(GTK_BOX(carbon_ring_picker), buttom_5C);
    GtkWidget* buttom_6C = gtk_button_new_with_label("6-C");
    gtk_box_append(GTK_BOX(carbon_ring_picker), buttom_6C);
    GtkWidget* buttom_6Arom = gtk_button_new_with_label("6-Arom");
    gtk_box_append(GTK_BOX(carbon_ring_picker), buttom_6Arom);
    GtkWidget* buttom_7C = gtk_button_new_with_label("7-C");
    gtk_box_append(GTK_BOX(carbon_ring_picker), buttom_7C);
    GtkWidget* buttom_8C = gtk_button_new_with_label("8-C");
    gtk_box_append(GTK_BOX(carbon_ring_picker), buttom_8C);
    GtkWidget* buttom_EnvResidues = gtk_button_new_with_label("Env. Residues");
    gtk_box_append(GTK_BOX(carbon_ring_picker), buttom_EnvResidues);
    GtkWidget* buttom_Key = gtk_button_new_with_label("Key");
    gtk_box_append(GTK_BOX(carbon_ring_picker), buttom_Key);

    GtkWidget* canvas_space = gtk_box_new(GTK_ORIENTATION_HORIZONTAL,0);
    gtk_box_append(GTK_BOX(mainbox), canvas_space);
    // Canvas space: chemical element picker
    GtkWidget* chem_element_picker = gtk_box_new(GTK_ORIENTATION_VERTICAL,5);
    gtk_box_append(GTK_BOX(canvas_space), chem_element_picker);

    GtkWidget* C_button = gtk_button_new_with_label("C");
    gtk_box_append(GTK_BOX(chem_element_picker), C_button);
    GtkWidget* N_button = gtk_button_new_with_label("N");
    gtk_box_append(GTK_BOX(chem_element_picker), N_button);
    GtkWidget* O_button = gtk_button_new_with_label("O");
    gtk_box_append(GTK_BOX(chem_element_picker), O_button);
    GtkWidget* S_button = gtk_button_new_with_label("S");
    gtk_box_append(GTK_BOX(chem_element_picker), S_button);
    GtkWidget* P_button = gtk_button_new_with_label("P");
    gtk_box_append(GTK_BOX(chem_element_picker), P_button);
    GtkWidget* H_button = gtk_button_new_with_label("H");
    gtk_box_append(GTK_BOX(chem_element_picker), H_button);
    GtkWidget* F_button = gtk_button_new_with_label("F");
    gtk_box_append(GTK_BOX(chem_element_picker), F_button);
    GtkWidget* Cl_button = gtk_button_new_with_label("Cl");
    gtk_box_append(GTK_BOX(chem_element_picker), Cl_button);
    GtkWidget* Br_button = gtk_button_new_with_label("Br");
    gtk_box_append(GTK_BOX(chem_element_picker), Br_button);
    GtkWidget* I_button = gtk_button_new_with_label("I");
    gtk_box_append(GTK_BOX(chem_element_picker), I_button);
    GtkWidget* X_button = gtk_button_new_with_label("X");
    gtk_box_append(GTK_BOX(chem_element_picker), X_button);
    // Canvas space: Canvas
    //GtkWidget* canvas = gtk_label_new("The widget will land here");
    auto* canvas = coot_ligand_editor_canvas_new();
    gtk_box_append(GTK_BOX(canvas_space), GTK_WIDGET(canvas));
    // Statusbar / the bottom
    GtkWidget* bottom_bar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 10);
    gtk_box_append(GTK_BOX(mainbox), bottom_bar);
    gtk_window_set_child(win, mainbox);
    
}

GMenu *build_menu() {
    GMenu *ret = g_menu_new();
    // g_menu_append(GMenu *menu, const gchar *label, const gchar
    // *detailed_action);

    // File
    GMenu *file = g_menu_new();
    g_menu_append(file,"Test","app.test");
    g_menu_append_section(ret, "File", G_MENU_MODEL(file));
    // Display
    GMenu *display = g_menu_new();
    g_menu_append_section(ret, "Display", G_MENU_MODEL(display));
    // Help
    GMenu *help = g_menu_new();
    g_menu_append_section(ret, "Help", G_MENU_MODEL(help));

    return ret;
}

int main() {
    gtk_init();
    
    GtkApplication* app = gtk_application_new("org.pemsley.NewLigandEditor",G_APPLICATION_DEFAULT_FLAGS);
    GError *error = NULL;
    g_application_register(G_APPLICATION(app), NULL, &error);

    g_signal_connect(app,"activate",G_CALLBACK(+[](GtkApplication* app, gpointer user_data){
        //GtkWindow* win = GTK_WINDOW(user_data);
        gtk_application_set_menubar(app, G_MENU_MODEL(build_menu()));
        GtkWidget* win = gtk_application_window_new(app);
        gtk_window_set_application(GTK_WINDOW(win),app);
        gtk_application_add_window(app,GTK_WINDOW(win));
        gtk_application_window_set_show_menubar(GTK_APPLICATION_WINDOW(win), TRUE);
        build_main_window(GTK_WINDOW(win));
        gtk_widget_show(win);

    }),NULL);


    return g_application_run(G_APPLICATION(app),0,0);
}