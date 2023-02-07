#include "ligand-builder.hpp"
#include "ligand_editor_canvas.hpp"
#include "ligand_editor_canvas/tools.hpp"
#include <gtk/gtk.h>
#include <string>



void build_main_window(GtkWindow* win, CootLigandEditorCanvas* canvas) {
    using coot::ligand_editor_canvas::ActiveTool;
    using coot::ligand_editor_canvas::ElementInsertion;
    using Element = coot::ligand_editor_canvas::ElementInsertion::Element;

    GtkWidget* mainbox = gtk_box_new(GTK_ORIENTATION_VERTICAL,10);

    gtk_window_set_child(win, mainbox);
    gtk_widget_set_margin_start(mainbox,10);
    gtk_widget_set_margin_end(mainbox,10);
    // Top toolbar
    GtkWidget* top_toolbar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL,5);
    gtk_box_append(GTK_BOX(mainbox), top_toolbar);

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
    gtk_widget_set_margin_top(chem_element_picker,10);

    gtk_box_append(GTK_BOX(canvas_space), chem_element_picker);

    GtkWidget* C_button = gtk_button_new_with_label("C");
    g_signal_connect(C_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::C)));
    }), canvas);
    gtk_box_append(GTK_BOX(chem_element_picker), C_button);
    GtkWidget* N_button = gtk_button_new_with_label("N");
    g_signal_connect(N_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::N)));
    }), canvas);
    gtk_box_append(GTK_BOX(chem_element_picker), N_button);
    GtkWidget* O_button = gtk_button_new_with_label("O");
    g_signal_connect(O_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::O)));
    }), canvas);
    gtk_box_append(GTK_BOX(chem_element_picker), O_button);
    GtkWidget* S_button = gtk_button_new_with_label("S");
    g_signal_connect(S_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::S)));
    }), canvas);
    gtk_box_append(GTK_BOX(chem_element_picker), S_button);
    GtkWidget* P_button = gtk_button_new_with_label("P");
    g_signal_connect(P_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::P)));
    }), canvas);
    gtk_box_append(GTK_BOX(chem_element_picker), P_button);
    GtkWidget* H_button = gtk_button_new_with_label("H");
    g_signal_connect(H_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::H)));
    }), canvas);
    gtk_box_append(GTK_BOX(chem_element_picker), H_button);
    GtkWidget* F_button = gtk_button_new_with_label("F");
    g_signal_connect(F_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::F)));
    }), canvas);
    gtk_box_append(GTK_BOX(chem_element_picker), F_button);
    GtkWidget* Cl_button = gtk_button_new_with_label("Cl");
    g_signal_connect(Cl_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::Cl)));
    }), canvas);
    gtk_box_append(GTK_BOX(chem_element_picker), Cl_button);
    GtkWidget* Br_button = gtk_button_new_with_label("Br");
    g_signal_connect(Br_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::Br)));
    }), canvas);
    gtk_box_append(GTK_BOX(chem_element_picker), Br_button);
    GtkWidget* I_button = gtk_button_new_with_label("I");
    g_signal_connect(I_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::I)));
    }), canvas);
    gtk_box_append(GTK_BOX(chem_element_picker), I_button);
    GtkWidget* X_button = gtk_button_new_with_label("X");
    g_signal_connect(X_button, "clicked", G_CALLBACK(+[](GtkButton* _btn, gpointer user_data){
        CootLigandEditorCanvas* canvas = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
        coot_ligand_editor_set_active_tool(canvas, std::make_unique<ActiveTool>(ElementInsertion(Element::X)));
        // todo: make this work
    }), canvas);
    gtk_box_append(GTK_BOX(chem_element_picker), X_button);
    // Canvas space: Canvas
    GtkWidget* canvas_viewport = gtk_scrolled_window_new();
    gtk_widget_set_hexpand(canvas_viewport,TRUE);
    gtk_widget_set_vexpand(canvas_viewport,TRUE);
    gtk_widget_set_margin_top(canvas_viewport,10);
    //gtk_widget_set_margin_end(canvas_viewport,10);
    gtk_widget_set_margin_start(canvas_viewport,10);
    gtk_widget_set_margin_bottom(canvas_viewport,10);
 
    gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(canvas_viewport),GTK_WIDGET(canvas));
    gtk_box_append(GTK_BOX(canvas_space), GTK_WIDGET(canvas_viewport));
    // Statusbar / the bottom
    GtkWidget* bottom_bar = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 10);
    gtk_widget_set_margin_start(bottom_bar, 10);
    gtk_widget_set_margin_end(bottom_bar, 10);
    gtk_box_append(GTK_BOX(mainbox), bottom_bar);

    GtkWidget* smiles_label = gtk_label_new("SMILES:");
    gtk_box_append(GTK_BOX(bottom_bar),smiles_label);
    GtkWidget* smiles_display_label = gtk_label_new("");
    gtk_box_append(GTK_BOX(bottom_bar),smiles_display_label);
    gtk_widget_set_hexpand(smiles_display_label, TRUE);
    GtkWidget* scale_label = gtk_label_new("Scale");
    gtk_box_append(GTK_BOX(bottom_bar),scale_label);
    gtk_widget_set_halign(scale_label,GTK_ALIGN_END);
    GtkAdjustment* adj = gtk_adjustment_new(1, 0.1, 20, 0.1, 1, 2);
    GtkWidget* scale_spin_button = gtk_spin_button_new(adj, 0.1, 1);
    gtk_box_append(GTK_BOX(bottom_bar),scale_spin_button);
    gtk_widget_set_halign(scale_spin_button,GTK_ALIGN_END);

    GtkWidget* show_alerts_checkbutton = gtk_check_button_new_with_label("Show Alerts");
    gtk_widget_set_halign(show_alerts_checkbutton,GTK_ALIGN_END);
    gtk_box_append(GTK_BOX(mainbox), show_alerts_checkbutton);
    gtk_widget_set_margin_end(show_alerts_checkbutton, 10);

    GtkWidget* button_box = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 5);
    gtk_box_append(GTK_BOX(mainbox), button_box);
    gtk_widget_set_halign(button_box,GTK_ALIGN_END);
    gtk_widget_set_margin_end(button_box, 10);
    gtk_widget_set_margin_start(button_box, 10);
    gtk_widget_set_margin_top(button_box, 10);
    gtk_widget_set_margin_bottom(button_box, 10);

    GtkWidget* apply_button = gtk_button_new_with_label("Apply");
    gtk_box_append(GTK_BOX(button_box),apply_button);
    GtkWidget* close_button = gtk_button_new_with_label("Close");
    gtk_box_append(GTK_BOX(button_box),close_button);
    g_signal_connect(close_button, "clicked", G_CALLBACK(+[](GtkWidget* button, gpointer user_data){
        gtk_window_close(GTK_WINDOW(user_data));
    }), win);
    
}

GMenu *build_menu(GtkApplication* app, CootLigandEditorCanvas* canvas, GtkWindow* win) {
    using ExportMode = coot::ligand_editor::LigandBuilderState::ExportMode;

    GMenu *ret = g_menu_new();
    
    // g_menu_append(GMenu *menu, const gchar *label, const gchar
    // *detailed_action);
    auto new_menu_item = [app](const char* label,const char* action_name,GCallback func, gpointer userdata = nullptr){
        std::string detailed_action_name = "app.";
        detailed_action_name += action_name;
        GMenuItem* item = g_menu_item_new(label,detailed_action_name.c_str());
        GSimpleAction* action = g_simple_action_new(action_name,nullptr);
        g_action_map_add_action(G_ACTION_MAP(app), G_ACTION(action));
        g_signal_connect(action, "activate", func, userdata);
        return item;
    };

    // File
    GMenu *file = g_menu_new();
    g_menu_append_submenu(ret, "File", G_MENU_MODEL(file));
    g_menu_append_item(file, new_menu_item("_New", "file_new", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        g_info("TODO: New");
    })));
    g_menu_append_item(file, new_menu_item("_Open", "file_open", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        g_info("TODO: Open");
    })));
    g_menu_append_item(file, new_menu_item("Import from SMILES", "import_from_smiles", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        coot::ligand_editor::global_instance->load_from_smiles();
    })));
    g_menu_append_item(file, new_menu_item("Import Molecule", "import_molecule", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        g_info("TODO: Import Molecule");
    })));
    g_menu_append_item(file, new_menu_item("Fetch Molecule", "fetch_molecule", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        g_info("TODO: Fetch Molecule");
    })));
    g_menu_append_item(file, new_menu_item("Save", "file_save", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        g_info("TODO: Save");
    })));
    g_menu_append_item(file, new_menu_item("Export as PDF", "export_pdf", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        coot::ligand_editor::global_instance->file_export(ExportMode::PDF);
    })));
    g_menu_append_item(file, new_menu_item("Export as PNG", "export_png", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        coot::ligand_editor::global_instance->file_export(ExportMode::PNG);
    })));
    g_menu_append_item(file, new_menu_item("Export as SVG", "export_svg", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        coot::ligand_editor::global_instance->file_export(ExportMode::SVG);
    })));
    g_menu_append_item(file, new_menu_item("_Exit", "exit", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        gtk_window_close(GTK_WINDOW(user_data));
    }),win));
    // Edit
    GMenu *edit = g_menu_new();
    g_menu_append_submenu(ret, "Edit", G_MENU_MODEL(edit));
    g_menu_append_item(edit, new_menu_item("Undo", "undo", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        g_info("TODO: Undo");
    })));
    g_menu_append_item(edit, new_menu_item("Redo", "redo", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        g_info("TODO: Redo");
    })));
    // Display
    GMenu *display = g_menu_new();
    g_menu_append_submenu(ret, "Display", G_MENU_MODEL(display));
    g_menu_append_item(display, new_menu_item("Standard", "display_standard", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        g_info("TODO: Standard");
    })));
    g_menu_append_item(display, new_menu_item("Atom Indices", "display_atom_indices", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        g_info("TODO: Atom Indices");
    })));
    g_menu_append_item(display, new_menu_item("Atom Names", "display_atom_names", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        g_info("TODO: Atom Names");
    })));
    // Help
    GMenu *help = g_menu_new();
    g_menu_append_submenu(ret, "Help", G_MENU_MODEL(help));
    g_menu_append_item(help, new_menu_item("About", "about", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        g_info("TODO: About");
    })));

    return ret;
}

int main() {
    gtk_init();
    
    GtkApplication* app = gtk_application_new("org.pemsley.NewLigandEditor",G_APPLICATION_DEFAULT_FLAGS);
    GError *error = NULL;
    g_application_register(G_APPLICATION(app), NULL, &error);

    g_signal_connect(app,"activate",G_CALLBACK(+[](GtkApplication* app, gpointer user_data){
        //GtkWindow* win = GTK_WINDOW(user_data);
        GtkWidget* win = gtk_application_window_new(app);
        gtk_window_set_title(GTK_WINDOW(win),"New Ligand Editor");
        gtk_application_window_set_show_menubar(GTK_APPLICATION_WINDOW(win), TRUE);
        gtk_window_set_application(GTK_WINDOW(win),app);
        auto* canvas = coot_ligand_editor_canvas_new();
        coot::ligand_editor::initialize_global_instance(canvas,GTK_WINDOW(win));
        gtk_application_set_menubar(app, G_MENU_MODEL(build_menu(app,canvas,GTK_WINDOW(win))));
        gtk_application_add_window(app,GTK_WINDOW(win));
        build_main_window(GTK_WINDOW(win),canvas);
        gtk_widget_show(win);

    }),NULL);


    auto ret = g_application_run(G_APPLICATION(app),0,0);
    g_info("Exiting...");
    delete coot::ligand_editor::global_instance;
    return ret;
}