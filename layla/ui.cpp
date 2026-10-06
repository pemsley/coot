/* layla/ui.cpp
 * 
 * Copyright 2023 by Global Phasing Ltd.
 * Author: Jakub Smulski
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include "ui.hpp"
#include "state.hpp"
#include "ligand_editor_canvas.hpp"
#include "ligand_editor_canvas/core.hpp"
#include "ligand_editor_canvas/model.hpp"
#include "ligand_editor_canvas/tools.hpp"
#include "utils/coot-utils.hh"
#include "qed.hpp"
#include <functional>
#include <cmath>

namespace {

   // Per-cell state for a QED desirability-curve mini-plot: which property it
   // shows, and the current molecule's raw value for that property (if any).
   struct desirability_curve_cell_t {
      coot::layla::RDKit::QED::QEDPropName prop;
      bool has_value = false;
      double value = 0.0;
   };

   // "Nice numbers" for axis ticks (Heckbert). round_it picks the nearest nice
   // value; otherwise the smallest nice value >= x.
   double nice_num(double x, bool round_it) {
      if (x <= 0.0) return 0.0;
      double expv = std::floor(std::log10(x));
      double f = x / std::pow(10.0, expv);
      double nf;
      if (round_it) {
         if (f < 1.5) nf = 1.0; else if (f < 3.0) nf = 2.0; else if (f < 7.0) nf = 5.0; else nf = 10.0;
      } else {
         if (f <= 1.0) nf = 1.0; else if (f <= 2.0) nf = 2.0; else if (f <= 5.0) nf = 5.0; else nf = 10.0;
      }
      return nf * std::pow(10.0, expv);
   }

   // Generate ~n_intervals "nice" tick positions spanning [lo, hi].
   std::vector<double> nice_ticks(double lo, double hi, int n_intervals) {
      std::vector<double> out;
      if (hi <= lo || n_intervals < 1) return out;
      double d = nice_num(nice_num(hi - lo, false) / n_intervals, true);
      if (d <= 0.0) return out;
      double graph_lo = std::floor(lo / d) * d;
      double graph_hi = std::ceil(hi / d) * d;
      for (double v = graph_lo; v <= graph_hi + 0.5 * d; v += d)
         out.push_back(std::fabs(v) < 1e-9 ? 0.0 : v); // avoid "-0"
      return out;
   }

   // GtkDrawingArea draw function: draws the property's ADS desirability curve
   // d(x) over its Fig.1 range, filled underneath, with a marker at the current
   // molecule's value. (Parameters live in layla/qed.cpp, from Bickerton et al.
   // 2012, Suppl. Table 2.)
   void draw_desirability_curve(GtkDrawingArea *area, cairo_t *cr,
                                int width, int height, gpointer user_data) {

      using QED = coot::layla::RDKit::QED;
      auto *cell = static_cast<desirability_curve_cell_t *>(user_data);
      if (! cell) return;

      // Adapt colours to the theme: read the widget's foreground colour and use
      // its luminance to decide light vs dark background (light text => dark bg).
      // The paper's dark blue reads poorly on a dark background, so brighten it.
      GdkRGBA fg = {0.0f, 0.0f, 0.0f, 1.0f}; // default: assume a light background
#if GTK_MINOR_VERSION >= 10
      gtk_widget_get_color(GTK_WIDGET(area), &fg);
#endif
      double fg_lum = 0.2126 * fg.red + 0.7152 * fg.green + 0.0722 * fg.blue;
      bool dark_background = fg_lum > 0.5;
      double line_r, line_g, line_b, fill_a;
      if (dark_background) {
         line_r = 0.45; line_g = 0.70; line_b = 1.00; fill_a = 0.22; // brighter blue
      } else {
         line_r = 0.169; line_g = 0.247; line_b = 0.667; fill_a = 0.25; // paper #2b3faa
      }

      const QED::ADSparameter &p = QED::get_ads_parameter(cell->prop);
      QED::plot_range_t range = QED::get_plot_range(cell->prop);
      double span = range.x_max - range.x_min;
      if (span <= 0.0) return;

      const double ml = 6.0, mr = 6.0, mt = 6.0, mb = 16.0; // mb leaves room for tick labels
      double pw = width  - ml - mr;
      double ph = height - mt - mb;
      if (pw <= 1.0 || ph <= 1.0) return;

      auto x_to_px = [&](double x){ return ml + (x - range.x_min) / span * pw; };
      auto d_to_py = [&](double d){ return mt + (1.0 - d) * ph; };
      auto clamp01 = [](double d){ return d < 0.0 ? 0.0 : (d > 1.0 ? 1.0 : d); };

      const int N = 96;

      // filled area under the curve
      cairo_new_path(cr);
      cairo_move_to(cr, x_to_px(range.x_min), d_to_py(0.0));
      for (int i=0; i<=N; i++) {
         double x = range.x_min + span * i / static_cast<double>(N);
         cairo_line_to(cr, x_to_px(x), d_to_py(clamp01(QED::ads(x, p))));
      }
      cairo_line_to(cr, x_to_px(range.x_max), d_to_py(0.0));
      cairo_close_path(cr);
      cairo_set_source_rgba(cr, line_r, line_g, line_b, fill_a); // faint fill under curve
      cairo_fill(cr);

      // the curve line
      cairo_new_path(cr);
      for (int i=0; i<=N; i++) {
         double x = range.x_min + span * i / static_cast<double>(N);
         double px = x_to_px(x), py = d_to_py(clamp01(QED::ads(x, p)));
         if (i == 0) cairo_move_to(cr, px, py); else cairo_line_to(cr, px, py);
      }
      cairo_set_line_width(cr, 1.5);
      cairo_set_source_rgb(cr, line_r, line_g, line_b); // theme-adapted blue
      cairo_stroke(cr);

      // marker at the current molecule's value
      if (cell->has_value) {
         double xv = cell->value;
         if (xv < range.x_min) xv = range.x_min;
         if (xv > range.x_max) xv = range.x_max;
         double d = clamp01(QED::ads(cell->value, p));
         double px = x_to_px(xv);
         cairo_set_line_width(cr, 1.0);
         cairo_set_source_rgba(cr, fg.red, fg.green, fg.blue, 0.55); // dropline (theme fg)
         cairo_move_to(cr, px, d_to_py(0.0));
         cairo_line_to(cr, px, d_to_py(d));
         cairo_stroke(cr);
         cairo_arc(cr, px, d_to_py(d), 3.0, 0.0, 2.0 * M_PI);
         cairo_set_source_rgb(cr, 0.85, 0.15, 0.15); // red dot
         cairo_fill(cr);
      }

      // x-axis baseline and "nice" tick labels
      double y0 = d_to_py(0.0);
      cairo_set_line_width(cr, 1.0);
      cairo_set_source_rgba(cr, fg.red, fg.green, fg.blue, 0.35);
      cairo_move_to(cr, x_to_px(range.x_min), y0);
      cairo_line_to(cr, x_to_px(range.x_max), y0);
      cairo_stroke(cr);

      cairo_set_font_size(cr, 8.0);
      cairo_set_source_rgba(cr, fg.red, fg.green, fg.blue, 0.75);
      for (double t : nice_ticks(range.x_min, range.x_max, 3)) {
         if (t < range.x_min - 1e-9 || t > range.x_max + 1e-9) continue;
         double px = x_to_px(t);
         cairo_move_to(cr, px, y0);
         cairo_line_to(cr, px, y0 + 3.0); // short tick mark
         cairo_stroke(cr);
         char lab[32];
         g_snprintf(lab, sizeof lab, "%g", t);
         cairo_text_extents_t ext;
         cairo_text_extents(cr, lab, &ext);
         double tx = px - ext.width / 2.0 - ext.x_bearing;   // centre under the tick
         if (tx < 1.0) tx = 1.0;                              // keep end labels on-screen
         if (tx + ext.width > width - 1.0) tx = width - 1.0 - ext.width;
         cairo_move_to(cr, tx, height - 3.0);
         cairo_show_text(cr, lab);
      }
   }
}

void setup_actions(coot::layla::LaylaState* state, GtkApplicationWindow* win, GtkBuilder* builder) {
    using namespace coot::layla;

    auto new_action = [win](const char* action_name, GCallback func, gpointer userdata = nullptr){
        std::string detailed_action_name = "win.";
        detailed_action_name += action_name;
        GSimpleAction* action = g_simple_action_new(action_name,nullptr);
        g_action_map_add_action(G_ACTION_MAP(win), G_ACTION(action));
        g_signal_connect(action, "activate", func, userdata);
        //return std::make_pair(detailed_action_name,action);
    };

    auto new_stateful_action = [win](const char* action_name,const GVariantType *state_type, GVariant* default_state, GCallback func, gpointer userdata = nullptr){
        std::string detailed_action_name = "win.";
        detailed_action_name += action_name;
        GSimpleAction* action = g_simple_action_new_stateful(action_name, state_type, default_state);
        g_action_map_add_action(G_ACTION_MAP(win), G_ACTION(action));
        g_signal_connect(action, "activate", func, userdata);
        //return std::make_pair(detailed_action_name,action);
    };

    using ExportMode = coot::layla::LaylaState::ExportMode;

    // File
    new_action("file_new", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LaylaState*)user_data)->file_new();
    }),state);
    new_action("file_open", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LaylaState*)user_data)->file_open();
    }),state);
    new_action("import_from_smiles", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LaylaState*)user_data)->load_from_smiles();
    }),state);
    new_action("import_molecule", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LaylaState*)user_data)->file_import_molecule();
    }),state);
    new_action("fetch_molecule", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LaylaState*)user_data)->file_fetch_molecule();
    }),state);
    new_action("file_save", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LaylaState*)user_data)->file_save();
    }),state);
    new_action("file_save_as", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LaylaState*)user_data)->file_save_as();
    }),state);
    new_action("export_pdf", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LaylaState*)user_data)->file_export(ExportMode::PDF);
    }),state);
    new_action("export_png", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LaylaState*)user_data)->file_export(ExportMode::PNG);
    }),state);
    new_action("export_svg", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LaylaState*)user_data)->file_export(ExportMode::SVG);
    }),state);
    new_action("file_exit", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LaylaState*)user_data)->file_exit();
    }),state);

    // Edit;
    new_action("undo", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LaylaState*)user_data)->edit_undo();
    }),state);
    new_action("redo", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        ((LaylaState*)user_data)->edit_redo();
    }),state);
    // Display

    using coot::ligand_editor_canvas::DisplayMode;
    GVariant* display_mode_action_defstate = g_variant_new("s",coot::ligand_editor_canvas::display_mode_to_string(DisplayMode::Standard));
    new_stateful_action(
        "switch_display_mode",
        G_VARIANT_TYPE_STRING,
        display_mode_action_defstate, 
        G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
            const gchar* mode_name = g_variant_get_string(parameter,nullptr);
            auto mode = coot::ligand_editor_canvas::display_mode_from_string(mode_name);
            if(mode.has_value()) {
                ((LaylaState*)user_data)->switch_display_mode(mode.value());
                g_simple_action_set_state(self, parameter);
            } else {
                g_error("Could not parse display mode from string!: '%s'",mode_name);
            }
        }
    ),state);

    // Help

    new_action("show_about_dialog", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        auto* about_dialog = GTK_WINDOW(user_data);
        gtk_window_present(GTK_WINDOW(about_dialog));
    }),gtk_builder_get_object(builder, "layla_about_dialog"));

    new_action("show_shortcuts_window", G_CALLBACK(+[](GSimpleAction* self, GVariant* parameter, gpointer user_data){
        auto* window = GTK_WINDOW(user_data);
        gtk_window_present(GTK_WINDOW(window));
    }),gtk_builder_get_object(builder, "layla_shortcuts_window"));

}

GtkApplicationWindow* coot::layla::setup_main_window(GtkApplication* app, GtkBuilder* builder) {

    GtkApplicationWindow* win = (GtkApplicationWindow*) gtk_builder_get_object(builder, "layla_window");
    gtk_window_set_application(GTK_WINDOW(win),app);
    GtkWidget* status_label = (GtkWidget*) gtk_builder_get_object(builder, "layla_status_label");
    GtkScrolledWindow* viewport = (GtkScrolledWindow*) gtk_builder_get_object(builder, "layla_canvas_viewport");
    auto* canvas = coot_ligand_editor_canvas_new();

    g_signal_connect(canvas, "status-updated", G_CALLBACK(+[](CootLigandEditorCanvas* canvas, const gchar* status_text, gpointer user_data){
        gtk_label_set_text(GTK_LABEL(user_data), status_text);
    }), status_label);

    // "Show Alerts": while enabled, hovering the canvas shows which structural
    // alerts were matched (the on-canvas circles show where). Read live so the
    // list stays current as the molecule is edited.
    gtk_widget_set_has_tooltip(GTK_WIDGET(canvas), TRUE);
    g_signal_connect(canvas, "query-tooltip", G_CALLBACK(+[](GtkWidget* widget, gint x, gint y, gboolean keyboard_mode, GtkTooltip* tooltip, gpointer user_data) -> gboolean {
        CootLigandEditorCanvas* c = (CootLigandEditorCanvas*) user_data;
        if(!coot_ligand_editor_canvas_get_show_alerts(c)) {
            return FALSE;
        }
        std::string summary = coot_ligand_editor_canvas_get_alert_summary(c);
        if(summary.empty()) {
            return FALSE;
        }
        std::string text = "Structural alerts: " + summary;
        gtk_tooltip_set_text(tooltip, text.c_str());
        return TRUE;
    }), canvas);

    GtkSpinButton* scale_spinbutton = (GtkSpinButton*) gtk_builder_get_object(builder, "layla_scale_spinbutton");
    g_signal_connect(canvas, "scale-changed", G_CALLBACK(+[](CootLigandEditorCanvas* canvas, float new_scale, gpointer user_data){
        GtkSpinButton* spinbutton = GTK_SPIN_BUTTON(user_data);
        gtk_spin_button_set_value(spinbutton, new_scale);
    }), scale_spinbutton);


    GtkGrid* smiles_display_grid = (GtkGrid*) gtk_builder_get_object(builder, "layla_smiles_display_grid");
    g_signal_connect(canvas, "smiles-changed", G_CALLBACK(+[](CootLigandEditorCanvas* self, gpointer user_data) {
        GtkGrid* display_grid = GTK_GRID(user_data);
        // Don't clear the widget all the time. It will create horrible mess
        // for(auto* i = gtk_widget_get_first_child(GTK_WIDGET(display_grid)); i != nullptr; i = gtk_widget_get_next_sibling(GTK_WIDGET(i))) {
        //     gtk_grid_remove(display_grid, i);
        // }
        struct WidgetsT {
            GtkWidget* smiles_label;
            GtkWidget* inchi_label;
        };
        auto get_widgets_for_mol_id = [display_grid](unsigned int id) -> std::optional<WidgetsT> {
            WidgetsT ret;
            ret.smiles_label = nullptr;
            ret.inchi_label = nullptr;
            for (auto* i = gtk_widget_get_first_child(GTK_WIDGET(display_grid)); i != nullptr; i = gtk_widget_get_next_sibling(GTK_WIDGET(i))) {
                bool is_editable = GTK_IS_EDITABLE_LABEL(i);
                gpointer inchi_label_gptr = g_object_get_data(G_OBJECT(i), "inchi_label");
                if (!is_editable && !inchi_label_gptr) {
                    // Skipping irrelevant labels / widgets
                    continue;
                }
                // g_info("Inspecting a widget for \"mol_id\"");
                gpointer mol_id_gptr = g_object_get_data(G_OBJECT(i), "mol_id");
                if (mol_id_gptr) {
                    // Differeniate between nullptr and zero
                    if (GPOINTER_TO_UINT(mol_id_gptr) - 1 == id) {
                        if (GTK_IS_EDITABLE_LABEL(i)) {
                            ret.smiles_label = i;
                        } else {
                            ret.inchi_label = i;
                        }
                    }
                }
            }
            // g_info("smiles_label %p inchi_label %p", ret.smiles_label, ret.inchi_label);
            if (!ret.inchi_label || !ret.smiles_label) {
                return std::nullopt;
            } else {
                return ret;
            }
        };

        auto smiles_map = coot_ligand_editor_canvas_get_smiles(self);

        auto do_inchi_lookup = [&](unsigned int mol_idx) -> std::string {
            std::string inchi_key = coot_ligand_editor_canvas_get_inchi_key_for_molecule(self, mol_idx);
            auto lookup_result = coot::layla::global_instance->lookup_inchi_key(inchi_key);
            if(lookup_result.has_value()) {
                auto lookup_result_value = lookup_result.value();
                return lookup_result_value.chemical_name + " " +  lookup_result_value.monomer_id;
            } else {
                return "";
            }
        };

        for (const auto& [mol_idx, smiles_code] : smiles_map) {
            auto widgets = get_widgets_for_mol_id(mol_idx);
            if (widgets.has_value()) {
                const auto m_widgets = widgets.value();
                if (!gtk_editable_label_get_editing(GTK_EDITABLE_LABEL(m_widgets.smiles_label))) {
                    g_info("Updating SMILES text for mol %u", mol_idx);
                    gtk_editable_set_text(GTK_EDITABLE(m_widgets.smiles_label), smiles_code.c_str());
                } else {
                    g_info("Not updating SMILES for mol %u, with SMILES text currently being edited.", mol_idx);
                }
                auto inchi_str = do_inchi_lookup(mol_idx);
                gtk_label_set_text(GTK_LABEL(m_widgets.inchi_label), inchi_str.c_str());
            } else {
                g_info("Creating SMILES row for mol %u ", mol_idx);
                // Create widgets
                auto l_str = std::to_string(mol_idx) + ":";
                GtkLabel* label = (GtkLabel*)gtk_label_new(l_str.c_str());
                // Differeniate between nullptr and zero
                g_object_set_data(G_OBJECT(label), "mol_id", GUINT_TO_POINTER(mol_idx + 1));
                gtk_grid_attach(display_grid, GTK_WIDGET(label), 0, mol_idx, 1, 1);

                GtkEditableLabel* smiles_label = (GtkEditableLabel*)gtk_editable_label_new(smiles_code.c_str());
                // Differeniate between nullptr and zero
                g_object_set_data(G_OBJECT(smiles_label), "mol_id", GUINT_TO_POINTER(mol_idx + 1));
                g_signal_connect(smiles_label, "changed", G_CALLBACK(+[](GtkEditable* smiles_label, gpointer user_data) {
                    if (!gtk_editable_label_get_editing(GTK_EDITABLE_LABEL(smiles_label))) {
                        // We don't want to update the mol from smiles if the text changed
                        // as a result of the molecule being modified on screen by the user
                        return;
                    }
                    CootLigandEditorCanvas* self = COOT_COOT_LIGAND_EDITOR_CANVAS(user_data);
                    std::string smiles_text = gtk_editable_get_text(smiles_label);
                    unsigned int mol_id = GPOINTER_TO_UINT(g_object_get_data(G_OBJECT(smiles_label), "mol_id")) - 1;
                    coot_ligand_editor_canvas_update_molecule_from_smiles(self, mol_id, smiles_text.c_str());
                }),
                self);

                GtkCssProvider* provider = gtk_css_provider_new();
                std::string smiles_label_font_string = ".smiles_label { font-family: monospace; }";
    #if (GTK_MAJOR_VERSION == 4) && (GTK_MINOR_VERSION > 11)  // available since 4.12
                gtk_css_provider_load_from_string(provider, smiles_label_font_string.c_str());
    #else
    gtk_css_provider_load_from_data(provider, smiles_label_font_string.c_str(), smiles_label_font_string.size());
    #endif
                GtkStyleContext* context = gtk_widget_get_style_context(GTK_WIDGET(smiles_label));
                gtk_style_context_add_provider(context, GTK_STYLE_PROVIDER(provider), GTK_STYLE_PROVIDER_PRIORITY_APPLICATION);
                gtk_style_context_add_class(context, "smiles_label");

                gtk_grid_attach(display_grid, GTK_WIDGET(smiles_label), 1, mol_idx, 1, 1);

                auto inchi_str = do_inchi_lookup(mol_idx);
                GtkLabel* inchi_data_label = (GtkLabel*) gtk_label_new(inchi_str.c_str());
                g_object_set_data(G_OBJECT(inchi_data_label), "mol_id", GUINT_TO_POINTER(mol_idx + 1));
                g_object_set_data(G_OBJECT(inchi_data_label), "inchi_label", GUINT_TO_POINTER(1));
                gtk_grid_attach(display_grid, GTK_WIDGET(inchi_data_label), 2, mol_idx, 1, 1);
            }
        }
    }),
    smiles_display_grid);

    g_signal_connect(canvas, "molecule-deleted", G_CALLBACK(+[](CootLigandEditorCanvas* self, unsigned int deleted_mol_idx, gpointer user_data){
        GtkGrid* display_grid = GTK_GRID(user_data);
        // Prevents iterator invalidation
        std::vector<GtkWidget*> to_be_removed(3);
        for(auto* i = gtk_widget_get_first_child(GTK_WIDGET(display_grid)); i != nullptr; i = gtk_widget_get_next_sibling(GTK_WIDGET(i))) {
            // if(g_object_get_data(G_OBJECT(i),"is_id_label")) {
            //     continue;
            // }
            gpointer mol_id_gptr = g_object_get_data(G_OBJECT(i), "mol_id");
            if(mol_id_gptr) {
                // Differeniate between nullptr and zero
                if(GPOINTER_TO_UINT(mol_id_gptr) - 1 == deleted_mol_idx) {
                    to_be_removed.push_back(i);
                }
            }
        }
        for(const auto& i: to_be_removed) {
            gtk_grid_remove(display_grid, i);
        }
    }), smiles_display_grid);

    GtkNotebook* qed_notebook = (GtkNotebook*) gtk_builder_get_object(builder, "layla_qed_notebook");

    auto qed_info_updated_handler = [] (CootLigandEditorCanvas* self, 
                                        unsigned int molecule_id,
                                        const ligand_editor_canvas::CanvasMolecule::QEDInfo *qed_info,
                                        gpointer user_data) {

        GtkNotebook* qed_notebook = GTK_NOTEBOOK(user_data);
        auto find_or_create_tab_for_mol_id = [qed_notebook](unsigned int molecule_id){
            auto no_pages = gtk_notebook_get_n_pages(qed_notebook);
            auto mol_id_as_str = std::to_string(molecule_id);
            for(int i = 0; i != no_pages; i++) {
                GtkWidget* tab = gtk_notebook_get_nth_page(qed_notebook, i);
                const gchar* label = gtk_notebook_get_tab_label_text(qed_notebook, tab);
                if(g_strcmp0(label, mol_id_as_str.c_str()) == 0) {
                    return tab;
                }
            }
            // No tab found. We have to create a new one.
            GtkWidget* n_label = gtk_label_new(mol_id_as_str.c_str());
            GtkWidget* qed_grid = gtk_grid_new();
            /// Setup contents
            gtk_grid_set_column_spacing(GTK_GRID(qed_grid), 15);
            gtk_grid_set_row_spacing(GTK_GRID(qed_grid), 5);
            gtk_widget_set_margin_top(qed_grid, 6);
            gtk_widget_set_margin_bottom(qed_grid, 6);


            auto build_progressbar_info_box = [] (const std::string &label /* range or anything else ??*/) {
                GtkWidget* ret = gtk_box_new(GTK_ORIENTATION_VERTICAL, 5);
                GtkWidget* gtk_label = gtk_label_new(label.c_str());
                gtk_box_append(GTK_BOX(ret), gtk_label);
                GtkWidget* progress_bar = gtk_progress_bar_new();
                gtk_progress_bar_set_show_text(GTK_PROGRESS_BAR(progress_bar), TRUE);
                gtk_box_append(GTK_BOX(ret), progress_bar);
                return ret;
            };

            // A property cell: caption label, the desirability curve (drawn with
            // Cairo, marked at this molecule's value) and a value caption below.
            auto build_curve_cell = [] (const std::string &label,
                                        coot::layla::RDKit::QED::QEDPropName prop) {
                GtkWidget* box = gtk_box_new(GTK_ORIENTATION_VERTICAL, 3);
                gtk_box_append(GTK_BOX(box), gtk_label_new(label.c_str()));
                GtkWidget* area = gtk_drawing_area_new();
                gtk_widget_set_size_request(area, 150, 96);
                gtk_widget_set_hexpand(area, TRUE);
                auto* cell = new desirability_curve_cell_t{prop, false, 0.0};
                g_object_set_data_full(G_OBJECT(area), "curve-cell", cell,
                    +[](gpointer d){ delete static_cast<desirability_curve_cell_t*>(d); });
                gtk_drawing_area_set_draw_func(GTK_DRAWING_AREA(area),
                                               draw_desirability_curve, cell, nullptr);
                gtk_box_append(GTK_BOX(box), area);
                gtk_box_append(GTK_BOX(box), gtk_label_new("")); // value caption
                return box;
            };

            using QEDPropName = coot::layla::RDKit::QED::QEDPropName;
            gtk_grid_attach(GTK_GRID(qed_grid), build_progressbar_info_box("QED"),             0, 0, 1, 1);
            gtk_grid_attach(GTK_GRID(qed_grid), build_curve_cell("MW",        QEDPropName::MW),     0, 1, 1, 1);
            gtk_grid_attach(GTK_GRID(qed_grid), build_curve_cell("PSA",       QEDPropName::PSA),    1, 1, 1, 1);
            gtk_grid_attach(GTK_GRID(qed_grid), build_curve_cell("cLogP",     QEDPropName::ALOGP),  2, 1, 1, 1);
            gtk_grid_attach(GTK_GRID(qed_grid), build_curve_cell("#HBA",      QEDPropName::HBA),    3, 1, 1, 1);
            gtk_grid_attach(GTK_GRID(qed_grid), build_curve_cell("#HBD",      QEDPropName::HBD),    0, 2, 1, 1);
            gtk_grid_attach(GTK_GRID(qed_grid), build_curve_cell("#RotBonds", QEDPropName::ROTB),   1, 2, 1, 1);
            gtk_grid_attach(GTK_GRID(qed_grid), build_curve_cell("#Arom",     QEDPropName::AROM),   2, 2, 1, 1);
            gtk_grid_attach(GTK_GRID(qed_grid), build_curve_cell("#Alerts",   QEDPropName::ALERTS), 3, 2, 1, 1);

            gtk_notebook_append_page(qed_notebook, qed_grid, n_label);
            return qed_grid;
        };

        enum class num_rep_t {FLOAT, INT};

        GtkWidget* tab = find_or_create_tab_for_mol_id(molecule_id);

        auto update_progressbar_info_box = [] (GtkWidget *info_box, num_rep_t t, double value, double progress_bar_value) {
            GtkWidget* label = gtk_widget_get_first_child(info_box);
            GtkWidget* progress_bar = gtk_widget_get_next_sibling(label);
            auto value_as_str = std::to_string(value);
            if (t == num_rep_t::INT) {
               int i = static_cast<int>(value);
               value_as_str = std::to_string(i);
            }
            gtk_progress_bar_set_text(GTK_PROGRESS_BAR(progress_bar), value_as_str.c_str());
            gtk_progress_bar_set_fraction(GTK_PROGRESS_BAR(progress_bar), progress_bar_value);
        };

        // Update a curve cell: set the marker value, redraw, and show "value d=0.xx".
        auto update_curve_cell = [] (GtkWidget *box, num_rep_t t, double value) {
            GtkWidget* caption = gtk_widget_get_first_child(box);
            GtkWidget* area    = gtk_widget_get_next_sibling(caption);
            GtkWidget* vlabel  = gtk_widget_get_next_sibling(area);
            auto* cell = static_cast<desirability_curve_cell_t*>(
                             g_object_get_data(G_OBJECT(area), "curve-cell"));
            if (! cell) return;
            cell->value = value;
            cell->has_value = true;
            double d = coot::layla::RDKit::QED::ads(
                          value, coot::layla::RDKit::QED::get_ads_parameter(cell->prop));
            if (d < 0.0) d = 0.0;
            if (d > 1.0) d = 1.0;
            char buf[64];
            if (t == num_rep_t::INT)
               g_snprintf(buf, sizeof buf, "%d   d=%.2f", static_cast<int>(value), d);
            else
               g_snprintf(buf, sizeof buf, "%.1f   d=%.2f", value, d);
            gtk_label_set_text(GTK_LABEL(vlabel), buf);
            gtk_widget_queue_draw(area);
        };

        // these are (carefully) accessed by grid location, not name:
        update_progressbar_info_box(gtk_grid_get_child_at(GTK_GRID(tab), 0, 0), num_rep_t::FLOAT, qed_info->qed_score, qed_info->qed_score);
        update_curve_cell(gtk_grid_get_child_at(GTK_GRID(tab), 0, 1), num_rep_t::FLOAT, qed_info->molecular_weight);
        update_curve_cell(gtk_grid_get_child_at(GTK_GRID(tab), 1, 1), num_rep_t::FLOAT, qed_info->molecular_polar_surface_area);
        update_curve_cell(gtk_grid_get_child_at(GTK_GRID(tab), 2, 1), num_rep_t::FLOAT, qed_info->alogp);
        update_curve_cell(gtk_grid_get_child_at(GTK_GRID(tab), 3, 1), num_rep_t::INT,   qed_info->number_of_hydrogen_bond_acceptors);
        update_curve_cell(gtk_grid_get_child_at(GTK_GRID(tab), 0, 2), num_rep_t::INT,   qed_info->number_of_hydrogen_bond_donors);
        update_curve_cell(gtk_grid_get_child_at(GTK_GRID(tab), 1, 2), num_rep_t::INT,   qed_info->number_of_rotatable_bonds);
        update_curve_cell(gtk_grid_get_child_at(GTK_GRID(tab), 2, 2), num_rep_t::INT,   qed_info->number_of_aromatic_rings);
        update_curve_cell(gtk_grid_get_child_at(GTK_GRID(tab), 3, 2), num_rep_t::INT,   qed_info->number_of_alerts);

    };
    g_signal_connect(canvas, "qed-info-updated", G_CALLBACK(+qed_info_updated_handler), qed_notebook);

    g_signal_connect(canvas, "molecule-deleted", G_CALLBACK(+[](CootLigandEditorCanvas* self, unsigned int deleted_mol_idx, gpointer user_data) {
        GtkNotebook* qed_notebook = GTK_NOTEBOOK(user_data);
        auto no_pages = gtk_notebook_get_n_pages(qed_notebook);
        auto mol_id_as_str = std::to_string(deleted_mol_idx);
        for(int i = 0; i != no_pages; i++) {
            GtkWidget* tab = gtk_notebook_get_nth_page(qed_notebook, i);
            const gchar* label = gtk_notebook_get_tab_label_text(qed_notebook, tab);
            if(g_strcmp0(label, mol_id_as_str.c_str()) == 0) {
                gtk_notebook_remove_page(qed_notebook, i);
            }
        }
    }), qed_notebook);

    gtk_scrolled_window_set_child(viewport, GTK_WIDGET(canvas));
    coot::layla::initialize_global_instance(canvas,GTK_WINDOW(win),GTK_LABEL(status_label));
    setup_actions(coot::layla::global_instance, win, builder);
    return win;
}

GtkBuilder* coot::layla::load_gtk_builder() {

        g_info("Loading Layla's UI...");

        std::string dir = coot::package_data_dir();
        // all ui files should live here:
        std::string dir_ui = coot::util::append_dir_dir(dir, "ui");
        std::string ui_file_name = "layla.ui";
        std::string ui_file_full = coot::util::append_dir_file(dir_ui, ui_file_name);
        // allow local override
        if(coot::file_exists(ui_file_name)) {
            ui_file_full = ui_file_name;
        }
        GError* error = NULL;
        GtkBuilder* builder = gtk_builder_new();
        gboolean status = gtk_builder_add_from_file(builder, ui_file_full.c_str(), &error);
        if (status == FALSE) {
            g_error("Failed to read or parse %s: %s", ui_file_full.c_str(), error->message);
        }

        return builder;
}
