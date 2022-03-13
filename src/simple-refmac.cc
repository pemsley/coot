
#ifdef USE_PYTHON
#include <Python.h>  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "globjects.h" //includes gtk/gtk.h
#include "interface.h"
#include "graphics-info.h"
#include "c-interface.h" // for is_valid_model_molecule()
#include "cc-interface.hh"
#include "c-interface-gtk-widgets.h"
#include "coot-fileselections.h"

#include "widget-from-builder.hh"

// should this be part of graphics_info, with a wrapper?

void
wrapped_create_simple_refmac_dialog() {

   // GtkWidget *w = create_simple_refmac_dialog();
   GtkWidget *dialog = widget_from_builder("simple_refmac_dialog");

   std::cout << "wrapped_create_simple_refmac_dialog() found dialog " << dialog << std::endl;
   graphics_info_t g;
   int imol_active = -1;
   GCallback callback_func = 0;
   // GtkWidget *combobox_coords = lookup_widget(w, "simple_refmac_coordinates_combobox");
   // GtkWidget *combobox_file   = lookup_widget(w, "simple_refmac_mtz_file_combobox");
   GtkWidget *combobox_coords = widget_from_builder("simple_refmac_coordinates_combobox");
   GtkWidget *combobox_file   = widget_from_builder("simple_refmac_mtz_file_combobox");

   g.fill_combobox_with_coordinates_options(combobox_coords, callback_func, imol_active);

   if (!g.mtz_file_for_refmac.empty()) {
      gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(combobox_file),
                                     g.mtz_file_for_refmac.c_str());
      gtk_combo_box_set_active(GTK_COMBO_BOX(combobox_file), 0);
   }
   gtk_widget_show(dialog);
}

#include "cc-interface-scripting.hh"

void
simple_refmac_run_refmac(GtkWidget *dialog) {

   // GtkWidget *combobox_coords = lookup_widget(dialog, "simple_refmac_coordinates_combobox");
   // GtkWidget *combobox_file   = lookup_widget(dialog, "simple_refmac_mtz_file_combobox");
   GtkWidget *combobox_coords = widget_from_builder("simple_refmac_coordinates_combobox");
   GtkWidget *combobox_file   = widget_from_builder("simple_refmac_mtz_file_combobox");

   graphics_info_t g;
   int imol_coords = g.combobox_get_imol(GTK_COMBO_BOX(combobox_coords));
   std::string mtz_in_filename = g.get_active_label_in_comboboxtext(GTK_COMBO_BOX_TEXT(combobox_file));

   if (! mtz_in_filename.empty())
      g.mtz_file_for_refmac = mtz_in_filename;

   if (is_valid_model_molecule(imol_coords)) {
      short int make_molecules_flag = 1; // not a sub-thread, (so do things
                                         // the normal/old way).
      std::string refmac_dir = coot::get_directory("coot-refmac");
      std::string pdb_in_filename  = coot::util::append_dir_file(refmac_dir, g.molecules[imol_coords].Refmac_in_name());
      std::string pdb_out_filename = coot::util::append_dir_file(refmac_dir, g.molecules[imol_coords].Refmac_out_name());
      std::string mtz_out_filename = coot::util::append_dir_file(refmac_dir, g.molecules[imol_coords].Refmac_mtz_out_name());
      std::string cif_lib_filename;
      std::string fobs_col;
      std::string sigfobs_col;
      std::string r_free_col;
      int sensible_r_free_col = false;
      int refmac_count = g.molecules[imol_coords].Refmac_count();
      std::string refmac_count_str = coot::util::int_to_string(refmac_count);
      int imol_map_refmac = -1;
      bool diff_map_flag = true;

      if (g.cif_dictionary_filename_vec->size() > 0)
         cif_lib_filename = g.cif_dictionary_filename_vec->at(0);

      int ierr = g.molecules[imol_coords].write_pdb_file(pdb_in_filename);
      if (! ierr) {
         safe_python_command("import refmac");
         execute_refmac_real(pdb_in_filename, pdb_out_filename,
                             mtz_in_filename, mtz_out_filename,
                             cif_lib_filename,
                             fobs_col, sigfobs_col, r_free_col, sensible_r_free_col,
                             make_molecules_flag,
                             refmac_count_str,
                             g.swap_pre_post_refmac_map_colours_flag,
                             imol_map_refmac,
                             diff_map_flag,
                             false, "", "", "");
      }
   }

}
