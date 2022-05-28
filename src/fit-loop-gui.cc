
#include "Python.h"  // for the fitting function

#include "positioned-widgets.h"
#include "widget-from-builder.hh"
#include "c-interface-gtk-widgets.h"
#include "cc-interface-scripting.hh"

#include "graphics-info.h"


void
new_fill_combobox_with_coordinates_options(GtkWidget *combobox_molecule,
                                           GCallback callback_func,
                                           int imol_active) {
   graphics_info_t g;
   g.new_fill_combobox_with_coordinates_options(combobox_molecule, callback_func, imol_active);
}

void
fill_mutate_sequence_dialog_gtkbuilder_version() {
   // c.f. wrapped_create_mutate_sequence_dialog()

   auto get_active_molecule_index = [] () {
                                int imol = -1;
                                for (int i=0; i<graphics_info_t::n_molecules(); i++) {
                                   if (graphics_info_t::molecules[i].has_model()) {
                                      imol = i;
                                      break;
                                   }
                                }
                                return imol;
                             };

   auto clear_box = [] (GtkWidget *box) {
                       GList *dlist = gtk_container_get_children(GTK_CONTAINER(box));
                       GList *free_list = dlist;
                       while (dlist) {
                          GtkWidget *list_item = (GtkWidget *) (dlist->data);
                          gtk_widget_destroy(GTK_WIDGET(list_item));
                          dlist = dlist->next;
                       }
                       g_list_free(free_list);
                    };

   GtkWidget *hbox_mol   = widget_from_builder("mutate_sequence_hbox_for_molecule_combobox");
   GtkWidget *hbox_chain = widget_from_builder("mutate_sequence_hbox_for_molecule_chain_combobox_text");

   // clear out hbox_mol and hbox_chain
   //
   clear_box(hbox_mol);
   clear_box(hbox_chain);

   GtkWidget *combobox_molecule = gtk_combo_box_new();   // number and name
   GtkWidget *combobox_chain = gtk_combo_box_text_new(); // just the chain id

#if (GTK_MAJOR_VERSION == 3 && GTK_MINOR_VERSION == 94) || (GTK_MAJOR_VERSION == 4)
      // 20220528-PE FIXME box packing
#else
   gtk_box_pack_start(GTK_BOX(hbox_mol),   combobox_molecule, FALSE, FALSE, 6);
   gtk_box_pack_start(GTK_BOX(hbox_chain), combobox_chain,    FALSE, FALSE, 6);
#endif
   gtk_widget_show(combobox_molecule);
   gtk_widget_show(combobox_chain);

   graphics_info_t g;
   int imol = get_active_molecule_index();
   g.mutate_sequence_imol = imol;
   new_fill_combobox_with_coordinates_options(combobox_molecule, NULL, imol);

   g.fill_combobox_with_chain_options(combobox_chain, imol, NULL);

}

GtkWidget *
create_fit_loop_rama_search_dialog_gtkbuilder_version() {

   GtkWidget *dialog             = widget_from_builder("mutate_sequence_dialog");
   GtkWidget *label              = widget_from_builder("function_for_molecule_label");
   GtkWidget *method_frame       = widget_from_builder("loop_fit_method_frame");
   GtkWidget *mutate_ok_button   = widget_from_builder("mutate_sequence_ok_button");
   GtkWidget *fit_loop_ok_button = widget_from_builder("fit_loop_ok_button");
   GtkWidget *checkbutton        = widget_from_builder("mutate_sequence_do_autofit_checkbutton");
   GtkWidget *rama_checkbutton   = widget_from_builder("mutate_sequence_use_ramachandran_restraints_checkbutton");

   set_transient_and_position(COOT_MUTATE_RESIDUE_RANGE_WINDOW, dialog);
   fill_mutate_sequence_dialog_gtkbuilder_version();

   gtk_label_set_text(GTK_LABEL(label), "\nFit loop in Molecule:\n");
   gtk_widget_hide(mutate_ok_button);
   gtk_widget_hide(checkbutton);
   gtk_widget_show(fit_loop_ok_button);
   gtk_widget_show(rama_checkbutton);
   gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(rama_checkbutton), TRUE);

   gtk_widget_show(method_frame);

   return dialog;
}


void
fit_loop_using_dialog()  {

   auto get_first_child = [] (GtkWidget *box) {
                             GtkWidget *first_child = 0;
                             GList *dlist = gtk_container_get_children(GTK_CONTAINER(box));
                             first_child = static_cast<GtkWidget *>(dlist->data);
                             g_list_free(dlist);
                             return first_child;
                          };

   auto get_sequence = [] () {
                          std::string seq;
                          GtkWidget *text = widget_from_builder("mutate_molecule_sequence_text");
                          GtkTextView *tv = GTK_TEXT_VIEW(text);
                          GtkTextBuffer* tb = gtk_text_view_get_buffer(tv);
                          GtkTextIter startiter;
                          GtkTextIter enditer;
                          gtk_text_buffer_get_iter_at_offset(tb, &startiter, 0);
                          gtk_text_buffer_get_iter_at_offset(tb, &enditer, -1);
                          char *txt = gtk_text_buffer_get_text(tb, &startiter, &enditer, 0);
                          if (txt) seq = std::string(txt);
                          return seq;
                       };

   std::cout << ":::::::::::::::::::::: read the gui, fit the loop! " << std::endl;

   GtkWidget *entry_1 = widget_from_builder("mutate_molecule_resno_1_entry");
   GtkWidget *entry_2 = widget_from_builder("mutate_molecule_resno_2_entry");

   const gchar *entry_1_text = gtk_entry_get_text(GTK_ENTRY(entry_1));
   const gchar *entry_2_text = gtk_entry_get_text(GTK_ENTRY(entry_2));

   GtkWidget *hbox_mol   = widget_from_builder("mutate_sequence_hbox_for_molecule_combobox");
   GtkWidget *hbox_chain = widget_from_builder("mutate_sequence_hbox_for_molecule_chain_combobox_text");

   try {
      graphics_info_t g;
      int resno_1 = coot::util::string_to_int(entry_1_text);
      int resno_2 = coot::util::string_to_int(entry_2_text);
      GtkWidget *molecule_combobox   = get_first_child(hbox_mol);
      GtkWidget *chain_combobox_text = get_first_child(hbox_chain);
      int imol = g.combobox_get_imol(GTK_COMBO_BOX(molecule_combobox));
      std::string chain_id = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(chain_combobox_text));
      GtkWidget *checkbutton_fit  = widget_from_builder("mutate_sequence_do_autofit_checkbutton");
      GtkWidget *checkbutton_rama = widget_from_builder("mutate_sequence_use_ramachandran_restraints_checkbutton");
      bool autofit_flag  = false;
      bool use_rama_flag = false;
      if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkbutton_fit)))  autofit_flag  = true;
      if (gtk_toggle_button_get_active(GTK_TOGGLE_BUTTON(checkbutton_rama))) use_rama_flag = true;
      std::string sequence = get_sequence();
      
      if (g.is_valid_model_molecule(imol)) {
         int res_number_delta = resno_2 - resno_1;
         int seq_length = sequence.length();
         if ((res_number_delta+1) == seq_length) {
            // good
         } else {
            sequence.clear();
            for (int i=0; i<(resno_2 - resno_1 + 1); i++)
               sequence += "A";
         }
         short int state_lang = coot::STATE_PYTHON;
         safe_python_command("import gap");
         std::vector<std::string> cmd_strings;
         cmd_strings.push_back("gap.fit_gap"); // was just "fit-gap" - safe_scheme_command will have to deal with that.
         cmd_strings.push_back(graphics_info_t::int_to_string(imol));
         cmd_strings.push_back(coot::util::single_quote(chain_id));
         cmd_strings.push_back(std::to_string(resno_1));
         cmd_strings.push_back(std::to_string(resno_2));
         cmd_strings.push_back(coot::util::single_quote(sequence));
         cmd_strings.push_back(std::to_string(static_cast<int>(use_rama_flag)));
         std::string cmd = g.state_command(cmd_strings, state_lang);
         safe_python_command(cmd);
      }
   }
   catch (const std::runtime_error & rte) {
      std::cout << "WARNING:: " << rte.what() << std::endl;
   }

}
