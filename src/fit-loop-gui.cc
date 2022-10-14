
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

#if (GTK_MAJOR_VERSION >= 4)
#else
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
#endif

   // GtkWidget *combobox_molecule = gtk_combo_box_new();   // number and name
   // GtkWidget *combobox_chain = gtk_combo_box_text_new(); // just the chain id

   GtkWidget *molecule_combobox   = widget_from_builder("mutate_sequence_molecule_combobox");
   GtkWidget *chain_combobox_text = widget_from_builder("mutate_sequence_chain_combobox_text");

   graphics_info_t g;
   int imol = get_active_molecule_index();
   std::cout << "debug::active index is " << imol << std::endl;
   g.mutate_sequence_imol = imol;
   g.new_fill_combobox_with_coordinates_options(molecule_combobox, NULL, imol);

   g.fill_combobox_with_chain_options(chain_combobox_text, imol, NULL);

}

GtkWidget *
create_fit_loop_rama_search_dialog_gtkbuilder_version() {

   GtkWidget *dialog              = widget_from_builder("mutate_sequence_dialog");
   GtkWidget *label               = widget_from_builder("function_for_molecule_label");
   GtkWidget *method_frame        = widget_from_builder("loop_fit_method_frame");
   GtkWidget *mutate_ok_button    = widget_from_builder("mutate_sequence_ok_button");
   GtkWidget *fit_loop_ok_button  = widget_from_builder("fit_loop_ok_button");
   GtkWidget *autofit_checkbutton = widget_from_builder("mutate_sequence_do_autofit_checkbutton");
   GtkWidget *rama_checkbutton    = widget_from_builder("mutate_sequence_use_ramachandran_restraints_checkbutton");

   set_transient_and_position(COOT_MUTATE_RESIDUE_RANGE_WINDOW, dialog);
   fill_mutate_sequence_dialog_gtkbuilder_version();

   gtk_label_set_text(GTK_LABEL(label), "\nFit loop in Molecule:\n");
   gtk_widget_hide(mutate_ok_button);
   gtk_widget_hide(autofit_checkbutton);
   gtk_widget_show(fit_loop_ok_button);
   gtk_widget_show(rama_checkbutton);
   gtk_check_button_set_active(GTK_CHECK_BUTTON(rama_checkbutton), TRUE);

   gtk_widget_show(method_frame);

   return dialog;
}


void
fit_loop_using_dialog()  {

   auto get_first_child = [] (GtkWidget *box) {
#if (GTK_MAJOR_VERSION >= 4)
      return gtk_widget_get_first_child(box);
#else
                             GtkWidget *first_child = 0;
                             GList *dlist = gtk_container_get_children(GTK_CONTAINER(box));
                             first_child = static_cast<GtkWidget *>(dlist->data);
                             g_list_free(dlist);
                             return first_child;
#endif
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

   std::cout << ":::::::::::::::::::::: fit_loop_using_dialog() read the gui, fit the loop! " << std::endl;

   GtkWidget *entry_1 = widget_from_builder("mutate_molecule_resno_1_entry");
   GtkWidget *entry_2 = widget_from_builder("mutate_molecule_resno_2_entry");

   const gchar *entry_1_text = gtk_editable_get_text(GTK_EDITABLE(entry_1));
   const gchar *entry_2_text = gtk_editable_get_text(GTK_EDITABLE(entry_2));

   // why do I want these?
   // GtkWidget *hbox_mol   = widget_from_builder("mutate_sequence_hbox_for_molecule_combobox");
   // GtkWidget *hbox_chain = widget_from_builder("mutate_sequence_hbox_for_molecule_chain_combobox_text");

   try {
      graphics_info_t g;
      int resno_1 = coot::util::string_to_int(entry_1_text);
      int resno_2 = coot::util::string_to_int(entry_2_text);

      GtkWidget *molecule_combobox   = widget_from_builder("mutate_sequence_molecule_combobox");
      GtkWidget *chain_combobox_text = widget_from_builder("mutate_sequence_chain_combobox_text");

      std::cout << "debug:: molecule_combobox: " << molecule_combobox << std::endl;
      std::cout << "debug:: chain_combobox text: " << chain_combobox_text << std::endl;

      if (! molecule_combobox)   { std::cout << "ERROR:: bad molecule_combobox lookup " << std::endl; return; }
      if (! chain_combobox_text) { std::cout << "ERROR:: bad chain combobox lookup "    << std::endl; return; }

      int imol = g.combobox_get_imol(GTK_COMBO_BOX(molecule_combobox));

      if (imol == -1) {
         std::cout << "ERROR:: bad imol " << imol << std::endl;
         return;
      }

      std::cout << "debug: imol " << imol << std::endl;

      std::string chain_id = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(chain_combobox_text));

      std::cout << "debug:: chain_id " << chain_id << std::endl;

      GtkWidget *checkbutton_rama = widget_from_builder("mutate_sequence_use_ramachandran_restraints_checkbutton");
      bool use_rama_flag = false;
      if (gtk_check_button_get_active(GTK_CHECK_BUTTON(checkbutton_rama))) use_rama_flag = true;
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

         // Old
         // safe_python_command("import gap");
         // std::vector<std::string> cmd_strings;
         // cmd_strings.push_back("gap.fit_gap"); // was just "fit-gap" - safe_scheme_command will have to deal with that.
         // cmd_strings.push_back(graphics_info_t::int_to_string(imol));
         // cmd_strings.push_back(coot::util::single_quote(chain_id));
         // cmd_strings.push_back(std::to_string(resno_1));
         // cmd_strings.push_back(std::to_string(resno_2));
         // cmd_strings.push_back(coot::util::single_quote(sequence));
         // cmd_strings.push_back(std::to_string(static_cast<int>(use_rama_flag)));
         // std::string cmd = g.state_command(cmd_strings, state_lang);
         // safe_python_command(cmd);

         // New 20221010-PE
         std::vector<coot::command_arg_t> args;
         args.push_back(coot::command_arg_t(imol));
         args.push_back(chain_id);
         args.push_back(resno_1);
         args.push_back(resno_2);
         args.push_back(sequence);
         args.push_back(static_cast<int>(use_rama_flag));
         std::string sc = g.state_command("gap", "fit_gap", args, state_lang); // coot.gap
         safe_python_command("import gap"); // should be coot.gap with correct Python namespacing FIXME

         std::cout << ":::::::::::::: " << sc << std::endl;
         safe_python_command(sc);

      }
   }
   catch (const std::runtime_error & rte) {
      std::cout << "WARNING:: " << rte.what() << std::endl;
   }

}
