/* src/dynaram-main.cc
 *
 * Author: Bernhard Lohkamp
 * Copyright 2015
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

#if !defined _MSC_VER
#include <unistd.h> // so that getopt declaration in /usr/include/x86_64-linux-gnu/bits/getopt_core.h
                    // comes before Coot's one
#endif

#ifdef HAVE_GOOCANVAS

#ifdef __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#else
#define __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#undef __GNU_LIBRARY__
#endif

#include <iostream>

#include <sys/stat.h>

#include <gtk/gtk.h>

#ifdef HAVE_GOOCANVAS
#include <goocanvas.h>
#endif

#include "utils/coot-utils.hh"
#include "rama_plot.hh"
#include "coords/mmdb-extras.h"
#include "coords/mmdb.hh"
#include <mmdb2/mmdb_manager.h>

// #include "coot-surface/rgbreps.h"


// Dummy definitions for stand alone version
extern "C" {
   void set_dynarama_is_displayed(GtkWidget *dynarama_widget, int imol) {}
   void set_go_to_atom_molecule(int imol) {}
   int set_go_to_atom_chain_residue_atom_name(const char *t1_chain_id, int iresno, const char *t3_atom_name)
   { return 0;}
   short int is_valid_model_molecule(int imol) { return 0;}
   void set_moving_atoms(double phi, double psi) {}
   void accept_phi_psi_moving_atoms() {}
   void clear_moving_atoms_object() {}
}


// void accept_phi_psi_moving_atoms() {}
// void set_dynarama_is_displayed(GtkWidget *dynarama_widget, int imol) {}


void
print_help(std::string cmd) {
   std::cout << "Usage: " << cmd << "\n"
             << "(--pdbin) pdb-in-filename\n"
             << "[--selection atom-selection-string]\n"
             << "[--chain chain-id]\n"
             << "[--chain2 chain-id2] (for kleywegt plot)\n"
             << "[--selection2 atom-selection-string] (for kleywegt plot)\n"
             << "[--pdbin2 pdb-in-filename2 (for kleywegt plot, otherwise assume pdbin)]\n"
             << "[--kleywegt (to make kleywegt, autoamtically for multiple selections, chains)]\n"
             << "[--edit (edit mode, currently debug only)]\n"
             << "[--psiaxis (change psi axis to -120 to 240)]\n"
             << "[--help (this help)]\n"
             << "\n";
   std::cout << "     where pdbin is the protein and pdbin2 a second one for a Kleywegt plot.\n";

}

int
main(int argc, char *argv[]) {


   std::string pdb_file_name;
   std::string pdb_file_name2;
   std::string selection;
   std::string selection2;
   std::string chain_id;
   std::string chain_id2;
   selection = "";
   selection2 = "";
   int index;
   int edit_res_no = -9999;
   int n_used_args = 0;
   int is_kleywegt_plot_flag = 0;
   int psi_axis_option = coot::rama_plot::PSI_CLASSIC;
   float block_size = 2;
   bool do_help = false;
   bool print_bg = false;

   const char *optstr = "i:s:j:t:c:d:e:b:kphx";
   struct option long_options[] = {
   {"pdbin", 1, 0, 0},
   {"selection", 1, 0, 0},
   {"pdbin2", 1, 0, 0},
   {"selection2", 1, 0, 0},
   {"chain", 1, 0, 0},
   {"chain2", 1, 0, 0},
   {"edit", 1, 0, 0},   // BL Note:: maybe there should be an edit selection
   {"blocksize", 1, 0, 0},
   {"kleywegt", 0, 0, 0},
   {"psiaxis", 0, 0, 0},
   {"help", 0, 0, 0},
   {0, 0, 0, 0}
   };

   int ch;
   int option_index = 0;
   while ( -1 !=
           (ch = coot_getopt_long(argc, argv, optstr, long_options, &option_index))) {
      switch(ch) {

      case 0:
         if (coot_optarg) {
            std::string arg_str = long_options[option_index].name;

            if (arg_str == "pdbin") {
               pdb_file_name = coot_optarg;
               n_used_args += 2;
            }
            if (arg_str == "pdbin2") {
               pdb_file_name2 = coot_optarg;
               is_kleywegt_plot_flag = 1;
               n_used_args += 2;
            }
            if (arg_str == "selection") {
               selection = coot_optarg;
               n_used_args += 2;
            }
            if (arg_str == "selection2") {
               selection2 = coot_optarg;
               is_kleywegt_plot_flag = 1;
               n_used_args += 2;
            }
            if (arg_str == "chain") {
               chain_id = coot_optarg;
               n_used_args += 2;
            }
            if (arg_str == "chain2") {
               chain_id2 = coot_optarg;
               is_kleywegt_plot_flag = 1;
               n_used_args += 2;
            }
            if (arg_str == "edit") {
               edit_res_no = coot::util::string_to_int(coot_optarg);
               n_used_args += 2;
            }
            if (arg_str == "blocksize") {
               block_size = coot::util::string_to_float(coot_optarg);
               n_used_args += 2;
            }
         } else {

            // options without arguments:

            // long argument without parameter:
            std::string arg_str(long_options[option_index].name);

            if (arg_str == "kleywegt") {
               is_kleywegt_plot_flag = 1;
               n_used_args++;
            }
            if (arg_str == "psiaxis") {
               psi_axis_option = coot::rama_plot::PSI_MINUS_120;
               n_used_args++;
            }
            if (arg_str == "help") {
               print_help(argv[0]);
               do_help = true;
               n_used_args++;
            }
         }
         break;

      case 'i':
         pdb_file_name = coot_optarg;
         n_used_args += 2;
         break;

      case 'j':
         pdb_file_name2 = coot_optarg;
         n_used_args += 2;
         break;

      case 's':
         selection = coot_optarg;
         n_used_args += 2;
         break;

      case 't':
         selection2 = coot_optarg;
         n_used_args += 2;
         break;

      case 'c':
         chain_id = coot_optarg;
         n_used_args += 2;
         break;

      case 'd':
         chain_id2 = coot_optarg;
         n_used_args += 2;
         break;

      case 'e':
         edit_res_no = coot::util::string_to_int(coot_optarg);
         n_used_args += 2;
         break;

      case 'b':
         block_size = coot::util::string_to_float(coot_optarg);
         n_used_args += 2;
         break;

      case 'k':
         is_kleywegt_plot_flag = 1;
         n_used_args++;
         break;

      case 'p':
         psi_axis_option = coot::rama_plot::PSI_MINUS_120;
         n_used_args++;
         break;

      case 'h':
         print_help(argv[0]);
         n_used_args++;
         do_help = true;
         break;

      // This is used to make background images
      case 'x':
         print_help(argv[0]);
         n_used_args++;
         print_bg = true;
         do_help = false;
         break;

      case '?':
         std::cout << "Unrecognised option: " << optopt << std::endl;
         break;

      default:
         std::cout << "default coot_optarg: " << coot_optarg << std::endl;
         break;
      }
   }

   // any unhandled arguments which could be a pdb?
   for (index = n_used_args+1; index < argc; index++) {
      if (coot::util::extension_is_for_coords(coot::util::file_name_extension(argv[index])))
         pdb_file_name = argv[index];
      else {
         std::cout <<"BL INFO:: have an unknown arg: " << argv[index] <<std::endl;
         print_help(argv[0]);
         do_help = true;
      }
   }


   if (not do_help) {

      mmdb::Manager *mol = new mmdb::Manager();
      mmdb::Manager *mol2 = new mmdb::Manager();

      if (pdb_file_name.length() == 0) {
         std::cout << "WARNING:: Missing input PDB file\n";
         //exit(1);

      } else {
         mol->ReadPDBASCII(pdb_file_name.c_str());
      }

      gtk_init (&argc, &argv);

      float level_prefered = 0.02;
      float level_allowed = 0.002;
      int imol = 0; // dummy for now
      int imol2 = 0;

      // edit plot?
      if (edit_res_no > -9999) {
         // make an edit plot
         coot::rama_plot *edit_phi_psi_plot = new coot::rama_plot;
         edit_phi_psi_plot->init("phi/psi-edit", psi_axis_option);
         edit_phi_psi_plot->set_stand_alone();

         if (chain_id.size() == 0)
            chain_id = "A";
         // make a mmdb res and the use
         mmdb::PResidue *SelResidue;
         int nRes;
         int SelHnd;
         std::vector <coot::util::phi_psi_t> vp;
         //g_print("BL DEBUG:: edit plot with chain %s and resno %i\n", chain_id.c_str(), edit_res_no);
         for (int resno=edit_res_no; resno<edit_res_no+2; resno++) {
            SelHnd = mol->NewSelection();
            mol->Select(SelHnd,
                        mmdb::STYPE_RESIDUE,  // select residues
                        0,              // any model
                        chain_id.c_str(),          // chains "A" and "B" only
                        resno-1,"*",resno+1,"*", // any residue in range 30 to 100,
                        // any insertion code
                        "*",            // any residue name
                        "*",            // any atom name
                        "*",           // any chemical element but sulphur
                        "*",            // any alternative location indicator
                        mmdb::SKEY_NEW);         // OR-selection
            mol->GetSelIndex(SelHnd, SelResidue, nRes);
            if (nRes == 3) {
               std::pair<bool, coot::util::phi_psi_t> phi_psi_all = coot::util::get_phi_psi(SelResidue);
               if (phi_psi_all.first)
                  vp.push_back(phi_psi_all.second);
            } else {
               std::cout<<"BL WARNING:: not 3 residues in selection"<<std::endl;
               break;
            }
            mol->DeleteSelection(SelHnd);
         }
         if (vp.size() == 2)
            edit_phi_psi_plot->draw_it(vp);
         else
            g_print("BL WARNING:: problem making phi psi for 2 residues %i and %i\n",
                    edit_res_no, edit_res_no+1);

      } else {
         coot::rama_plot *rama = new coot::rama_plot;
         // normal rama (or kleywy)
         if (is_kleywegt_plot_flag) {
            mol2 = mol;
            imol2 = imol + 1;
            if (chain_id.size() == 0 && chain_id2.size() == 0 &&
                selection.size() == 0 && selection2.size()== 0) {
               chain_id = "A";
               chain_id2 = "B";
            }
         }

         if (pdb_file_name2.length() > 0) {
            // can be same mol as well, usually...., so maybe another flag
            is_kleywegt_plot_flag = 1;
            mol2->ReadPDBASCII(pdb_file_name2.c_str());
            imol2 = imol + 1;
            // Make a double name
            pdb_file_name = pdb_file_name + " vs. \n" + pdb_file_name2;
            // what about selections? dont care if given ok,
            // if not then compare Rama of 2 structures
         }

         // Now select mols and pass the handle
         // Alternatively: pass the selection string
         // probably better to pass the Handle
         // OR even make a new mol and only work with these

         int selHnd = -1;
         int selHnd2 = -1;
         int nRes;
         mmdb::PResidue *SelResidue;
         mmdb::PChain *SelChain;

         if (selection.size() > 0) {
            selHnd = mol->NewSelection();
            //selection="//A";
            mol->Select(selHnd,
                        mmdb::STYPE_RESIDUE,
                        selection.c_str(),
                        mmdb::SKEY_NEW);
            mol->GetSelIndex(selHnd, SelResidue, nRes);
         } else {
            if (chain_id.size() > 0) {
               selHnd = mol->NewSelection();
               //selection="//A";
               mol->Select(selHnd,
                           mmdb::STYPE_RESIDUE,
                           0,              // any model
                           chain_id.c_str(),          // chains
                           mmdb::ANY_RES, "*",
                           mmdb::ANY_RES, "*", // all residues
                           // any insertion code
                           "*",            // any residue name
                           "*",            // any atom name
                           "*",           // any chemical element but sulphur
                           "*",            // any alternative location indicator
                           mmdb::SKEY_NEW);         // OR-selection
               mol->GetSelIndex(selHnd, SelResidue, nRes);
            }
         }
         if (selection2.size() > 0) {
            if (mol2) {
               selHnd2 = mol2->NewSelection();
               //selection="//A";
               mol2->Select(selHnd2,
                            mmdb::STYPE_RESIDUE,
                            selection2.c_str(),
                            mmdb::SKEY_NEW);
            } else {
               g_print("BL INFO:: no mol2, no 2nd selection.");
            }
         }

         rama->set_stand_alone();
         rama->init(imol,
                    pdb_file_name,
                    level_prefered,
                    level_allowed,
                    block_size,
                    is_kleywegt_plot_flag,
                    psi_axis_option);
         if (rama->dynawin) {
            if (is_kleywegt_plot_flag) {
               gtk_widget_hide(rama->selection_hbox);
               if (selHnd > -1 && selHnd2 > -1) {
                  rama->draw_it(imol, imol2,
                                mol, mol2,
                                selHnd, selHnd2);
               } else {
                  if (chain_id.size() != 0 && chain_id2.size() != 0)  {
                     rama->draw_it(imol, imol2,
                                   mol, mol2,
                                   chain_id, chain_id2);
                  } else {
                     if (mol != mol2) {
                        rama->draw_it(imol, imol2,
                                      mol, mol2);
                     } else {
                        g_print("BL INFO:: have different imols but same selection. Should not happen. No idea what to do");
                     }
                  }
               }
            } else {
               if (selHnd > -1) {
                  rama->draw_it(mol, selHnd, 1);
                  // FIXME:: maybe for chains too?! And for kleywegt at some point
                  if (selection.size() > 0) {
                     //rama->show_selection_widget(1);
                     gtk_toggle_button_set_active(GTK_TOGGLE_BUTTON(rama->selection_checkbutton),
                                                  TRUE);
                     gtk_entry_set_text(GTK_ENTRY(rama->selection_entry),
                                        selection.c_str());
                  }
               }
               else {
                  if (mol) {
                     rama->draw_it(mol);
                  } else {
                     g_print("BL INFO:: no mol and no selection, so no plot. Sorry.");
                  }
               }
            }

            if (rama) {
               if (print_bg) {
                  rama->make_bg_images();
               }
            }
         } else {
            g_print("BL WARNING:: something went wrong initialising the rama plot\n");
            return 1;
         }
      }

      gtk_main ();
      if (mol)
         delete mol;
      if (mol2)
         delete mol2;

      gtk_init(&argc, &argv);
   }

   std::cout << "BL DBEUG:: returning now n_used_args" << n_used_args <<std::endl;

   return 0;

}

#else

#include <gtk/gtk.h>
// Dummy definitions for stand alone version
extern "C" {
   void set_dynarama_is_displayed(GtkWidget *dynarama_widget, int imol) {}
   void set_go_to_atom_molecule(int imol) {}
   int set_go_to_atom_chain_residue_atom_name(const char *t1_chain_id, int iresno, const char *t3_atom_name)
   { return 0;}
   short int is_valid_model_molecule(int imol) { return 0;}
   void set_moving_atoms(double phi, double psi) {}
   void accept_phi_psi_moving_atoms() {}
   void clear_moving_atoms_object() {}
}

#include <iostream>

int
main(int argc, char *argv[]) {

   std::cout << "No goo canvas at compile-time, no dynarama " << std::endl;
   return 0;
}



#endif

