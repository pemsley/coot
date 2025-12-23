/* src/graphics-info.cc
 *
 * Copyright 2002, 2003, 2004, 2005 by The University of York
 * Copyright 2008  by The University of Oxford
 * Copyright 2016 by Medical Research Council
 * Author: Paul Emsley
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */

#include "validation-graphs/validation-graphs.hh"
#include "widget-from-builder.hh"
#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

#include "compat/coot-sysdep.h"


#ifndef HAVE_STRING
#define HAVE_STRING
#include <string>
#endif

#ifndef HAVE_VECTOR
#define HAVE_VECTOR
#include <vector>
#endif

#include <gtk/gtk.h>  // must come after mmdb_manager on MacOS X Darwin
// #include <GL/glut.h>  // for some reason...  // Eh?

#include <iostream>

#include <mmdb2/mmdb_manager.h>

#include "coords/mmdb-extras.hh"
#include "coords/mmdb.hh"
#include "coords/mmdb-crystal.hh"
#include "coords/Bond_lines.hh"

#include "clipper/core/map_utils.h" // Map_stats
#include "skeleton/graphical_skel.h"


#include "graphics-info.h"
#include "interface.h"

#include "molecule-class-info.h"
// #include "rama_plot.hh"
#include "skeleton/BuildCas.h"


#include "coot-utils/gl-matrix.h" // for baton rotation
#include "analysis/bfkurt.hh"

#include "globjects.h"
#ifdef USE_DUNBRACK_ROTAMERS
#include "ligand/dunbrack.hh"
#else
#include "ligand/richardson-rotamer.hh"
#endif
#include "ligand/ligand.hh"

#include "coot-utils/coot-map-utils.hh"
#include "geometry-graphs.hh"


void graphics_info_t::refresh_validation_graph_model_list() {

   // g_debug("refresh_validation_graph_model_list() called.");

   // std::cout << "-------------- refresh_validation_graph_model_list ------- " << std::endl;

   gtk_tree_model_foreach(
                          GTK_TREE_MODEL(validation_graph_model_list),
                          +[](GtkTreeModel* model, GtkTreePath* path, GtkTreeIter* iter, gpointer data) -> gboolean {
                             GtkListStore* list = GTK_LIST_STORE(model);
                             return ! gtk_list_store_remove(list,iter);
                          },
                          NULL
                          );
   int idx_active = -1; // use this to set active_validation_graph_model_idx
   for(int i=0; i<graphics_info_t::n_molecules(); i++) {
      if (graphics_info_t::molecules[i].has_model()) {
         std::string label = graphics_info_t::molecules[i].dotted_chopped_name();
         GtkTreeIter iter;
         // std::cout << "----- refresh_validation_graph_model_list adding label " << label << std::endl;
         gtk_list_store_append(validation_graph_model_list, &iter);
         gtk_list_store_set(validation_graph_model_list, &iter, 0, label.c_str(), 1, i, -1);
         if (idx_active  == -1)
            idx_active = i;
      }
   }
   if (idx_active != -1)
      active_validation_graph_model_idx = idx_active;

   // g_warning("refresh_validation_graph_model_list(): todo: Check if the active model ID is still on the list and react appropriately");
   // if (model is no longer on the list) {
   // 	// destroy all opened validation graphs (via calls to destroy_validation_graph())
   // }

   if (idx_active != -1) {
      if (!is_valid_model_molecule(active_validation_graph_model_idx)) {
         if (false)
            std::cout << "TODO:: in refresh_validation_graph_model_list() Destroy graphs for model "
                      << active_validation_graph_model_idx << " here..." << std::endl;
      }
   }
}

void graphics_info_t::update_active_validation_graph_model(int model_idx) {

   // this happens when the user changes the active model in the model combobox in the validation graph dialog

   // 1. Update the model active model variable
   active_validation_graph_model_idx = model_idx;
   std::cout << "update_active_validation_graph_model() active_validation graph model idx"
             << active_validation_graph_model_idx << std::endl;
   // 2. Handle chains
   g_warning("todo: update_active_validation_graph_model(): handle chains");
   // 3. Recompute all validation data of active validation graphs (by looking up widgets, not the data) and trigger a redraw
   for (const std::pair<const coot::validation_graph_type,GtkWidget*>& i : validation_graph_widgets) {
      g_warning("Todo: Display/rebuild validation graph data for: %s [model index changed to %i]",
                coot::validation_graph_type_to_human_name(i.first).c_str(), model_idx);
      coot::validation_graph_type graph_type = i.first;
      GtkWidget *graph = i.second;
      if (graph_type == coot::validation_graph_type::density_fit) { }
      if (graph_type == coot::validation_graph_type::omega) { }
      if (graph_type == coot::validation_graph_type::rama) { }
      if (graph_type == coot::validation_graph_type::rota) { }

   }
}

void graphics_info_t::change_validation_graph_chain(const std::string& chain_id) {

   // 20230527-PE It will be a while before this gets filled I think!
   g_debug("Todo: change_validation_graph_chain");
}


void graphics_info_t::refresh_ramachandran_plot_model_list() {

   // what is this - I mean, who calls it/when does it run? Is this an old method now that we have rama_plot_boxes?

   // noise
   // std::cout << "----------------------- refresh_ramachandran_plot_model_list --------- " << std::endl;

   auto fn = +[] (GtkTreeModel* model, GtkTreePath* path, GtkTreeIter* iter, gpointer data) {
      GtkListStore* list = GTK_LIST_STORE(model);
      return gboolean(!gtk_list_store_remove(list,iter));
   };

   gtk_tree_model_foreach(GTK_TREE_MODEL(ramachandran_plot_model_list), fn, NULL);

   for(int i=0; i<graphics_info_t::n_molecules(); i++) {
      if (graphics_info_t::molecules[i].has_model()) {
         std::string label = graphics_info_t::molecules[i].dotted_chopped_name();
         GtkTreeIter iter;
         // std::cout << "----- refresh_ramachandran_plot_model_list adding label " << label << std::endl;
         gtk_list_store_append(ramachandran_plot_model_list, &iter);
         gtk_list_store_set(ramachandran_plot_model_list, &iter, 0, label.c_str(), 1, i, -1);
      }
   }
   // std::cout << "----------------------- done refresh_ramachandran_plot_model_list --------- " << std::endl;
}

// TODO: we're not using tabs right now. Do we ever intend to do so?

// void create_tab_for_validation_graph(coot::validation_graph_type type, GtkWidget* the_graph) {
// 	GtkWidget* notebook = widget_from_builder("validation_graph_notebook");
// 	// we assume that when this function is called, there is no tab for the graph type
// 	GtkWidget* sw = gtk_scrolled_window_new();
// 	gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(sw), the_graph);
// 	gtk_notebook_append_page(GTK_NOTEBOOK(notebook), sw, gtk_label_new(coot::validation_graph_type_to_human_name(type).c_str()));
// }

// void destroy_tab_for_validation_graph(coot::validation_graph_type type) {
// 	GtkWidget* notebook = widget_from_builder("validation_graph_notebook");
// 	auto find_tab_idx = [notebook](coot::validation_graph_type graph_type) -> int {
// 		std::string target_label = coot::validation_graph_type_to_human_name(graph_type);
// 		for(int i = 0; i < gtk_notebook_get_n_pages(GTK_NOTEBOOK(notebook));i++) {
// 			const char* page_label = gtk_notebook_get_tab_label_text(GTK_NOTEBOOK(notebook),gtk_notebook_get_nth_page(GTK_NOTEBOOK(notebook),i));
// 			if (!page_label) {
// 				g_error("NULL page label");
// 			}
// 			if (page_label == target_label) {
// 				return i;
// 			}
// 		}
// 		return -1;
// 	};
// 	auto idx = find_tab_idx(type);
// 	if (idx == -1) {
// 		g_warning("Failed to find tab for graph type: %s",coot::validation_graph_type_to_human_name(type).c_str());
// 	} else {
// 		gtk_notebook_remove_page(GTK_NOTEBOOK(notebook),idx);
// 	}
// }

void insert_validation_graph(GtkWidget* graph) {

   GtkWidget* target_box = widget_from_builder("main_window_validation_graph_box");
   if(! gtk_widget_get_first_child(target_box)) {
      // Empty validation_graph_box means that we need to make the validation_graph_frame visible first
      GtkWidget* frame = widget_from_builder("main_window_validation_graph_frame");
      gtk_widget_set_visible(frame, TRUE);
   }
   //g_debug("Inserting %p to the validation graph box.",graph);
   gtk_box_append(GTK_BOX(target_box), graph);

}

void remove_validation_graph(GtkWidget* graph) {

   // 20230424-PE we move this to the box in the paned in the main window
   GtkWidget* target_box = widget_from_builder("main_window_validation_graph_box");
   //g_debug("Removing %p from the validation graph box.",graph);
   gtk_box_remove(GTK_BOX(target_box), graph);
   if (! gtk_widget_get_first_child(target_box)) {
      // If the validation_graph_box is empty now, we need to make the validation_graph_frame invisible
      GtkWidget* frame = widget_from_builder("main_window_validation_graph_frame");
      gtk_widget_set_visible(frame, FALSE);
   }
}

#include "validation-graphs/validation-information.hh"
#include "validation-graphs/validation-graph-widget.hh"

coot::validation_information_t
get_validation_data_for_density_fit_analysis(int imol) {

   graphics_info_t g;  // remove this when added to the class

   coot::validation_information_t r;
   r.name = "Density fit analysis";
   r.type = coot::graph_data_type::Score;

   int imol_map = g.Imol_Refinement_Map();
   if (! g.is_valid_model_molecule(imol))   return r;
   if (! g.is_valid_map_molecule(imol_map)) return r;

   const clipper::Xmap<float> &xmap = g.molecules[imol_map].xmap;

   mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol; // this can be removed when in place
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int n_res = chain_p->GetNumberOfResidues();
         for (int ires=0; ires<n_res; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            if (residue_p) {
               coot::residue_spec_t res_spec(residue_p);
               res_spec.int_user_data = imol; // this is used in the residue block click callback
               mmdb::PAtom *residue_atoms=0;
               int n_residue_atoms = 0;
               residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
               if (n_residue_atoms > 0) {
                  double residue_density_score = coot::util::map_score(residue_atoms, n_residue_atoms, xmap, 1);
                  std::string l = "Chain ID: "+ res_spec.chain_id + "     Residue number: "+ std::to_string(res_spec.res_no);
                  std::string atom_name = coot::util::intelligent_this_residue_mmdb_atom(residue_p)->GetAtomName();
                  const std::string &chain_id = res_spec.chain_id;
                  int this_resno = res_spec.res_no;
                  coot::atom_spec_t atom_spec(chain_id, this_resno, res_spec.ins_code, atom_name, "");
                  double score_per_residue = residue_density_score / static_cast<double>(n_residue_atoms);
                  coot::residue_validation_information_t rvi(res_spec, atom_spec, score_per_residue, l);
                  r.add_residue_validation_information(rvi, chain_id);
               }
            }
         }
      }
   }
   r.set_min_max();
   return r;
}

coot::validation_information_t
get_validation_data_for_density_correlation_analysis(int imol) {

   graphics_info_t g;  // remove this when added to the class

   coot::validation_information_t vi;
   vi.name = "Density correlation analysis";
   vi.type = coot::graph_data_type::Correlation;

   int imol_map = g.Imol_Refinement_Map();
   if (! g.is_valid_model_molecule(imol))   return vi;
   if (! g.is_valid_map_molecule(imol_map)) return vi;

   const clipper::Xmap<float> &xmap = g.molecules[imol_map].xmap;

   mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol; // this can be removed when in place
   unsigned short int atom_mask_mode = 0;
   float atom_radius = 2.0;

   std::vector<coot::residue_spec_t> residue_specs;
   std::vector<mmdb::Residue *> residues = coot::util::residues_in_molecule(mol);
   for (unsigned int i=0; i<residues.size(); i++)
      residue_specs.push_back(coot::residue_spec_t(residues[i]));

   std::vector<std::pair<coot::residue_spec_t, float> > correlations =
      coot::util::map_to_model_correlation_per_residue(mol,
                                                       residue_specs,
                                                       atom_mask_mode,
                                                       atom_radius, // for masking
                                                       xmap);

   std::vector<std::pair<coot::residue_spec_t, float> >::const_iterator it;
   for (it=correlations.begin(); it!=correlations.end(); ++it) {
      const auto &r_spec(it->first);
      const auto &correl(it->second);

      auto res_spec = r_spec;
      res_spec.int_user_data = imol;
      std::string atom_name = " CA ";
      coot::atom_spec_t atom_spec(r_spec.chain_id, r_spec.res_no, r_spec.ins_code, atom_name, "");
      std::string label = "Correl: ";
      coot::residue_validation_information_t rvi(res_spec, atom_spec, correl, label);
      vi.add_residue_validation_information(rvi, r_spec.chain_id);

   }
   vi.set_min_max();
   return vi;
}


#include "coords/ramachandran-validation.hh"

coot::validation_information_t
get_validation_data_for_ramachandran_analysis(int imol) {

   coot::validation_information_t vi;
   vi.name = "Ramachandran analysis";
   vi.type = coot::graph_data_type::Probability;

   graphics_info_t g;
   if (! g.is_valid_model_molecule(imol)) return vi;

   mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;

   const ramachandrans_container_t rc;
   std::vector<coot::phi_psi_prob_t> rv = coot::ramachandran_validation(mol, rc);
   for (unsigned int i=0; i<rv.size(); i++) {
      std::string chain_id = rv[i].phi_psi.chain_id;
      coot::residue_spec_t residue_spec(rv[i].phi_psi.chain_id, rv[i].phi_psi.residue_number, rv[i].phi_psi.ins_code);
      residue_spec.int_user_data = imol;
      double pr = rv[i].probability;
      std::string label = rv[i].phi_psi.chain_id + std::string(" ") + std::to_string(rv[i].phi_psi.residue_number);
      if (! rv[i].phi_psi.ins_code.empty())
         label += std::string(" ") + rv[i].phi_psi.ins_code;
      coot::atom_spec_t atom_spec(residue_spec.chain_id, residue_spec.res_no, residue_spec.ins_code, " CA ", "");
      coot::residue_validation_information_t rvi(residue_spec, atom_spec, pr, label);
      vi.add_residue_validation_information(rvi, chain_id);
   }
   vi.set_min_max();
   return vi;
}

coot::validation_information_t
get_validation_data_for_rotamer_analysis(int imol) {

   coot::validation_information_t vi;
   vi.name = "Rotamer analysis";
   vi.type = coot::graph_data_type::Probability;

   graphics_info_t g;
   if (! g.is_valid_model_molecule(imol)) return vi;

   mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;

   // fill these
   mmdb::PResidue *SelResidues = 0;
   int nSelResidues = 0;

   int selHnd = mol->NewSelection(); // yes, it's deleted.
   int imod = 1; // multiple models don't work on validation graphs

   mol->Select(selHnd, mmdb::STYPE_RESIDUE, imod,
                        "*", // chain_id
                        mmdb::ANY_RES, "*",
                        mmdb::ANY_RES, "*",
                        "*",  // residue name
                        "*",  // Residue must contain this atom name?
                        "*",  // Residue must contain this Element?
                        "*",  // altLocs
                        mmdb::SKEY_NEW // selection key
                        );
   mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

   if (nSelResidues > 2) {

      for (int ir=0; ir<nSelResidues; ir++) {
         mmdb::Residue *residue_p = SelResidues[ir];
         coot::residue_spec_t res_spec(residue_p);
         res_spec.int_user_data = imol;
         mmdb::PAtom *residue_atoms=0;
         int n_residue_atoms;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);

         // double residue_density_score = coot::util::map_score(residue_atoms, n_residue_atoms, xmap, 1);

         if (n_residue_atoms > 5) {

            std::string res_name = residue_p->GetResName();
            if (true) {

               coot::rotamer rot(residue_p);
               coot::rotamer_probability_info_t rpi = rot.probability_of_this_rotamer();
               double prob = rpi.probability * 0.01; // to range 0->1

               std::string l = "Chain ID: " + res_spec.chain_id+"     Residue number: " + std::to_string(res_spec.res_no);
               std::string atom_name = coot::util::intelligent_this_residue_mmdb_atom(residue_p)->GetAtomName();
               const std::string &chain_id = res_spec.chain_id;
               int this_resno = res_spec.res_no;
               coot::atom_spec_t atom_spec(chain_id, this_resno, res_spec.ins_code, atom_name, "");
               coot::residue_validation_information_t rvi(res_spec, atom_spec, prob, l);
               vi.add_residue_validation_information(rvi, chain_id);
            }
         }
      }
      mol->DeleteSelection(selHnd);
   }
   vi.set_min_max();
   return vi;
}

#include "ideal/simple-restraint.hh"

coot::validation_information_t
get_validation_data_for_peptide_omega_analysis(int imol) {

   coot::validation_information_t vi;
   vi.name = "Peptide Omega analysis";
   vi.type = coot::graph_data_type::UNSET; // should it have a type?

   graphics_info_t g;
   const coot::protein_geometry &geom = *g.Geom_p();
   if (! g.is_valid_model_molecule(imol)) return vi;

   mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
   int imodel = 1;
   mmdb::Model *model_p = mol->GetModel(imodel);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         // std::cout << "ichain: " << ichain << std::endl;
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         // std::cout << "ichain: " << ichain << " " << chain_p << std::endl;
         std::string chain_id(chain_p->GetChainID());
         coot::restraints_container_t rc(g.molecules[imol].atom_sel, chain_id, nullptr);
         coot::omega_distortion_info_container_t odi = rc.omega_trans_distortions(geom, true);
         if (false) {
            std::cout << "odi: chain_id "  << odi.chain_id << std::endl;
            std::cout << "odi: min_resno " << odi.min_resno << std::endl;
            std::cout << "odi: max_resno " << odi.max_resno << std::endl;
            std::cout << "odi: n omega_distortions " << odi.omega_distortions.size() << std::endl;
         }

         coot::chain_validation_information_t cvi(chain_id);
         for (const auto &od : odi.omega_distortions) {
            coot::residue_spec_t res_spec(chain_id, od.resno, "");
            res_spec.int_user_data = imol;
            coot::atom_spec_t atom_spec(chain_id, od.resno, "", " CA ", "");
            std::string label = od.info_string;
            coot::residue_validation_information_t rvi(res_spec, atom_spec, od.distortion, label);
            cvi.add_residue_validation_information(rvi);
         }
         vi.cviv.push_back(cvi);
      }
   }
   vi.set_min_max();
   return vi;
}


coot::validation_information_t
get_validation_data_for_temperature_factor_analysis(int imol) {

   coot::validation_information_t vi;
   vi.name = "Temperature Factor analysis";
   vi.type = coot::graph_data_type::UNSET; // should it have a type?


   graphics_info_t g;
   if (! g.is_valid_model_molecule(imol)) return vi;

   mmdb::Manager *mol = g.molecules[imol].atom_sel.mol;
   bool is_shelx_mol = g.molecules[imol].is_from_shelx_ins();
   coot_extras::b_factor_analysis bfa(mol, is_shelx_mol);
   std::vector<coot_extras::my_chain_of_stats_t> bfa_chain_info = bfa.chain_details(); // bfkurt.hh

   for (const auto &chain : bfa_chain_info) {
      coot::chain_validation_information_t cvi(chain.chain_id);
      for (const auto &res_prop : chain.residue_properties) {
         coot::residue_spec_t res_spec(chain.chain_id, res_prop.resno, res_prop.inscode);
         coot::atom_spec_t   atom_spec(chain.chain_id, res_prop.resno, res_prop.inscode, res_prop.atom_name, "");
         res_spec.int_user_data = imol;
         std::string label = "Label";
         coot::residue_validation_information_t rvi(res_spec, atom_spec, res_prop.mean, label);
         cvi.add_residue_validation_information(rvi);
      }
      vi.cviv.push_back(cvi);
   }

   vi.set_min_max();
   return vi;
}

coot::validation_information_t
graphics_info_t::get_validation_data_for_geometry_analysis(int imol) {

   coot::validation_information_t vi;
   vi.name = "Geometry Distortion analysis";
   vi.type = coot::graph_data_type::Distortion;

   if (! is_valid_model_molecule(imol)) return vi;

   bool with_nbcs = false; // 20230417-PE for now
   std::vector<coot::geometry_distortion_info_container_t> gd =
      geometric_distortions_from_mol(imol, molecules[imol].atom_sel, with_nbcs);

   for (const auto &chain : gd) {
      coot::chain_validation_information_t cvi(chain.chain_id);
      std::map<coot::residue_spec_t, double> residue_distortion_sum_map;
      // we have a vector of restraints - what restraints are they, I wonder.
      for (const auto &rest : chain.geometry_distortion) {
         std::map<coot::residue_spec_t, double>::iterator it = residue_distortion_sum_map.find(rest.residue_spec);
         if (it == residue_distortion_sum_map.end()) {
            residue_distortion_sum_map[rest.residue_spec] = rest.distortion_score;
         } else {
            it->second += rest.distortion_score;
         }
      }

      std::map<coot::residue_spec_t, double>::iterator it;
      for (it=residue_distortion_sum_map.begin(); it != residue_distortion_sum_map.end(); ++it) {
         coot::atom_spec_t atom_spec; // not set (yet) - what is it used for?
         auto res_spec = it->first;
         res_spec.int_user_data = imol;
         std::string label = "Label";
         double distortion_sum = it->second;
         coot::residue_validation_information_t rvi(res_spec, atom_spec, distortion_sum, label);
         cvi.add_residue_validation_information(rvi);
      }
      vi.cviv.push_back(cvi);
   }
   // vi.set_min_max(); // this is auto-scaling, we don't want that.
   vi.min_max = coot::validation_information_min_max_t(0.0, 200.0);
   return vi;
}

coot::validation_information_t
get_validation_data(int imol, coot::validation_graph_type type) {

   graphics_info_t g;

   // types in validation-graphs.hh

   coot::validation_information_t vi;
   if (type == coot::validation_graph_type::density_fit)
      vi = get_validation_data_for_density_fit_analysis(imol);
   if (type == coot::validation_graph_type::density_correlation)
      vi = get_validation_data_for_density_correlation_analysis(imol);
   if (type == coot::validation_graph_type::rama)
      vi = get_validation_data_for_ramachandran_analysis(imol);
   if (type == coot::validation_graph_type::rota)
      vi = get_validation_data_for_rotamer_analysis(imol);
   if (type == coot::validation_graph_type::temp_factor)
      vi = get_validation_data_for_temperature_factor_analysis(imol);
   if (type == coot::validation_graph_type::omega)
      vi = get_validation_data_for_peptide_omega_analysis(imol);
   if (type == coot::validation_graph_type::geometry)
      vi = g.get_validation_data_for_geometry_analysis(imol);

   return vi;

}


// pass active_validation_graph_model_idx as imol
// static
void
graphics_info_t::create_validation_graph(int imol, coot::validation_graph_type type) {

   std::cout << "Yes! create_validation_graph() for " << imol << " type: "
             << coot::validation_graph_type_to_human_name(type) << std::endl;

   if (imol != -1) {
      // 1. instantiate the validation graph
      CootValidationGraph *cvg = coot_validation_graph_new();
      GtkWidget *this_is_the_graph = GTK_WIDGET(cvg);
      // 3. Compute data
      coot::validation_information_t vi = get_validation_data(imol, type);
      validation_graph_widgets[type] = this_is_the_graph;
      // 4. Store the data in std::maps
      std::shared_ptr<coot::validation_information_t> vip = std::make_shared<coot::validation_information_t>(vi);
      validation_graph_data[type] = vip;
      // 5. Set the data for the graph
      coot_validation_graph_set_validation_information(cvg, vip);
      // 6. Show the graph
      insert_validation_graph(this_is_the_graph);

      auto callback = +[] (G_GNUC_UNUSED CootValidationGraph* self,
                           const coot::residue_validation_information_t* residue_vip,
                           G_GNUC_UNUSED gpointer user_data) {

         std::cout << "residue-clicked handler " << residue_vip->label << " " << residue_vip->residue_spec << std::endl;
         int imol = residue_vip->residue_spec.int_user_data; // set by constructor of the validation information
         graphics_info_t g;
         g.go_to_residue(imol, residue_vip->residue_spec);
      };

      g_signal_connect(cvg, "residue-clicked", G_CALLBACK(callback), nullptr);

   } else {
      g_warning("graphics_info_t::create_validation_graph(): There is no valid active validation graph model.");
   }
}

// see update_validation_graphs(imol) below


void
graphics_info_t::destroy_validation_graph(coot::validation_graph_type type) {

   // 1. Remove the graph and its data from std::maps
   auto* widget = validation_graph_widgets[type];
   validation_graph_widgets.erase(type);
   validation_graph_data.erase(type);
   // 2. Destroy the graph widget
   remove_validation_graph(widget);

}

// Validation stuff	    //


// this should be a wrapper - and the real function be in graphics_info_t
//
#ifdef HAVE_GOOCANVAS
void
coot::set_validation_graph(int imol, coot::geometry_graph_type type, GtkWidget *dialog) {

   if (graphics_info_t::is_valid_model_molecule(imol)) {

      bool found = 0;
      if (type == coot::GEOMETRY_GRAPH_GEOMETRY) {
	 found = 1;
	 graphics_info_t::molecules[imol].validation_graphs.geometry_graph = dialog;
      }
      if (type == coot::GEOMETRY_GRAPH_B_FACTOR) {
	 found = 1;
	 graphics_info_t::molecules[imol].validation_graphs.b_factor_variance_graph = dialog;
      }
////B
      if (type == coot::GEOMETRY_GRAPH_CALC_B_FACTOR) {
	 found = 1;
	 graphics_info_t::molecules[imol].validation_graphs.b_factor_graph = dialog;
      }
////E
      if (type == coot::GEOMETRY_GRAPH_DENSITY_FIT) {
	 found = 1;
	 graphics_info_t::molecules[imol].validation_graphs.residue_density_fit_graph = dialog;
      }
      if (type == coot::GEOMETRY_GRAPH_OMEGA_DISTORTION) {
	 found = 1;
	 graphics_info_t::molecules[imol].validation_graphs.omega_distortion_graph = dialog;
      }
      if (type == coot::GEOMETRY_GRAPH_ROTAMER) {
	 found = 1;
	 graphics_info_t::molecules[imol].validation_graphs.rotamer_graph = dialog;
      }
      if (type == coot::GEOMETRY_GRAPH_NCS_DIFFS) {
	 found = 1;
	 graphics_info_t::molecules[imol].validation_graphs.ncs_diffs_graph = dialog;
      }
      if (type == coot::SEQUENCE_VIEW) {
	 found = 1;
	 graphics_info_t::molecules[imol].validation_graphs.sequence_view_is_displayed = dialog;
      }
      if (type == coot::RAMACHANDRAN_PLOT) {
	 found = 1;
	 graphics_info_t::molecules[imol].validation_graphs.dynarama_is_displayed = dialog;
      }


      if (!found) {
	 std::cout << "ERROR:: graph type " << type << " not found " << std::endl;
      }

   } else {
      std::cout << "WARNING:: set_validation_graph no valid molecule for imol = "
		<< imol << std::endl;
   }
}
#endif


#ifdef HAVE_GOOCANVAS
GtkWidget *
coot::get_validation_graph(int imol, coot::geometry_graph_type type) {

   GtkWidget *w = 0;
   if (graphics_info_t::is_valid_model_molecule(imol)) {
	switch(type) {
	case coot::GEOMETRY_GRAPH_GEOMETRY:
	   w = graphics_info_t::molecules[imol].validation_graphs.geometry_graph;
	   break;
	case coot::GEOMETRY_GRAPH_B_FACTOR:
	   w = graphics_info_t::molecules[imol].validation_graphs.b_factor_variance_graph;
	   break;
	case coot::GEOMETRY_GRAPH_CALC_B_FACTOR:
	   w = graphics_info_t::molecules[imol].validation_graphs.b_factor_graph;
	   break;
	case coot::GEOMETRY_GRAPH_DENSITY_FIT:
	   w = graphics_info_t::molecules[imol].validation_graphs.residue_density_fit_graph;
	   break;
	case coot::GEOMETRY_GRAPH_OMEGA_DISTORTION:
	   w = graphics_info_t::molecules[imol].validation_graphs.omega_distortion_graph;
	   break;
	case coot::GEOMETRY_GRAPH_ROTAMER:
	   w = graphics_info_t::molecules[imol].validation_graphs.rotamer_graph;
	   break;
	case coot::GEOMETRY_GRAPH_NCS_DIFFS:
	   w = graphics_info_t::molecules[imol].validation_graphs.ncs_diffs_graph;
	   break;
	case coot::SEQUENCE_VIEW:
	   w = graphics_info_t::molecules[imol].validation_graphs.sequence_view_is_displayed;
	   break;
	case coot::RAMACHANDRAN_PLOT:
	   w = graphics_info_t::molecules[imol].validation_graphs.dynarama_is_displayed; // terrible name for a widget
	   break;
	default:
	   break;
	}
   }
   return w;
}
#endif


// convenience function
void
graphics_info_t::update_geometry_graphs(int imol) {

   update_geometry_graphs(molecules[imol].atom_sel, imol);
}


// imol map is passed in case that the density fit graph was displayed.
// If there is no imol_map then the geometry graph could not have been displayed.  You can
// pass -1 for the map in that case.
void
graphics_info_t::update_geometry_graphs(mmdb::PResidue *SelResidues, int nSelResidues, int imol, int imol_map) { // searching for update_validation_graphs? Check the next function also

#ifdef HAVE_GOOCANVAS
   GtkWidget *graph = coot::get_validation_graph(imol, coot::GEOMETRY_GRAPH_ROTAMER);
   if (graph) {
      coot::geometry_graphs *gr = geometry_graph_dialog_to_object(graph);
      if (!gr) {
	 std::cout << "ERROR:: failed to get rotamer_graph from dialog\n";
      } else {
	 std::vector<coot::geometry_graph_block_info_generic> dv =
	    rotamers_from_residue_selection(SelResidues, nSelResidues, imol);
	 gr->update_residue_blocks(dv);
      }
   }
#endif

   update_validation(imol);
}

#include "nsv.hh"

// The molecule-based version of the above.
void
graphics_info_t::update_geometry_graphs(const atom_selection_container_t &moving_atoms_asc_local,  // searching for update_validation_graphs?
					int imol_moving_atoms) {

   update_validation(imol_moving_atoms);

#ifdef HAVE_GOOCANVAS
   GtkWidget *graph = coot::get_validation_graph(imol_moving_atoms, coot::GEOMETRY_GRAPH_GEOMETRY);
   if (graph) {
      // get deviations and replace those positions in the graph:
      coot::geometry_graphs *gr = geometry_graph_dialog_to_object(graph);
      if (!gr) {
	 std::cout << "ERROR:: failed to get geometry_graph from dialog\n";
      } else {
	 bool with_nbcs = false;
	 std::vector<coot::geometry_distortion_info_container_t> dv =
	    geometric_distortions_from_mol(imol_moving_atoms, moving_atoms_asc_local, with_nbcs);
	 for(unsigned int ich=0; ich<dv.size(); ich++)
// 	    std::cout << "       ich " << ich << " residue blocks for updating:\n"
// 		      << dv[ich] << std::endl;
	 for(unsigned int ich=0; ich<dv.size(); ich++)
	    gr->update_residue_blocks(dv[ich]);
      }
   }

   graph = coot::get_validation_graph(imol_moving_atoms, coot::GEOMETRY_GRAPH_DENSITY_FIT);

   if (graph) {
      coot::geometry_graphs *gr = geometry_graph_dialog_to_object(graph);
      if (!gr) {
	 std::cout << "ERROR:: failed to get residue_density_fit_graph from dialog\n";
      } else {
	 // Imol_Refinement_Map should be set by the time we get
	 // here.  There may be a pathological case where it has
	 // been closed when we get here.  density_fit_from_mol checks that.
	 std::vector<coot::geometry_graph_block_info_generic> dv =
	    density_fit_from_mol(moving_atoms_asc_local,
				 imol_moving_atoms,
				 Imol_Refinement_Map());
	 gr->update_residue_blocks(dv);
      }
   }

   graph = coot::get_validation_graph(imol_moving_atoms, coot::GEOMETRY_GRAPH_ROTAMER);
   if (graph) {
      coot::geometry_graphs *gr = geometry_graph_dialog_to_object(graph);
      if (!gr) {
	 std::cout << "ERROR:: failed to get rotamer_graph from dialog\n";
      } else {
	 std::vector<coot::geometry_graph_block_info_generic> dv =
	    rotamers_from_mol(moving_atoms_asc_local, imol_moving_atoms);

	 gr->update_residue_blocks(dv);
      }
   }

   graph = coot::get_validation_graph(imol_moving_atoms, coot::GEOMETRY_GRAPH_NCS_DIFFS);
   if (graph) {
      coot::geometry_graphs *gr = geometry_graph_dialog_to_object(graph);
      if (!gr) {
	 std::cout << "ERROR:: failed to get rotamer_graph from dialog\n";
      } else {
	 std::vector<coot::geometry_graph_block_info_generic> dv =
	    ncs_diffs_from_mol(imol_moving_atoms); // update everything
	 gr->update_residue_blocks(dv);
      }
   }

   graph = coot::get_validation_graph(imol_moving_atoms, coot::GEOMETRY_GRAPH_OMEGA_DISTORTION);
   if (graph) {
      coot::geometry_graphs *gr = geometry_graph_dialog_to_object(graph);
      if (!gr) {
	 std::cout << "ERROR:: failed to get omega_graph from dialog\n";
      } else {

	 if (! moving_atoms_asc_local.empty()) {

	    // We do this long handedly (c.f. above) because here we use
	    // render_omega_blocks() which needs the offset (which is a
	    // per-chain variable:
	    //
	    int n_models = moving_atoms_asc_local.mol->GetNumberOfModels();
	    for (int imodel = 1; imodel <= n_models; imodel++) {
	       mmdb::Model *model_p = moving_atoms_asc_local.mol->GetModel(imodel);
	       mmdb::Chain *chain_p;
	       const char *chain_id;
	       int n_chains = model_p->GetNumberOfChains();

	       for (int ich=0; ich<n_chains; ich++) {
		  chain_p = model_p->GetChain(ich);
		  chain_id = chain_p->GetChainID();
		  std::pair<short int, int> m = coot::util::min_resno_in_chain(chain_p);
		  if (m.first) {
		     // not used:
		     // int offset = m.second - 1; // min resno = 1 -> offset = 0

		     coot::omega_distortion_info_container_t om_dist =
			omega_distortions_from_mol(moving_atoms_asc_local, chain_id);

		     if (0)
			std::cout << "DEBUG:: update omega dist graph chain "
				  << om_dist.chain_id << " " << om_dist.omega_distortions.size()
				  << " blocks" << std::endl;

		     gr->update_omega_blocks(om_dist, ich, std::string(chain_id));
		  }
	       }
	    }
	 }
      }
   }

   graph = coot::get_validation_graph(imol_moving_atoms, coot::SEQUENCE_VIEW);
   if (graph) {

      exptl::nsv *sequence_view = static_cast<exptl::nsv *>(g_object_get_data(G_OBJECT(graph), "nsv"));

      if (sequence_view) {
	 mmdb::Manager *mol = molecules[imol_moving_atoms].atom_sel.mol;
	 sequence_view->regenerate(mol);
      }
   }

   // and now ramachandran also

   // 20211201-PE Hmm... this is blank - I wonder why...
   // Let's add an explicit function to update the rama plot (taking an imol)

   // update_ramachandran_plot(imol_moving_atoms);
#endif // HAVE_GOOCANVAS

}

void
graphics_info_t::update_ramachandran_plot(int imol) {

   std::cout << "update_ramachandran_plot() " << imol << " FIXME? " << std::endl;
}

#include "dynamic-validation.hh"

void
graphics_info_t::update_validation(int imol_changed_model) {

   // 20230910-PE, well, we only want to update the validation if it was already being displayed.
   // 20250202-PE so the caller sets the visibility of validation_boxes_vbox to true.

   if (! use_graphics_interface_flag) return;

   bool do_dynamic_validation = false;
   GtkWidget* vbox_dv = widget_from_builder("validation_boxes_vbox");
   if (gtk_widget_get_visible(vbox_dv)) do_dynamic_validation = true;

   update_validation_graphs(imol_changed_model);
   update_ramachandran_plot(imol_changed_model);
   if (do_dynamic_validation)
      update_dynamic_validation_for_molecule(imol_changed_model); // maybe.

   if (coot_all_atom_contact_dots_are_begin_displayed_for(imol_changed_model)) {
      mmdb::Manager *mol = molecules[imol_changed_model].atom_sel.mol;
      coot_all_atom_contact_dots_instanced(mol, imol_changed_model);
   }
}

void
graphics_info_t::update_validation_graphs(int imol_changed_model) {

   if (! use_graphics_interface_flag) return;

   // imol has change (e.g. a rotamer or RSR) and now I want to update the graphs
   // for that molecule if they are displayed.

   if (false) {
      g_debug("update_validation() called");
      g_warning("Reimplement update_validation(). "
                "The function should iterate over the std::map holding validation data for each active graph "
                "and recompute it, then trigger a redraw.");
   }

   // 20230527-PE maybe this is the right wasy to do it. But I think that I have done it a different way for now.
   //
   // update_ramachandran_plot(imol);

   if (active_validation_graph_model_idx == imol_changed_model) {
      for (const std::pair<const coot::validation_graph_type,GtkWidget*>& i : validation_graph_widgets) {
         coot::validation_graph_type graph_type = i.first;
         GtkWidget *graph = i.second;
         coot::validation_information_t vi = get_validation_data(imol_changed_model, graph_type);
         std::shared_ptr<coot::validation_information_t> vip = std::make_shared<coot::validation_information_t>(vi);
         CootValidationGraph *cvg = (CootValidationGraph *)(graph);
         coot_validation_graph_set_validation_information(cvg, vip);
      }
   }
}



void
graphics_info_t::delete_residue_from_geometry_graphs(int imol, coot::residue_spec_t res_spec) {

   update_validation(imol); // 20230528-PE we are not so clever (to be specific about what gets updated) now
}

void
graphics_info_t::delete_residues_from_geometry_graphs(int imol, const std::vector<coot::residue_spec_t> &res_specs) {

   update_validation(imol); // 20230528-PE again we are not so clever (to be specific about what gets updated) now

#ifdef HAVE_GOOCANVAS

   std::vector<coot::geometry_graph_type> graph_types;
   graph_types.push_back(coot::GEOMETRY_GRAPH_DENSITY_FIT);
   graph_types.push_back(coot::GEOMETRY_GRAPH_GEOMETRY);
   graph_types.push_back(coot::GEOMETRY_GRAPH_B_FACTOR);
   graph_types.push_back(coot::GEOMETRY_GRAPH_DENSITY_FIT);
   graph_types.push_back(coot::GEOMETRY_GRAPH_OMEGA_DISTORTION);
   graph_types.push_back(coot::GEOMETRY_GRAPH_ROTAMER);
   graph_types.push_back(coot::GEOMETRY_GRAPH_NCS_DIFFS);

   for (unsigned int igt=0; igt<graph_types.size(); igt++) {
      GtkWidget *graph =
	 coot::get_validation_graph(imol_moving_atoms, graph_types[igt]);
      if (graph) {
	 coot::geometry_graphs *gr = geometry_graph_dialog_to_object(graph);
	 if (gr) {

	    for (std::size_t ires=0; ires<res_specs.size(); ires++) {
	       const coot::residue_spec_t &res_spec = res_specs[ires];
	       gr->delete_block(res_spec.chain_id, res_spec.res_no);
	    }
	 }
      }
   }
#endif
}

void
graphics_info_t::delete_chain_from_geometry_graphs(int imol, const std::string &chain_id) {

#ifdef HAVE_GOOCANVAS

   std::vector<coot::geometry_graph_type> graph_types;
   graph_types.push_back(coot::GEOMETRY_GRAPH_DENSITY_FIT);
   graph_types.push_back(coot::GEOMETRY_GRAPH_GEOMETRY);
   graph_types.push_back(coot::GEOMETRY_GRAPH_B_FACTOR);
   graph_types.push_back(coot::GEOMETRY_GRAPH_DENSITY_FIT);
   graph_types.push_back(coot::GEOMETRY_GRAPH_OMEGA_DISTORTION);
   graph_types.push_back(coot::GEOMETRY_GRAPH_ROTAMER);
   graph_types.push_back(coot::GEOMETRY_GRAPH_NCS_DIFFS);

   for (unsigned int igt=0; igt<graph_types.size(); igt++) {
      GtkWidget *graph =
	 coot::get_validation_graph(imol_moving_atoms, graph_types[igt]);
      if (graph) {
	 coot::geometry_graphs *gr = geometry_graph_dialog_to_object(graph);
	 if (gr) {
	    // for res_spec in residues of chain that has bee deleted...
	    // gr->delete_block(res_spec.chain_id, res_spec.res_no);
	 }
      }
   }
#endif
}


void
graphics_info_t::geometric_distortion(int imol) {

#ifdef HAVE_GOOCANVAS

   // we need to assign these
//    int resno_1;
//    int resno_2;
   std::string chain_id_1;

   //    short int irest = 0;  // make 1 if restraints were found

   // make the selection and build a new molecule inside restraints.

   // short int have_flanking_residue_at_start = 0;
   // short int have_flanking_residue_at_end = 0;
   // short int have_disulfide_residues = 0;  // other residues are included in the
   // residues_mol for disphide
   // restraints.

   // 9 Sept 2003: The atom selection goes mad if residue with seqnum
   // iend_res+1 does not exist, but is not at the end of the chain.

   // Therefore we will set 2 flags, which tell us if istart_res-1 and
   // iend_res+1 exist.  And we do that by trying to select atoms from
   // them - if they exist, the number of selected atoms will be more
   // than 0.

//    istart_minus_flag = 0;  // from simple restraint code
//    iend_plus_flag    = 0;

   mmdb::Manager *mol = molecules[imol].atom_sel.mol; // short-hand usage

   // This dcv used to be inside the mol test, but that tickled what I
   // believe a compiler bug (on deconstructing this vector).  So now
   // it is here and coot doesn't crash when doing geometry analysis
   // of NMR model(s).  (In minimal testing).

   std::vector<coot::geometry_distortion_info_container_t> dcv;

   if (mol) {

      // big copy:
      bool with_nbcs = false;
      dcv = geometric_distortions_from_mol(imol, molecules[imol].atom_sel, with_nbcs);

      int max_chain_length = coot::util::max_min_max_residue_range(mol);
      if (max_chain_length <= 0) {
	 std::cout << "WARNING:: Funny coords - no graphs for this molecule" << std::endl;
      } else {
	 int nchains = coot::util::number_of_chains(mol);

	 std::string name = graphics_info_t::molecules[imol].name_for_display_manager();
	 coot::geometry_graphs *graphs = new coot::geometry_graphs(coot::GEOMETRY_GRAPH_GEOMETRY,
								   imol, name, nchains, max_chain_length);
	 coot::set_validation_graph(imol, coot::GEOMETRY_GRAPH_GEOMETRY, graphs->dialog()); // store for potential updates

	 // debugging.  Problem was negative occupancies.
// 	 for(unsigned int i=0; i<dcv.size(); i++) {
// 	    std::cout << i << " chain: " << dcv[i].chain_id << std::endl;
// 	    for(unsigned int j=0; j<dcv[i].geometry_distortion.size(); j++) {
// 	       std::cout << j << " " << dcv[i].geometry_distortion[j].distortion_score
// 			 << std::endl;
// 	    }
// 	 }

	 for(unsigned int i=0; i<dcv.size(); i++) {
	    graphs->render_to_canvas(dcv[i], i);
	 }
      }
   }
#endif
}

coot::geometry_distortion_info_container_t
graphics_info_t::geometric_distortions(int imol, mmdb::Residue *residue_p, bool with_nbcs) {

   coot::geometry_distortion_info_container_t gdc(NULL, 0, "");

   if (residue_p) {
      mmdb::Manager *mol = coot::util::create_mmdbmanager_from_residue(residue_p);
      if (mol) {
	 atom_selection_container_t asc = make_asc(mol);
	 std::vector<coot::geometry_distortion_info_container_t> v =
	    geometric_distortions_from_mol(imol, asc, with_nbcs);
	 if (v.size() == 1) {
	    if (v[0].geometry_distortion.size() > 1) {
	       gdc = v[0];
	    }
	 }
	 asc.clear_up();
      }
   }
   return gdc;
}

std::vector<coot::geometry_distortion_info_container_t>
graphics_info_t::geometric_distortions_from_mol(int imol, const atom_selection_container_t &asc,
						bool with_nbcs) {

   std::vector<coot::geometry_distortion_info_container_t> dcv;
   std::string altconf("");  // use this (e.g. "A") or "".

   if (! asc.mol)
      return dcv;

   int n_models = asc.mol->GetNumberOfModels();

   if (n_models > 0) {

      // 20100629, crash!  Ooops, we can't run over many models
      // because geometry_graphs (and its data member `blocks' are
      // sized to the number of chains).  If we run over all models,
      // then there are too many chains for the indexing of `blocks`
      // -> crash.  So we just use the first model.

      // for (int imod=1; imod<=n_models; imod++) {
      int imod=1;
      {

	 mmdb::Model *model_p = asc.mol->GetModel(imod);
	 mmdb::Chain *chain_p;
	 int nchains = model_p->GetNumberOfChains();
	 const char *chain_id;


	 for (int ichain=0; ichain<nchains; ichain++) {

	    chain_p = model_p->GetChain(ichain);

	    if (! chain_p->isSolventChain()) {
	       chain_id = chain_p->GetChainID();

	       // First make an atom selection of the residues selected to regularize.
	       //
	       int selHnd = asc.mol->NewSelection(); // yes, it's deleted.
	       int nSelResidues;
	       mmdb::PResidue *SelResidues = NULL;

	       // Consider as the altconf the altconf of one of the residues (we
	       // must test that the altlocs of the selected atoms to be either
	       // the same as each other (A = A) or one of them is "".  We need to
	       // know the mmdb syntax for "either".  Well, now I know that's ",A"
	       // (for either blank or "A").
	       //
	       //
	       //
	       asc.mol->Select(selHnd, mmdb::STYPE_RESIDUE, imod,
			       chain_id,
			       mmdb::ANY_RES, "*",
			       mmdb::ANY_RES, "*",
			       "*",  // residue name
			       "*",  // Residue must contain this atom name?
			       "*",  // Residue must contain this Element?
			       "*",  // altLocs
			       mmdb::SKEY_NEW // selection key
			       );
	       asc.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
	       std::pair<int, std::vector<std::string> > icheck =
		  check_dictionary_for_residue_restraints(imol, SelResidues, nSelResidues);

	       if (icheck.first == 0) {
		  for (unsigned int icheck_res=0; icheck_res<icheck.second.size(); icheck_res++) {
		     std::cout << "WARNING:: Failed to find restraints for "
			       << icheck.second[icheck_res] << std::endl;
		  }
	       }

	       std::cout << "INFO:: " << nSelResidues
			 << " residues selected for geometry checking object" << std::endl;

	       if (nSelResidues <= 0) {

		  std::cout << "ERROR:: No Residues!!   This should never happen:" << std::endl;
		  std::cout << "  in create_regularized_graphical_object" << std::endl;

	       } else { // normal

		  std::vector<mmdb::Atom *> fixed_atoms;
		  std::vector<coot::atom_spec_t> fixed_atom_specs;

		  // Notice that we have to make 2 atom selections, one, which includes
		  // flanking (and disulphide eventually) residues that is used for the
		  // restraints (restraints_container_t constructor) and one that is the
		  // moving atoms (which does not have flanking atoms).
		  //
		  // The restraints_container_t moves the atom of the mol that is passes to
		  // it.  This must be the same mol as the moving atoms mol so that the
		  // changed atom positions can be seen.  However (as I said) the moving
		  // atom mol should not have flanking residues shown.  So we make an asc
		  // that has the same mol as that passed to the restraints but a different
		  // atom selection (it is the atom selection that is used in the bond
		  // generation).
		  //

//	          20100210 try vector
// 		  coot::restraints_container_t restraints(SelResidues, nSelResidues,
// 							  std::string(chain_id),
// 							  asc.mol);
		  std::vector<std::pair<bool,mmdb::Residue *> > residue_vec;
		  for (int ires=0; ires<nSelResidues; ires++)
		     residue_vec.push_back(std::pair<bool, mmdb::Residue *> (0, SelResidues[ires]));

		  std::vector<mmdb::Link> links;
		  clipper::Xmap<float> dummy_xmap;

		  coot::restraints_container_t restraints(residue_vec,
							  links,
							  *Geom_p(),
							  asc.mol,
							  fixed_atom_specs,
							  &dummy_xmap);

		  // coot::restraint_usage_Flags flags = coot::BONDS;
		  // coot::restraint_usage_Flags flags = coot::BONDS_AND_ANGLES;
		  // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_AND_PLANES;
		  // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_TORSIONS_AND_PLANES;
		  coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_AND_NON_BONDED;
		  flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
		  flags = coot::BONDS_ANGLES_AND_PLANES;
		  flags = coot::BONDS_ANGLES_PLANES_AND_CHIRALS;

		  if (with_nbcs)
		     flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;

		  unsigned int n_threads = coot::get_max_number_of_threads();
		  if (n_threads > 0)
		     restraints.thread_pool(&static_thread_pool, n_threads);
		  short int do_residue_internal_torsions = 0;

		  // 	       if (do_torsion_restraints) {
		  // 		  do_residue_internal_torsions = 1;
		  // 		  flags = coot::BONDS_ANGLES_TORSIONS_PLANES_AND_NON_BONDED;
		  // 	       }

		  // 	       if (do_peptide_torsion_restraints)
		  // 		  do_link_torsions = 1;

		  coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
		  bool do_trans_peptide_restraints = false;
		  int nrestraints =
		     restraints.make_restraints(imol, *geom_p,
						flags,
						do_residue_internal_torsions,
						do_trans_peptide_restraints,
						0.0, 0, false, false, false,
						pseudos);

		  if (nrestraints > 0) {

// 		     std::cout << "DEBUG:: model " << imod << " pushing back " << nrestraints
// 			       << " restraints" << std::endl;

		     dcv.push_back(restraints.geometric_distortions());

		  } else {

		     // don't give this annoying dialog if restraints
		     // have been read for the residues in this chain.
		     // e.g. a single CLs residues in a chain.
		     std::vector<std::string> res_types = coot::util::residue_types_in_chain(chain_p);
		     bool hd = geom_p->have_restraints_dictionary_for_residue_types(res_types, imol,
                                                                                    cif_dictionary_read_number);
		     cif_dictionary_read_number += res_types.size();
		     if (! hd) {
			if (use_graphics_interface_flag) {
			   // GtkWidget *widget = create_no_restraints_info_dialog();
			   GtkWidget *widget = widget_from_builder("no_restraints_info_dialog");
			   gtk_widget_set_visible(widget, TRUE);
			} else {
			   std::cout << "WARNING:: No dictionary for some residue types " << std::endl;
			}
		     }
		  }
	       }
	       asc.mol->DeleteSelection(selHnd);
	    }
	 }
      }
   }
   // print_geometry_distortion(dcv);
   return dcv;
}

void
graphics_info_t::print_geometry_distortion(const std::vector<coot::geometry_distortion_info_container_t> &v) const {
   for (unsigned int i=0; i<v.size(); i++) {
      std::cout << v[i] << "\n";
   }
}


////B
void
graphics_info_t::calc_b_factor_graphs(int imol) {

#ifdef HAVE_GOOCANVAS

   if (imol<n_molecules()) {
      if (imol >= 0) {
	 if (molecules[imol].has_model()) {
	    mmdb::Manager *mol = molecules[imol].atom_sel.mol;
	    bool is_shelx_mol = molecules[imol].is_from_shelx_ins();

	    coot_extras::b_factor_analysis bfa(mol, is_shelx_mol);
	    std::vector<coot_extras::my_chain_of_stats_t> bfa_chain_info =
	       bfa.chain_details();

	    int n_models = mol->GetNumberOfModels();
	    int max_chain_length = coot::util::max_min_max_residue_range(mol);
	    if (max_chain_length <= 0) {
	       std::cout << "WARNING:: Funny coords - no graphs for this molecule"
			 << std::endl;
	    } else {
	       for (int imodel = 1; imodel <= n_models; imodel++) {
		  mmdb::Model *model_p = mol->GetModel(imodel);
		  mmdb::Chain *chain_p;
		  unsigned int n_chains = model_p->GetNumberOfChains();
		  coot::geometry_graphs *graphs =
		     new coot::geometry_graphs(coot::GEOMETRY_GRAPH_CALC_B_FACTOR,
					       imol,
					       graphics_info_t::molecules[imol].name_for_display_manager(),
					       n_chains, max_chain_length);
		  // b_factor_variance_graph[imol] = graphs->dialog();
		  set_validation_graph(imol, coot::GEOMETRY_GRAPH_CALC_B_FACTOR, graphs->dialog());

		  coot::b_factor_block_info_t bfi[3];
		  float std_dev;
		  int offset;

		  for (unsigned int ich=0; ich<n_chains; ich++) {

		     if (ich < bfa_chain_info.size()) {
			chain_p = model_p->GetChain(ich);
			std::pair<short int, int> m = coot::util::min_resno_in_chain(chain_p);

			if (m.first) {
			   std::vector<coot::b_factor_block_info_t> bfiv;
			   offset = m.second - 1;
			  for (unsigned int ires=0; ires<bfa_chain_info[ich].residue_properties.size(); ires++) {

				double variance[2]={0.0,0.0};
////	B FACTOR CALC FOR MAIN AND SIDECHAIN
				mmdb::PPAtom residue_atoms;
				mmdb::PResidue residue_p = chain_p->GetResidue(ires);
   				int nResidueAtoms=0;
				double running_sum[4]	= {0.0,0.0,0.0,0.0};
 				double mean[3]		= {0.0,0.0,0.0};
				double std_dev[2]	= {0.0,0.0};
				double bf  	= 0.0;
				double occ 	= 0.0;
				double bfo 	= 0.0;
 				double div[2] 	= {0.0,0.0};
   				int    bMC 	= 0;

				residue_p->GetAtomTable(residue_atoms, nResidueAtoms);
				if (nResidueAtoms > 0) {
					for (int i=0; i<nResidueAtoms; i++) {
						bMC = ( coot::is_main_chain_p(residue_atoms[i]) )?1:0;
						std::string ele = residue_atoms[i]->element;
						if ((ele != " H") && (ele != " D")) {
							bf  = residue_atoms[i]->tempFactor;
							occ = residue_atoms[i]->occupancy;
							if ( ((bf > 0.0) && (occ >= 0.0) && (occ <= 1.0)) ||
							(is_shelx_mol && (occ < 11.001) && (occ > 10.999))) {
								if (is_shelx_mol)
									occ = 1.0;
								div[bMC] 		+= occ;
								bfo			 = bf*occ;
								running_sum[0+2*bMC] 	+= bfo;
								running_sum[1+2*bMC] 	+= bfo*bfo;
							}
						}
					}
					if ( div[0] > 0 || div[1] > 0 ) {
						mean[0]     = (div[0]>0)?(running_sum[0]/div[0]):(0.0);		// notMC
						mean[1]     = (running_sum[0]+running_sum[2])/(div[0]+div[1]);
					}
				}else{ //SO THIS SHOULD NEVER HAPPEN WILL JUST YIELD ZERO PLOTS
					std::cout << "ERROR::  IN B FACTOR CALCULATION, EMPTY RESIDUE" << std::endl;
				}
////	END AVERAGE B FACTOR
////	SIDECHAIN <B>_sc
			      bfi[0].resno = bfa_chain_info[ich].residue_properties[ires].resno;
			      bfi[0].b_factor_var = mean[0]+mean[1]; // IN ORDER TO PLOT SC ON THE TOP
			      bfi[0].info_string  = int_to_string(bfi[0].resno);
			      bfi[0].info_string += chain_p->GetChainID();
			      bfi[0].info_string += " ";
			      bfi[0].info_string += bfa_chain_info[ich].residue_properties[ires].resname;
			      bfi[0].info_string += ":SC: ";
			      bfi[0].info_string += float_to_string(bfi[0].b_factor_var - mean[1]); // USER SEES AVERAGE SC B FACTOR
			      bfi[0].atom_name = bfa_chain_info[ich].residue_properties[ires].atom_name;
////	TOTAL <B>_tot
		 	      bfi[1].resno = bfa_chain_info[ich].residue_properties[ires].resno;
			      bfi[1].b_factor_var = mean[1];
			      bfi[1].info_string  = int_to_string(bfi[1].resno);
			      bfi[1].info_string += chain_p->GetChainID();
			      bfi[1].info_string += " ";
			      bfi[1].info_string += bfa_chain_info[ich].residue_properties[ires].resname;
			      bfi[1].info_string += ":TOT: ";
			      bfi[1].info_string += float_to_string(bfi[1].b_factor_var); // USER SEES THE ACTUAL VALUE SINCE IT IS IN THE FRONT
			      bfi[1].atom_name = bfa_chain_info[ich].residue_properties[ires].atom_name;

			      bfiv.push_back(bfi[0]); //SC
			      bfiv.push_back(bfi[1]); //TOT
			 }
			 graphs->render_b_factor_blocks(imol, ich, bfa_chain_info[ich].chain_id,
							  offset, bfiv);
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
#endif
}
////E

void
graphics_info_t::b_factor_graphs(int imol) {

#ifdef HAVE_GOOCANVAS

   if (imol<n_molecules()) {
      if (imol >= 0) {
	 if (molecules[imol].has_model()) {
	    mmdb::Manager *mol = molecules[imol].atom_sel.mol;
	    bool is_shelx_mol = molecules[imol].is_from_shelx_ins();

	    coot_extras::b_factor_analysis bfa(mol, is_shelx_mol);
	    std::vector<coot_extras::my_chain_of_stats_t> bfa_chain_info =
	       bfa.chain_details();

	    int n_models = mol->GetNumberOfModels();
	    int max_chain_length = coot::util::max_min_max_residue_range(mol);
	    if (max_chain_length <= 0) {
	       std::cout << "WARNING:: Funny coords - no graphs for this molecule"
			 << std::endl;
	    } else {
	       for (int imodel = 1; imodel <= n_models; imodel++) {
		  mmdb::Model *model_p = mol->GetModel(imodel);
		  mmdb::Chain *chain_p;
		  unsigned int n_chains = model_p->GetNumberOfChains();
		  coot::geometry_graphs *graphs =
		     new coot::geometry_graphs(coot::GEOMETRY_GRAPH_B_FACTOR,
					       imol,
					       graphics_info_t::molecules[imol].name_for_display_manager(),
					       n_chains, max_chain_length);
		  // b_factor_variance_graph[imol] = graphs->dialog();
		  set_validation_graph(imol, coot::GEOMETRY_GRAPH_B_FACTOR, graphs->dialog());

		  coot::b_factor_block_info_t bfi;
		  float std_dev;
		  int offset;

		  for (unsigned int ich=0; ich<n_chains; ich++) {

		     if (ich < bfa_chain_info.size()) {
			chain_p = model_p->GetChain(ich);
			std::pair<short int, int> m = coot::util::min_resno_in_chain(chain_p);

			if (m.first) {
			   std::vector<coot::b_factor_block_info_t> bfiv;
			   offset = m.second - 1;

			   for (unsigned int ires=0; ires<bfa_chain_info[ich].residue_properties.size(); ires++) {
			      bfi.resno = bfa_chain_info[ich].residue_properties[ires].resno;
			      std_dev = bfa_chain_info[ich].residue_properties[ires].std_dev;
			      bfi.b_factor_var = std_dev * std_dev;
			      bfi.info_string  = int_to_string(bfi.resno);
			      bfi.info_string += chain_p->GetChainID();
			      bfi.info_string += " ";
			      bfi.info_string += bfa_chain_info[ich].residue_properties[ires].resname;
			      bfi.info_string += ": ";
			      bfi.info_string += float_to_string(bfi.b_factor_var);
			      bfi.atom_name = bfa_chain_info[ich].residue_properties[ires].atom_name;
			      bfiv.push_back(bfi);
			   }
			   graphs->render_b_factor_blocks(imol, ich, bfa_chain_info[ich].chain_id,
							  offset, bfiv);
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
#endif
}

void
graphics_info_t::omega_graphs(int imol) {

#ifdef HAVE_GOOCANVAS

   if (imol >= 0) {
      if (imol < n_molecules()) {
	 if (molecules[imol].has_model()) {
	    mmdb::Manager *mol = molecules[imol].atom_sel.mol;
	    int n_models = mol->GetNumberOfModels();
	    int max_chain_length = coot::util::max_min_max_residue_range(mol);
	    if (max_chain_length <= 0) {
	       std::cout << "WARNING:: Funny coords - no graphs for this molecule" << std::endl;
	    } else {
	       for (int imodel = 1; imodel <= n_models; imodel++) {
		  mmdb::Model *model_p = mol->GetModel(imodel);
		  mmdb::Chain *chain_p;
		  const char *chain_id;
		  int n_chains = model_p->GetNumberOfChains();
		  coot::geometry_graphs *graphs =
		     new coot::geometry_graphs(coot::GEOMETRY_GRAPH_OMEGA_DISTORTION,
					       imol,
					       graphics_info_t::molecules[imol].name_for_display_manager(),
					       n_chains, max_chain_length);
		  // omega_distortion_graph[imol] = graphs->dialog();
		  coot::set_validation_graph(imol, coot::GEOMETRY_GRAPH_OMEGA_DISTORTION, graphs->dialog());

		  for (int ich=0; ich<n_chains; ich++) {
		     chain_p = model_p->GetChain(ich);
		     chain_id = chain_p->GetChainID();
		     std::pair<short int, int> m = coot::util::min_resno_in_chain(chain_p);
		     if (m.first) {
			int offset = m.second - 1; // min resno = 1 -> offset = 0
			int selHnd = mol->NewSelection();
			mmdb::PResidue *SelResidues = NULL;
			int nSelResidues;

			mol->Select(selHnd, mmdb::STYPE_RESIDUE, 0,
				    chain_id,
				    mmdb::ANY_RES, "*",
				    mmdb::ANY_RES, "*",
				    "*",  // residue name
				    "*",  // Residue must contain this atom name?
				    "*",  // Residue must contain this Element?
				    "*",  // altLocs
				    mmdb::SKEY_NEW // selection key
				    );
			mol->GetSelIndex(selHnd, SelResidues, nSelResidues);
			clipper::Xmap<float> dummy_xmap;

			if (nSelResidues > 0) {
			   coot::restraints_container_t restraints(molecules[imol].atom_sel,
								   std::string(chain_id),
								   &dummy_xmap);

			   coot::omega_distortion_info_container_t om_dist =
			      restraints.omega_trans_distortions(*geom_p,
								 mark_cis_peptides_as_bad_flag);
			   // std::cout << "DEBUG: got om_dist." << std::endl;

			   graphs->render_omega_blocks(om_dist, ich, std::string(chain_id),
						       offset);
			}
			mol->DeleteSelection(selHnd);
		     }
		  }
	       }
	    }
	 }
      }
   }
#endif // HAVE_GOOCANVAS
}

coot::omega_distortion_info_container_t
graphics_info_t::omega_distortions_from_mol(const atom_selection_container_t &asc,
					    const std::string &chain_id) {

   clipper::Xmap<float> dummy_xmap;
   coot::restraints_container_t restraints(asc, chain_id, &dummy_xmap);
   coot::omega_distortion_info_container_t om_dist =
      restraints.omega_trans_distortions(*geom_p, mark_cis_peptides_as_bad_flag);
   return om_dist;
}

// 20240420-PE this return value can be exported to python rather than displayed internally
coot::rotamer_graphs_info_t
graphics_info_t::rotamer_graphs(int imol) {

   coot::rotamer_graphs_info_t info;

   if (is_valid_model_molecule(imol)) {
      mmdb::Manager *mol = molecules[imol].atom_sel.mol;

      for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
         mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               int n_res = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<n_res; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     std::string alt_conf = ""; // fix this lazyness one day
                     coot::rotamer_probability_info_t d_score =
                        get_rotamer_probability(residue_p, alt_conf, mol, rotamer_lowest_probability, 1);
                     int this_resno = residue_p->GetSeqNum();
                     std::string this_inscode = residue_p->GetInsCode();
                     std::string chain_id = residue_p->GetChainID();
                     std::string str = std::to_string(this_resno);
                     str += chain_id;
                     str += " ";
                     str += residue_p->name;
                     str += " ";
                     double distortion = 0.0;
                     switch (d_score.state) {

                     case 1: {
                        // On reflection we don't want found
                        // rotamers with low probabilities to be
                        // marked (nearly as or more) bad than not
                        // found rotamers:
                        distortion = 100.0 * rotamer_distortion_scale/d_score.probability;
                        distortion = distortion > 100.0 ? 100.0 : distortion;
                        str += " ";
                        str += d_score.rotamer_name;
                        str += " ";
                        str += "Probability: ";
                        str += float_to_string(d_score.probability);
                        str += "%";
                        coot::graph_rotamer_info_t ri(chain_id, this_resno, this_inscode, d_score.probability, d_score.rotamer_name);
                        info.info.push_back(ri);
                        break;
                     }

                     case 0: {
                        distortion = 100.0;
                        str += "Missing Atoms";
                        coot::graph_rotamer_info_t ri(chain_id, this_resno, this_inscode, 0.0, "Missing Atoms");
                        info.info.push_back(ri);
                        break;
                     }

                     case -1: {
                        distortion = 100.0;
                        coot::graph_rotamer_info_t ri(chain_id, this_resno, this_inscode, 0.0, "Rotamer not recognised");
                        info.info.push_back(ri);
                        str += "Rotamer not recognised";
                        break;
                     }

                     case -2: {  // don't plot a block for this one.
                        distortion = 0.0;
                        break;
                     }
                     }
                  }
               }
            }
         }
      }
   }
   return info;
}


#ifdef HAVE_GOOCANVAS

coot::rotamer_graphs_info_t
graphics_info_t::old_rotamer_graphs(int imol) {

   coot::rotamer_graphs_info_t info; // return value

   if (imol >= 0) {
      if (imol < n_molecules()) {
	 if (molecules[imol].has_model()) {
	    mmdb::Manager *mol = molecules[imol].atom_sel.mol;
	    int n_models = mol->GetNumberOfModels();
	    int max_chain_length = coot::util::max_min_max_residue_range(mol);
	    if (max_chain_length <= 0) {
	       std::cout << "WARNING:: Funny coords - no graphs for this molecule" << std::endl;
	    } else {
	       for (int imodel = 1; imodel <= n_models; imodel++) {
		  mmdb::Model *model_p = mol->GetModel(imodel);
		  mmdb::Chain *chain_p;
		  const char *chain_id;
		  int n_chains = model_p->GetNumberOfChains();
		  coot::geometry_graphs *graphs = 0;
		  if (use_graphics_interface_flag) {
		     std::string mol_name = graphics_info_t::molecules[imol].name_for_display_manager();
		     graphs = new coot::geometry_graphs(coot::GEOMETRY_GRAPH_ROTAMER,
							imol, mol_name,
							n_chains, max_chain_length);

		     // rotamer_graph[imol] = graphs->dialog();
		     coot::set_validation_graph(imol, coot::GEOMETRY_GRAPH_ROTAMER, graphs->dialog());
		  }
		  for (int ich=0; ich<n_chains; ich++) {
		     chain_p = model_p->GetChain(ich);
		     if (! chain_p->isSolventChain()) {
			chain_id = chain_p->GetChainID();
			std::pair<short int, int> m = coot::util::min_resno_in_chain(chain_p);
			if (m.first) {
			   int offset = m.second - 1;
			   int selHnd = mol->NewSelection();
			   mmdb::PResidue *SelResidues = NULL;
			   int nSelResidues;

			   mol->Select(selHnd, mmdb::STYPE_RESIDUE, 0,
				       chain_id,
				       mmdb::ANY_RES, "*",
				       mmdb::ANY_RES, "*",
				       "*",  // residue name
				       "*",  // Residue must contain this atom name?
				       "*",  // Residue must contain this Element?
				       "*",  // altLocs
				       mmdb::SKEY_NEW // selection key
				       );
			   mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

			   if (nSelResidues > 0) {
			      int max_resno = -9999;

			      std::string altconf = ""; // fixme, I guess
			      std::vector<coot::geometry_graph_block_info_generic> v;
			      for (int ir=0; ir<nSelResidues; ir++) {
				 int this_resno = SelResidues[ir]->GetSeqNum();
				 std::string res_name = SelResidues[ir]->GetResName();
				 if (coot::util::is_standard_amino_acid_name(res_name)) {
				 // if (res_name != "HOH") {
				    std::string this_inscode = SelResidues[ir]->GetInsCode();
				    if (this_resno > max_resno)
				       max_resno = this_resno;
				    coot::rotamer_probability_info_t d_score =
				       get_rotamer_probability(SelResidues[ir], altconf, mol,
							       rotamer_lowest_probability, 1);

				    double distortion = 0.0;
				    std::string str = int_to_string(this_resno);
				    str += chain_id;
				    str += " ";
				    str += SelResidues[ir]->name;
				    str += " ";
				    coot::atom_spec_t atom_spec(chain_id, this_resno,
								"", " CA ", "");
				    coot::graph_rotamer_info_t ri("", -9999, "", -99.9, "unset");
				    switch (d_score.state) {

				    case 1:
				       // On reflection we don't want found
				       // rotamers with low probabilities to be
				       // marked (nearly as or more) bad than not
				       // found rotamers:
				       distortion = 100.0 * rotamer_distortion_scale/d_score.probability;
				       distortion = distortion > 100.0 ? 100.0 : distortion;
				       str += " ";
				       str += d_score.rotamer_name;
				       str += " ";
				       str += "Probability: ";
				       str += float_to_string(d_score.probability);
				       str += "%";
				       ri = coot::graph_rotamer_info_t (chain_id, this_resno,
									this_inscode,
									d_score.probability,
									d_score.rotamer_name);
				       info.info.push_back(ri);
				       break;

				    case 0:
				       distortion = 100.0;
				       str += "Missing Atoms";
				       atom_spec.string_user_data = "Missing Atoms";
				       ri = coot::graph_rotamer_info_t (chain_id, this_resno,
									this_inscode,
									0.0, "Missing Atoms");
				       info.info.push_back(ri);
				       break;

				    case -1:
				       distortion = 100.0;
				       ri = coot::graph_rotamer_info_t (chain_id, this_resno,
									this_inscode,
									0.0, "Rotamer not recognised");
				       info.info.push_back(ri);
				       str += "Rotamer not recognised";
				       break;

				    case -2:   // don't plot a block for this one.
				       distortion = 0.0;
				       break;
				    }
				    if (d_score.state != -2) {
				       v.push_back(coot::geometry_graph_block_info_generic(imol, this_resno, atom_spec, distortion, str));
				    }
				 }
			      } // end residue for loop
			      // done residue loop:
// 			      std::cout << "render_to_canvas: (rotamer) chain: " << ich << " min_resno: "
// 					<< m.second << " max_resno: " << max_resno
// 					<< " offset: " << offset << std::endl;
			      if (use_graphics_interface_flag && max_resno > 0)
				 graphs->render_to_canvas(v, ich, std::string(chain_id),
							  max_resno, m.second, offset);

			   }

			   mol->DeleteSelection(selHnd);
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   return info;
}
#endif // HAVE_GOOCANVAS

#ifdef HAVE_GOOCANVAS
std::vector<coot::geometry_graph_block_info_generic>
graphics_info_t::rotamers_from_mol(const atom_selection_container_t &asc,
				  int imol_moving_atoms) {

   // this does not use the provided atom_selection_container_t asc

   std::vector<coot::geometry_graph_block_info_generic> dv;

   mmdb::Manager *mol = molecules[imol_moving_atoms].atom_sel.mol;
   // int n_models = mol->GetNumberOfModels();
   // for (int imodel = 1; imodel <= n_models; imodel++) {
   int imodel = 1;
   mmdb::Model *model_p = mol->GetModel(imodel);
   mmdb::Chain *chain_p;
   const char *chain_id;
   int n_chains = model_p->GetNumberOfChains();

   for (int ich=0; ich<n_chains; ich++) {
      chain_p = model_p->GetChain(ich);
      if (! chain_p->isSolventChain()) {
	 chain_id = chain_p->GetChainID();
	 std::pair<short int, int> m = coot::util::min_resno_in_chain(chain_p);
	 if (m.first) {
	    // int offset = m.second - 1;
	    int selHnd = mol->NewSelection();
	    mmdb::PResidue *SelResidues = NULL;
	    int nSelResidues;

	    mol->Select(selHnd, mmdb::STYPE_RESIDUE, 0,
			chain_id,
			mmdb::ANY_RES, "*",
			mmdb::ANY_RES, "*",
			"*",  // residue name
			"*",  // Residue must contain this atom name?
			"*",  // Residue must contain this Element?
			"*",  // altLocs
			mmdb::SKEY_NEW // selection key
			);
	    mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

	    if (nSelResidues > 0) {
	       int max_resno = -9999;

	       for (int ir=0; ir<nSelResidues; ir++) {
		  int this_resno = SelResidues[ir]->GetSeqNum();
		  std::string res_name = SelResidues[ir]->GetResName();
		  if (res_name != "HOH") {
		     if (this_resno > max_resno)
			max_resno = this_resno;
		     std::string altconf = "";
		     coot::rotamer_probability_info_t d_score =
			get_rotamer_probability(SelResidues[ir], altconf, mol,
						rotamer_lowest_probability, 1);
		     double distortion = 0.0;
		     std::string str = int_to_string(this_resno);
		     str += chain_id;
		     str += " ";
		     str += SelResidues[ir]->name;
		     str += " ";
		     coot::atom_spec_t atom_spec(chain_id, this_resno,
						 "", " CA ", "");
		     switch (d_score.state) {

		     case 1:
			// On reflection we don't want found
			// rotamers with low probabilities to be
			// marked (nearly as or more) bad than not
			// found rotamers:
			distortion = 100.0 * rotamer_distortion_scale/d_score.probability;
			distortion = distortion > 100.0 ? 100.0 : distortion;
			str += " ";
			str += d_score.rotamer_name;
			str += " ";
			str += "Probability: ";
			str += float_to_string(d_score.probability);
			str += "%";
			break;

		     case 0:
			distortion = 100.0;
			str += "Missing Atoms";
			atom_spec.string_user_data = "Missing Atoms";
			break;

		     case -1:
			distortion = 100.0;
			str += "Rotamer not recognised";
			break;

		     case -2:   // don't plot a block for this one.
			distortion = 0.0;
			break;
		     }
		     if (d_score.state != -2) {
			dv.push_back(coot::geometry_graph_block_info_generic(imol_moving_atoms, this_resno, atom_spec, distortion, str));
		     }
		  }
	       }
	    }
	    mol->DeleteSelection(selHnd);
	 }
      }
   }
   return dv;
}
#endif

#ifdef HAVE_GOOCANVAS
std::vector<coot::geometry_graph_block_info_generic>
graphics_info_t::rotamers_from_residue_selection(mmdb::PResidue *SelResidues,
						 int nSelResidues, int imol) {

   std::vector<coot::geometry_graph_block_info_generic> v;
   for (int ires=0; ires<nSelResidues; ires++) {
      std::string res_name = SelResidues[ires]->GetResName();
      if (res_name != "HOH") {
	 int this_resno = SelResidues[ires]->GetSeqNum();
	 std::string chain_id = SelResidues[ires]->GetChainID();
	 std::string alt_conf = ""; // fixme?
	 coot::rotamer_probability_info_t d_score =
	    get_rotamer_probability(SelResidues[ires], alt_conf, 0, // 0? woo..
				    rotamer_lowest_probability, 1);
	 double distortion = 0.0;
	 std::string str = int_to_string(this_resno);
	 str += chain_id;
	 str += " ";
	 str += SelResidues[ires]->name;
	 str += " ";
	 switch (d_score.state) {

	 case 1:
	    // On reflection we don't want found
	    // rotamers with low probabilities to be
	    // marked (nearly as or more) bad than not
	    // found rotamers:
	    distortion = 100.0 * rotamer_distortion_scale/d_score.probability;
	    distortion = distortion > 100.0 ? 100.0 : distortion;
	    str += " ";
	    str += d_score.rotamer_name;
	    str += " ";
	    str += "Probability: ";
	    str += float_to_string(d_score.probability);
	    str += "%";
	    break;

	 case 0:
	    distortion = 100.0;
	    str += "Missing Atoms";
	    break;

	 case -1:
	    distortion = 100.0;
	    str += "Rotamer not recognised";
	    break;

	 case -2:   // don't plot a block for this one.
	    distortion = 0.0;
	    break;
	 }
	 if (d_score.state != -2) {
	    coot::atom_spec_t atom_spec(chain_id, this_resno,
					"", " CA ", "");
	    v.push_back(coot::geometry_graph_block_info_generic(imol_moving_atoms,
								this_resno,
								atom_spec,
								distortion,
								str));
	 }
      }
   }
   return v;
}
#endif

void
graphics_info_t::density_fit_graphs(int imol) {

#ifdef HAVE_GOOCANVAS

   if (imol >= 0) {
      if (imol < n_molecules()) {
	 if (molecules[imol].has_model()) {

	    int imol_for_map = Imol_Refinement_Map();
	    if (imol_for_map == -1)
	       show_select_map_dialog();
            // maybe we have it now?!
            imol_for_map = Imol_Refinement_Map();
	    if (imol_for_map > -1) {
	       // std::cout << "DEBUG:: starting" << std::endl;
	       mmdb::Manager *mol = molecules[imol].atom_sel.mol;
	       // double max_grid_factor = coot::util::max_gridding(molecules[imol_for_map].xmap_list[0]);
	       // std::cout << "DEBUG:: max_grid_factor: " << max_grid_factor << std::endl;
	       // double sq_max_grid_fac = max_grid_factor * max_grid_factor;
	       // 1.8A data -> max_grid_factor = 0.6
	       // 3.0A data -> max_grid_factor = 1.0
	       // squared max_grid_factor: 0.36: 4.0 is good
	       // squared max_grid_factor: 1.0 : 4.0 is 4 times too much
	       // so we want 1/squared(max_grid_factor) not 4.0
	       int n_models = mol->GetNumberOfModels();
	       int max_chain_length = coot::util::max_min_max_residue_range(mol);
	       if (max_chain_length <= 0) {
		  std::cout << "WARNING:: Funny coords - no graphs\n";
	       } else {
		  for (int imodel = 1; imodel <= n_models; imodel++) {
		     mmdb::Model *model_p = mol->GetModel(imodel);
                     if (! model_p) continue;
		     int n_chains = model_p->GetNumberOfChains();
                     std::string mol_name = graphics_info_t::molecules[imol].name_for_display_manager();
		     coot::geometry_graphs *graphs =
			new coot::geometry_graphs(coot::GEOMETRY_GRAPH_DENSITY_FIT,
						  imol, mol_name,
						  n_chains, max_chain_length);

		     // residue_density_fit_graph[imol] = graphs->dialog();
		     set_validation_graph(imol, coot::GEOMETRY_GRAPH_DENSITY_FIT, graphs->dialog());

		     for (int ich=0; ich<n_chains; ich++) {
                        mmdb::Chain *chain_p = model_p->GetChain(ich);
			const char *chain_id = chain_p->GetChainID();
			std::pair<short int, int> m = coot::util::min_resno_in_chain(chain_p);
			if (m.first) {
			   int offset = m.second - 1;
			   int selHnd = mol->NewSelection(); // d
			   mmdb::PResidue *SelResidues = NULL;
			   int nSelResidues;

			   mol->Select(selHnd, mmdb::STYPE_RESIDUE, 0,
				       chain_id,
				       mmdb::ANY_RES, "*",
				       mmdb::ANY_RES, "*",
				       "*",  // residue name
				       "*",  // Residue must contain this atom name?
				       "*",  // Residue must contain this Element?
				       "*",  // altLocs
				       mmdb::SKEY_NEW // selection key
				       );
			   mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

			   std::vector<coot::geometry_graph_block_info_generic> v =
			      graphics_info_t::density_fit_from_residues(SelResidues, nSelResidues,
									 imol,
									 imol_for_map);

			   if (nSelResidues > 0) {
			      int max_resno = -9999;
			      for (int ires=0; ires<nSelResidues; ires++)
				 if (SelResidues[ires]->GetSeqNum() > max_resno)
				    max_resno = SelResidues[ires]->GetSeqNum();

			      graphs->render_to_canvas(v, ich, std::string(chain_id),
						       max_resno, m.second, offset);
			   }
			   mol->DeleteSelection(selHnd);
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
#endif // HAVE_GOOCANVAS
}


// Use this to update the molecule only.  (Because we call a graph
// function that doesn't setup the canvas and draw the ticks (it only
// updates the blocks that it has been given).
//
// We pass imol_moving_atoms because we will be updating a graph and
// we want to know which graph to update.
//
#ifdef HAVE_GOOCANVAS
std::vector<coot::geometry_graph_block_info_generic>
graphics_info_t::density_fit_from_mol(const atom_selection_container_t &asc,
				      int imol_moving_atoms,
				      int imol_map) {


   std::vector<coot::geometry_graph_block_info_generic> drv;
   std::string altconf("");  // use this (e.g. "A") or "".

   if (! asc.mol)
      return drv;

   if (imol_map < n_molecules() && graphics_info_t::molecules[imol_map].has_xmap()) {
      int n_models = asc.mol->GetNumberOfModels();

      if (n_models > 0) {

	 for (int imod=1; imod<=n_models; imod++) {

	    mmdb::Model *model_p = asc.mol->GetModel(imod);
	    mmdb::Chain *chain_p;
	    int nchains = model_p->GetNumberOfChains();

	    for (int ichain=0; ichain<nchains; ichain++) {

	       chain_p = model_p->GetChain(ichain);
               const char *chain_id = chain_p->GetChainID();

	       // Maybe we could do a chain->GetResidueTable() here
	       // instead of a selection.

	       int selHnd = asc.mol->NewSelection();
	       int nSelResidues;
	       mmdb::PResidue *SelResidues = NULL;

	       asc.mol->Select(selHnd, mmdb::STYPE_RESIDUE, 0,
			       chain_id,
			       mmdb::ANY_RES, "*",
			       mmdb::ANY_RES, "*",
			       "*",  // residue name
			       "*",  // Residue must contain this atom name?
			       "*",  // Residue must contain this Element?
			       "*",  // altLocs
			       mmdb::SKEY_NEW // selection key
			       );
	       asc.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

	       std::vector<coot::geometry_graph_block_info_generic> v =
		  graphics_info_t::density_fit_from_residues(SelResidues, nSelResidues,
							     imol_moving_atoms, imol_map);

	       for (unsigned int i=0; i<v.size(); i++)
		  drv.push_back(v[i]);

	       // the graph update is done in the function that calls this one.

	       asc.mol->DeleteSelection(selHnd);

	    }
	 }
      }
   }
   return drv;
}
#endif


// To be called for each chain in the molecule (or atom selection).
//
#ifdef HAVE_GOOCANVAS
std::vector<coot::geometry_graph_block_info_generic>
graphics_info_t::density_fit_from_residues(mmdb::PResidue *SelResidues, int nSelResidues,
					   int imol,
					   int imol_for_map) const {

   std::vector<coot::geometry_graph_block_info_generic> v;
   if (nSelResidues > 0) {
      int max_resno = -9999;
      double max_grid_factor = coot::util::max_gridding(molecules[imol_for_map].xmap);

      for (int ir=0; ir<nSelResidues; ir++) {
         mmdb::Residue *residue_p = SelResidues[ir];
	 int this_resno = residue_p->GetSeqNum();
	 if (this_resno > max_resno)
	    max_resno = this_resno;

	 mmdb::PAtom *residue_atoms=0;
	 int n_residue_atoms;

	 SelResidues[ir]->GetAtomTable(residue_atoms, n_residue_atoms);

	 double residue_density_score =
	    coot::util::map_score(residue_atoms,
				  n_residue_atoms,
				  molecules[imol_for_map].xmap, 1);
	 double occ_sum = coot::util::occupancy_sum(residue_atoms, n_residue_atoms);
	 if (occ_sum > 0) {
            float distortion_max_abs = 132.0;
            distortion_max_abs = 140; // 20220115-PE let's say
            float distortion_max = distortion_max_abs;
	    residue_density_score /= occ_sum;
	    std::string str = int_to_string(this_resno);
	    str += residue_p->GetChainID();
	    str += " ";
	    str += residue_p->name;
	    str += " ";
	    str += float_to_string(residue_density_score);

	    if (residue_density_score < 0.0001)
	       residue_density_score = 0.0001;

	    // std::cout << "DEBUG::          max_grid_factor " << max_grid_factor
	    // << " score " << residue_density_score << std::endl;
	    double sf = residue_density_fit_scale_factor * 0.72; // 20220115-PE  was 1.25;
	    // high resolution maps have high grid factors (say 0.5) and high
	    // residue_density_ scores (say 2.0)
	    double distortion =  sf/(pow(max_grid_factor,3) * residue_density_score);
	    distortion =  sf/(pow(max_grid_factor,4) * residue_density_score); // seems reasonable!

	    // distortion *= distortion; // non-linear, provides distinction.

	    if (distortion > distortion_max)
	       distortion = distortion_max;
	    // use intelligent atom name here
	    std::string chain_id = residue_p->GetChainID();
            std::string atom_name = coot::util::intelligent_this_residue_mmdb_atom(residue_p)->GetAtomName();
            coot::atom_spec_t atom_spec(chain_id, this_resno, "", atom_name, "");
            // std::cout << "creating block with distortion " << distortion << std::endl;
	    v.push_back(coot::geometry_graph_block_info_generic(imol, this_resno, atom_spec, distortion, str));
	 }
      }
      //      graphs->render_to_canvas(v, ich, std::string(chain_id),
      //                               max_resno, m.second, offset);
   }
   return v;
}
#endif

#ifdef HAVE_GOOCANVAS
std::vector<coot::geometry_graph_block_info_generic>
graphics_info_t::ncs_diffs(int imol, const coot::ncs_chain_difference_t &d) {
   std::vector<coot::geometry_graph_block_info_generic> v;

   std::cout << "peer chain id in ncs_diffs: " << d.peer_chain_id << std::endl;
   for (unsigned int ires=0; ires<d.residue_info.size(); ires++) {
      if (d.residue_info[ires].filled) {
	 float distance = d.residue_info[ires].mean_diff;
	 mmdb::Atom *at = molecules[imol].atom_intelligent(d.peer_chain_id,
						      d.residue_info[ires].resno,
						      d.residue_info[ires].inscode);
	 std::string atom_name = " CA ";
	 std::string altconf = "";
	 if (at) {
	    atom_name = at->name;
	    altconf = at->altLoc;
	 }

	 coot::atom_spec_t as(d.peer_chain_id,
			      d.residue_info[ires].resno,
			      d.residue_info[ires].inscode,
			      atom_name, altconf);
	 std::string str = coot::util::int_to_string(d.residue_info[ires].resno);
	 str += d.peer_chain_id;
	 str += " mean d = ";
	 str += coot::util::float_to_string(distance);
	 str += "A";
	 coot::geometry_graph_block_info_generic block(imol, d.residue_info[ires].resno, as,
						       30.0*distance, str);
	 v.push_back(block);
      }
   }
   return v;
}
#endif

#ifdef HAVE_GOOCANVAS
std::vector<coot::geometry_graph_block_info_generic>
graphics_info_t::ncs_diffs_from_mol(int imol) {

   int imodel = 1;
   std::vector<coot::geometry_graph_block_info_generic> drv;
   std::string altconf("");  // use this (e.g. "A") or "".

   if (is_valid_model_molecule(imol)) {
      std::pair<bool, std::string> master_info = graphics_info_t::molecules[imol].first_ncs_master_chain_id();
      if (master_info.first) {
	 std::string master = master_info.second; // a target_chain_id
	 float w = 1.0; // main chain weight
	 coot::ncs_differences_t diff = graphics_info_t::molecules[imol].ncs_chain_differences(master, w);

	 mmdb::Manager *mol = molecules[imol].atom_sel.mol;
	 int n_chains = diff.diffs.size();
	 int max_chain_length = coot::util::max_min_max_residue_range(mol);
	 coot::geometry_graphs *graphs =
	    new coot::geometry_graphs(coot::GEOMETRY_GRAPH_NCS_DIFFS, imol,
				      graphics_info_t::molecules[imol].name_for_display_manager(),
				      n_chains, max_chain_length);

	 // ncs_diffs_graph[imol] = graphs->dialog();
	 set_validation_graph(imol, coot::GEOMETRY_GRAPH_NCS_DIFFS, graphs->dialog());
	 for (unsigned int incs_set=0; incs_set<diff.diffs.size(); incs_set++) {

	    // do this for each chain
	    int min_resno =  99999;
	    int max_resno = -99999;
	    int offset = 0;

	    // diffs.diffs is a vector of chain differences (vector of ncs_chain_differences_t)
	    // A ncs_chain_differences_t contains a vector of residue difference infos.

	    for (unsigned int ires=0; ires<diff.diffs[incs_set].residue_info.size(); ires++) {
	       if (0)
		  std::cout << "DEBUG:: resno for diffs: "
			    << diff.diffs[incs_set].residue_info[ires].resno
			    << std::endl;

	       if (diff.diffs[incs_set].residue_info[ires].resno < min_resno) {
		  min_resno = diff.diffs[incs_set].residue_info[ires].resno;
	       }
	       if (diff.diffs[incs_set].residue_info[ires].resno > max_resno) {
		  max_resno = diff.diffs[incs_set].residue_info[ires].resno;
	       }
	    }
	    offset = min_resno - 1;
	    //       std::cout << "max_resno, min_resno " << max_resno << " "
	    // 		<< min_resno << std::endl;

	    std::vector<coot::geometry_graph_block_info_generic> v =
	       graphics_info_t::ncs_diffs(imol, diff.diffs[incs_set]);
	    graphs->render_to_canvas(v, incs_set, diff.diffs[incs_set].peer_chain_id,
				     max_resno, min_resno, offset);
	 }
      }
   }
   return drv;
}
#endif
