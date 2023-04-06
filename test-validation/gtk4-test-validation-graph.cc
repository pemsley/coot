
// use: ./gtk4-test-validation-graph tutorial-modern.pdb rnasa-1.8-all_refmac1.mtz
// there files are in the data directory.

#include "coot-utils/atom-selection-container.hh"
#include "coot-utils/coot-map-utils.hh"
#include "residue-validation-information.hh"
#include "validation-information.hh"
#include "validation-graph-widget.hh"
#include "ligand/rotamer.hh"
#include <gtk/gtk.h>

#include <map>
#include <vector>

#include "coords/ramachandran-validation.hh"

bool
read_mtz(const std::string &file_name,
         const std::string &f, const std::string &phi, const std::string &weight,
         bool use_weight, bool is_a_difference_map,
         clipper::Xmap<float> *xmap_p) {

   bool status = coot::util::map_fill_from_mtz(xmap_p, file_name, f, phi, weight, use_weight, is_a_difference_map);
   return status;
}



coot::validation_information_t
density_fit_analysis(const std::string &pdb_file_name, const std::string &mtz_file_name) {

   coot::validation_information_t r;
   r.name = "Density fit analysis";

   // fill these
   mmdb::PResidue *SelResidues = 0;
   int nSelResidues = 0;

   auto atom_sel = get_atom_selection(pdb_file_name, true, false, false);
   int selHnd = atom_sel.mol->NewSelection(); // yes, it's deleted.
   int imod = 1; // multiple models don't work on validation graphs

   atom_sel.mol->Select(selHnd, mmdb::STYPE_RESIDUE, imod,
                        "*", // chain_id
                        mmdb::ANY_RES, "*",
                        mmdb::ANY_RES, "*",
                        "*",  // residue name
                        "*",  // Residue must contain this atom name?
                        "*",  // Residue must contain this Element?
                        "*",  // altLocs
                        mmdb::SKEY_NEW // selection key
                        );
   atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);

   clipper::Xmap<float> xmap;
   bool status = read_mtz(mtz_file_name, "FWT", "PHWT", "W", 0, 0, &xmap);

   if (! status) {

      std::cout << "Bad mtz file read " << mtz_file_name << std::endl;

   } else {

      for (int ir=0; ir<nSelResidues; ir++) {
         mmdb::Residue *residue_p = SelResidues[ir];
         coot::residue_spec_t res_spec(residue_p);
         mmdb::PAtom *residue_atoms=0;
         int n_residue_atoms;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
         double residue_density_score =
            coot::util::map_score(residue_atoms, n_residue_atoms, xmap, 1);
         //std::string l = res_spec.label();
         std::string l = "Chain ID: "+res_spec.chain_id+"     Residue number: "+std::to_string(res_spec.res_no);
         std::string atom_name = coot::util::intelligent_this_residue_mmdb_atom(residue_p)->GetAtomName();
         const std::string &chain_id = res_spec.chain_id;
         int this_resno = res_spec.res_no;
         coot::atom_spec_t atom_spec(chain_id, this_resno, res_spec.ins_code, atom_name, "");
         coot::residue_validation_information_t rvi(res_spec, atom_spec, residue_density_score, l);
         r.add_residue_validation_information(rvi, chain_id);
      }
      atom_sel.mol->DeleteSelection(selHnd);
   }
   return r;
}

coot::validation_information_t
density_correlation(const std::string &pdb_file_name, const std::string &mtz_file_name) {

   coot::validation_information_t r;

   auto atom_sel = get_atom_selection(pdb_file_name, true, false, false);
   if (atom_sel.read_success) {
      clipper::Xmap<float> xmap;
      bool status = read_mtz(mtz_file_name, "FWT", "PHWT", "W", 0, 0, &xmap);
      if (status) {

         // signature:
         //
         // std::vector<std::pair<coot::residue_spec_t, float> >
         //    coot::util::map_to_model_correlation_per_residue(mmdb::Manager *mol,
         //                                                     const std::vector<coot::residue_spec_t> &specs,
         //                                                     unsigned short int atom_mask_mode,
         //                                                     float atom_radius, // for masking
         //                                                     const clipper::Xmap<float> &reference_map)

         // atom_mask_mode:
         //
         // 0: all-atoms
         // 1: main-chain atoms if is standard amino-acid, else all atoms
         // 2: side-chain atoms if is standard amino-acid, else all atoms
         // 3: side-chain atoms-exclusing CB if is standard amino-acid, else all atoms

         unsigned short int atom_mask_mode = 0;
         float atom_radius = 2.0;

         std::vector<coot::residue_spec_t> residue_specs;
         std::vector<mmdb::Residue *> residues = coot::util::residues_in_molecule(atom_sel.mol);
         for (unsigned int i=0; i<residues.size(); i++)
            residue_specs.push_back(coot::residue_spec_t(residues[i]));

         std::vector<std::pair<coot::residue_spec_t, float> > correlations =  
            coot::util::map_to_model_correlation_per_residue(atom_sel.mol,
                                                             residue_specs,
                                                             atom_mask_mode,
                                                             atom_radius, // for masking
                                                             xmap);

         std::vector<std::pair<coot::residue_spec_t, float> >::const_iterator it;
         for (it=correlations.begin(); it!=correlations.end(); ++it) {
            const auto &r_spec(it->first);
            const auto &correl(it->second);

            std::string atom_name = " CA ";
            coot::atom_spec_t atom_spec(r_spec.chain_id, r_spec.res_no, r_spec.ins_code, atom_name, "");
            std::string label = "Correl: ";
            coot::residue_validation_information_t rvi(r_spec, atom_spec, correl, label);
            r.add_residue_validation_information(rvi, r_spec.chain_id);
            
         }
      } else {
         std::cout << "Bad mtz file read " << mtz_file_name << std::endl;
      }
   } else {
      std::cout << "Bad read for pdb file " << pdb_file_name << std::endl;
   }
   return r;
}



#include "ligand/rotamer.hh"

coot::validation_information_t
rotamer_analysis(const std::string &pdb_file_name) {

   coot::validation_information_t r;
   r.name = "Rotamer analysis";

   // fill these
   mmdb::PResidue *SelResidues = 0;
   int nSelResidues = 0;

   auto atom_sel = get_atom_selection(pdb_file_name, true, false, false);
   int selHnd = atom_sel.mol->NewSelection(); // yes, it's deleted.
   int imod = 1; // multiple models don't work on validation graphs

   atom_sel.mol->Select(selHnd, mmdb::STYPE_RESIDUE, imod,
                        "*", // chain_id
                        mmdb::ANY_RES, "*",
                        mmdb::ANY_RES, "*",
                        "*",  // residue name
                        "*",  // Residue must contain this atom name?
                        "*",  // Residue must contain this Element?
                        "*",  // altLocs
                        mmdb::SKEY_NEW // selection key
                        );
   atom_sel.mol->GetSelIndex(selHnd, SelResidues, nSelResidues);


   if (! atom_sel.read_success) {

      std::cout << "Bad pdb file read " << pdb_file_name << std::endl;

   } else {

      for (int ir=0; ir<nSelResidues; ir++) {
         mmdb::Residue *residue_p = SelResidues[ir];
         coot::residue_spec_t res_spec(residue_p);
         mmdb::PAtom *residue_atoms=0;
         int n_residue_atoms;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);

         // double residue_density_score = coot::util::map_score(residue_atoms, n_residue_atoms, xmap, 1);

         if (n_residue_atoms > 5) {

            std::string res_name = residue_p->GetResName();
            if (true) {

               coot::rotamer rot(residue_p);
               coot::rotamer_probability_info_t rpi = rot.probability_of_this_rotamer();
               double prob = rpi.probability;

               std::string l = "Chain ID: "+res_spec.chain_id+"     Residue number: "+std::to_string(res_spec.res_no);
               std::string atom_name = coot::util::intelligent_this_residue_mmdb_atom(residue_p)->GetAtomName();
               const std::string &chain_id = res_spec.chain_id;
               int this_resno = res_spec.res_no;
               coot::atom_spec_t atom_spec(chain_id, this_resno, res_spec.ins_code, atom_name, "");
               coot::residue_validation_information_t rvi(res_spec, atom_spec, prob, l);
               r.add_residue_validation_information(rvi, chain_id);
            }
         }
      }
      atom_sel.mol->DeleteSelection(selHnd);
   }
   return r;
}

struct graphs_shipment_t {
   /// we use tab label as key here
   std::map<std::string,CootValidationGraph*> graphs_for_tabs;
   std::map<std::string,CootValidationGraph*> compacted_graphs;
   GtkComboBoxText* chain_selector;

   graphs_shipment_t(GtkComboBoxText* chain_selector) {
      this->chain_selector = chain_selector;
   }

   void push_graph(coot::validation_information_t&& data) {
      auto data_ptr = std::make_shared<coot::validation_information_t>(data);
      auto* graph_ptr = coot_validation_graph_new();
      std::string name = data_ptr->name;

      graphs_for_tabs[name] = graph_ptr;
      coot_validation_graph_set_validation_information(graph_ptr,data_ptr);
      gtk_widget_set_margin_bottom(GTK_WIDGET(graph_ptr),10);
      gtk_widget_set_margin_start(GTK_WIDGET(graph_ptr),10);
      gtk_widget_set_margin_end(GTK_WIDGET(graph_ptr),10);
      gtk_widget_set_margin_top(GTK_WIDGET(graph_ptr),10);

      graph_ptr = coot_validation_graph_new();
      compacted_graphs[name] = graph_ptr;
      coot_validation_graph_set_validation_information(graph_ptr,data_ptr);
   }
};

GtkWidget* build_graph_vbox(CootValidationGraph* validation_graph) {
   GtkWidget* vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL,10);
   gtk_widget_set_margin_bottom(vbox,10);
   gtk_widget_set_margin_top(vbox,10);
   gtk_widget_set_margin_start(vbox,10);
   gtk_widget_set_margin_end(vbox,10);

   GtkWidget* host_frame = gtk_frame_new(NULL);

   gtk_widget_set_margin_bottom(host_frame,10);
   gtk_widget_set_margin_start(host_frame,10);
   gtk_widget_set_margin_end(host_frame,10);
   gtk_widget_set_margin_top(host_frame,10);

   gtk_frame_set_child(GTK_FRAME(host_frame),GTK_WIDGET(validation_graph));

   GtkWidget* host_scrolled_window = gtk_scrolled_window_new();
   gtk_widget_set_hexpand(host_scrolled_window,TRUE);
   gtk_widget_set_vexpand(host_scrolled_window,TRUE);
   gtk_widget_set_size_request(host_scrolled_window,720,400);
   gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(host_scrolled_window),GTK_WIDGET(host_frame));

   GtkWidget* outer_frame = gtk_frame_new("Container for the experimental Validation Graph Widget");
   gtk_frame_set_child(GTK_FRAME(outer_frame),host_scrolled_window);

   gtk_box_append(GTK_BOX(vbox),outer_frame);
   GtkWidget* target_label = gtk_label_new("");
   gtk_box_append(GTK_BOX(vbox),target_label);

   g_signal_connect(validation_graph,"residue-clicked",
      G_CALLBACK(+[](CootValidationGraph* self, const coot::residue_validation_information_t* residue, gpointer userdata){
         GtkLabel* label = GTK_LABEL(userdata);
         gtk_label_set_text(label,residue->label.c_str());
         g_debug("Inside 'residue-clicked' handler: %s",residue->label.c_str());
      }),
   target_label);

   GtkWidget* scale = gtk_scale_new_with_range(GTK_ORIENTATION_HORIZONTAL, 0.1f, 5.f, 0.1f);
   gtk_box_append(GTK_BOX(vbox), scale);
   gtk_scale_set_draw_value(GTK_SCALE(scale), TRUE);
   gtk_range_set_value(GTK_RANGE(scale),1.f);
   g_signal_connect(scale, "value-changed", G_CALLBACK(+[](GtkScale* scale, gpointer user_data){
      CootValidationGraph* graph = COOT_COOT_VALIDATION_GRAPH(user_data);
      coot_validation_graph_set_horizontal_zoom_scale(graph, gtk_range_get_value(GTK_RANGE(scale)));
   }), validation_graph);

   return vbox;
}

GtkWidget* build_graph_stack(graphs_shipment_t* graphs) {
   GtkWidget* vbox = gtk_box_new(GTK_ORIENTATION_VERTICAL,10);
   gtk_widget_set_margin_bottom(vbox,10);
   gtk_widget_set_margin_top(vbox,10);
   gtk_widget_set_margin_start(vbox,10);
   gtk_widget_set_margin_end(vbox,10);

   GtkWidget* host_frame = gtk_frame_new(NULL);

   gtk_widget_set_margin_bottom(host_frame,10);
   gtk_widget_set_margin_start(host_frame,10);
   gtk_widget_set_margin_end(host_frame,10);
   gtk_widget_set_margin_top(host_frame,10);

   GtkWidget* vbox_inner = gtk_box_new(GTK_ORIENTATION_VERTICAL,0);
   gtk_widget_set_margin_bottom(vbox_inner,10);
   gtk_widget_set_margin_start(vbox_inner,10);
   gtk_widget_set_margin_end(vbox_inner,10);
   gtk_widget_set_margin_top(vbox_inner,10);
   gtk_frame_set_child(GTK_FRAME(host_frame),GTK_WIDGET(vbox_inner));

   GtkWidget* target_label = gtk_label_new("");
   GtkWidget* scale = gtk_scale_new_with_range(GTK_ORIENTATION_HORIZONTAL, 0.1f, 5.f, 0.1f);

   for(auto el : graphs->compacted_graphs) {
      gtk_box_append(GTK_BOX(vbox_inner),GTK_WIDGET(el.second));
      coot_validation_graph_set_single_chain_mode(el.second, "A");
      g_signal_connect(el.second,"residue-clicked",
         G_CALLBACK(+[](CootValidationGraph* self, const coot::residue_validation_information_t* residue, gpointer userdata){
            GtkLabel* label = GTK_LABEL(userdata);
            gtk_label_set_text(label,residue->label.c_str());
            g_debug("Inside 'residue-clicked' handler: %s",residue->label.c_str());
         }),
      target_label);

      g_signal_connect(scale, "value-changed", G_CALLBACK(+[](GtkScale* scale, gpointer user_data){
         CootValidationGraph* graph = COOT_COOT_VALIDATION_GRAPH(user_data);
         coot_validation_graph_set_horizontal_zoom_scale(graph, gtk_range_get_value(GTK_RANGE(scale)));
      }), el.second);

      g_signal_connect(graphs->chain_selector, "changed",G_CALLBACK(+[](GtkComboBoxText* selector, gpointer user_data){
         CootValidationGraph* graph = COOT_COOT_VALIDATION_GRAPH(user_data);
         coot_validation_graph_set_single_chain_mode(graph, gtk_combo_box_text_get_active_text(selector));
      }), el.second);
   }

   GtkWidget* host_scrolled_window = gtk_scrolled_window_new();
   gtk_widget_set_hexpand(host_scrolled_window,TRUE);
   gtk_widget_set_vexpand(host_scrolled_window,TRUE);
   gtk_widget_set_size_request(host_scrolled_window,720,400);
   gtk_scrolled_window_set_child(GTK_SCROLLED_WINDOW(host_scrolled_window),GTK_WIDGET(host_frame));

   GtkWidget* outer_frame = gtk_frame_new("Container for the experimental Validation Graph Widget [stacked display]");
   gtk_frame_set_child(GTK_FRAME(outer_frame),host_scrolled_window);

   gtk_box_append(GTK_BOX(vbox),outer_frame);
   
   gtk_box_append(GTK_BOX(vbox),target_label);

   gtk_box_append(GTK_BOX(vbox), scale);
   gtk_scale_set_draw_value(GTK_SCALE(scale), TRUE);
   gtk_range_set_value(GTK_RANGE(scale),1.f);


   gtk_box_append(GTK_BOX(vbox), GTK_WIDGET(graphs->chain_selector));

   return vbox;
}

void build_main_window(GtkWindow* main_window, graphs_shipment_t* graphs) {
   
   GtkWidget* graph_notebook = gtk_notebook_new();
   gtk_window_set_child(main_window,graph_notebook);

   for(auto el : graphs->graphs_for_tabs) {
      gtk_notebook_append_page(GTK_NOTEBOOK(graph_notebook), build_graph_vbox(el.second), gtk_label_new(el.first.c_str()));
   }
   gtk_notebook_append_page(GTK_NOTEBOOK(graph_notebook), build_graph_stack(graphs), gtk_label_new("Stacked view"));

}

coot::validation_information_t
ramachandran_analysis(const std::string &pdb_file_name) {

   // internals copied from the function of the same name in molecules_container.cc

   coot::validation_information_t vi;
   vi.name = "Ramachandran analysis";

   auto atom_sel = get_atom_selection(pdb_file_name, true, false, false);
   mmdb::Manager *mol = atom_sel.mol;

   const ramachandrans_container_t rc;
   std::vector<coot::phi_psi_prob_t> rv = coot::ramachandran_validation(mol, rc);
   for (unsigned int i=0; i<rv.size(); i++) {
      std::string chain_id = rv[i].phi_psi.chain_id;
      coot::residue_spec_t residue_spec(rv[i].phi_psi.chain_id, rv[i].phi_psi.residue_number, rv[i].phi_psi.ins_code);
      double pr = rv[i].probability;
      std::string label = rv[i].phi_psi.chain_id + std::string(" ") + std::to_string(rv[i].phi_psi.residue_number);
      if (! rv[i].phi_psi.ins_code.empty())
         label += std::string(" ") + rv[i].phi_psi.ins_code;
      coot::atom_spec_t atom_spec(residue_spec.chain_id, residue_spec.res_no, residue_spec.ins_code, " CA ", "");
      coot::residue_validation_information_t rvi(residue_spec, atom_spec, pr, label);
      if (false)
         std::cout << "         " << residue_spec << " " << rv[i].phi_psi.phi() << " " << rv[i].phi_psi.psi()
                   << " pr " << pr << " " << std::endl;
      vi.add_residue_validation_information(rvi, chain_id);
   }
   vi.set_min_max();

   return vi;
}

#include "ideal/simple-restraint.hh"

coot::validation_information_t
peptide_omega_analysis(const std::string &pdb_file_name) {

   coot::validation_information_t vi;
   vi.name = "Peptide Omega Analysis";
   coot::protein_geometry geom;
   auto atom_sel = get_atom_selection(pdb_file_name, true, false, false);
   if (! atom_sel.read_success) return vi;
   mmdb::Manager *mol = atom_sel.mol;
   int imodel = 1;
   mmdb::Model *model_p = mol->GetModel(imodel);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         std::cout << "ichain: " << ichain << std::endl;
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         std::cout << "ichain: " << ichain << " " << chain_p << std::endl;
         std::string chain_id(chain_p->GetChainID());
         coot::restraints_container_t rc(atom_sel, chain_id, nullptr);
         coot::omega_distortion_info_container_t odi = rc.omega_trans_distortions(geom, true);
         std::cout << "odi: chain_id "  << odi.chain_id << std::endl;
         std::cout << "odi: min_resno " << odi.min_resno << std::endl;
         std::cout << "odi: max_resno " << odi.max_resno << std::endl;
         std::cout << "odi: n omega_distortions " << odi.omega_distortions.size() << std::endl;

         coot::chain_validation_information_t cvi(chain_id);
         for (const auto &od : odi.omega_distortions) {
            coot::residue_spec_t res_spec(chain_id, od.resno, "");
            coot::atom_spec_t atom_spec(chain_id, od.resno, "", " CA ", "");
            std::string label = od.info_string;
            coot::residue_validation_information_t rvi(res_spec, atom_spec, od.distortion, label);
            cvi.add_residue_validation_information(rvi);
         }
         vi.cviv.push_back(cvi);
      }
   }
   return vi;
}

int main(int argc, char **argv) {

   if (argc > 2) {
      std::string pdb_file_name = argv[1];
      std::string mtz_file_name = argv[2];
      coot::validation_information_t vid = density_fit_analysis(pdb_file_name, mtz_file_name);
      coot::validation_information_t vir = rotamer_analysis(pdb_file_name);
      coot::validation_information_t vit = ramachandran_analysis(pdb_file_name);
      coot::validation_information_t vio; // = peptide_omega_analysis(pdb_file_name); smashes the stack.
      // for now, so that the graphs display correctly
      vio.name = "Peptide Omega Analysis";

      // now do something (i.e. make a pretty interactive graph) with vid and vir.

      for (const auto &cvi : vid.cviv) {
         std::cout << "Chain " << cvi.chain_id << std::endl;
         for (const auto &ri : cvi.rviv) {
            std::cout << "[Density validation] Residue " << ri.residue_spec << " " << ri.function_value << std::endl;
         }
      }

      for (const auto &cvi : vir.cviv) {
         std::cout << "Chain " << cvi.chain_id << std::endl;
         for (const auto &ri : cvi.rviv) {
            std::cout << "[Rotamer validation] Residue " << ri.residue_spec << " " << ri.function_value << std::endl;
         }
      }

      for (const auto &cvi : vit.cviv) {
         std::cout << "Chain " << cvi.chain_id << std::endl;
         for (const auto &ri : cvi.rviv) {
            std::cout << " [Ramachandran Validation] Residue " << ri.residue_spec << " " << ri.function_value << std::endl;
         }
      }

      for (const auto &cvi : vio.cviv) {
         std::cout << "Chain " << cvi.chain_id << std::endl;
         for (const auto &ri : cvi.rviv) {
            std::cout << " [Peptide Omega Validation] Residue " << ri.residue_spec << " " << ri.function_value << std::endl;
         }
      }


      gtk_init();
      
      GtkApplication* app = gtk_application_new("org.pemsley.Test-validation-graphs",G_APPLICATION_FLAGS_NONE);
      GError *error = NULL;
      g_application_register(G_APPLICATION(app), NULL, &error);


      GtkWidget* chain_selector = gtk_combo_box_text_new();
      for(const auto& chain: vid.cviv) {
         gtk_combo_box_text_append_text(GTK_COMBO_BOX_TEXT(chain_selector), chain.chain_id.c_str());
      }

      graphs_shipment_t* gs = new graphs_shipment_t(GTK_COMBO_BOX_TEXT(chain_selector));

      // for the sake of testing
      vid.type = coot::graph_data_type::Energy;
      vir.type = coot::graph_data_type::Energy;
      vit.type = coot::graph_data_type::Probability;
      // vit.type = coot::graph_data_type::LogProbability;

      gs->push_graph(std::move(vid));
      gs->push_graph(std::move(vir));
      gs->push_graph(std::move(vit));
      gs->push_graph(std::move(vio));

      g_signal_connect(app,"activate",G_CALLBACK(+[](GtkApplication* app, gpointer user_data) {
         //GtkWindow* win = GTK_WINDOW(user_data);
         GtkWidget* win = gtk_application_window_new(app);
         gtk_application_add_window(app,GTK_WINDOW(win));
         gtk_window_set_application(GTK_WINDOW(win),app);
         graphs_shipment_t* gs = static_cast<graphs_shipment_t*>(user_data);
         build_main_window(GTK_WINDOW(win),gs);
         gtk_widget_show(win);

         delete gs;
      }),gs);


      return g_application_run(G_APPLICATION(app),0,0);
   } else {
     std::cout << "Two commandline args needed.\n";
     return 2;
   }
}
