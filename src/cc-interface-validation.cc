/*
 * src/cc-interface-validation.cc
 *
 * Copyright 2025 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
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
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#include "cc-interface.hh"
#include "graphics-info.h"
#include "coot-utils/json.hpp"
#include "geometry/residue-and-atom-specs.hh"

using json = nlohmann::json;


class interesting_position_button_t {
public:
   coot::Cartesian position;
   std::string label;
   int button_index;
   std::vector<std::string> actions;
   interesting_position_button_t() : button_index(-1) {}
   interesting_position_button_t(const coot::Cartesian &p,
                                 const std::string &l,
                                 int button_idx) : position(p), label(l), button_index(button_idx) {}
};

class atom_spec_pair_t {
public:
   coot::atom_spec_t spec_1;
   coot::atom_spec_t spec_2;
   std::string label;
   int button_index;
   atom_spec_pair_t() {}
   atom_spec_pair_t(const coot::atom_spec_t &a1, const coot::atom_spec_t &a2) : spec_1(a1), spec_2(a2) {}
};

void show_interesting_positions_dialog(int imol,
                                       const std::string &title,
                                       const std::vector<coot::atom_spec_t> &atom_specs,
                                       const std::vector<atom_spec_pair_t> &atom_pair_specs,
                                       const std::vector<coot::residue_spec_t> &residue_specs,
                                       const std::vector<interesting_position_button_t> &positions) {

   auto atom_spec_to_position = [] (int imol, const coot::atom_spec_t &atom_spec) {
      bool status = false;
      clipper::Coord_orth co(0,0,0);
      mmdb::Atom *at = graphics_info_t::molecules[imol].get_atom(atom_spec);
      if (at) {
         status = true;
         co = clipper::Coord_orth(at->x, at->y, at->z);
      }
      return std::make_pair(status, co);
   };

   class button_info_t : public labelled_button_info_t {
   public:
      GtkWidget *button;
      button_info_t(const std::string &l, const clipper::Coord_orth &co, GtkWidget *b) :
         labelled_button_info_t(l, co), button(b) {}
   };

   // This is a bit of a strange "encoding" - the button index is inside the atom and residue spes
   // and positions - I will try to keep the button order

   auto button_callback = +[] (GtkButton *button, gpointer user_data) {
      coot::Cartesian *cc = static_cast<coot::Cartesian *>(user_data);
      bool force_jump = true; // No slide, so that the following updates are done
                              // when we reach the new centre.
      graphics_info_t g;
      bool have_jumped = g.setRotationCentre(*cc, force_jump);
      bool do_updates_now = have_jumped;
      if (do_updates_now) {
         g.update_things_on_move();
      }
      g.graphics_draw();
   };

   auto make_button_label = [] (const std::string &l, float badness) {
      if (badness < 0) {
         // unassigned
         return l;
      } else {
         std::string col = "red";
         if (badness < 0.7) col = "orange";
         if (badness < 0.6) col = "yellow";
         if (badness < 0.4) col = "yellowgreen";
         if (badness < 0.2) col = "greenyellow";
         std::string ll = l + "   <span foreground='";
         ll += col;
         ll +="'>▆</span>";
         return ll;
      }
   };

   auto make_button = [make_button_label, button_callback] (const std::string &l, float badness,
                                                            const clipper::Coord_orth &co) {
      GtkWidget *button = gtk_button_new();
      std::string ll = make_button_label(l, badness);
      GtkWidget *label = gtk_label_new(ll.c_str());
      gtk_label_set_use_markup(GTK_LABEL(label), TRUE);
      gtk_label_set_markup(GTK_LABEL(label), ll.c_str());
      gtk_widget_set_halign(label, GTK_ALIGN_START);
      gtk_button_set_child(GTK_BUTTON(button), label);
      gtk_widget_set_margin_start (button, 4);
      gtk_widget_set_margin_end   (button, 4);
      gtk_widget_set_margin_top   (button, 2);
      gtk_widget_set_margin_bottom(button, 2);
      coot::Cartesian cc(co.x(), co.y(), co.z());
      coot::Cartesian *ccp = new coot::Cartesian(cc);
      g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(button_callback), ccp);
      return button;
   };

   // c.f. void overlaps_peptides_cbeta_ramas_and_rotas_internal(int imol)

   std::vector<button_info_t> buttons;
   if (! graphics_info_t::is_valid_model_molecule(imol)) return;
   mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;
   if (mol == nullptr) return;

   {

      for (unsigned int i=0; i<atom_specs.size(); i++) {
         int idx_button_for_this_atom_spec = atom_specs[i].int_user_data;
         if (true) {
            const coot::atom_spec_t &atom_spec = atom_specs[i];
            const std::string &l = atom_spec.string_user_data;
            std::pair<bool, clipper::Coord_orth> co = atom_spec_to_position(imol, atom_spec);
            if (co.first) {
               float badness =  -1; // unassigned
               if (atom_spec.int_user_data == 1) {
                  badness = atom_spec.float_user_data;
               }
               GtkWidget *button = make_button(l, badness, co.second);
               button_info_t bi(l, co.second, button);
               buttons.push_back(bi);
            }
         }
      }
      for (unsigned int i=0; i<atom_pair_specs.size(); i++) {
         int idx_button_for_this_atom_spec = atom_pair_specs[i].button_index;
         if (true) {
            const coot::atom_spec_t &atom_spec_1 = atom_pair_specs[i].spec_1;
            const coot::atom_spec_t &atom_spec_2 = atom_pair_specs[i].spec_2;
            const std::string &l = atom_spec_1.string_user_data;
            std::cout << "Button for atom spec pair B  " << " " << atom_spec_1 << " "  << atom_spec_2
                      << " " << l << std::endl;
            std::pair<bool, clipper::Coord_orth> co_1 = atom_spec_to_position(imol, atom_spec_1);
            std::pair<bool, clipper::Coord_orth> co_2 = atom_spec_to_position(imol, atom_spec_2);
            if (co_1.first) {
               if (co_2.first) {
                  float badness = -1; // as yet unassigned
                  clipper::Coord_orth combined = co_1.second + co_2.second;
                  clipper::Coord_orth mid_point(combined * 0.5);
                  GtkWidget *button = make_button(l, badness, mid_point);
                  button_info_t bi(l, mid_point, button);
                  buttons.push_back(bi);
               }
            }
         }
      }
      for (unsigned int i=0; i<residue_specs.size(); i++) {
         if (true) {
            const coot::residue_spec_t &res_spec = residue_specs[i];
            std::cout << "handling residue spec B " << i << " " << res_spec << std::endl;
            const std::string &l = res_spec.string_user_data;
            std::pair<bool, clipper::Coord_orth> co = coot::util::get_residue_mid_point(mol, res_spec);
            std::cout << "debug:: res_spec " << res_spec << " midpoint " << co.first << " " << co.second.format()
                      << std::endl;
            if (co.first) {
               std::cout << "handling residue spec C " << std::endl;
               float badness = -1; // as yet unassigned
               if (res_spec.int_user_data == 1) {
                  badness = res_spec.float_user_data;
               }
               GtkWidget *button = make_button(l, badness, co.second);
               button_info_t bi(l, co.second, button);
               buttons.push_back(bi);
            }
         }
      }
      for (unsigned int i=0; i<positions.size(); i++) {
         if (true) {
            const std::string &l = positions[i].label;
            const coot::Cartesian cc = positions[i].position;
            clipper::Coord_orth pos(cc.x(), cc.y(), cc.z());
            float badness = -1; // as yet unassigned
            GtkWidget *button = make_button(l, badness, pos);
            button_info_t bi(l, pos, button);
            buttons.push_back(bi);
         }
      }
   }

   graphics_info_t g;

   GtkWidget *mwravp  = widget_from_builder("main_window_ramchandran_and_validation_pane");
   GtkWidget *mwgrvgp = widget_from_builder("main_window_graphics_rama_vs_graphics_pane");
   GtkWidget* vbox_vbox = widget_from_builder("validation_boxes_vbox");
   GtkWidget* vvf = widget_from_builder("main_window_vertical_validation_frame");
   gtk_widget_set_visible(vbox_vbox, TRUE);

   GtkWidget *vbox_outer  = widget_from_builder("interesting_places_outer_vbox");
   GtkWidget *vbox  = widget_from_builder("interesting_places_vbox");
   GtkWidget *label = widget_from_builder("interesting_places_label");

   int pos = gtk_paned_get_position(GTK_PANED(mwgrvgp));
   if (pos < 200)
      gtk_paned_set_position(GTK_PANED(mwgrvgp), 200);

   gtk_label_set_text(GTK_LABEL(label), title.c_str());
   gtk_widget_set_visible(mwravp,     TRUE);
   gtk_widget_set_visible(mwgrvgp,    TRUE);
   gtk_widget_set_visible(vvf,        TRUE);
   gtk_widget_set_visible(vbox_outer, TRUE);
   gtk_widget_set_visible(vbox,       TRUE);
   gtk_widget_set_visible(label,      TRUE);

   for (auto &button : buttons) {
      // std::cout << "appending button " << button.button << std::endl;
      gtk_box_append(GTK_BOX(vbox), button.button);
   }


}

void read_interesting_places_json_file(const std::string &file_name) {

   bool debug = true;
   if (coot::file_exists(file_name)) {

      graphics_info_t g;
      std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = g.active_atom_spec_simple();
      if (! pp.first) {
         // log
         std::cout << "WARNING:: no molecule found" << std::endl;
      } else {
         int imol = pp.second.first;

         std::vector<coot::atom_spec_t> atom_specs;
         std::vector<atom_spec_pair_t> atom_spec_pairs;
         std::vector<coot::residue_spec_t> residue_specs;
         std::vector<interesting_position_button_t> positions;

         auto get_atom_spec = [] (const json &j) {
            coot::atom_spec_t spec;
            if (j.size() == 5) {
               const json &chain_id_item  = j[0];
               const json &res_no_item    = j[1];
               const json &ins_code_item  = j[2];
               const json &atom_name_item = j[3];
               const json &alt_conf_item  = j[4];
               std::string chain_id  =  chain_id_item.get<std::string>();
               int res_no            =    res_no_item.get<int>();
               std::string ins_code  =  ins_code_item.get<std::string>();
               std::string atom_name = atom_name_item.get<std::string>();
               std::string alt_conf  =  alt_conf_item.get<std::string>();
               spec = coot::atom_spec_t(chain_id, res_no, ins_code, atom_name, alt_conf);
            }
            return spec;
         };

         auto get_residue_spec = [] (const json &j) {

            // std::cout << "Here in get_residue_spec() --- start ---" << std::endl;
            coot::residue_spec_t spec;
            // std::cout << "Here in get_residue_spec() j.size " << j.size() << std::endl;
            if (j.size() == 3) {
               const json &chain_id_item  = j[0];
               const json &res_no_item    = j[1];
               const json &ins_code_item  = j[2];
               std::string chain_id  =  chain_id_item.get<std::string>();
               int res_no            =    res_no_item.get<int>();
               std::string ins_code  =  ins_code_item.get<std::string>();
               spec = coot::residue_spec_t(chain_id, res_no, ins_code);
            }
            if (j.size() == 5) {
               const json &chain_id_item  = j[0];
               const json &res_no_item    = j[1];
               const json &ins_code_item  = j[2];
               std::string chain_id  =  chain_id_item.get<std::string>();
               int res_no            =    res_no_item.get<int>();
               std::string ins_code  =  ins_code_item.get<std::string>();
               // anomaly in the input file
               if (ins_code == " ") ins_code = "";
               spec = coot::residue_spec_t(chain_id, res_no, ins_code);
            }
            return spec;
         };

         std::string s;
         std::fstream f(file_name);
         f.seekg(0, std::ios::end);
         s.reserve(f.tellg());
         f.seekg(0, std::ios::beg);
         s.assign((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());

         std::string title = "<Title>";

         try {
            json j = json::parse(s);
            unsigned int n_outer = j.size();
            if (debug)
               std::cout << "Found " << n_outer << " parts in the json" << std::endl;
            json::const_iterator j_title = j.find("title");
            if (j_title != j.end()) {
               title = j_title.value();
            }
            json j_sections = j["sections"];
            unsigned int n_sections = j_sections.size();
            if (true) {
               if (debug)
                  std::cout << "here A n_sections " << n_sections << std::endl;
               for (std::size_t i=0; i<n_sections; i++) {
                  const json &j_section = j_sections[i];
                  std::string section_title;
                  json::const_iterator it_s = j_section.find(std::string("title"));
                  if (it_s != j_section.end()) section_title = it_s.value();

                  // now iterate through the items of a section

                  json j_items = j_section["items"];
                  unsigned int n_items = j_items.size();

                  for (std::size_t i=0; i<n_items; i++) {
                     if (debug)
                        std::cout << "item " << i << " of " << n_items << std::endl;
                     const json &j_item = j_items[i];
                     json::const_iterator it_1 = j_item.find(std::string("position-type"));
                     if (it_1 != j_item.end()) {
                        std::string position_type = it_1.value();
                        std::string label;
                        coot::Cartesian position(0,0,0);
                        json::const_iterator it_2 = j_item.find(std::string("label"));
                        if (it_2 != j_item.end()) {
                           if (debug)
                              std::cout << "found the label" << std::endl;
                           label = it_2.value();
                        } else {
                           std::cout << "DEBUG:: " << i << " didn't find the label" << std::endl;
                        }

                        if (position_type == "by-atom-spec") {
                           json::const_iterator it_3 = j_item.find(std::string("atom-spec"));
                           if (it_3 != j_item.end()) {
                              const json &j_atom_spec = *it_3;
                              coot::atom_spec_t atom_spec = get_atom_spec(j_atom_spec);
                              if (atom_spec.chain_id != "unset") {
                                 atom_spec.string_user_data = label;
                                 atom_spec.int_user_data = i;
                                 it_3 = j_item.find(std::string("badness")); // optional
                                 if (it_3 != j_item.end()) {
                                    float badness = it_3->get<float>();
                                    atom_spec.float_user_data = badness;
                                    atom_spec.int_user_data = 1; // float user data was set
                                 }
                                 atom_specs.push_back(atom_spec);
                              }
                           }
                        }

                        if (position_type == "by-atom-spec-pair") {
                           // std::cout << "--- found a by-atom-spec-pair" << std::endl;
                           coot::atom_spec_t atom_1_spec;
                           coot::atom_spec_t atom_2_spec;
                           std::string label;
                           json::const_iterator it_2 = j_item.find(std::string("label"));
                           if (it_2 != j_item.end()) {
                              label = it_2.value();
                           }
                           json::const_iterator it_3 = j_item.find(std::string("atom-1-spec"));
                           if (it_3 != j_item.end()) {
                              const json &j_atom_spec = *it_3;
                              atom_1_spec = get_atom_spec(j_atom_spec);
                              if (atom_1_spec.chain_id != "unset") {
                                 atom_1_spec.string_user_data = label;
                                 atom_1_spec.int_user_data = i;
                              }
                           }
                           it_3 = j_item.find(std::string("atom-2-spec"));
                           if (it_3 != j_item.end()) {
                              const json &j_atom_spec = *it_3;
                              atom_2_spec = get_atom_spec(j_atom_spec);
                              if (atom_2_spec.chain_id != "unset") {
                                 atom_2_spec.string_user_data = label;
                                 atom_2_spec.int_user_data = i;
                              }
                           }
                           if (atom_1_spec.chain_id != "unset") {
                              if (atom_2_spec.chain_id != "unset") {
                                 atom_spec_pair_t asp(atom_1_spec, atom_2_spec);
                                 asp.label = label;
                                 asp.button_index = i;
                                 atom_spec_pairs.push_back(asp);
                              }
                           }
                        }

                        if (position_type == "by-residue-spec") {
                           // std::cout << "Here I --- found a by-residue-spec" << std::endl;
                           json::const_iterator it_3 = j_item.find(std::string("residue-spec"));
                           if (it_3 != j_item.end()) {
                              if (debug)
                                 std::cout << "found residue-spec" << std::endl;
                              const json &j_residue_spec = *it_3;
                              coot::residue_spec_t residue_spec = get_residue_spec(j_residue_spec);
                              if (residue_spec.chain_id != "unset") {
                                 residue_spec.string_user_data = label;
                                 residue_spec.int_user_data = i;
                                 it_3 = j_item.find(std::string("badness")); // optional
                                 if (it_3 != j_item.end()) {
                                    float badness = it_3->get<float>();
                                    residue_spec.float_user_data = badness;
                                    residue_spec.int_user_data = 1; // float user data was set
                                 }
                                 residue_specs.push_back(residue_spec);
                              }
                           } else {
                              if (debug)
                                 std::cout << "DEBUG:: didn't find residue spec" << std::endl;
                           }
                        }

                        if (position_type == "by-atom-spec") {
                           if (debug)
                              std::cout << "Here J --- found a by-atom-spec" << std::endl;
                           json::const_iterator it_3 = j_item.find(std::string("atom-spec"));
                           if (it_3 != j_item.end()) {
                              const json &j_atom_spec = *it_3;
                              coot::atom_spec_t atom_spec = get_atom_spec(j_atom_spec);
                              if (atom_spec.chain_id != "unset") {
                                 atom_spec.string_user_data = label;
                                 atom_spec.int_user_data = i;
                                 it_3 = j_item.find(std::string("badness")); // optional
                                 if (it_3 != j_item.end()) {
                                    float badness = it_3->get<float>();
                                    atom_spec.float_user_data = badness;
                                    atom_spec.int_user_data = 1; // float user data was set
                                 }
                                 atom_specs.push_back(atom_spec);
                              }
                           }
                        }

                        if (position_type == "by-coordinates") {
                           // std::cout << "--- found a by-coordinates" << std::endl;
                           json::const_iterator it_3 = j_item.find(std::string("position"));
                           if (it_3 != j_item.end()) {
                              const json &j_pos = *it_3;
                              unsigned int l = j_pos.size();
                              if (l == 3) {
                                 std::cout << "Found a position" << std::endl;
                                 const json &x_item = j_pos[0];
                                 const json &y_item = j_pos[1];
                                 const json &z_item = j_pos[2];
                                 float x = x_item.get<float>();
                                 float y = y_item.get<float>();
                                 float z = z_item.get<float>();
                                 coot::Cartesian c(x,y,z);
                                 // try to find "badness" here?
                                 interesting_position_button_t ipb(c, label, i);
                                 positions.push_back(ipb);
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
         catch(const nlohmann::detail::type_error &e) {
            std::cout << "ERROR:: " << e.what() << std::endl;
         }
         catch(const nlohmann::detail::parse_error &e) {
            std::cout << "ERROR:: " << e.what() << std::endl;
         }


         std::cout << "debug:: ------ " << std::endl;
         std::cout << "debug:: atom_specs " << atom_specs.size() << std::endl;
         std::cout << "debug:: atoms-spec-pairs " << atom_spec_pairs.size() << std::endl;
         std::cout << "debug:: residue_specs " << residue_specs.size() << std::endl;
         std::cout << "debug:: positions " << positions.size() << std::endl;
         if (! atom_specs.empty() || ! residue_specs.empty() || positions.empty()) {
            show_interesting_positions_dialog(imol, title, atom_specs, atom_spec_pairs, residue_specs, positions);
         }
      }
   } else {
      std::cout << "File does not exist " << file_name << std::endl;
   }
}

// not really a validation function - a new file cc-interface-analysis.cc is needed?
std::vector<std::string> get_types_in_molecule(int imol) {

   std::vector<std::string> v;
   graphics_info_t g;
   if (g.is_valid_model_molecule(imol)) {
      v = graphics_info_t::molecules[imol].get_types_in_molecule();
   }
   return v;
}
