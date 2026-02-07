/*
 * src/dynamic-validation.cc
 *
 * Copyright 2023 by Medical Research Council
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

#include "coot-utils/coot-coord-utils.hh"
#include "graphics-info.h"
#include "widget-from-builder.hh"
#include "dynamic-validation.hh"
#include "coot-utils/c-beta-deviations.hh"

void dynamic_validation_internal(int imol, int imol_map) {

   overlaps_peptides_cbeta_ramas_and_rotas_internal(imol); // for now

   // map tools here.
}

void update_dynamic_validation() {

   graphics_info_t g;
   std::pair<bool, std::pair<int, coot::atom_spec_t> > pp = g.active_atom_spec_simple();
   if (pp.first) {
      int imol = pp.second.first;
      update_dynamic_validation_for_molecule(imol);
   }
}

void update_dynamic_validation_for_molecule(int imol) {

   // 20250202-PE as of today, this function makes validation_boxes_vbox visible.
   // It is up to the calling function to decide if this function should be called
   // and hence overlaps_peptides_cbeta_ramas_and_rotas_internal()
   //

   // only update if dynamic validation was being displayed already.

   // does main_window_ramchandran_and_validation_pane exist today? I can't find it.
   // If not, what does this function do?
   //
   GtkWidget* pane_1 = widget_from_builder("main_window_ramchandran_and_validation_pane");
   GtkWidget* vbox_2 = widget_from_builder("validation_boxes_vbox");
   GtkWidget *vbox_1  = widget_from_builder("dynamic_validation_outliers_vbox");

   gtk_widget_set_visible(vbox_2, TRUE);
   if (gtk_widget_get_visible(vbox_1)) {
      if (gtk_widget_get_visible(pane_1)) {
         overlaps_peptides_cbeta_ramas_and_rotas_internal(imol);
      } else {
         std::cout << "ERROR:: pane main_window_ramchandran_and_validation_pane not found " << std::endl;
      }
   }
}


void overlaps_peptides_cbeta_ramas_and_rotas_internal(int imol) {

   // add "missing side-chains? " button

   auto to_label = [] (mmdb:: Atom *at) {
      std::string s = std::string(at->GetChainID());
      s += "";
      s += std::to_string(at->GetSeqNum());
      std::string ins_code(at->GetInsCode());
      if (!ins_code.empty()) s += " " + ins_code;
      s += " ";
      s += coot::util::remove_whitespace(std::string(at->GetAtomName()));
      std::string alt_conf(at->altLoc);
      if (!alt_conf.empty()) s += " " + alt_conf;
      return s;
   };

   auto overlap_button_callback = +[] (GtkButton *button, G_GNUC_UNUSED gpointer user_data) {
      std::string target_position_x_str(static_cast<const char *>(g_object_get_data(G_OBJECT(button), "target-position-x")));
      std::string target_position_y_str(static_cast<const char *>(g_object_get_data(G_OBJECT(button), "target-position-y")));
      std::string target_position_z_str(static_cast<const char *>(g_object_get_data(G_OBJECT(button), "target-position-z")));
      try {
         float x = coot::util::string_to_float(target_position_x_str);
         float y = coot::util::string_to_float(target_position_y_str);
         float z = coot::util::string_to_float(target_position_z_str);
         clipper::Coord_orth p(x,y,z);
         graphics_info_t::set_rotation_centre(p);
      }
      catch (const std::runtime_error &e) {
         std::cout << "WARNING::" << e.what() << std::endl;
      }
      graphics_info_t g;
      g.graphics_grab_focus();
   };


   // How do I go to a position from a residue elsewhere in the code!? I don't recall!
   //
   auto rota_button_callback = +[] (GtkButton *button, G_GNUC_UNUSED gpointer user_data) {
      std::string target_position_x_str(static_cast<const char *>(g_object_get_data(G_OBJECT(button), "target-position-x")));
      std::string target_position_y_str(static_cast<const char *>(g_object_get_data(G_OBJECT(button), "target-position-y")));
      std::string target_position_z_str(static_cast<const char *>(g_object_get_data(G_OBJECT(button), "target-position-z")));
      try {
         float x = coot::util::string_to_float(target_position_x_str);
         float y = coot::util::string_to_float(target_position_y_str);
         float z = coot::util::string_to_float(target_position_z_str);
         clipper::Coord_orth p(x,y,z);
         graphics_info_t::set_rotation_centre(p);
      }
      catch (const std::runtime_error &e) {
         std::cout << "WARNING::" << e.what() << std::endl;
      }
      graphics_info_t g;
      g.graphics_grab_focus();
   };

   auto set_target_position_data = [] (GtkWidget *button, const clipper::Coord_orth &p) {
      // std::to_string(p.x()).c_str()
      char *x = new char[10];
      char *y = new char[10];
      char *z = new char[10];
      for (unsigned int i=0; i<10; i++) x[i] = 0;
      for (unsigned int i=0; i<10; i++) y[i] = 0;
      for (unsigned int i=0; i<10; i++) z[i] = 0;
      strncpy(x, std::to_string(p.x()).c_str(), 9);
      strncpy(y, std::to_string(p.y()).c_str(), 9);
      strncpy(z, std::to_string(p.z()).c_str(), 9);
      g_object_set_data(G_OBJECT(button), "target-position-x", x);
      g_object_set_data(G_OBJECT(button), "target-position-y", y);
      g_object_set_data(G_OBJECT(button), "target-position-z", z);
   };


   auto set_target_position_from_residue = [set_target_position_data] (mmdb::Residue *residue_p, GtkWidget *button) {

      std::pair<bool, clipper::Coord_orth> ptcb = coot::util::get_CB_position_in_residue(residue_p);
      clipper::Coord_orth p(0,0,0);
      if (ptcb.first) {
         p = ptcb.second;
      } else {
         // crazy
         std::pair<bool, clipper::Coord_orth> pt = coot::util::get_residue_centre(residue_p);
         if (pt.first) {
            p = pt.second;
         }
      }
      set_target_position_data(button, p);
   };

   auto mmdb_to_clipper = [] (mmdb::Atom *at) {
      return clipper::Coord_orth(at->x, at->y, at->z);
   };

   // using atom overlap
   auto get_target_position = [mmdb_to_clipper] (const coot::atom_overlap_t &ao) {
      clipper::Coord_orth p1 = mmdb_to_clipper(ao.atom_1);
      clipper::Coord_orth p2 = mmdb_to_clipper(ao.atom_2);
      return 0.5 * (p1 + p2);
   };

   GtkWidget *pane_to_show  = widget_from_builder("main_window_ramchandran_and_validation_pane");
   gtk_widget_set_visible(pane_to_show,  TRUE);
   GtkWidget *pane  = widget_from_builder("main_window_graphics_rama_vs_graphics_pane");

   GtkWidget *outer_vbox  = widget_from_builder("dynamic_validation_vbox");
   gtk_widget_set_visible(outer_vbox,  TRUE);

   GtkWidget *vbox  = widget_from_builder("dynamic_validation_outliers_vbox");
   GtkWidget *label = widget_from_builder("dynamic_validation_outliers_label");

   auto make_overlap_buttons = [to_label, get_target_position, set_target_position_data,
                                overlap_button_callback] (mmdb::Manager *mol) {

      double vol_crit = 1.2;
      std::vector<std::pair<coot::residue_spec_t, GtkWidget *> > buttons;
      bool ignore_waters = false;
      coot::atom_overlaps_container_t overlaps(mol, graphics_info_t::Geom_p(), ignore_waters, 0.5, 0.25);
      overlaps.make_all_atom_overlaps();
      std::vector<coot::atom_overlap_t> olv = overlaps.overlaps;

      for (unsigned int ii=0; ii<olv.size(); ii++) {
         const auto &o = olv[ii];
         if (o.overlap_volume < vol_crit) continue;
	 if (false) // debug
	    std::cout << "Overlap " << ii << " "
		      << coot::atom_spec_t(o.atom_1) << " "
		      << coot::atom_spec_t(o.atom_2) << " overlap-vol "
		      << o.overlap_volume << " r_1 "
		      << o.r_1 << " r_2 " << o.r_2 << std::endl;
         std::string lab = "Atom Overlap " + to_label(o.atom_1) + " - " + to_label(o.atom_2);
         lab += " OV: " + coot::util::float_to_string_using_dec_pl(o.overlap_volume, 2);
         lab += "Å³";
         // GtkWidget *button = gtk_button_new_with_label(lab.c_str());
         GtkWidget *button = gtk_button_new();
         GtkWidget *label = gtk_label_new(lab.c_str());
         gtk_widget_set_halign(label, GTK_ALIGN_START);
         gtk_button_set_child(GTK_BUTTON(button), label);
         gtk_widget_set_margin_start (button, 4);
         gtk_widget_set_margin_end   (button, 4);
         gtk_widget_set_margin_top   (button, 2);
         gtk_widget_set_margin_bottom(button, 2);
         clipper::Coord_orth target_position = get_target_position(o);
         set_target_position_data(button, target_position);
         g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(overlap_button_callback), nullptr);
         auto p = std::make_pair(coot::residue_spec_t(o.atom_1->residue), button);
         buttons.push_back(p);
      }
      return buttons;
   };

   auto make_rota_buttons = [imol, set_target_position_from_residue, rota_button_callback] (mmdb::Manager *mol) {
      float prob_critical = 1.0; // from 0 to 1.0;
      std::vector<std::pair<coot::residue_spec_t, GtkWidget *> > buttons;
      std::vector<mmdb::Residue *> residues =  graphics_info_t::molecules[imol].get_all_protein_residues();
      graphics_info_t g;
      float lp = graphics_info_t::rotamer_lowest_probability;
      for(auto &residue_p : residues) {
         std::string alt_conf;
         coot::rotamer_probability_info_t d_score = g.get_rotamer_probability(residue_p, alt_conf, mol, lp, 1);
         if (d_score.state == coot::rotamer_probability_info_t::OK) {
            float f = d_score.probability; // 0 to 100
            // std::cout << "   " << coot::residue_spec_t(residue_p) << " " << f << std::endl;
            if (f < prob_critical) {
               std::string lab = "Rotamer ";
               lab += std::string(residue_p->GetChainID()) + std::string(" ") + std::to_string(residue_p->GetSeqNum());
               lab += " ";
               lab += residue_p->GetResName();
               lab += " ";
               lab += coot::util::float_to_string_using_dec_pl(f, 2);
               lab += "%";
               // GtkWidget *button = gtk_button_new_with_label(lab.c_str());
               GtkWidget *button = gtk_button_new();
               GtkWidget *label = gtk_label_new(lab.c_str());
               gtk_widget_set_halign(label, GTK_ALIGN_START);
               gtk_button_set_child(GTK_BUTTON(button), label);
               set_target_position_from_residue(residue_p, button);
               g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(rota_button_callback), nullptr);
               gtk_widget_set_margin_start (button, 4);
               gtk_widget_set_margin_end   (button, 4);
               gtk_widget_set_margin_top   (button, 2);
               gtk_widget_set_margin_bottom(button, 2);
               auto p = std::make_pair(coot::residue_spec_t(residue_p), button);
               buttons.push_back(p);
            }
         }
      }
      return buttons;
   };

   auto chiral_volume_button_clicked_callback = +[] (GtkButton *button, gpointer user_data) {
      coot::atom_spec_t *atom_spec_p = static_cast<coot::atom_spec_t *>(user_data);
      int imol = atom_spec_p->int_user_data;
      const auto &atom_spec(*atom_spec_p);
      graphics_info_t g;
      graphics_info_t::set_go_to_atom(imol, atom_spec);
      g.try_centre_from_new_go_to_atom(); // oh dear!
      g.graphics_grab_focus();
   };

   auto make_chiral_volume_buttons = [chiral_volume_button_clicked_callback] (int imol) {
      std::vector<std::pair<coot::residue_spec_t, GtkWidget *> > buttons;
      // the first is types with no dictionary
      double chiral_volume_distortion_limit = 6.0;
      std::pair<std::vector<std::string>, std::vector<std::pair<coot::atom_spec_t, double> > > v =
         graphics_info_t::molecules[imol].distorted_chiral_volumes(chiral_volume_distortion_limit);

      for (unsigned int i=0; i<v.second.size(); i++) {
         const auto &spec = v.second[i].first;
         float distortion = v.second[i].second;
         coot::atom_spec_t *spec_p = new coot::atom_spec_t(spec);
         spec_p->int_user_data = imol;
         std::string lab = "Chiral Volume Outlier ";
         coot::residue_spec_t res_spec(*spec_p);
         mmdb::Residue *residue_p = graphics_info_t::molecules[imol].get_residue(res_spec);
         if (residue_p) {
            lab += std::string(residue_p->GetChainID()) + std::string(" ") + std::to_string(residue_p->GetSeqNum());
            lab += " ";
            lab += coot::util::remove_whitespace(std::string(spec_p->atom_name));
            lab += " ";
            lab += residue_p->GetResName();
            lab += " ";
            lab += coot::util::float_to_string_using_dec_pl(distortion, 1);
            // GtkWidget *button = gtk_button_new_with_label(lab.c_str());
            GtkWidget *button = gtk_button_new();
            GtkWidget *label = gtk_label_new(lab.c_str());
            gtk_widget_set_halign(label, GTK_ALIGN_START);
            gtk_button_set_child(GTK_BUTTON(button), label);
            gtk_widget_set_margin_start (button, 4);
            gtk_widget_set_margin_end   (button, 4);
            gtk_widget_set_margin_top   (button, 2);
            gtk_widget_set_margin_bottom(button, 2);
            g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(chiral_volume_button_clicked_callback), spec_p);
            auto p = std::make_pair(res_spec, button);
            buttons.push_back(p);
         }
      }
      return buttons;
   };

   // these are mostly the same and can be consolidated.
   auto rama_button_clicked_callback = +[] (GtkButton *button, gpointer user_data) {
      coot::atom_spec_t *atom_spec_p = static_cast<coot::atom_spec_t *>(user_data);
      int imol = atom_spec_p->int_user_data;
      const auto &atom_spec(*atom_spec_p);
      graphics_info_t g;
      graphics_info_t::set_go_to_atom(imol, atom_spec);
      g.try_centre_from_new_go_to_atom(); // oh dear!
      g.graphics_grab_focus();
   };

   auto make_rama_buttons = [rama_button_clicked_callback] (int imol) {

      double prob_crit = 0.003; // rama_level_allowed is 0.002
      std::vector<std::pair<coot::residue_spec_t, GtkWidget *> > buttons;

      const auto &m = graphics_info_t::molecules[imol];
      coot::rama_score_t rama_score = m.get_all_molecule_rama_score();

      std::vector<coot::rama_score_t::scored_phi_psi_t>::const_iterator it;
      for (it=rama_score.scores.begin(); it!=rama_score.scores.end(); ++it) {
         const auto &spp = *it;
         if (spp.score < prob_crit) {
            if (false)
               std::cout << "debug:: rama  " << spp.res_spec << "  " << spp.score << "  "
                         << spp.residue_prev << " " << spp.residue_this << " " << spp.residue_next << std::endl;
            if (! spp.residue_this) continue;

            std::string lab = "Ramachandran ";
            lab += spp.res_spec.chain_id;
            lab += " ";
            lab += std::to_string(spp.res_spec.res_no);
            lab += " ";
            lab += std::string(spp.residue_this->GetResName());
            lab += " pr: ";
            lab += coot::util::float_to_string_using_dec_pl(spp.score, 3);
            GtkWidget *button = gtk_button_new();
            GtkWidget *label = gtk_label_new(lab.c_str());
            gtk_widget_set_halign(label, GTK_ALIGN_START);
            gtk_button_set_child(GTK_BUTTON(button), label);
            gtk_widget_set_margin_start (button, 4);
            gtk_widget_set_margin_end   (button, 4);
            gtk_widget_set_margin_top   (button, 2);
            gtk_widget_set_margin_bottom(button, 2);

            coot::atom_spec_t atom_spec(spp.res_spec.chain_id, spp.res_spec.res_no, spp.res_spec.ins_code, " CA ", "");
            coot::atom_spec_t *spec_p = new coot::atom_spec_t(atom_spec);
            spec_p->int_user_data = imol;
            g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(rama_button_clicked_callback), spec_p);
            auto p = std::make_pair(spp.res_spec, button);
            buttons.push_back(p);
         }
      }
      return buttons;
   };

   auto non_pro_cis_button_clicked_callback = +[] (GtkButton *button, gpointer user_data) {
      coot::atom_spec_t *atom_spec_p = static_cast<coot::atom_spec_t *>(user_data);
      int imol = atom_spec_p->int_user_data;
      const auto &atom_spec(*atom_spec_p);
      graphics_info_t g;
      graphics_info_t::set_go_to_atom(imol, atom_spec);
      g.try_centre_from_new_go_to_atom(); // oh dear!
      g.graphics_grab_focus();
   };

   auto make_non_pro_cis_peptide_buttons = [non_pro_cis_button_clicked_callback] (int imol, mmdb::Manager *mol) {
      std::vector<std::pair<coot::residue_spec_t, GtkWidget *> > buttons;

      std::vector<coot::util::cis_peptide_info_t> v = coot::util::cis_peptides_info_from_coords(mol);
      for (unsigned int i=0; i<v.size(); i++) {
         const auto &cpi = v[i];
         coot::residue_spec_t res_spec(cpi.chain_id_1, cpi.resno_1, cpi.ins_code_1);
         mmdb::Residue *residue_p = coot::util::get_residue(res_spec, mol);
         if (residue_p) {
            if (cpi.residue_name_2 != "PRO") {
               std::string lab = "Non-PRO <i>cis</i> Peptide ";
               lab += cpi.chain_id_1;
               lab += " ";
               lab += std::to_string(cpi.resno_1);
               lab += " ω: ";
               lab += coot::util::float_to_string_using_dec_pl(cpi.omega_torsion_angle, 2);
               lab += "°";
               // GtkWidget *button = gtk_button_new_with_label(lab.c_str());
               GtkWidget *button = gtk_button_new();
               GtkWidget *label = gtk_label_new(lab.c_str());
               gtk_label_set_use_markup(GTK_LABEL(label), TRUE);
               gtk_widget_set_halign(label, GTK_ALIGN_START);
               gtk_button_set_child(GTK_BUTTON(button), label);
               gtk_widget_set_margin_start (button, 4);
               gtk_widget_set_margin_end   (button, 4);
               gtk_widget_set_margin_top   (button, 2);
               gtk_widget_set_margin_bottom(button, 2);
               coot::atom_spec_t spec(cpi.chain_id_1, cpi.resno_1, cpi.ins_code_1, " CA ", "");
               coot::atom_spec_t *spec_p = new coot::atom_spec_t(spec);
               spec_p->int_user_data = imol;
               g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(non_pro_cis_button_clicked_callback), spec_p);
               auto p = std::make_pair(res_spec, button);
               buttons.push_back(p);
            }
         }
      }
      return buttons;
   };

   auto twisted_trans_button_clicked_callback = +[] (GtkButton *button, gpointer user_data) {
      coot::atom_spec_t *atom_spec_p = static_cast<coot::atom_spec_t *>(user_data);
      int imol = atom_spec_p->int_user_data;
      const auto &atom_spec(*atom_spec_p);
      graphics_info_t g;
      graphics_info_t::set_go_to_atom(imol, atom_spec);
      g.try_centre_from_new_go_to_atom(); // oh dear!
      g.graphics_grab_focus();
   };

   auto make_twisted_trans_buttons = [twisted_trans_button_clicked_callback] (int imol, mmdb::Manager *mol) {
      std::vector<std::pair<coot::residue_spec_t, GtkWidget *> > buttons;

      graphics_info_t g;
      std::vector<coot::util::cis_peptide_quad_info_t> v = coot::cis_peptide_quads_from_coords(mol, 0, g.Geom_p(), false);
      for (unsigned int i=0; i<v.size(); i++) {
	 if (v[i].type == coot::util::cis_peptide_quad_info_t::TWISTED_TRANS) {
            try {
	       coot::residue_spec_t r1(v[i].quad.atom_1->GetResidue());
               coot::residue_spec_t r2(v[i].quad.atom_4->GetResidue());
               float omega = v[i].quad.torsion();
               std::string lab = "Twisted <i>trans</i> ";
               lab += std::string(v[i].quad.atom_2->GetChainID());
               lab += " ";
               lab += std::to_string(v[i].quad.atom_2->GetSeqNum());
               lab += " ω: ";
               lab += coot::util::float_to_string_using_dec_pl(omega, 2);
               lab += "°";
               // GtkWidget *button = gtk_button_new_with_label(lab.c_str());
               GtkWidget *button = gtk_button_new();
               GtkWidget *label = gtk_label_new(lab.c_str());
               gtk_label_set_use_markup(GTK_LABEL(label), TRUE);
               gtk_widget_set_halign(label, GTK_ALIGN_START);
               gtk_button_set_child(GTK_BUTTON(button), label);
               gtk_widget_set_margin_start (button, 4);
               gtk_widget_set_margin_end   (button, 4);
               gtk_widget_set_margin_top   (button, 2);
               gtk_widget_set_margin_bottom(button, 2);
               coot::residue_spec_t res_spec(v[i].quad.atom_2->residue);
               coot::atom_spec_t *spec_p = new coot::atom_spec_t(v[i].quad.atom_2);
               spec_p->int_user_data = imol;
               g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(twisted_trans_button_clicked_callback), spec_p);
               auto p = std::make_pair(res_spec, button);
               buttons.push_back(p);
            }
            catch (const std::runtime_error &e) {
               std::cout << "WARNING::" << e.what() << std::endl;
            }
         }
      }
      return buttons;
   };

   auto c_beta_devi_button_clicked_callback = +[] (GtkButton *button, gpointer user_data) {
      coot::atom_spec_t *atom_spec_p = static_cast<coot::atom_spec_t *>(user_data);
      int imol = atom_spec_p->int_user_data;
      const auto &atom_spec(*atom_spec_p);
      graphics_info_t g;
      graphics_info_t::set_go_to_atom(imol, atom_spec);
      g.try_centre_from_new_go_to_atom(); // oh dear!
      g.graphics_grab_focus();
   };

   auto make_cbeta_devi_buttons = [c_beta_devi_button_clicked_callback] (int imol, mmdb::Manager *mol) {
      double dist_crit = 0.25;
      std::vector<std::pair<coot::residue_spec_t, GtkWidget *> > buttons;
      // the key of the map is the alt conf.
      std::map<mmdb::Residue *, std::map<std::string, coot::c_beta_deviation_t> >
         residue_c_beta_map = coot::get_c_beta_deviations(mol);
      std::map<mmdb::Residue *, std::map<std::string, coot::c_beta_deviation_t> >::const_iterator it;
      for (it=residue_c_beta_map.begin(); it!=residue_c_beta_map.end(); ++it) {
         // multiple alt confs
         mmdb::Residue *res_key = it->first;
         const std::map<std::string, coot::c_beta_deviation_t> &value_map = it->second;
         std::map<std::string, coot::c_beta_deviation_t>::const_iterator it_inner;
         for (it_inner=value_map.begin(); it_inner!=value_map.end(); ++it_inner) {
            const std::string alt_conf_key = it_inner->first;
            const coot::c_beta_deviation_t &cbd = it_inner->second;
            if (cbd.dist > dist_crit) {
               coot::residue_spec_t res_spec(res_key);
               std::string lab = "Cβ Deviation ";
               lab += " ";
               lab += cbd.at->GetChainID();
               lab += " ";
               lab += std::to_string(cbd.at->GetSeqNum());
               lab += " ";
               lab += std::string(cbd.at->residue->GetResName());
               lab += " d: ";
               lab += coot::util::float_to_string_using_dec_pl(cbd.dist, 2);
               if (! alt_conf_key.empty())
                  lab += std::string(" ") + alt_conf_key;
               GtkWidget *button = gtk_button_new();
               GtkWidget *label = gtk_label_new(lab.c_str());
               gtk_widget_set_halign(label, GTK_ALIGN_START);
               gtk_button_set_child(GTK_BUTTON(button), label);
               gtk_widget_set_margin_start (button, 4);
               gtk_widget_set_margin_end   (button, 4);
               gtk_widget_set_margin_top   (button, 2);
               gtk_widget_set_margin_bottom(button, 2);
               coot::atom_spec_t spec(cbd.at);
               coot::atom_spec_t *spec_p = new coot::atom_spec_t(spec);
               spec_p->int_user_data = imol;
               g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(c_beta_devi_button_clicked_callback), spec_p);
               auto p = std::make_pair(res_spec, button);
               buttons.push_back(p);
            }
         }
      }
      return buttons;
   };

   auto make_template_baddies = [chiral_volume_button_clicked_callback] (mmdb::Manager *mol) {
      std::vector<std::pair<coot::residue_spec_t, GtkWidget *> > buttons;
      std::vector<coot::util::cis_peptide_info_t> v; // = something
      for (unsigned int i=0; i<v.size(); i++) {
         const auto &cpi = v[i];
         mmdb::Residue *residue_p = nullptr;
         std::string lab = "xx";
         lab += " ";
         GtkWidget *button = gtk_button_new();
         GtkWidget *label = gtk_label_new(lab.c_str());
         gtk_widget_set_halign(label, GTK_ALIGN_START);
         gtk_button_set_child(GTK_BUTTON(button), label);
         gtk_widget_set_margin_start (button, 4);
         gtk_widget_set_margin_end   (button, 4);
         gtk_widget_set_margin_top   (button, 2);
         gtk_widget_set_margin_bottom(button, 2);
         coot::residue_spec_t res_spec(residue_p);
         // coot::atom_spec_t *spec_p = new coot::atom_spec_t(spec);
         // spec_p->int_user_data = imol;
         g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(chiral_volume_button_clicked_callback), nullptr);
         auto p = std::make_pair(res_spec, button);
         buttons.push_back(p);
      }
      return buttons;
   };

   auto sorter = [] (const std::pair<coot::residue_spec_t, GtkWidget *> &r1,
                     const std::pair<coot::residue_spec_t, GtkWidget *> &r2) {
      return (r1.first<r2.first);
   };

   auto sort_buttons = [sorter] (const std::vector<std::pair<coot::residue_spec_t, GtkWidget *> > &buttons_in) {
      auto buttons = buttons_in;
      std::sort(buttons.begin(), buttons.end(), sorter);
      return buttons;
   };

   int pos = gtk_paned_get_position(GTK_PANED(pane));
   // // std::cout << "here in overlaps_peptides_cbeta_ramas_and_rotas_internal(): with pos " << pos << std::endl;
   if (pos < 320)
      gtk_paned_set_position(GTK_PANED(pane), 320);

   graphics_info_t::clear_out_container(vbox);

   mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;

   std::vector<std::pair<coot::residue_spec_t, GtkWidget *> > buttons;
   auto overlap_buttons = make_overlap_buttons(mol);
   buttons.insert(buttons.end(), overlap_buttons.begin(), overlap_buttons.end());
   auto rota_buttons = make_rota_buttons(mol);
   buttons.insert(buttons.end(), rota_buttons.begin(), rota_buttons.end());
   auto chiral_buttons = make_chiral_volume_buttons(imol);
   buttons.insert(buttons.end(), chiral_buttons.begin(), chiral_buttons.end());
   auto rama_buttons = make_rama_buttons(imol);
   buttons.insert(buttons.end(), rama_buttons.begin(), rama_buttons.end());
   auto npcp_buttons = make_non_pro_cis_peptide_buttons(imol, mol);
   buttons.insert(buttons.end(), npcp_buttons.begin(), npcp_buttons.end());
   auto twisted_trans_buttons = make_twisted_trans_buttons(imol, mol);
   buttons.insert(buttons.end(), twisted_trans_buttons.begin(), twisted_trans_buttons.end());
   auto cbeta_devi_buttons = make_cbeta_devi_buttons(imol, mol);
   buttons.insert(buttons.end(), cbeta_devi_buttons.begin(), cbeta_devi_buttons.end());

   if (label) {
      unsigned int n_baddies = buttons.size();
      std::string label_txt = "Outlier Count: ";
      label_txt += std::to_string(n_baddies);
      gtk_label_set_text(GTK_LABEL(label), label_txt.c_str());
   }

   auto sorted_buttons = sort_buttons(buttons); // sort by chain-id and residue number

   for(auto button_info : sorted_buttons)
      gtk_box_append(GTK_BOX(vbox), button_info.second);

}

void phenix_geo_validation_buttons(int imol,
                                   const coot::phenix_geo::phenix_geometry &pg,
                                   double residual_cutoff) {

   auto atom_spec_to_label = [] (const coot::atom_spec_t &spec) {
      std::string s = spec.chain_id;
      s += " ";
      s += std::to_string(spec.res_no);
      if (!spec.ins_code.empty()) s += " " + spec.ins_code;
      s += " ";
      s += coot::util::remove_whitespace(spec.atom_name);
      if (!spec.alt_conf.empty()) s += " " + spec.alt_conf;
      return s;
   };

   auto button_callback = +[] (GtkButton *button, G_GNUC_UNUSED gpointer user_data) {
      std::string target_position_x_str(static_cast<const char *>(g_object_get_data(G_OBJECT(button), "target-position-x")));
      std::string target_position_y_str(static_cast<const char *>(g_object_get_data(G_OBJECT(button), "target-position-y")));
      std::string target_position_z_str(static_cast<const char *>(g_object_get_data(G_OBJECT(button), "target-position-z")));
      try {
         float x = coot::util::string_to_float(target_position_x_str);
         float y = coot::util::string_to_float(target_position_y_str);
         float z = coot::util::string_to_float(target_position_z_str);
         clipper::Coord_orth p(x, y, z);
         graphics_info_t::set_rotation_centre(p);
      }
      catch (const std::runtime_error &e) {
         std::cout << "WARNING::" << e.what() << std::endl;
      }
      graphics_info_t g;
      g.graphics_grab_focus();
   };

   auto set_target_position_data = [] (GtkWidget *button, const clipper::Coord_orth &p) {
      char *x = new char[10];
      char *y = new char[10];
      char *z = new char[10];
      for (unsigned int i=0; i<10; i++) x[i] = 0;
      for (unsigned int i=0; i<10; i++) y[i] = 0;
      for (unsigned int i=0; i<10; i++) z[i] = 0;
      strncpy(x, std::to_string(p.x()).c_str(), 9);
      strncpy(y, std::to_string(p.y()).c_str(), 9);
      strncpy(z, std::to_string(p.z()).c_str(), 9);
      g_object_set_data(G_OBJECT(button), "target-position-x", x);
      g_object_set_data(G_OBJECT(button), "target-position-y", y);
      g_object_set_data(G_OBJECT(button), "target-position-z", z);
   };

   auto sorter = [] (const std::pair<coot::residue_spec_t, GtkWidget *> &r1,
                     const std::pair<coot::residue_spec_t, GtkWidget *> &r2) {
      return (r1.first < r2.first);
   };

   GtkWidget *pane_to_show = widget_from_builder("main_window_ramchandran_and_validation_pane");
   gtk_widget_set_visible(pane_to_show, TRUE);
   GtkWidget *pane = widget_from_builder("main_window_graphics_rama_vs_graphics_pane");

   GtkWidget *vbox_2 = widget_from_builder("validation_boxes_vbox");
   gtk_widget_set_visible(vbox_2, TRUE);

   GtkWidget *outer_vbox = widget_from_builder("dynamic_validation_vbox");
   gtk_widget_set_visible(outer_vbox, TRUE);

   GtkWidget *vbox  = widget_from_builder("dynamic_validation_outliers_vbox");
   GtkWidget *label = widget_from_builder("dynamic_validation_outliers_label");

   int pos = gtk_paned_get_position(GTK_PANED(pane));
   if (pos < 320)
      gtk_paned_set_position(GTK_PANED(pane), 320);

   graphics_info_t::clear_out_container(vbox);

   mmdb::Manager *mol = graphics_info_t::molecules[imol].atom_sel.mol;

   std::vector<std::pair<coot::residue_spec_t, GtkWidget *> > buttons;

   // Bond outliers
   for (unsigned int i = 0; i < pg.geo_bonds.size(); i++) {
      const coot::phenix_geo::phenix_geo_bond &gb = pg.geo_bonds[i];
      if (gb.residual < residual_cutoff) continue;

      mmdb::Atom *at_1 = coot::util::get_atom_using_fuzzy_search(gb.atom_1, mol);
      mmdb::Atom *at_2 = coot::util::get_atom_using_fuzzy_search(gb.atom_2, mol);
      if (!at_1 || !at_2) continue;

      clipper::Coord_orth p1(at_1->x, at_1->y, at_1->z);
      clipper::Coord_orth p2(at_2->x, at_2->y, at_2->z);
      clipper::Coord_orth midpoint = 0.5 * (p1 + p2);

      std::string lab = "Bond " + atom_spec_to_label(gb.atom_1) + " - " + atom_spec_to_label(gb.atom_2);
      lab += " Δ: " + coot::util::float_to_string_using_dec_pl(gb.delta, 3) + "Å";
      lab += " (" + coot::util::float_to_string_using_dec_pl(std::sqrt(gb.residual), 1) + "σ)";

      GtkWidget *button = gtk_button_new();
      GtkWidget *button_label = gtk_label_new(lab.c_str());
      gtk_widget_set_halign(button_label, GTK_ALIGN_START);
      gtk_button_set_child(GTK_BUTTON(button), button_label);
      gtk_widget_set_margin_start (button, 4);
      gtk_widget_set_margin_end   (button, 4);
      gtk_widget_set_margin_top   (button, 2);
      gtk_widget_set_margin_bottom(button, 2);

      set_target_position_data(button, midpoint);
      g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(button_callback), nullptr);

      coot::residue_spec_t res_spec(gb.atom_1);
      buttons.push_back(std::make_pair(res_spec, button));
   }


   // Angle outliers
   for (unsigned int i = 0; i < pg.geo_angles.size(); i++) {
      const coot::phenix_geo::phenix_geo_angle &ga = pg.geo_angles.angles[i];
      if (ga.residual < residual_cutoff) continue;

      mmdb::Atom *at_2 = coot::util::get_atom_using_fuzzy_search(ga.atom_2, mol); // middle atom of angle
      if (!at_2) {
         std::cout << "DEBUG:: Angle outliers: failed to get ga.atom_2 " << ga.atom_2 << std::endl;
         continue;
      }

      clipper::Coord_orth target(at_2->x, at_2->y, at_2->z);

      std::string lab = "Angle " + atom_spec_to_label(ga.atom_1) + " - ";
      lab += atom_spec_to_label(ga.atom_2) + " - " + atom_spec_to_label(ga.atom_3);
      lab += " Δ: " + coot::util::float_to_string_using_dec_pl(ga.delta, 1) + "°";
      lab += " (" + coot::util::float_to_string_using_dec_pl(std::sqrt(ga.residual), 1) + "σ)";

      GtkWidget *button = gtk_button_new();
      GtkWidget *button_label = gtk_label_new(lab.c_str());
      gtk_widget_set_halign(button_label, GTK_ALIGN_START);
      gtk_button_set_child(GTK_BUTTON(button), button_label);
      gtk_widget_set_margin_start (button, 4);
      gtk_widget_set_margin_end   (button, 4);
      gtk_widget_set_margin_top   (button, 2);
      gtk_widget_set_margin_bottom(button, 2);

      set_target_position_data(button, target);
      g_signal_connect(G_OBJECT(button), "clicked", G_CALLBACK(button_callback), nullptr);

      coot::residue_spec_t res_spec(ga.atom_2);
      buttons.push_back(std::make_pair(res_spec, button));
   }


   if (label) {
      unsigned int n_outliers = buttons.size();
      std::string label_txt = "Phenix Geo Outliers: ";
      label_txt += std::to_string(n_outliers);
      gtk_label_set_text(GTK_LABEL(label), label_txt.c_str());
   }

   std::sort(buttons.begin(), buttons.end(), sorter);

   for (const auto &button_info : buttons) {
      gtk_widget_set_visible(button_info.second, TRUE);
      gtk_box_append(GTK_BOX(vbox), button_info.second);
   }
   gtk_widget_set_visible(vbox, TRUE);
}
