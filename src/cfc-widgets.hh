/* src/cfc-widgets.cc
 * 
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
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 */


#ifndef CFC_WIDGETS_HH
#define CFC_WIDGETS_HH

#ifdef MAKE_ENHANCED_LIGAND_TOOLS

#include <gtk/gtk.h>

#include "cfc.hh"

namespace cfc {
   
   GtkWidget *wrapped_create_cfc_dialog(extracted_cluster_info_from_python extracted_cluster_info);
   void on_cfc_site_button_clicked(GtkButton *button, gpointer user_data);
   void on_cfc_water_cluster_button_clicked(GtkButton *button, gpointer user_data);
   void on_cfc_water_cluster_structure_button_clicked(GtkButton *button, gpointer user_data);

   void on_cfc_pharmacophore_cluster_button_clicked(GtkButton *button, gpointer user_data);
   void on_cfc_pharmacophore_cluster_structure_button_clicked(GtkButton *button, gpointer user_data);
   void cfc_table_show_hide(std::string show_this_one_name, GtkWidget *vbox);
   

   // pass by value because we need to use map[] operator (or do we?)
   // 
   void cfc_dialog_add_waters(unsigned int site_number,
			      extracted_cluster_info_from_python extracted_cluster_info,
			      GtkWidget *cfc_dialog);
   void cfc_dialog_add_pharmacophores(unsigned int site_number,
				      extracted_cluster_info_from_python extracted_cluster_info,
				      GtkWidget *cfc_dialog);

   class pharm_button_set {
   public:
      GtkWidget *left_hand_button;
      GtkWidget *structure_buttons_hbox;
      double conservedness;
      std::string type;
      int type_rank;
      pharm_button_set(std::string type_in,
		       GtkWidget *left_hand_button_in,
		       double conservedness_in) {
	 left_hand_button = left_hand_button_in;
	 conservedness = conservedness_in;
	 type = type_in;
	 type_rank = 6;
	 if (type == "Donor")
	    type_rank = 1;
	 if (type == "Acceptor")
	    type_rank = 2;
	 if (type == "Aromatic")
	    type_rank = 3;
	 if (type == "Hydrophobic")
	    type_rank = 4;
      }
      static bool sorter(const pharm_button_set &p1,
			 const pharm_button_set &p2) {
	 if (p2.conservedness < p1.conservedness) {
	    return true;
	 } else {
	    if (p2.conservedness > p1.conservedness) {
	       return false;
	    } else {
	       return (p1.type_rank < p2.type_rank);
	    }
	 }
      }
      void add_structure_buttons_hbox(GtkWidget *b) {
	 structure_buttons_hbox = b;
      }
   };

   // on clicked, go a cluster mean, show all structures that
   // contribute to this cluster, and the chemical features of
   // those ligands that match the type of the cluster.
   // 
   // We need cluster mean, and the contributors to this cluster.
   //
   class on_pharmacophore_click_info_t {
   public:
      clipper::Coord_orth pos;
      std::vector<std::pair<int, coot::residue_spec_t> > imol_residue_specs;
      on_pharmacophore_click_info_t(const clipper::Coord_orth &pos_in,
				    const std::vector<std::pair<int, coot::residue_spec_t> > &imol_residue_specs_in) {
	 pos = pos_in;
	 imol_residue_specs = imol_residue_specs_in;
      }
      on_pharmacophore_click_info_t(clipper::Coord_orth &pos_in) {
	 pos = pos_in;
      }
      void add_imol_spec(int imol, const coot::residue_spec_t &r) {
	 std::pair<int, coot::residue_spec_t> p(imol, r);
	 imol_residue_specs.push_back(p);
      }
   };

   // on clicked, go to this pharmacophore position,
   // show only this structure and this ligand and the pharmacophores
   // for this ligand. We need imol, res_spec, and pharmacophore pos.
   //
   class on_pharmacophore_structure_click_info_t {
   public:
      clipper::Coord_orth pos;
      int imol;
      coot::residue_spec_t res_spec;
      on_pharmacophore_structure_click_info_t(const clipper::Coord_orth &pos_in,
					      int imol_in,
					      const coot::residue_spec_t &spec_in) {
	 pos = pos_in;
	 imol = imol_in;
	 res_spec = spec_in;
      }
   };

   void chemical_feature_clusters_add_site_info(unsigned int site_number,
						const extracted_cluster_info_from_python &eci,
						GtkWidget *cfc_dialog);

   void cfc_dialog_add_site_info(unsigned int site_number,
				 const extracted_cluster_info_from_python &eci);
   

}

#endif // MAKE_ENHANCED_LIGAND_TOOLS

#endif // CFC_WIDGETS_HH

