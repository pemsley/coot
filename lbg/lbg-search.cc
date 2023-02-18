/* lbg/lbg-search.cc
 * 
 * Copyright 2010, 2011, 2012 by The University of Oxford
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

#ifdef EMSCRIPTEN_THING

#include <fstream>

#include <stdlib.h> // for getenv()

#include "clipper/core/coords.h"

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include "geometry/srs-interface.hh"
#ifdef RDKIT_HAS_CAIRO_SUPPORT
#include <cairo.h>
// used to be installed in the wrong directory
#include <GraphMol/MolDraw2D/MolDraw2DCairo.h>
#include "lidia-core/rdkit-interface.hh"
#else
#include "lidia-core/rdkit-interface.hh"
#endif // RDKIT_HAS_CAIRO_SUPPORT
#endif // MAKE_ENHANCED_LIGAND_TOOLS

#include "lbg.hh"

mmdb::math::Graph *
get_graph_2() {

   mmdb::math::Graph *graph = new mmdb::math::Graph;

   std::vector<std::pair<std::string, int> > v;
   v.push_back(std::pair<std::string, int>("N",   7));
   v.push_back(std::pair<std::string, int>("CA",  6));
   v.push_back(std::pair<std::string, int>("CB",  6));
   v.push_back(std::pair<std::string, int>("CG",  6));
   v.push_back(std::pair<std::string, int>("CD1", 6));
   v.push_back(std::pair<std::string, int>("NE1", 7));
   v.push_back(std::pair<std::string, int>("CE2", 6));
   v.push_back(std::pair<std::string, int>("CD2", 6));
   v.push_back(std::pair<std::string, int>("CE3", 6));
   v.push_back(std::pair<std::string, int>("CZ3", 6));
   v.push_back(std::pair<std::string, int>("CH2", 6));
   v.push_back(std::pair<std::string, int>("CZ2", 6));
   v.push_back(std::pair<std::string, int>("C",   6));
   v.push_back(std::pair<std::string, int>("O",   8));

   // mmdb::math::Vertex *vert = new mmdb::math::Vertex("", "");
   // graph->AddVertex(vert);
   for (unsigned int i=0; i<v.size(); i++) {
      std::string ele = "C";
      if (v[i].second == 7) ele = "N";
      if (v[i].second == 8) ele = "0";
      std::string name = v[i].first;
      mmdb::math::Vertex *vert = new mmdb::math::Vertex(ele.c_str(), name.c_str());
      graph->AddVertex(vert);
   }

   std::vector<std::pair<std::pair<int, int>, int> > e;
   e.push_back(std::pair<std::pair<int, int>, int> (std::pair<int,int>(1,  2), 1));
   e.push_back(std::pair<std::pair<int, int>, int> (std::pair<int,int>(2,  3), 1));
   e.push_back(std::pair<std::pair<int, int>, int> (std::pair<int,int>(3,  4), 1));
   e.push_back(std::pair<std::pair<int, int>, int> (std::pair<int,int>(4,  8), 1));
   e.push_back(std::pair<std::pair<int, int>, int> (std::pair<int,int>(4,  5), 2));
   e.push_back(std::pair<std::pair<int, int>, int> (std::pair<int,int>(5,  6), 1));
   e.push_back(std::pair<std::pair<int, int>, int> (std::pair<int,int>(6,  7), 1));
   e.push_back(std::pair<std::pair<int, int>, int> (std::pair<int,int>(7, 12), 1));
   e.push_back(std::pair<std::pair<int, int>, int> (std::pair<int,int>(7,  8), 2));
   e.push_back(std::pair<std::pair<int, int>, int> (std::pair<int,int>(8,  9), 1));
   e.push_back(std::pair<std::pair<int, int>, int> (std::pair<int,int>(9, 10), 2));
   e.push_back(std::pair<std::pair<int, int>, int> (std::pair<int,int>(10, 11), 1));
   e.push_back(std::pair<std::pair<int, int>, int> (std::pair<int,int>(11, 12), 2));
   e.push_back(std::pair<std::pair<int, int>, int> (std::pair<int,int>( 2, 13), 1));
   e.push_back(std::pair<std::pair<int, int>, int> (std::pair<int,int>(13, 14), 1));

   for (unsigned int i=0; i<e.size(); i++) {
      mmdb::math::Edge *ed = new mmdb::math::Edge(e[i].first.first, e[i].first.second,
						  e[i].second);
      graph->AddEdge(ed);
   }

   return graph;
}


#ifdef HAVE_CCP4SRS
// not const because try_dynamic_add() can be called
void
lbg_info_t::search() {

   double local_search_similarity = get_search_similarity();

   mmdb::Residue *res = 0;
   
   mmdb::math::Graph *graph = new mmdb::math::Graph;
   int n_atoms = 0;

   // mol atom indexing -> graph vertex indexing
   int vertex_indexing[mol.atoms.size()];

   // atoms are vertices
   for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {
      if (! mol.atoms[iat].is_closed()) {
	 std::string ele = mol.atoms[iat].element;
	 std::string name = mol.atoms[iat].get_atom_id();
	 // std::cout << "adding atom with ele " << ele << " and name \"" << name << "\"" << std::endl;
	 mmdb::math::Vertex *v = new mmdb::math::Vertex(ele.c_str(), name.c_str());
	 vertex_indexing[n_atoms] = iat;
	 v->SetUserID(iat);
	 graph->AddVertex(v);
	 n_atoms++;
      }
   }

   // bonds are edges
   int n_bonds = 0;
   for (unsigned int ib=0; ib<mol.bonds.size(); ib++) {
      if (! mol.bonds[ib].is_closed()) {
	 int mmdb_bond_type = mol.bonds[ib].mmdb_bond_type();
	 if (mmdb_bond_type != UNASSIGNED_INDEX) {
	    int ind_1 = mol.bonds[ib].get_atom_1_index();
	    int ind_2 = mol.bonds[ib].get_atom_2_index();
 	    mmdb::math::Edge *e = new mmdb::math::Edge(vertex_indexing[ind_1]+1,  // 1-indexed
						       vertex_indexing[ind_2]+1,
						       mmdb_bond_type);
	    graph->AddEdge(e);
	    n_bonds++;
	 }
      }
   }

   // hack for testing
   // graph = get_graph_2();
   // n_atoms = 14;
   
   graph->SetName ("Coot-LBG-Query");
   graph->MakeVertexIDs();
   int build_result = graph->Build(true);

   if (build_result != 0) {

      std::cout << "Bad graph build result" << std::endl;

   } else {

      if (stand_alone_flag) {

	 // non-const protein_geometry *
	 coot::protein_geometry *geom_p_local = new coot::protein_geometry;

	 const char *d1 = getenv(MONOMER_DIR_STR); // "COOT_CCP4SRS_DIR"
	 std::string srs_dir = PKGDATADIR;
	 if (d1)
	    srs_dir = d1;
	 
	 std::cout << "------------ geom_p_local init_ccp4srs with srs_dir " << srs_dir << std::endl;

	 geom_p_local->init_ccp4srs(srs_dir);

	 graph->MakeSymmetryRelief(true);
	 graph->Print();
	 std::cout << "graph search using similarity  " << local_search_similarity << std::endl;
	 std::cout << "graph build returns: " << build_result << std::endl;
	 bool fill_graph_match = true;
	 std::vector<coot::match_results_t> v =
	    geom_p_local->compare_vs_ccp4srs(graph, local_search_similarity, n_atoms,
					     -1, -1, fill_graph_match);
	 delete graph;
	 std::cout << "found " << v.size() << " close matches" << std::endl;
	 display_search_results(v);

      } else {

	 if (geom_p) { 
	    graph->MakeSymmetryRelief(true);
	    graph->Print();
	    std::cout << "graph search using similarity  " << local_search_similarity << std::endl;
	    std::cout << "graph build returns: " << build_result << std::endl;
	    bool fill_graph_match = true;
	    std::vector<coot::match_results_t> v =
	       geom_p->compare_vs_ccp4srs(graph, local_search_similarity, n_atoms, -1, -1,
					  fill_graph_match);
	    delete graph;
	    display_search_results(v);
	 } else {
	    std::cout << "WARNING:: No geometry library for ligand search() " << std::endl;
	 }
      }
   }
}


// allow access of the Search button callback to the search
// similarity combox box text.
double
lbg_info_t::get_search_similarity() const {

   double r = search_similarity;

   gchar *txt = gtk_combo_box_text_get_active_text(GTK_COMBO_BOX_TEXT(lbg_search_combobox));
   if (txt) {
      try { 
	 r = lig_build::string_to_float(txt);
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "WARNING:: Bad number " << txt << " " << rte.what() << std::endl;
	 std::cout << "WARNING:: Using " << r << std::endl;
      }
   } 
   return r;
}
#endif



#ifdef HAVE_CCP4SRS
coot::match_results_t
lbg_info_t::residue_from_best_match(mmdb::math::Graph &graph_1, mmdb::math::Graph &graph_2,
				    mmdb::math::GraphMatch &match, int n_match, 
				    ccp4srs::Monomer *monomer_p) const {

   coot::match_results_t r("","",NULL);
   if (geom_p)
      r = geom_p->residue_from_best_match(graph_1, graph_2, match, n_match, monomer_p);
   return r;
}
#endif // HAVE_CCP4SRS

// not const because try_dynamic_add() can be called (to make images)
//
void
lbg_info_t::display_search_results(const std::vector<coot::match_results_t> &v) {

   bool new_dialog = false;

   if (new_dialog) {

      gtk_widget_show(lbg_sbase_search_results_dialog);
      GtkSettings *default_settings = gtk_settings_get_default();
      g_object_set(default_settings, "gtk-button-images", TRUE, NULL);

      // clear GTK_BOX(lbg_sbase_search_results_vbox) here
      // 
      GList* glist = gtk_container_get_children(GTK_CONTAINER(lbg_sbase_search_results_vbox));
      while (glist) {
	 GtkWidget *w = GTK_WIDGET(glist->data);
	 gtk_widget_destroy(w); // clear existing children
	 glist = glist->next;
      }

      for (unsigned int i=0; i<v.size(); i++) {
	 std::string lab = v[i].comp_id;
	 lab += ":  ";
	 lab += v[i].name;

            // monomer search does it like this (it has left-aligned text)
            //
            // GtkWidget *button = gtk_button_new();
            // GtkWidget *label  = gtk_label_new(l.c_str());
            // GtkWidget *button_hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
            // gtk_container_add(GTK_CONTAINER(button), button_hbox);
            // gtk_box_pack_start(GTK_BOX(button_hbox), label, FALSE, FALSE, 0);

	 GtkWidget *button = gtk_button_new_with_label(lab.c_str());
	 gtk_box_pack_start(GTK_BOX(lbg_sbase_search_results_vbox),
			    GTK_WIDGET(button), FALSE, FALSE, 3);
	 std::string *comp_id = new std::string(v[i].comp_id);
	 g_signal_connect(GTK_WIDGET(button), "clicked",
			  G_CALLBACK(on_sbase_search_result_button_clicked),
			  (gpointer) (comp_id));

	 // gtk_button_set_alignment(GTK_BUTTON(button), 0, 0.5); // deprecated  - but how to FIX-IT?

	 g_object_set_data(G_OBJECT(button), "lbg", (gpointer) this);
	 gtk_widget_show(button);
      }
   } else {

      // put the results into lbg_srs_search_results_scrolledwindow
      gtk_widget_show(lbg_srs_search_results_scrolledwindow);
      if (lbg_srs_search_results_vbox) {

#ifdef HAVE_CCP4SRS
	 std::string srs_dir = coot::get_srs_dir(); // use environment variables or ccp4/prefix-dir fall-back
	 ccp4srs::Manager *srs_manager = new ccp4srs::Manager;
	 srs_manager->loadIndex(srs_dir.c_str());
#endif // HAVE_CCP4SRS
	 gtk_widget_show(lbg_srs_search_results_vbox);

	 for (unsigned int i=0; i<v.size(); i++) {
	    std::string lab = v[i].comp_id;
	    lab += ":  ";
	    lab += v[i].name;
	    // GtkWidget *button = gtk_button_new_with_label(lab.c_str());
	    GtkWidget *button = gtk_button_new();
	    gtk_box_pack_start(GTK_BOX(lbg_srs_search_results_vbox),
			       GTK_WIDGET(button), FALSE, FALSE, 3);
	    std::string *comp_id = new std::string(v[i].comp_id);
	    g_signal_connect(GTK_WIDGET(button), "clicked",
			     G_CALLBACK(on_sbase_search_result_button_clicked),
			     (gpointer) (comp_id));
	    gtk_button_set_alignment(GTK_BUTTON(button), 0, 0.5);
	    g_object_set_data(G_OBJECT(button), "lbg", (gpointer) this);

	    GtkWidget *label  = gtk_label_new(lab.c_str());
	    // GtkWidget *button_hbox = gtk_hbox_new(FALSE, 0);
            GtkWidget *button_hbox = gtk_box_new(GTK_ORIENTATION_HORIZONTAL, 0);
	    gtk_container_add(GTK_CONTAINER(button), button_hbox);
	    int imol = 0; // dummy

#ifdef HAVE_CCP4SRS
#ifdef RDKIT_HAS_CAIRO_SUPPORT
	    GtkWidget *wp = get_image_widget_for_comp_id(v[i].comp_id, imol, srs_manager);
	    if (wp) {
	       gtk_widget_show(wp);
	       // std::cout << "adding image " << wp << std::endl;
	       // std::cout << "gtk_box_pack_start() " << button_hbox << " " << wp << std::endl;
	       gtk_box_pack_start(GTK_BOX(button_hbox), wp, FALSE, FALSE, 0);
	    } else {
	       std::cout << "Null image for " << v[i].comp_id << std::endl;
	    }
#endif
#endif
	    gtk_box_pack_start(GTK_BOX(button_hbox), label, FALSE, FALSE, 0);
	    gtk_widget_show(label);
	    gtk_widget_show(button);
	    gtk_widget_show(button_hbox);
	 }
      } else {
	 std::cout << "ERROR:: Null lbg_srs_search_results_vbox" << std::endl;
      }
   }
}


#ifdef HAVE_CCP4SRS
// because now we don't construct rdkit molecules by looking up the residue in the dictionary,
// we get the dictionary (directly) from the srs::Monomer.  This works around a crash that
// happens when we access the dictionary after using monomer functions.  Not tested in main coot yet.
//
GtkWidget *
lbg_info_t::get_image_widget_for_comp_id(const std::string &comp_id, int imol, ccp4srs::Manager *srs_manager) {

   GtkWidget *r = 0;

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#ifdef RDKIT_HAS_CAIRO_SUPPORT

   if (false) {
      // pass this
      int cif_dictionary_read_number = 53;
      std::cout << "calling try_dynamic_add() on " << comp_id << " " << cif_dictionary_read_number
		<< std::endl;
      if (geom_p)
	 geom_p->try_dynamic_add(comp_id, cif_dictionary_read_number++);
   }

   // std::string srs_dir = coot::get_srs_dir(); // use environment variables or ccp4/prefix-dir fall-back
   // = new ccp4srs::Manager;
   // srs_manager->loadIndex(srs_dir.c_str());
   coot::dictionary_residue_restraints_t dict;
   bool status = dict.fill_using_ccp4srs(srs_manager, comp_id);

   if (status) {

      try {
	 RDKit::RWMol rdk_m = coot::rdkit_mol(dict);
	 coot::assign_formal_charges(&rdk_m);
	 coot::rdkit_mol_sanitize(rdk_m);
	 RDKit::RWMol rdk_mol_with_no_Hs = coot::remove_Hs_and_clean(rdk_m);

	 int iconf_2d = RDDepict::compute2DCoords(rdk_mol_with_no_Hs);
	 WedgeMolBonds(rdk_mol_with_no_Hs, &(rdk_mol_with_no_Hs.getConformer(iconf_2d)));

	 std::string smb = RDKit::MolToMolBlock(rdk_mol_with_no_Hs, true, iconf_2d);
	 if (false) { // debug
	    std::string fn = "test-" + comp_id + ".mol";
	    std::ofstream f(fn.c_str());
	    if (f)
	       f << smb << std::endl;
	    f.close();
	 }

	 int n_conf = rdk_mol_with_no_Hs.getNumConformers();
	 // std::cout << "n_conf for " << comp_id << " is " << n_conf << std::endl;
	 if (n_conf > 0) {
	    {
	       RDKit::MolDraw2DCairo drawer(150, 150);
	       drawer.drawMolecule(rdk_mol_with_no_Hs);
	       drawer.finishDrawing();
	       std::string dt = drawer.getDrawingText();

	       if (false) { // debugging.  Mac build has 0 bytes in dt
		  std::string png_file_name = "image-" + comp_id + ".png";
		  std::cout << "writing png " << png_file_name << std::endl;
		  drawer.writeDrawingText(png_file_name.c_str());
		  std::string dt_filename = "dt-" + comp_id + ".png";
		  std::cout << "writing " << dt.length() << " chars to " << dt_filename
			    << std::endl;
		  std::ofstream f(dt_filename.c_str());
		  f << dt;
		  f.close();
	       }

	       // now convert dt to a Pixbuf
	       GError *error = NULL;
	       GdkPixbufLoader *pbl = gdk_pixbuf_loader_new_with_type ("png", &error);
	       guchar tmp[dt.length()];
	       for (unsigned int i=0; i<dt.size(); i++)
		  tmp[i] = dt[i];
	       gboolean write_status = gdk_pixbuf_loader_write(pbl, tmp, dt.length(), &error);
	       gdk_pixbuf_loader_close(pbl, &error);
	       GdkPixbuf *pixbuf = gdk_pixbuf_loader_get_pixbuf(pbl);
	       r = gtk_image_new_from_pixbuf(pixbuf);
	    }
	 }
      }
      catch (...) {
	 std::cout << "WARNING:: hack caught a ... exception " << std::endl;
      }
   } else {
      std::cout << "No dictionary for rdkit_mol from " << comp_id << std::endl;
   }

#endif   // RDKIT_HAS_CAIRO_SUPPORT
#endif   // MAKE_ENHANCED_LIGAND_TOOLS

   return r;
}
#endif // HAVE_CCP4SRS

// static
void
lbg_info_t::on_sbase_search_result_button_clicked (GtkButton *button,
						   gpointer   user_data) {

   std::string *ud = static_cast<std::string *> (user_data);
   if (! ud) {
      std::cout << "Null user_data" << std::endl;

   } else { 
      std::string comp_id = *ud;
      // std::cout << "Do something with " << comp_id << std::endl;

      // 20120110 new style, call an import function, using a pointer.
      lbg_info_t *lbg = (lbg_info_t *) g_object_get_data(G_OBJECT(button), "lbg");
      if (!lbg) { 
	 std::cout << "ERROR NULL lbg in on_sbase_search_result_button_clicked() " << std::endl;
      } else {
// 	 if (! lbg->sbase_import_func_ptr) {
// 	    std::cout << "WARNING:: null SBase import function " << std::endl;
// 	 } else {
// 	    lbg->sbase_import_func_ptr(comp_id);
// 	 }
	 // new style
	 lbg->import_srs_monomer(comp_id);
      }
      
      // 20120110, old-style write out the comp-id to a special file name
      // 
      // std::string file_name = ".sbase-to-coot-comp-id";
      // std::ofstream of(file_name.c_str());
      // /if (of) {
      // of << comp_id;
      // }

   }
}

#endif // EMSCRIPTEN
