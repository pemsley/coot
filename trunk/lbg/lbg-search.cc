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

#ifdef HAVE_GOOCANVAS

#ifdef USE_PYTHON
#include <Python.h>
#endif

#include <fstream>

#include <stdlib.h> // for getenv()

#include "clipper/core/coords.h"


#include "lbg.hh"

#ifdef HAVE_CCP4SRS
void
lbg_info_t::search() const {

   double local_search_similarity = get_search_similarity();

   // mol.debug();

   CResidue *res = 0;
   
   CGraph *graph = new CGraph;
   int n_atoms = 0;

   // mol atom indexing -> graph vertex indexing
   int vertex_indexing[mol.atoms.size()];

   // atoms are vertices
   for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {
      if (! mol.atoms[iat].is_closed()) {
	 std::string ele = mol.atoms[iat].element;
	 std::string name = mol.atoms[iat].get_atom_id();

	 // kludge: 
	 // name = "C" + coot::util::int_to_string(iat+1);
	    
	 CVertex *v = new CVertex(ele.c_str(), name.c_str());
	 vertex_indexing[n_atoms] = iat;
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
// 	    CEdge *e = new CEdge(mol.bonds[ib].get_atom_1_index()+1, // 1-indexed
// 				 mol.bonds[ib].get_atom_2_index()+1,
// 				 mmdb_bond_type);
	    int ind_1 = mol.bonds[ib].get_atom_1_index();  
	    int ind_2 = mol.bonds[ib].get_atom_2_index();
 	    CEdge *e = new CEdge(vertex_indexing[ind_1] + 1,  // 1-indexed
 				 vertex_indexing[ind_2] + 1,
 				 mmdb_bond_type);
	    graph->AddEdge(e);
	    n_bonds++;
	 }
      }
   }
   graph->SetName ("Coot-LBG-Query");
   graph->MakeVertexIDs();
   
   int build_result = graph->Build(False);

   if (build_result != 0) {

      std::cout << "Bad graph build result" << std::endl;

   } else { 

      if (geom_p) { 
	 graph->MakeSymmetryRelief(False);
	 graph->Print();
	 std::cout << "graph search using similarity  " << local_search_similarity << std::endl;
	 std::cout << "graph build returns: " << build_result << std::endl;
	 std::vector<coot::match_results_t> v =
	    geom_p->compare_vs_ccp4srs(graph, local_search_similarity, n_atoms);
	 delete graph;
	 display_search_results(v);
      } else {
	 std::cout << "WARNING:: No geometry in search() " << std::endl;
      } 
   }
}


// allow access of the Search button callback to the search
// similarity combox box text.
double
lbg_info_t::get_search_similarity() const {

   double r = search_similarity;

   gchar *txt = gtk_combo_box_get_active_text(GTK_COMBO_BOX(lbg_search_combobox));
   if (txt) {
      try { 
	 r = lig_build::string_to_float(txt);
      }
      catch (std::runtime_error rte) {
	 std::cout << "WARNING:: Bad number " << txt << " " << rte.what() << std::endl;
	 std::cout << "WARNING:: Using " << r << std::endl;
      }
   } 
   return r;
}
#endif 




#ifdef HAVE_CCP4SRS
PCGraph
lbg_info_t::makeTestQueryGraph() const {
   
   // benzene
   PCGraph  G;

  G = new CGraph();

  G->AddVertex ( new CVertex(" C","C1") );
  G->AddVertex ( new CVertex(" C","C2") );
  G->AddVertex ( new CVertex(" C","C3") );
  G->AddVertex ( new CVertex(" C","C4") );
  G->AddVertex ( new CVertex(" C","C5") );
  G->AddVertex ( new CVertex(" C","C6") );

  G->AddEdge ( new CEdge(1,2,2) );
  G->AddEdge ( new CEdge(2,3,1) );
  G->AddEdge ( new CEdge(3,4,2) );
  G->AddEdge ( new CEdge(4,5,1) );
  G->AddEdge ( new CEdge(5,6,2) );
  G->AddEdge ( new CEdge(6,1,1) );

  G->MakeVertexIDs();
  G->MakeSymmetryRelief ( False );
  // G->Build ( True ); // This makes the search results go crazy.
  G->Build ( False ); // no bond orders!

  return G;
}
#endif 


#ifdef HAVE_CCP4SRS   
coot::match_results_t
lbg_info_t::residue_from_best_match(CGraph &graph_1, CGraph &graph_2,
				    CGraphMatch &match, int n_match, 
				    CCP4SRSMonomer *monomer_p) const {

   coot::match_results_t r("","",NULL);
   if (geom_p)
      r = geom_p->residue_from_best_match(graph_1, graph_2, match, n_match, monomer_p);
   return r;
} 
#endif // HAVE_CCP4SRS

void
lbg_info_t::display_search_results(const std::vector<coot::match_results_t> &v) const {

   gtk_widget_show(lbg_sbase_search_results_dialog);
   
   // clear GTK_BOX(lbg_sbase_search_results_vbox) here
   // 
   GList* glist = gtk_container_get_children(GTK_CONTAINER(lbg_sbase_search_results_vbox));
   while (glist) {
      GtkWidget *w = GTK_WIDGET(glist->data);
      gtk_widget_destroy(w);
      glist = glist->next;
   }

   for (int i=0; i<v.size(); i++) {
      std::string lab = v[i].comp_id;
      lab += ":  ";
      lab += v[i].name;
      GtkWidget *button = gtk_button_new_with_label(lab.c_str());
      gtk_box_pack_start(GTK_BOX(lbg_sbase_search_results_vbox),
			 GTK_WIDGET(button), FALSE, FALSE, 3);
      std::string *comp_id = new std::string(v[i].comp_id);
      g_signal_connect(GTK_WIDGET(button), "clicked",
		       GTK_SIGNAL_FUNC(on_sbase_search_result_button_clicked),
		       (gpointer) (comp_id));
      gtk_object_set_data(GTK_OBJECT(button), "lbg", (gpointer) this);
      gtk_widget_show(button);
   }
}

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
      lbg_info_t *lbg = (lbg_info_t *) gtk_object_get_data(GTK_OBJECT(button), "lbg");
      if (!lbg) { 
	 std::cout << "ERROR NULL lbg in on_sbase_search_result_button_clicked() " << std::endl;
      } else { 
	 if (! lbg->sbase_import_func_ptr) {
	    std::cout << "WARNING:: null SBase import function " << std::endl;
	 } else {
	    lbg->sbase_import_func_ptr(comp_id);
	 }
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

#endif // HAVE_GOOCANVAS
