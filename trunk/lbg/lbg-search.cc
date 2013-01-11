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
   
      graph->MakeSymmetryRelief(False);
      graph->Print();
      std::cout << "graph search using similarity  " << local_search_similarity << std::endl;
      std::cout << "graph build returns: " << build_result << std::endl;
      std::vector<lbg_info_t::match_results_t> v =
	 compare_vs_sbase(graph, local_search_similarity, n_atoms);
      delete graph;
      display_search_results(v);
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



std::vector<lbg_info_t::match_results_t>
lbg_info_t::compare_vs_sbase(CGraph *graph_1, float similarity, int n_vertices) const {

   std::vector<lbg_info_t::match_results_t> v;
   if (! SBase) {
      std::cout << "SBase not initialized" << std::endl;
   } else {
      //  4.  Run the query through all databsae

      //  There are several methods for retrieving graphs
      //  from the sbase, here we use one most convenient
      //  for serial extractions.
      PCFile graphFile = SBase->GetGraphFile();
      if (!graphFile)  {
	 printf ( "\n SBASE graph file not found.\n" );
      }
  
      int exclude_H_flag = 1;  // neglect hydrogens
      CGraph *graph_2 = NULL;
      int min_match = get_min_match(n_vertices, similarity);

      if (0) // debug
	 std::cout << "min_match " << min_match
		   << " n_vertices: " << n_vertices << " "
		   <<     (similarity * float(n_vertices)) << " " 
		   << int (similarity * float(n_vertices)) << " " 
		   << std::endl;

      std::cout << "Close fragments must match " << min_match << " atoms of "
		<< n_vertices << std::endl;
      

      int nStructures = SBase->GetNofStructures();
      int n_match = 0;
      std::cout << "searching " << nStructures << " SBase structures\n";
      for (int is=0; is<nStructures; is++)  {
	 int rc = SBase->GetGraph(graphFile, graph_2, exclude_H_flag);
	 if (graph_2 == NULL) {
	    std::cout << "bad status on get graph " << is << std::endl;
	 } else {
	    // std::cout << "graph check on structure " << is << std::endl;
	    int n2 = graph_2->GetNofVertices();

	    if ((n2 >= int (double(similarity) * double(n_vertices))) &&
		(n2 < (2.0 - double(similarity)) * double(n_vertices))) { 

	       graph_2->MakeVertexIDs();
	       graph_2->Build(False); // 20100608 was True

	       CGraphMatch *match  = new CGraphMatch();
	       if (min_match > 0) { 

		  match->MatchGraphs(graph_1, graph_2, min_match);
		  int nMatches = match->GetNofMatches();
		  if (nMatches > 0) {
		     if (0)
			std::cout << "found " << nMatches << " match(es) for query in structure "
				  << is << " " << graph_1->GetName() << " vs " << graph_2->GetName()
				  << std::endl;
		     CFile *sf = SBase->GetStructFile();
		     CSBStructure *SBS = SBase->GetStructure(is, sf);
		     if (SBS) {
			if (SBS->Name) {
			   std::cout << "    " << n_match << " " << graph_2->GetName() << " : "
				     << SBS->Name << "\n";
			   lbg_info_t::match_results_t res =
			      residue_from_best_match(*graph_1, *graph_2, *match, nMatches, SBS);
			   v.push_back(res);
			   n_match++;
			}
		     }
		  }
	       }
	       delete match;
	    }
	 }
      }
      std::cout << "Search complete" << std::endl;
  
      graphFile->shut();
      delete graphFile;

      delete graph_2;
   }
   return v;
}



lbg_info_t::match_results_t
lbg_info_t::residue_from_best_match(CGraph &graph1, CGraph &graph2,
				    CGraphMatch &match, int n_match,
				    CSBStructure *SBS) const {

   lbg_info_t::match_results_t r("", "", NULL);
   int best_match = UNASSIGNED_INDEX;
   for (int imatch=0; imatch<n_match; imatch++) {
      r.success = 1;
      r.name = SBS->Name;
      r.comp_id = graph2.GetName();
      int n;
      realtype p1, p2;
      ivector FV1, FV2;
      match.GetMatch(imatch, FV1, FV2, n, p1, p2); // n p1 p2 set
      if (0)
	 std::cout << "   match " << imatch << " " << " set n pairs " << n << std::endl;
      int n_type_match = 0;
      for (int ipair=1; ipair<=n; ipair++) {
	 PCVertex V1 = graph1.GetVertex ( FV1[ipair] );
	 PCVertex V2 = graph2.GetVertex ( FV2[ipair] );
	 if ((!V1) || (!V2))  {
	    std::cout << "Can't get vertices for match " << ipair << std::endl;
	 } else {
	    int type_1 = V1->GetType();
	    int type_2 = V2->GetType();
	    if (type_1 == type_2) {
	       // std::cout << "   type match on " << type_1 << std::endl;
	       n_type_match++;
	    }
	 }
      }
      if (0)
	 std::cout << "This match matches " << n_type_match << " types out of " << n << std::endl;
   }
   return r;
}


// return mmdb sbase return codes
//
// Try to use the MONOMER_DIR_STR, ie. COOT_SBASE_DIR first, if that
// fails then use the fallback directory sbase_monomer_dir_in
// 
int
lbg_info_t::init_sbase(const std::string &sbase_monomer_dir_in) {
      
   int RC = SBASE_FileNotFound; // initial status.
   const char *monomer_dir = getenv(MONOMER_DIR_STR);
   
   if (!monomer_dir) {
      if (coot::is_directory_p(sbase_monomer_dir_in)) {
	 monomer_dir = sbase_monomer_dir_in.c_str();
      } else { 
	 RC = SBASE_FileNotFound; // fail
      }
   }

   if (monomer_dir) { 

      std::cout << "sbase monomer dir: " << monomer_dir << std::endl;
      SBase = new CSBase;
      RC = SBase->LoadIndex(monomer_dir);

      if (RC != SBASE_Ok) {
         std::cout << "sbase files not found in " << monomer_dir << std::endl;
	 delete SBase;
	 SBase = NULL;
      } else { 
         // std::cout << "sbase files found" << std::endl; 
      }
   }
   return RC; 
}

void
lbg_info_t::display_search_results(const std::vector<lbg_info_t::match_results_t> v) const {

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
