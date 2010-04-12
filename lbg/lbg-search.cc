
#include <fstream>

#include <stdlib.h> // for getenv()

#include "clipper/core/coords.h"


#include "lbg.hh"

void
lbg_info_t::search() const {

   // mol.debug();

   CResidue *res = 0;
   
   CGraph *graph = new CGraph;
   int n_atoms = 0;

   // mol atom indexing -> graph vertex indexing
   int vertex_indexing[mol.atoms.size()];

   // atoms are vertices
   for (unsigned int iat=0; iat<mol.atoms.size(); iat++) {
      if (! mol.atoms[iat].is_closed()) {
	 CVertex *v = new CVertex(mol.atoms[iat].element.c_str(), mol.atoms[iat].name.c_str());
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
   graph->SetName ("Liebig-Query");
   int build_result = graph->Build(True);
   graph->MakeSymmetryRelief(False);
   graph->Print();
   std::cout << "graph build returns: " << build_result << std::endl;
   std::vector<lbg_info_t::match_results_t> v =
      compare_vs_sbase(graph, search_similarity, n_atoms);
   delete graph;
   display_search_results(v);
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
      int min_match = get_min_match(n_vertices);

      std::cout << "Close fragments must match " << min_match << " atoms of "
		<< n_vertices << std::endl;
      

      int nStructures = SBase->GetNofStructures();
      std::cout << "searching " << nStructures << " SBase structures\n";
      for (int is=0; is<nStructures; is++)  {
	 int rc = SBase->GetGraph(graphFile, graph_2, exclude_H_flag);
	 if (graph_2 == NULL) {
	    std::cout << "bad status on get graph " << is << std::endl;
	 } else {
	    // std::cout << "graph check on structure " << is << std::endl;
	    int n2 = graph_2->GetNofVertices();

	    graph_2->MakeVertexIDs();
	    graph_2->Build(True);

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
			std::cout << "    " << graph_2->GetName() << " : " << SBS->Name << "\n";
			lbg_info_t::match_results_t res =
			   residue_from_best_match(*graph_1, *graph_2, *match, nMatches, SBS);
			v.push_back(res);
		     }
		  }
	       }
	    }
	    delete match;
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
      InitMatType();

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
      std::string file_name = ".sbase-to-coot-comp-id";
      std::ofstream of(file_name.c_str());
      if (of) {
	 of << comp_id;
      }
   }
}
