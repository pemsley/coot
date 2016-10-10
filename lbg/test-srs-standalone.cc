
#include <iostream>
#include <vector>
#include <cstdlib>
#include <mmdb2/mmdb_manager.h>
#include <mmdb2/mmdb_math_graph.h>
#include <ccp4srs/ccp4srs_manager.h>

// a container for the results of the comparison vs CCP4SRS graph matching.
//
class match_results_t {
public:
   bool success;
   std::string name;
   std::string comp_id;
   mmdb::Residue *res;
   match_results_t(const std::string &comp_id_in, const std::string &name_in, mmdb::Residue *res_in) {
      name = name_in;
      comp_id = comp_id_in;
      res = res_in;
      if (res_in)
	 success = true;
      else
	 success = false;
   }
};
      


// #define HAVE_CCP4SRS

#ifdef TEST_WITH_GEOMETRY
mmdb::math::Graph *
get_graph_1() {
   coot::protein_geometry *geom_p = new coot::protein_geometry;
   const char *d1 = getenv(MONOMER_DIR_STR); // "COOT_CCP4SRS_DIR"
   std::string srs_dir = PKGDATADIR;
   if (d1)
      srs_dir = d1;
   std::cout << "---- geom_p init_ccp4srs with srs_dir " << srs_dir
	     << std::endl;
   geom_p->init_ccp4srs(srs_dir);
   geom_p->try_dynamic_add(monomer_type, true);
   std::pair<bool, coot::dictionary_residue_restraints_t> rest =
      geom_p->get_monomer_restraints(monomer_type);
   mmdb::math::Graph *graph = rest.second.make_graph(inc_Hs);
   return graph;
}
#endif // TEST_WITH_GEOMETRY


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

void
init_ccp4srs(ccp4srs::Manager *ccp4srs) {

   const char *d1 = getenv("COOT_CCP4SRS_DIR");
   std::string srs_dir;
   if (d1)
      srs_dir = d1;
   int RC = ccp4srs->loadIndex(srs_dir.c_str());
   if (RC != ccp4srs::CCP4SRS_Ok) {
      std::cout << "CCP4SRS init problem." << std::endl;
   }
}

#ifdef HAVE_CCP4SRS
std::vector<match_results_t>
compare_vs_ccp4srs(mmdb::math::Graph *graph_1, float similarity, int n_vertices) {

   std::vector<match_results_t> v;

   int minMatch = int(similarity * n_vertices);
   std::cout << "INFO:: match.MatchGraphs must match at least "
	     << minMatch << " atoms from " << similarity << " * " << n_vertices << std::endl;

   ccp4srs::Manager *ccp4srs = new ccp4srs::Manager;
   init_ccp4srs(ccp4srs);

   if (! ccp4srs) {
      std::cout << "WARNING:: CCP4SRS is not initialized" << std::endl;
   } else {
      int l = ccp4srs->n_entries();
      std::cout << "INFO:: compare_vs_ccp4srs(): found " << l << " entries in CCP4 SRS" << std::endl;
      mmdb::math::Graph  *graph_2 = NULL;
      int rc = 0;
      for (int i=1; i<l; i++)  {
	 ccp4srs::Monomer *Monomer = ccp4srs->getMonomer(i, NULL);
	 if (Monomer)  {
	    std::string id = Monomer->ID();
	    // std::cout << "i " << i <<  " monomer id  " << id << std::endl;
	    if (id.length()) {
	       graph_2 = Monomer->getGraph(&rc);
	       if (graph_2) {
		  graph_2->Build(true);
		  graph_2->MakeSymmetryRelief(true);

		  if (rc < 10000) {
		     mmdb::math::GraphMatch match;
		     match.SetTimeLimit(2); // seconds

		     mmdb::math::VERTEX_EXT_TYPE vertex_ext=mmdb::math::EXTTYPE_Equal; // mmdb default
		     bool vertext_type = true;
		     match.MatchGraphs(graph_1, graph_2, minMatch, vertext_type, vertex_ext);
		     int n_match = match.GetNofMatches();
		     if (n_match > 0) {
			std::cout << "INFO:: " << id
				  << " match NumberofMatches (similar graphs): " << n_match << std::endl;

			bool really_match = false;
			for (int imatch=0; imatch<n_match; imatch++) {
			   int n;
			   mmdb::realtype p1, p2;
			   mmdb::ivector FV1, FV2;
			   match.GetMatch(imatch, FV1, FV2, n, p1, p2); // n p1 p2 set
			   std::cout << "   match " << imatch << " matched " << n << " atoms"
				     << std::endl;
			   if (n >= minMatch) really_match = true;
			   for (int ipair=1; ipair<=n; ipair++) {
			      mmdb::math::Vertex *V1 = graph_1->GetVertex(FV1[ipair]);
			      mmdb::math::Vertex *V2 = graph_2->GetVertex(FV2[ipair]);
			      if ((!V1) || (!V2))  {
				 std::cout << "Can't get vertices for match " << ipair << std::endl;
			      } else  {
				 std::cout << "   " << V1->GetUserID() << " " << V2->GetUserID()
					   << std::endl;
			      }
			   }
			   if (really_match)
			      break;
			}

			if (really_match) {
			   mmdb::Residue *residue_p = NULL; // for now
			   std::string name = Monomer->chem_name();
			   match_results_t mr(id, name, residue_p);
			   v.push_back(mr);
			}
		     }
		  }
	       }
	    }
	 }
      }
   }
   return v;
}
#endif // HAVE_CCP4SRS
   

int test_ccp4srs_graph_search() {

   int r = 0;

#ifdef HAVE_CCP4SRS

   double search_similarity = 0.85;
   bool inc_Hs = false;
   std::string monomer_type = "TRP";
   unsigned int n_atoms = 14;

   mmdb::math::Graph *graph = get_graph_2();

   graph->SetName ("Coot-LBG-Query");
   graph->MakeVertexIDs();
   int build_result = graph->Build(true);
   if (build_result != 0) {
      std::cout << "Bad graph build result" << std::endl;
   } else {
      std::cout << "graph build returns status: " << build_result << std::endl;
      graph->MakeSymmetryRelief(true);
      graph->Print();
      std::cout << "graph search using similarity  " << search_similarity << std::endl;
      std::vector<match_results_t> v = compare_vs_ccp4srs(graph, search_similarity, n_atoms);
      delete graph;
      std::cout << "found " << v.size() << " close matches" << std::endl;
      for (unsigned int i=0; i<v.size(); i++) {
	 std::cout << "   " << i << " " << v[i].comp_id << " " << v[i].name << std::endl;
      }
   }
#endif // HAVE_CCP4SRS
   return r;
}


int main(int argc, char **argv) {

   int status = 0;

   test_ccp4srs_graph_search();

   return status;
}
