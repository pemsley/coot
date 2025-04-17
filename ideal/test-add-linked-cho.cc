
#include "clipper/core/xmap.h"
#include "clipper/core/map_utils.h"
#include "clipper/ccp4/ccp4_map_io.h"
#include "clipper/ccp4/ccp4_mtz_io.h"
#include "clipper/core/hkl_compute.h"
#include "clipper/core/map_utils.h" // Map_stats

#include "coot-utils/atom-selection-container.hh"
#include "add-linked-cho.hh"


// copied from mini-rsr. Consider making this a util function
std::pair<bool, clipper::Xmap<float> >
map_from_mtz(std::string mtz_file_name,
	     std::string f_col,
	     std::string phi_col,
	     std::string weight_col,
	     int use_weights,
	     bool is_debug_mode) {


   bool status = 0; // not filled

   clipper::HKL_info myhkl;
   clipper::MTZdataset myset;
   clipper::MTZcrystal myxtl;
   clipper::Xmap<float> xmap;

   try {
      std::cout << "reading mtz file..." << std::endl;
      clipper::CCP4MTZfile mtzin;
      mtzin.open_read( mtz_file_name );       // open new file
      mtzin.import_hkl_info( myhkl );         // read sg, cell, reso, hkls
      clipper::HKL_data< clipper::datatypes::F_sigF<float> >   f_sigf_data(myhkl, myxtl);
      clipper::HKL_data< clipper::datatypes::Phi_fom<float> > phi_fom_data(myhkl, myxtl);
      clipper::HKL_data< clipper::datatypes::F_phi<float> >       fphidata(myhkl, myxtl);

      std::string mol_name = mtz_file_name + " ";
      mol_name += f_col;
      mol_name += " ";
      mol_name += phi_col;

      if (use_weights) {
	 mol_name += " ";
	 mol_name += weight_col;
      }

      if (use_weights) {
	 clipper::String dataname = "/*/*/[" + f_col + " " + f_col + "]";
	 std::cout << dataname << "\n";
	 mtzin.import_hkl_data(  f_sigf_data, myset, myxtl, dataname );
	 dataname = "/*/*/[" + phi_col + " " + weight_col + "]";
	 std::cout << dataname << "\n";
	 mtzin.import_hkl_data( phi_fom_data, myset, myxtl, dataname );
	 mtzin.close_read();
	 std::cout << "We should use the weights: " << weight_col << std::endl;
	 // it seems to me that we should make 2 data types, an F_sigF and a phi fom
	 // and then combine them using a Convert_fsigf_phifom_to_fphi();

	 fphidata.compute(f_sigf_data, phi_fom_data,
			  clipper::datatypes::Compute_fphi_from_fsigf_phifom<float>());

      } else {
	 clipper::String dataname = "/*/*/[" + f_col + " " + phi_col + "]";
	 mtzin.import_hkl_data(     fphidata, myset, myxtl, dataname );
	 mtzin.close_read();
      }
      std::cout << "Number of reflections: " << myhkl.num_reflections() << "\n";
      xmap.init( myhkl.spacegroup(), myhkl.cell(),
		 clipper::Grid_sampling( myhkl.spacegroup(),
					 myhkl.cell(),
					 myhkl.resolution()) );
      std::cout << "Grid..." << xmap.grid_sampling().format() << "\n";
      std::cout << "doing fft..." << std::endl;

      if (is_debug_mode) {
	 int count = 0;
	 clipper::HKL_info::HKL_reference_index hri;
	 for (hri=fphidata.first(); !hri.last(); hri.next()) {
	    if (count == 10)
	       break;
	    std::cout << "sample data " << " "
		      << hri.hkl().h() << " "
		      << hri.hkl().k() << " "
		      << hri.hkl().l() << " : "
		      << fphidata[hri].f() << " " << fphidata[hri].phi()*180/M_PI << std::endl;
	    count++;
	 }
      }
      xmap.fft_from(fphidata);                  // generate map
      std::cout << "done fft..." << std::endl;
      status = 1;
   }
   catch (const clipper::Message_base &exc) {  // "exception" is a protected word, it seems.
      std::cout << "Failed to read mtz file " << mtz_file_name << std::endl;
   }

   return std::pair<bool, clipper::Xmap<float> > (status, xmap);
}

#include <memory>
#include <variant>

struct Node; // Forward declaration

struct Edge {
   std::shared_ptr<Node> target;
   std::string link_type;
};

struct Node {
   std::string res_type;
   std::vector<Edge> edges;
   coot::residue_spec_t spec;
   unsigned int level;
};

void printNodeInfo(const Node& node) {

   std::cout << "Handle Node with Info: " << node.res_type << " spec: " << node.spec << std::endl;
}

void printEdgeInfo(const Edge& edge) {
   std::cout << "  Edge Info: " << edge.link_type << std::endl;
}

void traverseTree(const Node& node) {
    printNodeInfo(node);
    for (const auto& edge : node.edges) {
        printEdgeInfo(edge);
        traverseTree(*edge.target); // Recursive traversal
    }
}

void build_onto_node(const Node &node,
                     mmdb::Manager *mol,
                     int imol,
                     coot::protein_geometry &geom,
                     const clipper::Xmap<float> *xmap,
                     float new_atoms_b_factor) {

   float map_weight = 400.0;

   for (const auto& edge : node.edges) {
      unsigned int new_level = node.level + 1;;
      std::pair<std::string, std::string> res_pair(edge.link_type, edge.target->res_type);
      coot::residue_spec_t new_res_spec =
         coot::cho::add_linked_residue_add_cho_function(mol, imol, node.spec, res_pair, new_level,
                                                        new_atoms_b_factor, geom, xmap, map_weight);
      edge.target->spec = new_res_spec;
   }
}

void traverse_tree_and_build(const Node& node,
                             mmdb::Manager *mol,
                             int imol,
                             coot::protein_geometry &geom,
                             const clipper::Xmap<float> *xmap,
                             float new_atoms_b_factor) {
    printNodeInfo(node);
    build_onto_node(node, mol, imol, geom, xmap, new_atoms_b_factor);
    for (const auto& edge : node.edges) {
        printEdgeInfo(edge);
        traverse_tree_and_build(*edge.target, mol, imol, geom, xmap, new_atoms_b_factor);
    }
}

int main(int argc, char **argv) {

   int  status = 0;

   auto root    = std::make_shared<Node>();   root->res_type  = "ASN";   root->level = 0;
   auto child1  = std::make_shared<Node>(); child1->res_type  = "NAG"; child1->level = 1;
   auto child2  = std::make_shared<Node>(); child2->res_type  = "NAG"; child2->level = 2;
   auto child3  = std::make_shared<Node>(); child3->res_type  = "BMA"; child3->level = 3;
   auto child4  = std::make_shared<Node>(); child4->res_type  = "MAN"; child4->level = 4;
   auto child5  = std::make_shared<Node>(); child5->res_type  = "MAN"; child5->level = 5;
   auto child6  = std::make_shared<Node>(); child6->res_type  = "MAN"; child6->level = 6;
   auto child7  = std::make_shared<Node>(); child7->res_type  = "MAN"; child7->level = 4;
   auto child8  = std::make_shared<Node>(); child8->res_type  = "MAN"; child8->level = 5;
   auto child9  = std::make_shared<Node>(); child9->res_type  = "MAN"; child9->level = 6;
   auto child10 = std::make_shared<Node>(); child10->res_type = "MAN"; child10->level = 5;

   // high mannose

   root->edges.push_back({child1,    "pyr-ASN"});
   child1->edges.push_back({child2,  "BETA1-4"});
   child2->edges.push_back({child3,  "BETA1-4"});
   child3->edges.push_back({child4,  "ALPHA1-3"});
   child4->edges.push_back({child5,  "ALPHA1-2"});
   child5->edges.push_back({child6,  "ALPHA1-2"});

   child3->edges.push_back({child7,  "ALPHA1-6"});
   child7->edges.push_back({child8,  "ALPHA1-6"});
   child8->edges.push_back({child9,  "ALPHA1-2"});
   child7->edges.push_back({child10, "ALPHA1-3"});

   // 8zwp human galactosylltransferase

   root   = std::make_shared<Node>();   root->res_type = "ASN";   root->level = 0;
   child1 = std::make_shared<Node>(); child1->res_type = "NAG"; child1->level = 1;
   child2 = std::make_shared<Node>(); child2->res_type = "NAG"; child2->level = 2;
   child3 = std::make_shared<Node>(); child3->res_type = "BMA"; child3->level = 3;
   child4 = std::make_shared<Node>(); child4->res_type = "MAN"; child4->level = 4;
   child5 = std::make_shared<Node>(); child5->res_type = "MAN"; child5->level = 4;
   child6 = std::make_shared<Node>(); child6->res_type = "FUC"; child6->level = 2;

   root->edges.push_back({child1,    "pyr-ASN"});
   child1->edges.push_back({child2,  "BETA1-4"});
   child2->edges.push_back({child3,  "BETA1-4"});
   child3->edges.push_back({child4,  "ALPHA1-3"});
   child3->edges.push_back({child5,  "ALPHA1-6"});

   child1->edges.push_back({child6,  "ALPHA1-6"});

   traverseTree(*root);

   std::string pdb_file_name = "2qc1-sans-cho.pdb";
   std::string mtz_file_name = "2qc1_map.mtz";
   std::string asn_chain_id = "B";
   int asn_res_no = 141;

   pdb_file_name = "pdb8zwp-sans-cho.pdb";
   mtz_file_name = "8zwp_map.mtz";
   asn_res_no = 174;


   int imol = 0;
   bool use_gemmi = false;
   atom_selection_container_t asc = get_atom_selection(pdb_file_name, use_gemmi);
   if (asc.read_success) {

      std::pair<bool, clipper::Xmap<float> > xmap_pair =
         map_from_mtz(mtz_file_name, "FWT", "PHWT", "W", false, false);
      if (xmap_pair.first) {
         const clipper::Xmap<float> &xmap = xmap_pair.second;
         int cif_read_number = 60;
         coot::protein_geometry geom;
         geom.init_standard();
         std::vector<std::string> av1 = { " C1 ", " C2 ", " C4 ", " C5 "};
         std::vector<std::string> av2 = { " C2 ", " C3 ", " C5 ", " O5 "};
         std::vector<std::string> av3 = { " C3 ", " C4 ", " O5 ", " C1 "};
         std::vector<std::string> rns = {"NAG", "MAN", "BMA", "FUC", "GLC", "GAL", "XYL", "SIA"};
         for (const auto &rn : rns)
            geom.use_unimodal_ring_torsion_restraints(imol, rn, cif_read_number++);
         for (const auto &rn : rns) {
            geom.add_pyranose_pseudo_ring_plane_restraints(rn, imol, "pseudo-plane-1", av1, 0.01);
            geom.add_pyranose_pseudo_ring_plane_restraints(rn, imol, "pseudo-plane-2", av2, 0.01);
            geom.add_pyranose_pseudo_ring_plane_restraints(rn, imol, "pseudo-plane-3", av3, 0.01);
         }

         coot::residue_spec_t parent(asn_chain_id, asn_res_no, "");
         float new_atoms_b_factor = 30.0;

         // std::pair<std::string, std::string> res_pair("pyr-ASN", "NAG");
         // coot::residue_spec_t new_res_spec =
         // coot::cho::add_linked_residue_add_cho_function(asc.mol, imol,
         // parent, res_pair,
         // new_atoms_b_factor,
         // geom, &xmap);

         root->spec = parent;
         traverse_tree_and_build(*root, asc.mol, imol, geom, &xmap, new_atoms_b_factor);
      }
   }
   asc.mol->WritePDBASCII("done.pdb");
   return status;
}
