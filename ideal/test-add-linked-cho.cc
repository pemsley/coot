
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

   for (const auto& edge : node.edges) {
      std::pair<std::string, std::string> res_pair(edge.link_type, edge.target->res_type);
      coot::residue_spec_t new_res_spec =
         coot::cho::add_linked_residue_add_cho_function(mol, imol, node.spec, res_pair,
                                                        new_atoms_b_factor, geom, xmap);
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

   auto root   = std::make_shared<Node>();   root->res_type = "ASN";
   auto child1 = std::make_shared<Node>(); child1->res_type = "NAG";
   auto child2 = std::make_shared<Node>(); child2->res_type = "NAG";
   auto child3 = std::make_shared<Node>(); child3->res_type = "BMA";
   auto child4 = std::make_shared<Node>(); child4->res_type = "MAN";
   auto child5 = std::make_shared<Node>(); child5->res_type = "MAN";

   root->edges.push_back({child1, "pyr-ASN"});
   child1->edges.push_back({child2, "BETA1-4"});
   child2->edges.push_back({child3, "BETA1-4"});
   child3->edges.push_back({child4, "ALPHA1-6"});
   child3->edges.push_back({child5, "ALPHA1-3"});

   traverseTree(*root);


   int imol = 0;
   bool use_gemmi = false;
   atom_selection_container_t asc = get_atom_selection("2qc1-sans-cho.pdb", use_gemmi);
   if (asc.read_success) {

      std::pair<bool, clipper::Xmap<float> > xmap_pair =
         map_from_mtz("2qc1_map.mtz", "FWT", "PHWT", "W", false, false);
      if (xmap_pair.first) {
         const clipper::Xmap<float> &xmap = xmap_pair.second;
         coot::protein_geometry geom;
         geom.init_standard();

         coot::residue_spec_t parent("B", 141, "");
         std::pair<std::string, std::string> res_pair("pyr-ASN", "NAG");
         float new_atoms_b_factor = 30.0;

         // coot::residue_spec_t new_res_spec =
         // coot::cho::add_linked_residue_add_cho_function(asc.mol, imol,
         // parent, res_pair,
         // new_atoms_b_factor,
         // geom, &xmap);

         root->spec = parent;
         traverse_tree_and_build(*root, asc.mol, imol, geom, &xmap, new_atoms_b_factor);
      }
   }
   return status;
}
