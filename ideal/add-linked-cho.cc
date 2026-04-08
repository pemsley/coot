
#include <fstream>
#include <memory>
#include <variant>

#include "utils/logging.hh"

#include "coot-utils/atom-selection-container.hh"
#include "extra-restraints.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "coot-utils/coot-map-utils.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "coot-utils/atom-overlaps.hh"
#include "coot-utils/glyco-tree.hh"
#include "coot-utils/glyco-torsions.hh"
#include "coot-utils/dict-link-info.hh"
#include "torsion-bonds.hh"
#include "simple-restraint.hh"

#include "add-linked-cho.hh"
#include "utils/coot-utils.hh"

extern logging logger;

namespace coot {
   namespace cho {
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
         logger.log(log_t::INFO, {logging::ltw("Handle Node with Info: "), logging::ltw(node.res_type),
                                  logging::ltw("spec"), logging::ltw(node.spec.format())});
      }

      void printEdgeInfo(const Edge& edge) {
         // std::cout << "  Edge Info: " << edge.link_type << std::endl;
         logger.log(log_t::INFO, logging::ltw("    Edge Info: "), logging::ltw(edge.link_type));
      }

      void traverseTree(const Node& node) {
         printNodeInfo(node);
         for (const auto& edge : node.edges) {
            printEdgeInfo(edge);
            traverseTree(*edge.target); // Recursive traversal
         }
      }
      void build_onto_node(const coot::cho::Node &node,
                           atom_selection_container_t *asc,
                           int imol,
                           coot::protein_geometry &geom,
                           const clipper::Xmap<float> *xmap,
                           float new_atoms_b_factor);

      void traverse_tree_and_build(const coot::cho::Node& node,
                                   atom_selection_container_t *asc,
                                   int imol,
                                   coot::protein_geometry &geom,
                                   const clipper::Xmap<float> *xmap,
                                   float new_atoms_b_factor);
   }
}

void coot::cho::build_onto_node(const coot::cho::Node &node,
                                atom_selection_container_t *asc,
                                int imol,
                                coot::protein_geometry &geom,
                                const clipper::Xmap<float> *xmap,
                                float new_atoms_b_factor) {

   float map_weight = 40.0; // should be passed

   for (const auto& edge : node.edges) {
      unsigned int new_level = node.level + 1;;
      std::pair<std::string, std::string> res_pair(edge.link_type, edge.target->res_type);
      coot::residue_spec_t new_res_spec =
         coot::cho::add_linked_residue_add_cho_function(asc, imol, node.spec, res_pair, new_level,
                                                        new_atoms_b_factor, geom, xmap, map_weight);
      edge.target->spec = new_res_spec;
   }
}

void
coot::cho::traverse_tree_and_build(const coot::cho::Node& node,
                                   atom_selection_container_t *asc,
                                   int imol,
                                   coot::protein_geometry &geom,
                                   const clipper::Xmap<float> *xmap,
                                   float new_atoms_b_factor) {
    printNodeInfo(node);
    // debug_atom_selection_container(*asc);
    build_onto_node(node, asc, imol, geom, xmap, new_atoms_b_factor);
    for (const auto& edge : node.edges) {
       { // debugging
          asc->user_data++;
          // std::cout << "------------------------------------------ debug_write_pdb() " << asc->user_data << std::endl;
          // asc->debug_write_pdb();
       }
       // std::cout << "=============================== edge info with count " << asc->user_data << std::endl;
       printEdgeInfo(edge);
       traverse_tree_and_build(*edge.target, asc, imol, geom, xmap, new_atoms_b_factor);
    }
}

// "NAG-NAG-BMA" or "high-mannose" or "hybrid" or "mammalian-biantennary" or "plant-biantennary"
//
void
coot::cho::add_named_glyco_tree(const std::string &glycoylation_name, atom_selection_container_t *asc, int imol,
                                float new_atoms_b_factor,
                                const clipper::Xmap<float> &xmap, coot::protein_geometry *geom,
                                std::string asn_chain_id, int asn_res_no) {

   auto make_high_mannose_tree = [] () {
      // high mannose

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

      return root;
   };

   // e.g 8zwp human galactosylltransferase

   auto make_hybrid = [] () {

      auto root   = std::make_shared<Node>();   root->res_type = "ASN";   root->level = 0;
      auto child1 = std::make_shared<Node>(); child1->res_type = "NAG"; child1->level = 1;
      auto child2 = std::make_shared<Node>(); child2->res_type = "NAG"; child2->level = 2;
      auto child3 = std::make_shared<Node>(); child3->res_type = "BMA"; child3->level = 3;
      auto child4 = std::make_shared<Node>(); child4->res_type = "MAN"; child4->level = 4;
      auto child5 = std::make_shared<Node>(); child5->res_type = "MAN"; child5->level = 4;
      auto child6 = std::make_shared<Node>(); child6->res_type = "FUC"; child6->level = 2;

      root->edges.push_back({child1,    "pyr-ASN"});
      child1->edges.push_back({child2,  "BETA1-4"});
      child2->edges.push_back({child3,  "BETA1-4"});
      child3->edges.push_back({child4,  "ALPHA1-3"});
      child3->edges.push_back({child5,  "ALPHA1-6"});
      child1->edges.push_back({child6,  "ALPHA1-6"});

      return root;
   };

   //"mammalian-biantennary" or "plant-biantennary"

   auto make_mammalian_biantennary = [] () {

      auto root   = std::make_shared<Node>();   root->res_type = "ASN";   root->level = 0;
      auto child1 = std::make_shared<Node>(); child1->res_type = "NAG"; child1->level = 1;
      auto child2 = std::make_shared<Node>(); child2->res_type = "NAG"; child2->level = 2;
      auto child3 = std::make_shared<Node>(); child3->res_type = "BMA"; child3->level = 3;
      auto child4 = std::make_shared<Node>(); child4->res_type = "MAN"; child4->level = 4;
      auto child5 = std::make_shared<Node>(); child5->res_type = "MAN"; child5->level = 4;

      root->edges.push_back({child1,    "pyr-ASN"});
      child1->edges.push_back({child2,  "BETA1-4"});
      child2->edges.push_back({child3,  "BETA1-4"});
      child3->edges.push_back({child4,  "ALPHA1-3"});
      child3->edges.push_back({child5,  "ALPHA1-6"});

      return root;
   };

   auto make_plant_biantennary = [] () {

      auto root   = std::make_shared<Node>();   root->res_type = "ASN";   root->level = 0;
      auto child1 = std::make_shared<Node>(); child1->res_type = "NAG"; child1->level = 1;
      auto child2 = std::make_shared<Node>(); child2->res_type = "NAG"; child2->level = 2;
      auto child3 = std::make_shared<Node>(); child3->res_type = "BMA"; child3->level = 3;
      auto child4 = std::make_shared<Node>(); child4->res_type = "MAN"; child4->level = 4;
      auto child5 = std::make_shared<Node>(); child5->res_type = "MAN"; child5->level = 4;

      root->edges.push_back({child1,    "pyr-ASN"});
      child1->edges.push_back({child2,  "BETA1-4"});
      child2->edges.push_back({child3,  "BETA1-4"});
      child3->edges.push_back({child4,  "ALPHA1-3"});
      child3->edges.push_back({child5,  "ALPHA1-6"});

      return root;
   };

   auto make_NAG_NAG_BMA = [] () {

      auto root   = std::make_shared<Node>();   root->res_type = "ASN";   root->level = 0;
      auto child1 = std::make_shared<Node>(); child1->res_type = "NAG"; child1->level = 1;
      auto child2 = std::make_shared<Node>(); child2->res_type = "NAG"; child2->level = 2;
      auto child3 = std::make_shared<Node>(); child3->res_type = "BMA"; child3->level = 3;

      root->edges.push_back({child1,    "pyr-ASN"});
      child1->edges.push_back({child2,  "BETA1-4"});
      child2->edges.push_back({child3,  "BETA1-4"});

      return root;
   };

   // "NAG-NAG-BMA" or "high-mannose" or "hybrid" or "mammalian-biantennary" or "plant-biantennary"
   //
   std::shared_ptr<Node> root;
   if (glycoylation_name == "NAG-NAG-BMA")           root = make_NAG_NAG_BMA();
   if (glycoylation_name == "high-mannose")          root = make_high_mannose_tree();
   if (glycoylation_name == "hybrid")                root = make_hybrid();
   if (glycoylation_name == "mammalian-biantennary") root = make_mammalian_biantennary();
   if (glycoylation_name == "plant-biantennary")     root = make_plant_biantennary();

   if (root) {

      traverseTree(*root); // testing

      std::vector<std::string> av1 = { " C1 ", " C2 ", " C4 ", " C5 "};
      std::vector<std::string> av2 = { " C2 ", " C3 ", " C5 ", " O5 "};
      std::vector<std::string> av3 = { " C3 ", " C4 ", " O5 ", " C1 "};
      std::vector<std::string> rns = {"NAG", "MAN", "BMA", "FUC", "GLC", "GAL", "XYL", "SIA"};
      int cif_read_number = 60;
      for (const auto &rn : rns)
         geom->use_unimodal_ring_torsion_restraints(imol, rn, cif_read_number++);
      for (const auto &rn : rns) {
         geom->add_pyranose_pseudo_ring_plane_restraints(rn, imol, "pseudo-plane-1", av1, 0.01);
         geom->add_pyranose_pseudo_ring_plane_restraints(rn, imol, "pseudo-plane-2", av2, 0.01);
         geom->add_pyranose_pseudo_ring_plane_restraints(rn, imol, "pseudo-plane-3", av3, 0.01);
      }

      coot::residue_spec_t parent(asn_chain_id, asn_res_no, "");

      // std::pair<std::string, std::string> res_pair("pyr-ASN", "NAG");
      // coot::residue_spec_t new_res_spec =
      // coot::cho::add_linked_residue_add_cho_function(asc.mol, imol,
      // parent, res_pair,
      // new_atoms_b_factor,
      // geom, &xmap);

      root->spec = parent;

      traverse_tree_and_build(*root, asc, imol, *geom, &xmap, new_atoms_b_factor);
   }

 }


int
coot::cho::clashes_with_symmetry(mmdb::Manager *mol, const coot::residue_spec_t &res_spec, float clash_dist,
                                 const coot::protein_geometry &geom) {

   int r = -1;
   mmdb::Residue *residue_p = util::get_residue(res_spec, mol);
   if (mol) {
      if (residue_p) {
         std::vector<mmdb::Residue *> dummy; // neighbours
         atom_overlaps_container_t ao(residue_p, dummy, mol, &geom);
         std::vector<coot::atom_overlap_t> v = ao.symmetry_contacts(clash_dist);
         if (v.empty())
            r = 0;
         else
            r = 1;
      }
   }
   return r;
}

bool
coot::cho::is_well_fitting(mmdb::Residue *residue_p,
                           mmdb::Manager *mol,
                           const clipper::Xmap<float> &xmap,
                           const coot::protein_geometry &geom) {

   float add_linked_residue_tree_correlation_cut_off = 0.50;
   float clash_dist = 2.0;
   float atom_radius = 1.6;

   bool status = false;
   float radius = 4.0;
   residue_spec_t res_spec(residue_p);
   std::vector<mmdb::Residue *> neighbours = residues_near_residue(residue_p, mol, radius);
   std::vector<residue_spec_t> residues_for_masking;
   for(mmdb::Residue *r : neighbours)
      residues_for_masking.push_back(residue_spec_t(r));
   std::vector<residue_spec_t> residues_for_cc = { res_spec };
   unsigned short int atom_mask_mode = 0; // all atom

   float c = util::map_to_model_correlation(mol, residues_for_cc, residues_for_masking, atom_mask_mode, atom_radius, xmap);
   // std::cout << "debug:: is_well_fitting: " << coot::residue_spec_t(residue_p) << " correllation: " << c << std::endl;
   logger.log(log_t::DEBUG, std::string("is_well_fitting:"), coot::residue_spec_t(residue_p).format(),
              std::string("correlation:"), c);
   logger.show_last();
   if (c > add_linked_residue_tree_correlation_cut_off) {
      int symm_clash = clashes_with_symmetry(mol, res_spec, clash_dist, geom);
      if (symm_clash == 0) {
         status = true;
      }
   }
   return status;
}



bool
coot::cho::is_het_residue(mmdb::Residue *residue_p) {

   bool status = false;

   if (residue_p) {
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for(int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            if (at->Het) {
               status =  true;
               break;
            }
         }
      }
   }
   return status;
}


// return state, max_resno + 1, or 0, 1 of no residues in chain.
//
// new_res_no_by_hundreds is default false
std::pair<short int, int>
coot::cho::next_residue_number_in_chain(mmdb::Chain *w,
                                        bool new_res_no_by_hundreds) {

   std::pair<short int, int> p(0,1);
   int max_res_no = -9999;

   if (w) {
      int nres = w->GetNumberOfResidues();
      mmdb::Residue *residue_p;
      if (nres > 0) {
         for (int ires=nres-1; ires>=0; ires--) {
            residue_p = w->GetResidue(ires);
            if (residue_p->seqNum > max_res_no) {
               max_res_no = residue_p->seqNum;
               bool is_het_residue_flag = is_het_residue(residue_p);
               if (is_het_residue_flag) {
                  p = std::pair<short int, int>(1, residue_p->seqNum+1);
               } else {
                  if (new_res_no_by_hundreds) {
                     if (max_res_no < 9999) {
                        int res_no = coot::util::round_up_by_hundreds(max_res_no+1);
                        p = std::pair<short int, int>(1, res_no+1);
                     }
                  } else {
                     if (max_res_no < 9999) {
                        p = std::pair<short int, int>(1, max_res_no+1);
                     }
                  }
               }
            }
         }
         if (! p.first) {
            //  first the first space starting from the front
            int test_resno_start = 1001;
            bool is_clear = false;
            while (! is_clear) {
               is_clear = true;
               for (int iser=0; iser<nres; iser++) {
                  int resno_res = w->GetResidue(iser)->seqNum;
                  if (resno_res >= test_resno_start) {
                     if (resno_res <= (test_resno_start+10)) {
                        is_clear = false;
                     }
                  }
                  if (! is_clear)
                     break;
               }
               test_resno_start += 100;
            }
            p = std::pair<short int, int> (1, test_resno_start);
         }
      }
   }
   return p;
}


// This doesn't do a backup or finalise model.
mmdb::Residue *
coot::cho::copy_and_add_residue_to_chain(mmdb::Manager *mol,
                                         mmdb::Chain *this_model_chain,
                                         mmdb::Residue *add_model_residue,
                                         bool new_resno_by_hundreds_flag) {

   mmdb::Residue *res_copied = NULL;
   if (add_model_residue) {
      bool whole_res_flag = true;
      int udd_atom_index_handle = 1; // does this matter?
      bool add_this = true;
      // check for overlapping water (could be generalised for same residue type?!
      std::vector<mmdb::Residue *> close_residues;
      close_residues = coot::residues_near_residue(add_model_residue, mol, 0.05);
      for (unsigned int i=0; i<close_residues.size(); i++) {
         if (close_residues[i]->isSolvent() && add_model_residue->isSolvent()) {
            add_this = false;
            // std::cout << "INFO:: not adding water because of overlap" << std::endl;
            logger.log(log_t::INFO, std::string("not adding water because of overlap"));
            break;
         }
      }
      if (add_this) {

         /* No - this does an implicit embed-in-chain - that is not what we want
         mmdb::Residue *residue_copy = coot::deep_copy_this_residue(add_model_residue,
                                                                    "",
                                                                    whole_res_flag,
                                                                    udd_atom_index_handle);
         */
         mmdb::Residue *residue_copy = coot::util::deep_copy_this_residue(add_model_residue);

         if (residue_copy) {
            std::pair<short int, int> res_info =
               next_residue_number_in_chain(this_model_chain, new_resno_by_hundreds_flag);
            int new_res_resno = 9999;
            if (res_info.first)
               new_res_resno = res_info.second;
            residue_copy->seqNum = new_res_resno; // try changing the seqNum before AddResidue().
            this_model_chain->AddResidue(residue_copy);
            res_copied = residue_copy;
         }
      }
   }
   return res_copied;
}

void
coot::cho::asn_hydrogen_position_swap(std::vector<std::pair<bool, mmdb::Residue *> > residues) {

   if (residues[0].second) {
      if (residues[1].second) {
	 std::string rn0(residues[0].second->GetResName());
	 std::string rn1(residues[1].second->GetResName());
	 mmdb::Residue *r_0 = 0;
	 mmdb::Residue *r_1 = 0;
	 if (rn0 == "ASN") {
	    if (rn1 == "NAG") {
	       r_0 = residues[0].second;
	       r_1 = residues[1].second;
	    }
	 }
	 if (rn1 == "ASN") {
	    if (rn0 == "NAG") {
	       r_1 = residues[0].second;
	       r_0 = residues[1].second;
	    }
	 }

	 if (r_1 && r_0) {
	    mmdb::Atom *at_hd21 = 0;
	    mmdb::Atom *at_hd22 = 0;
	    mmdb::Atom **residue_atoms_0 = 0;
	    int n_residue_atoms_0;
	    r_0->GetAtomTable(residue_atoms_0, n_residue_atoms_0);
	    for (int iat=0; iat<n_residue_atoms_0; iat++) {
	       mmdb::Atom *at = residue_atoms_0[iat];
	       std::string atom_name(at->GetAtomName());
	       if (atom_name == "HD21") at_hd21 = at;
	       if (atom_name == "HD22") at_hd22 = at;
	    }
	    if (at_hd21 && at_hd22) {
	       clipper::Coord_orth co21 = coot::co(at_hd21);
	       clipper::Coord_orth co22 = coot::co(at_hd22);
	       at_hd21->x = co22.x();
	       at_hd21->y = co22.y();
	       at_hd21->z = co22.z();
	       at_hd22->x = co21.x(); // this atom will be deleted.
	       at_hd22->y = co21.y();
	       at_hd22->z = co21.z();
	    }
	 }
      }
   }
}

#include <string.h>

void
coot::cho::make_link(mmdb::Manager *mol, const coot::atom_spec_t &spec_1,
                          const coot::atom_spec_t &spec_2,
                          const std::string &link_name, float length,
                          const coot::protein_geometry &geom) {

   // 2014: link_name and length are not part curently of a mmdb::Link.
   // Perhaps they should not be passed then?

   mmdb::Atom *at_1 = util::get_atom(spec_1, mol);
   mmdb::Atom *at_2 = util::get_atom(spec_2, mol);

   if (! at_1) {
      // std::cout << "WARNING:: atom " << spec_1 << " not found - abandoning LINK addition " << std::endl;
      logger.log(log_t::WARNING, std::string("atom"), spec_1.format(), "not found - abandoning LINK addtion");
      logger.show_last();
   } else {
      if (! at_2) {
	 // std::cout << "WARNING:: atom " << spec_1 << " not found - abandoning LINK addition " << std::endl;
         logger.log(log_t::WARNING, std::string("atom"), spec_2.format(), "not found - abandoning LINK addtion");
         logger.show_last();
      } else {

	 mmdb::Model *model_1 = at_1->GetModel();
	 mmdb::Model *model_2 = at_1->GetModel();

	 if (model_1 != model_2) {

	    std::cout << "WARNING:: specified atoms have mismatching models - abandoning LINK addition"
		      << std::endl;

	 } else {

            mmdb::Link *link = new mmdb::Link; // sym ids default to 1555 1555

	    strncpy(link->atName1,  at_1->GetAtomName(), 20);
	    strncpy(link->aloc1,    at_1->altLoc, 20);
	    strncpy(link->resName1, at_1->GetResName(), 19);
	    strncpy(link->chainID1, at_1->GetChainID(), 9);
	    strncpy(link->insCode1, at_1->GetInsCode(), 9);
	    link->seqNum1         = at_1->GetSeqNum();

	    strncpy(link->atName2,  at_2->GetAtomName(), 20);
	    strncpy(link->aloc2,    at_2->altLoc, 20);
	    strncpy(link->resName2, at_2->GetResName(), 19);
	    strncpy(link->chainID2, at_2->GetChainID(), 9);
	    strncpy(link->insCode2, at_2->GetInsCode(), 9);
	    link->seqNum2         = at_2->GetSeqNum();

	    model_1->AddLink(link);
	    mol->FinishStructEdit();

	    // now, do we need to do any deletions to the model that
	    // are defined in the dictionary?
	    //
	    std::vector<std::pair<bool, mmdb::Residue *> > residues(2);
	    residues[0] = std::pair<bool, mmdb::Residue *> (0, at_1->residue);
	    residues[1] = std::pair<bool, mmdb::Residue *> (0, at_2->residue);
	    std::vector<coot::atom_spec_t> dummy_fixed_atom_specs;

	    // convert to restraints_container_t interface
	    mmdb::Link local_link;
	    local_link.Copy(link);
	    std::vector<mmdb::Link> links(1);
	    links[0] = local_link;
	    clipper::Xmap<float> dummy_xmap;

	    restraints_container_t restraints(residues, links, geom, mol,
                                              dummy_fixed_atom_specs, &dummy_xmap);
	    bonded_pair_container_t bpc = restraints.bonded_residues_from_res_vec(geom);

	    asn_hydrogen_position_swap(residues); // HD21 and HD22 (HD22 will be deleted)
	    bpc.apply_chem_mods(geom);
	    mol->FinishStructEdit();
	 }
      }
   }
}



std::pair<bool, mmdb::Residue *>
coot::cho::add_residue(mmdb::Manager *mol,
                       mmdb::Residue *new_res,
                       const std::string &chain_id_in) {

   bool status = false;
   mmdb::Residue *res_copied = NULL;
   bool new_resno_by_hundreds_flag =  true;
   int imod = 1;
   if (new_res) {
      if (mol) {
	 mmdb::Model *model_p = mol->GetModel(imod);
	 mmdb::Chain *chain_p;
	 if (model_p) {
	    int n_chains = model_p->GetNumberOfChains();
	    for (int ichain=0; ichain<n_chains; ichain++) {
	       chain_p = model_p->GetChain(ichain);
	       std::string chain_id(chain_p->GetChainID());
	       if (chain_id == chain_id_in) {
		  res_copied = copy_and_add_residue_to_chain(mol, chain_p, new_res, new_resno_by_hundreds_flag);
		  status = true;
		  mol->FinishStructEdit();
		  break;
	       }
	    }
	 }
      }
   }
   return std::pair<bool, mmdb::Residue *> (status, res_copied);
}


coot::residue_spec_t
coot::cho::add_linked_residue_by_atom_torsions(mmdb::Manager *mol,
                                               const residue_spec_t &parent,
                                               const std::pair<std::string, std::string> &new_link_types,
                                               protein_geometry &geom,
                                               float new_atoms_b_factor) {

   coot::residue_spec_t new_residue_spec;
   mmdb::Residue *residue_ref = util::get_residue(parent, mol);
   if (residue_ref) {
      try {
         const std::string &new_link_type = new_link_types.first;
         const std::string &new_res_type  = new_link_types.second;
         link_by_torsion_t l(new_link_type, new_res_type);
         l.set_temperature_factor(new_atoms_b_factor);
         mmdb::Residue *result_residue = l.make_residue(residue_ref);
         mol->FinishStructEdit();
         std::pair<bool, mmdb::Residue *> status_pair = add_residue(mol, result_residue, parent.chain_id);
         if (status_pair.first) {
            mmdb::Residue *residue_new = status_pair.second;
            new_residue_spec = residue_spec_t(status_pair.second);
            dict_link_info_t link_info(residue_ref, residue_new, new_link_type, geom);
            make_link(mol, link_info.spec_ref, link_info.spec_new, new_link_type, link_info.dist, geom);
         }
      }
      catch (const std::runtime_error &rte) {
         logger.log(log_t::WARNING, rte.what());
         logger.show_last();
      }
   }
   return new_residue_spec;
}

// from backrub-rotamer.cc
void
coot::cho::replace_coords(mmdb::Manager *fragment_mol, mmdb::Manager *mol) {
   if (mol) {
      int imod = 1;
      mmdb::Model *model_p = fragment_mol->GetModel(imod);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            int nres = chain_p->GetNumberOfResidues();
            for (int ires=0; ires<nres; ires++) {
               mmdb::Residue *residue_p = chain_p->GetResidue(ires);
               int n_atoms = residue_p->GetNumberOfAtoms();
               for (int iat=0; iat<n_atoms; iat++) {
                  mmdb::Atom *at = residue_p->GetAtom(iat);
                  atom_spec_t spec(at);
                  mmdb::Atom *at_mol = util::get_atom(spec, mol);
                  if (at_mol) {
                     at_mol->x = at->x;
                     at_mol->y = at->y;
                     at_mol->z = at->z;
                  }
               }
            }
         }
      }
   }
}

#include  <filesystem>

coot::residue_spec_t
coot::cho::add_linked_residue(atom_selection_container_t *asc,
                              int imol,
                              const residue_spec_t &parent,
                              const std::pair<std::string, std::string> &new_link_types,
                              unsigned int level,
                              float new_atoms_b_factor,
                              int mode,
                              coot::protein_geometry &geom,
                              const clipper::Xmap<float> *xmap,
                              float map_weight) {


   auto make_refinement_residues = [](const std::vector<residue_spec_t> &residue_specs,
                                      mmdb::Manager *mol) {

      std::vector<std::pair<bool,mmdb::Residue *> > residues;
      for (const auto &r : residue_specs) {
         mmdb::Residue *residue_p = util::get_residue(r, mol);
         if (residue_p) {
            residues.push_back(std::make_pair(false, residue_p));
         }
      }
      return residues;
   };

   auto debug_coords = [] (mmdb::Manager *mol) {
      int icount = 0;
      for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
         mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               int n_res = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<n_res; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     int n_atoms = residue_p->GetNumberOfAtoms();
                     for (int iat=0; iat<n_atoms; iat++) {
                        mmdb::Atom *at = residue_p->GetAtom(iat);
                        if (! at->isTer()) {
                           if (icount < 10) {
                              std::cout << "atom: " << icount << " "
                                        << at->GetChainID() << " "
                                        << at->residue->GetSeqNum() << " "
                                        << at->GetAtomName() << " "
                                        << at->x << " " << at->y << " " << at->z << " "
                                        << std::endl;
                           } else {
                              break;
                           }
                           icount++;
                        }
                     }
                  }
               }
            }
         }
      }
   };

   auto get_glyco_dir_path = [] () -> std::filesystem::path {
      // look in package_data_dir ... share/coot/data/cho-models
      std::string pkg_data_dir = package_data_dir();
      std::filesystem::path pkg_data_path = std::filesystem::path(pkg_data_dir);
      std::filesystem::path cho_models_path = pkg_data_path / "data" / "cho-models";
      // std::cout << "debug:: pkd_data_dir is " << pkg_data_dir << std::endl;
      return cho_models_path;
   };

   class atom_pair_restraints_t {
   public:
      std::string atom_name_1;
      std::string atom_name_2;
      residue_spec_t parent_residue_spec;
      residue_spec_t new_residue_spec;
      double target_dist;
      double sigma;
      atom_pair_restraints_t(const std::string &a1, const std::string &a2,
                             const residue_spec_t p, const residue_spec_t &n,
                             double t, double s) : atom_name_1(a1), atom_name_2(a2),
                                                   parent_residue_spec(p), new_residue_spec(n),
                                                   target_dist(t), sigma(s) {}
   };

   auto add_extra_restraints = [get_glyco_dir_path] (coot::restraints_container_t *restraints_p,
                                                     int imol,
                                                     unsigned int level,
                                                     const residue_spec_t &parent_spec,
                                                     const std::string &parent_type,
                                                     const residue_spec_t &new_residue_spec,
                                                     const std::string &new_type,
                                                     const std::string &link_type,
                                                     const protein_geometry &geom) {

      std::vector<atom_pair_restraints_t> extra_restraints;

      //  new-type + link + parent-type
      std::string fn = "model-level-" + std::to_string(level) + "-" + new_type + "-" + link_type + "-" + new_type + ".tab";
      std::filesystem::path cho_models_path = get_glyco_dir_path();
      std::filesystem::path model_path = cho_models_path / fn;
      // std::cout << "model file name: " << model_path.string() << std::endl;
      logger.log(log_t::DEBUG, "model file name", model_path.string());
      if (std::filesystem::exists(model_path)) {
         std::ifstream f(model_path);
         if (f) {
            std::vector<std::string> lines;
            std::string line;
            while (std::getline(f, line)) {
               lines.push_back(line);
            }
            if (! lines.empty()) {
               for (const auto &line : lines) {
                  std::vector<std::string> parts = coot::util::split_string(line, " ");
                  if (parts.size() == 7) {
                     const std::string atom_name_1 = parts[0];
                     const std::string atom_name_2 = parts[1];
                     double target_dist = coot::util::string_to_double(parts[2]);
                     double rmsd        = coot::util::string_to_double(parts[3]);
                     unsigned int n     = coot::util::string_to_int(parts[4]);
                     double d           = coot::util::string_to_double(parts[5]);
                     // parts 6 is mod-sarle
                     if (n > 20) {
                        if (d < 0.42) {
                           double s = 0.05;
                           atom_pair_restraints_t apr(atom_name_1, atom_name_2, parent_spec, new_residue_spec,
                                                      target_dist, s);
                           extra_restraints.push_back(apr);
                        }
                     }
                  }
               }
            }
            f.close();
         }
      }

      // let's add them outside the loop
      if (!extra_restraints.empty()) {
         extra_restraints_t er;
         for (const auto &r : extra_restraints) {
            coot::atom_spec_t atom_spec_1(r.parent_residue_spec.chain_id, r.parent_residue_spec.res_no,
                                          r.parent_residue_spec.ins_code, r.atom_name_1, "");
            coot::atom_spec_t atom_spec_2(r.new_residue_spec.chain_id, r.new_residue_spec.res_no,
                                          r.new_residue_spec.ins_code, r.atom_name_2, "");
            extra_restraints_t::extra_bond_restraint_t bond(atom_spec_1, atom_spec_2, r.target_dist, r.sigma);
            er.bond_restraints.push_back(bond);
         }
         std::string descr("inter-residue CHO distance retraints");
         restraints_p->add_extra_restraints(imol, descr, er, geom);
      }
   };

   auto refine_direct = [add_extra_restraints] (const std::vector<std::pair<bool,mmdb::Residue *> > &residues,
                                                mmdb::Manager *mol, int imol_no,
                                                const std::vector<mmdb::Link> &links,
                                                unsigned int level, // for extra restraints
                                                const residue_spec_t &parent_spec,
                                                const std::string &parent_type,
                                                const residue_spec_t &new_residue_spec,
                                                const std::string &new_type,
                                                const std::string &link_type,
                                                const protein_geometry &geom,
                                                const std::vector<coot::atom_spec_t> &fixed_atom_specs,
                                                const clipper::Xmap<float> &xmap,
                                                double map_weight,
                                                double torsion_restraints_weight) {

      unsigned int number_of_threads = 4;
      unsigned int n_cycles = 1000;
      coot::restraints_container_t restraints(residues, links, geom, mol, fixed_atom_specs, &xmap);

      if (false)
         restraints.set_quiet_reporting();
      restraints.set_torsion_restraints_weight(torsion_restraints_weight);
      restraints.add_map(map_weight);
      restraint_usage_Flags flags = TYPICAL_RESTRAINTS_WITH_TORSIONS;
      pseudo_restraint_bond_type pseudos = NO_PSEUDO_BONDS;

      ctpl::thread_pool local_thread_pool(number_of_threads);
      restraints.thread_pool(&local_thread_pool, number_of_threads);
      int imol = imol_no;
      bool do_auto_helix_restraints = true;
      bool do_auto_strand_restraints = false;
      bool do_h_bond_restraints = false;
      bool do_residue_internal_torsions = true;
      bool make_trans_peptide_restraints = true;
      double rama_plot_weight = 1.0;
      bool do_rama_plot_restraints = false;
      bool print_chi_sq_flag = true;
      restraints.make_restraints(imol, geom, flags,
                                 do_residue_internal_torsions,
                                 make_trans_peptide_restraints,
                                 rama_plot_weight, do_rama_plot_restraints,
                                 do_auto_helix_restraints,
                                 do_auto_strand_restraints,
                                 do_h_bond_restraints,
                                 pseudos);
      add_extra_restraints(&restraints, imol_no, level, parent_spec, parent_type,
                           new_residue_spec, new_type, link_type, geom); // modify restraints
      restraints.minimize(flags, n_cycles, print_chi_sq_flag);

   };

   int n_trials = 500;
   int cif_dictionary_read_number = 60; // mergh. Pass a reference to this?

   const std::string &new_link     = new_link_types.first;
   const std::string &new_res_type = new_link_types.second;
   if (geom.have_dictionary_for_residue_type_no_dynamic_add(new_res_type, imol)) {
   } else {
      int status = geom.try_dynamic_add(new_res_type, cif_dictionary_read_number);
      if (status == 0) { // fail
         // std::cout << "WARNING:: failed to add dictionary for type \"" << new_res_type << "\"" << std::endl;
         // std::cout << "WARNING:: add_linked_residue() stops here " << std::endl;
         logger.log(log_t::WARNING, "failed to add dictionary for type", new_res_type,
                    "add_linked_residue() stops here");
         logger.show_last();
         return residue_spec_t();
      }
   }
   cif_dictionary_read_number++;
   mmdb::Manager *mol = asc->mol;
   residue_spec_t new_res_spec = add_linked_residue_by_atom_torsions(mol, parent, new_link_types, geom, new_atoms_b_factor);
   asc->regen_atom_selection();
   // delete extra restraints for new_res_spec - hmm
   if (! new_res_spec.unset_p()) {
      if (mode == 2 || mode == 3) {
         if (xmap) {

            std::vector<residue_spec_t > residue_specs;
            residue_specs.push_back(parent);
            residue_specs.push_back(new_res_spec);
            std::string parent_type;
            mmdb::Residue *parent_residue_p = util::get_residue(parent, mol);
            if (parent_residue_p)
               parent_type = parent_residue_p->GetResName();

            if (false) {
               std::vector<mmdb::Residue *> rv = residues_near_residue(parent_residue_p, mol, 7.5);
               for (auto &r : rv) {
                  bool fs = coot::cho::is_well_fitting(r, mol, *xmap, geom);
                  std::cout << "debug:: well-fitting: " << coot::residue_spec_t(r) << std::endl;
               }
            }

            coot::residue_spec_t parent_spec(parent_residue_p);
            mmdb::Manager *moving_mol = util::create_mmdbmanager_from_residue_specs(residue_specs, mol);

            int n_rounds_of_fit_and_refine = 1;
            std::vector<std::pair<bool, clipper::Coord_orth> > avoid_these_atoms;
            for (int ii=0; ii<n_rounds_of_fit_and_refine; ii++) {
               // std::string file_name = "mrtf-pre-" + std::to_string(ii) + ".pdb";
               // moving_mol->WritePDBASCII(file_name.c_str());
               multi_residue_torsion_fit_map(imol, moving_mol, *xmap, avoid_these_atoms, n_trials, &geom);
               // file_name = "mrtf-post-" + std::to_string(ii) + ".pdb";
               // moving_mol->WritePDBASCII(file_name.c_str());

               if (mode == 3) {

                  std::vector<std::pair<bool,mmdb::Residue *> > residues =
                     make_refinement_residues(residue_specs, moving_mol);
                  std::vector<mmdb::Link> links;
                  std::vector<atom_spec_t> fixed_atom_specs;
                  double torsions_weight = 6.0;

                  // debug_coords(moving_mol);
                  refine_direct(residues, moving_mol, imol, links, level,
                                parent_spec, parent_type, new_res_spec, new_res_type, new_link,
                                geom, fixed_atom_specs,
                                *xmap, map_weight, torsions_weight);

                  // debug_coords(moving_mol);
                  if (false) {
                     std::string file_name = "post-refine-" + std::to_string(ii) + ".pdb";
                     moving_mol->WritePDBASCII(file_name.c_str());
                  }

                  mmdb::Residue *new_residue_p = util::get_residue(new_res_spec, mol);
                  if (new_residue_p) {
                     bool well_fitting = is_well_fitting(new_residue_p, mol, *xmap, geom);
                     // std::cout << "INFO:: new residue " << new_res_spec << " is-well-fitting-status: "
                     //           << well_fitting << std::endl;
                     logger.log(log_t::INFO, "new residue", new_res_spec.format(),
                                "is-well_fitting-status:", well_fitting);
                  }

                  replace_coords(moving_mol, mol);
               }
            }
            // OK, we we want to add the atoms of moving_mol into mol
            replace_coords(moving_mol, mol);
            delete moving_mol;
         }
      }
   }
   return new_res_spec;
}

//! do the thing.
//! res-pair is new-link-type and new-res-type
coot::residue_spec_t
coot::cho::add_linked_residue_add_cho_function(atom_selection_container_t *asc,
                                               int imol,
                                               const coot::residue_spec_t &parent,
                                               const std::pair<std::string, std::string> &new_link_types,
                                               unsigned int level,
                                               float new_atoms_b_factor,
                                               coot::protein_geometry &geom,
                                               const clipper::Xmap<float> *xmap,
                                               float map_weight) { // can be null

   int mode = 3; // mode is either 1: add  2: add and fit  3: add, fit and refine

   const std::string &new_link     = new_link_types.first;
   const std::string &new_res_type = new_link_types.second;

   mmdb::Residue *residue_p = util::get_residue(parent, asc->mol);
   coot::glyco_tree_t t(residue_p, asc->mol, &geom); // geom is not const
   std::vector<mmdb::Residue *> tree_residues = t.residues(parent);
   coot::residue_spec_t new_res_spec = add_linked_residue(asc, imol, parent, new_link_types, level,
                                                          new_atoms_b_factor, mode, geom, xmap, map_weight);

   if (false) {
      std::cout << "===== add_linked_residue() returned " << new_res_spec << std::endl;
      // if mode is 3, then this will be a post-refined model.
      std::string fn = "added-" + new_link + "-" + new_res_type + ".pdb";
      asc->mol->WritePDBASCII(fn.c_str());
   }

   return new_res_spec;
}


