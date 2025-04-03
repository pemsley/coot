
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
                           clipper::Xmap<float> &xmap,
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
            std::cout<<"INFO:: not adding water because of overlap\n"<<std::endl;
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
      std::cout << "WARNING:: atom " << spec_1 << " not found - abandoning LINK addition " << std::endl;
   } else {
      if (! at_2) {
	 std::cout << "WARNING:: atom " << spec_1 << " not found - abandoning LINK addition " << std::endl;
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
         std::cout << "WARNING:: " << rte.what() << std::endl;
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
                     if (false)
                        std::cout << "--------------- replacing coords for atom "
                                  << atom_spec_t(at) << " "
                                  << at_mol->x << " " << at_mol->y << " " << at_mol->z << " "
                                  <<     at->x << " " <<     at->y << " " <<     at->z << " "
                                  << std::endl;
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



coot::residue_spec_t
coot::cho::add_linked_residue(mmdb::Manager *mol,
                              int imol,
                              const residue_spec_t &parent,
                              const std::pair<std::string, std::string> &new_link_types,
                              float new_atoms_b_factor,
                              int mode,
                              coot::protein_geometry &geom,
                              const clipper::Xmap<float> *xmap) {

   int n_trials = 5000;
   int cif_dictionary_read_number = 60; // mergh. Pass a reference to this?

   const std::string &new_link     = new_link_types.first;
   const std::string &new_res_type = new_link_types.second;
   if (geom.have_dictionary_for_residue_type_no_dynamic_add(new_res_type, imol)) {
   } else {
      int status = geom.try_dynamic_add(new_res_type, cif_dictionary_read_number);
      if (status == 0) { // fail
         std::cout << "WARNING:: failed to add dictionary for " << new_res_type << std::endl;
         std::cout << "WARNING:: add_linked_residue() stops here " << std::endl;
         return residue_spec_t();
      }
   }
   cif_dictionary_read_number++;
   residue_spec_t new_res_spec = add_linked_residue_by_atom_torsions(mol, parent, new_link_types, geom, new_atoms_b_factor);
   // delete extra restraints for new_res_spec - hmm
   if (! new_res_spec.unset_p()) {
      if (mode == 2 || mode == 3) {
         if (xmap) {

            std::vector<residue_spec_t > residue_specs;
            residue_specs.push_back(parent);
            residue_specs.push_back(new_res_spec);
            mmdb::Manager *moving_mol = util::create_mmdbmanager_from_residue_specs(residue_specs, mol);

            int n_rounds_of_fit_and_refine = 3;
            std::vector<std::pair<bool, clipper::Coord_orth> > avoid_these_atoms;
            for (int ii=0; ii<n_rounds_of_fit_and_refine; ii++) {
               std::string file_name = "mrtf-pre-" + std::to_string(ii) + ".pdb";
               moving_mol->WritePDBASCII(file_name.c_str());
               multi_residue_torsion_fit_map(imol, moving_mol, *xmap, avoid_these_atoms, n_trials, &geom);
               file_name = "mrtf-post-" + std::to_string(ii) + ".pdb";
               moving_mol->WritePDBASCII(file_name.c_str());

               if (mode == 3) {
                  // refine the residue_specs. Not sure that I can do that here.
               }
            }
            // OK, we we want to add the atoms of moving_mol into mol
            replace_coords(moving_mol, mol);
         }
      }
   }
   return new_res_spec;
}

//! do the thing.
//! res-pair is new-link-type and new-res-type
coot::residue_spec_t
coot::cho::add_linked_residue_add_cho_function(mmdb::Manager *mol, int imol,
                                               const coot::residue_spec_t &parent,
                                               const std::pair<std::string, std::string> &new_link_types,
                                               float new_atoms_b_factor,
                                               coot::protein_geometry &geom,
                                               const clipper::Xmap<float> *xmap) { // can be null

   int mode = 2; // mode is either 1: add  2: add and fit  3: add, fit and refine

   const std::string &new_link     = new_link_types.first;
   const std::string &new_res_type = new_link_types.second;

   mmdb::Residue *residue_p = util::get_residue(parent, mol);
   coot::glyco_tree_t t(residue_p, mol, &geom); // geom is not const
   std::vector<mmdb::Residue *> tree_residues = t.residues(parent);
   coot::residue_spec_t new_res_spec = add_linked_residue(mol, imol, parent, new_link_types,
                                                          new_atoms_b_factor, mode, geom, xmap);

   if (true) {
      std::string fn = "added-" + new_link + "-" + new_res_type + ".pdb";
      mol->WritePDBASCII(fn.c_str());
   }

   return new_res_spec;
}


