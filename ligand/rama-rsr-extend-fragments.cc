
#include "utils/split-indices.hh"
#include "coot-utils/coot-map-utils.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "rama-rsr-extend-fragments.hh"
#include "residue_by_phi_psi.hh"
#include "ideal/simple-restraint.hh"
#include "new-residue-by-3-phi-psi.hh"

std::pair<int, float> get_sum_of_density_for_residue(mmdb::Residue *residue_p, const clipper::Xmap<float> &xmap) {

   mmdb::Atom **residue_atoms = 0;
   int n_residue_atoms = 0;
   float density_sum = 0.0f;
   int n_atoms = 0;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int iat=0; iat<n_residue_atoms; iat++) {
      mmdb::Atom *at = residue_atoms[iat];
      if (! at->isTer()) {
         clipper::Coord_orth co(coot::co(at));
         float d = coot::util::density_at_point(xmap, co);
         n_atoms++;
         density_sum += d;
      }
   }
   return std::make_pair(density_sum, n_atoms);
}

float get_average_density_per_atom(mmdb::Manager *mol, const clipper::Xmap<float> &xmap) {

   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      float sum = 0.0f;
      unsigned int n_atoms = 0;
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int n_res = chain_p->GetNumberOfResidues();
         for (int ires=0; ires<n_res; ires++) {
            mmdb::Residue *residue_p = chain_p->GetResidue(ires);
            if (residue_p) {
               auto dp = get_sum_of_density_for_residue(residue_p, xmap);
               sum     += dp.first;
               n_atoms += dp.second;
            }
         }
      }
      if (n_atoms > 0) {
         float av = sum/static_cast<float>(n_atoms);
         return av;
      }
   }
   std::cout << "ERROR:: no atoms in get_average_density_per_atom() " << std::endl;
   return 0.0f;
   
}



void
rama_rsr_extend_fragments(mmdb::Manager *mol, const clipper::Xmap<float> &xmap, ctpl::thread_pool  *thread_pool_p, unsigned int n_threads,
                          float weight, unsigned int n_phi_psi_trials, const coot::protein_geometry &geom,
                          unsigned int *update_count) {

   // return the sum and the count
   auto get_density_sum_for_new_residues = [&xmap] (const coot::minimol::molecule &m) { // capture xmap by reference

                                              float sum = 0.0;
                                              unsigned int n_atoms = 0;

                                              for(unsigned int ifrag=0; ifrag<m.fragments.size(); ifrag++) {
                                                 for(int ires=m[ifrag].min_res_no(); ires<=m[ifrag].max_residue_number(); ires++) {
                                                    for (unsigned int iat=0; iat<m[ifrag][ires].atoms.size(); iat++) {
                                                       const clipper::Coord_orth atom_pos(m[ifrag][ires][iat].pos);
                                                       sum += coot::util::density_at_point(xmap, atom_pos);
                                                       n_atoms++;
                                                    }
                                                 }
                                              }
                                              return std::make_pair(sum, n_atoms);
                                           };

   
   // return the sum and the count - this is the the 3 residue case - with occupancy weighting
   auto get_density_score_for_new_residues_frag = [&xmap] (const coot::minimol::fragment &frag) { // capture xmap by reference

                                              float sum = 0.0;
                                              float sum_w = 0.0;
                                              unsigned int n_atoms = 0;

                                              for(int ires=frag.min_res_no(); ires<=frag.max_residue_number(); ires++) {
                                                 for (unsigned int iat=0; iat<frag[ires].atoms.size(); iat++) {
                                                    const clipper::Coord_orth atom_pos(frag[ires][iat].pos);
                                                    float f = coot::util::density_at_point(xmap, atom_pos);
                                                    float w = frag[ires][iat].occupancy;
                                                    sum += f * w;
                                                    sum_w += w;
                                                    n_atoms++;
                                                 }
                                              }
                                              return std::make_pair(sum/sum_w, n_atoms);
                                           };


   auto add_residues_to_chain = [] (const coot::minimol::molecule &m, mmdb::Chain *chain_p, int seq_num_of_anchor_residue,
                                    const std::string &terminus_type, unsigned int n_residues_to_add=2) {

                                   bool status = true;

                                   for(unsigned int ifrag=0; ifrag<m.fragments.size(); ifrag++) {
                                      std::string chain_id = m[ifrag].fragment_id;
                                      if (terminus_type == "N") {
                                         if (n_residues_to_add == 2) {
                                            mmdb::Residue *residue_minus_1_p = m[ifrag][seq_num_of_anchor_residue-1].make_residue();
                                            mmdb::Residue *residue_minus_2_p = m[ifrag][seq_num_of_anchor_residue-2].make_residue();
                                            if (residue_minus_1_p && residue_minus_2_p) {
                                               chain_p->InsResidue(residue_minus_1_p, 0);
                                               chain_p->InsResidue(residue_minus_2_p, 0);
                                            } else {
                                               status = false;
                                               std::cout << "in add_residues_to_chain() N failed to get residues " << residue_minus_1_p << "  " << residue_minus_2_p
                                                         << std::endl;
                                            }
                                         }

                                         if (n_residues_to_add == 3) {
                                            mmdb::Residue *residue_minus_1_p = m[ifrag][seq_num_of_anchor_residue-1].make_residue();
                                            mmdb::Residue *residue_minus_2_p = m[ifrag][seq_num_of_anchor_residue-2].make_residue();
                                            mmdb::Residue *residue_minus_3_p = m[ifrag][seq_num_of_anchor_residue-3].make_residue();
                                            if (residue_minus_1_p && residue_minus_2_p && residue_minus_3_p) {
                                               chain_p->InsResidue(residue_minus_1_p, 0);
                                               chain_p->InsResidue(residue_minus_2_p, 0);
                                               chain_p->InsResidue(residue_minus_3_p, 0);
                                            } else {
                                               status = false;
                                               std::cout << "in add_residues_to_chain() N failed to get residues " << residue_minus_1_p << "  "
                                                         << residue_minus_2_p << " " << residue_minus_3_p << std::endl;
                                            }
                                         }
                                      }

                                      if (terminus_type == "C") {
                                         if (n_residues_to_add == 2) {
                                            mmdb::Residue *residue_plus_1_p = m[ifrag][seq_num_of_anchor_residue+1].make_residue();
                                            mmdb::Residue *residue_plus_2_p = m[ifrag][seq_num_of_anchor_residue+2].make_residue();
                                            if (residue_plus_1_p && residue_plus_2_p) {
                                               chain_p->AddResidue(residue_plus_1_p);
                                               chain_p->AddResidue(residue_plus_2_p);
                                               std::cout << "debug:: seq_num_of_anchor_residue " << seq_num_of_anchor_residue
                                                         << " add_residue to chain " << coot::residue_spec_t(residue_plus_1_p) << std::endl;
                                            } else {
                                               status = false;
                                               std::cout << "in add_residues_to_chain() C failed to get residues " << residue_plus_1_p << "  " << residue_plus_2_p
                                                         << std::endl;
                                            }
                                         }

                                         if (n_residues_to_add == 3) {
                                            mmdb::Residue *residue_plus_1_p = m[ifrag][seq_num_of_anchor_residue+1].make_residue();
                                            mmdb::Residue *residue_plus_2_p = m[ifrag][seq_num_of_anchor_residue+2].make_residue();
                                            mmdb::Residue *residue_plus_3_p = m[ifrag][seq_num_of_anchor_residue+3].make_residue();
                                            if (residue_plus_1_p && residue_plus_2_p && residue_plus_3_p) {
                                               chain_p->AddResidue(residue_plus_1_p);
                                               chain_p->AddResidue(residue_plus_2_p);
                                               chain_p->AddResidue(residue_plus_3_p);
                                               if (false)
                                                  std::cout << "debug:: add_residues_to_chain() seq_num_of_anchor_residue " << seq_num_of_anchor_residue
                                                            << " add_residues "
                                                            << coot::residue_spec_t(residue_plus_1_p) << " "
                                                            << coot::residue_spec_t(residue_plus_2_p) << " "
                                                            << coot::residue_spec_t(residue_plus_3_p) << std::endl;
                                            } else {
                                               status = false;
                                               std::cout << "in add_residues_to_chain() C failed to get residues " << residue_plus_1_p << "  " << residue_plus_2_p
                                                         << std::endl;
                                            }
                                         }
                                      }
                                   }
                                   return status;
                                };


   // this is a copy of res-tracer make_CB_ideal_pos. You might consider consolidating them
   //
   // pass the chain_id for debugging.
   auto calculate_CB_ideal_pos = [] (const coot::minimol::residue &res, const std::string &chain_id) {
                                    bool found = false;
                                    clipper::Coord_orth pos;
                                    // don't add one if there is one already there.
                                    std::pair<bool, coot::minimol::atom> CB = res.get_atom(" CB ");
                                    if (! CB.first) {
                                       auto CA = res.get_atom(" CA ");
                                       auto C  = res.get_atom(" C  ");
                                       auto N  = res.get_atom(" N  ");
                                       if (CA.first) {
                                          if (C.first) {
                                             if (N.first) {
                                                clipper::Coord_orth C_to_N = N.second.pos - C.second.pos;
                                                clipper::Coord_orth C_to_N_mid_point(0.5 * (N.second.pos + C.second.pos));
                                                clipper::Coord_orth CA_to_CN_mid_point = C_to_N_mid_point - CA.second.pos;
                                                clipper::Coord_orth CA_to_CN_mid_point_uv(CA_to_CN_mid_point.unit());
                                                clipper::Coord_orth perp(clipper::Coord_orth::cross(C_to_N, CA_to_CN_mid_point));
                                                clipper::Coord_orth perp_uv(perp.unit());
                                                // guess and fiddle these - good enough
                                                clipper::Coord_orth CB_pos(CA.second.pos + 1.21 * perp_uv - 0.95 * CA_to_CN_mid_point_uv);
                                                pos = CB_pos;
                                                found = true;
                                             } else {
                                                std::cout << "INFO:: calculate_CB_ideal_pos(): sad residue " << res << " in chain " << chain_id << ", has no N " << std::endl;
                                             }
                                          } else {
                                             std::cout << "INFO:: calculate_CB_ideal_pos(): sad residue " << res << " in chain " << chain_id << ", has no C " << std::endl;
                                          }
                                       } else {
                                          std::cout << "INFO:: calculate_CB_ideal_pos(): sad residue " << res << " in chain " << chain_id << ", has no CA " << std::endl;
                                       }
                                    } else {
                                       std::cout << "INFO:: calculate_CB_ideal_pos(): residue " << res << " in chain " << chain_id << " " << res.seqnum
                                                 << ", already has a CB " << std::endl;
                                    }
                                    return std::make_pair(found, pos);
                                 };

   auto add_CB_to_residue_maybe = [calculate_CB_ideal_pos]
      (coot::minimol::molecule *mmol_p, const clipper::Xmap<float> &xmap, float average_density_per_atom_for_molecule) {

                                     float crit = 0.0;
                                     coot::minimol::molecule &m(*mmol_p);
                                     for(unsigned int ifrag=0; ifrag<m.fragments.size(); ifrag++) {
                                        const std::string &chain_id(m[ifrag].fragment_id);
                                        for(int ires=m[ifrag].min_res_no(); ires<=m[ifrag].max_residue_number(); ires++) {
                                           coot::minimol::residue &residue = m[ifrag][ires];
                                           if (! residue.is_empty()) {
                                              std::pair<bool, clipper::Coord_orth > pt = calculate_CB_ideal_pos(residue, chain_id);
                                              if (pt.first) {
                                                 float dv = coot::util::density_at_point(xmap, pt.second);
                                                 std::cout << "debug:: in add_CB_to_residue_maybe() for " << m[ifrag].fragment_id << " " << ires
                                                           << " dv is " << dv << std::endl;
                                                 if (dv > crit * average_density_per_atom_for_molecule) {
                                                    // add it
                                                    coot::minimol::atom cb(" CB ", " C", pt.second, "", 20.0f);
                                                    residue.addatom(cb);
                                                 } else {
                                                    std::cout << "INFO:: Density too wispy to add CB for chain-id " << chain_id << " ires " << ires << " "
                                                              << " av for molecule " << average_density_per_atom_for_molecule
                                                              << " cut-off " << 0.7 * average_density_per_atom_for_molecule
                                                              << " this CB-density " << dv << std::endl;
                                                 }
                                              } else {
                                                 std::cout << "WARNING:: in add_CB_to_residue_maybe() failed to calculate ideal CB pos "
                                                           << chain_id << " " << ires << std::endl;
                                              }
                                           } else {
                                              // std::cout << "DEBUG:: residue in " << chain_id << " ires " << ires << " residue.seqnum: " << residue.seqnum
                                              // << " is empty" << std::endl;
                                           }
                                        }
                                     }
                                  };

   auto refine_triple = [] (int seqnum, mmdb::Chain *chain_p, mmdb::Manager *mol, const std::string &terminus_type, const clipper::Xmap<float> &xmap,
                            const coot::protein_geometry &geom, ctpl::thread_pool *thread_pool_p, unsigned int n_threads, float weight) {

                           std::vector<std::pair<bool, mmdb::Residue *> > residues;
                           std::vector<int> residue_seqnums = { seqnum-2, seqnum-1, seqnum};
                           if (terminus_type == "C")
                              residue_seqnums = { seqnum, seqnum+1, seqnum+2};
                           std::string chain_id(chain_p->GetChainID());
                           for (unsigned int i=0; i<residue_seqnums.size(); i++) {
                              mmdb::Residue *residue_p = chain_p->GetResidue(residue_seqnums[i], ""); //  1st arg is for seqno
                              if (residue_p) {
                                 residues.push_back(std::make_pair(false, residue_p));
                              }
                           }
                           if (residues.size() == 3) {

                              std::cout << "debug:: in refine_triple() seqno is " << seqnum << " and residue numbers "
                                        << residues[0].second->GetSeqNum() << " "
                                        << residues[1].second->GetSeqNum() << " "
                                        << residues[2].second->GetSeqNum() << std::endl;
                              std::vector<mmdb::Link> links;
                              std::vector<coot::atom_spec_t> fixed_atom_specs;
                              coot::restraint_usage_Flags flags = coot::TYPICAL_RESTRAINTS;
                              coot::restraints_container_t restraints(residues, links, geom, mol, fixed_atom_specs, &xmap);
                              restraints.thread_pool(thread_pool_p, n_threads);
                              restraints.set_quiet_reporting();
                              coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
                              bool do_internal_torsions = false;
                              restraints.add_map(weight);
                              int imol = 0;
                              restraints.make_restraints(imol, geom, flags, do_internal_torsions, false, 0, 0, true, true, false, pseudos);
                              restraints.minimize(flags);
                           }
                        };

   auto refine_quad = [] (int seqnum, mmdb::Chain *chain_p, mmdb::Manager *mol, const std::string &terminus_type, const clipper::Xmap<float> &xmap,
                          const coot::protein_geometry &geom, ctpl::thread_pool *thread_pool_p, unsigned int n_threads, float weight) {

                           std::vector<std::pair<bool, mmdb::Residue *> > residues;
                           std::vector<int> residue_seqnums = { seqnum-3, seqnum-2, seqnum-1, seqnum};
                           if (terminus_type == "C")
                              residue_seqnums = { seqnum, seqnum+1, seqnum+2, seqnum+3};
                           std::string chain_id(chain_p->GetChainID());
                           for (unsigned int i=0; i<residue_seqnums.size(); i++) {
                              mmdb::Residue *residue_p = chain_p->GetResidue(residue_seqnums[i], ""); //  1st arg is for seqno
                              if (residue_p) {
                                 residues.push_back(std::make_pair(false, residue_p));
                              }
                           }
                           if (residues.size() == 4) {

                              std::cout << "debug:: in refine_quad() seqno is " << seqnum << " and residue numbers "
                                        << residues[0].second->GetSeqNum() << " "
                                        << residues[1].second->GetSeqNum() << " "
                                        << residues[2].second->GetSeqNum() << " "
                                        << residues[3].second->GetSeqNum() << std::endl;
                              auto tp_0 = std::chrono::high_resolution_clock::now();
                              std::vector<mmdb::Link> links;
                              std::vector<coot::atom_spec_t> fixed_atom_specs;
                              coot::restraint_usage_Flags flags = coot::TYPICAL_RESTRAINTS;
                              coot::restraints_container_t restraints(residues, links, geom, mol, fixed_atom_specs, &xmap);
                              restraints.thread_pool(thread_pool_p, n_threads);
                              restraints.set_quiet_reporting();
                              coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
                              bool do_internal_torsions = false;
                              restraints.add_map(weight);
                              int imol = 0;
                              restraints.make_restraints(imol, geom, flags, do_internal_torsions, false, 0, 0, true, true, false, pseudos);
                              restraints.minimize(flags);
                              auto tp_1 = std::chrono::high_resolution_clock::now();
                              auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
                              // std::cout << "Timings: refine_quad(): " << chain_id << " " << seqnum << " " << d10 << " milliseconds" << std::endl;
                              // std::cout << "INFO:: in refine_quad() refinement finished for " << residues.size() << " residues " << std::endl;
                              if (false) {
                                 std::string file_name = "refine-quad-term-" + terminus_type + "-chain-" + chain_id + "-" + std::to_string(seqnum) + ".pdb";
                                 mol->WritePDBASCII(file_name.c_str());
                              }
                           } else {
                              std::cout << "ERROR:: in refine_quad() residues.size() was " << residues.size() << std::endl;
                           }
                        };

   // refine this quad of residues without considering other residues (for non-bonded contacts)
   auto refine_isolated_quad = [] (int seqnum, mmdb::Chain *chain_p, mmdb::Manager *mol, const std::string &terminus_type, const clipper::Xmap<float> &xmap,
                                   const coot::protein_geometry &geom, ctpl::thread_pool *thread_pool_p, unsigned int n_threads, float weight) {

                                  // 1: Make a new molecule from the specs
                                  // 2: Refine that chain
                                  // 3: copy back the positions of the refined atoms into the residues of the original (input) chain

                                  // Maybe I will need to take account of chains that don't have 4 residues?
                                  // Maybe I willl need to add a non-refine "anchor" residue

                                  std::vector<coot::residue_spec_t> specs;
                                  std::vector<mmdb::Residue *> original_residues;
                                  std::vector<int> residue_seqnums = { seqnum-3, seqnum-2, seqnum-1, seqnum};
                                  std::string chain_id(chain_p->GetChainID());
                                  if (terminus_type == "C") residue_seqnums = { seqnum, seqnum+1, seqnum+2, seqnum+3};
                                  for (unsigned int i=0; i<residue_seqnums.size(); i++) {
                                     mmdb::Residue *residue_p = chain_p->GetResidue(residue_seqnums[i], ""); //  1st arg is for seqno
                                     if (residue_p) {
                                        specs.push_back(coot::residue_spec_t(residue_p));
                                        original_residues.push_back(residue_p);
                                     }
                                  }
                                  if (specs.size() == 4) {
                                     mmdb::Manager *refmol = coot::util::create_mmdbmanager_from_residue_specs(specs, mol);
                                     std::vector<std::pair<bool, mmdb::Residue *> > ref_residues;
                                     for (unsigned int i=0; i<4; i++) {
                                        mmdb::Residue *residue_p = coot::util::get_residue(specs[i], refmol);
                                        if (residue_p) {
                                           ref_residues.push_back(std::make_pair(false, residue_p));
                                        } else {
                                           std::cout << "ERROR:: in refine_isolated_quad() failed to extract residue " << specs[i] << std::endl;
                                        }
                                     }
                                     if (ref_residues.size() == 4) {
                                        auto tp_0 = std::chrono::high_resolution_clock::now();
                                        std::vector<mmdb::Link> links;
                                        std::vector<coot::atom_spec_t> fixed_atom_specs;
                                        coot::restraint_usage_Flags flags = coot::TYPICAL_RESTRAINTS;
                                        coot::restraints_container_t restraints(ref_residues, links, geom, refmol, fixed_atom_specs, &xmap);
                                        restraints.thread_pool(thread_pool_p, n_threads);
                                        restraints.set_quiet_reporting();
                                        coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
                                        bool do_internal_torsions = false;
                                        restraints.add_map(weight);
                                        int imol = 0;
                                        restraints.make_restraints(imol, geom, flags, do_internal_torsions, false, 0, 0, true, true, false, pseudos);
                                        restraints.minimize(flags, 500); // default 1000 steps
                                        auto tp_1 = std::chrono::high_resolution_clock::now();
                                        auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
                                        // std::cout << "Timings: refine_isolated_quad(): " << chain_id << " " << seqnum << " " << d10 << " milliseconds" << std::endl;
                                        // std::cout << "INFO:: in refine_isolated_quad() refinement finished for " << ref_residues.size() << " residues " << std::endl;

                                        // copy atoms back into original chain
                                        for (unsigned int i=0; i<4; i++) {
                                           mmdb::Residue *original_residue_p = original_residues[i];
                                           mmdb::Residue *refined_residue_p  = ref_residues[i].second;
                                           int n_atoms_o =  original_residue_p->GetNumberOfAtoms();
                                           int n_atoms_r =   refined_residue_p->GetNumberOfAtoms();
                                           if (n_atoms_o == n_atoms_r) {
                                              for (int iat=0; iat<n_atoms_o; iat++) {
                                                 mmdb::Atom *atom_o = original_residue_p->GetAtom(iat);
                                                 mmdb::Atom *atom_r =  refined_residue_p->GetAtom(iat);
                                                 atom_o->x = atom_r->x;
                                                 atom_o->y = atom_r->y;
                                                 atom_o->z = atom_r->z;
                                              }
                                           }
                                        }
                                        
                                     } else {
                                        std::cout << "ERROR:: in refine_isolated_quad() ref residues size " << ref_residues.size() << std::endl;
                                     }
                                     delete refmol;
                                  } else {
                                     std::cout << "ERROR:: in refine_isolated_quad(): specs size " << specs.size() << std::endl;
                                  }
                               };

   // delete the added residue (because it was bad)
   auto delete_extra_residue = [] (int seqnum, mmdb::Chain *chain_p, mmdb::Manager *mol, const std::string &terminus_type) {
                                  // Dangerous beacuse no checking?
                                  if (terminus_type == "N") {
                                     // mmdb::Residue *residue_p = chain_p->GetResidue(0);
                                     // std::cout << "debug:: in delete_extra_residue() deleting residue " << coot::residue_spec_t(residue_p)  << std::endl;
                                     chain_p->DeleteResidue(0);
                                     mol->FinishStructEdit();
                                  }
                                  if (terminus_type == "C") {
                                     int n_res = chain_p->GetNumberOfResidues();
                                     int residue_index = n_res -1;
                                     chain_p->DeleteResidue(residue_index);
                                     mol->FinishStructEdit();
                                  }
                               };

   // delete the 2 added residue beyond the first one (build 3, keep 1)
   auto delete_extra_residues = [] (int seqnum, mmdb::Chain *chain_p, mmdb::Manager *mol, const std::string &terminus_type) {
                                   // Dangerous beacuse no checking?
                                   if (terminus_type == "N") {
                                      // std::cout << "debug:: in delete_extra_residue() deleting residue " << coot::residue_spec_t(residue_p)  << std::endl;
                                      chain_p->DeleteResidue(1);
                                      chain_p->DeleteResidue(0);
                                      mol->FinishStructEdit();
                                   }
                                   if (terminus_type == "C") {
                                      int n_res = chain_p->GetNumberOfResidues();
                                      int residue_index = n_res -1;
                                      chain_p->DeleteResidue(residue_index);
                                      chain_p->DeleteResidue(residue_index-1);
                                      mol->FinishStructEdit();
                                   }
                                };

   // This function is in res-tracer.cc
   //
   // twisted trans or cis
   auto peptide_is_twisted = [] (mmdb::Residue *residue_with_CO, mmdb::Residue *residue_with_N, double deformation_limit_deg = 20.0) {

                                bool status = false;
                                mmdb::Atom *CA_1 = residue_with_CO->GetAtom(" CA ");
                                mmdb::Atom *C_1  = residue_with_CO->GetAtom(" C  ");
                                mmdb::Atom *N_2  = residue_with_N->GetAtom(" N  ");
                                mmdb::Atom *CA_2 = residue_with_N->GetAtom(" CA ");
                                if (CA_1 && C_1 && N_2 && CA_2) {
                                   clipper::Coord_orth ca_1_pt = coot::co(CA_1);
                                   clipper::Coord_orth c_1_pt  = coot::co(C_1);
                                   clipper::Coord_orth n_2_pt  = coot::co(N_2);
                                   clipper::Coord_orth ca_2_pt = coot::co(CA_2);
                                   double torsion = clipper::Coord_orth::torsion(ca_1_pt, c_1_pt, n_2_pt, ca_2_pt);
                                   double torsion_deg = clipper::Util::rad2d(torsion);
                                   if (torsion_deg > (-180.0 + deformation_limit_deg))
                                      if (torsion_deg < (180.0 - deformation_limit_deg))
                                         status = true;
                                   if (false)
                                      if (status)
                                         std::cout << "debug:: twisted peptide " << coot::residue_spec_t(residue_with_CO) << " "
                                                   << coot::residue_spec_t(residue_with_N) << " " << torsion_deg << std::endl;
                                } else {
                                   std::cout << "ERROR:: peptide_is_twisted(): missing atoms torsion " << std::endl;
                                }
                                return status;
                             };

   auto new_residue_is_deformed = [peptide_is_twisted] (mmdb::Residue *residue_p, const std::string &terminus_type, mmdb::Manager *mol) {

                                     double deformation_limit_deg = 55.0;

                                     bool status = false;
                                     coot::residue_spec_t anchor_res_spec(residue_p);
                                     coot::residue_spec_t new_residue_spec = anchor_res_spec.next();
                                     if (terminus_type == "N")
                                        new_residue_spec = anchor_res_spec.previous();
                                     mmdb::Residue *added_residue_p = coot::util::get_residue(new_residue_spec, mol);
                                     if (added_residue_p) {
                                        if (terminus_type == "N") {
                                           bool twisted = peptide_is_twisted(added_residue_p, residue_p, deformation_limit_deg);
                                           if (twisted) status = true;
                                        } else {
                                           bool twisted = peptide_is_twisted(residue_p, added_residue_p, deformation_limit_deg);
                                           if (twisted) status = true;
                                        }
                                     } else {
                                        std::cout << "DEBUG:: new_residue_is_deformed() added_residue_p not found, hence failure "
                                                  << "added_residue_spec: " << new_residue_spec << std::endl;
                                        status = true; // can't find it - that can't be good!
                                     }
                                     return status;
                                  };

   auto new_residue_crashing_into_other_chain = [] (mmdb::Residue *residue_p, const std::string &terminus_type, mmdb::Manager *mol) {
                                                   bool status = false;

                                                   coot::residue_spec_t anchor_res_spec(residue_p);
                                                   coot::residue_spec_t new_residue_spec = anchor_res_spec.next();
                                                   if (terminus_type == "N")
                                                      new_residue_spec = anchor_res_spec.previous();
                                                   mmdb::Residue *added_residue_p = coot::util::get_residue(new_residue_spec, mol);
                                                   if (added_residue_p) {

                                                      mmdb::realtype local_dist_max = 3.5; // hmmm... tricky
                                                      mmdb::Contact *pscontact = NULL;
                                                      int n_contacts = 0;
                                                      long i_contact_group = 1;
                                                      mmdb::mat44 my_matt;
                                                      for (int i=0; i<4; i++)
                                                         for (int j=0; j<4; j++)
                                                            my_matt[i][j] = 0.0;
                                                      for (int i=0; i<4; i++) my_matt[i][i] = 1.0;
                                                      //

                                                      int SelHnd_1 = mol->NewSelection(); // d
                                                      int SelHnd_2 = mol->NewSelection(); // d

                                                      mol->SelectAtoms(SelHnd_2, 0, "*", mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", " CA ", " C", "");
                                                      new_residue_spec.select_atoms(mol, SelHnd_1, mmdb::SKEY_NEW);
                                                      mmdb::Atom **atom_selection_1 = 0;
                                                      mmdb::Atom **atom_selection_2 = 0;
                                                      int n_selected_atoms_1 = 0;
                                                      int n_selected_atoms_2 = 0;
                                                      mol->GetSelIndex(SelHnd_1, atom_selection_1, n_selected_atoms_1);
                                                      mol->GetSelIndex(SelHnd_2, atom_selection_2, n_selected_atoms_2);

                                                      mol->SeekContacts(atom_selection_1, n_selected_atoms_1,
                                                                        atom_selection_2, n_selected_atoms_2,
                                                                        0.0, local_dist_max, // min, max distances
                                                                        1,        // seqDist 0 -> in same res also
                                                                        pscontact, n_contacts,
                                                                        0, &my_matt, i_contact_group);

                                                      if (n_contacts > 0) {
                                                         // std::cout << "debug:: Selection had " << n_contacts << " CA-CA contacts" << std::endl;
                                                         if (pscontact) {
                                                            unsigned int n_clash = 0;
                                                            for (int i=0; i<n_contacts; i++) {
                                                               mmdb::Atom *at_1 = atom_selection_1[pscontact[i].id1];
                                                               mmdb::Atom *at_2 = atom_selection_2[pscontact[i].id2];
                                                               if (at_1->GetResidue() == at_2->GetResidue()) {
                                                                  // ignore
                                                               } else {
                                                                  if (at_1->GetChain() == at_2->GetChain()) {
                                                                     int res_no_delta = at_1->GetResidue()->GetSeqNum() - at_2->GetResidue()->GetSeqNum();
                                                                     if (res_no_delta < 2) {
                                                                        // ignore interactions to bonded neighbour
                                                                     } else {
                                                                        n_clash++;
                                                                     }
                                                                  } else {
                                                                     // normal clashing case
                                                                     n_clash++;
                                                                  }
                                                               }
                                                            }
                                                            if (n_clash > 2) // will need validation
                                                               status = true;

                                                            if (false)
                                                               std::cout << "debug:: new_residue_crashing_into_other_chain() " << coot::residue_spec_t(residue_p)
                                                                         << " n_clash " << n_clash << std::endl;
                                                         }
                                                      }
                                                      mol->DeleteSelection(SelHnd_1);
                                                      mol->DeleteSelection(SelHnd_2);
                                                   }
                                                   // return status; // reinstate the correct return value at some stage.
                                                   return false;
                                                };

   // return a bool - for success status (residue added)
   auto extend_this_chain = [get_density_sum_for_new_residues, add_residues_to_chain, add_CB_to_residue_maybe,
                             refine_triple, delete_extra_residue, new_residue_is_deformed, new_residue_crashing_into_other_chain,
                             thread_pool_p, n_threads, xmap] (mmdb::Chain *chain_p, mmdb::Manager *mol,
                                                              const std::string &terminus_type,
                                                              float average_density_per_atom_for_molecule,
                                                              float weight, unsigned int n_phi_psi_trials,
                                                              const coot::protein_geometry &geom) {

                               bool status = false; // initially no residue added.
                               int n_res = chain_p->GetNumberOfResidues();
                               if (n_res > 0) {
                                  std::string chain_id(chain_p->GetChainID());
                                  int residue_index = 0;
                                  int n_residues_in_chain =  chain_p->GetNumberOfResidues();
                                  if (terminus_type == "C") {
                                     residue_index = n_residues_in_chain -1;
                                  }
                                  mmdb::Residue *residue_p = chain_p->GetResidue(residue_index);
                                  std::cout << "debug:: in extend_this_chain " << chain_id << " terminus type " << terminus_type
                                            << " n_residues_in_chain " << n_residues_in_chain <<  " residue_p with index : "
                                            << residue_index << " " << coot::residue_spec_t(residue_p) << std::endl;
                                  if (residue_p) {

                                     std::string  residue_type = "ALA";
                                     int residue_seqnum = residue_p->GetSeqNum();

                                     float b_factor = 20.0;
                                     coot::residue_by_phi_psi rphipsi(terminus_type, residue_p, chain_id, residue_type, b_factor);
                                     rphipsi.thread_pool(thread_pool_p, n_threads);
                                     coot::minimol::molecule m = rphipsi.best_fit_phi_psi(n_phi_psi_trials, false, true, xmap);
                                     if (false) {
                                        std::string file_name = "rama-trial-" + terminus_type + "-" + chain_id + ".pdb";
                                        m.write_file(file_name, b_factor);
                                     }
                                     auto density_sum_for_new_residues = get_density_sum_for_new_residues(m);
                                     unsigned int n_atoms_in_new_residues = density_sum_for_new_residues.second;
                                     if (n_atoms_in_new_residues < 4) {
                                        std::cout << "ERROR:: extend_this_chain(): too few atoms in residue " << n_atoms_in_new_residues << std::endl;
                                     } else {
                                        float av = density_sum_for_new_residues.first/static_cast<float>(density_sum_for_new_residues.second);
                                        std::cout << "debug:: in extend_this_chain " << chain_id << " terminus type " << terminus_type
                                                  << " residue_p with index : " << residue_index << " " << coot::residue_spec_t(residue_p)
                                                  << " densities " << average_density_per_atom_for_molecule
                                                  << " crit: " << 0.6 * average_density_per_atom_for_molecule << " this_residue_av: " << av << std::endl;
                                        if (av > 0.6 * average_density_per_atom_for_molecule) {
                                           add_CB_to_residue_maybe(&m, xmap, average_density_per_atom_for_molecule);
                                           int n_residues_pre = chain_p->GetNumberOfResidues();
                                           add_residues_to_chain(m, chain_p, residue_seqnum, terminus_type);
                                           mol->FinishStructEdit(); // because mol is not passed to add_residue_to_chain. Perhaps it should be.
                                           coot::util::pdbcleanup_serial_residue_numbers(mol);
                                           int n_residues_post = chain_p->GetNumberOfResidues();
                                           if (n_residues_post == n_residues_pre) {
                                              std::cout << "ERROR:: extend_this_chain(): residue was not added to chain " << chain_id << std::endl;
                                              status = false;
                                           } else {
                                              // happy path
                                              // refine 3, keep 2: i.e. adding 2 new residues, keeping 1.
                                              refine_triple(residue_seqnum, chain_p, mol, terminus_type, xmap, geom, thread_pool_p, n_threads, weight);
                                              delete_extra_residue(residue_seqnum, chain_p, mol, terminus_type);
                                              n_residues_post = chain_p->GetNumberOfResidues();
                                              if (n_residues_post == n_residues_pre)
                                                 status = false;
                                              if (new_residue_is_deformed(residue_p, terminus_type, mol)) {
                                                 std::cout << "INFO:: extend_this_chain(): stop trace because new residue is deformed" << std::endl;
                                                 delete_extra_residue(residue_seqnum, chain_p, mol, terminus_type);
                                                 status = false;
                                              } else {
                                                 if (new_residue_crashing_into_other_chain(residue_p, terminus_type, mol)) {
                                                    std::cout << "INFO:: extend_this_chain(): stop trace because new residue overlaps other chain" << std::endl;
                                                    delete_extra_residue(residue_seqnum, chain_p, mol, terminus_type);
                                                    status = false;
                                                 } else {
                                                    status = true; // keep going!
                                                 }
                                              }
                                           }
                                        }
                                     }
                                  }
                               }
                               return status;
                            };

   auto add_udd_atom_index_to_molecule = [] (mmdb::Manager *mol) {

                                            int imod = 1;
                                            mmdb::Model *model_p = mol->GetModel(imod);
                                            if (model_p) {
                                               int udd_atom_index_handle = mol->RegisterUDInteger(mmdb::UDR_ATOM, "atom index");

                                               if (false) {
                                                  int udd_atom_index_handle_test = mol->GetUDDHandle(mmdb::UDR_ATOM, "atom index");
                                                  std::cout << "debug:: add_udd_atom_index_to_molecule(): " << mol << " udd_atom_index_handle_test "
                                                            << udd_atom_index_handle_test << std::endl;
                                               }

                                               int n_chains = model_p->GetNumberOfChains();
                                               int atom_index = 0;
                                               for (int ichain=0; ichain<n_chains; ichain++) {
                                                  mmdb::Chain *chain_p = model_p->GetChain(ichain);
                                                  int n_res = chain_p->GetNumberOfResidues();
                                                  for (int ires=0; ires<n_res; ires++) {
                                                     mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                                                     if (residue_p) {
                                                        int n_atoms = residue_p->GetNumberOfAtoms();
                                                        for (int iat=0; iat<n_atoms; iat++) {
                                                           mmdb::Atom *at = residue_p->GetAtom(iat);
                                                           int ierr = at->PutUDData(udd_atom_index_handle, atom_index);
                                                           if (ierr != mmdb::UDDATA_Ok) {
                                                              std::cout << "ERROR:: giving udd data atom index to atom residue " << coot::atom_spec_t(at) << std::endl;
                                                           }
                                                           ierr = at->GetUDData(udd_atom_index_handle, atom_index);
                                                           if (false) { // debugging
                                                              if (ierr == mmdb::UDDATA_Ok) {
                                                                 std::cout << "debug:: add_udd_atom_index_to_molecule() atom " << at << " " << coot::atom_spec_t(at)
                                                                           << " had index " << atom_index << std::endl;
                                                              } else {
                                                                 std::cout << "debug:: add_udd_atom_index_to_molecule() atom " << at << " failed to get uddata" << std::endl;
                                                              }
                                                           }
                                                           atom_index++;
                                                        }
                                                     }
                                                  }
                                               }
                                            }
                                         };


   // This its the 3-residue version of extend_this_chain
   //
   // return a bool - for success status (residue added).
   auto build_3_keep_1 = [get_density_score_for_new_residues_frag, add_residues_to_chain, add_udd_atom_index_to_molecule,
                          refine_isolated_quad, delete_extra_residue, delete_extra_residues, new_residue_is_deformed,
                          thread_pool_p, n_threads, xmap] (mmdb::Chain *chain_p, mmdb::Manager *mol,
                                                           const std::string &terminus_type,
                                                           float average_density_per_atom_for_molecule,
                                                           float weight, unsigned int n_phi_psi_trials,
                                                           const coot::protein_geometry &geom) {

                            float crit_sf = 0.345; // the critcal ratio between the density for the new fragment and the average density for the model so far
                            bool status = false; // initially no residue added.
                            float density_level_crit_for_main_chain = 0.5 * average_density_per_atom_for_molecule;
                            // std::cout << "debug density_level_crit_for_main_chain: " << density_level_crit_for_main_chain << std::endl;
                            int n_res = chain_p->GetNumberOfResidues();
                            if (n_res > 0) {
                               std::string chain_id(chain_p->GetChainID());
                               int residue_index = 0;
                               int n_residues_in_chain =  chain_p->GetNumberOfResidues();
                               if (terminus_type == "C") {
                                  residue_index = n_residues_in_chain -1;
                               }
                               mmdb::Residue *residue_p  = chain_p->GetResidue(residue_index);
                               mmdb::Residue *res_prev_p = chain_p->GetResidue(residue_index-1);
                               if (false)
                                  std::cout << "DEBUG:: in build_3_keep_1(): for chain-id " << chain_id << " terminus type " << terminus_type
                                            << " current n_residues_in_chain " << n_residues_in_chain <<  " residue_p with index : "
                                            << residue_index << " " << coot::residue_spec_t(residue_p) << std::endl;
                               if (!residue_p) {

                                  std::cout << "DEBUG:: in build_3_keep_1(): NULL residue_p" << std::endl;

                               } else {
                                  // happy path
                                  int residue_seqnum = residue_p->GetSeqNum();

                                  coot::new_residue_by_3_phi_psi nr3phipsi(terminus_type, residue_p, chain_id);
                                  if (terminus_type == "C")
                                     nr3phipsi.set_upstream_neighbour(res_prev_p);
                                  if (terminus_type == "N") {
                                     mmdb::Residue *res_next_p = chain_p->GetResidue(1);
                                     nr3phipsi.set_downstream_neighbour(res_next_p);
                                  }
                                  nr3phipsi.add_thread_pool(thread_pool_p, n_threads);
                                  // auto tp_0 = std::chrono::high_resolution_clock::now();
                                  float min_density_level_for_connecting_atom = 0.5 * average_density_per_atom_for_molecule;
                                  coot::minimol::fragment frag = nr3phipsi.best_fit_phi_psi(n_phi_psi_trials, xmap, min_density_level_for_connecting_atom);
                                  // auto tp_1 = std::chrono::high_resolution_clock::now();
                                  //auto d10 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_1 - tp_0).count();
                                  // std::cout << "Timings: from build_3_keep_1(): best_fit_phi_psi() " << d10 << " milliseconds" << std::endl;

                                  if (false) {
                                     std::string file_name = "rama-trial-b3k1-" + terminus_type + "-" + chain_id + ".pdb";
                                     frag.write_file(file_name);
                                  }

                                  auto density_score_for_new_residues = get_density_score_for_new_residues_frag(frag);
                                  unsigned int n_atoms_in_new_residues = density_score_for_new_residues.second;
                                  if (n_atoms_in_new_residues < 4) {
                                     std::cout << "ERROR:: build_3_keep_1(): too few atoms in residue: " << n_atoms_in_new_residues << std::endl;
                                  } else {
                                     float av = density_score_for_new_residues.first;
                                     coot::minimol::molecule m(frag); // because I don't want to rewrite add_CB_to_residue_maybe() and add_residues_to_chain()

                                     if (true) {
                                        std::cout << "DEBUG:: in build_3_keep_1(): chain " << chain_id << " terminus type " << terminus_type
                                                  << " residue_p with index : " << residue_index << " " << coot::residue_spec_t(residue_p)
                                                  << " densities " << average_density_per_atom_for_molecule
                                                  << " crit: " << crit_sf * average_density_per_atom_for_molecule << " this_residue_av: " << av << " ";
                                        if (av < crit_sf * average_density_per_atom_for_molecule)
                                           std::cout << "Low density" << std::endl;
                                        else
                                           std::cout << std::endl;
                                     }

                                     bool density_is_good_enough_to_continue = true;
                                     if (av < crit_sf * average_density_per_atom_for_molecule)
                                        density_is_good_enough_to_continue = false;

                                     if (terminus_type == "C") {
                                        // N of next residue must be in density at least 1rmsd!
                                        std::pair<bool, coot::minimol::atom> N_pair = frag[residue_seqnum+1].get_atom(" N  ");
                                        if (N_pair.first) {
                                           const clipper::Coord_orth &pos = N_pair.second.pos;
                                           float f = coot::util::density_at_point(xmap, pos);
                                           if (f < density_level_crit_for_main_chain) {
                                              density_is_good_enough_to_continue = false;
                                              std::cout << "Terminate trace with Bad N density for residue " << residue_seqnum+1 << std::endl;
                                           }
                                        } else {
                                           density_is_good_enough_to_continue = false; // couldn't find the atom
                                        }
                                     }
                                     if (terminus_type == "N") {
                                        // C of previous residue must be in density at least 1rmsd!
                                        std::pair<bool, coot::minimol::atom> C_pair = frag[residue_seqnum-1].get_atom(" C  ");
                                        if (C_pair.first) {
                                           const clipper::Coord_orth &pos = C_pair.second.pos;
                                           float f = coot::util::density_at_point(xmap, pos);
                                           if (f < density_level_crit_for_main_chain) {
                                              density_is_good_enough_to_continue = false;
                                              std::cout << "Terminate trace with Bad C density for residue " << residue_seqnum-1 << std::endl;
                                           }
                                        } else {
                                           density_is_good_enough_to_continue = false; // couldn't find the atom
                                        }
                                     }

                                     if (! density_is_good_enough_to_continue) {
                                        // Sadge.
                                     } else {
                                        // add_CB_to_residue_maybe(&m, xmap, average_density_per_atom_for_molecule); // already has one now
                                        int n_residues_pre = chain_p->GetNumberOfResidues();
                                        add_residues_to_chain(m, chain_p, residue_seqnum, terminus_type, 3);
                                        mol->FinishStructEdit(); // because mol is not passed to add_residue_to_chain. Perhaps it should be.
                                        coot::util::pdbcleanup_serial_residue_numbers(mol);
                                        add_udd_atom_index_to_molecule(mol); // so that coot::util::create_mmdbmanager_from_residue_vector(const std::vector<mmdb::Residue *> &res_vec,
                                                                             // doesn't complain about it missing.
                                        int n_residues_post = chain_p->GetNumberOfResidues();
                                        if (n_residues_post == n_residues_pre) {
                                           std::cout << "ERROR:: in build_3_keep_1(): residue was not added to chain " << chain_id << std::endl;
                                           status = false;
                                        } else {
                                           // happy path
                                           // refine 3, keep 2: i.e. adding 2 new residues, keeping 1.
                                           unsigned int n_threads_for_quad_refine = 10;
                                           refine_isolated_quad(residue_seqnum, chain_p, mol, terminus_type, xmap, geom, thread_pool_p, n_threads_for_quad_refine, weight);
                                           delete_extra_residues(residue_seqnum, chain_p, mol, terminus_type);
                                           n_residues_post = chain_p->GetNumberOfResidues();
                                           if (n_residues_post == n_residues_pre)
                                              status = false;
                                           if (new_residue_is_deformed(residue_p, terminus_type, mol)) {
                                              std::cout << "INFO:: in build_3_keep_1(): stop trace because new residue is deformed "
                                                        << coot::residue_spec_t(residue_p) << std::endl;
                                              delete_extra_residue(residue_seqnum, chain_p, mol, terminus_type);
                                              status = false;
                                           } else {
                                              status = true; // keep going!
                                              if (false)
                                                 std::cout << "DEBUG:: in build_3_keep_1(): chain " << chain_id << " terminus type " << terminus_type
                                                           << " residue_p with index : " << residue_index << " " << coot::residue_spec_t(residue_p)
                                                           << " status: keep going! " << " " << status << std::endl;
                                           }
                                        }
                                     }
                                  }
                               }
                            }

                            // status = false; // stop after the first one
                            return status;
                         };
   


   // ------------------------------------------------------------------------------------------------------
   // ------------------------------------------------------------------------------------------------------
   // rama_rsr_extend_fragments()
   // ------------------------------------------------------------------------------------------------------
   // ------------------------------------------------------------------------------------------------------

   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);

   auto extend_chain_func = [build_3_keep_1] (mmdb::Chain *chain_p, mmdb::Manager *mol,
                                              float average_density_per_atom,
                                              float weight,
                                              unsigned int n_phi_psi_trials,
                                              const coot::protein_geometry &geom,
                                              unsigned int *update_count) {
                      
                               bool status = false;
                               std::vector<std::string> terminus_types = {"N", "C"};
                               std::pair<mmdb::Chain *, mmdb::Manager *> chain_mol_pair = coot::util::copy_chain(chain_p); // add it into a molecule hierarchy

                               for (auto it = terminus_types.begin(); it != terminus_types.end(); ++it) {
                                  const std::string &terminus_type(*it);
                                  do {
                                     status = build_3_keep_1(chain_mol_pair.first, chain_mol_pair.second, terminus_type,
                                                             average_density_per_atom, weight, n_phi_psi_trials, geom);
                                     if (status) (*update_count)++;
                                  } while (status);
                               }
                               // now copy the atoms (including any new residues and atoms) of chain_mol.first into chain_p.
                               // Using mol to do a FinishStructEdit() here looks dangerous. Perhaps do it at the end of this function (not this lambda)
                               coot::util::replace_chain_contents_with_atoms_from_chain(chain_p, mol, chain_mol_pair.first); // (to_chain, from_chain)
                            };

   if (model_p) {

      float average_density_per_atom_for_molecule = get_average_density_per_atom(mol, xmap);
      int n_chains = model_p->GetNumberOfChains();
#if 1 // threaded.
      unsigned int n_rounds = 50;
      std::vector<std::pair<unsigned int, unsigned int> > cir = coot::atom_index_ranges(n_chains, n_rounds);

      for (unsigned int i=0; i<cir.size(); i++) {
         const auto &chain_index_pair = cir[i];
         std::cout << "New chain batch: (chain index range " << chain_index_pair.first << " " << chain_index_pair.second << ")"
                   << std::endl;
         for (unsigned int ich=chain_index_pair.first; ich<chain_index_pair.second; ich++) {
            mmdb::Chain *chain_p = model_p->GetChain(ich);
            std::cout << " " << chain_p->GetChainID();
         }
         std::cout << std::endl;

         std::vector<std::thread> threads;
         for (unsigned int ichain=cir[i].first; ichain<cir[i].second; ichain++) {
            mmdb::Chain *chain_p = model_p->GetChain(ichain);
            threads.push_back(std::thread(extend_chain_func, chain_p, mol, average_density_per_atom_for_molecule, weight, n_phi_psi_trials,
                                          std::cref(geom), update_count));
         }
         for (unsigned int ithr=0; ithr<threads.size(); ithr++)
            threads[ithr].join();
      }
#endif
#if 0
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         extend_chain_func(chain_p, mol, average_density_per_atom_for_molecule, weight, n_phi_psi_trials, geom, update_count);
      }
#endif
   }

}
