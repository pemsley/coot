
#include <algorithm>

#include "coot-utils/coot-coord-utils.hh"
#include "merge-atom-selections.hh"
#include "protein_db/protein_db_utils.h"


// use_symmetry default true
// using_missing_loop_fit default true
void
coot::merge_C_and_N_termii(mmdb::Manager *mol, bool use_symmetry, bool using_missing_loop_fit) {

   enum close_type { NONE, C_AND_N, N_AND_C };

   auto residues_are_close = [] (close_type check_this_close_type, mmdb::Residue *r_1, mmdb::Residue *r_2, double dist_crit_for_close_atoms=5.0) {

                                bool status = false;
                                mmdb::Atom *r_1_C = r_1->GetAtom(" C  ");
                                mmdb::Atom *r_1_N = r_1->GetAtom(" N  ");
                                mmdb::Atom *r_2_C = r_2->GetAtom(" C  ");
                                mmdb::Atom *r_2_N = r_2->GetAtom(" N  ");
                                if (r_1_C && r_1_N && r_2_C && r_2_N) {
                                   clipper::Coord_orth p_1 = co(r_1_C);
                                   clipper::Coord_orth p_2 = co(r_2_N);
                                   clipper::Coord_orth p_3 = co(r_1_N);
                                   clipper::Coord_orth p_4 = co(r_2_C);
                                   double dd_CN = (p_2-p_1).lengthsq();
                                   double dd_NC = (p_3-p_4).lengthsq();
                                   if (check_this_close_type == C_AND_N) {
                                      if (dd_CN < dd_NC) { // this has to be true for residues going in the same direction
                                         double d = std::sqrt(dd_CN);
                                         if (d < dist_crit_for_close_atoms)
                                            status = true;
                                      }
                                   }
                                   if (check_this_close_type == N_AND_C) {
                                      if (dd_NC < dd_CN) {
                                         double d = std::sqrt(dd_NC);
                                         if (d < dist_crit_for_close_atoms)
                                            status = true;
                                      }
                                   }
                                }
                                return status;
                             };

   auto have_close_terminii = [residues_are_close] (mmdb::Chain *i_chain_p, mmdb::Chain *j_chain_p) {
                                 bool status = false;
                                 close_type ct = NONE;
                                 int n_residues_i_chain = i_chain_p->GetNumberOfResidues();
                                 int n_residues_j_chain = j_chain_p->GetNumberOfResidues();
                                 if (n_residues_i_chain > 2) {
                                    if (n_residues_j_chain > 2) {
                                       mmdb::Residue *i_chain_N_terminus = i_chain_p->GetResidue(0);
                                       mmdb::Residue *j_chain_N_terminus = j_chain_p->GetResidue(0);
                                       mmdb::Residue *i_chain_C_terminus = i_chain_p->GetResidue(n_residues_i_chain-1);
                                       mmdb::Residue *j_chain_C_terminus = j_chain_p->GetResidue(n_residues_j_chain-1);

                                       double dist_crit_for_close_atoms = 5.0; // for 0 or 1
                                       // dist_crit_for_close_atoms = 8.0 for 2 residue insertion, i.e. more than this is 2 far
                                       // appart to be melded by adding 2 residues in the gap.

                                       if (residues_are_close(N_AND_C, i_chain_N_terminus, j_chain_C_terminus, dist_crit_for_close_atoms)) {
                                          status = true;
                                          ct = N_AND_C;
                                       }
                                       if (residues_are_close(C_AND_N, i_chain_C_terminus, j_chain_N_terminus, dist_crit_for_close_atoms)) {
                                          status = true;
                                          ct = C_AND_N;
                                       }
                                    }
                                 }
                                 return std::make_pair(status, ct);
                              };

   auto mergeable_with_0_residues_insertion = [residues_are_close] (mmdb::Chain *i_chain_p, mmdb::Chain *j_chain_p, close_type ct) {
                                                 bool status = false;
                                                 int n_residues_i_chain = i_chain_p->GetNumberOfResidues();
                                                 int n_residues_j_chain = j_chain_p->GetNumberOfResidues();
                                                 mmdb::Residue *i_chain_N_terminus = i_chain_p->GetResidue(0);
                                                 mmdb::Residue *j_chain_N_terminus = j_chain_p->GetResidue(0);
                                                 mmdb::Residue *i_chain_C_terminus = i_chain_p->GetResidue(n_residues_i_chain-1);
                                                 mmdb::Residue *j_chain_C_terminus = j_chain_p->GetResidue(n_residues_j_chain-1);
                                                 double dist_crit = 4.0;
                                                 if (ct == C_AND_N) {
                                                    status = residues_are_close(ct, i_chain_C_terminus, j_chain_N_terminus, dist_crit);
                                                 }
                                                 if (ct == N_AND_C) {
                                                    status = residues_are_close(ct, i_chain_N_terminus, j_chain_C_terminus, dist_crit); // did I swap correctly?
                                                 }
                                                 return status;
                                              };

   auto mergeable_with_1_residue_insertion = [] (mmdb::Chain *i_chain_p, mmdb::Chain *j_chain_p, close_type ct) {
                                                bool status = false;
                                                return status;
                                             };

   auto mergeable_with_2_residues_insertion = [] (mmdb::Chain *i_chain_p, mmdb::Chain *j_chain_p, close_type ct) {
                                                 bool status = false;
                                                 return status;
                                              };

   // 1 copy all of the residues of the j_chain into the i_chain
   // 2 delete the j_chain
   // 3 finishstructuraledit
   auto merge_chains_N_to_C = [] (mmdb::Chain *chain_with_CO_p, mmdb::Chain *chain_with_N_p, mmdb::Manager *mol, int gap_size) {
                                // 1 copy all of the residues of the j_chain into the i_chain
                                // 2 delete the j_chain
                                // 3 finishstructuraledit
                                int n_residues_i_chain = chain_with_CO_p->GetNumberOfResidues();
                                int n_residues_j_chain = chain_with_N_p->GetNumberOfResidues();
                                mmdb::Residue *i_chain_C_terminus = chain_with_CO_p->GetResidue(n_residues_i_chain-1);
                                int rn_base = i_chain_C_terminus->GetSeqNum() + 1;  // or should I use a rn-delta?
                                for (int ires=0; ires<n_residues_j_chain; ires++) {
                                   mmdb::Residue *residue_p = chain_with_N_p->GetResidue(ires);
                                   if (residue_p) {
                                      mmdb::Residue *new_residue = util::deep_copy_this_residue(residue_p);
                                      new_residue->seqNum = rn_base + ires + gap_size;
                                      chain_with_CO_p->AddResidue(new_residue);
                                   }
                                }
                                mol->FinishStructEdit();
                                util::pdbcleanup_serial_residue_numbers(mol);
                                delete chain_with_N_p;
                                mol->FinishStructEdit();
                                util::pdbcleanup_serial_residue_numbers(mol);
                              };

   
   auto make_fragment_chain = [] (const std::vector<coot::residue_spec_t> &residue_specs, mmdb::Manager *mol) {
                      
                                 ProteinDB::Chain chain;

                                 std::vector<coot::residue_spec_t> local_specs = residue_specs;
                                 std::map<std::string, int> chain_id_map;

                                 std::vector<coot::residue_spec_t>::const_iterator it;
                                 for (it=local_specs.begin(); it!= local_specs.end(); ++it)
                                    chain_id_map[it->chain_id]++;

                                 if (chain_id_map.size() != 1) {
                                    std::cout << "WARNING:: all residues need to be in the same chain. Aborted loop selection"
                                              << std::endl;
                                 } else {

                                    std::sort(local_specs.begin(), local_specs.end());

                                    int n_loop_residues = local_specs.back().res_no - local_specs.begin()->res_no + 1;
      
                                    for (int i_loop_res=0; i_loop_res<n_loop_residues; i_loop_res++) {

                                       coot::residue_spec_t test_spec(*local_specs.begin());
                                       test_spec.res_no += i_loop_res;

                                       if (std::find(local_specs.begin(), local_specs.end(), test_spec) == local_specs.end()) {
                                          // Add a null residue
                                          std::cout << "Added a null for " << test_spec << std::endl;
                                          ProteinDB::Residue residue;
                                          chain.add_residue(residue);
                                       } else { 
                                          mmdb::Residue *residue_p = coot::util::get_residue(test_spec, mol);
                                          if (! residue_p) {
                                             std::cout << "oops - missing residue " << test_spec << std::endl;
                                             std::cout << "Added a null for " << test_spec << std::endl;
                                             ProteinDB::Residue residue;
                                             chain.add_residue(residue);
                                          } else {
                                             std::string type = residue_p->GetResName();
                                             clipper::Coord_orth  n_pos;
                                             clipper::Coord_orth ca_pos;
                                             clipper::Coord_orth  c_pos;

                                             mmdb::Atom *at_n  = residue_p->GetAtom(" N  ");
                                             mmdb::Atom *at_ca = residue_p->GetAtom(" CA ");
                                             mmdb::Atom *at_c  = residue_p->GetAtom(" C  ");
                                             if (at_n)
                                                n_pos = clipper::Coord_orth( at_n->x,  at_n->y,  at_n->z);
                                             if (at_ca)
                                                ca_pos = clipper::Coord_orth(at_ca->x, at_ca->y, at_ca->z);
                                             if (at_c)
                                                n_pos = clipper::Coord_orth( at_c->x,  at_c->y,  at_c->z);


                                             if (at_n && at_ca && at_c) { 
                                                ProteinDB::Residue residue(n_pos, ca_pos, c_pos, type);
                                                chain.add_residue(residue);
                                                // std::cout << "Added a real residue for " << test_spec << std::endl;
                                             }
                                          }
                                       }
                                    }
                                 }
                                 return chain;
                              };


   auto merge_chains = [merge_chains_N_to_C] (mmdb::Chain *i_chain_p, mmdb::Chain *j_chain_p, mmdb::Manager *mol, close_type ct,
                                              int number_of_residues_for_insertion) {

                          if (number_of_residues_for_insertion == 0) {
                             if (ct == C_AND_N) merge_chains_N_to_C(i_chain_p, j_chain_p, mol, 0);
                             if (ct == N_AND_C) merge_chains_N_to_C(j_chain_p, i_chain_p, mol, 0);
                          }

                          if (number_of_residues_for_insertion == 1) {
                             int res_no_of_residue_with_CO = -1;
                             mmdb::Chain *chain_p = 0;
                             if (ct == C_AND_N) chain_p = i_chain_p;
                             if (ct == N_AND_C) chain_p = j_chain_p;
                             if (ct == C_AND_N) res_no_of_residue_with_CO = i_chain_p->GetResidue(i_chain_p->GetNumberOfResidues()-1)->GetSeqNum();
                             if (ct == N_AND_C) res_no_of_residue_with_CO = j_chain_p->GetResidue(j_chain_p->GetNumberOfResidues()-1)->GetSeqNum();
                             if (ct == C_AND_N) merge_chains_N_to_C(i_chain_p, j_chain_p, mol, 1);
                             if (ct == N_AND_C) merge_chains_N_to_C(j_chain_p, i_chain_p, mol, 1);
                             if (chain_p) {
                                // db loop things
                             }
                             // now the residues are in same chain and have a gap of 1 residue between them
                          }
                       };

   if (mol) {
      int imodel = 1;
      mmdb::Model *model_p = mol->GetModel(imodel);
      if (model_p) {
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<(n_chains-1); ichain++) {
            mmdb::Chain *i_chain_p = model_p->GetChain(ichain);
            for (int jchain=(ichain+1); jchain<n_chains; jchain++) {
               mmdb::Chain *j_chain_p = model_p->GetChain(jchain);
               auto ct = have_close_terminii(i_chain_p, j_chain_p);
                  std::cout << "DEBUG:: chain " << i_chain_p->GetChainID() << " " << j_chain_p->GetChainID()
                            << " have_close_terminii: " << ct.first << " " << ct.second << std::endl;
               if (ct.first) {
                  // should I also check that the density between the C and the N is sensible?
                  bool status = mergeable_with_0_residues_insertion(i_chain_p, j_chain_p, ct.second);
                  std::cout << "debug:: chain " << i_chain_p->GetChainID() << " " << j_chain_p->GetChainID()
                            << " merge status for 0 residue insertion " << status << std::endl;
                  if (status) {
                     merge_chains(i_chain_p, j_chain_p, mol, ct.second, 0);
                  } else {
                     status = mergeable_with_1_residue_insertion(i_chain_p, j_chain_p, ct.second);
                     if (status) {
                        merge_chains(i_chain_p, j_chain_p, mol, ct.second, 1);
                     } else {
                        status = mergeable_with_2_residues_insertion(i_chain_p, j_chain_p, ct.second);
                        if (status) {
                           merge_chains(i_chain_p, j_chain_p, mol, ct.second, 2);
                        }
                     }
                  }
               }
            }
         }
      }
   }

}
