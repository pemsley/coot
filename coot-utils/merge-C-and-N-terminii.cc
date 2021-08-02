
#include <algorithm>

#define _USE_MATH_DEFINES
#include <math.h> // come on C++, get this sorted

#include "utils/coot-utils.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "merge-C-and-N-terminii.hh"
#include "protein_db/protein_db_utils.h"

void
coot::merge_C_and_N_terminii_0_gap(mmdb::Manager *mol) {

   clipper::Xmap<float> xmap_dummy;
   // ideally, use_symmetry would be true
   merge_C_and_N_terminii(mol, xmap_dummy, false, false);

}


// use_symmetry default true - but what does it do?
// using_missing_loop_fit optional argument default true
void
coot::merge_C_and_N_terminii(mmdb::Manager *mol,
                             const clipper::Xmap<float> &xmap,
                             bool use_symmetry,
                             bool using_missing_loop_fit) {

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
                                                 double dist_crit = 4.0;
                                                 bool status = false;
                                                 int n_residues_i_chain = i_chain_p->GetNumberOfResidues();
                                                 int n_residues_j_chain = j_chain_p->GetNumberOfResidues();
                                                 mmdb::Residue *i_chain_N_terminus = i_chain_p->GetResidue(0);
                                                 mmdb::Residue *j_chain_N_terminus = j_chain_p->GetResidue(0);
                                                 mmdb::Residue *i_chain_C_terminus = i_chain_p->GetResidue(n_residues_i_chain-1);
                                                 mmdb::Residue *j_chain_C_terminus = j_chain_p->GetResidue(n_residues_j_chain-1);

                                                 // Also
                                                 // the resulting CA-CA-CA angles must make sense
                                                 // the resulting CA-prev-CA-next distances must make sense

                                                 if (ct == C_AND_N) {
                                                    bool rac = residues_are_close(ct, i_chain_C_terminus, j_chain_N_terminus, dist_crit);
                                                    if (rac) {
                                                       mmdb::Residue *i_chain_prev_res_p = i_chain_p->GetResidue(n_residues_i_chain-2);
                                                       mmdb::Residue *j_chain_next_res_p = j_chain_p->GetResidue(1);
                                                       if (i_chain_prev_res_p && j_chain_next_res_p) {
                                                          mmdb::Atom *at_1 = i_chain_prev_res_p->GetAtom(" CA "); // PDBv3 FIXME
                                                          mmdb::Atom *at_2 = i_chain_C_terminus->GetAtom(" CA ");
                                                          mmdb::Atom *at_3 = j_chain_N_terminus->GetAtom(" CA ");
                                                          mmdb::Atom *at_4 = j_chain_next_res_p->GetAtom(" CA ");
                                                          if (at_1 && at_2 && at_3 && at_4) {
                                                             clipper::Coord_orth p1(coot::co(at_1));
                                                             clipper::Coord_orth p2(coot::co(at_2));
                                                             clipper::Coord_orth p3(coot::co(at_3));
                                                             clipper::Coord_orth p4(coot::co(at_4));
                                                             double angle_1 = clipper::Coord_orth::angle(p1, p2, p3);
                                                             double angle_2 = clipper::Coord_orth::angle(p2, p3, p4);
                                                             double d1 = std::sqrt((p3-p1).lengthsq());
                                                             double d2 = std::sqrt((p4-p2).lengthsq());
                                                             if (angle_1 > 0.49 * M_PI)
                                                                if (angle_2 > 0.49 * M_PI)
                                                                   if (d1 > 5.1)
                                                                      if (d2 > 5.1)
                                                                         status = true;
                                                          }
                                                       }
                                                    }
                                                 }
                                                 if (ct == N_AND_C) {
                                                    bool rac = residues_are_close(ct, i_chain_N_terminus, j_chain_C_terminus, dist_crit);
                                                    if (rac) {
                                                       mmdb::Residue *i_chain_next_res_p = i_chain_p->GetResidue(1);
                                                       mmdb::Residue *j_chain_prev_res_p = j_chain_p->GetResidue(n_residues_i_chain-2);
                                                       if (i_chain_next_res_p && j_chain_prev_res_p) {
                                                          mmdb::Atom *at_1 = j_chain_prev_res_p->GetAtom(" CA "); // PDBv3 FIXME
                                                          mmdb::Atom *at_2 = j_chain_C_terminus->GetAtom(" CA ");
                                                          mmdb::Atom *at_3 = i_chain_N_terminus->GetAtom(" CA ");
                                                          mmdb::Atom *at_4 = i_chain_next_res_p->GetAtom(" CA ");
                                                          if (at_1 && at_2 && at_3 && at_4) {
                                                             clipper::Coord_orth p1(coot::co(at_1));
                                                             clipper::Coord_orth p2(coot::co(at_2));
                                                             clipper::Coord_orth p3(coot::co(at_3));
                                                             clipper::Coord_orth p4(coot::co(at_4));
                                                             double angle_1 = clipper::Coord_orth::angle(p1, p2, p3);
                                                             double angle_2 = clipper::Coord_orth::angle(p2, p3, p4);
                                                             double d1 = std::sqrt((p3-p1).lengthsq());
                                                             double d2 = std::sqrt((p4-p2).lengthsq());
                                                             if (angle_1 > 0.49 * M_PI)
                                                                if (angle_2 > 0.49 * M_PI)
                                                                   if (d1 > 5.1)
                                                                      if (d2 > 5.1)
                                                                         status = true;
                                                          }
                                                       }
                                                    }
                                                 }
                                                 return status;
                                              };

   auto mergeable_with_1_residue_insertion = [residues_are_close] (mmdb::Chain *i_chain_p, mmdb::Chain *j_chain_p, close_type ct) {
                                                bool status = false;
                                                double dist_crit = 5.0;
                                                int n_residues_i_chain = i_chain_p->GetNumberOfResidues();
                                                int n_residues_j_chain = j_chain_p->GetNumberOfResidues();
                                                mmdb::Residue *i_chain_N_terminus = i_chain_p->GetResidue(0);
                                                mmdb::Residue *j_chain_N_terminus = j_chain_p->GetResidue(0);
                                                mmdb::Residue *i_chain_C_terminus = i_chain_p->GetResidue(n_residues_i_chain-1);
                                                mmdb::Residue *j_chain_C_terminus = j_chain_p->GetResidue(n_residues_j_chain-1);
                                                if (ct == C_AND_N) {
                                                   status = residues_are_close(ct, i_chain_C_terminus, j_chain_N_terminus, dist_crit);
                                                }
                                                if (ct == N_AND_C) {
                                                   status = residues_are_close(ct, i_chain_N_terminus, j_chain_C_terminus, dist_crit); // did I swap correctly?
                                                }
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

                                 ProteinDB::Chain chain; /// return this

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
                                             std::cout << "WARNING:: make_fragment_chain(): oops - missing residue " << test_spec << std::endl;
                                             std::cout << "WARNIGN:: make_fragment_chain(): Added a null for " << test_spec << std::endl;
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

   // std::vector<ProteinDB::Chain>
   auto protein_db_loops = [make_fragment_chain, &xmap] (const std::vector<coot::residue_spec_t> &residue_specs,
                                                         mmdb::Manager *mol, int nfrags) {

                              std::string pkg_data_dir = coot::package_data_dir();
                              std::string dir = coot::util::append_dir_dir(pkg_data_dir, "protein_db");
                              std::string file_name = coot::util::append_dir_file(dir, "protein.db");

                              std::vector<clipper::Coord_orth> clash_coords;

                              ProteinDB::Chain chain = make_fragment_chain(residue_specs, mol);

                              ProteinDB::ProteinDBSearch protein_db_search(file_name);
                              std::vector<ProteinDB::Chain> chains = protein_db_search.search(chain, nfrags, xmap, clash_coords);

                              return chains;
                           };

   auto add_chain_to_molecule = [] (const ProteinDB::Chain &chain, const std::string &chain_id,
                                    int first_res_no, mmdb::Manager *mol) {

                                   std::vector<mmdb::Residue *> needs_cb_and_o;

                                   mmdb::Model *model_p = new mmdb::Model;
                                   mmdb::Chain *chain_p = new mmdb::Chain;
                                   mol->AddModel(model_p);
                                   model_p->AddChain(chain_p);
                                   chain_p->SetChainID(chain_id.c_str());
                                   for (int ires=0; ires<chain.size(); ires++) {
                                      if (chain[ires].flag() == ProteinDB::Residue::NONE) {
                                         // do nothing
                                      } else {
                                         // we have a CA at least, and possibly a N and C too.
                                         mmdb::Residue *residue_p = new mmdb::Residue;
                                         chain_p->AddResidue(residue_p);
                                         residue_p->seqNum = ires+first_res_no;
                                         mmdb::Atom *at_p = new mmdb::Atom;
                                         residue_p->AddAtom(at_p);
                                         residue_p->SetResName("UNK");    // Or ALA.
                                         clipper::Coord_orth ca_pos = chain[ires].coord_ca();
                                         at_p->SetCoordinates(ca_pos.x(), ca_pos.y(), ca_pos.z(), 1.0, 30.0);
                                         at_p->SetElementName(" C");
                                         at_p->SetAtomName(" CA ");
                                         if (chain[ires].flag() == ProteinDB::Residue::NORMAL) {
                                            clipper::Coord_orth n_pos = chain[ires].coord_n();
                                            clipper::Coord_orth c_pos = chain[ires].coord_c();
                                            mmdb::Atom *at_p_1 = new mmdb::Atom;
                                            mmdb::Atom *at_p_2 = new mmdb::Atom;
                                            residue_p->AddAtom(at_p_1);
                                            residue_p->AddAtom(at_p_2);
                                            at_p_1->SetCoordinates(n_pos.x(), n_pos.y(), n_pos.z(), 1.0, 30.0);
                                            at_p_2->SetCoordinates(c_pos.x(), c_pos.y(), c_pos.z(), 1.0, 30.0);
                                            at_p_1->SetElementName(" N");
                                            at_p_2->SetElementName(" C");
                                            at_p_1->SetAtomName(" N  ");
                                            at_p_2->SetAtomName(" C  ");
                                            needs_cb_and_o.push_back(residue_p);
                                         }
                                      }
                                   }
                                   mol->FinishStructEdit();
                                   return needs_cb_and_o;
                                };

   auto make_mol = [add_chain_to_molecule] (ProteinDB::Chain chain, const std::string &chain_id, int first_res_no) {
                      mmdb::Manager *mol = new mmdb::Manager;
                      std::vector<mmdb::Residue *> residues = add_chain_to_molecule(chain, chain_id, first_res_no, mol);
                      // add_cbs_and_os(residues, mol);
                      return mol;
                   };

   // chain mol is the model with the extra residue
   // rs is the residue spec of the inserted residue
   // chain_p is the chain in mol when needs to have the residue inserted.
   // clumsy variables names I now think.
   auto merge_loop_residue_into_main_molecule = [] (mmdb::Manager *chain_mol, mmdb::Manager *mol, mmdb::Chain *chain_p,
                                                    const coot::residue_spec_t &rs) {
                                                   mmdb::Residue *residue_p = util::get_residue(rs, chain_mol);
                                                   if (residue_p) {
                                                      mmdb::Residue *residue_copy_p = util::deep_copy_this_residue(residue_p);
                                                      if (residue_copy_p) {
                                                         int insertion_point_res_no = rs.res_no -1;
                                                         // now find that
                                                         int n_residues = chain_p->GetNumberOfResidues();
                                                         for (int ires=0; ires<n_residues; ires++) {
                                                            mmdb::Residue *ch_residue_p = chain_p->GetResidue(ires);
                                                            if (ch_residue_p) {
                                                               int ch_res_no = ch_residue_p->GetSeqNum();
                                                               if (ch_res_no == insertion_point_res_no) {
                                                                  // found it
                                                                  int insertion_index = ires;
                                                                  chain_p->InsResidue(residue_copy_p, insertion_index);
                                                                  mol->FinishStructEdit();
                                                                  break;
                                                               }
                                                            }
                                                         }
                                                      }
                                                   }
                                                };

   auto merge_chains = [merge_chains_N_to_C, protein_db_loops, make_mol, merge_loop_residue_into_main_molecule]
      (mmdb::Chain *i_chain_p, mmdb::Chain *j_chain_p, mmdb::Manager *mol, close_type ct,
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
                                std::string chain_id(chain_p->GetChainID());
                                std::vector<coot::residue_spec_t> residue_specs;
                                int offset_from_merge_point = 3;
                                int merge_point_res_no = res_no_of_residue_with_CO + 1;
                                int first_res_no = res_no_of_residue_with_CO - offset_from_merge_point + 1;
                                int last_res_no  = res_no_of_residue_with_CO + 1 + offset_from_merge_point;
                                for (int ires=first_res_no; ires<=last_res_no; ires++) {
                                   mmdb::Residue *residue_p = chain_p->GetResidue(ires, ""); // use SeqNum API
                                   if (residue_p)
                                      residue_specs.push_back(residue_spec_t(residue_p));
                                }
                                int n_frags = 1;
                                std::vector<ProteinDB::Chain> chains = protein_db_loops(residue_specs, mol, n_frags);
                                if (!chains.empty()) {
                                   for(unsigned int ich=0; ich<chains.size(); ich++) {
                                      mmdb::Manager *chain_mol = make_mol(chains[ich], chain_id, first_res_no);
                                      residue_spec_t rs(chain_id, merge_point_res_no, "");
                                      merge_loop_residue_into_main_molecule(chain_mol, mol, chain_p, rs);
                                   }
                                }
                             }
                             // now the residues are in same chain and have a gap of 1 residue between them
                          }
                       };

   if (mol) {
      int imodel = 1;
      mmdb::Model *model_p = mol->GetModel(imodel);

      if (model_p) {
         bool keep_looping = false; // unless something was merged
         do {
            keep_looping = false;

            int n_chains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<(n_chains-1); ichain++) {
               mmdb::Chain *i_chain_p = model_p->GetChain(ichain);
               if (i_chain_p) {
                  for (int jchain=(ichain+1); jchain<n_chains; jchain++) {
                     mmdb::Chain *j_chain_p = model_p->GetChain(jchain);
                     if (j_chain_p) {
                        auto ct = have_close_terminii(i_chain_p, j_chain_p);
                        if (false)
                           std::cout << "DEBUG:: chain " << i_chain_p->GetChainID() << " " << j_chain_p->GetChainID()
                                     << " have_close_terminii: " << ct.first << " " << ct.second << std::endl;
                        if (ct.first) {
                           // should I also check that the density between the C and the N is sensible?
                           bool status = mergeable_with_0_residues_insertion(i_chain_p, j_chain_p, ct.second);
                           if (true)
                              std::cout << "DEBUG:: merge_C_and_N_terminii(): chain " << i_chain_p->GetChainID() << " "
                                        << j_chain_p->GetChainID() << " merge status for 0 residue insertion "
                                        << status << std::endl;
                           if (status) {
                              merge_chains(i_chain_p, j_chain_p, mol, ct.second, 0);
                           } else {
                              if (using_missing_loop_fit) {
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
                           keep_looping = status; // because something changed.
                        }
                     }
                  }
               }
            }
         } while (keep_looping);
      }
   }
}
