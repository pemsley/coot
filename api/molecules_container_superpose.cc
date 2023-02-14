
#include <iostream>
#include <iomanip>
#include "molecules_container.hh"

// 20230214-PE -------------- all of these function copied from graphics-info-superpose.cc.
//                            I think that they should be in the coot-utils library.

//! superposition (using SSM)
void
molecules_container_t::SSM_superpose(int imol_ref, const std::string &chain_id_ref,
                                     int imol_mov, const std::string &chain_id_mov) {

   int imol_new = -1;
   if (is_valid_model_molecule(imol_ref)) {
      if (is_valid_model_molecule(imol_mov)) {
         atom_selection_container_t asc_ref = molecules[imol_ref].atom_sel;
         atom_selection_container_t asc_mov = molecules[imol_mov].atom_sel;
         // now replace the atom selections
         int save_sel_handle_ref = asc_ref.SelectionHandle;
         int save_sel_handle_mov = asc_mov.SelectionHandle;
         asc_ref.SelectionHandle = asc_ref.mol->NewSelection();
         asc_mov.SelectionHandle = asc_mov.mol->NewSelection();

         asc_ref.mol->SelectAtoms(asc_ref.SelectionHandle, 0, chain_id_ref.c_str(), mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");
         asc_mov.mol->SelectAtoms(asc_mov.SelectionHandle, 0, chain_id_mov.c_str(), mmdb::ANY_RES, "*", mmdb::ANY_RES, "*", "*", "*", "*", "*");

         std::string name = "Chain";
         std::string ref_name = "Chain";
         // superpose in place
         imol_new = superpose_with_atom_selection(asc_ref, asc_mov, imol_mov, name, ref_name, false);

         asc_ref.delete_atom_selection();
         asc_mov.delete_atom_selection();

         asc_ref.SelectionHandle = save_sel_handle_ref;
         asc_mov.SelectionHandle = save_sel_handle_mov;

      }
   }
   // because we replace the position of imol_mov the return value is not interesting.
}


int
molecules_container_t::superpose_with_atom_selection(atom_selection_container_t asc_ref,
                                                     atom_selection_container_t asc_mov,
                                                     int imol_mov,
                                                     std::string moving_mol_name,
                                                     std::string reference_mol_name,
                                                     bool move_copy_of_imol2_flag) {


   int imodel_return = -1;

#ifdef HAVE_SSMLIB

   ssm::PRECISION precision = ssm::PREC_Normal;
   ssm::CONNECTIVITY connectivity = ssm::CONNECT_Flexible;

   // probably not necessary, not sure:
   ssm::SetConnectivityCheck(ssm::CONNECT_Flexible);
   ssm::SetMatchPrecision(ssm::PREC_Normal);

   if (asc_ref.n_selected_atoms > 0) {
      if (asc_mov.n_selected_atoms > 0) {

         mmdb::Manager *mol1 = asc_ref.mol;
         mmdb::Manager *mol2 = asc_mov.mol;

         // debug atom selections:
         mmdb::PAtom *debug_atom_selection1 = NULL;
         mmdb::PAtom *debug_atom_selection2 = NULL;
         int n_selected_atoms_1_debug, n_selected_atoms_2_debug;
         asc_mov.mol->GetSelIndex(asc_mov.SelectionHandle,
                                  debug_atom_selection1,
                                  n_selected_atoms_1_debug);
         asc_ref.mol->GetSelIndex(asc_ref.SelectionHandle,
                                  debug_atom_selection2,
                                  n_selected_atoms_2_debug);

         std::string name1 = reference_mol_name;
         std::string name2 = moving_mol_name;
         std::cout << "superposing..." << std::endl;

         // Remove the pointer one day.
         //
         ssm::Align *SSMAlign = new ssm::Align();
         int rc = SSMAlign->AlignSelectedMatch(mol2, mol1, precision, connectivity,
                                               asc_mov.SelectionHandle,
                                               asc_ref.SelectionHandle);

         if (rc)  {
            std::string ws;
            switch (rc)  {
            case ssm::RC_NoHits :
               std::cout << " *** secondary structure does not match.\n";
               ws = "secondary structure does not match";
               break;
            case ssm::RC_NoSuperposition :
               std::cout << " *** structures are too remote.\n";
               ws = "structures are too remote";
               break;

            // It seems to me that these error codes from SSM are the wrong
            // way round.  This code is adjusted to compensate.

            case ssm::RC_NoGraph :
               std::cout << " *** can't make graph for " << name2 << "\n";
               ws = "can't make graph for " + name2;
               ws += " structure";
               break;
            case ssm::RC_NoVertices :
               std::cout << " *** empty graph for " << name2 << "\n";
               ws = "empty graph for " + name2;
               break;
            case ssm::RC_NoGraph2 :
               std::cout << " *** can't make graph for " << name1 << "\n";
               ws = "can't make graph for " + name1;
               break;
            case ssm::RC_NoVertices2 :
               std::cout << " *** empty graph for " << name1 << "\n";
               ws = "empty graph for " + name1;
               break;
            default :
               std::cout << " *** undocumented return code: " << rc << "\n";
            }
         } else  {

            if (move_copy_of_imol2_flag) {
               mol2 = new mmdb::Manager;
               mol2->Copy(asc_mov.mol, mmdb::MMDBFCM_All);
               std::string name = "Copy_of_";
               name += moving_mol_name;
               // int imol2_new = g.create_molecule();
               // g.molecules[imol2_new].install_model(imol2_new, make_asc(mol2), g.Geom_p(), name, 1);
               // imol2 = imol2_new;
               int imol_new = molecules.size();
               atom_selection_container_t asc_new = make_asc(mol2);
               coot::molecule_t mol_new(asc_new, imol_new, name);
               imol_mov = imol_new;
               imodel_return = install_model(mol_new);
            }

            // OK, let's get a consistent naming system:  1 is moving: 2 is reference

            mmdb::PAtom *atom_selection1 = NULL;
            mmdb::PAtom *atom_selection2 = NULL;
            int n_selected_atoms_1, n_selected_atoms_2;
            asc_mov.mol->GetSelIndex(SSMAlign->selHndCa1, atom_selection1, n_selected_atoms_1);
            asc_ref.mol->GetSelIndex(SSMAlign->selHndCa2, atom_selection2, n_selected_atoms_2);

            std::cout << "number of Ca atoms in selections: " << n_selected_atoms_1
                      << " (moving) and " << n_selected_atoms_2 << " (reference)"
                      << std::endl;

            // mmdb::realtype *rv = SSMAlign->dist1;
            std::string res_type1;
            std::string res_type2;
            std::string res_no_str_1;
            std::string res_no_str_2;
            std::string reference_seq, moving_seq;

            // another go at the sequence alignment used in the superposition.
            //
            // correct: number of residues in reference structure: " << SSMAlign->nres2

            print_ssm_sequence_alignment(SSMAlign,
                                         asc_ref, asc_mov,
                                         atom_selection1, atom_selection2,
                                         n_selected_atoms_1, n_selected_atoms_2,
                                         move_copy_of_imol2_flag);

            make_and_print_horizontal_ssm_sequence_alignment(SSMAlign, asc_ref, asc_mov,
                                                             atom_selection1, atom_selection2,
                                                             n_selected_atoms_1, n_selected_atoms_2);

            // map assignment of secondary structure for Ana
            //
            if (false)
               map_secondary_structure_headers(SSMAlign, asc_ref, asc_mov,
                                              atom_selection1, atom_selection2,
                                              n_selected_atoms_1, n_selected_atoms_2);

            molecules[imol_mov].transform_by(SSMAlign->TMatrix);
            molecules[imol_mov].set_show_symmetry(0);

            //
            // mmdb::ivector Ca1, Ca2;
            // rvector dist1;
            // int nCa1, nCa2;
            // mmdb::mat44 tmat;
            // mmdb::realtype rmsdAchieved, seqIdentity, nCombs;
            // int nAligned, nGaps, nMisD;
            // CSuperpose Superpose;  type in class CSSMAlign
            //
            // Superpose->GetSuperposition(Ca1, dist1, nCa1, Ca2,
            // nCa2, tmat, rmsdAchieved, nAligned, nGaps, seqIdentity,
            // nMisD, nCombs);
            //
            std::cout << "INFO: core rmsd achieved: " << SSMAlign->rmsd << " Angstroems\n"
                      << "      number of residues in reference structure: " << SSMAlign->nres2 << "\n"
                      << "      number of residues in moving structure:    " << SSMAlign->nres1 << "\n"
                      << "      number of residues in aligned sections (reference):  " << SSMAlign->nsel2 << "\n"
                      << "      number of residues in aligned sections (moving):     " << SSMAlign->nsel1 << "\n"
                      << "      number of aligned residues:  " << SSMAlign->nalgn << "\n"
                      << "      number of gaps:              " << SSMAlign->ngaps << "\n"
                      << "      number of misdirections:     " << SSMAlign->nmd << "\n"
                      << "      number of SSE combinations:  " << SSMAlign->ncombs << "\n"
                      << "      sequence identity:           " << SSMAlign->seqIdentity*100.0 << "%\n";

         }
         delete SSMAlign;
      } else {
         std::cout << "WARNING:: Molecule moving has no model atoms\n";
      }
   } else {
      std::cout << "WARNING:: Molecule reference has no model atoms\n";
   }

#endif // HAVE_SSMLIB

   return imodel_return;

}

#ifdef HAVE_SSMLIB
void
molecules_container_t::make_and_print_horizontal_ssm_sequence_alignment(ssm::Align *SSMAlign,
                                                                        atom_selection_container_t asc_ref,
                                                                        atom_selection_container_t asc_mov,
                                                                        mmdb::PAtom *atom_selection1, mmdb::PAtom *atom_selection2,
                                                                        int n_selected_atoms_1, int n_selected_atoms_2) const {

   std::pair<std::string, std::string> aligned_sequences =
      get_horizontal_ssm_sequence_alignment(SSMAlign, asc_ref, asc_mov,
                                            atom_selection1, atom_selection2,
                                            n_selected_atoms_1, n_selected_atoms_2);

   print_horizontal_ssm_sequence_alignment(aligned_sequences);
}
#endif // HAVE_SSMLIB

#ifdef HAVE_SSMLIB
void
molecules_container_t::map_secondary_structure_headers(ssm::Align *SSMAlign,
                                                       atom_selection_container_t asc_ref,
                                                       atom_selection_container_t asc_mov,
                                                       mmdb::PAtom *atom_selection1,
                                                       mmdb::PAtom *atom_selection2,
                                                       int n_selected_atoms_1, int n_selected_atoms_2) const {

   bool debug = true;

   mmdb::Model *model_p = asc_ref.mol->GetModel(1);
   if (model_p) {
      int nhelix = model_p->GetNumberOfHelices();
      int nsheet = model_p->GetNumberOfSheets();
      std::cout << "INFO:: From header, there are " << nhelix << " helices and "
                << nsheet << " sheets\n";
      std::cout << "               Sheet info: " << std::endl;
      std::cout << "------------------------------------------------\n";
      for (int is=1; is<=nsheet; is++) {
         mmdb::Sheet *sheet_p = model_p->GetSheet(is);

         int nstrand = sheet_p->nStrands;
         for (int istrand=0; istrand<nstrand; istrand++) {
            mmdb::Strand *strand_p = sheet_p->strand[istrand];
            if (strand_p) {
               std::cout << strand_p->sheetID << " " << strand_p->strandNo << " "
                         << strand_p->initChainID << " " << strand_p->initSeqNum
                         << " " << strand_p->endChainID << " " << strand_p->endSeqNum
                         << std::endl;
            }
         }
      }
      std::cout << "------------------------------------------------\n";

      std::map<std::string, std::vector<coot::saved_strand_info_t> > saved_strands;
      for (int is=1; is<=nsheet; is++) {
         mmdb::Sheet *sheet_p = model_p->GetSheet(is);
         std::cout << "------- sheet " << sheet_p->sheetID << std::endl;
         int nstrand = sheet_p->nStrands;
         for (int istrand=0; istrand<nstrand; istrand++) {
            mmdb::Strand *strand_p = sheet_p->strand[istrand];
            if (strand_p) {
               bool start_was_found = false;
               bool   end_was_found = false;
               coot::residue_spec_t save_matched_start_mov_res;
               coot::residue_spec_t ref_start_res(strand_p->initChainID, strand_p->initSeqNum, "");
               coot::residue_spec_t   ref_end_res(strand_p->endChainID,  strand_p->endSeqNum,  "");

               for (int i1=0; i1<SSMAlign->nsel1; i1++) {
                  int t_index = SSMAlign->Ca1[i1];
                  // was t_index sensible?
                  if (t_index < SSMAlign->nsel2 && t_index >= 0) {
                     int s_index = SSMAlign->Ca2[t_index];
                     if (s_index == i1) {
                        coot::residue_spec_t matched_atom_res_ref(atom_selection1[i1]->GetResidue()); // SSM match, that is
                        coot::residue_spec_t matched_atom_res_mov(atom_selection2[t_index]->GetResidue());
                        // if we find it in mov, save the ref!  Weird...
                        if (ref_start_res == matched_atom_res_mov) {
                           std::cout << "found start " << ref_start_res << " -> " << matched_atom_res_ref 
                                     << " " << strand_p->sheetID << " " << strand_p->strandNo << std::endl;
                           save_matched_start_mov_res = matched_atom_res_ref;
                           // save_matched_start_mov_res.string_user_data = atom_selection2[t_index]->GetResName();
                           start_was_found = true;
                        }
                        if (ref_end_res == matched_atom_res_mov) {
                           end_was_found = true;
                           std::cout << "found end   " << ref_end_res << " -> " << matched_atom_res_ref
                                     << " " << strand_p->sheetID << " " << strand_p->strandNo << std::endl;
                           if (start_was_found) {
                              coot::saved_strand_info_t ss(save_matched_start_mov_res, matched_atom_res_ref, strand_p->strandNo);
                              saved_strands[strand_p->sheetID].push_back(ss);
                           }
                        }
                     }
                  }
               }
               if (! start_was_found)
                  std::cout << "missing start - no match found for residue " << ref_start_res << std::endl;
               if (! end_was_found)
                  std::cout << "missing end   - no match found for residue " << ref_end_res << std::endl;
            }
         }
      }
      std::cout << "... got " << saved_strands.size() << " saved strands" << std::endl;
      std::map<std::string, std::vector<coot::saved_strand_info_t> >::const_iterator it;
      for (it = saved_strands.begin(); it != saved_strands.end(); it++) {
         std::string key = it->first;
         std::cout << " ---- Sheet " << key << std::endl;
         for (unsigned int ii=0; ii<it->second.size(); ii++) {
            std::cout << "   " << it->second[ii].start << " " << it->second[ii].end << " " << it->second[ii].strand_idx
                      << std::endl;
         }
      }
      if (saved_strands.size() > 0) {
         mmdb::Sheet *sheet_p = new mmdb::Sheet;
         for (it = saved_strands.begin(); it != saved_strands.end(); it++) {
            std::string key = it->first;
            for (unsigned int ii=0; ii<it->second.size(); ii++) {
               mmdb::Strand *strand_p = new mmdb::Strand;
               strcpy(strand_p->sheetID, key.c_str());
               // strcpy(strand_p->initResName, it->second[ii].start.string_user_data.c_str());
               strcpy(strand_p->initChainID, it->second[ii].start.chain_id.c_str());
               strcpy(strand_p->initICode,   it->second[ii].start.ins_code.c_str());
               // strcpy(strand_p->endResName,  it->second[ii].end.string_user_data.c_str());
               strcpy(strand_p->endChainID,  it->second[ii].end.chain_id.c_str());
               strcpy(strand_p->endICode,    it->second[ii].end.ins_code.c_str());
            }
         }
      }
   }

}
#endif // HAVE_SSMLIB


#ifdef HAVE_SSMLIB
//
// To make a GUI dialog with the alignment, we need this function to
// generate a string, rather than print to the screen.  Easily converted.
//
void
molecules_container_t::print_horizontal_ssm_sequence_alignment(std::pair<std::string, std::string> aligned_sequences) const {

   bool debug = 0;

   if (debug) {
      std::cout << "DEBUG:: moving :" << aligned_sequences.first  << ":" << std::endl;
      std::cout << "DEBUG:: target :" << aligned_sequences.second << ":" << std::endl;
   }

   int chars_per_line = 70;
   int lf = aligned_sequences.first.length();
   int ls = aligned_sequences.second.length();
   int l = lf;
   if (ls > l)
      l = ls;
   int n_lines = 1 + l/chars_per_line;
   // std::cout << "DEUBG:: n_lines: " << n_lines << " " << lf << " " << ls << std::endl;

   for (int i=0; i<n_lines; i++) {
      int f_start = i*chars_per_line;
      int f_end = chars_per_line;
      if (f_end > lf)
         f_end = lf - f_start;
      // std::cout << "DEBUG:: comparing  first " << f_start << " with " << lf << std::endl;
      if (f_start < lf)
         std::cout << " Moving: " << aligned_sequences.first.substr(f_start, f_end) << std::endl;

      int s_start = i*chars_per_line;
      int s_end = chars_per_line;
      if (s_end > ls)
         s_end = ls - s_start;
      // std::cout << "DEBUG:: comparing second " << s_start << " with " << ls << std::endl;
      if (s_start < ls)
         std::cout << " Target: " << aligned_sequences.second.substr(s_start, s_end) << std::endl;
      std::cout << std::endl; // for neatness
   }
}
#endif // HAVE_SSMLIB


#ifdef HAVE_SSMLIB
std::pair<std::string, std::string>
molecules_container_t::get_horizontal_ssm_sequence_alignment(ssm::Align *SSMAlign,
                                                             atom_selection_container_t asc_ref,
                                                             atom_selection_container_t asc_mov,
                                                             mmdb::PAtom *atom_selection1, mmdb::PAtom *atom_selection2,
                                                             int n_selected_atoms_1, int n_selected_atoms_2) const {

   std::string s;
   std::string t;
   int previous_t_index = -1;
   int previous_s_index = -1;
   bool debug = false;

   // Talk to Eugene about this: turn on debugging and run Alice
   // Dawson test.
   //
   if (debug) {
      std::cout << "DEBUG:: t_indexes: ";
      for (int i1=0; i1<SSMAlign->nsel1; i1++) {
         std::cout << "[" << i1 << " " << SSMAlign->Ca1[i1] << "] " ;
      }
      std::cout << std::endl;
      std::cout << "DEBUG:: s_indexes: ";
      for (int i2=0; i2<SSMAlign->nsel2; i2++) {
         std::cout << "[" << i2 << " " << SSMAlign->Ca2[i2] << "] " ;
      }
      std::cout << std::endl;
   }

   for (int i1=0; i1<SSMAlign->nsel1; i1++) {
      if (debug)
         std::cout << "getting t_indexes from Ca1 index i1 " << i1 << " range " << SSMAlign->nsel1
                   << std::endl;
      int t_index = SSMAlign->Ca1[i1];
      if (debug)
         std::cout << "debug i1: " << i1 << " t_index is " << t_index << " and range is "
                   << SSMAlign->nsel2 << std::endl;
      if (t_index == -1) {
         s += coot::util::three_letter_to_one_letter(atom_selection1[i1]->GetResName());
         t += "-";
      } else {
         // was t_index sensible?
         if (t_index < SSMAlign->nsel2 && t_index >= 0) {
            int s_index = SSMAlign->Ca2[t_index];

            if (s_index == i1) {
               s += coot::util::three_letter_to_one_letter(atom_selection1[i1]->GetResName());
               t += coot::util::three_letter_to_one_letter(atom_selection2[t_index]->GetResName());
            }
         } else {
            t += "^";
            s += "^";
         }

         // for next round
         previous_t_index = t_index;
      }
   }
   std::cout << std::endl;
   return std::pair<std::string, std::string> (s,t);
}
#endif // HAVE_SSMLIB

#ifdef HAVE_SSMLIB
void
molecules_container_t::print_ssm_sequence_alignment(ssm::Align *SSMAlign,
                                                    atom_selection_container_t asc_ref,
                                                    atom_selection_container_t asc_mov,
                                                    mmdb::PAtom *atom_selection1, mmdb::PAtom *atom_selection2,
                                                    int n_selected_atoms_1, int n_selected_atoms_2,
                                                    bool move_copy_of_imol2_flag) {

   // std::cout << "Another Go...\n\n";

   mmdb::Chain *moving_chain_p = 0;
   mmdb::Chain *reference_chain_p = 0;
   std::string mov_chain_id = std::string(atom_selection1[0]->GetChainID());
   std::string ref_chain_id = std::string(atom_selection2[0]->GetChainID());
   std::string slc_1, slc_2; // single letter code.

   int nchains_ref = asc_ref.mol->GetNumberOfChains(1);
   for (int ich=0; ich<nchains_ref; ich++) {
      mmdb::Chain *chain_p = asc_ref.mol->GetChain(1, ich);
      std::string mol_chain_id(chain_p->GetChainID());
      if (mol_chain_id == std::string(ref_chain_id)) {
         reference_chain_p = chain_p;
         break;
      }
   }
   int nchains_mov = asc_mov.mol->GetNumberOfChains(1);
   for (int ich=0; ich<nchains_mov; ich++) {
      mmdb::Chain *chain_p = asc_mov.mol->GetChain(1, ich);
      std::string mol_chain_id(chain_p->GetChainID());
      if (mol_chain_id == std::string(mov_chain_id)) {
         moving_chain_p = chain_p;
         break;
      }
   }

   if (moving_chain_p && reference_chain_p) {

      // print alignment (distance) table (not sequence)
      //
      if (n_selected_atoms_1 > 0) {
          clipper::RTop_orth ssm_matrix = coot::util::matrix_convert(SSMAlign->TMatrix);
         std::cout << "     Moving      Reference   Distance(/A)" << std::endl;
         for (int ires=0; ires<n_selected_atoms_1; ires++) {
            if (ires < SSMAlign->nalgn) {
               mmdb::Atom *mov_at = atom_selection1[ires];
               std::string ins_code_mov(mov_at->GetInsCode());

               int mov_index = SSMAlign->Ca1[ires];
               std::cout << "      " << mov_at->GetChainID() << " "
                         << std::setw(3) << mov_at->GetSeqNum() << ins_code_mov;
               if ((mov_index > -1) && (mov_index < n_selected_atoms_1)) {
                  mmdb::Atom *ref_at = atom_selection2[mov_index];
                  if (ref_at) {
                     clipper::Coord_orth pos1 = coot::co(mov_at);
                     clipper::Coord_orth pos2 = coot::co(ref_at);
                     clipper::Coord_orth pos3 = pos1.transform(ssm_matrix);
                     double d = clipper::Coord_orth::length(pos3, pos2);
                     std::string ins_code_ref(ref_at->GetInsCode());
                     std::cout << "  <--->  " << ref_at->GetChainID() << " "
                               << std::setw(3) << ref_at->GetSeqNum() << ins_code_ref << "  :  "
                               << std::right << std::setprecision(4) << std::fixed
                               << d << "\n";
                  }
               } else {
                  std::cout << "\n";
               }
            }
         }
      }
   } else {
      std::cout << "ERROR:: Failed to get moving or reference_chain pointer\n";
   }
}
#endif
