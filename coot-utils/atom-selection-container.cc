/* coords/atom-selection-container.hh
 * -*-c++-*-  
 * 
 * Copyright 2005 by The University of York
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include <sys/types.h>  // stating
#include <sys/stat.h>

#include <string.h>
#include "utils/coot-utils.hh"
#include "atom-selection-container.hh"
#include "read-sm-cif.hh"
#include "coot-shelx.hh"
#include "geometry/residue-and-atom-specs.hh"
#include "compat/coot-sysdep.h"
#include "lidia-core/lig-build.hh"
#include "lidia-core/lbg-molfile.hh"
#include "lidia-core-functions.hh"
#include "analysis/stats.hh"

#ifdef USE_GEMMI
#include "gemmi/mmread.hpp"
#include "gemmi/mmdb.hpp"
#endif

#include "utils/logging.hh"
extern logging logger;

mmdb::Residue *
atom_selection_container_t::get_next(mmdb::Residue *residue_in) const {

   mmdb::Residue *r = NULL;

   mmdb::Chain *chain = residue_in->GetChain();
   int this_res_no = residue_in->GetSeqNum();
   int res_no_next = this_res_no + 1;
   for (int i=0; i<n_selected_atoms; i++) {
      if (atom_selection[i]->GetChain() == chain) {
         // for rigor we should do some testing for insertion codes here abouts
         // std::cout << "get_next(): comparing " << atom_selection[i]->GetSeqNum() << " "
         // << res_no_next << std::endl;
         if (atom_selection[i]->GetSeqNum() == res_no_next) {
            r = atom_selection[i]->GetResidue();
            break;
         }
      }
   }
   return r;
}

mmdb::Residue *
atom_selection_container_t::get_previous(mmdb::Residue *residue_in) const {

   mmdb::Residue *r = NULL;

   mmdb::Chain *chain = residue_in->GetChain();
   int this_res_no = residue_in->GetSeqNum();
   int res_no_prev = this_res_no - 1;
   for (int i=0; i<n_selected_atoms; i++) {
      if (atom_selection[i]->GetChain() == chain) {
         // for rigor we should do some testing for insertion codes here abouts
         if (atom_selection[i]->GetSeqNum() == res_no_prev) {
            r = atom_selection[i]->GetResidue();
            break;
         }
      }
   }
   return r;
}

//! clear the atom selection of all pointers
void
atom_selection_container_t::clear_up() {

   if (read_success)
      if (SelectionHandle)
         if (mol)
            mol->DeleteSelection(SelectionHandle);
   delete mol;
   atom_selection = 0;
   mol = 0;
   read_success = 0;
}



// This is used for pick_test  (a function that returns
// an atom selection from a pdb_file name string is not
// generally so useful).
//
atom_selection_container_t
get_atom_selection(std::string pdb_name,
                   bool use_gemmi,
                   bool allow_duplseqnum,
                   bool verbose_mode) {

   auto atom_name_fix_ups = [] (mmdb::Manager *mol) {

      // c.f. fix_wrapped_names() - should these functions be combined?

      // Currently fixes "HH  " in "TYR" // mmdb2 parse mmcif with H atoms

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
                     if (strncmp(residue_p->GetResName(), "TYR", 3) == 0) {
                        int n_atoms = residue_p->GetNumberOfAtoms();
                        for (int iat=0; iat<n_atoms; iat++) {
                           mmdb::Atom *at = residue_p->GetAtom(iat);
                           if (strncmp(at->GetAtomName(), "HH  ", 4) == 0) {
                              at->SetAtomName(" HH ");
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   };

   auto debug_atom_names = [] (mmdb::Manager *mol) {

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
                           std::cout << "       " << coot::atom_spec_t(at) << std::endl;
                        }
                     }
                  }
               }
            }
         }
      }
   };

   auto add_links_to_models = [] (mmdb::Manager *mol, const std::vector<mmdb::Link> &mmdb_links) {

      for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
         mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
            for(const auto &ml : mmdb_links) {
               mmdb::Link *l = new mmdb::Link(ml);
               model_p->AddLink(l);
            }
         }
      }
      if (! mmdb_links.empty())
         mol->FinishStructEdit();
   };

#ifdef USE_GEMMI
   auto transfer_links = [add_links_to_models] (const gemmi::Structure &st, mmdb::Manager *mol) {

      std::vector<mmdb::Link> mmdb_links;
      for (const gemmi::Connection &con : st.connections) {
         if (false) // debugging
            // std::cout << "INFO:: gemmi connection: " << con.name
            //           << " " << con.partner1.str() << " " << con.partner2.str() << std::endl;
            logger.log(log_t::INFO, "gemmi connection:", con.name, con.partner1.str(), con.partner2.str());
         mmdb::Link l;
         std::string atom_name = con.partner1.atom_name;
         bool atom_1_is_metal = false;
         if (atom_name == "ZN") atom_1_is_metal = true;
         if (atom_name == "MG") atom_1_is_metal = true;
         if (atom_name == "FE") atom_1_is_metal = true;
         if (atom_name == "CO") atom_1_is_metal = true;
         if (atom_name == "NI") atom_1_is_metal = true;
         int ll = atom_name.length();
         std::string new_atom_name;
         if (ll == 1) new_atom_name = std::string(" ") + atom_name + std::string("  ");
         if (ll == 2) {
            if (atom_1_is_metal) {
               new_atom_name = atom_name + std::string("  "); // form for metal
            } else {
               new_atom_name = std::string(" ") + atom_name + std::string(" ");
            }
         }
         if (ll == 3) new_atom_name = atom_name + std::string(" ");
         strcpy(l.atName1, new_atom_name.c_str());
         l.aloc1[0] = con.partner1.altloc;
         l.aloc1[1] = 0;
         strcpy(l.resName1, con.partner1.res_id.name.c_str());
         strcpy(l.chainID1, con.partner1.chain_name.c_str());
         l.insCode1[0] = con.partner1.res_id.seqid.icode;
         if (con.partner1.res_id.seqid.icode == ' ') l.insCode1[0] = 0;
         l.insCode1[1] = 0;
         atom_name = con.partner2.atom_name;
         bool atom_2_is_metal = false;
         if (atom_name == "ZN") atom_2_is_metal = true;
         if (atom_name == "MG") atom_2_is_metal = true;
         if (atom_name == "FE") atom_2_is_metal = true;
         if (atom_name == "CO") atom_2_is_metal = true;
         if (atom_name == "NI") atom_2_is_metal = true;
         ll = atom_name.length();
         new_atom_name = atom_name;
         if (ll == 1) new_atom_name = std::string(" ") + atom_name + std::string("  ");
         if (ll == 2) {
            if (atom_2_is_metal) {
               new_atom_name = atom_name + std::string("  "); // form for metal
            } else {
               new_atom_name = std::string(" ") + atom_name + std::string(" ");
            }
         }
         if (ll == 3) new_atom_name = atom_name + std::string(" ");
         strcpy(l.atName2, new_atom_name.c_str());
         l.aloc2[0] = con.partner2.altloc;
         l.aloc2[1] = 0;
         strcpy(l.resName2, con.partner2.res_id.name.c_str());
         strcpy(l.chainID2, con.partner2.chain_name.c_str());
         l.insCode2[0] = con.partner2.res_id.seqid.icode;
         if (con.partner2.res_id.seqid.icode == ' ') l.insCode2[0] = 0;
         l.insCode2[1] = 0;
         if (con.partner1.res_id.seqid.num.has_value()) {
            if (con.partner2.res_id.seqid.num.has_value()) {
               // std::cout << "   debug:: on pushing back link: at-Name-1 :" << l.atName1 << ":" << std::endl;
               // std::cout << "   debug:: on pushing back link: at-Name-2 :" << l.atName2 << ":" << std::endl;
               l.seqNum1 = con.partner1.res_id.seqid.num.value;
               l.seqNum2 = con.partner2.res_id.seqid.num.value;
               if (false)
                  std::cout << "Here's the fresh link: "
                            << l.chainID1 << " seqNum: " << l.seqNum1 << " ins-code: \"" << l.insCode1 << "\" \"" << l.atName1 << "\" to "
                            << l.chainID2 << " seqNum: " << l.seqNum2 << " ins-code: \"" << l.insCode2 << "\" \"" << l.atName2 << "\"" << std::endl;
               mmdb_links.push_back(l);
            }
         }
      }
      add_links_to_models(mol, mmdb_links);
   };
#endif // USE_GEMMI

#ifdef USE_GEMMI
   auto file_name_to_manager_via_gemmi = [transfer_links] (const std::string &pdb_name) {
      mmdb::Manager *mol = nullptr;
      try {
         gemmi::Structure st = gemmi::read_structure_file(pdb_name);
         if (! st.models.empty()) {
            st.merge_chain_parts();
            mol = new mmdb::Manager;
            gemmi::copy_to_mmdb(st, mol);
            transfer_links(st, mol);

            // now fix (cootify) the atom names
            //
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
                              std::string atom_name = at->GetAtomName();
                              int l = atom_name.length();
                              std::string new_atom_name;
                              if (l == 1) new_atom_name = std::string(" ") + atom_name + std::string("  ");
                              if (l == 2) new_atom_name = atom_name + std::string("  ");
                              if (l == 3) new_atom_name = atom_name + std::string(" ");
                              if (l == 1 || l == 2 || l == 3) {
                                 at->SetAtomName(new_atom_name.c_str());
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }
      catch (const std::runtime_error &e) {
         std::cout << "WARNING::" << e.what() << std::endl;
      }
      return mol;
   };
#endif // USE_GEMMI

   if (false) // too noisy
      std::cout << "DEBUG:: get_atom_selection() with file \"" << pdb_name << "\""
                << " use_gemmi: " << use_gemmi << std::endl;

   logger.log(log_t::DEBUG, logging::function_name_t(__FUNCTION__),
              {logging::ltw("with file"), pdb_name, logging::ltw("use_gemmi"), use_gemmi});

   mmdb::ERROR_CODE err;
   mmdb::Manager* MMDBManager;

   // Needed for the error message printing:
   // MMDBManager->GetInputBuffer(S, lcount);
   // Used by reference and as a pointer.  Grimness indeed.
   int  error_count;
   char error_buf[500];

   //   Make routine initializations
   //
   mmdb::InitMatType();

   atom_selection_container_t asc;

   std::string extension = coot::util::file_name_extension(pdb_name);

   // returns e.g. ".ins"

   if (coot::util::extension_is_for_mdl_mol_or_mol2_coords(extension)) {

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
       asc = coot::mol_to_asc_rdkit(pdb_name); // (not a PDB file of course)
       // OK, if that failed, maybe it was an MDL mol file format.
       // Use my parser for that for now.
       if (! asc.read_success) {
          lig_build::molfile_molecule_t m;
          m.read(pdb_name);
          asc = coot::mdl_mol_to_asc(m);
       }
#else
       lig_build::molfile_molecule_t m;
       m.read(pdb_name);
       asc = coot::mdl_mol_to_asc(m);
#endif

    } else {

       if (coot::util::extension_is_for_shelx_coords(extension)) {

          coot::ShelxIns s;
          coot::shelx_read_file_info_t srf = s.read_file(pdb_name);
          // atom_selection_container_t.mol is of type Mymmdb::Manager *
          // currently.
          asc = make_asc(srf.mol);
          MMDBManager = asc.mol;
          if (asc.mol)
             asc.read_success = 1;  // a good idea?

       } else {

          if (use_gemmi) {

#ifdef USE_GEMMI
             MMDBManager = file_name_to_manager_via_gemmi(pdb_name); // new mmdb::Manager;
             if (MMDBManager) {
                asc.read_success = 1;
                asc.mol = MMDBManager;
             }
#else
             MMDBManager = nullptr;
             std::cout << "No GEMMI - sad times " << std::endl;
#endif
          } else {

             MMDBManager = new mmdb::Manager;

             // For mmdb version 1.0.8 and beyond:

             if (allow_duplseqnum)
                MMDBManager->SetFlag ( mmdb::MMDBF_IgnoreBlankLines |
                                       mmdb::MMDBF_IgnoreDuplSeqNum |
                                       mmdb::MMDBF_IgnoreNonCoorPDBErrors |
                                       mmdb::MMDBF_IgnoreHash |
                                       mmdb::MMDBF_IgnoreRemarks);
             else
                MMDBManager->SetFlag ( mmdb::MMDBF_IgnoreBlankLines |
                                       mmdb::MMDBF_IgnoreNonCoorPDBErrors |
                                       mmdb::MMDBF_IgnoreHash |
                                       mmdb::MMDBF_IgnoreRemarks);

#if 0 // 20231020-PE debugging atom names
             for(int imod = 1; imod<=MMDBManager->GetNumberOfModels(); imod++) {
                mmdb::Model *model_p = MMDBManager->GetModel(imod);
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
                                  std::cout << "       " << coot::atom_spec_t(at) << std::endl;
                               }
                            }
                         }
                      }
                   }
                }
             }
#endif
             MMDBManager->PDBCleanup(mmdb::PDBCLEAN_ELEMENT);

             // std::cout << "INFO:: Reading coordinate file: " << pdb_name.c_str() << "\n";
             logger.log(log_t::INFO, "Reading coordinate file:", pdb_name);

             err = MMDBManager->ReadCoorFile(pdb_name.c_str());

             if (err) {

                // try Small-molecule cif
                coot::smcif sm;
                mmdb::Manager *mol = sm.read_sm_cif(pdb_name);

                if (mol) {

                   delete MMDBManager;
                   MMDBManager = mol;
                   err = mmdb::ERROR_CODE(0); // success

                } else {

                   // We also failed to read a small molecule cif, but
                   // write the mmCIF error message.

                   std::cout << "There was an error reading " << pdb_name.c_str() << ". \n";
                   std::cout << "ERROR " << err << " READ: "
                             << mmdb::GetErrorDescription(err) << std::endl;
                   //
                   MMDBManager->GetInputBuffer(error_buf, error_count);
                   if (error_count >= 0) {
                      std::cout << "         LINE #" << error_count << "\n     "
                                << error_buf << std::endl << std::endl;
                   } else {
                      if (error_count == -1) {
                         std::cout << "       CIF ITEM: " << error_buf << std::endl << std::endl;
                      }
                   }
                   asc.read_success = 0; // FAIL
                   asc.read_error_message = error_buf;
                }

             }

             if (! err) {

                atom_name_fix_ups(MMDBManager);

                MMDBManager->PDBCleanup(mmdb::PDBCLEAN_ELEMENT);

                if (verbose_mode) {
                   // std::cout << "INFO:: file " << pdb_name.c_str() << " has been read.\n";
                   logger.log(log_t::INFO, "File", pdb_name, "has been read");
                }
                asc.read_success = 1; // TRUE

                // atom_selection_container.read_error_message = NULL; // its a string
                asc.mol = MMDBManager;
             } else {

#ifdef USE_GEMMI
                mmdb::Manager *mol = file_name_to_manager_via_gemmi(pdb_name); // new mmdb::Manager;
                if (mol) {
                   asc.read_success = 1;
                   asc.mol = mol;
                }
#else
             std::cout << "No GEMMI fallback - sad times " << std::endl;
#endif
             }
          }
       }

       if (MMDBManager) {
          char *str = MMDBManager->GetSpaceGroup();
          if (str) {
             if (verbose_mode) {
                std::string sgrp(str);
                // std::cout << "Spacegroup: " << sgrp << "\n";
                logger.log(log_t::INFO, "File", pdb_name, "has spacegroup", sgrp);
             }
          } else {
             if (verbose_mode) {
                // std::cout << "No Spacegroup found for this PDB file\n";
                logger.log(log_t::INFO, "File", pdb_name, "no spacegroup found");
             }
          }
       }

       // Make handle_read_draw_molecule use make_asc which add the
       // UDD "atom index".
       //
       if (asc.read_success) {
          asc = make_asc(asc.mol);

          // debug atom names
          if (false) {
             for (int i=0; i<asc.n_selected_atoms; i++) {
                std::cout << i << " "
                          << asc.atom_selection[i]->GetChainID() << " "
                          << asc.atom_selection[i]->GetSeqNum() << " :"
                          << asc.atom_selection[i]->name << ":" <<std::endl;
             }
          }

          fix_element_name_lengths(asc.mol); // should not be needed with new mmdb
          fix_away_atoms(asc);
          // fix_wrapped_names(asc); // 20240302-PE remove this. Surely it's no longer needed
                                     // (and it has a memory leak)
       }
    }

    // debug_atom_selection_container(asc);
    return asc;
}

void
atom_selection_container_t::regen_atom_selection() {

   SelectionHandle = mol->NewSelection();
   mol->SelectAtoms (SelectionHandle, 0, "*",
                     mmdb::ANY_RES, // starting resno, an int
                     "*", // any insertion code
                     mmdb::ANY_RES, // ending resno
                     "*", // ending insertion code
                     "*", // any residue name
                     "*", // atom name
                     "*", // elements
                     "*"  // alt loc.
                     );
   mol->GetSelIndex(SelectionHandle, atom_selection, n_selected_atoms);
   UDDAtomIndexHandle = mol->GetUDDHandle(mmdb::UDR_ATOM, "atom index");
   for (int i=0; i<n_selected_atoms; i++)
      atom_selection[i]->PutUDData(UDDAtomIndexHandle, i);
   UDDOldAtomIndexHandle = -1;
}

atom_selection_container_t
coot::mdl_mol_to_asc(const lig_build::molfile_molecule_t &m) {

   return mdl_mol_to_asc(m, 20.0);
}


atom_selection_container_t
coot::mdl_mol_to_asc(const lig_build::molfile_molecule_t &m, float b_factor) {

   atom_selection_container_t asc;

   // set these in the constuctor?
   asc.mol = 0;
   asc.n_selected_atoms = 0;

   if (m.atoms.size()) {
      mmdb::Residue *residue_p = new mmdb::Residue;
      for (unsigned int iat=0; iat<m.atoms.size(); iat++) {
         mmdb::Atom *at = new mmdb::Atom;
         at->SetCoordinates(m.atoms[iat].atom_position.x(),
                            m.atoms[iat].atom_position.y(),
                            m.atoms[iat].atom_position.z(),
                            1.0, b_factor);
         at->SetAtomName(m.atoms[iat].name.c_str());
         at->SetElementName(m.atoms[iat].element.c_str());
         residue_p->AddAtom(at);
      }

      mmdb::Chain *chain_p = new mmdb::Chain;
      mmdb::Model *model_p = new mmdb::Model;

      chain_p->SetChainID("A");
      residue_p->SetResID("UNL", 1, ""); // insertion code of blank

      chain_p->AddResidue(residue_p);
      model_p->AddChain(chain_p);
      mmdb::Manager *mol = new mmdb::Manager;
      mol->AddModel(model_p);
      asc = make_asc(mol);
   }
   return asc;
}

void
fix_element_name_lengths(mmdb::Manager *mol) {

   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (model_p) {
         mmdb::Chain *chain_p;
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            chain_p = model_p->GetChain(ichain);
            if (chain_p) {
               int nres = chain_p->GetNumberOfResidues();
               mmdb::Residue *residue_p;
               mmdb::Atom *at;
               for (int ires=0; ires<nres; ires++) {
                  residue_p = chain_p->GetResidue(ires);
                  if (residue_p) {
                     int n_atoms = residue_p->GetNumberOfAtoms();
                     for (int iat=0; iat<n_atoms; iat++) {
                        at = residue_p->GetAtom(iat);
                        std::string ele(at->element);
                        if (ele.length() == 1) {
                           ele = " " + ele;
                           at->SetElementName(ele.c_str());
                        }
                     }
                  }
               }
            }
         }
      }
   }
}


// Return the number of residue names changed.
//
// Tinker with asc as necessary.
int
fix_nucleic_acid_residue_names(atom_selection_container_t asc) {

   int istat = 0;

   if (asc.n_selected_atoms > 0) {

      int n_models = asc.mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) {

         mmdb::Model *model_p = asc.mol->GetModel(imod);
         // model can legitimately be null if that particular model
         // number was not in the PDB file.
         if (model_p) {
            mmdb::Chain *chain_p;
            // run over chains of the existing mol
            int nchains = model_p->GetNumberOfChains();
            if (nchains <= 0) {
               std::cout << "bad nchains in molecule " << nchains
                         << std::endl;
            } else {
               for (int ichain=0; ichain<nchains; ichain++) {
                  chain_p = model_p->GetChain(ichain);
                  if (chain_p == NULL) {
                     // This should not be necessary. It seem to be a
                     // result of mmdb corruption elsewhere - possibly
                     // DeleteChain in update_molecule_to().
                     std::cout << "NULL chain in ... " << std::endl;
                  } else {
                     int nres = chain_p->GetNumberOfResidues();
                     mmdb::PResidue residue_p;
                     for (int ires=0; ires<nres; ires++) {
                        residue_p = chain_p->GetResidue(ires);
                        std::string residue_name(residue_p->name);

                        if (residue_name == "T" ||
                            residue_name == "U" ||
                            residue_name == "A" ||
                            residue_name == "C" ||
                            residue_name == "G" ||
                            residue_name == "DA" ||
                            residue_name == "DG" ||
                            residue_name == "DT" ||
                            residue_name == "DC") {

                           istat += fix_nucleic_acid_residue_name(residue_p);
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return istat;
}

int fix_nucleic_acid_residue_name(mmdb::Residue *r) {

   int istat=0;

   mmdb::PAtom *residue_atoms;
   int n_residue_atoms;
   bool found_o2_star = 0;

   r->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int i=0; i<n_residue_atoms; i++) {
      std::string atom_name(residue_atoms[i]->name);
      if (atom_name == " O2*") {
         found_o2_star = 1;
         break;
      }
      // I suppose we should call bases with O2' RNA too (not that it
      // does much good because the dictionary will not match the atom
      // names).
      if (atom_name == " O2'") {
         found_o2_star = 1;
         break;
      }
   }

   convert_to_old_nucleotide_atom_names(r);

   std::string res_name = r->name;
   std::string new_name_stub = res_name.substr(0,1);
   if (res_name == "DA" || res_name == "DT" ||
       res_name == "DC" || res_name == "DG")
      new_name_stub = res_name.substr(1,1);

   if (n_residue_atoms > 0)
      istat = 1;

   if (istat == 1) {
      if (found_o2_star) {
         new_name_stub += "r";
      } else {
         new_name_stub += "d";
      }
      r->SetResName(new_name_stub.c_str());
   }
   return istat;
}

// " H5'" -> "H5*1"
// " H5'" -> "H5*1"
// "H5''" -> "H5*2"
void
convert_to_old_nucleotide_atom_names(mmdb::Residue *r) {

   mmdb::PAtom *residue_atoms;
   int n_residue_atoms;
   r->GetAtomTable(residue_atoms, n_residue_atoms);
   for (int i=0; i<n_residue_atoms; i++) {
      std::string atom_name(residue_atoms[i]->name);
      std::string name_orig = atom_name;
      std::string ele(residue_atoms[i]->element);
      char c3 = atom_name[2]; // 3rd char
      char c4 = atom_name[3]; // 4th char
      if (coot::is_hydrogen(ele)) {
         if (c3 == '\'') {
            atom_name[2] = '*';
            if (c4 == '\'')
               atom_name[3] = '2';
            else
               atom_name[3] = '1';
         } else {
            if (c4 == '\'') {
               if (atom_name == " H5'")
                  atom_name = "H5*1";
               else
                  atom_name[3] = '*';
            }
         }
         strncpy(residue_atoms[i]->name, atom_name.c_str(),5);
      } else {
         // if it is not a hydrogen, simply change the prime to a star
         if (c4 == '\'') {
            atom_name[3] = '*';
            strncpy(residue_atoms[i]->name, atom_name.c_str(),5);
         }

         if (atom_name == " OP1") {
            atom_name = " O1P";
            strncpy(residue_atoms[i]->name, atom_name.c_str(),5);
         }
         if (atom_name == " OP2") {
            atom_name = " O2P";
            strncpy(residue_atoms[i]->name, atom_name.c_str(),5);
         }
      }
      // debug
      // std::cout << "from :" << name_orig << ": to :"
      // << atom_name << ":" << std::endl;
   }
}


int
fix_away_atoms(atom_selection_container_t asc) {

   int nat = 0;
   for (int i=0; i<asc.n_selected_atoms; i++) {
      if ((asc.atom_selection[i]->x > 9998.0) &&
          (asc.atom_selection[i]->y > 9998.0) &&
          (asc.atom_selection[i]->z > 9998.0)) {
         asc.atom_selection[i]->x =  0.0;
         asc.atom_selection[i]->y =  0.0;
         asc.atom_selection[i]->z =  0.0;
         nat++;
      }
   }
   return nat;
}

// Return the number of residue names changed.
//
// Tinker with asc as necessary.
int
fix_wrapped_names(atom_selection_container_t asc) {

   int n_changed = 0;
   int uddHnd_old =
      asc.mol->RegisterUDString(mmdb::UDR_ATOM , "initial hydrogen name");
   int uddHnd_new =
      asc.mol->RegisterUDString(mmdb::UDR_ATOM , "new hydrogen name");
//    std::cout << "udd_old: create time " << uddHnd_old << std::endl;
//    std::cout << "udd_new: create time " << uddHnd_new << std::endl;

   // e.g. "3HB " -> " HB3", and "2HG2" -> "HG22"
   for (int i=0; i<asc.n_selected_atoms; i++) {
      // std::string ele(asc.atom_selection[i]->element);

      if (1) {
         std::string atom_name(asc.atom_selection[i]->name);
         if (atom_name[0] == '1' ||
             atom_name[0] == '2' ||
             atom_name[0] == '3' ||
             atom_name[0] == '4' ||
             atom_name[0] == '*') {
            // switch it.
            std::string new_atom_name = atom_name.substr(1,3) + atom_name[0];
            if (atom_name[3] != ' ') {
               if (atom_name[3] == ' ') {
                  new_atom_name = atom_name.substr(1,2) + atom_name[0];
                  new_atom_name += ' ';
               }
               if (atom_name[2] == ' ') {
                  new_atom_name = atom_name.substr(1,1) + atom_name[0];
                  new_atom_name += ' ';
                  new_atom_name += ' ';
               }
            } else {
               // atom_name length is 3 presumably
               new_atom_name = ' ';
               new_atom_name += atom_name.substr(1,2) + atom_name[0];
            }
//               std::cout << "DEBUG:: atom_name switch :" <<  atom_name << ": -> :"
//                         << new_atom_name << ":\n";
            if (uddHnd_old >= 0)
               asc.atom_selection[i]->PutUDData(uddHnd_old,
                                                asc.atom_selection[i]->name);
            if (uddHnd_new >= 0)
               asc.atom_selection[i]->PutUDData(uddHnd_new,
                                                new_atom_name.c_str());
            asc.atom_selection[i]->SetAtomName(new_atom_name.c_str());
            n_changed++;;
          } else {
            // refmac calls it " H "
            if (atom_name == " H0 ") {
               std::string new_atom_name = " H  ";
               if (uddHnd_old >= 0)
                  asc.atom_selection[i]->PutUDData(uddHnd_old,
                                                   asc.atom_selection[i]->name);
               if (uddHnd_new >= 0)
                  asc.atom_selection[i]->PutUDData(uddHnd_new,
                                                   (char *) new_atom_name.c_str());
               asc.atom_selection[i]->SetAtomName(new_atom_name.c_str());
               n_changed++;
            }
         }
      }
   }
   // std::cout << "done hydrogen names " << n_changed << std::endl;
   return n_changed;
}

bool
coot::is_hydrogen(const std::string &ele) {
   if (ele == " H" || ele == " D")
      return true;
   else
      return false;
}

bool
coot::is_deuterium(const std::string &ele) {
   if (ele == " D")
      return true;
   else
      return false;
}

void
debug_atom_selection_container(atom_selection_container_t asc) {

   bool all_atoms = false;

   std::cout << "DEBUG: asc " << "mol=" << asc.mol << std::endl;
   std::cout << "DEBUG: asc " << "n_selected_atoms=" << asc.n_selected_atoms << std::endl;
   std::cout << "DEBUG: asc " << "atom_selection=" << asc.atom_selection << std::endl;
   std::cout << "DEBUG: asc " << "read_error_message=" << asc.read_error_message << std::endl;
   std::cout << "DEBUG: asc " << "read_success=" << asc.read_success << std::endl;

//    cout << "DEBUG: asc " << "cell="
//         << asc.mol->get_cell_p()->a << " "
//         << asc.mol->get_cell_p()->b << " "
//         << asc.mol->get_cell_p()->c << " "
//         << asc.mol->get_cell_p()->alpha << " "
//         << asc.mol->get_cell_p()->beta << " "
//         << asc.mol->get_cell_p()->gamma << std::endl;

//    cout << "DEBUG: asc " << "spacegroup=" << asc.mol->get_cell_p()->spaceGroup
//         << std::endl;

   if (asc.n_selected_atoms > 10) {
      std::cout << "DEBUG start 10 atoms: " << std::endl;
      for (int ii = 0; ii< 10; ii++) {
         std::cout << ii << " " << asc.atom_selection[ii] << " " ;
         mmdb:: Atom *ap = asc.atom_selection[ii];
         std::cout << coot::atom_spec_t(ap) << std::endl;
      }
      std::cout << "DEBUG end 10 atoms: " << std::endl;
      for (int ii = asc.n_selected_atoms - 10; ii< asc.n_selected_atoms; ii++) {
         std::cout << ii << " " << asc.atom_selection[ii] << " " ;
         mmdb:: Atom *ap = asc.atom_selection[ii];
         std::cout << coot::atom_spec_t(ap) << std::endl;
      }
   }

   if (all_atoms) {
      for (int ii = 0; ii< asc.n_selected_atoms; ii++) {
         std::cout << ii << " " << asc.atom_selection[ii] << " " ;
         mmdb:: Atom *ap = asc.atom_selection[ii];
         std::cout << coot::atom_spec_t(ap) << std::endl;
      }
   }
}

atom_selection_container_t
make_asc(mmdb::Manager *mol, bool transfer_atom_index_flag) {

   atom_selection_container_t asc;
   asc.mol = mol;

   asc.SelectionHandle = mol->NewSelection();
   asc.mol->SelectAtoms (asc.SelectionHandle, 0, "*",
                     mmdb::ANY_RES, // starting resno, an int
                     "*", // any insertion code
                     mmdb::ANY_RES, // ending resno
                     "*", // ending insertion code
                     "*", // any residue name
                     "*", // atom name
                     "*", // elements
                     "*"  // alt loc.
                     );

   asc.mol->GetSelIndex(asc.SelectionHandle, asc.atom_selection, asc.n_selected_atoms);

   int uddHnd = mol->RegisterUDInteger(mmdb::UDR_ATOM, "atom index");
   if (false)
      std::cout << "debug:: in make_asc(): uddHnd " << uddHnd << " for 'atom index' for mol "
                << mol << std::endl;
   if (uddHnd < 0) {
      std::cout << "ERROR:: ----------------- atom index registration failed.\n";
   } else {
      // std::cout << "in make_asc() saving UDDAtomIndexHandle " << uddHnd << std::endl;
      asc.UDDAtomIndexHandle = uddHnd;
      for (int i=0; i<asc.n_selected_atoms; i++) {
         int status = asc.atom_selection[i]->PutUDData(uddHnd,i);
      }
   }
   asc.read_error_message = "No error";
   asc.read_success = 1;
   asc.UDDOldAtomIndexHandle = -1;

   if (transfer_atom_index_flag) {
      int udd_atom_index_handle = mol->GetUDDHandle(mmdb::UDR_ATOM, "atom index");
      asc.UDDOldAtomIndexHandle = udd_atom_index_handle;
   }

   return asc;
}

void
atom_selection_container_t::add_old_atom_indices() {

   if (mol) {
      UDDOldAtomIndexHandle = mol->RegisterUDInteger(mmdb::UDR_ATOM, "old atom index");
      for (int i=0; i<n_selected_atoms; i++)
         atom_selection[i]->PutUDData(UDDOldAtomIndexHandle, i);
   }

}


void
atom_selection_container_t::fill_links(mmdb::Manager *mol_other) {

   if (mol_other) {
      mmdb::Model *model_p = mol_other->GetModel(1);
      if (model_p) {
         unsigned int n_links = model_p->GetNumberOfLinks();
         for (unsigned int i=1; i<=n_links; i++) {
            mmdb::Link *ref_link = model_p->GetLink(i);
            if (ref_link) {
               mmdb::Link l(*ref_link);
               links.push_back(l);
            }
         }
      }
   }
}



atom_selection_container_t read_standard_residues() {

   std::string standard_env_dir = "COOT_STANDARD_RESIDUES";
   atom_selection_container_t standard_residues_asc;

   const char *filename = getenv(standard_env_dir.c_str());
   if (! filename) {

      // std::string standard_file_name = PKGDATADIR;
      std::string standard_file_name = coot::package_data_dir();
      standard_file_name += "/";
      standard_file_name += "standard-residues.pdb";

      struct stat buf;
      int status = stat(standard_file_name.c_str(), &buf);
      if (status != 0) { // standard-residues file was not found in
                         // default location either...
         std::cout << "WARNING: environment variable for standard residues ";
         std::cout << standard_env_dir << "\n";
         std::cout << "         is not set.";
         std::cout << " Mutations will not be possible\n";
         // mark as not read then:
         standard_residues_asc.read_success = 0;
         // std::cout << "DEBUG:: standard_residues_asc marked as
         // empty" << std::endl;
      } else {
         // stat success:
         standard_residues_asc = get_atom_selection(standard_file_name, false, true, false);
      }
   } else {
      standard_residues_asc = get_atom_selection(filename, false, true, false);
   }

   return standard_residues_asc;
}

void
atom_selection_container_t::debug_write_pdb() {

   if (mol) {
      std::string fn = "debug-" + std::to_string(user_data) + ".pdb";
      mol->WritePDBASCII(fn.c_str());
   }

}


// return an estimate of the molecule diameter
float
coot::get_molecule_diameter(const atom_selection_container_t &asc) {

   // c.f. radius_of_gyration

   float f = -1;

   int n_max = asc.n_selected_atoms;

   stats::single s;
   for (unsigned int i=0; i<1000; i++) {
      float f1 = coot::util::random_f();
      float f2 = coot::util::random_f();
      int ff_1 = static_cast<float>(n_max) * f1;
      int ff_2 = static_cast<float>(n_max) * f2;
      // std::cout << "f1 " << f1 << " f2 " << f2 << " ff_1 " << ff_1 << " ff_2 " << ff_2 << std::endl;
      if (f1 < 1.0) {
         if (f2 < 1.0) {
            int idx_1 = static_cast<int>(ff_1);
            int idx_2 = static_cast<int>(ff_2);
            if (idx_1 != idx_2) {
               mmdb:: Atom *at_1 = asc.atom_selection[idx_1];
               mmdb:: Atom *at_2 = asc.atom_selection[idx_2];
               float dx = at_2->x - at_1->x;
               float dy = at_2->y - at_1->y;
               float dz = at_2->z - at_1->z;
               float dd = dx*dx + dy*dy + dz*dz;
               float d = std::sqrt(dd);
               s.add(d);
            }
         }
      }
   }

   if (s.size() > 10) {
      // std::cout << "info:: idx 10 highest " << s.get_ith_highest(10) << " idx 10 lowest " << s.get_ith_lowest(10) << std::endl;
      f = s.get_ith_highest(10);
   }

   return f;

}


