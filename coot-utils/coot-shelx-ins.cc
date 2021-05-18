/*
 * 
 * Copyright 2005, 2006, 2007 by The University of York
 * Copyright 2008 by The University of Oxford
 * Copyright 2013, 2014, 2015 by Medical Research Council
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */


#include <string.h>
#include <string>

#include <mmdb2/mmdb_manager.h>
#include "clipper/core/spacegroup.h"
#include "clipper/core/coords.h"

#include "utils/coot-utils.hh"
#include "coot-coord-utils.hh"
#include "coot-shelx.hh"
#include <iostream>
#include <iomanip>
#include <fstream>

#include "compat/coot-sysdep.h"

enum { NONE=0, ONE_HALF=6, ONE_THIRD=4, ONE_QUARTER=3, ONE_SIXTH=2, TWO_THIRDS=8,
       THREE_QUARTERS=9, FIVE_SIXTHS=10, MINUS_ONE_HALF= -6, MINUS_ONE_THIRD = -4,
       MINUS_ONE_QUARTER = -3, MINUS_ONE_SIXTH = -2, MINUS_TWO_THIRDS = -8,
       MINUS_THREE_QUARTERS = -9 , MINUS_FIVE_SIXTHS  = -10};
       
coot::ShelxIns::ShelxIns(const std::string &filename) {

   init();
   coot::shelx_read_file_info_t p =
      read_file(filename); // sets filled_flag and udd_afix_handle

}


// Shelx Free Variables
//
// Free Variables are specified on the FVAR line.  The FVAR lines come
// just before the atom lines.
// 
// The free variables are used typically for occupancy refinement of
// disordered side-chains, and can be used in an any order (which is
// useful).
//
// So, if we want to split a side chain of a shelx molecule, how do we
// do that?
// 
// If we decode that occ from 11 (1 + 1) Fixed at 1.
// 
// When we split that we need to introduce a new free variable
// (default 0.5? - or whatever was set in the occ chooser?).  Say that
// free variable is number 23.
// 
// When it comes to writting the PARTs we need to set the occ to 23.0
// for PART A and -23.0 for PART B.
//
// We need to know that that residue has been split 



// Return a pair: status (0: bad), udd_afix_handle (-1 bad)
// 
coot::shelx_read_file_info_t
coot::ShelxIns::read_file(const std::string &filename) {

   coot::shelx_card_info_t card;
   int istate = 0;
   udd_afix_handle = -1;
   mmdb::Manager *mol = NULL;
   std::ifstream f(filename.c_str());
   int latt = 0;  // special not-set value
   int resi_count = 0; // use to determine if this is a protein

   if (f) {

      istate = 1; // success
      filled_flag = 1;
      short int in_residue_flag = 0;
      std::string current_res_name;
      int current_res_no = 0; // shelx default 
      clipper::Spacegroup space_group;
      mmdb::InitMatType();
      mol = new mmdb::Manager;
      mmdb::Model *model = new mmdb::Model;
      mmdb::Chain *chain = new mmdb::Chain; // possibly ends up empty?
      std::map<std::string, mmdb::Chain *> chain_map;
      chain->SetChainID("");
      chain_map[""] = chain;
      model->AddChain(chain);
      std::vector<mmdb::Atom *> atom_vector;
      std::vector<std::string> symm_vec;
      std::vector<double> cell_local;
      std::string altconf;
      short int encountered_atoms_flag = 0;
      short int post_atoms_flag = 0; // we have got past atoms?, old, not useful?
      short int post_END_flag = 0;   // END marks the end of atoms, I hope.
      int udd_non_riding_atom_flag_handle = mol->RegisterUDInteger(mmdb::UDR_ATOM, "non_riding_atom");

      // class variables
      udd_afix_handle = mol->RegisterUDInteger(mmdb::UDR_ATOM, "shelx afix");
      udd_riding_atom_negative_u_value_handle = mol->RegisterUDReal(mmdb::UDR_ATOM, "riding_atom_negative_u");
      
      bool have_udd_atoms = false;
      if (udd_afix_handle>=0)  {
         have_udd_atoms = true;
      }
      // float u_to_b = 8.0 * M_PI * M_PI;  // perhaps this should be a function
      int current_afix = -1; // undefined initially.  AFIX is
                             // independent of RESI and should not be
                             // reset for every residue (as had been
                             // the practice 20060102).
      int nlines = 0;


      // GMS says:

      // FVAR is the last card to come before atom cards.  There can
      // be any number of FVAR cards and there can be any number of
      // variables on the FVAR cards.

      // The end of the atoms is marked by END or HKLF cards.  The Q
      // values might be of interest - e.g. they are peaks that might
      // be waters.
         
      while (!f.eof()) {
         card = read_card(f);
         nlines++;
         // std::cout << card.card << std::endl;
         
         if (card.words.size() > 0) {
            std::string card_word_0 = util::upcase(card.words[0]);
            if (card_word_0 == "RESI") {
               encountered_atoms_flag = 1;
               resi_count++;
               if (atom_vector.size() > 0) { 
                  mmdb::Residue *residue =
                     add_shelx_residue(atom_vector,
                                       current_res_name,
                                       current_res_no);
                  chain->AddResidue(residue);
                  // std::cout << "adding residue to chain " << residue_spec_t(residue) << std::endl;
               }
               // now update
               if (card.words.size() > 1) {

                  // res_no_string can be "1001" or "A:1001"
                  std::string res_no_string = card.words[1];
                  // first try a simple string -> int
                  try {
                     current_res_no = util::string_to_int(res_no_string);
                  }
                  catch (const std::runtime_error &rte) {
                     // OK, was it "A:1001" style?
                     // find the colon and try string -> int on the rest of the string
                     std::size_t p = res_no_string.find(':');
                     if (p != std::string::npos) {
                        std::string res_chain_id = res_no_string.substr(0,p);
                        std::string r_string = res_no_string.substr(p+1);
                        std::map<std::string, mmdb::Chain *>::iterator it = chain_map.find(res_chain_id);
                        if (it == chain_map.end()) {
                           chain = new mmdb::Chain;
                           chain->SetChainID(res_chain_id.c_str());
                           model->AddChain(chain);
                           chain_map[res_chain_id] = chain;
                           std::cout << "debug:: ## new chain wiith chain-id " << res_chain_id << std::endl;
                        } else {
                           chain = it->second;
                        }
                        try {
                           current_res_no = util::string_to_int(r_string);
                        }
                        catch (const std::runtime_error &rte) {
                           std::cout << "WARNING failed to parse residue from line " << card.card << std::endl;
                        }
                     } else {
                        std::cout << "WARNING failed to parse residue from line " << card.card << std::endl;
                     }
                  }
               }
               if (card.words.size() > 2) 
                  current_res_name = card.words[2];
               else
                  current_res_name = "";
               atom_vector.clear();
            } else {

               if (card_word_0 == "AFIX") {
                  if (card.words.size() > 1) {
                     current_afix = atoi(card.words[1].c_str());
                  }
               } else {
                  if (card_word_0 == "PART") {
                     // 20081008 allow for a site occupancy factor on the PART line
                     // but then throw it away.  Not good. Should fix at some stage.
                     if (card.words.size() == 2 || card.words.size() == 3) {
                        if (card.words[1] == "0") {
                           altconf = "";
                        } else {
                           if (card.words[1] == "1") {
                              altconf = "A";
                           } else {
                              if (card.words[1] == "2") {
                                 altconf = "B";
                              } else {
                                 if (card.words[1] == "3") {
                                    altconf = "C";
                                 } else {
                                    if (card.words[1] == "4") {
                                       altconf = "D";
                                    } else {
                                       if (card.words[1] == "5") {
                                          altconf = "E";
                                       } else {
                                          if (card.words[1] == "6") {
                                             altconf = "F";
                                          } else {
                                             if (card.words[1] == "7") {
                                                altconf = "G";
                                             } else {
                                                if (card.words[1] == "8") {
                                                   altconf = "H";
                                                } else {
                                                   
                                                   // negative PARTS: shelxl
                                                   // generation of special position
                                                   // constraints suppressed.
                                                   
                                                   if (card.words[1] == "-1") {
                                                      altconf = "a";
                                                   } else { 
                                                      if (card.words[1] == "-2") {
                                                         altconf = "b";
                                                      } else { 
                                                         if (card.words[1] == "-3") {
                                                            altconf = "c";
                                                         } else {
                                                            if (card.words[1] == "-4") {
                                                               altconf = "d";
                                                            } else {
                                                               if (card.words[1] == "-5") {
                                                                  altconf = "e";
                                                               } else {
                                                                  if (card.words[1] == "-6") {
                                                                     altconf = "f";
                                                                  } else {
                                                                     if (card.words[1] == "-7") {
                                                                        altconf = "g";
                                                                     } else {
                                                                        if (card.words[1] == "-8") {
                                                                           altconf = "h";
                                                                        }
                                                                     }
                                                                  }
                                                               }
                                                            }
                                                         }
                                                      }
                                                   }
                                                }
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }
                  } else {
                        
                     // Not in a residue, e.g. TITL, CELL, SYMM, LATT
                     // etc.
                     // Hmmm, but small molecule atoms are not in a
                     // residue either.

                     // It was a non-atom keyword after atoms?

                     if (encountered_atoms_flag) {
                        if (post_atoms_flag) {
                           post_atom_lines.push_back(card.card);
                        }
                     } 
                     // Push back all pre-atom lines, except FVARs
                     // which are handled differently.
                     bool is_fvar_line = false;
                     bool is_sfac_line = false;
                     bool is_disp_line = false;
                     if (card.words.size() > 0)
                        if (card_word_0 == "FVAR")
                           is_fvar_line = 1; // flag for later

                     if (card_word_0 == "TITL") {
                        title = card.card;
                     } else {
                        if (card_word_0 == "CELL") {
                           if (card.words.size() > 7) {
                              have_cell_flag = 1;
                              cell_local.resize(6);
                              cell_local[0] = atof(card.words[2].c_str());
                              cell_local[1] = atof(card.words[3].c_str());
                              cell_local[2] = atof(card.words[4].c_str());
                              cell_local[3] = atof(card.words[5].c_str());
                              cell_local[4] = atof(card.words[6].c_str());
                              cell_local[5] = atof(card.words[7].c_str());
                              clipper::Cell_descr cell_d(cell_local[0], cell_local[1], cell_local[2], 
                                                         clipper::Util::d2rad(cell_local[3]),
                                                         clipper::Util::d2rad(cell_local[4]),
                                                         clipper::Util::d2rad(cell_local[5]));
                              cell.init(cell_d); // cell is class member
                           }
                        } else { 
                           if (card_word_0 == "UNIT") {
                              for (unsigned int i=1; i<card.words.size(); i++)
                                 unit.push_back(atoi(card.words[i].c_str())); // number of atoms or each sfac type
                           } else {
                              if (card_word_0 == "SFAC") {
                                 std::cout << "SFAC LINE: " << card.card << std::endl;
                                 for (unsigned int i=1; i<card.words.size(); i++)
                                    sfac.push_back(card.words[i]); // atoms
                                 // std::cout << "DEBUG:: sfac is now of size " << sfac.size() << std::endl;
                                 is_sfac_line = 1;
                              } else {
                                 if (card_word_0 == "SYMM") {
                                    symm_cards.push_back(card.card); // save for output
                                    std::string s;
                                    for (unsigned int i=1; i<card.words.size(); i++) { 
                                       s += card.words[i];
                                       s += " ";
                                    }
                                    symm_vec.push_back(s);

                                 } else {

                                    if (card_word_0 == "DISP") {
                                       disp_cards.push_back(card.card);
                                       is_disp_line = true;

                                    } else {
                                       if (card_word_0 == "HKLF") {
                                          post_atoms_flag = 1;
                                          post_atom_lines.push_back(card.card); // special case, first post atom line
                                          mmdb::Residue *residue = add_shelx_residue(atom_vector,
                                                                                     current_res_name,
                                                                                     current_res_no);
                                          chain->AddResidue(residue);
                                          atom_vector.clear();
                                       } else {
                                          if (card_word_0.substr(0, 4) == "WGHT") {
                                             // handle WGHT
                                             // post_atoms_flag = 1; // Definately not, e.g. T.RES
                                          } else {
                                             if (card_word_0.substr(0, 3) == "END") {
                                                // handle END
                                                mmdb::Residue *residue = add_shelx_residue(atom_vector,
                                                                                           current_res_name,
                                                                                           current_res_no);
                                                chain->AddResidue(residue);
                                                atom_vector.clear();
                                                post_atoms_flag = 1;
                                                post_END_flag = 1;

                                             } else {

                                                if (card_word_0.substr(0, 4) == "ZERR") { // zerror 
                                             
                                                } else {

                                                   if (card_word_0.substr(0, 3) == "REM") { // REM
                                                      // Special case:
                                                      // we have a REM after atoms but before HKL 4 line
                                                      // e.g. insulin.res from GMS
                                                      if (encountered_atoms_flag)
                                                         if (!post_atoms_flag)
                                                            post_atom_lines.push_back(card.card);
                                                
                                                   } else {
                                                      if (card_word_0.substr(0, 4) == "FVAR") {

                                                         save_fvars(card); 

                                                      } else {

                                                         if ( (card_word_0.substr(0, 4) == "DEFS") || // DEFS and others
                                                              (card_word_0.substr(0, 4) == "CGLS") ||
                                                              (card_word_0.substr(0, 4) == "SHEL") ||
                                                              (card_word_0.substr(0, 4) == "FMAP") ||
                                                              (card_word_0.substr(0, 4) == "SIZE") ||
                                                              (card_word_0.substr(0, 4) == "STIR") ||
                                                              (card_word_0.substr(0, 4) == "TEMP") ||
                                                              (card_word_0.substr(0, 4) == "BLOC") ||
                                                              (card_word_0.substr(0, 4) == "SADI") ||
                                                              // (card_word_0.substr(0, 4) == "DISP") ||
                                                              (card_word_0.substr(0, 4) == "SPEC") ||
                                                              (card_word_0.substr(0, 4) == "PLAN") ||
                                                              (card_word_0.substr(0, 4) == "LIST") ||
                                                              (card_word_0.substr(0, 4) == "FREE") ||
                                                              (card_word_0.substr(0, 4) == "HTAB") ||
                                                              (card_word_0.substr(0, 4) == "DELU") ||
                                                              (card_word_0.substr(0, 4) == "SIMU") ||
                                                              (card_word_0.substr(0, 4) == "CONN") ||
                                                              (card_word_0.substr(0, 4) == "EQIV") ||
                                                              (card_word_0.substr(0, 4) == "BUMP") ||
                                                              (card_word_0.substr(0, 4) == "MORE") ||
                                                              (card_word_0.substr(0, 4) == "ISOR") ||
                                                              (card_word_0.substr(0, 4) == "MERG") ||
                                                              (card_word_0.substr(0, 4) == "TREF") ||
                                                              (card_word_0.substr(0, 4) == "ACTA") ||
                                                              (card_word_0.substr(0, 4) == "TWIN") ||
                                                              (card_word_0.substr(0, 4) == "OMIT") ||
                                                              (card_word_0.substr(0, 4) == "SWAT") ||
                                                              (card_word_0.substr(0, 4) == "ANIS") ||
                                                              (card_word_0.substr(0, 4) == "BASF") ||
                                                              (card_word_0.substr(0, 4) == "SAME") ||
                                                              (card_word_0.substr(0, 4) == "MOLE") ||
                                                              (card_word_0.substr(0, 4) == "BIND") ||
                                                              (card_word_0.substr(0, 4) == "L.S.") ||
                                                              (card_word_0.substr(0, 4) == "SUMP") ||
                                                              (card_word_0.substr(0, 4) == "BOND") ||
                                                              (card_word_0.substr(0, 4) == "RIGU") ||
                                                              (card_word_0.substr(0, 4) == "CONF") ||
                                                              (card_word_0.substr(0, 4) == "MPLA") ||
                                                              (card_word_0.substr(0, 4) == "HOPE") ||
                                                              (card_word_0.substr(0, 4) == "EXTI") ||
                                                              (card_word_0.substr(0, 4) == "XNPD") ||
                                                              (card_word_0.substr(0, 4) == "WPDB")) { 
                                                         } else {
                                                            if (card_word_0.substr(0, 4) == "LATT") {
                                                               std::cout << "LATT LINE: " << card.card << std::endl;
                                                               latt = atoi(card.words[1].c_str());  // potential crash here
                                                            } else { 
                                                               if ( (card_word_0.substr(0, 4) == "DFIX" ) || 
                                                                    (card_word_0.substr(0, 5) == "DFIX_")) {
                                                               } else {
                                                                  if ((card_word_0.substr(0, 4) == "FLAT") ||
                                                                      (card_word_0.substr(0, 5) == "SADI_")) {
                                                               
                                                                  } else {
                                                                     if (card_word_0.substr(0, 5) == "CHIV_") {
                                                                     } else {
                                                                        if ( (card_word_0.substr(0, 4) == "DANG" ) ||
                                                                             (card_word_0.substr(0, 5) == "DANG_") ||
                                                                             (card_word_0.substr(0, 4) == "RTAB")  ||
                                                                             (card_word_0.substr(0, 5) == "RTAB_") ) {
                                                                        } else {

                                                                           if (card.words.size() <= 4) {
                                                                              //
                                                                              // A dos file has ^Ms for blank lines.  We don't want to
                                                                              // say that those are bad atoms
                                                                              short int handled_dos = 0;
                                                                              if (card.words.size() == 1) {
                                                                                 if (card_word_0.length() == 1) {
                                                                                    handled_dos = 1;
                                                                                 }
                                                                              }
                                                                              if (!handled_dos) 
                                                                                 std::cout << "WARNING:: BAD ATOM line " << nlines << " " 
                                                                                           << card.words.size()
                                                                                           << " field(s) :" << card.card << ":"
                                                                                           << std::endl;
                                                                           
                                                                           } else {
                                                                              if (! post_END_flag) { 
                                                                                 // it's an atom
                                                                                 // std::cout << "it's an atom! " << std::endl;
                                                                                 mmdb::Atom *at = make_atom(card, altconf,
                                                                                                            udd_afix_handle,
                                                                                                            udd_non_riding_atom_flag_handle,
                                                                                                            udd_riding_atom_negative_u_value_handle,
                                                                                                            have_udd_atoms, current_afix,
                                                                                                            cell, atom_vector);
                                                                                 if (at)
                                                                                    atom_vector.push_back(at);
                                                                                 else
                                                                                    std::cout << "WARNING:: BAD ATOM on line #" << nlines << " " 
                                                                                              << " #fields: " << card.words.size()
                                                                                              << " field(s) :" << card.card << ":"
                                                                                              << std::endl;
                                                                                 encountered_atoms_flag = 1; // stop adding to pre_atom lines
                                                                              }
                                                                           }
                                                                        }
                                                                     }
                                                                  }
                                                               }
                                                            }
                                                         }
                                                      }
                                                   }
                                                }
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
                           }
                        }
                     }
                     if ((!is_fvar_line) && (!encountered_atoms_flag) && (!is_sfac_line) && (!is_disp_line))
                        pre_atom_lines.push_back(card.card);
                  }
               }
            }
         } else {
            if (encountered_atoms_flag) {
               if (post_atoms_flag) {
                  post_atom_lines.push_back(card.card);
               }
            } else {
               pre_atom_lines.push_back(card.card);
            }
            in_residue_flag = 0;
         }
      }
      // we get garbage on the last line from the read for some
      // reason.  Let's pop it off:
      if (!post_atom_lines.empty()) 
         post_atom_lines.pop_back();

      if (cell_local.size() == 6) { 
         clipper::Cell_descr cell_d(cell_local[0], cell_local[1], cell_local[2], 
                                    clipper::Util::d2rad(cell_local[3]),
                                    clipper::Util::d2rad(cell_local[4]),
                                    clipper::Util::d2rad(cell_local[5]));
         cell.init(cell_d);
         mol->SetCell(cell_local[0], cell_local[1], cell_local[2],
                      cell_local[3], cell_local[4], cell_local[5]);
         std::cout << "INFO:: CELL set to "
                   << cell_local[0] << " "
                   << cell_local[1] << " "
                   << cell_local[2] << " "
                   << cell_local[3] << " "
                   << cell_local[4] << " "
                   << cell_local[5] << "\n";
         have_cell_flag = 1;
      } else {
         std::cout << "Cell not set" << std::endl;
      }


      std::string symmetry_ops("");

      // There was no LATT card, so the default is centro-symmetric
      // 
      // We need to multiply each of the symm ops with -1
      // 
      // What a pain.
      //
      // There is something else that is going on here, but I have
      // forgotten what it was.  GMS does it "in shorthand", but the
      // symop.lib has expanded form using symops.  Possibly something
      // to do with the lattice (yes, I think that's it).
      //
      // The lattice/cenop/symop expansion is not to be done here.
      // It should be a coot-utils function.
      //
      // If the LATT is negative, that means non-centrosymmetric
      // If the LATT is positive, that's centrosymmetric.
      // If there is no LATT given, the default is 1

      if (latt == 0)
         latt = 1; // P1 bar
      std::vector<std::string> cs = clipper_symm_strings(symm_vec, latt);

//       for (unsigned int ics=0; ics<cs.size(); ics++)
//           std::cout << "DEBUG:: clipper cs " << ics << " " << cs[ics] << std::endl;

// there were 2 "X,Y,Z"s and clipper was failing to init the space group
//       if (symm_vec.size() > 0)
//          symmetry_ops = "X,Y,Z ; ";

         
      for (unsigned int isym=0; isym<cs.size(); isym++) { 
          symmetry_ops += cs[isym];
         symmetry_ops += " ; ";
      }
      if (symmetry_ops != "") {

         bool spacegroup_ok = true;
         try {
            space_group.init(clipper::Spgr_descr(symmetry_ops, clipper::Spgr_descr::Symops));
         } catch (const clipper::Message_base &exc) {
            std::cout << "Oops, trouble.  No such spacegroup\n";
            spacegroup_ok = false;
         }

         if (spacegroup_ok) {
            // std::cout << "INFO:: set space group to \"" << space_group.descr().symbol_xhm() << "\""
            // << std::endl;
            mmdb::cpstr sg = space_group.descr().symbol_xhm().c_str();
            mol->SetSpaceGroup(sg);
            char *spg = mol->GetSpaceGroup();
            if (spg) { 
               std::string sgrp(spg);
               std::cout << "READ-INS:: Spacegroup: \"" << sgrp << "\"\n";
            } else {
               std::cout << "READ-INS:: No Spacegroup found in res file\n";
               // Is that because mmdb doesn't know about SYMINFO?
               char *si = getenv("SYMINFO");
               if (! si) {
                  std::cout << "   possible cause: SYMINFO not set.\n";
               }
            }
         }
      }

      std::cout << "INFO:: read_file() chain with chain id " << chain->GetChainID() << " has "
                << chain->GetNumberOfResidues() << " residues" << std::endl;
      mol->AddModel(model);

      if (cell_local.size() != 6) { // found a cell?
         std::cout << "WARNING:: no cell found in shelx file\n";
      } else {
         int n_models = mol->GetNumberOfModels();
         for (int imod=1; imod<=n_models; imod++) { 
      
            mmdb::Model *model_p = mol->GetModel(imod);
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
                     std::cout << "NULL chain in read_file() " << imod << std::endl;
                  } else { 
                     int nres = chain_p->GetNumberOfResidues();
                     mmdb::PResidue residue_p;
                     mmdb::Atom *at;
                     for (int ires=0; ires<nres; ires++) { 
                        residue_p = chain_p->GetResidue(ires);
                        int n_atoms = residue_p->GetNumberOfAtoms();

                        for (int iat=0; iat<n_atoms; iat++) {
                           at = residue_p->GetAtom(iat);
                           if (! at) {
                              std::cout << "Trapped a null atom for atom "
                                        << iat << " of residue " << coot::residue_spec_t(residue_p)
                                        << std::endl;
                           } else {
                              if (false) // debug
                                 std::cout << "Mol Hierarchy atom: " << iat << " "
                                           << " " << at->name << " "
                                           << at->GetResName() << " " << at->GetSeqNum() << " "
                                           << at->x << " " <<  at->y << " " << at->z << std::endl;
                              clipper::Coord_frac pf(at->x, at->y, at->z);
                              clipper::Coord_orth po = pf.coord_orth(cell);
                              at->x = po.x();
                              at->y = po.y();
                              at->z = po.z();
                           }
                        }
                     }
                  }
               }
            }
         }
      }
      mol->FinishStructEdit();
      mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      // mol->WritePDBASCII("testing-shelx-read-file.pdb");
      // write_ins_file(mol, "just-read-this.res");
   }

   // SHELX ins files now can have chain ids attached to the residue numbers
   // If that was the case, the molecule is already split into chains and we
   // don't need to unshelx.  Which did we read? If there are more than 1 chains
   // or the chain-id is not blank then we read new style
   //
   mmdb::Manager *shelx_mol = 0;
   if (mol) { // maybe file not found?

      // shall we turn off needs_unshelx?
      bool needs_unshelx = mol_needs_shelx_transfer(mol);

      if (needs_unshelx) 
         shelx_mol = unshelx(mol);
      else
         shelx_mol = mol;
   }

   if (shelx_mol) {
      shelx_read_file_info_t ri(istate, udd_afix_handle, shelx_mol);
      if (resi_count > 10)
         ri.is_protein_flag = 1;

      if (false) { // debugging block
         std::cout << "debug:: shelx read_file() returns mol: " << std::endl;
         int imod = 1;
         mmdb::Model *model_p = mol->GetModel(imod);
         if (! model_p) {
            std::cout << "debug:: shelx read_file() No model for 1 " << std::endl;
         } else {
            int n_chains = model_p->GetNumberOfChains();
            std::cout << "debug:: shelx read_file() n_chains: " << n_chains << std::endl;
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               int nres = chain_p->GetNumberOfResidues();
               for (int ires=0; ires<nres; ires++) {
                  mmdb::Residue *residue_p = chain_p->GetResidue(ires);
                  std::cout << "debug:: shelx read_file() " << residue_spec_t(residue_p) << std::endl;
               }
            }
         }
      }

      return ri;
   } else {
      return coot::shelx_read_file_info_t(istate, udd_afix_handle, mol);
   }
}


mmdb::Atom *
coot::ShelxIns::make_atom(const coot::shelx_card_info_t &card, const std::string &altconf,
                          int udd_afix_handle_in,
                          int udd_non_riding_atom_flag_handle_in,
                          int udd_riding_atom_negative_u_value_handle_in,
                          bool have_udd_atoms, int current_afix,
                          clipper::Cell &cell_in, const std::vector<mmdb::Atom *> &atom_vector) const {

   mmdb::Atom *at = new mmdb::Atom;  // returned

   // local
   int sfac_index;
   float u_to_b = 8.0 * M_PI * M_PI;  // perhaps this should be a function
   
   // std::cout << "DEBUG:: new atom for " <<  card.card << std::endl;
   
   sfac_index = atoi(card.words[1].c_str());
   std::string element = make_atom_element(card.words[0].c_str(), sfac_index);
   if (element == "ERROR-in-SFAC") {
      std::cout << "WARNING:: problem making element - ignoring atom" << std::endl;
      delete at;
      at = 0;
   } else { 
   
      at->SetAtomName(make_atom_name(card.words[0].c_str(), element).c_str());
      at->x = atof(card.words[2].c_str());
      at->y = atof(card.words[3].c_str());
      at->z = atof(card.words[4].c_str());
      float occupancy = 1.0;
      float b_synth= 10.0;

      try {
         if (card.words.size() > 5)
            occupancy = util::string_to_float(card.words[5]); // 11.0
         
         at->SetCoordinates(util::string_to_float(card.words[2].c_str()),
                            util::string_to_float(card.words[3].c_str()),
                            util::string_to_float(card.words[4].c_str()),
                            occupancy, b_synth);
         at->SetElementName(element.c_str());
         strncpy(at->altLoc, altconf.c_str(), 2);
      }
      catch (const std::runtime_error &rte) {
         // do nothing
      }
      if (card.words.size() >= 5) {
         
//          for (unsigned int iword=0; iword<card.words.size(); iword++) {
//             std::cout << "    " << iword << " " << card.words[iword] << std::endl;
//          }

         // What makes an atom non-riding?
         // 
         // I'm guessing that there is an anisotropic U or sing U > 0
         
         if (card.words.size() > 6) {
            if (card.words.size() < 8) {
               // isotropic temperature factor
               mmdb::realtype u_factor_from_card = atof(card.words[6].c_str());
               if (u_factor_from_card > 0.0 ) { 
                  at->tempFactor = u_to_b * u_factor_from_card;
                  at->WhatIsSet = at->WhatIsSet | 4; // is isotropic
                  at->PutUDData(udd_non_riding_atom_flag_handle_in, 1);
               } else {
                  // negative U:
                  // 
                  // If U is between -0.5 and -5.0 then apply riding U rule
                  // 
                  if ((u_factor_from_card <= -0.5) && (u_factor_from_card >= -5.0)) {
                     // Find previous non-riding atom and use that to determine b for this atom

                     mmdb::Atom *prev = previous_non_riding_atom(atom_vector, udd_non_riding_atom_flag_handle_in);
                     if (prev) {
                        int status = at->PutUDData(udd_riding_atom_negative_u_value_handle_in, u_factor_from_card);
                        at->tempFactor = prev->tempFactor * -u_factor_from_card;
                     } else {
                        // Don't know what to do.  Does this ever happen?
                        at->tempFactor = u_factor_from_card;
                     } 
                     // 
                  } else {
                     // Don't know what to do.  Does this ever happen?
                     at->tempFactor = u_factor_from_card;
                  }
               }
            } else {
               if (card.words.size() > 11) {
                  // anisotropic temperature factor
                  at->u11 = atof(card.words[ 6].c_str());
                  at->u22 = atof(card.words[ 7].c_str());
                  at->u33 = atof(card.words[ 8].c_str());
                  at->u23 = atof(card.words[ 9].c_str());
                  at->u13 = atof(card.words[10].c_str());
                  at->u12 = atof(card.words[11].c_str());

                  double a = cell_in.a();
                  double b = cell_in.b();
                  double c = cell_in.c();

                  // Now othogonalize the U values:
                  // clipper::U_aniso_frac ocaf(at->u11, at->u22, at->u33,
                  //                            at->u12, at->u13, at->u23);
                  clipper::U_aniso_frac caf(at->u11/(a*a), at->u22/(b*b), at->u33/(c*c),
                                            at->u12/(a*b), at->u13/(a*c), at->u23/(b*c));
                  clipper::U_aniso_orth cao = caf.u_aniso_orth(cell_in);

                  at->u11 = cao(0,0);
                  at->u22 = cao(1,1);
                  at->u33 = cao(2,2);
                  at->u12 = cao(0,1);
                  at->u13 = cao(0,2);
                  at->u23 = cao(1,2);

//                   std::cout << "DEBUG::  pre-orthog:\n" << ocaf.format() << std::endl;
//                   std::cout << "DEBUG:: post-orthog:\n" <<  cao.format() << std::endl;

                  at->WhatIsSet |= mmdb::ASET_Anis_tFac; // is anisotropic
                  float u_synth = (at->u11 + at->u22 + at->u33)/3.0;
                  at->WhatIsSet |= mmdb::ASET_tempFactor; // has synthetic B factor
                  at->tempFactor = 8.0 * M_PI * M_PI * u_synth;
                  at->PutUDData(udd_non_riding_atom_flag_handle_in, 1);
               }
            }
         } else {
            // An atom with minimal description. Let's make up a
            // temperature factor:
            at->tempFactor = 1.0;
            at->WhatIsSet = at->WhatIsSet | 4; // is isotropic
         } 
//          std::cout << "on setting WhatIsSet is "
//                    << at->WhatIsSet << "\n";
         if (have_udd_atoms) { 
            at->PutUDData(udd_afix_handle_in, current_afix);
         }
      
      } else {
         std::cout << "(make_atom) bad atom: " << card.card << std::endl;
         delete at;
         at = 0;
      }
   }
   return at;
}

void
coot::ShelxIns::save_fvars(const shelx_card_info_t &card) {

   for (unsigned int i=1; i<card.words.size(); i++) {
//       std::cout << "DEBUG:: save_fvars: pushing back fvar " << i << " "
//                 << atof(card.words[i].c_str()) << std::endl;
      fvars.push_back(atof(card.words[i].c_str()));
   }
}


coot::shelx_card_info_t
coot::ShelxIns::read_line(std::ifstream &f) {

   coot::shelx_card_info_t shelx_card_info;
   std::string s;
   std::vector<std::string> sv;

   char c;
   std::string running_word;
   while (!f.eof()) {
      c = f.get();
      int c_as_int = c;
      // BL says:: we shall include tab here too (9)
      if ((c_as_int >= 32) || (c_as_int == 10) || (c_as_int == 9)) { // George Sheldrick suggestion 20070415
         if (c == '\n') {
            if (running_word.length() > 0)
               sv.push_back(running_word);
            break;
         } else {
            // BL says:: let's include tabs as separator of words too
            if (c == ' ' || c == '\t') {
               if (running_word.length() > 0)
                  sv.push_back(running_word);
               running_word = "";
            } else {
               running_word += c;
            }
            s += c;
         }
      }
   }
//    std::cout << "---- FOUND " << sv.size() << " words in " << s << std::endl;
   shelx_card_info.card = s;
   shelx_card_info.words = sv;

   // was that line starting with one or more spaces?  If so, then we
   // add info to the card, saying so, because (in the calling
   // routine) we don't want this line if it was starting with spaces
   // (unless the previous line was a continuation).

   if (s.length() > 0) {
      if (s[0] == ' ') {
         shelx_card_info.spaced_start = 1;
      }
      if (s[0] == '\t') {
         shelx_card_info.spaced_start = 1;
      }
   }

   
//    for (int i=0; i<shelx_card_info.words.size(); i++) {
//       std::cout << ":" << shelx_card_info.words[i] << ": ";
//    }
//   std::cout << std::endl;
   return shelx_card_info; 
}

mmdb::Residue *
coot::ShelxIns::add_shelx_residue(std::vector<mmdb::Atom *> &atom_vector,
                                  const std::string &current_res_name,
                                  int &current_res_no) const {
   
   mmdb::Residue *residue = new mmdb::Residue;
   residue->SetResName(current_res_name.c_str());
   residue->seqNum = current_res_no;
   bool srn = util::is_standard_residue_name(current_res_name);
   for (unsigned int i=0; i<atom_vector.size(); i++) {
      if (! srn)
         atom_vector[i]->Het = 1;
      residue->AddAtom(atom_vector[i]);
   } 
   return residue;
}


coot::shelx_card_info_t
coot::ShelxIns::read_card(std::ifstream &f) {

   // OK, we have to deal with a commented line (using !) (that it
   // ends with with a "=" doesn't mean that it is continued).
   //
   // How shall we do that?  shelx_card_info returns raw strings, so
   // here we shall carve off the commented lines

   coot::shelx_card_info_t shelx_card_info = read_line(f);
   shelx_card_info.strip_post_bang();
   if (shelx_card_info.words.size() > 0) {
      if (shelx_card_info.spaced_start == 0) {
         if (shelx_card_info.last_word_is_equal_symbol()) {
            shelx_card_info.add_card(read_card_extended(f));
//             std::cout << "DEBUG:: extending" << shelx_card_info.card << "\n";
//             for (int i=0; i<shelx_card_info.words.size(); i++)
//                std::cout << "DEBUG:: extending... " << i << " :" << shelx_card_info.words[i]
//                          << ":\n";
         }
      }
   }

   // we want to return empty if spaced_start was true, unless the
   // previous line ended with a "=".
   if (shelx_card_info.spaced_start) {
      shelx_card_info.empty_yorself();
   } 
   return shelx_card_info;
}



// Read a extension card, a card like
// "    1.1 2.2 3.3"
// that follows a card ending in a "="
//
coot::shelx_card_info_t
coot::ShelxIns::read_card_extended(std::ifstream &f) {

   coot::shelx_card_info_t shelx_card_info = read_line(f);
   shelx_card_info.strip_post_bang();
   if (shelx_card_info.words.size() > 0) {
      if (shelx_card_info.spaced_start == 0) {
         //          if (shelx_card_info.words.back() == "=") {
         if (shelx_card_info.last_word_is_equal_symbol()) {
            shelx_card_info.add_card(read_card_extended(f));
            // std::cout << "DEBUG::" << shelx_card_info.card << "\n";
            // for (int i=0; i<shelx_card_info.words.size(); i++)
            // std::cout << "DEBUG:: " << i << " :" << shelx_card_info.words[i] << ":\n";
         }
      }
   }

   return shelx_card_info;
}


bool
coot::shelx_card_info_t::last_word_is_equal_symbol() const {

   bool r = 0;
   if (words.size() > 0) { 
      std::string s = words.back();
      if (s == "=") {
         r = 1;
      } else {
         if (s.length() == 2) {
            if (s[0] == 61) { // = symbol
               if (s[1] == 13) {
                  std::cout << "windows =" << std::endl;
                  r = 1;
               }
            }
         }
      }
   }
   return r;

}


// modify yourself so that words after a ! word (and the ! word) are
// removed.
void
coot::shelx_card_info_t::strip_post_bang() {

   if (bang_index() != -1) {
      std::vector<std::string> new_words;
      for (unsigned int i=0; i<words.size(); i++) { 
         if (words[i][0] == '!') {
            words = new_words;
            break;
         } else {
            new_words.push_back(words[i]);
         }
      }
   }
}


// return -1 on no bang, else return the index of the ! (bang).
int
coot::shelx_card_info_t::bang_index() const {

   int index = -1; // no bang
   for (unsigned int i=0; i<words.size(); i++) { 
      if (words[i][0] == '!') {
         index = i;
         break;
      }
   }
   return index; 
}

// This is not used, I think.
std::string
coot::ShelxIns::make_atom_name(const std::string &atom_name_in,
                               const int &atomic_weight) const {

   std::string s;
   int len = atom_name_in.length();
   if (len == 4)
      s = atom_name_in;
   if (len == 3)
      s = std::string(" ") + atom_name_in;
   if (len == 2)
      s = std::string("  ") + atom_name_in;
   if (len == 1)
      s = std::string("   ") + atom_name_in;

//    std::cout << "debug:: make_atom_name converts :" << atom_name_in
//              << ": to :" << s << ":\n";
   return s;
}

std::string
coot::ShelxIns::make_atom_name(const std::string &atom_name_in,
                               const std::string &element) const {

   std::string s = coot::pad_atom_name(atom_name_in, element);
//    std::cout << "debug:: make_atom_name converts :" << atom_name_in
//               << ": to :" << s << ":\n";
   return s;
}



std::string
coot::ShelxIns::make_atom_element(const std::string &atom_name_in,
                                  const int &sfac_index) const {

   // sfac_index is fortran indexing
   std::string s("ERROR-in-SFAC");

   int vind = sfac_index - 1; 
   if ( vind < int(sfac.size())) {
      if ( vind >= 0) {
         s = sfac[sfac_index-1];
         if (s.length() == 1)
            s = std::string(" ") + s;
      } else {
         std::cout << "ERROR:: Bad vind! " << vind << " sfac index limit: " << sfac.size() << "\n";
         std::cout << "        sfac_index: " << sfac_index << " for atom name :" << atom_name_in
                   << ":" << std::endl;
      }
   } else {
      std::cout << "ERROR:: Bad vind! over end: " << vind << " sfac index limit: " << sfac.size() << "\n";
   }
   return s;
}


bool
coot::ShelxIns::try_assign_cell(mmdb::Manager *mol) {

   if (!have_cell_flag) {
      mmdb::mat44 my_matt;
      int err = mol->GetTMatrix(my_matt, 0, 0, 0, 0);
      if (err != 0) {
         std::cout << "!! Warning:: No symmetry available for this template molecule"
                   << std::endl;
      } else {
         mmdb::realtype a[6];
         mmdb::realtype vol;
         int orthcode;
         mol->GetCell(a[0], a[1], a[2], a[3], a[4], a[5], vol, orthcode);
         
         clipper::Cell_descr cdr(a[0], a[1], a[2],
                                 clipper::Util::d2rad(a[3]),
                                 clipper::Util::d2rad(a[4]),
                                 clipper::Util::d2rad(a[5]));
         cell = clipper::Cell(cdr);
//          std::cout << "try_assign_cell assigned cell:"
//                    << cell.format() << std::endl;
         have_cell_flag = 1;
      }
   }
   return have_cell_flag;
}

// mol_is_from_shelx_ins is default arg with default value true
//
// This is the API from code using this class (e.g. in molecule-class-info)
//
std::pair<int, std::string>
coot::ShelxIns::write_ins_file(mmdb::Manager *mol_in,
                               const std::string &filename,
                               bool mol_is_from_shelx_ins) {

   std::pair<int, std::string> r(-1,"");

   if (!have_cell_flag) { // Need cell for orth->frac convertion for atoms
      have_cell_flag = try_assign_cell(mol_in);
   }
   if (mol_is_from_shelx_ins) {
      bool mol_needs_reshelx = mol_needs_shelx_transfer(mol_in);
      if (mol_needs_reshelx) {
         mmdb::Manager *mol = reshelx(mol_in);
         r = write_ins_file_internal(mol, filename, true);
         delete mol;
      } else {
         r = write_ins_file_internal(mol_in, filename, true);
      }
   } else {
      r = write_ins_file_internal(mol_in, filename, false);
   }
   return r;
}

bool
coot::ShelxIns::mol_needs_shelx_transfer(mmdb::Manager *mol) const {

   // the test is the same for unshelx or re-shelx

   bool needs_unshelx = true;
   if (! mol) {
      std::cout << "   ERROR:: mol_needs_shelx_transfer() was passed a null mol " << std::endl;
   } else {
      int n_models = mol->GetNumberOfModels();
      // std::cout << "   debug:: mol_needs_shelx_transfer() mol has " << n_models
      // << " models " << std::endl;
      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      if (! model_p) {
         std::cout << "   ERROR:: shelx read_file() No model for 1 " << std::endl;
      } else {
         int n_chains = model_p->GetNumberOfChains();
         if (n_chains > 1) {
            needs_unshelx = false;
         } else {
            for (int ichain=0; ichain<n_chains; ichain++) {
               mmdb::Chain *chain_p = model_p->GetChain(ichain);
               std::string chain_id = chain_p->GetChainID();
               if (false)
                  std::cout << "   DEBUG:: in mol_needs_shelx_transfer() chain_p "
                            << chain_p << " has " << chain_p->GetNumberOfResidues()
                            << " and chain id :" << chain_p->GetChainID() << ":" << std::endl;
               if (chain_id.length() > 0)
                  needs_unshelx = false;
            }
         }
      }
   }
   // std::cout << "   DEBUG:: mol_needs_shelx_transfer() returns " << needs_unshelx << std::endl;
   return needs_unshelx;
}

void
coot::ShelxIns::write_orthodox_pre_atom_lines(std::ofstream &f) const {

   bool sfac_done = false;
   for (unsigned int i=0; i<pre_atom_lines.size(); i++) {
      if (is_unit_line(pre_atom_lines[i])) {
         write_sfac_line(f);
         write_disp_lines(f);
         sfac_done = true;
         f << pre_atom_lines[i];
         if (sfac.size() >= unit.size()) {
            // std::cout << "INFO :: padding UNIT card from size "
            // << unit.size() << " to " << sfac.size() << std::endl;
            for (unsigned int iextra=0; iextra<(sfac.size() - unit.size()); iextra++) {
               f << " 0"; // pad with 0s, like GMS suggests, 20080601
            }
         }
         f << "\n";
      } else {
         f << pre_atom_lines[i];
         f << "\n";
      } 
   }
         
   // SFAC line
   if (! sfac_done)
      write_sfac_line(f);

   // FVAR lines
   int fvar_count = 0;
   for (unsigned int i=0; i<fvars.size(); i++) {
      if (fvar_count == 0)
         f << "FVAR ";
      f.precision(5); // Notes from Doug Kuntz, need 5 not 8, else
      // shelxl bombs.
      f << "  " << fvars[i];
      fvar_count++;
      if (fvar_count == 7) { 
         f << "\n";
         fvar_count = 0;
      }
   }
   f << "\n\n";
}


// not const because we add atoms to sfac.
void
coot::ShelxIns::write_synthetic_pre_atom_lines(mmdb::Manager *mol,
                                               std::ofstream &f) {

   // TITL
   // CELL
   // LATT
   // SYMM
   // SFAC ??
   // UNIT

   f << "TITL PDB->ins\n";

   if (have_cell_flag) { 
      mmdb::cpstr sg = mol->GetSpaceGroup();
      if (sg) {
         std::pair<clipper::Cell, clipper::Spacegroup> cell_sg = util::get_cell_symm(mol);
         f << "CELL 1.54178  ";
         f << std::right << std::setprecision(4) << std::fixed
           << cell.descr().a() << " "
           << cell.descr().b() << " "
           << cell.descr().c() << " " 
           << clipper::Util::rad2d(cell.descr().alpha()) << " " 
           << clipper::Util::rad2d(cell.descr().beta())  << " " 
           << clipper::Util::rad2d(cell.descr().gamma()) << "\n";
         int n_symm = cell_sg.second.num_symops();
         f << "ZERR " << n_symm << "         "
           << 0.001 * cell.descr().a() << "  " 
           << 0.001 * cell.descr().b() << "  " 
           << 0.001 * cell.descr().c() << "  ";
      
         // This (ZERR for angles) is a terrible hack - we need to
         // determine the variable angles from the space group.
         // 
         double zerr_angle[3] = {0.05, 0.05, 0.05};
         if (util::close_double_p(clipper::Util::rad2d(cell.descr().alpha()),  90.0))
            zerr_angle[0] = 0.0;
         if (util::close_double_p(clipper::Util::rad2d(cell.descr().alpha()), 120.0))
            zerr_angle[0] = 0.0;
         if (util::close_double_p(clipper::Util::rad2d(cell.descr().beta()),   90.0))
            zerr_angle[1] = 0.0;
         if (util::close_double_p(clipper::Util::rad2d(cell.descr().beta()),  120.0))
            zerr_angle[1] = 0.0;
         if (util::close_double_p(clipper::Util::rad2d(cell.descr().gamma()),  90.0))
            zerr_angle[2] = 0.0;
         if (util::close_double_p(clipper::Util::rad2d(cell.descr().gamma()), 120.0))
            zerr_angle[2] = 0.0;
         f << zerr_angle[0] << "  "
           << zerr_angle[1] << "  "
           << zerr_angle[2] << "\n";

         std::string spg(sg);
         if (spg.length() > 1) {
            int latt = 1;
            char lat_char = spg[0];
            if (lat_char == 'P')
               latt = 1;
            if (lat_char == 'I')
               latt = 2;
            if (lat_char == 'R')
               latt = 3;
            if (lat_char == 'F')
               latt = 4;
            if (lat_char == 'A')
               latt = 5;
            if (lat_char == 'B')
               latt = 6;
            if (lat_char == 'C')
               latt = 7;

            if (! cell_sg.second.is_null()) {
               // P31 (and others presumably) has 1 inversion symop, it seems.
               //
               // make the latt negative if non-centrosymmetric.
               if (cell_sg.second.num_inversion_symops() <= 1)
                  latt = -latt;
            }
            f << "LATT " << latt << "\n";
            for (int isym=1; isym<cell_sg.second.num_primitive_symops(); isym++) { 
               f << "SYMM " << util::Upper(cell_sg.second.primitive_symop(isym).format()) << "\n";
            }
            f << "\n";
         }

         // SFAC & UNIT
         std::map<std::string, unsigned int> atomic_contents = get_atomic_contents(mol);
         std::map<std::string, unsigned int>::const_iterator it;
         if (atomic_contents.size()) {
            f << "SFAC ";
            for (it=atomic_contents.begin(); it!=atomic_contents.end(); it++)
               f << " " << it->first << " ";
            f << "\n";
            f << "UNIT " ;
            for (it=atomic_contents.begin(); it!=atomic_contents.end(); it++)
               f << it->second * n_symm << " ";
            f << "\n";
            for (it=atomic_contents.begin(); it!=atomic_contents.end(); it++) { 
               add_sfac(it->first);
            }
         }
      }
   }

   // f << "DEFS 0.02 0.1 0.01 0.04\n";
   f << "CGLS 30 -1\n"; // 30 cycles, use-R-free
   f << "SHEL 10 0.1\n";
   f << "FMAP 2\n";
   f << "PLAN 200 2.3\n";
   f << "LIST 6\n";
   f << "WPDB 2\n";

   bool is_aniso = mol_is_anisotropic(mol);

   // f << "ISOR 0.1 O_3001 > LAST\n";
   // f << "CONN 0 O_3001 > LAST\n";

   // --- ISOR and CONN ------
   //
   std::vector<hetatom_range> hetatom_ranges;

   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   mmdb::Chain *chain_p;
   int n_chains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<n_chains; ichain++) {
      int resno_offset = 0;
      resno_offset = ichain * 1000;
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      mmdb::Residue *residue_p = 0;
      mmdb::Atom *at = 0;
      mmdb::Atom *running_hetatm = 0;
      hetatom_range current_range;
      for (int ires=0; ires<nres; ires++) { 
         residue_p = chain_p->GetResidue(ires);
         int n_atoms = residue_p->GetNumberOfAtoms();
         for (int iat=0; iat<n_atoms; iat++) {
            at = residue_p->GetAtom(iat);
            if (! at->isTer()) {
               if (at->Het) {
                  running_hetatm = at;
                  if (current_range.range_first == 0) {
                     current_range.range_first = at;
                     current_range.resno_offset = resno_offset;
                  }
               } else {
                  if (current_range.range_first) {
                     // Bizarro-world
                     current_range.range_last = running_hetatm;
                     hetatom_ranges.push_back(current_range);
                  }
               } 
            }
         }
      }
      if (residue_p) {
         if (running_hetatm) {
            // this is the last residue in the chain, so set last value
            // if it has not been set.
            if (! current_range.range_last) {
               current_range.range_last = running_hetatm;
               hetatom_ranges.push_back(current_range);
            }
         }
      }
   }

   if (is_aniso) {
      for (unsigned int ir=0; ir<hetatom_ranges.size(); ir++) {
         const hetatom_range &hr = hetatom_ranges[ir];
         f << "ISOR 0.1 "
           << util::remove_leading_spaces(hr.range_first->element) << "_"
           << hr.range_first->GetSeqNum() + hr.resno_offset << " > " 
           << util::remove_leading_spaces(hr.range_last->element) << "_" 
           << hr.range_last->GetSeqNum() + hr.resno_offset << "\n";
      }
   }
   for (unsigned int ir=0; ir<hetatom_ranges.size(); ir++) {
      const hetatom_range &hr = hetatom_ranges[ir];
      f << "CONN 0 "
        << util::remove_leading_spaces(hr.range_first->element) << "_"
        << hr.range_first->GetSeqNum() + hr.resno_offset << " > " 
        << util::remove_leading_spaces(hr.range_last->element) << "_" 
        << hr.range_last->GetSeqNum() + hr.resno_offset << "\n";
        // << " > LAST\n";
   }
   if (is_aniso) {
      std::map<std::string, unsigned int> atomic_contents = get_atomic_contents(mol);
      std::map<std::string, unsigned int>::const_iterator it;
      if (atomic_contents.size()) {
         f << "DELU ";
         for (it=atomic_contents.begin(); it!=atomic_contents.end(); it++)
            f << " %" << it->first << "_* ";
         f << "\n";
      }
      f << "RIGU\n";
   }
   
   f << "BUMP\n";

} 


// for SFAC & UNIT
//
std::map<std::string, unsigned int>
coot::ShelxIns::get_atomic_contents(mmdb::Manager *mol) const {

   std::map<std::string, unsigned int> m;

   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   mmdb::Chain *chain_p;
   int n_chains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<n_chains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      mmdb::Residue *residue_p;
      mmdb::Atom *at;
      for (int ires=0; ires<nres; ires++) { 
         residue_p = chain_p->GetResidue(ires);
         int n_atoms = residue_p->GetNumberOfAtoms();
         for (int iat=0; iat<n_atoms; iat++) {
            at = residue_p->GetAtom(iat);
            if (! at->isTer()) { 
               std::string ele(at->element);
               if (! ele.empty())  // now we added the isTer() test this probably
                                   //  won't catch anything.
                  m[ele]++; // initial/default value is 0.
            }
         }
      }
   }
   return m;
} 


// write_synthetic_pre_atom_lines changes sfac so this can't be const (boo).
// 
// mol_is_from_shelx_ins is a default argument (false).
std::pair<int, std::string>
coot::ShelxIns::write_ins_file_internal(mmdb::Manager *mol,
                                        const std::string &filename,
                                        bool mol_is_from_shelx_ins) {

   int istat = 0;
   std::string message;

   int udd_riding_atom_negative_u_value_handle_local = mol->GetUDDHandle(mmdb::UDR_ATOM, "riding_atom_negative_u");

   float u_to_b = 8.0 * M_PI * M_PI;  // perhaps this should be a function

   if (have_cell_flag) { // Need cell for orth->frac convertion for atoms
      
      std::ofstream f(filename.c_str());
      double a = cell.descr().a();
      double b = cell.descr().b();
      double c = cell.descr().c();
      if (f) {

         // pre-atom lines.
         //
         // THE SFAC line has to be inserted into the pre-atom lines -
         // before the UNIT card (AFAICS - maybe after SYMM).  The
         // UNIT line need to be padded with 0s to the length of the
         // SFAC line.
         //

         if (pre_atom_lines.size() == 0) {

            // we are trying to write a shelx ins file from a PDB file
            // 
            write_synthetic_pre_atom_lines(mol, f);

         } else {
            write_orthodox_pre_atom_lines(f);
         }
      
         std::string current_altloc = "";
         int n_models = mol->GetNumberOfModels();
         int current_afix = -1; // AFIX is independent of residue, it
                                // should not be reset for every
                                // residue, put it here (outside the
                                // residue atoms loop).
         std::vector<std::string> afix_failure_vec;
         for (int imod=1; imod<=n_models; imod++) { 
      
            mmdb::Model *model_p = mol->GetModel(imod);
            mmdb::Chain *chain_p;
            // run over chains of the existing mol
            int nchains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<nchains; ichain++) {
               int resno_offset = 0;
               if (! mol_is_from_shelx_ins)
                  resno_offset = ichain * 1000;  // fake an offset
               chain_p = model_p->GetChain(ichain);
               if (chain_p == NULL) {  
                  // This should not be necessary. It seem to be a
                  // result of mmdb corruption elsewhere - possibly
                  // DeleteChain in update_molecule_to().
                  std::cout << "NULL chain in ... " << std::endl;
               } else { 
                  int nres = chain_p->GetNumberOfResidues();
                  mmdb::PResidue residue_p;
                  mmdb::Atom *at;
                  for (int ires=0; ires<nres; ires++) { 
                     residue_p = chain_p->GetResidue(ires);
                     
                     int n_atoms = residue_p->GetNumberOfAtoms();
                     if (n_atoms > 0) {
                        try { 
                           f <<  "RESI ";
                           // format to look like shelx resi:
                           bool output_old_style = false;
                           int resno = residue_p->GetSeqNum();
                           int resno_out = resno + resno_offset;
                           if (output_old_style) {
                              if (resno_out < 1000)
                                 f << " ";
                              if (resno_out < 100)
                                 f << " ";
                              if (resno_out < 10)
                                 f << " ";
                              f << resno_out << "   " << residue_p->GetResName() << "\n";
                           } else {
                              if (resno_out < 1000)
                                 f << " ";
                              if (resno_out < 100)
                                 f << " ";
                              if (resno_out < 10)
                                 f << " ";
                              std::string chain_id = residue_p->GetChainID();
                              if (! chain_id.empty()) {
                                 f << chain_id <<  ":";
                              } else {
                                 std::cout << "in write_ins_file_internal() residue " << residue_spec_t(residue_p)
                                           << " had empty chain id " << std::endl;
                              }
                              f << resno_out << " " << residue_p->GetResName() << "\n";
                           }
                        }
                        catch (const std::ios::failure &e) { 
                           std::cout << "WARNING:: IOS exception caught on RESI start "
                                     << e.what() << std::endl;
                        }
                     }
                        
                     int ic;
                     
                     for (int iat=0; iat<n_atoms; iat++) {
                        try { 
                           at = residue_p->GetAtom(iat);
                           float site_occ_factor = at->occupancy;
                           // reset occupancy to shelx standard occ
                           // (FVAR 1 is implied if not set in the .ins file)
                           if (! mol_is_from_shelx_ins)
                              site_occ_factor = 11.000;

                           clipper::Coord_orth co(at->x, at->y, at->z);
                           clipper::Coord_frac cf = co.coord_frac(cell);
                           int sfac_index = get_sfac_index(at->element);

                           // AFIX comes before PART
                           if (at->GetUDData(udd_afix_handle, ic) == mmdb::UDDATA_Ok) {
                              if (ic != current_afix) { 
                                 f << "AFIX ";
                                 // a bit of formatting to match SHELX format for AFIX
                                 if (ic < 100)
                                    f << " ";
                                 if (ic < 10)
                                    f << " ";
                                 f << ic << "\n";
                              }
                              current_afix = ic;
                           } else {
                              std::string s = message_for_atom("WARNING:: failed to get AFIX handle for ", at);
                              afix_failure_vec.push_back(s);
                           }
                        
                           // std::cout << "on writting WhatIsSet is " << at->WhatIsSet
                           // << "\n";
                        
                           std::string this_altloc = at->altLoc;
                           if (this_altloc != current_altloc) {
                              int ipart = altloc_to_part_no(std::string(at->altLoc));
                              f << "PART    " << ipart << "\n";
                           }
                           if (at->WhatIsSet & mmdb::ASET_Anis_tFac) {
                              // Anisotropic

                              clipper::U_aniso_orth cao(at->u11, at->u22, at->u33,
                                                        at->u12, at->u13, at->u23);
                              clipper::U_aniso_frac caf = cao.u_aniso_frac(cell);

                              std::string at_name(at->name);
                              f.setf(std::ios::fixed);
                              f.precision(9);
                              f << coot::util::remove_leading_spaces(at_name)
                                << "   " << sfac_index << "  ";
                              f << cf.u() << "  " << cf.v() << "   " << cf.w() << " "
                                << site_occ_factor << "    =\n     ";
                              f.precision(5);
                              // f << at->u11 << "  " << at->u22 << "  " << at->u33 << "  "
                              //   << at->u23 << "  " << at->u13 << "  " << at->u12
                              //   << "\n";
                              f << caf(0,0)*a*a << "  " << caf(1,1)*b*b << "  " << caf(2,2)*c*c
                                << "  "
                                << caf(1,2)*b*c << "  " << caf(0,2)*a*c << "  " << caf(0,1)*a*b
                                << "\n";
                           } else { 
                              if (at->WhatIsSet & mmdb::ASET_tempFactor) {
                                 // Isotropic B factor
                                 std::string at_name(at->name);
                                 float b_factor = at->tempFactor;
                                 f.setf(std::ios::fixed);
                                 f.precision(7);

                                 mmdb::realtype negative_u;
                                 int status_2 = at->GetUDData(udd_riding_atom_negative_u_value_handle_local,
                                                              negative_u);

                                 if (status_2 == mmdb::UDDATA_Ok) {
                                    f << coot::util::remove_leading_spaces(at_name)
                                      << "   " << sfac_index << "  "
                                      << cf.u() << "  " << cf.v() << "   " << cf.w() << " "
                                      << site_occ_factor << "    "
                                      << negative_u << "\n";
                                 } else {
                                    f.setf(std::ios::fixed);
                                    f.precision(7);
                                    if (b_factor > 0.0) // (riding?) hydrogens are not (-1.2)
                                       b_factor /= u_to_b;
                                    f << coot::util::remove_leading_spaces(at_name)
                                      << "   " << sfac_index << "  "
                                      << cf.u() << "  " << cf.v() << "   " << cf.w() << " "
                                      << site_occ_factor << "    "
                                      << b_factor << "\n";
                                 }
                              }
                           }
                           current_altloc = this_altloc;
                        }
                        catch (const std::ios::failure &e) { 
                           std::cout << "WARNING:: IOS exception caught: " << e.what() << std::endl;
                        }
                     }
                     try { 
                        f << " \n"; // end of a RESI
                     } 
                     catch (const std::ios::failure &e) { 
                        std::cout << "WARNING:: IOS exception caught on end of a RESI " << e.what() << std::endl;
                     }
                  }
               }
            }
         }

         // print out 10 first AFIX misses
         int n_afix_failures = afix_failure_vec.size();
         if (n_afix_failures > 10)
            n_afix_failures = 10;
         for (int i=0; i<n_afix_failures; i++)
            std::cout << afix_failure_vec[i] << std::endl;
         if (n_afix_failures > 10)
            std::cout << "WARNING:: and " << afix_failure_vec.size() - 10
                      << " more AFIX failures" << std::endl;


         if (post_atom_lines.size() > 0) { 
            for (unsigned int i=0; i<post_atom_lines.size(); i++) {
               try { 
                  f << post_atom_lines[i] << "\n";
               }
               catch (const std::ios::failure &e) { 
                  std::cout << "WARNING:: IOS exception caught in post atom lines " << e.what() << std::endl;
               }
            }
         } else {

            // Synthetic post-atom lines
            f << "HKLF 3\n";
            f << "END\n";
            
         } 
      }
      f.close();
      message =  "INFO:: SHELXL file ";
      message += filename;
      message += " written.";
      istat = 1;
      
   } else {
      std::cout << "WARNING:: no cell available... failure to write ins file."
                << std::endl;
      message = "WARNING:: no cell available... failure to write ins file.";
   }
   return std::pair<int, std::string>(istat, message);
}

std::string
coot::ShelxIns::message_for_atom(const std::string &in_string, mmdb::Atom *at) const {

   std::string s = in_string;
   s += "\""; 
   s += at->GetChainID();
   s += "\""; 
   s += " "; 
   s += coot::util::int_to_string(at->GetSeqNum());
   s += " ";
   s += "\""; 
   s += at->GetResName();
   s += "\""; 
   s += " ";
   s += "\""; 
   s += at->GetAtomName();
   s += "\"";
   if (std::string(at->altLoc).length()) {
      s += " ,";
      s += "\""; 
      s += at->altLoc;
      s += "\"";
   }
   return s;
}


void
coot::ShelxIns::write_sfac_line(std::ostream &f) const {

   if (sfac.size() > 0) { 
      f << "SFAC";
      for (unsigned int i=0; i<sfac.size(); i++) {
         f << "  " << sfac[i];
      }
      f << "\n";
   }
}

void
coot::ShelxIns::write_disp_lines(std::ostream &f) const {
   for (std::size_t i=0; i<disp_cards.size(); i++)
      f << disp_cards[i] << "\n";
}


// return -1 on index not found
int
coot::ShelxIns::get_sfac_index(const std::string &element) const {

   std::string ele = element;
   if (ele[0] == ' ') {
      ele = element.substr(1,1);
   }

   int indx = -1;
   // std::cout << "in get_sfac_index() for :" << ele << ": sfac.size() is " << sfac.size() << std::endl;
   for (unsigned int i=0; i<sfac.size(); i++) {
      // std::cout << "comparing :" << ele << ": :" << sfac[i] << ":\n";
      if (ele == sfac[i]) {
         indx = i+1;
         break;
      }
   }
   return indx;
}

// do redundancy checking.
void
coot::ShelxIns::add_sfac(const std::string &ele_in) {

   std::string ele=util::remove_leading_spaces(ele_in);

   bool found = false;
   for (unsigned int i=0; i<sfac.size(); i++) {
      if (sfac[i] == ele) {
         found = true;
         break;
      }
   }

   if (! found) {
      sfac.push_back(ele);
   }
}



int
coot::ShelxIns::altloc_to_part_no(const std::string &altloc) const {

   int ipart = 0;
   if (altloc == "")
      return 0;
   if (altloc == "A")
      return 1;
   if (altloc == "B")
      return 2;
   if (altloc == "C")
      return 3;
   if (altloc == "D")
      return 4;
   if (altloc == "E")
      return 5;
   if (altloc == "F")
      return 6;
   if (altloc == "G")
      return 7;
   if (altloc == "H")
      return 8;
   if (altloc == "a")
      return -1;
   if (altloc == "b")
      return -2;
   if (altloc == "c")
      return -3;
   if (altloc == "d")
      return -4;
   if (altloc == "e")
      return -5;
   if (altloc == "f")
      return -6;
   if (altloc == "g")
      return -7;
   if (altloc == "h")
      return -8;
   return ipart;

}

int
coot::ShelxIns::add_fvar(float f) {

   // FVAR 1 is not written to SHELX file
   //
   // 20060315 BUT! the overall scale (FO to FC?) is written in the
   // .res file as if it was free variable number one.

   // so if shelx ins file seems to have 1
   // FVAR value, then we've just created shelx FVAR 2.

   fvars.push_back(f);
   return fvars.size();  // fvars is increased in size now, of course.
}

void
coot::ShelxIns::set_fvar(int fvar_no, float f) {

   if (int(fvars.size()) >= (fvar_no-1))
      if (fvar_no >= 0)
         fvars[fvar_no - 1] = f;  // i.e. free varable no 2 is at index 1
}


void
coot::ShelxIns::debug() const {

   std::cout << "DEBUG ShelxIns title: " << title << std::endl;
   std::cout << "DEBUG ShelxIns filled_flag: " << filled_flag << std::endl;
   std::cout << "DEBUG ShelxIns : have_cell_flag: " << have_cell_flag<< std::endl;

   std::cout << "DEBUG ShelxIns : cell " << cell.format() << std::endl;
   std::cout << "DEBUG ShelxIns : sfac size " << sfac.size() << std::endl;
   std::cout << "DEBUG ShelxIns : unit size " << unit.size() << std::endl;
   std::cout << "DEBUG ShelxIns : defs size " <<  defs.size() << std::endl;
   std::cout << "DEBUG ShelxIns : fvars size " <<  fvars.size() << std::endl;
   std::cout << "DEBUG ShelxIns : pre_atom_lines.size() " <<  pre_atom_lines.size() << std::endl;
   std::cout << "DEBUG ShelxIns : post_atom_lines.size() " <<  post_atom_lines.size() << std::endl;
   
}

// Do lattice expansion and possible centro-symmetric expansion
// and add in X,Y,Z which shelx ins file does not require.
// 
std::vector<std::string>
coot::clipper_symm_strings(const std::vector<std::string> &symm_vec,
                           int shelx_latt) {

   std::vector<std::string> v;
   v.push_back("X,Y,Z");
   std::vector<std::string> r;
   
   for (unsigned int i=0; i<symm_vec.size(); i++)
      v.push_back(symm_vec[i]);

//     for (unsigned int i=0; i<v.size(); i++)
//        std::cout << "DEBUG:: v in clipper_symm_strings: "
//                  << i << " " << v[i] << std::endl;

   for (unsigned int i=0; i<v.size(); i++) { 
      symm_card_composition_t sc(v[i]);
      // std::cout << "symm_card on :" << v[i] << ": gives " << sc.symm_card() << std::endl;
      std::vector<std::string> cards = sc.symm_cards_from_lat(shelx_latt);
//       std::cout << "INFO:: There are " << cards.size()
//                 << " elements in cards" << std::endl;
      for (unsigned int isc=0; isc<cards.size(); isc++) {
         r.push_back(cards[isc]);
      }
   }
   return r;
} 



// constructor
// 
// should be able to cope with 0.25-x, 0.333+y, 0.1667+z
//
coot::symm_card_composition_t::symm_card_composition_t(const std::string &symm_card) {

   // So first we split into 3 strings on comma.
   //
   // That gives us elements like: "x+2/3" say, or "-x+y" or "x-0.3333"
   //
   // Next we need to split the x/y/z component from the numerical
   // component
   //
   // -> [x],2/3   or   [-x,y],0  or  [x],-0.333
   //
   // How do we do that? by stripping out
   // "+x" "-x" then "x" ... "+y" "-y" then "y" ... "+z" "-z" then "z".
   // If we find those components then flag in the x_element array.
   // What we have left is the numerical component.
   //
   // Then we interpret the numerical component.
   // Does it contain a slash?  If yes, just leave it.
   // If no, then turn it into a number and multiply it by 12.
   // 
   // Then see if it is 2, 3, 4, 6 or - those numbers,
   // if it isn't then throw hands up in air and grown loadly.
   // (These days, I guess, an exception should be thrown).

//    std::cout << "\nconstructing symm_card_composition_t from "
//              << symm_card << std::endl << std::endl;

   for (int i=0; i<3; i++) {
      x_element[i] = 0;
      y_element[i] = 0;
      z_element[i] = 0;
      frac_trans[i] = 0;
   }

   // now split on comma
   std::vector<std::string> ele(3, "");

   // Find the first comma:
   std::string::size_type icomma_1 = symm_card.find(",");
   if (icomma_1 != std::string::npos) {
      std::string a = symm_card.substr(0,icomma_1);
      std::string b = symm_card.substr(icomma_1+1);
      std::string::size_type icomma_2 = b.find(",");
      // Find the second comma
      if (icomma_2 != std::string::npos) { 
         std::string c = b.substr(0, icomma_2);
         std::string d = b.substr(icomma_2+1);
         ele[0] = a;
         ele[1] = c;
         ele[2] = d;
      }
   }

   // we have split on comma into ele.  now examine the ele for x,y,z
   // and translational part.
   // 
   for (int iele=0; iele<3; iele++) {
      // these need case insensitivity

      std::string this_ele = coot::util::upcase(ele[iele]);
      
      std::string::size_type ixp = this_ele.find("+X"); 
      std::string::size_type ixm = this_ele.find("-X");
      std::string::size_type ix  = this_ele.find( "X");
      std::string::size_type iyp = this_ele.find("+Y");
      std::string::size_type iym = this_ele.find("-Y");
      std::string::size_type iy  = this_ele.find( "Y");
      std::string::size_type izp = this_ele.find("+Z");
      std::string::size_type izm = this_ele.find("-Z");
      std::string::size_type iz  = this_ele.find( "Z");

      if (ixp == std::string::npos)
         ixp = this_ele.find("+ X"); 
      if (ixm == std::string::npos)
         ixm = this_ele.find("- X"); 

      if (iyp == std::string::npos)
         iyp = this_ele.find("+ Y"); 
      if (iym == std::string::npos)
         iym = this_ele.find("- Y"); 

      if (izp == std::string::npos)
         izp = this_ele.find("+ Z"); 
      if (izm == std::string::npos)
         izm = this_ele.find("- Z"); 

//       std::cout << "looking in iele=" << iele << " " << this_ele
//                 << std::endl;
      
      std::string substituted = this_ele;
      if (ixp != std::string::npos) {
         substituted[ixp] = ' ';
         substituted[ixp+1] = ' ';
      }
      if (ixm != std::string::npos) {
         substituted[ixm] = ' ';
         substituted[ixm+1] = ' ';
      }
      if ((ixp == std::string::npos) &&
          (ixm == std::string::npos)) { 
         if (ix != std::string::npos) {
            substituted[ix] = ' ';
         }
      }
      if (ixp != std::string::npos)
         x_element[iele] = 1;
      if (ixm != std::string::npos)
         x_element[iele] = -1;

      if ((ixp == std::string::npos) &&
          (ixm == std::string::npos))
         if (ix != std::string::npos) {
            x_element[iele] = 1;
         }
      
      if (iyp != std::string::npos) {
         substituted[iyp] = ' ';
         substituted[iyp+1] = ' ';
      }
      if (iym != std::string::npos) {
         substituted[iym] = ' ';
         substituted[iym+1] = ' ';
      }
      if ((iyp == std::string::npos) &&
          (iym == std::string::npos)) { 
         if (iy != std::string::npos) {
            substituted[iy] = ' ';
         }
      }
      if (iyp != std::string::npos)
         y_element[iele] = 1;
      if (iym != std::string::npos)
         y_element[iele] = -1;
      if ((iyp == std::string::npos) &&
          (iym == std::string::npos))
         if (iy != std::string::npos)
            y_element[iele] = 1;

      if (izp != std::string::npos) {
         substituted[izp] = ' ';
         substituted[izp+1] = ' ';
      }
      if (izm != std::string::npos) {
         substituted[izm] = ' ';
         substituted[izm+1] = ' ';
      }
      if ((izp == std::string::npos) &&
          (izm == std::string::npos)) { 
         if (iz != std::string::npos) {
            substituted[iz] = ' ';
         }
      }
      if (izp != std::string::npos)
         z_element[iele] = 1;
      if (izm != std::string::npos)
         z_element[iele] = -1;
      if ((izp == std::string::npos) &&
          (izm == std::string::npos))
         if (iz != std::string::npos)
            z_element[iele] = 1;

      // debugging stuff:
      if (0)
         std::cout << "  after setting: x/y/z eles: (["
                   << x_element[0] << " "
                   << x_element[1] << " "
                   << x_element[2] << "] "
                   << ")  (["
                   << y_element[0] << " "
                   << y_element[1] << " "
                   << y_element[2] << " ] "
                   << ")  (["
                   << z_element[0] << " "
                   << z_element[1] << " "
                   << z_element[2] << "] "
                   << ") " << std::endl;

   
      // so now substituted should be stripped of all X, Y, Z
      // commonents and now we need to examine the trans component.
      // 
      // set whether substituted is blank 
      short int substituted_is_blank = 1;
      for (unsigned int il=0; il<substituted.length(); il++) {
         if (substituted[il] != ' ') {
            substituted_is_blank = 0;
            break;
         }
      }

      if (substituted_is_blank) {
         // no trans components.
         frac_trans[iele] = 0;

      } else { 

         std::string::size_type islash = substituted.find("/");
         if (islash != std::string::npos) {
            std::string::size_type ifind;
            ifind = substituted.find("-1/2");
            if (ifind != std::string::npos) { 
               frac_trans[iele] = MINUS_ONE_HALF;
            } else { 
               ifind = substituted.find("1/2");
               if (ifind != std::string::npos)
                  frac_trans[iele] = ONE_HALF;
            }
            ifind = substituted.find("-1/3");
            if (ifind != std::string::npos) { 
               frac_trans[iele] = MINUS_ONE_THIRD;
            } else { 
               ifind = substituted.find("1/3");
               if (ifind != std::string::npos)
                  frac_trans[iele] = ONE_THIRD;
            }
            ifind = substituted.find("-1/4");
            if (ifind != std::string::npos) { 
               frac_trans[iele] = MINUS_ONE_QUARTER;
            } else { 
               ifind = substituted.find("1/4");
               if (ifind != std::string::npos)
                  frac_trans[iele] = ONE_QUARTER;
            }
            ifind = substituted.find("-1/6");
            if (ifind != std::string::npos) { 
               frac_trans[iele] = MINUS_ONE_SIXTH;
            } else { 
               ifind = substituted.find("1/6");
               if (ifind != std::string::npos)
                  frac_trans[iele] = ONE_SIXTH;
            }
            ifind = substituted.find("-2/3");
            if (ifind != std::string::npos) { 
               frac_trans[iele] = MINUS_TWO_THIRDS;
            } else { 
               ifind = substituted.find("2/3");
               if (ifind != std::string::npos)
                  frac_trans[iele] = TWO_THIRDS;
            }
            ifind = substituted.find("-3/4");
            if (ifind != std::string::npos) { 
               frac_trans[iele] = MINUS_THREE_QUARTERS;
            } else { 
               ifind = substituted.find("3/4");
               if (ifind != std::string::npos)
                  frac_trans[iele] = THREE_QUARTERS;
            }
            ifind = substituted.find("-5/6");
            if (ifind != std::string::npos) { 
               frac_trans[iele] = MINUS_FIVE_SIXTHS;
            } else { 
               ifind = substituted.find("5/6");
               if (ifind != std::string::npos)
                  frac_trans[iele] = FIVE_SIXTHS;
            }

         } else {

            // no slash
            float f = atof(substituted.c_str());
            double ff = f*12;
            frac_trans[iele] = int(round(ff));
            if (0) { 
               std::cout << "DEBUG:: no slash substitued string "
                         << substituted << ": -> " << f << std::endl;
               std::cout << "floating trans: " << f << " " << frac_trans[iele]
                         << std::endl;
            }
         }
      }
   }

   // debugging stuff:
   if (0)
      std::cout << "  x/y/z eles: (["
                << x_element[0] << " "
                << x_element[1] << " "
                << x_element[2] << "] "
                << frac_trans[0] << ")  (["
                << y_element[0] << " "
                << y_element[1] << " "
                << y_element[2] << " ] "
                << frac_trans[1] << ")  (["
                << z_element[0] << " "
                << z_element[1] << " "
                << z_element[2] << "] "
                << frac_trans[2] << ") " << std::endl;
} 


// The centring is:
//  5: A centring additional x, 1/2+y, 1/2+z
//  6: B centring additional x+1/2, y, 1/2+z
//  7: C centring additional 1/2+x, 1/2+y, z
//  4: F centring additional x, 1/2+y, 1/2+z
//                           x+1/2, y, 1/2+z
//                           1/2+x, 1/2+y, z
//
//  3: H centring additional x+2/3,y+1/3,z+1/3
//                           x+1/3,y+2/3,z+2/3
//
//  2: I centring additional x+1/2,y+1/2,z+1/2
// 
std::vector<std::string>
coot::symm_card_composition_t::symm_cards_from_lat(int latt) {

   // This function has been unshadowed - perhaps not in the most intelligent way

   std::vector<std::string> r;
   std::vector<std::vector<int> > latt_additions;
   
   int abs_latt = abs(latt);

   std::vector<int> a(3);
   a[0] = NONE;  a[1] = NONE; a[2] = NONE;
   latt_additions.push_back(a);

   if (abs_latt == 2) {
      std::vector<int> aa(3);
      aa[0] = ONE_HALF;  aa[1] = ONE_HALF; aa[2] = ONE_HALF;
      latt_additions.push_back(aa);
   }
   if (abs_latt == 3) { 
      std::vector<int> aa(3);
      aa[0] = TWO_THIRDS;  aa[1] = ONE_THIRD; aa[2] = ONE_THIRD;
      latt_additions.push_back(aa);
      aa[0] = ONE_THIRD;  aa[1] = TWO_THIRDS; aa[2] = TWO_THIRDS;
      latt_additions.push_back(aa);
   }
   if (abs_latt == 4) { 
      std::vector<int> aa(3);
      aa[0] = 0;  aa[1] = ONE_HALF; aa[2] = ONE_HALF;
      latt_additions.push_back(aa);
      aa[0] = 0;  aa[1] = ONE_HALF; aa[2] = ONE_HALF;
      latt_additions.push_back(aa);
      aa[0] = ONE_HALF;  aa[1] = ONE_HALF; aa[2] = NONE;
      latt_additions.push_back(aa);
   }
   if (abs_latt == 5) { 
      std::vector<int> aa(3);
      aa[0] = 0;  aa[1] = ONE_HALF; aa[2] = ONE_HALF;
      latt_additions.push_back(aa);
   }
   if (abs_latt == 6) {
      std::vector<int> aa(3);
      aa[0] = ONE_HALF;  aa[1] = NONE; aa[2] = ONE_HALF;
      latt_additions.push_back(aa);
   }
   if (abs_latt == 7) {
      std::vector<int> aa(3);
      aa[0] = ONE_HALF;  aa[1] = ONE_HALF; aa[2] = NONE;
      latt_additions.push_back(aa);
   }

   for (unsigned int i=0; i<latt_additions.size(); i++) {
      symm_card_composition_t t = *this;
      t.add_centering_frac(latt_additions[i][0],
                           latt_additions[i][1],
                           latt_additions[i][2]);
      r.push_back(t.symm_card());
      if (latt>0) {
         symm_card_composition_t it = *this;
         it.invert();
         r.push_back(it.symm_card());
      }
   }
   return r;
}


void
coot::symm_card_composition_t::invert() {

   frac_trans[0] = -frac_trans[0];
   frac_trans[1] = -frac_trans[1];
   frac_trans[2] = -frac_trans[2];

   for (int i=0; i<3; i++) {
      x_element[i] = -x_element[i];
      y_element[i] = -y_element[i];
      z_element[i] = -z_element[i];
   }
}


void
coot::symm_card_composition_t::add_centering_frac(int ix, int iy, int iz) {

   frac_trans[0] += ix;
   frac_trans[1] += iy;
   frac_trans[2] += iz;

   for (int i=0; i<3; i++) {
      if (frac_trans[i] < -12)
         frac_trans[i] += 12;
      if (frac_trans[i] > 12)
         frac_trans[i] -= 12;
   }
} 


std::string
coot::symm_card_composition_t::symm_card() const {

   std::string r;

   if (x_element[0] == 1)
      r+= "X";
   if (x_element[0] == -1)
      r+= "-X";

   if (y_element[0] == 1)
      r+= "+Y";
   if (y_element[0] == -1)
      r+= "-Y";
   if (z_element[0] == 1)
      r+= "+Z";
   if (z_element[0] == -1)
      r+= "-Z";
   if (frac_trans[0]) {
      r += fract_trans_to_str(frac_trans[0]);
   }
   r += ",";


   if (x_element[1] == 1)
      r+= "+X";
   if (x_element[1] == -1)
      r+= "-X";
   if (y_element[1] == 1)
      r+= "+Y";
   if (y_element[1] == -1)
      r+= "-Y";
   if (z_element[1] == 1)
      r+= "+Z";
   if (z_element[1] == -1)
      r+= "-Z";
   if (frac_trans[1]) {
      r += fract_trans_to_str(frac_trans[1]);
   }
   r += ",";



   if (x_element[2] == 1)
      r+= "+X";
   if (x_element[2] == -1)
      r+= "-X";
   if (y_element[2] == 1)
      r+= "+Y";
   if (y_element[2] == -1)
      r+= "-Y";
   if (z_element[2] == 1)
      r+= "+Z";
   if (z_element[2] == -1)
      r+= "-Z";
   if (frac_trans[2]) {
      r += fract_trans_to_str(frac_trans[2]);
   }
   
   return r;
}


std::string
coot::symm_card_composition_t::fract_trans_to_str(int itrans_frac) const {

   std::string r;
   if (itrans_frac == NONE)
      r = "";
   if (itrans_frac == ONE_HALF)
      r = "+1/2";
   if (itrans_frac == ONE_THIRD)
      r = "+1/3";
   if (itrans_frac == ONE_QUARTER)
      r = "+1/4";
   if (itrans_frac == ONE_SIXTH)
      r = "+1/6";
   if (itrans_frac == TWO_THIRDS)
      r = "+2/3";
   if (itrans_frac == THREE_QUARTERS)
      r = "+3/4";
   if (itrans_frac == FIVE_SIXTHS)
      r = "+5/6";

   if (itrans_frac == MINUS_ONE_HALF)
      r = "-1/2";
   if (itrans_frac == MINUS_ONE_THIRD)
      r = "-1/3";
   if (itrans_frac == MINUS_ONE_QUARTER)
      r = "-1/4";
   if (itrans_frac == MINUS_ONE_SIXTH)
      r = "-1/6";
   if (itrans_frac == MINUS_TWO_THIRDS)
      r = "-2/3";
   if (itrans_frac == MINUS_THREE_QUARTERS)
      r = "-3/4";
   if (itrans_frac == MINUS_FIVE_SIXTHS)
      r = "-5/6";
   

   return r;
}



// Return -1 on a problem
// 
// static
int
coot::ShelxIns::shelx_occ_to_fvar(float shelx_occ) {

   // e.g. return 18 if shelx_occ 181.00

   int r = -1;

   double aocc = fabs(shelx_occ);
   if (aocc > 10.0)
      r = int(aocc/10.0);

   return r;
}


// Convert the single chained mol into a mol with multiple chains.
// 
// return null on no conversion.
mmdb::Manager *
coot::unshelx(mmdb::Manager *shelx_mol) {

   int skip_chain_step = 21;
   mmdb::Manager *mol = 0;

   if (!shelx_mol) {
      std::cout << "ERROR:: Null shelx_mol" << std::endl;
      return mol;
   }
   
   int imod = 1;
   mmdb::Model *shelx_model_p = shelx_mol->GetModel(imod);
   if (! shelx_model_p) {
      std::cout << "ERROR: unshelx() no model 1 in molecule " << std::endl;
      return NULL;
   }
   mmdb::Chain *chain_p = NULL;
   std::string r("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
   int r_index = 0;
   // run over chains of the existing mol
   int nchains = shelx_model_p->GetNumberOfChains();
   if (nchains != 1) {
      std::cout << "Opps.  Don't know what to do. There are "
                << nchains << " chains and there should be just 1 "
                << std::endl;
   } else {
      mol = new mmdb::Manager;
      int udd_afix_handle_shelx = shelx_mol->GetUDDHandle(mmdb::UDR_ATOM, "shelx afix");
      int udd_afix_handle = mol->RegisterUDInteger(mmdb::UDR_ATOM, "shelx afix");
      int udd_riding_atom_negative_u_value_handle_shelx = shelx_mol->GetUDDHandle(mmdb::UDR_ATOM, "riding_atom_negative_u");
      int udd_riding_atom_negative_u_value_handle_local =  mol->RegisterUDInteger(mmdb::UDR_ATOM, "riding_atom_negative_u");
      mmdb::Model *model_p = new mmdb::Model;
      mol->AddModel(model_p);
      mmdb::Chain *shelx_chain_p = shelx_model_p->GetChain(0);
      int nres = shelx_chain_p->GetNumberOfResidues();
      bool need_new_chain = 1;
      int ires_prev = -1000;
      bool made_afix_transfer_message = false;
      for (int ires=0; ires<nres; ires++) {

         mmdb::Residue *shelx_residue_p = shelx_chain_p->GetResidue(ires);
         int resno = shelx_residue_p->GetSeqNum();
         if (resno > (ires_prev + skip_chain_step))
            // this is a new chain
            need_new_chain = 1;
         
         if (need_new_chain) {
            chain_p = new mmdb::Chain;
            std::string new_chain_id = r.substr(r_index, 1);
            r_index++;
            chain_p->SetChainID(new_chain_id.c_str());
            model_p->AddChain(chain_p);
            need_new_chain = 0;
         }

         mmdb::Residue *copy_residue_p = coot::util::deep_copy_this_residue(shelx_residue_p);
         chain_p->AddResidue(copy_residue_p);

         // apply the shelx afix numbers:
         int shelx_natoms;
         mmdb::PAtom *shelx_residue_atoms = NULL;
         shelx_residue_p->GetAtomTable(shelx_residue_atoms, shelx_natoms);

         int copy_natoms;
         mmdb::PAtom *copy_residue_atoms = NULL;
         copy_residue_p->GetAtomTable(copy_residue_atoms, copy_natoms);

         if (shelx_natoms == copy_natoms) { 
            for (int iat=0; iat<copy_natoms; iat++) {
               int afix;
               int istatus = shelx_residue_atoms[iat]->GetUDData(udd_afix_handle_shelx, afix);
               // istatus can fail and that's OK...
               if (istatus == mmdb::UDDATA_Ok) { 
                  copy_residue_atoms[iat]->PutUDData(udd_afix_handle, afix);
               }
               mmdb::realtype negative_u_value;
               int istatus_2 = shelx_residue_atoms[iat]->GetUDData(udd_riding_atom_negative_u_value_handle_shelx,
                                                                   negative_u_value);
               // transfer of negative_u_value often fails (it's missing from original molecule) and that's OK
               if (istatus_2 == mmdb::UDDATA_Ok) { 
                  int istatus_3 = copy_residue_atoms[iat]->PutUDData(udd_riding_atom_negative_u_value_handle_local,
                                                                     negative_u_value);
               }
            }
         } else {
            std::cout << "ERROR transfering afix: bad copy number of atoms "
                      << shelx_natoms << " " << copy_natoms << std::endl;
         }

         ires_prev = shelx_residue_p->GetSeqNum(); // set up for next round
      }

      // Fix the index of the residues
      // run over chains of the existing mol
      int nchains_local = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains_local; ichain++) {
         chain_p = model_p->GetChain(ichain);
         chain_p->TrimResidueTable();
         mmdb::PResidue residue_p;
         for (int ires=0; ires<nres; ires++) { 
            residue_p = chain_p->GetResidue(ires);
            if (residue_p) 
               residue_p->index = ires;
         }
      }
      mol->FinishStructEdit();
      mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
      

      // need to copy over cell and symmetry info:
      mmdb::realtype a[6];
      mmdb::realtype vol;
      int orthcode;
      shelx_mol->GetCell(a[0], a[1], a[2], a[3], a[4], a[5], vol, orthcode);
      mol->SetCell(a[0], a[1], a[2], a[3], a[4], a[5]);
      mmdb::cpstr sg = shelx_mol->GetSpaceGroup();
      // Don't mess around with copying the string here!
      // SetSpaceGroup() copies the string.
       if (sg)                   
          mol->SetSpaceGroup(sg); 
   }
   return mol;
}

// If the residues are simply put in a different chain when they are
// imported (i.e. there is no renumbering of the residues) then we
// don't need to change the residue numbers back when re-export to
// shelx format.  Which means that the ShelxIns ins_info is not used.
// Hmmm... does that work OK?
mmdb::Manager *
coot::reshelx(mmdb::Manager *mol) {

   mmdb::Manager *shelx_mol = new mmdb::Manager;

   int imod = 1;
   mmdb::Model *shelx_model_p = new mmdb::Model;
   shelx_mol->AddModel(shelx_model_p);
   mmdb::Chain *shelx_chain_p = new mmdb::Chain;
   shelx_model_p->AddChain(shelx_chain_p);
   bool made_afix_transfer_message = false;
   
   int udd_afix_handle = mol->GetUDDHandle(mmdb::UDR_ATOM, "shelx afix");
   int udd_riding_atom_negative_u_value_handle = mol->GetUDDHandle(mmdb::UDR_ATOM, "riding_atom_negative_u");

   // local register for shelx_mol
   int udd_riding_atom_negative_u_value_handle_local =
       shelx_mol->RegisterUDInteger(mmdb::UDR_ATOM, "riding_atom_negative_u");
   int udd_afix_handle_local = shelx_mol->RegisterUDInteger(mmdb::UDR_ATOM, "shelx afix");
   
   // run over chains of the existing mol
   mmdb::Model *model_p = mol->GetModel(imod);
   mmdb::Chain *chain_p;
   int nchains_local = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains_local; ichain++) {
      // int residue_offset = ichain*ins_info.new_chain_offset;
      int residue_offset = 0; // no need to mess with the residue numbers, I think.
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      mmdb::PResidue residue_p;
      for (int ires=0; ires<nres; ires++) { 
         residue_p = chain_p->GetResidue(ires);
         mmdb::Residue *copy_residue_p = coot::util::deep_copy_this_residue(residue_p);
         copy_residue_p->seqNum = residue_p->GetSeqNum() + residue_offset;
         shelx_chain_p->AddResidue(copy_residue_p);

         // apply the shelx afix numbers:
         int unshelxed_natoms;
         mmdb::PAtom *unshelxed_residue_atoms = NULL;
         residue_p->GetAtomTable(unshelxed_residue_atoms, unshelxed_natoms);

         int copy_natoms;
         mmdb::PAtom *copy_residue_atoms = NULL;
         copy_residue_p->GetAtomTable(copy_residue_atoms, copy_natoms);

         if (unshelxed_natoms == copy_natoms) { 
            for (int iat=0; iat<copy_natoms; iat++) {
               int afix;
               int istatus = unshelxed_residue_atoms[iat]->GetUDData(udd_afix_handle, afix);
               if (istatus == mmdb::UDDATA_Ok) { 
                  copy_residue_atoms[iat]->PutUDData(udd_afix_handle_local, afix);
               } else {
                  if (!made_afix_transfer_message) { 
                     std::cout << "ERROR transfering AFIX back" << std::endl;
                     made_afix_transfer_message = 1;
                  }
               }
               mmdb::realtype negative_u;
               int istatus_2 =
                  unshelxed_residue_atoms[iat]->GetUDData(udd_riding_atom_negative_u_value_handle, negative_u);
               if (istatus_2 == mmdb::UDDATA_Ok) {
                  int istatus_3 =
                     copy_residue_atoms[iat]->PutUDData(udd_riding_atom_negative_u_value_handle_local, negative_u);
               } 
            }
         } else {
            std::cout << "ERROR transfering afix back: bad copy number of atoms "
                      << unshelxed_natoms << " " << copy_natoms << std::endl;
         } 
      }
   }

   // need to copy over cell and symmetry info:
   mmdb::realtype a[6];
   mmdb::realtype vol;
   int orthcode;
   mol->GetCell(a[0], a[1], a[2], a[3], a[4], a[5], vol, orthcode);
   shelx_mol->SetCell(a[0], a[1], a[2], a[3], a[4], a[5]);
   mmdb::cpstr sg = mol->GetSpaceGroup();

   if (sg)  
      shelx_mol->SetSpaceGroup(sg); 
   
   shelx_mol->FinishStructEdit();
   shelx_mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);
   return shelx_mol;
}


std::ostream &
coot::operator<<(std::ostream &s, const symm_card_composition_t &sc) {

   s << sc.x_element[0] << " " << sc.y_element[0] << " " << sc.z_element[0] << "\n"
     << sc.x_element[1] << " " << sc.y_element[1] << " " << sc.z_element[1] << "\n"
     << sc.x_element[2] << " " << sc.y_element[2] << " " << sc.z_element[2] << "\n"
     << "translations: "
     << sc.trans_frac(0) << " "
     << sc.trans_frac(1) << " "
     << sc.trans_frac(2) << std::endl;

   return s;
} 
