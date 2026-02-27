/* coot-utils/coot-coord-extras.cc
 * 
 * Copyright 2006, by The University of York
 * Copyright 2009 by The University of Oxford
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

#include "utils/coot-utils.hh"
#include "coot-coord-utils.hh"

#include "utils/logging.hh"
extern logging logger;

// LSQing
//
std::pair<short int, clipper::RTop_orth>
coot::util::get_lsq_matrix(mmdb::Manager *mol1,
                           mmdb::Manager *mol2,
                           const std::vector<coot::lsq_range_match_info_t> &matches,
                           int every_nth,
                           bool summary_to_screen) {

   short int istat = 0;
   clipper::RTop_orth rtop(clipper::Mat33<double>(0,0,0,0,0,0,0,0,0),
                           clipper::Coord_orth(0,0,0));
   int SelHnd1 = mol1->NewSelection(); // d
   int SelHnd2 = mol2->NewSelection(); // d

   std::vector<clipper::Coord_orth> co1v;
   std::vector<clipper::Coord_orth> co2v;
   for (unsigned int i=0; i<matches.size(); i++) {
      std::pair<std::vector<clipper::Coord_orth>, std::vector<clipper::Coord_orth> > p =
        get_matching_indices(mol1, mol2, matches[i], every_nth);
      if ((p.first.size() > 0) && (p.first.size() == p.second.size())) {
         if (false) { // debugging
            if (p.first.size() == 1)
               std::cout << "    " << matches[i] << " " << p.first[0].format() << " "
                         << p.second[0].format() << std::endl;
         }
         for (unsigned int j=0; j<p.first.size(); j++) {
            co1v.push_back(p.first[j]);
            co2v.push_back(p.second[j]);
         }
      }
   }

   // Now convert v1 and v2 to Coord_orth vectors and call Clipper's
   // LSQ function.
   //
   if (co1v.size() > 0) {
      if (co1v.size() > 2) {
         if (co2v.size() > 2) {
            logger.log(log_t::INFO, "LSQ matched", co1v.size(), "atoms");
            rtop = clipper::RTop_orth(co2v, co1v);
            double sum_dist = 0.0;
            double sum_dist2 = 0.0;
            double mind =  999999999.9;
            double maxd = -999999999.9;
            double d;
            for (unsigned int i=0; i<co2v.size(); i++) {
               d = clipper::Coord_orth::length(co1v[i],
                                               clipper::Coord_orth(co2v[i].transform(rtop)));
               sum_dist  += d;
               sum_dist2 += d*d;
               if (d>maxd)
                  maxd = d;
               if (d<mind)
                  mind = d;
            }
            double mean = sum_dist/double(co2v.size());
            // not variance about mean, variance from 0.
            //double var  = sum_dist2/double(co2v.size()) - mean*mean;
            double v    = sum_dist2/double(co2v.size());

            logger.log(log_t::INFO, co1v.size(), "matched atoms with");
            logger.log(log_t::INFO, "   mean devi", mean);
            logger.log(log_t::INFO, "    rms devi", sqrt(v));
            logger.log(log_t::INFO, "    max devi", maxd);
            logger.log(log_t::INFO, "    min devi", mind);

            if (summary_to_screen) {
               // std::cout << "INFO:: " << co1v.size() << " matched atoms had: \n"
               //           << "   mean devi: " << mean << "\n"
               //           << "    rms devi: " << sqrt(v) << "\n"
               //           << "    max devi: " << maxd << "\n"
               //           << "    min devi: " << mind << std::endl;
               logger.log(log_t::INFO, co1v.size(), "matched atoms had:");
               logger.log(log_t::INFO, "  mean devi:", mean);
               logger.log(log_t::INFO, "   rms devi:", sqrt(v));
               logger.log(log_t::INFO, "   max devi:", maxd);
               logger.log(log_t::INFO, "   min devi:", mind);
            }
            istat = 1;
         } else {
            std::cout << "WARNING:: not enough points to do matching (matching)"
                      << std::endl;
         }
      } else {
         std::cout << "WARNING:: not enough points to do matching (reference)"
                   << std::endl;
      }
   } else {
      std::cout << "WARNING:: no points to do matching" << std::endl;
   }
   mol1->DeleteSelection(SelHnd1);
   mol2->DeleteSelection(SelHnd2);
   return std::pair<short int, clipper::RTop_orth> (istat, rtop);
}

// On useful return, first.size() == second.size() and first.size() > 0.
//
std::pair<std::vector<clipper::Coord_orth>, std::vector<clipper::Coord_orth> >
coot::util::get_matching_indices(mmdb::Manager *mol1,
                                 mmdb::Manager *mol2,
                                 const coot::lsq_range_match_info_t &match,
                                 int every_nth) {

   std::vector<clipper::Coord_orth> v1;
   std::vector<clipper::Coord_orth> v2;

   // general main chain atom names
   std::vector<std::string> mc_at_names;
   // Set up the amino acid main chain atom names
   std::vector<std::string> amc_at_names;
   amc_at_names.push_back(" CA ");
   amc_at_names.push_back(" N  ");
   amc_at_names.push_back(" O  ");
   amc_at_names.push_back(" C  ");

   std::vector<std::string> aa_core_mc_at_names; // no O (for LSQ fitting of single residues),
   aa_core_mc_at_names.push_back(" CA ");
   aa_core_mc_at_names.push_back(" N  ");
   aa_core_mc_at_names.push_back(" C  ");

   // Set up the nucleotide main chain atom names
   std::vector<std::string> nmc_at_names;
   nmc_at_names.push_back(" P  ");
   nmc_at_names.push_back(" O5*");
   nmc_at_names.push_back(" O5'");
   nmc_at_names.push_back(" C5*");
   nmc_at_names.push_back(" C5'");
   nmc_at_names.push_back(" C4*");
   nmc_at_names.push_back(" C4'");
   nmc_at_names.push_back(" O4*");
   nmc_at_names.push_back(" O4'");
   nmc_at_names.push_back(" C1*");
   nmc_at_names.push_back(" C1'");
   nmc_at_names.push_back(" C2*");
   nmc_at_names.push_back(" C2'");
   nmc_at_names.push_back(" O2*");
   nmc_at_names.push_back(" O2'");
   nmc_at_names.push_back(" C3*");
   nmc_at_names.push_back(" C3'");
   nmc_at_names.push_back(" O3*");
   nmc_at_names.push_back(" O3'");
   nmc_at_names.push_back(" O3T");

   std::vector<std::string> warnings;

   if (false) { //debugging atom selection
      int imod = 1;
      mmdb::Model *model_p = mol2->GetModel(imod);
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
                        std::cout << "  mol1 " << coot::atom_spec_t(at) << std::endl;
                     }
                  }
               }
            }
         }
      }
   }

//   for (int ires=match.to_reference_start_resno; ires<=match.to_reference_end_resno; ires++) {
   if (every_nth < 1 || every_nth > 10)
     every_nth = 1;   // reset for nonsense values
   for (int ires=match.to_reference_start_resno; ires<=match.to_reference_end_resno; ires+=every_nth) {
      int ires_matcher = ires - match.to_reference_start_resno + match.from_matcher_start_resno;
      int SelHnd_res1 = mol1->NewSelection();
      int SelHnd_res2 = mol2->NewSelection();
      mmdb::PResidue *SelResidue_1 = NULL;
      mmdb::PResidue *SelResidue_2 = NULL;
      int nSelResidues_1, nSelResidues_2;

      if (false) {
         std::cout << "Searching for residue number " << ires << " "
                   << match.reference_chain_id << " in reference molecule" << std::endl;
         std::cout << "Searching for residue number " << ires_matcher << " "
                   << match.matcher_chain_id << " in matcher molecule" << std::endl;
      }

      mol1->Select (SelHnd_res1, mmdb::STYPE_RESIDUE,
                    match.model_number_reference,
                    match.reference_chain_id.c_str(), // Chain(s)
                    ires, "*",  // starting res
                    ires, "*",  // ending res
                    "*",  // residue name
                    "*",  // Residue must contain this atom name?
                    "*",  // Residue must contain this Element?
                    "*",  // altLocs
                    mmdb::SKEY_NEW // selection key
                    );
      mol2->Select (SelHnd_res2, mmdb::STYPE_RESIDUE, // .. TYPE,
                    match.model_number_matcher,
                    match.matcher_chain_id.c_str(), // Chain(s)
                    ires_matcher, "*",  // starting res
                    ires_matcher, "*",  // ending res
                    "*",  // residue name
                    "*",  // Residue must contain this atom name?
                    "*",  // Residue must contain this Element?
                    "*",  // altLocs
                    mmdb::SKEY_NEW // selection key
                    );


      mol1->GetSelIndex(SelHnd_res1, SelResidue_1, nSelResidues_1);
      mol2->GetSelIndex(SelHnd_res2, SelResidue_2, nSelResidues_2);

      if (false) {
         std::cout << "in get_matching_indices() mol1 " << mol1 << std::endl;
         std::cout << "in get_matching_indices() mol2 " << mol2 << std::endl;
         std::cout << "mol1 uses " << ires << std::endl;
         std::cout << "mol2 uses " << ires_matcher << std::endl;
         std::cout << ":::::::::: nSelResidues_1 " << nSelResidues_1
                   << " nSelResidues_2 " << nSelResidues_2 << std::endl;
      }

      if (nSelResidues_1 == 0 || nSelResidues_2 == 0) {

         if (nSelResidues_1 == 0) {
            std::string s = "WARNING:: - no residue for reference molecule residue number ";
            s += util::int_to_string(ires);
            s += " for reference chain-id: \"";
            s += match.reference_chain_id;
            s += "\"";
            warnings.push_back(s);
         }
         if (nSelResidues_2 == 0) {
            std::string s = "WARNING:: - no residue for moving    molecule residue number ";
            s += util::int_to_string(ires_matcher);
            s += " for reference chain-id: \"";
            s += match.reference_chain_id;
            s += "\"";
            warnings.push_back(s);
         }
      } else {

         // So we have 2 good residues

         // Let's get their residue type:
         std::string res_type_1(SelResidue_1[0]->GetResName());
         std::string res_type_2(SelResidue_2[0]->GetResName());

         // ---------------- CA ----------------------------

         if (match.match_type_flag == lsq_t::CA) {

            // CA/P names
            std::string ca_name;

            mmdb::Atom *at1 = NULL;
            mmdb::Atom *at2 = NULL;
            if (SelResidue_1[0]->isAminoacid()) {
              ca_name = " CA ";
            } else {
              if (SelResidue_1[0]->isNucleotide()) {
                 ca_name = " P  ";
              } else {
                ca_name = " P  ";
                std::cout << "WARNING:: residue is not amino acid or nucleotide! "
                          << "Assuming non-standard nucleotide." << std::endl;
              }
            }
            at1 = SelResidue_1[0]->GetAtom(ca_name.c_str());
            at2 = SelResidue_2[0]->GetAtom(ca_name.c_str());

            if (at1 == NULL) {
              std::cout << "WARNING:: no " << ca_name << " in this reference residue " << ires
                        << std::endl;
            }
            if (at2 == NULL) {
              std::cout << "WARNING:: no " << ca_name << " in this reference residue " << ires
                        << std::endl;
            }
            if (at1 && at2) {
               v1.push_back(clipper::Coord_orth(at1->x, at1->y, at1->z));
               v2.push_back(clipper::Coord_orth(at2->x, at2->y, at2->z));
            }
         }

         // ----------------- Mainchain ---------------------
         //
         if ((match.match_type_flag == lsq_t::MAIN) ||
             (match.match_type_flag == lsq_t::NCAC) ||
             (match.match_type_flag == lsq_t::ALL && res_type_1 != res_type_2)) {

            if (SelResidue_1[0]->isNucleotide()) {
               mc_at_names = nmc_at_names;
            } else {
               mc_at_names = amc_at_names;
               if (match.match_type_flag == lsq_t::NCAC)
                  mc_at_names = aa_core_mc_at_names;
            }
            if (! match.is_single_atom_match) {
               for (unsigned int iat=0; iat<mc_at_names.size(); iat++) {
                  mmdb::Atom *at1 = SelResidue_1[0]->GetAtom(mc_at_names[iat].c_str());
                  mmdb::Atom *at2 = SelResidue_2[0]->GetAtom(mc_at_names[iat].c_str());
                  if (at1) {
                     if (!at2) {
                        std::cout << "WARNING:: no " << mc_at_names[iat]
                                  << " in this moving residue " << ires_matcher
                                  << std::endl;
                     } else {
                        if (false)
                           std::cout << "Found match "
                                     << atom_spec_t(at1) << " to "
                                     << atom_spec_t(at2) << std::endl;
                        v1.push_back(clipper::Coord_orth(at1->x, at1->y, at1->z));
                        v2.push_back(clipper::Coord_orth(at2->x, at2->y, at2->z));
                     }
                  }
               }
            }
         }

         // NCACBC
         if (match.match_type_flag == lsq_t::NCACBC) {
            if (! match.is_single_atom_match) {
               std::vector<std::string> atom_names = {" N  ", " CA ", " CB ", " C  "}; // PDBV3-FIXME
               for (unsigned int iat=0; iat<atom_names.size(); iat++) {
                  const std::string &atom_name = atom_names[iat];
                  mmdb::Atom *at1 = SelResidue_1[0]->GetAtom(atom_name.c_str());
                  mmdb::Atom *at2 = SelResidue_2[0]->GetAtom(atom_name.c_str());
                  if (at1) {
                     if (at2) {
                        if (! at1->isTer()) {
                           if (! at2->isTer()) {
                              v1.push_back(clipper::Coord_orth(at1->x, at1->y, at1->z));
                              v2.push_back(clipper::Coord_orth(at2->x, at2->y, at2->z));
                           }
                        }
                     }
                  }
               }
            }
         }

         // ----------------- All Atom ---------------------
         //
         if (match.match_type_flag == lsq_t::ALL) {

            if (! match.is_single_atom_match) {
               mmdb::PAtom *residue_atoms1 = NULL;
               mmdb::PAtom *residue_atoms2 = NULL;
               int n_residue_atoms1;
               int n_residue_atoms2;
               SelResidue_1[0]->GetAtomTable(residue_atoms1, n_residue_atoms1);
               SelResidue_2[0]->GetAtomTable(residue_atoms2, n_residue_atoms2);
               for (int iat=0; iat<n_residue_atoms1; iat++) {
                  mmdb::Atom *at1 = residue_atoms1[iat];
                  std::string at1_name(at1->name);
                  std::string at1_altconf(at1->altLoc);
                  for (int jat=0; jat<n_residue_atoms2; jat++) {
                     mmdb::Atom *at2 = residue_atoms2[jat];
                     std::string at2_name(at2->name);
                     std::string at2_altconf(at2->altLoc);

                     if (at1_name == at2_name) {
                        if (at1_altconf == at2_altconf) {
                           v1.push_back(clipper::Coord_orth(at1->x, at1->y, at1->z));
                           v2.push_back(clipper::Coord_orth(at2->x, at2->y, at2->z));
                           break;
                        }
                     }
                  }
               }
            } else {

               // is single atom match
               mmdb::PAtom *residue_atoms1 = NULL;
               mmdb::PAtom *residue_atoms2 = NULL;
               int n_residue_atoms1;
               int n_residue_atoms2;
               SelResidue_1[0]->GetAtomTable(residue_atoms1, n_residue_atoms1);
               SelResidue_2[0]->GetAtomTable(residue_atoms2, n_residue_atoms2);
//                std::cout << "single atom match: residue: ref: "
//                          << coot::residue_spec_t(SelResidue_1[0])
//                          << " to "
//                          << coot::residue_spec_t(SelResidue_2[0])
//                          << std::endl;
               for (int iat=0; iat<n_residue_atoms1; iat++) {
                  mmdb::Atom *at1 = residue_atoms1[iat];
                  std::string at1_name(at1->name);
                  std::string at1_altconf(at1->altLoc);
                  if (at1_name == match.reference_atom_name) {
                     if (at1_altconf == match.reference_alt_conf) {
                        for (int jat=0; jat<n_residue_atoms2; jat++) {
                           mmdb::Atom *at2 = residue_atoms2[jat];
                           std::string at2_name(at2->name);
                           std::string at2_altconf(at2->altLoc);
                           if (at2_name == match.matcher_atom_name) {
                              if (at2_altconf == match.matcher_alt_conf) {
                                 v1.push_back(clipper::Coord_orth(at1->x, at1->y, at1->z));
                                 v2.push_back(clipper::Coord_orth(at2->x, at2->y, at2->z));
                              }
                           }
                        }
                     }
                  }
               }
            }
         }
      }

      // tidy up
      mol1->DeleteSelection(SelHnd_res1);
      mol2->DeleteSelection(SelHnd_res2);
   }

   if (! warnings.empty()) {
      unsigned int n_warnings = 10;
      if (warnings.size() < n_warnings) n_warnings = warnings.size();
      for (unsigned int ii=0; ii<n_warnings; ii++)
         std::cout << warnings[ii] << "\n";
      if (warnings.size() > n_warnings)
         std::cout << "WARNING:: ... and others...\n";
   }

   return std::pair<std::vector<clipper::Coord_orth>, std::vector<clipper::Coord_orth> > (v1, v2);
}
