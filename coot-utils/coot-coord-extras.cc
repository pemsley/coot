/* coot-utils/coot-coord-extras.cc
 *
 * Copyright 2004, 2005, 2006, 2007 by The University of York
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


#include <cmath>
#include <stdexcept>
#include <iomanip>
#include <algorithm> // for std::find
#include <queue>
#include <string>
#include <functional>

#include "acedrg-types-for-residue.hh"
#include "string.h"

#include "compat/coot-sysdep.h"
#include "utils/coot-utils.hh"
#include "geometry/mol-utils.hh"
#include "coot-coord-utils.hh"
#include "coot-coord-extras.hh"
#include "atom-tree.hh"
#include "contact-info.hh"

#include "utils/logging.hh"
extern logging logger;


// Return 0 if any of the residues don't have a dictionary entry
// geom_p gets updated to include the residue restraints if necessary
//
std::pair<int, std::vector<std::string> >
coot::util::check_dictionary_for_residues(mmdb::PResidue *SelResidues, int nSelResidues,
                                          coot::protein_geometry *geom_p,
                                          int read_number) {

   int imol_enc = protein_geometry::IMOL_ENC_ANY; // pass this?

   std::pair<int, std::vector<std::string> > r;

   int status;
   int fail = 0; // not fail initially.

   for (int ires=0; ires<nSelResidues; ires++) {
      std::string resname(SelResidues[ires]->name);
      status = geom_p->have_dictionary_for_residue_type(resname, imol_enc, read_number);
      // This bit is redundant now that try_dynamic_add has been added
      // to have_dictionary_for_residue_type():
      if (status == 0) {
         status = geom_p->try_dynamic_add(resname, read_number);
         if (status == 0) {
            fail = 1; // we failed to find it then.
            r.second.push_back(resname);
         }
      }
   }

   if (fail)
      r.first = 0;
   return r;
}


// For use with wiggly ligands, constructed from a minimol residue,
// the get_contact_indices_from_restraints() needs a mmdb::Residue *.
// Caller disposes.
mmdb::Residue *
coot::GetResidue(const minimol::residue &res_in) {

   mmdb::Residue *res = new mmdb::Residue;

   std::string residue_name = res_in.name;
   int seqnum = res_in.seqnum;
   std::string ins_code = res_in.ins_code;
   res->SetResID(residue_name.c_str(),  seqnum, ins_code.c_str());

   for (unsigned int i=0; i<res_in.atoms.size(); i++) {
      coot::minimol::atom mat = res_in.atoms[i];
      mmdb::Atom *at = new mmdb::Atom;
      at->SetAtomName(mat.name.c_str());
      at->SetElementName(mat.element.c_str());
      at->SetCoordinates(mat.pos.x(), mat.pos.y(), mat.pos.z(),
                         mat.occupancy, mat.temperature_factor);
      unsigned int new_length = mat.altLoc.length() +1;
      char *new_alt_loc = new char [new_length];
      // reset new_alt_loc
      for (unsigned int ic=0; ic<new_length; ic++)
         new_alt_loc[ic] = 0;
      strncpy(at->altLoc, mat.altLoc.c_str(), new_length);
      res->AddAtom(at);
   }

   return res;
}




// 200900905 These days, the following may not be needed.
//
// We also now pass regular_residue_flag so that the indexing of the
// contacts is inverted in the case of not regular residue.  I don't
// know why this is necessary, but I have stared at it for hours, this
// is a quick (ugly hack) fix that works.  I suspect that there is
// some atom order dependency in mgtree that I don't understand.
// Please fix (remove the necessity of depending on
// regular_residue_flag) if you know how.
//
std::vector<std::vector<int> >
coot::util::get_contact_indices_from_restraints(mmdb::Residue *residue,
                                                coot::protein_geometry *geom_p,
                                                bool regular_residue_flag,
                                                bool add_reverse_contacts) {

   int nResidueAtoms = residue->GetNumberOfAtoms();
   std::vector<std::vector<int> > contact_indices(nResidueAtoms);
   std::string restype(residue->name);

   int n_monomers = geom_p->size();

   auto residue_has_deuterium_atoms = [] (mmdb::Residue *residue) {
                                         int nResidueAtoms = residue->GetNumberOfAtoms();
                                         bool has_deuterium_atoms = false;
                                         for (int iat=0; iat<nResidueAtoms; iat++) {
                                            mmdb::Atom *atom_p = residue->GetAtom(iat);
                                            if (! atom_p->isTer()) {
                                               std::string atom_ele(atom_p->element);
                                               if (atom_ele == " D") {
                                                  has_deuterium_atoms = true;
                                                  break;
                                               }
                                            }
                                         }
                                         return has_deuterium_atoms;
                                      };

   // this is a horrible and unconventional method of getting the restraints

   bool has_deuterium_atoms = residue_has_deuterium_atoms(residue);

   for (int icomp=0; icomp<n_monomers; icomp++) {
      const dictionary_residue_restraints_t&dict = (*geom_p)[icomp].second;
      if (dict.residue_info.comp_id == restype) {
         for (unsigned int ibr=0; ibr< dict.bond_restraint.size(); ibr++) {
            for (int iat=0; iat<nResidueAtoms; iat++) {
               mmdb::Atom *atom_p = residue->GetAtom(iat);
               std::string at_name(atom_p->GetAtomName());

               if (dict.bond_restraint[ibr].atom_id_1_4c() == at_name) {
                  int ibond_to = -1;  // initially unassigned.
                  std::string at_name_2;
                  for (int iat2=0; iat2<nResidueAtoms; iat2++) {
                     atom_p = residue->GetAtom(iat2);
                     at_name_2 = atom_p->GetAtomName();
                     if (dict.bond_restraint[ibr].atom_id_2_4c() == at_name_2) {
                        ibond_to = iat2;
                        break;
                     }
                  }
                  if (ibond_to > -1 ) {
                     if (add_reverse_contacts == 0) {
                        if (regular_residue_flag) {
                           contact_indices[iat].push_back(ibond_to);  // for ALA etc
                        } else {
                           contact_indices[ibond_to].push_back(iat);  // ligands
                           // contact_indices[iat].push_back(ibond_to);  // ALA etc
                        }
                     } else {
                        // add reverse contacts.
                        contact_indices[ibond_to].push_back(iat);
                        contact_indices[iat].push_back(ibond_to);
                     }
                  }
//                 This spits out the names of Hydrogens, often.
//                  else
//                       std::cout << "failed to find bonded atom "
//                                 << (*geom_p)[icomp].bond_restraint[ibr].atom_id_2_4c()
//                                 << std::endl;
               }


               if (has_deuterium_atoms) {
                  // same again, but change the dictionary atom names on the fly

                  std::string dict_atom_name_1 = dict.bond_restraint[ibr].atom_id_1_4c();
                  if (dict_atom_name_1[0] == 'H') dict_atom_name_1[0] = 'D';
                  if (dict_atom_name_1[1] == 'H') dict_atom_name_1[1] = 'D';
                  if (dict_atom_name_1 == at_name) {
                     int ibond_to = -1;  // initially unassigned.
                     for (int iat2=0; iat2<nResidueAtoms; iat2++) {
                        atom_p = residue->GetAtom(iat2);
                        std::string at_name_2 = atom_p->GetAtomName();
                        std::string dict_atom_name_2 = dict.bond_restraint[ibr].atom_id_2_4c();
                        if (dict_atom_name_2[0] == 'H') dict_atom_name_2[0] = 'D';
                        if (dict_atom_name_2[1] == 'H') dict_atom_name_2[1] = 'D';
                        if (dict_atom_name_2 == at_name_2) {
                           ibond_to = iat2;
                           break;
                        }
                     }
                     if (ibond_to > -1 ) {
                        if (add_reverse_contacts == 0) {
                           if (regular_residue_flag) {
                              contact_indices[iat].push_back(ibond_to);  // for ALA etc
                           } else {
                              contact_indices[ibond_to].push_back(iat);  // ligands
                              // contact_indices[iat].push_back(ibond_to);  // ALA etc
                           }
                        } else {
                           // add reverse contacts.
                           contact_indices[ibond_to].push_back(iat);
                           contact_indices[iat].push_back(ibond_to);
                        }
                     }
                  }
               }
            }
         }
      }
   }
   return contact_indices;
}

std::vector<std::vector<int> >
coot::util::get_contact_indices_from_restraints(mmdb::Residue *residue,
                                                const coot::dictionary_residue_restraints_t &restraints,
                                                bool regular_residue_flag,
                                                bool add_reverse_contacts) {

   int nResidueAtoms = residue->GetNumberOfAtoms();
   std::vector<std::vector<int> > contact_indices(nResidueAtoms);
   mmdb::Atom *atom_p;

   for (unsigned int ibr=0; ibr< restraints.bond_restraint.size(); ibr++) {
      for (int iat=0; iat<nResidueAtoms; iat++) {
         atom_p = residue->GetAtom(iat);
         std::string at_name(atom_p->GetAtomName());
         if (restraints.bond_restraint[ibr].atom_id_1_4c() == at_name ) {
            int ibond_to = -1;  // initially unassigned.
            std::string at_name_2;
            for (int iat2=0; iat2<nResidueAtoms; iat2++) {
               atom_p = residue->GetAtom(iat2);
               at_name_2 = atom_p->GetAtomName();
               if (restraints.bond_restraint[ibr].atom_id_2_4c() == at_name_2 ) {
                  ibond_to = iat2;
                  break;
               }
            }
            if (ibond_to > -1 ) {
               if (add_reverse_contacts == 0) {
                  if (regular_residue_flag) {
                     contact_indices[iat].push_back(ibond_to);  // for ALA etc
                  } else {
                     contact_indices[ibond_to].push_back(iat);  // ligands
                     // contact_indices[iat].push_back(ibond_to);  // ALA etc
                  }
               } else {
                  // add reverse contacts.
                  contact_indices[ibond_to].push_back(iat);
                  contact_indices[iat].push_back(ibond_to);
               }
            }
         }
      }
   }
   return contact_indices;
}


// The atoms of residue_atoms are in the "right" order for not making
// a tree along the main chain.
//
std::vector<std::vector<int> >
coot::util::get_contact_indices_for_PRO_residue(mmdb::PPAtom residue_atoms,
                                                int nResidueAtoms,
                                                coot::protein_geometry *geom_p) {

   std::vector<std::vector<int> > contact_indices(nResidueAtoms);
   mmdb::Atom *atom_p;

   std::pair<bool, coot::dictionary_residue_restraints_t> rest =
      geom_p->get_monomer_restraints("PRO", protein_geometry::IMOL_ENC_ANY);

   if (rest.first) {
      for (unsigned int ibr=0; ibr< rest.second.bond_restraint.size(); ibr++) {
         for (int iat=0; iat<nResidueAtoms; iat++) {
            atom_p = residue_atoms[iat];
            std::string at_name(atom_p->GetAtomName());
            if (rest.second.bond_restraint[ibr].atom_id_1_4c() == at_name ) {
               int ibond_to = -1;  // initially unassigned.
               std::string at_name_2;
               for (int iat2=0; iat2<nResidueAtoms; iat2++) {
                  atom_p = residue_atoms[iat2];
                  at_name_2 = atom_p->GetAtomName();
                  if (rest.second.bond_restraint[ibr].atom_id_2_4c() == at_name_2 ) {
                     ibond_to = iat2;
                     break;
                  }
               }
               if (ibond_to != -1)
                  contact_indices[iat].push_back(ibond_to);
            }
         }
      }
   }
   return contact_indices;
}


coot::util::dict_residue_atom_info_t::dict_residue_atom_info_t(const std::string &residue_name_in,
                                                               coot::protein_geometry *geom_p) {

   residue_name = residue_name_in;

   std::pair<short int, dictionary_residue_restraints_t> p =
      geom_p->get_monomer_restraints(residue_name, protein_geometry::IMOL_ENC_ANY);

   if (p.first) {
      for (unsigned int iat=0; iat<p.second.atom_info.size(); iat++) {
         std::string atom_name = p.second.atom_info[iat].atom_id_4c;
         short int isHydrogen = 0;
         if (p.second.atom_info[iat].type_symbol == "H" ||
             p.second.atom_info[iat].type_symbol == "D") {
            isHydrogen = 1;
         }
         atom_info.push_back(util::dict_atom_info_t(atom_name, isHydrogen));
      }
   }

}

// This one we can do a dynamic add.
//
bool
coot::util::is_nucleotide_by_dict_dynamic_add(mmdb::Residue *residue_p, coot::protein_geometry *geom_p) {

   bool is_nuc = 0;
   std::string residue_name = residue_p->GetResName();

   std::pair<short int, dictionary_residue_restraints_t> p =
      geom_p->get_monomer_restraints(residue_name, protein_geometry::IMOL_ENC_ANY);

   if (p.first) {
      if (p.second.residue_info.group == "RNA" ||
          p.second.residue_info.group == "DNA" ) {
         is_nuc = 1;
      }
   } else {
      int read_number = 40;
      int status = geom_p->try_dynamic_add(residue_name, read_number);
      if (status != 0) {
         // we successfully added it, let's try to run this function
         // again.  Or we could just test the last entry in
         // geom_p->dict_res_restraints(), but it is not public, so
         // it's messy.
         //
         p = geom_p->get_monomer_restraints(residue_name, protein_geometry::IMOL_ENC_ANY);
         if (p.first) {
            if (p.second.residue_info.group == "RNA" ||
                p.second.residue_info.group == "DNA" ) {
               is_nuc = 1;
            }
         }
      }
   }
   return is_nuc;
}


// This one we can NOT do a dynamic add.
//
bool
coot::util::is_nucleotide_by_dict(mmdb::Residue *residue_p, const coot::protein_geometry &geom) {

   bool is_nuc = 0;
   std::string residue_name = residue_p->GetResName();

   std::pair<short int, dictionary_residue_restraints_t> p =
      geom.get_monomer_restraints(residue_name, protein_geometry::IMOL_ENC_ANY);
   if (p.second.residue_info.group == "RNA" ||
       p.second.residue_info.group == "DNA" ) {
      is_nuc = 1;
   }
   return is_nuc;
}



// Move the atoms in res_moving.
//
// Return the number of rotated torsions.
//
coot::match_torsions::match_torsions(mmdb::Residue *res_moving_in, mmdb::Residue *res_ref_in,
                                     const coot::dictionary_residue_restraints_t &rest) {

   res_moving = res_moving_in;
   res_ref = res_ref_in;
   moving_residue_restraints = rest;

}

// Move the atoms in res_moving.
//
// Return the number of rotated torsions.
//
// tr_moving is not used - we use an atom name match to find the
// torsions in the moving residue.
//
int
coot::match_torsions::match(const std::vector <coot::dict_torsion_restraint_t>  &tr_moving,
                            const std::vector <coot::dict_torsion_restraint_t>  &tr_ref) {

   int n_matched = 0; // return value.

   coot::graph_match_info_t match_info = coot::graph_match(res_moving, res_ref, 0, 0);

   if (! match_info.success) {
      // std::cout << "WARNING:: Failed to match graphs " << std::endl;
      logger.log(log_t::WARNING, logging::function_name_t("match"), "Failed to match graphs");
   } else {
      std::string alt_conf = ""; // kludge it in

      // std::map<std::pair<std::string, std::string>, std::pair<std::string, std::string> > atom_name_map;

      // match_info.matching_atom_names have moving molecule first and
      // reference atoms second.  We want to map from reference atom
      // names to moving atom names.
      //
      std::map<std::string, std::string> atom_name_map;
      for (unsigned int i=0; i<match_info.matching_atom_names.size(); i++) {
         atom_name_map[match_info.matching_atom_names[i].second.first] =
            match_info.matching_atom_names[i].first.first;
         if (false)
            std::cout << "      name map construction  :"
                      << match_info.matching_atom_names[i].second.first
                      << ": -> :"
                      << match_info.matching_atom_names[i].first.first
                      << ":\n";
      }

      // for debugging
      std::vector<std::pair<coot::atom_name_quad, double> > check_quads;
      std::vector<double> starting_quad_torsions;

      for (unsigned int itr=0; itr<tr_ref.size(); itr++) {

         coot::atom_name_quad quad_ref(tr_ref[itr].atom_id_1_4c(),
                                       tr_ref[itr].atom_id_2_4c(),
                                       tr_ref[itr].atom_id_3_4c(),
                                       tr_ref[itr].atom_id_4_4c());

         coot::atom_name_quad quad_moving(atom_name_map[tr_ref[itr].atom_id_1_4c()],
                                          atom_name_map[tr_ref[itr].atom_id_2_4c()],
                                          atom_name_map[tr_ref[itr].atom_id_3_4c()],
                                          atom_name_map[tr_ref[itr].atom_id_4_4c()]);

         // test if quad_moving/tr_moving (whatever that is) is a ring torsion
         //
         // if (moving_residue_restraints.is_ring_torsion(quad_moving))


         if (moving_residue_restraints.is_ring_torsion(quad_moving)) {
            // std::cout << "    ignore this ring torsion " << quad_moving << std::endl;
         } else {
            // std::cout << "    OK moving torsion " << quad_moving << " is not a ring torsion"
            // << std::endl;


            if (quad_ref.all_non_blank()) {
               if (quad_moving.all_non_blank()) {

                  if (false)
                     std::cout << "  Reference torsion: "
                               << ":" << tr_ref[itr].format() << " maps to "
                               << quad_moving << std::endl;

                  double starting_quad_tor = quad_moving.torsion(res_moving); // debugging
                  std::pair<bool, double> result = apply_torsion(quad_moving, quad_ref, alt_conf);
                  if (! result.first) {
                     // no tree in restraints? Try without
                     if (false)
                        std::cout << "No tree in match (torsions) - try without" << std::endl;
                     result = apply_torsion_by_contacts(quad_moving, quad_ref, alt_conf);
                  }

                  if (result.first) {
                     n_matched++;
                     // result.second is in radians
                     std::pair<coot::atom_name_quad, double> cq (quad_moving, result.second);
                     check_quads.push_back(cq);
                     starting_quad_torsions.push_back(starting_quad_tor);
                  }
               } else {
                  logger.log(log_t::WARNING, logging::function_name_t("match_torsions::match"),
                             "quad_moving not all non-blank");
                  // let's diagnose that:
                  std::cout << "quad-ref:" << std::endl;
                  std::cout << "   " << tr_ref[itr].atom_id_1_4c() << std::endl;
                  std::cout << "   " << tr_ref[itr].atom_id_2_4c() << std::endl;
                  std::cout << "   " << tr_ref[itr].atom_id_3_4c() << std::endl;
                  std::cout << "   " << tr_ref[itr].atom_id_4_4c() << std::endl;
                  std::cout << "quad-moving:" << std::endl;
                  std::cout << "   " << atom_name_map[tr_ref[itr].atom_id_1_4c()] << std::endl;
                  std::cout << "   " << atom_name_map[tr_ref[itr].atom_id_2_4c()] << std::endl;
                  std::cout << "   " << atom_name_map[tr_ref[itr].atom_id_3_4c()] << std::endl;
                  std::cout << "   " << atom_name_map[tr_ref[itr].atom_id_4_4c()] << std::endl;
               }
            } else {
               std::cout << "WARNING:: in torsion match() quad ref not all non-blank" << std::endl;
               // let's diagnose that:
               std::cout << "quad-ref:" << std::endl;
               std::cout << "   " << tr_ref[itr].atom_id_1_4c() << std::endl;
               std::cout << "   " << tr_ref[itr].atom_id_2_4c() << std::endl;
               std::cout << "   " << tr_ref[itr].atom_id_3_4c() << std::endl;
               std::cout << "   " << tr_ref[itr].atom_id_4_4c() << std::endl;
            }
         }
      }

      // std::cout << "------ after matching, check the torsions " << std::endl;
      logger.log(log_t::INFO, "---- after matching, check the torsions");
      // after matching, check the torsions:
      for (unsigned int iquad=0; iquad<check_quads.size(); iquad++) {
         std::pair<bool, double> mtr = get_torsion(coot::match_torsions::MOVING_TORSION,
                                                   check_quads[iquad].first);
         if (mtr.first) {
            // std::cout << "   torsion check:  "
            //           << check_quads[iquad].first << "  was "
            //           << std::fixed << std::setw(7) << std::setprecision(2)
            //           << starting_quad_torsions[iquad] << " "
            //           << " should be " << std::fixed << std::setw(7) << std::setprecision(2)
            //           << check_quads[iquad].second * 180/M_PI
            //           << " and is "  << std::fixed << std::setw(7) << std::setprecision(2)
            //           << mtr.second * 180/M_PI;

            std::vector<logging::ltw> args = {"torsion check:",
                                              check_quads[iquad].first.format(),
                                              "was",
                                              starting_quad_torsions[iquad],
                                              "should be",
                                              check_quads[iquad].second * 180.0/M_PI,
                                              "and is",
                                              mtr.second * 180.0/M_PI};
            if (fabs(check_quads[iquad].second - mtr.second) > M_PI/180.0)
               args.push_back("  ----- WRONG!!!! ");
            logger.log(log_t::INFO, logging::function_name_t("match_torsions::match"), args);
         }
      }
   }
   return n_matched;
}

// return in radians
std::pair<bool, double>
coot::match_torsions::get_torsion(int torsion_type,
                                  const coot::atom_name_quad &quad) const {

   switch (torsion_type) {

   case coot::match_torsions::REFERENCE_TORSION:
      return get_torsion(res_ref, quad);

   case coot::match_torsions::MOVING_TORSION:
      return get_torsion(res_moving, quad);

   default:
      return std::pair<bool, double> (0,0);
   }
}

std::pair<bool, double>
coot::match_torsions::get_torsion(mmdb::Residue *res, const coot::atom_name_quad &quad) const {

   bool status = 0;
   double tors = 0;
   std::vector<mmdb::Atom *> atoms(4, static_cast<mmdb::Atom *> (NULL));
   atoms[0] = res->GetAtom(quad.atom_name(0).c_str());
   atoms[1] = res->GetAtom(quad.atom_name(1).c_str());
   atoms[2] = res->GetAtom(quad.atom_name(2).c_str());
   atoms[3] = res->GetAtom(quad.atom_name(3).c_str());

   if (atoms[0] && atoms[1] && atoms[2] && atoms[3]) {
      clipper::Coord_orth pts[4];
      for (unsigned int i=0; i<4; i++)
         pts[i] = clipper::Coord_orth(atoms[i]->x, atoms[i]->y, atoms[i]->z);
      tors = clipper::Coord_orth::torsion(pts[0], pts[1], pts[2], pts[3]); // radians
      status = 1;
   }
   return std::pair<bool, double> (status, tors);
}

// Move the atoms of res_moving to match the torsion of res_ref - and
// the torsion of res_ref is determined from the
// torsion_restraint_reference atom names.
//
// The alt conf is the alt conf of the moving residue.
//
// return the torsion which we applied (in radians).
//
std::pair<bool, double>
coot::match_torsions::apply_torsion(const coot::atom_name_quad &moving_quad,
                                    const coot::atom_name_quad &reference_quad,
                                    const std::string &alt_conf) {

   bool status = 0;
   double new_angle = 0;
   std::pair<bool, double> tors = get_torsion(res_ref, reference_quad);
   if (tors.first) {
      try {
         coot::atom_tree_t tree(moving_residue_restraints, res_moving, alt_conf);

         new_angle = tree.set_dihedral(moving_quad.atom_name(0), moving_quad.atom_name(1),
                                       moving_quad.atom_name(2), moving_quad.atom_name(3),
                                       tors.second * 180/M_PI);
         status = 1; // may not happen if set_dihedral() throws an exception
      }
      catch (const std::runtime_error &rte) {
         // std::cout << "WARNING tree-based setting dihedral failed, " << rte.what() << std::endl;
      }
   }
   return std::pair<bool, double> (status, new_angle * M_PI/180.0);
}


std::pair<bool, double>
coot::match_torsions::apply_torsion_by_contacts(const coot::atom_name_quad &moving_quad,
                                                const coot::atom_name_quad &reference_quad,
                                                const std::string &alt_conf) {

   bool status = false;
   double new_angle = 0.0;

   try {
      bool add_reverse_contacts = true;
      std::vector<std::vector<int> > contact_indices =
         coot::util::get_contact_indices_from_restraints(res_moving, moving_residue_restraints, 1, add_reverse_contacts);
      std::pair<bool, double> tors = get_torsion(res_ref, reference_quad);

      int base_atom_index = 0; // hopefully this will work
      coot::minimol::residue ligand_residue(res_moving);
      coot::atom_tree_t tree(moving_residue_restraints, contact_indices, base_atom_index, ligand_residue, alt_conf);
      new_angle = tree.set_dihedral(moving_quad.atom_name(0), moving_quad.atom_name(1),
                                    moving_quad.atom_name(2), moving_quad.atom_name(3),
                                    tors.second * 180/M_PI);

      coot::minimol::residue wiggled_ligand_residue = tree.GetResidue();

      // transfer from the ligand_residue to res_moving
      mmdb::PPAtom residue_atoms = 0;
      int n_residue_atoms;
      int n_transfered = 0;
      res_moving->GetAtomTable(residue_atoms, n_residue_atoms);
      if (int(wiggled_ligand_residue.atoms.size()) <= n_residue_atoms) {
         for (unsigned int iat=0; iat<wiggled_ligand_residue.atoms.size(); iat++) {
            mmdb::Atom *at = res_moving->GetAtom(wiggled_ligand_residue.atoms[iat].name.c_str(), NULL, alt_conf.c_str());
            if (at) {
               if (0)
                  std::cout << "transfering coords was "
                            << at->z << " " << at->y << " " << at->z << " to "
                            << ligand_residue.atoms[iat] << std::endl;
               at->x = wiggled_ligand_residue.atoms[iat].pos.x();
               at->y = wiggled_ligand_residue.atoms[iat].pos.y();
               at->z = wiggled_ligand_residue.atoms[iat].pos.z();
               n_transfered++;
            }
         }
      }
      if (0) {
         std::cout << "-------------------------------- n_transfered " << n_transfered << "------------------- "
                   << std::endl;
         std::cout << "in apply_torsion_by_contacts() new_angle is " << new_angle << std::endl;
      }
      status = 1;
   }
   catch (const std::runtime_error &rte) {
      std::cout << "WARNING:: " << rte.what() << std::endl;
   }
   return std::pair<bool, double> (status, new_angle * M_PI/180.0);
}


// Don't return any hydrogen torsions - perhaps we should make that a
// passed parameter.
//
std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> >
coot::torsionable_bonds_monomer_internal(mmdb::Residue *residue_p,
                                         mmdb::PPAtom atom_selection, int n_selected_atoms,
                                         bool include_pyranose_ring_torsions_flag,
                                         coot::protein_geometry *geom_p) {

   std::vector<std::pair<mmdb::Atom *, mmdb::Atom *> > v;

   bool hydrogen_torsions = false;
   std::string rn = residue_p->GetResName();
   std::vector <dict_torsion_restraint_t> tors_restraints =
      geom_p->get_monomer_torsions_from_geometry(rn, hydrogen_torsions);
   bool is_pyranose = false; // reset maybe
   std::string group = geom_p->get_group(residue_p);
   // CCD dictionaries are marked at D-SACCHARIDE or SACCHARIDE (FUC)
   // SIA (sialic acid) is marked as NON-POLYMER.  Hmm..
   if (group == "pyranose" || group == "D-pyranose" || group == "L-pyranose" ||
       group == "D-SACCHARIDE" || group == "SACCHARIDE")
      is_pyranose = true;

   if (tors_restraints.size()) {
      for (unsigned int itor=0; itor<tors_restraints.size(); itor++) {

         if (! tors_restraints[itor].is_const()) {
            std::string tr_atom_name_2 = tors_restraints[itor].atom_id_2_4c();
            std::string tr_atom_name_3 = tors_restraints[itor].atom_id_3_4c();

            for (int iat1=0; iat1<n_selected_atoms; iat1++) {
               mmdb::Residue *res_1 = atom_selection[iat1]->residue;
               std::string atom_name_1 = atom_selection[iat1]->name;
               for (int iat2=0; iat2<n_selected_atoms; iat2++) {
                  if (iat1 != iat2) {
                     mmdb::Residue *res_2 = atom_selection[iat2]->residue;
                     if (res_1 == res_2) {
                        std::string atom_name_2 = atom_selection[iat2]->name;
                        if (atom_name_1 == tr_atom_name_2) {
                           if (atom_name_2 == tr_atom_name_3) {

                              if ((include_pyranose_ring_torsions_flag == 1) ||
                                  (is_pyranose && !tors_restraints[itor].is_pyranose_ring_torsion(rn)) ||
                                  (! is_pyranose)) {

                                 std::pair<mmdb::Atom *, mmdb::Atom *> p(atom_selection[iat1],
                                                               atom_selection[iat2]);
                                 v.push_back(p);
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
   return v;
}

// The quad version of this (for actually setting torsions)
//
// Don't return any hydrogen torsions - perhaps we should make that a
// passed parameter.
//
std::vector<coot::torsion_atom_quad>
coot::torsionable_bonds_monomer_internal_quads(mmdb::Residue *residue_p,
                                               mmdb::PPAtom atom_selection, int n_selected_atoms,
                                               bool include_pyranose_ring_torsions_flag,
                                               coot::protein_geometry *geom_p) {
   std::vector<coot::torsion_atom_quad> quads;
   bool hydrogen_torsions = false;
   std::string rn = residue_p->GetResName();
   std::vector <dict_torsion_restraint_t> tors_restraints =
      geom_p->get_monomer_torsions_from_geometry(rn, hydrogen_torsions);
   bool is_pyranose = false;
   std::string group = geom_p->get_group(residue_p);
   if (group == "pyranose" || group == "D-pyranose" || group == "L-pyranose" ||
       group == "D-SACCHARIDE" || group == "SACCHARIDE")
      is_pyranose = true;
   std::vector<std::string> residue_alt_confs = coot::util::get_residue_alt_confs(residue_p);
   for (unsigned int itor=0; itor<tors_restraints.size(); itor++) {
      if (! tors_restraints[itor].is_const()) {
         std::string tor_atom_name[5];
         std::vector<mmdb::Atom *> ats(5, static_cast<mmdb::Atom *>(NULL));
         tor_atom_name[1] = tors_restraints[itor].atom_id_1_4c();
         tor_atom_name[2] = tors_restraints[itor].atom_id_2_4c();
         tor_atom_name[3] = tors_restraints[itor].atom_id_3_4c();
         tor_atom_name[4] = tors_restraints[itor].atom_id_4_4c();
         if ((include_pyranose_ring_torsions_flag == 1) ||
             (is_pyranose && !tors_restraints[itor].is_pyranose_ring_torsion(rn)) ||
             (! is_pyranose)) {
            for (unsigned int ialt=0; ialt<residue_alt_confs.size(); ialt++) {
               for (int iat=0; iat<n_selected_atoms; iat++) {
                  std::string atom_name = atom_selection[iat]->name;
                  std::string alt_conf  = atom_selection[iat]->altLoc;
                  if (alt_conf == residue_alt_confs[ialt]) {
                     for (unsigned int jtor=1; jtor<5; jtor++) {
                        if (atom_name == tor_atom_name[jtor])
                           ats[jtor] = atom_selection[iat];
                     }
                  }
               }
               // yes we got for atoms (of matching alt confs)
               if (ats[1] && ats[2] && ats[3] && ats[4]) {
                  coot::torsion_atom_quad q(ats[1],ats[2],ats[3],ats[4],
                                            tors_restraints[itor].angle(),
                                            tors_restraints[itor].esd(),
                                            tors_restraints[itor].periodicity());
                  q.name = tors_restraints[itor].id();
                  q.residue_name = rn;
                  quads.push_back(q);

               }
            }
         }
      }
   }
   return quads;
}


coot::bonded_pair_container_t
coot::linkrs_in_atom_selection(mmdb::Manager *mol, mmdb::PPAtom atom_selection, int n_selected_atoms,
                               protein_geometry *geom_p) {
   coot::bonded_pair_container_t bpc;
#ifdef MMDB_WITHOUT_LINKR
#else
   // normal case
   std::vector<mmdb::Residue *> residues;
   for (int i=0; i<n_selected_atoms; i++) {
      mmdb::Residue *r = atom_selection[i]->residue;
      if (std::find(residues.begin(), residues.end(), r) == residues.end())
         residues.push_back(r);
   }

   bool found = false;
   mmdb::Model *model_p = mol->GetModel(1);
   int n_linkrs = model_p->GetNumberOfLinkRs();
   std::cout << "model has " << n_linkrs << " LINKR records"
             << " and " << model_p->GetNumberOfLinks() << " LINK records"
             << std::endl;
   for (int ilink=1; ilink<=n_linkrs; ilink++) {
      mmdb::PLinkR linkr = model_p->GetLinkR(ilink);
      coot::residue_spec_t link_spec_1(linkr->chainID1,
                                       linkr->seqNum1,
                                       linkr->insCode1);
      coot::residue_spec_t link_spec_2(linkr->chainID2,
                                       linkr->seqNum2,
                                       linkr->insCode2);
      for (unsigned int i=0; i<residues.size(); i++) {
         coot::residue_spec_t spec_1(residues[i]);
         if (spec_1 == link_spec_1) {
            for (unsigned int j=0; j<residues.size(); j++) {
               if (i != j) {
                  coot::residue_spec_t spec_2(residues[j]);
                  if (spec_2 == link_spec_2) {
                     found = true;
                     coot::bonded_pair_t pair(residues[i], residues[j], 0, 0, linkr->linkRID);
                     break;
                  }
               }
            }
         }
         if (found)
            break;
      }
      if (found)
         break;
   }

#endif
   return bpc;
}


// use residues-near-residue to find linked residues
std::vector<mmdb::Residue *>
coot::simple_residue_tree(mmdb::Residue *residue_centre, mmdb::Manager *mol, float close_dist_max) {

   double dist_crit = close_dist_max;
   std::vector<mmdb::Residue *> v;
   std::set<mmdb::Residue *> s;

   std::queue<mmdb::Residue *> q; // what is dequeue? (double-ended)

   q.push(residue_centre);
   s.insert(residue_centre);

   while (q.size()) {
      mmdb::Residue *test_residue = q.front();
      s.insert(test_residue);
      q.pop();

      // OK, what new ones shall we add?
      // Don't add residues that are already in the set. Everything that is in the queue
      // is in the set also.
      std::vector<mmdb::Residue *> residues = residues_near_residue(test_residue, mol, dist_crit);
      for (unsigned int ires=0; ires<residues.size(); ires++) {
         mmdb::Residue *rnr = residues[ires];
         std::set<mmdb::Residue *>::const_iterator it = s.find(rnr);
         if (it == s.end()) {
            q.push(rnr);
            s.insert(rnr);
         }
      }
   }

   std::set<mmdb::Residue *>::const_iterator its;
   for (its=s.begin(); its!=s.end(); its++)
      v.push_back(*its);

   return v;
}

// Is this the wrong file for this util function?
std::vector<std::pair<mmdb::Residue *, double> >
coot::util::CO_orientations(mmdb::Manager *mol) {

   std::vector<std::pair<mmdb::Residue *, double> > scores;
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   if (model_p) {
      int n_chains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<n_chains; ichain++) {
         mmdb::Chain *chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         if (nres > 2) {
            int res_idx_min = 1;
            int res_idx_max = nres-2;
            for (int ires=res_idx_min; ires<res_idx_max; ires++) {
               mmdb::Residue *prev_p = chain_p->GetResidue(ires-1);
               mmdb::Residue *this_p = chain_p->GetResidue(ires);
               mmdb::Residue *next_p = chain_p->GetResidue(ires+1);
               int delta_1 = this_p->GetSeqNum() - prev_p->GetSeqNum();
               int delta_2 = next_p->GetSeqNum() - this_p->GetSeqNum();
               if (delta_1 != 1) continue;
               if (delta_2 != 1) continue;
               if (prev_p && this_p && next_p) {
                  int n_atoms_prev = prev_p->GetNumberOfAtoms();
                  int n_atoms_this = this_p->GetNumberOfAtoms();
                  int n_atoms_next = next_p->GetNumberOfAtoms();
                  mmdb::Atom *prev_O = 0;
                  mmdb::Atom *prev_C = 0;
                  mmdb::Atom *this_O = 0;
                  mmdb::Atom *this_C = 0;
                  mmdb::Atom *next_O = 0;
                  mmdb::Atom *next_C = 0;
                  for (int iat=0; iat<n_atoms_prev; iat++) {
                     mmdb::Atom *at = prev_p->GetAtom(iat);
                     std::string atom_name(at->GetAtomName());
                     std::string alt_conf(at->altLoc);
                     if (alt_conf == "") {
                        if (atom_name == " C  ") prev_C = at;
                        if (atom_name == " O  ") prev_O = at;
                     }
                  }
                  if (! prev_C) continue;
                  if (! prev_O) continue;
                  for (int iat=0; iat<n_atoms_this; iat++) {
                     mmdb::Atom *at = this_p->GetAtom(iat);
                     std::string atom_name(at->GetAtomName());
                     std::string alt_conf(at->altLoc);
                     if (alt_conf == "") {
                        if (atom_name == " C  ") this_C = at;
                        if (atom_name == " O  ") this_O = at;
                     }
                  }
                  if (! this_C) continue;
                  if (! this_O) continue;
                  for (int iat=0; iat<n_atoms_next; iat++) {
                     mmdb::Atom *at = next_p->GetAtom(iat);
                     std::string atom_name(at->GetAtomName());
                     std::string alt_conf(at->altLoc);
                     if (alt_conf == "") {
                        if (atom_name == " C  ") next_C = at;
                        if (atom_name == " O  ") next_O = at;
                     }
                  }
                  if (! next_C) continue;
                  if (! next_O) continue;

                  // OK!
                  clipper::Coord_orth v1(co(prev_O) - co(prev_C));
                  clipper::Coord_orth v2(co(this_O) - co(this_C));
                  clipper::Coord_orth v3(co(next_O) - co(next_C));

                  clipper::Coord_orth v1n(v1.unit());
                  clipper::Coord_orth v2n(v2.unit());
                  clipper::Coord_orth v3n(v3.unit());

                  double dp_1 = clipper::Coord_orth::dot(v1n, v2n);
                  double dp_2 = clipper::Coord_orth::dot(v2n, v3n);
                  double sum = dp_1; // + dp_2;
                  // std::cout << "dp_1 " << dp_1 << " dp_2 " << dp_2 << "\n";
                  std::pair<mmdb::Residue *, double> s(this_p, sum);
                  scores.push_back(s);
               }
            }
         }
      }
   }
   return scores;
}

#include "atom-selection-container.hh"
void
coot::util::parse_prosmart_log_and_gen_CO_plot(const std::string &prosmart_log_file_helix,
                                               const std::string &prosmart_log_file_strand,
                                               const std::string &data_points_file_name,
                                               const std::string &pdb_file_name,
                                               const std::string &chain_id) {

   atom_selection_container_t asc = get_atom_selection(pdb_file_name, false, true, false);
   if (asc.read_success) {
      std::vector<std::pair<mmdb::Residue *, double> > co_scores = CO_orientations(asc.mol);
      std::map<residue_spec_t, double> scores_map;
      std::vector<std::pair<mmdb::Residue *, double> >::const_iterator it;
      for (it=co_scores.begin(); it!=co_scores.end(); ++it) {
         residue_spec_t spec(it->first);
         scores_map[spec] = it->second;
      }
      if (co_scores.size() > 0) {
         std::ifstream f_helix(prosmart_log_file_helix.c_str());
         std::ifstream f_strand(prosmart_log_file_strand.c_str());
         std::ofstream fo(data_points_file_name.c_str());
         std::string pdb_fn = file_name_non_directory(pdb_file_name);
         std::map<residue_spec_t, double> helix_scores;
         std::map<residue_spec_t, double> strand_scores;
         if (f_helix) {
            std::string line;
            while (std::getline(f_helix, line)) {
               std::vector<std::string> bits = util::split_string_on_whitespace_no_blanks(line);
               // std::cout << "read " << line << " from " << prosmart_log_file << " "
               //          << bits.size() << std::endl;
               if (bits.size() == 12) {
                  if (bits[3] == "ALA") { // does this work for strand test also?
                     try {
                        int res_no = string_to_int(bits[0]);
                        float deviation_flexi      = string_to_float(bits[6]);
                        float deviation_procrustes = string_to_float(bits[7]);
                        residue_spec_t res_spec(chain_id, res_no, "");
                        helix_scores[res_spec] = deviation_procrustes;
                     }
                     catch (const std::runtime_error &rte) {
                        // residue number was not a number - oh well
                        std::cout << "something bad parsing " << line  << " " << rte.what()
                                  << std::endl;
                     }
                  }
               }
            }
         }

         if (f_strand) {
            std::string line;
            while (std::getline(f_strand, line)) {
               std::vector<std::string> bits = util::split_string_on_whitespace_no_blanks(line);
               // std::cout << "read " << line << " from " << prosmart_log_file << " "
               //          << bits.size() << std::endl;
               if (bits.size() == 12) {
                  if (bits[3] == "ALA") { // does this work for strand test also?
                     try {
                        int res_no = string_to_int(bits[0]);
                        float deviation_flexi      = string_to_float(bits[6]);
                        float deviation_procrustes = string_to_float(bits[7]);
                        residue_spec_t res_spec(chain_id, res_no, "");
                        strand_scores[res_spec] = deviation_procrustes;
                     }
                     catch (const std::runtime_error &rte) {
                        // residue number was not a number - oh well
                        std::cout << "something bad parsing " << line  << " " << rte.what()
                                  << " " << prosmart_log_file_strand << std::endl;
                     }
                  }
               }
            }
         }

         for (it=co_scores.begin(); it!=co_scores.end(); ++it) {
            residue_spec_t res_spec(it->first);
            const double &co_dp(it->second);
            std::map<residue_spec_t, double>::const_iterator it_helix;
            std::map<residue_spec_t, double>::const_iterator it_strand;
            it_strand = strand_scores.find(res_spec);
            it_helix  = helix_scores.find(res_spec);
            if (it_helix != helix_scores.end()) {
               if (it_strand != strand_scores.end()) {
                  const double &helix_score(it_helix->second);
                  const double &strand_score(it_strand->second);
                  fo << pdb_fn << " " << chain_id << " " << res_spec.res_no
                     << " CO-dp: " << co_dp
                     << " helix: " << helix_score
                     << " strand: " << strand_score
                     << "\n";
               }
            } else {
               if (false) {
                  std::cout << "debug:: failed to find residue " << res_spec
                            << " in " << " helix map of size " << helix_scores.size() << std::endl;
                  std::map<residue_spec_t, double>::const_iterator it;
                  for (it=helix_scores.begin(); it!=helix_scores.end(); ++it) {
                     std::cout << "   " << it->first << " " << it->second << std::endl;
                  }
               }
            }
         }
      }
   }
}

void
coot::util::multi_parse_prosmart_log_and_gen_CO_plot() {
   std::string dir = "pdb";
   std::vector<std::string> sub_dirs = glob_files(dir, "*");
   for (unsigned int i=0; i<sub_dirs.size(); i++) {
      const std::string &sub_dir = sub_dirs[i];
      std::vector<std::string> sub_dir_files = glob_files(sub_dir, "*.pdb");
      for (unsigned int j=0; j<sub_dir_files.size(); j++) {
         const std::string &pdb_file = sub_dir_files[j];
         // std::cout << "pdb_file: " << pdb_file << std::endl;
         std::string code_pdb = util::file_name_non_directory(pdb_file);
         if (code_pdb.length() > 4) code_pdb = code_pdb.substr(0,4);
         std::string code = util::name_sans_extension(code_pdb);
         std::string code_star = code + "_*";
         // std::cout << "pdb file name: " << pdb_file << " code: " << code << "\n";

         std::string prosmart_1 = "../src/ProSMART_Output/Output_Files/Residue_Alignment_Scores";
         std::vector<std::string> chain_files = glob_files(prosmart_1, code_star);

         for (unsigned int k=0; k<chain_files.size(); k++) {
            std::string chain_file = chain_files[k];
            // ha! now I need to strip the directory
            std::string chain_file_file = util::file_name_non_directory(chain_file);
            // std::cout << "chain file file:  " << chain_file_file << std::endl;
            int cff_len = chain_file_file.length();
            if (cff_len == 6 || cff_len == 14) {
               std::string chain_id(chain_file_file.substr(5,1));
               if (cff_len == 14) chain_id = chain_file_file.substr(13,1);
               // std::cout << "chain file " << chain_file << " chain-id: " << chain_id << std::endl;
               std::string data_file_name = chain_file_file + ".data";
               std::string fn_helix  = chain_file_file + "_helix_A.txt";
               std::string fn_strand = chain_file_file + "_strand_A.txt";
               std::string log_file_name_helix  = append_dir_file(chain_file, fn_helix);
               std::string log_file_name_strand = append_dir_file(chain_file, fn_strand);

               parse_prosmart_log_and_gen_CO_plot(log_file_name_helix, log_file_name_strand,
                                                  data_file_name, pdb_file, chain_id);
            }
         }
      }
   }
}


coot::util::missing_atom_info
coot::util::missing_atoms(mmdb::Manager *mol,
                          bool do_missing_hydrogen_atoms_flag,
                          protein_geometry *geom_p) {

   bool ignore_missing_OXT = true;
   bool ignore_missing_OP3 = true;

   std::vector<mmdb::Residue *> residues_with_missing_atoms;
   std::vector<std::string> residues_no_dictionary;
   std::map<mmdb::Residue *, std::vector<std::string> > residue_missing_atom_names_map;
   // and these atoms will need to be deleted when we auto-complete the residue.
   std::vector<std::pair<mmdb::Residue *, std::vector<mmdb::Atom *> > > atoms_in_coords_but_not_in_dict;

   if (mol) {

      // residue_atoms is a vector of residues names (monomer comp_id)
      // together with a vector of atoms.  On each check in this list
      // for the residue type, associated with each atom name is a
      // flag which says if this is a hydrogen or not.
      std::vector<coot::util::dict_residue_atom_info_t> residue_atoms;

      int imod = 1;
      mmdb::Model *model_p = mol->GetModel(imod);
      mmdb::Chain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      for (int ichain=0; ichain<nchains; ichain++) {
         chain_p = model_p->GetChain(ichain);
         int nres = chain_p->GetNumberOfResidues();
         mmdb::PResidue residue_p;
         mmdb::Atom *at;
         for (int ires=0; ires<nres; ires++) {
            residue_p = chain_p->GetResidue(ires);
            std::string residue_name(residue_p->GetResName());
            short int found_dict = 0;
            std::vector<coot::util::dict_atom_info_t> residue_dict_atoms;
            for (unsigned int idict=0; idict<residue_atoms.size(); idict++) {
               std::string tmp_str = residue_atoms[idict].residue_name;
               if (residue_name == tmp_str) {
                  residue_dict_atoms = residue_atoms[idict].atom_info;
                  found_dict = 1;
               }
            }
            if (!found_dict) {
               // we need to try to make it from the dict_info restraint info.
               // (dymamic_add will happen automatically).
               //
               // set found_dict if we get it and make it.
               coot::util::dict_residue_atom_info_t residue_atoms_for_a_residue(residue_name,
                                                                                geom_p);
               if (!residue_atoms_for_a_residue.is_empty_p()) {
                  residue_atoms.push_back(residue_atoms_for_a_residue);
                  residue_dict_atoms = residue_atoms_for_a_residue.atom_info;
                  found_dict = 1;
               }
            }

            if (!found_dict) {

               // No dictionary available, even after we try to make
               // it (dynamic add failed, presumably).  So add this
               // residue type to the list of residues for which we
               // can't find the restraint info:
               residues_no_dictionary.push_back(residue_name);

            } else {

               // OK, we have a dictionary. Deal with hydrogens, by
               // assigning them as present if we don't care if they
               // don't exist.
               //
               // Note: this kludges the isHydrogen? flag.  We are
               // testing if it was set (existed) or not!
               //
               std::vector<coot::util::dict_atom_info_t> dict_atom_names_pairs;
               for (unsigned int iat=0; iat<residue_dict_atoms.size(); iat++) {
                  // Let's not add hydrogens to the dictionary atoms
                  // if we don't care if the don't exist.
                  if ((do_missing_hydrogen_atoms_flag == false) &&
                      (residue_dict_atoms[iat].is_Hydrogen_flag == 1)) {
                     // do nothing
                  } else {
                     // put it in the list and initially mark it as not found.
                     coot::util::dict_atom_info_t p(residue_dict_atoms[iat].name, 0);
                     // PDBv3 FIXME
                     bool really_missing = true;

                     if (ignore_missing_OXT) {
                        if (residue_dict_atoms[iat].name == " OXT")
                           really_missing = false;
                     }
                     if (ignore_missing_OP3) {
                        if (residue_dict_atoms[iat].name == " OP3")
                           really_missing = false;
                     }

                     if (really_missing) {
                        dict_atom_names_pairs.push_back(p);
                     }
                  }
               }

               // OK for every atom in the PDB residue, was it in the dictionary?
               //
               int n_atoms = residue_p->GetNumberOfAtoms();
               for (int iat=0; iat<n_atoms; iat++) {
                  at = residue_p->GetAtom(iat);
                  std::string atom_name(at->name);
                  // check against each atom in the dictionary:
                  for (unsigned int idictat=0; idictat<dict_atom_names_pairs.size(); idictat++) {
                     if (atom_name == dict_atom_names_pairs[idictat].name) {
                        // kludge the is_Hydrogen_flag! Use as a marker
                        dict_atom_names_pairs[idictat].is_Hydrogen_flag = 1; // mark as found
                        break;
                     }
                  }
               }

               // OK, so we have run through all atoms in the PDB
               // residue.  How many atoms in the dictionary list for
               // this residue were not found? Counterintuitive: the
               // is_Hydrogen_flag is used as a marker of being found!
               //
               std::vector<std::string> missing_atom_names;
               for (unsigned int idictat=0; idictat<dict_atom_names_pairs.size(); idictat++) {
                  if (! dict_atom_names_pairs[idictat].is_Hydrogen_flag) {
                     missing_atom_names.push_back(dict_atom_names_pairs[idictat].name);
                  }
               }

               if (! missing_atom_names.empty()) {
                  residues_with_missing_atoms.push_back(residue_p);
                  residue_missing_atom_names_map[residue_p] = missing_atom_names;
               }
            }
         }
      }
   }

   coot::util::missing_atom_info mai(residues_no_dictionary,
                                     residues_with_missing_atoms,
                                     atoms_in_coords_but_not_in_dict);
   mai.residue_missing_atom_names_map = residue_missing_atom_names_map; // a bit kludgy.

   return mai;

}


// mark up things that have omega > 210 or omega < 150. i.e, 180 +/- 30.
//
// strictly_cis_flag is false by default.  If strictly_cis_flag we catch twisted trans too (where twisted
// means that (delta omega) is more than 30 degrees from 180 trans).
//
std::vector<coot::util::cis_peptide_quad_info_t>
coot::cis_peptide_quads_from_coords(mmdb::Manager *mol,
                                    int model_number,
                                    const coot::protein_geometry *geom_p,
                                    bool strictly_cis_flag) {

   std::vector<util::cis_peptide_quad_info_t> v;

   if (!mol) {
      std::cout << "ERROR:: in cis_peptide_quads_from_coords() null mol " << std::endl;
      return v;
   }

   int mol_atom_index_handle = mol->GetUDDHandle(mmdb::UDR_ATOM, "atom index");

   int n_models = mol->GetNumberOfModels();
   if (n_models == 0)
      return v;

   for (int imod=1; imod<=n_models; imod++) {
      if (model_number == 0 || model_number == imod) {
         mmdb::Model *model_p = mol->GetModel(imod);
         if (model_p) {
            mmdb::Chain *chain_p;
            int nchains = model_p->GetNumberOfChains();
            for (int ichain=0; ichain<nchains; ichain++) {
               chain_p = model_p->GetChain(ichain);
               if (chain_p) {
                  int nres = chain_p->GetNumberOfResidues();
                  mmdb::Residue *residue_p_1 = 0;
                  mmdb::Residue *residue_p_2 = 0;
                  mmdb::Atom *at_1 = 0;
                  mmdb::Atom *at_2 = 0;
                  for (int ires=0; ires<(nres-1); ires++) {

                     mmdb::Atom *ca_first = NULL, *c_first = NULL, *n_next = NULL, *ca_next = NULL;
                     residue_p_1 = chain_p->GetResidue(ires);
                     residue_p_2 = chain_p->GetResidue(ires+1);

                     // if (residue_p_2->GetSeqNum() == (residue_p_1->GetSeqNum() + 1)) {
                     if (residue_p_1 && residue_p_2) {

                        bool is_pre_pro = false;

                        std::string res_name_1 = residue_p_1->GetResName();
                        std::string res_name_2 = residue_p_2->GetResName();

                        int n_atoms_1 = residue_p_1->GetNumberOfAtoms();
                        int n_atoms_2 = residue_p_2->GetNumberOfAtoms();

                        if (res_name_2 == "PRO")
                           is_pre_pro = true;

                        for (int iat=0; iat<n_atoms_1; iat++) {
                           at_1 = residue_p_1->GetAtom(iat);
                           if (std::string(at_1->GetAtomName()) == " CA ")
                              ca_first = at_1;
                           if (std::string(at_1->GetAtomName()) == " C  ")
                              c_first = at_1;
                        }

                        for (int iat=0; iat<n_atoms_2; iat++) {
                           at_2 = residue_p_2->GetAtom(iat);
                           if (std::string(at_2->GetAtomName()) == " CA ")
                              ca_next = at_2;
                           if (std::string(at_2->GetAtomName()) == " N  ")
                              n_next = at_2;
                        }

                        if (ca_first && c_first && n_next && ca_next) {

                           // Don't have peptide planes with mixed alt-confs.
                           std::set<std::string> alt_confs;
                           std::string ac[4];
                           ac[0] = ca_first->altLoc;
                           ac[1] = c_first->altLoc;
                           ac[2] = n_next->altLoc;
                           ac[3] = ca_next->altLoc;
                           for (int i=0; i<4; i++)
                              if (!ac[i].empty())
                                 alt_confs.insert(ac[i]);
                           if (alt_confs.size() < 2) {

                              // we don't want to include CISPEPs for residues that
                              // have a TER card between them.
                              //
                              bool is_ter = false;
                              for (int iat=0; iat<n_atoms_1; iat++) {
                                 mmdb::Atom *at = residue_p_1->GetAtom(iat);
                                 if (at->isTer()) {
                                    is_ter = true;
                                    break;
                                 }
                              }
                              if (! is_ter) {
                                 clipper::Coord_orth caf(ca_first->x, ca_first->y, ca_first->z);
                                 clipper::Coord_orth  cf( c_first->x,  c_first->y,  c_first->z);
                                 clipper::Coord_orth can( ca_next->x,  ca_next->y,  ca_next->z);
                                 clipper::Coord_orth  nn(  n_next->x,   n_next->y,   n_next->z);
                                 double tors = clipper::Coord_orth::torsion(caf, cf, nn, can);
                                 double torsion = clipper::Util::rad2d(tors);

                                 // no flags for C-N distances that are more than 2A apart:

                                 double dist = clipper::Coord_orth::length(nn, cf);

                                 if (dist <= 2.0) {

                                    // 20230929-PE is this an expensive test?
                                    // If the geom_p was not set by the caller,
                                    // then this test should pass.

                                    std::string res_1_group;
                                    std::string res_2_group;

                                    if (geom_p) {
                                       try {
                                          // "PRO" group is "P-peptide" (not "peptide")
                                          res_1_group = geom_p->get_group(res_name_1);
                                          res_2_group = geom_p->get_group(res_name_2);
                                       }
                                       catch(const std::runtime_error &rte) {
                                          std::cout << "WARNING:: " << rte.what() << std::endl;
                                       }
                                    }

                                    // std::cout << ":::::::::: geom_p " << geom_p
                                    //           << " residue_p_1 " << coot::residue_spec_t(residue_p_1)
                                    //           << " res_1_group " << res_1_group
                                    //           << " res_2_group " << res_2_group << " is_pre_pro " << is_pre_pro << std::endl;

                                    bool res_1_group_ok = (res_1_group == "peptide");
                                    bool res_2_group_ok = (res_2_group == "peptide");
                                    if (res_2_group == "P-peptide") res_2_group_ok = true;

                                    if (! geom_p || (res_1_group_ok && res_2_group_ok)) {

                                       // put torsion in the range -180 -> + 180
                                       //
                                       if (torsion > 180.0) torsion -= 360.0;
                                       double d = sqrt((cf - nn).lengthsq());
                                       if (d<3.0) { // the residues were close in space, not just close in sequence

                                          util::cis_peptide_quad_info_t::type_t type = util::cis_peptide_quad_info_t::UNSET_TYPE;

                                          double tors_crit = 90.0;
                                          // cis baddies: -90 to +90
                                          if ( (torsion > -tors_crit) && (torsion < tors_crit)) {
                                             if (is_pre_pro)
                                                type = util::cis_peptide_quad_info_t::PRE_PRO_CIS;
                                             else
                                                type = util::cis_peptide_quad_info_t::CIS;
                                          } else {

                                             if (! strictly_cis_flag) {

                                                double tors_twist_delta_max = 30.0; // degrees
                                                // baddies: -150 to +150
                                                if ((torsion > (-180+tors_twist_delta_max)) && (torsion < (180-tors_twist_delta_max)))
                                                   type = util::cis_peptide_quad_info_t::TWISTED_TRANS;

                                             }
                                          }

                                          if (type != util::cis_peptide_quad_info_t::UNSET_TYPE) {
                                             atom_quad q(ca_first, c_first, n_next, ca_next);
                                             int i1 = -1, i2 = -1, i3 = -1, i4 = -1;
                                             ca_first->GetUDData(mol_atom_index_handle, i1);
                                             c_first->GetUDData( mol_atom_index_handle, i2);
                                             n_next->GetUDData(  mol_atom_index_handle, i3);
                                             ca_next->GetUDData( mol_atom_index_handle, i4);
                                             atom_index_quad iq(i1, i2, i3, i4);
                                             util::cis_peptide_quad_info_t qi(q, iq, type);
                                             v.push_back(qi);
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
   return v;
}


#include "atom-tree.hh"

std::vector<mmdb::Residue *>
coot::util::get_dictionary_conformers(const dictionary_residue_restraints_t &restraints,
                                      bool remove_internal_clash_conformers) {

   std::vector<mmdb::Residue *> rv;
   std::string comp_id = restraints.residue_info.comp_id;
   bool include_hydrogen_torsions_flag = false;
   std::vector <coot::dict_torsion_restraint_t> torsion_restraints =
      restraints.get_non_const_torsions(include_hydrogen_torsions_flag);

   std::vector<unsigned int> conformers_per_torsion;
   unsigned int n_conformers = 1; // this gets multipled, not added to
   std::vector <coot::dict_torsion_restraint_t> rotatable_torsions; // fill this
   for (const auto &torsion : torsion_restraints) {
      if (torsion.is_pyranose_ring_torsion(comp_id)) {
         // pass
      } else {
         std::vector<std::vector<std::string> > ring_atoms_sets;
         if (false) { // test is_ring_torsion here
         } else {
            if (torsion.periodicity() > 1) {
               if (! torsion.is_peptide_torsion()) {
                  std::cout << "***************** " << torsion << " is not peptide torsion" << std::endl;
                  rotatable_torsions.push_back(torsion);
                  conformers_per_torsion.push_back(torsion.periodicity());
                  n_conformers *= torsion.periodicity();
               } else {
                  std::cout << "***************** " << torsion << " IS peptide torsion" << std::endl;
               }
            }
         }
      }
   }

   if (false) { // debug
      std::cout << "debug:: in get_dictionary_conformers(): here with rotatable_torsions size "
                << rotatable_torsions.size() << std::endl;
      for (unsigned int i_tor=0; i_tor<rotatable_torsions.size(); i_tor++) {
         std::cout << "   i_tor " << i_tor << " " << rotatable_torsions[i_tor] << std::endl;
      }
   }

   auto debug_torsion_angles = [] (const std::vector<std::vector<double> > &torsion_angles) {
      for (unsigned int i_conf=0; i_conf<torsion_angles.size(); i_conf++) {
         const auto &conformer_set = torsion_angles[i_conf];
         std::cout << " " << std::setw(2) << i_conf << " [" << conformer_set.size() << "] : ";
         for (unsigned int i_tor=0; i_tor<conformer_set.size(); i_tor++)
            std::cout << std::setw(3) << conformer_set[i_tor] << "  ";
         std::cout << std::endl;
      }
   };

   auto periods_to_torsions = [] (const std::vector<int> &torsions_periods,
                                  const std::vector <coot::dict_torsion_restraint_t> &rotatable_torsions) {

      std::vector<double> angles(torsions_periods.size());
      for (unsigned int i=0; i<torsions_periods.size(); i++) {
         const auto &torsion_restraint = rotatable_torsions[i];
         double variant = torsion_restraint.angle() + torsions_periods[i] * 360.0 / static_cast<double>(torsion_restraint.periodicity());
         angles[i] = variant;
      }
      return angles;
   };

   auto make_bond_or_angle_related_pairs = [] (const dictionary_residue_restraints_t &restraints,
                                               mmdb::Residue *residue_p) {

      std::vector<std::pair<int, int> > bond_or_angle_related_pairs;

      mmdb::Atom **residue_atoms_1 = 0;
      int n_residue_atoms_1 = 0;
      residue_p->GetAtomTable(residue_atoms_1, n_residue_atoms_1);

      // bonds
      for (unsigned int i=0; i<restraints.bond_restraint.size(); i++) {
         const auto &bond_restraint = restraints.bond_restraint[i];
         std::string atom_id_1 = bond_restraint.atom_id_1_4c();
         std::string atom_id_2 = bond_restraint.atom_id_2_4c();
         int atom_id_1_idx = -1;
         int atom_id_2_idx = -1;
         bool found_bond = false;
         for (int iat=0; iat<n_residue_atoms_1; iat++) {
            mmdb::Atom *at_1 = residue_atoms_1[iat];
            if (! at_1->isTer()) {
               std::string atom_name_1(at_1->GetAtomName());
               if (atom_name_1 == atom_id_1) {
                  mmdb::Atom **residue_atoms_2 = 0;
                  int n_residue_atoms_2 = 0;
                  residue_p->GetAtomTable(residue_atoms_2, n_residue_atoms_2);
                  for (int jat=0; jat<n_residue_atoms_2; jat++) {
                     mmdb::Atom *at_2 = residue_atoms_2[jat];
                     if (! at_2->isTer()) {
                        std::string atom_name_2(at_2->GetAtomName());
                        if (atom_name_2 == atom_id_2) {
                           atom_id_1_idx = iat;
                           atom_id_2_idx = jat;
                           bond_or_angle_related_pairs.push_back(std::make_pair(atom_id_1_idx, atom_id_2_idx));
                           found_bond = true;
                        }
                     }
                     if (found_bond) break;
                  }
               }
            }
            if (found_bond) break;
         }
      }

      // angles
      for (unsigned int i=0; i<restraints.angle_restraint.size(); i++) {
         const auto &angle_restraint = restraints.angle_restraint[i];
         std::string atom_id_1 = angle_restraint.atom_id_1_4c();
         std::string atom_id_3 = angle_restraint.atom_id_3_4c();
         int atom_id_1_idx = -1;
         int atom_id_3_idx = -1;
         bool found_angle = false;
         for (int iat=0; iat<n_residue_atoms_1; iat++) {
            mmdb::Atom *at_1 = residue_atoms_1[iat];
            if (! at_1->isTer()) {
               std::string atom_name_1(at_1->GetAtomName());
               if (atom_name_1 == atom_id_1) {
                  mmdb::Atom **residue_atoms_2 = 0;
                  int n_residue_atoms_2 = 0;
                  residue_p->GetAtomTable(residue_atoms_2, n_residue_atoms_2);
                  for (int jat=0; jat<n_residue_atoms_2; jat++) {
                     mmdb::Atom *at_2 = residue_atoms_2[jat];
                     if (! at_2->isTer()) {
                        std::string atom_name_2(at_2->GetAtomName());
                        if (atom_name_2 == atom_id_3) {
                           atom_id_1_idx = iat;
                           atom_id_3_idx = jat;
                           bond_or_angle_related_pairs.push_back(std::make_pair(atom_id_1_idx, atom_id_3_idx));
                           found_angle = true;
                        }
                     }
                  }
               }
            }
         }
      }

      if (false) { // debugging
         std::cout << "These " << bond_or_angle_related_pairs.size()  << " pairs are related by bond or angle" << std::endl;
         for (unsigned int ii=0; ii<bond_or_angle_related_pairs.size(); ii++) {
            int idx_1 = bond_or_angle_related_pairs[ii].first;
            int idx_2 = bond_or_angle_related_pairs[ii].second;
            std::cout << std::setw(2) << ii << " " << std::setw(2) << idx_1 << " : " << std::setw(2) << idx_2 << " "
                      << coot::atom_spec_t(residue_atoms_1[idx_1]) << " "
                      << coot::atom_spec_t(residue_atoms_1[idx_2]) << " "
                      << std::endl;
         }
      }

      return bond_or_angle_related_pairs;
   };

   // count *down*, not up.
   auto get_previous_period_set = [&conformers_per_torsion] (std::vector<int> period_set) {
      int period_set_size = period_set.size();
      bool done = false;
      for (int i=period_set.size()-1; i >= 0; i--) {
         if (period_set[i] > 0) {
            period_set[i] -= 1;
            done = true;
         } else {
            // the heart - let's find a previous period that isn't zero.
            int i_prev = i - 1;
            for (int j=i_prev; j >= 0; j--) {
               if (period_set[j] > 0) {
                  period_set[j] -= 1;
                  // and put the next torsions back to max index
                  for (int jj=0; jj<period_set_size; jj++) {
                     if (jj > j)
                        period_set[jj] = conformers_per_torsion[jj] -1;
                  }
                  done = true;
               }
               if (done) break;
            }
         }
         if (done) break;
      }
      return period_set;
   };

   std::function<std::vector<std::vector<int> >(std::vector<int>)> func = [&func, get_previous_period_set] (const std::vector<int> &period_set) {
      bool all_zeros = true;
      for (unsigned int i=0; i<period_set.size(); i++) {
         if (period_set[i] != 0) {
            all_zeros = false;
            break;
         }
      }
      if (all_zeros) {
         // end case
         return std::vector<std::vector<int> > {period_set};
      } else {
         std::vector<int> next_period = get_previous_period_set(period_set);
         std::vector<std::vector<int> > r = func(next_period);
         r.push_back(period_set);
         return r;
      }
   };

   auto get_self_clash = [] (mmdb::Residue *residue_p, const std::vector<std::pair<int, int> > &bond_or_angle_related_pairs) {

      // check non-bonded contacts are not too close. I guess I setup a refinement for this
      // and look through the non-bonded contacts. So the restraints should be passed to this function.

      // 20240817-PE ah, but if I am going to set up a restraints_container_t or a reduced_angle_info_container_t
      // then I will need to move this functionality into the ideal directory. I don't want to do that, so I will
      // make a function here that makes bond-or-angle-related pairs

      bool status = false; // return this

      // This 1-4 distance in a PHE is 2.79A, so we need to be a bit less than that
      float dist_crit = 2.65;

      // check is made both ways around - I guess that could be changed for speed if needed.
      //
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at_1 = residue_atoms[iat];
         std::string ele_1(at_1->element);
         if (ele_1 == " H") continue;
         if (! at_1->isTer()) {
            for (int jat=0; jat<n_residue_atoms; jat++) {
               if (iat == jat) continue;
               mmdb::Atom *at_2 = residue_atoms[jat];
               std::string ele_2(at_2->element);
               if (ele_2 == " H") continue;
               if (! at_2->isTer()) {
                  float dx = at_2->x - at_1->x;
                  float dy = at_2->y - at_1->y;
                  float dz = at_2->z - at_1->z;
                  if (fabsf(dx) < dist_crit) {
                     if (fabsf(dy) < dist_crit) {
                        if (fabsf(dz) < dist_crit) {
                           float dd = dx * dx + dy * dy + dz * dz;
                           float d = sqrtf(dd);
                           if (d < dist_crit) {
                              // 20240817-PE OK now ignore out all "O" and "OXT" clashes
                              // (this is something that might be best under user-control)
                              // but for now - just jam it in.
                              std::string atom_name_1(at_1->GetAtomName());
                              std::string atom_name_2(at_2->GetAtomName());
                              if (atom_name_1 == " O  ") continue;
                              if (atom_name_2 == " O  ") continue;
                              if (atom_name_1 == " OXT") continue;
                              if (atom_name_2 == " OXT") continue;
                              // OK, now check if this is OK because bond or angle
                              bool found = false;
                              for (unsigned int ii=0; ii<bond_or_angle_related_pairs.size(); ii++) {
                                 const int &i_1 = bond_or_angle_related_pairs[ii].first;
                                 const int &i_2 = bond_or_angle_related_pairs[ii].second;
                                 if (iat == i_1) { if (jat == i_2) { found = true; } }
                                 if (iat == i_2) { if (jat == i_1) { found = true; } }
                                 if (found) break;
                              }
                              if (! found) {
                                 status = true;
                                 if (false)
                                    std::cout << "clash! " << coot::atom_spec_t(at_1) << " " << coot::atom_spec_t(at_2)
                                              << " " << d << std::endl;
                              }
                           }
                        }
                     }
                  }
               }
               if (status) break;
            }
         }
         if (status) break;
      }
      return status;
   };

   // rotate_residue_about_torsions(r, rotatable_torsions, t);

   auto rotate_residue_about_torsions = [] (mmdb::Residue *residue_p,
                                            const coot::dictionary_residue_restraints_t &rest,
                                            const std::vector <coot::dict_torsion_restraint_t> &rotatable_torsions,
                                            const std::vector<double> &torsion_angles) {

      if (rotatable_torsions.size() == torsion_angles.size()) {
         for (unsigned int i=0; i<rotatable_torsions.size(); i++) {
            double torsion_angle = torsion_angles[i];
            const auto &torsion_restraint = rotatable_torsions[i];
            std::string atom_name_1 = torsion_restraint.atom_id_1_4c();
            std::string atom_name_2 = torsion_restraint.atom_id_2_4c();
            std::string atom_name_3 = torsion_restraint.atom_id_3_4c();
            std::string atom_name_4 = torsion_restraint.atom_id_4_4c();
            mmdb::Atom *at_1 = nullptr;
            mmdb::Atom *at_2 = nullptr;
            mmdb::Atom *at_3 = nullptr;
            mmdb::Atom *at_4 = nullptr;
            mmdb::Atom **residue_atoms = 0;
            int n_residue_atoms = 0;
            residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
            for (int iat=0; iat<n_residue_atoms; iat++) {
               mmdb::Atom *at = residue_atoms[iat];
               if (! at->isTer()) {
                  std::string atom_name(at->GetAtomName());
                  if (atom_name == atom_name_1) at_1 = at;
                  if (atom_name == atom_name_2) at_2 = at;
                  if (atom_name == atom_name_3) at_3 = at;
                  if (atom_name == atom_name_4) at_4 = at;
               }
            }
            if (at_1 && at_2 && at_3 && at_4) {
               try {

                  // 20240819-PE atom_tree_t is constructed from a restraints that has a tree.
                  // that seems not to be the case for restraints these days
                  coot::atom_quad quad(at_1, at_2, at_3, at_4);
                  if (! rest.tree.empty()) {
                     coot::atom_tree_t tree(rest, residue_p, "");
                     tree.set_dihedral(quad, torsion_angle, false);
                  }
               }
               catch (const std::runtime_error &e) {
                  std::cout << "WARNING::" << e.what() << std::endl;
               }
            } else {
               std::cout << "WARNING:: rotate_residue_about_torsions(): Missing atoms "
                         << at_1 << " " << at_2 << " " << at_3 << " " << at_4 << std::endl;
            }
         }
      }
   };

   auto transfer_coordinates = [] (mmdb::Residue *from_p, mmdb::Residue *to_p) {

      bool debug = false;
      mmdb::Atom **from_residue_atoms = 0;
      int n_from_residue_atoms = 0;
      from_p->GetAtomTable(from_residue_atoms, n_from_residue_atoms);
      for (int iat=0; iat<n_from_residue_atoms; iat++) {
         mmdb::Atom *at_from = from_residue_atoms[iat];
         if (! at_from->isTer()) {

            std::string atom_name_from = at_from->GetAtomName();
            std::string alt_conf_from  = at_from->altLoc;

            mmdb::Atom **to_residue_atoms = 0;
            int n_to_residue_atoms = 0;
            to_p->GetAtomTable(to_residue_atoms, n_to_residue_atoms);
            for (int iat=0; iat<n_to_residue_atoms; iat++) {
               mmdb::Atom *at_to = to_residue_atoms[iat];
               if (! at_to->isTer()) {

                  std::string atom_name_to = at_to->GetAtomName();
                  std::string alt_conf_to  = at_to->altLoc;

                  if (atom_name_from == atom_name_to) {
                     if (alt_conf_from == alt_conf_to) {

                        if (debug) {
                           std::vector<mmdb::realtype> was = {at_to->x, at_to->y, at_to->z};
                           std::cout << "transfered " << coot::atom_spec_t(at_to) << " "
                                     << std::setw(8) << was[0] << " "
                                     << std::setw(8) << was[1] << " "
                                     << std::setw(8) << was[2] << "  now "
                                     << std::setw(8) << at_from->x << " "
                                     << std::setw(8) << at_from->y << " "
                                     << std::setw(8) << at_from->z << " "
                                     << std::endl;
                        }

                        at_to->x = at_from->x;
                        at_to->y = at_from->y;
                        at_to->z = at_from->z;

                        break;
                     }
                  }
               }
            }
         }
      }
   };

   auto rotate_residue_about_torsions_sans_tree = [transfer_coordinates] (mmdb::Residue *residue_p,
                                                      const coot::dictionary_residue_restraints_t &rest,
                                                      const std::vector <coot::dict_torsion_restraint_t> &rotatable_torsions,
                                                      const std::vector<double> &torsion_angles) {


      bool debug =  false;

      // I could use multi-torsion here. Not sure that it's worth it.
      if (rotatable_torsions.size() == torsion_angles.size()) {

         mmdb::Manager *mol = create_mmdbmanager_from_residue(residue_p); // copies residue
         mmdb::Residue *copied_residue_p = get_first_residue(mol);
         std::vector<std::vector<int> > contact_indices =
            get_contact_indices_from_restraints(copied_residue_p, rest, true, true);
         int base_atom_index = 0;
         int SelHnd = mol->NewSelection();
         mol->Select(SelHnd, mmdb::STYPE_ATOM,
                     0, residue_p->GetChainID(),
                     residue_p->GetSeqNum(),  // starting resno, an int
                     residue_p->GetInsCode(), // any insertion code
                     residue_p->GetSeqNum(),  // starting resno, an int
                     residue_p->GetInsCode(), // any insertion code
                     "*", // any residue name
                     "*", // atom name
                     "*", // elements
                     "*", // alt loc.
                     mmdb::SKEY_OR);

         for (unsigned int i=0; i<rotatable_torsions.size(); i++) {
            double torsion_angle = torsion_angles[i];
            const auto &torsion_restraint = rotatable_torsions[i];
            // std::cout << "\ntorsion_angle " << i << " of " << rotatable_torsions.size() << " " << torsion_restraint << std::endl;
            std::string atom_name_1 = torsion_restraint.atom_id_1_4c();
            std::string atom_name_2 = torsion_restraint.atom_id_2_4c();
            std::string atom_name_3 = torsion_restraint.atom_id_3_4c();
            std::string atom_name_4 = torsion_restraint.atom_id_4_4c();
            mmdb::Atom *at_1 = nullptr;
            mmdb::Atom *at_2 = nullptr;
            mmdb::Atom *at_3 = nullptr;
            mmdb::Atom *at_4 = nullptr;
            mmdb::Atom **residue_atoms = 0;
            int n_residue_atoms = 0;
            copied_residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
            for (int iat=0; iat<n_residue_atoms; iat++) {
               mmdb::Atom *at = residue_atoms[iat];
               if (! at->isTer()) {
                  std::string atom_name(at->GetAtomName());
                  // std::cout << "looking for \"" << atom_name_1 << "\" found \"" << atom_name << "\"" << std::endl;
                  if (atom_name == atom_name_1) at_1 = at;
                  if (atom_name == atom_name_2) at_2 = at;
                  if (atom_name == atom_name_3) at_3 = at;
                  if (atom_name == atom_name_4) at_4 = at;
               }
               // std::cout << "debug:: here with at_1 " << at_1 << std::endl;
            }

            if (at_1 && at_2 && at_3 && at_4) {

               if (debug) {
                  std::cout << "in rotate_residue_about_torsions() lambda: " << std::endl;
                  mmdb::Atom **atom_selection = 0;
                  int n_selected_atoms = 0;
                  mol->GetSelIndex(SelHnd, atom_selection, n_selected_atoms);
                  for (int iat=0; iat<n_selected_atoms; iat++) {
                     mmdb::Atom *at = atom_selection[iat];
                     if (! at->isTer()) {
                        std::cout << "     " << iat << " " << atom_selection[iat] << " "
                                  << coot::atom_spec_t(atom_selection[iat]) << std::endl;
                     }
                  }
               }

               coot::atom_quad quad(at_1, at_2, at_3, at_4);
               coot::atom_tree_t tree(contact_indices, base_atom_index, mol, SelHnd);
               // std::cout << "quad: " << quad << std::endl;
               tree.set_dihedral(quad, torsion_angle, false);
               // so we have changed the atoms of copied_residue_p and the caller expects
               // the atoms of residue_p to be moved - so transfer the coordinates
               transfer_coordinates(copied_residue_p, residue_p);

            } else {
               std::cout << "WARNING:: rotate_residue_about_torsions_sans_tree(): Missing atoms "
                         << at_1 << " " << at_2 << " " << at_3 << " " << at_4 << std::endl;
            }
         }
         mol->DeleteSelection(SelHnd);
         delete mol;
      }
   };

   auto print_atom_positions = [] (mmdb::Residue *residue_p, const std::string &lab) {

      mmdb::Atom **residue_atoms = nullptr;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            std::cout << "   " << lab << "  " << iat << " " << coot::atom_spec_t(at) << " "
                      << at->x << " " << at->y << " " << at->z << std::endl;
         }
      }
   };

   // here find which atom index pairs are related by bond or angles.

   std::vector<int> period_set(conformers_per_torsion.size());
   for (unsigned int i=0; i<conformers_per_torsion.size(); i++)
      period_set[i] = conformers_per_torsion[i] - 1;

   std::vector<std::vector<int> > torsions_periods = func(period_set);
   std::vector<std::vector<double> > torsion_angles; // outer index is conformer index
                                                     // inner index is i_tor
   for (unsigned int i=0; i<torsions_periods.size(); i++) {
      const auto &torsion_periods = torsions_periods[i];
      std::vector<double> angles = periods_to_torsions(torsion_periods, rotatable_torsions);
      torsion_angles.push_back(angles);
   }

   // debug_torsion_angles(torsion_angles);

   mmdb::Residue *residue_p = restraints.GetResidue(false, 10.0f);
   std::vector<std::pair<int, int> > bond_or_angle_related_pairs = make_bond_or_angle_related_pairs(restraints, residue_p);

   for (unsigned int i=0; i<torsion_angles.size(); i++) {
      const std::vector<double> &t = torsion_angles[i];
      mmdb::Residue *r = deep_copy_this_residue(residue_p);

      if (! restraints.tree.empty()) {
         // print_atom_positions(r, "pre              ");
         rotate_residue_about_torsions(r, restraints, rotatable_torsions, t);
         // print_atom_positions(r, "post with tree   ");
      } else {
         // print_atom_positions(r, "pre              ");
         rotate_residue_about_torsions_sans_tree(r, restraints, rotatable_torsions, t);
         // print_atom_positions(r, "post without tree");
      }
      bool is_clashing = get_self_clash(r, bond_or_angle_related_pairs);
      if (! is_clashing)
         rv.push_back(r);
      else
         if (false)
            std::cout << "self clash torsion-set " << i << std::endl;
   }
   delete residue_p;
   return rv;
}


// 20240817-PE old scripting function is moved into libcootapi core
int
coot::util::mutate_by_overlap(mmdb::Residue *residue_p, mmdb::Manager *mol,
                              const dictionary_residue_restraints_t &restraints_current_type,
                              const dictionary_residue_restraints_t &restraints_new_type) {

   auto is_in_residue = [] (mmdb::Residue *residue_p, const std::string &atom_name) {

      bool status = false;
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            std::string res_atom_name(at->GetAtomName());
            if (res_atom_name == atom_name) {
               status = true;
               break;
            }
         }
      }
      return status;
   };

   auto reposition_copy_or_delete_atoms = [mol, is_in_residue] (mmdb::Residue *res_mutable,
                                                                mmdb::Residue *residue_ref,
                                                                bool move_O_atom,
                                                                bool is_nucleotide) {
      // first, delete the atoms of res_mutable that are not in residue_ref;
      std::vector<std::string> keep_atoms;
      std::set<mmdb::Atom *> delete_atoms;

      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_ref->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            std::string atom_name(at->GetAtomName());
            keep_atoms.push_back(atom_name);
         }
      }

      mmdb::Atom **residue_atoms_mutable = 0;
      int n_residue_atoms_mutable = 0;
      res_mutable->GetAtomTable(residue_atoms_mutable, n_residue_atoms_mutable);
      for (int iat=0; iat<n_residue_atoms_mutable; iat++) {
         mmdb::Atom *at = residue_atoms_mutable[iat];
         if (! at->isTer()) {
            std::string atom_name(at->GetAtomName());
            if (std::find(keep_atoms.begin(), keep_atoms.end(), atom_name) == keep_atoms.end()) {
               // not found
               delete_atoms.insert(at);
            }
         }
      }

      // delete_atoms.clear();

      for (auto atom : delete_atoms) {
         delete atom;
         atom = NULL;
         res_mutable->TrimAtomTable();
      }

      if (! delete_atoms.empty())
         res_mutable->TrimAtomTable();

      mol->FinishStructEdit();

      if (false) {
         mmdb::Atom **residue_atoms_mutable = 0;
         int n_residue_atoms_mutable = 0;
         res_mutable->GetAtomTable(residue_atoms_mutable, n_residue_atoms_mutable);
         for (int iat=0; iat<n_residue_atoms_mutable; iat++) {
            mmdb::Atom *at = residue_atoms_mutable[iat];
            if (! at->isTer()) {
               std::cout << "debug atom " << iat << " of " << n_residue_atoms_mutable
                         << " " << at << std::endl;
            }
         }
      }

      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            std::string atom_name(at->GetAtomName());

            for (int jat=0; jat<n_residue_atoms_mutable; jat++) {
               mmdb::Atom *at_mutable = residue_atoms_mutable[jat];
               if (! at_mutable->isTer()) {
                  std::string atom_name_mutable(at_mutable->GetAtomName());

                  if (atom_name == atom_name_mutable) {

                     if (atom_name != " O   " || move_O_atom) {

                        at_mutable->x = at->x;
                        at_mutable->y = at->y;
                        at_mutable->z = at->z;
                        if (false)
                           std::cout << "moved atom " << coot::atom_spec_t(at_mutable)
                                     << " to " << at->x << " " << at->y << " " << at->z << std::endl;
                     }
                  }
               }
            }
         }
      }

      // add new
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         if (! at->isTer()) {
            std::string atom_name(at->GetAtomName());
            if (! is_in_residue(res_mutable, atom_name)) {
               mmdb::Atom *at_copy = new mmdb::Atom;
               std::string at_name = at->GetAtomName();
               if (atom_name != " OXT") { // extra atom in an amino acid
                  if (! (atom_name != " OP3" && is_nucleotide)) {  // extra atom in a nucleic acid
                     std::string ele = at->element;
                     if (ele == " H") continue;
                     at_copy->Copy(at);
                     res_mutable->AddAtom(at_copy);
                  }
               }
            }
         }
      }
      mol->FinishStructEdit();
      mol->PDBCleanup(mmdb::PDBCLEAN_SERIAL|mmdb::PDBCLEAN_INDEX);

   };

   auto convert_to_hetatoms = [] (mmdb::Residue *residue_p) {

      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb::Atom *at = residue_atoms[iat];
         at->Het = 1;
      }
   };

   auto both_have_CB = [] (mmdb::Residue *restraints_residue_p, mmdb::Residue *residue_p) {

      bool in_first = false;
      bool in_second = false;
      mmdb::Atom **residue_atoms = 0;
      int n_residue_atoms = 0;
      residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb:: Atom *at = residue_atoms[iat];
         std::string name(at->GetAtomName());
         if (name == " CB ") {
            in_first = true;
            break;
         }
      }
      restraints_residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
      for (int iat=0; iat<n_residue_atoms; iat++) {
         mmdb:: Atom *at = residue_atoms[iat];
         std::string name(at->GetAtomName());
         if (name == " CB ") {
            in_second = true;
            break;
         }
      }
      if (in_first && in_second) return true;
      return false;
   };

   // --- main line

   // note to self: match_ligand_torsions should be here somewhere

   int status = 0;

   bool is_nucl = is_nucleotide(residue_p);
   bool is_aa   = residue_p->isAminoacid();

   bool debug = false;

   mmdb::Residue *restraints_residue_p = restraints_new_type.GetResidue(false, 10.0f);
   if (restraints_residue_p) {
      mmdb::Manager *mol_from_restraints_residue = create_mmdbmanager_from_residue(restraints_residue_p);
      mmdb::Residue *rr = get_first_residue(mol_from_restraints_residue);
      if (rr) {
         if (is_aa) {
            std::vector<lsq_range_match_info_t> lsq_matchers;
            std::vector<std::string> atom_names = {" N  ", " CA ", " C  "};
            if (both_have_CB(restraints_residue_p, residue_p))
               atom_names.push_back(" CB ");

            std::string ref_chain_id = residue_p->GetChainID();
            int ref_res_no = residue_p->GetSeqNum();
            std::string  ref_ins_code = residue_p->GetInsCode();
            std::string matcher_chain_id = restraints_residue_p->GetChainID();
            int matcher_res_no = restraints_residue_p->GetSeqNum();
            std::string matcher_ins_code = restraints_residue_p->GetInsCode();
            for (const auto &atom_name : atom_names) {
               std::string alt_conf;
               lsq_range_match_info_t m(ref_chain_id, ref_res_no, ref_ins_code, atom_name, alt_conf,
                                        matcher_chain_id, matcher_res_no, matcher_ins_code,
                                        atom_name, alt_conf);
               lsq_matchers.push_back(m);
            }

            std::pair<short int, clipper::RTop_orth> rtop_info =
               get_lsq_matrix(mol, mol_from_restraints_residue, lsq_matchers, 1, true);

            if (debug) {
               std::cout << "get_lsq_matrix() returned " << rtop_info.first << std::endl;
               std::cout << "get_lsq_matrix() returned\n" << rtop_info.second.format()
                         << std::endl;
            }

            if (rtop_info.first)
               transform_atoms(restraints_residue_p, rtop_info.second);

            if (debug) { // debugging
               // where is restraints_residue_p now? What is its orientation?
               mmdb::Manager *mmm = create_mmdbmanager_from_residue(restraints_residue_p);
               mmm->WriteCIFASCII("transformed-restraints-residue.pdb");
            }

            // moving and reference
            const auto &tr_ligand  = restraints_new_type.torsion_restraint;
            const auto &tr_res_ref = restraints_current_type.torsion_restraint;
            match_torsions mt(restraints_residue_p, residue_p, restraints_new_type);
            int n_torsions_moved = mt.match(tr_ligand, tr_res_ref);

            if (debug) { // debugging
               std::cout << "debug:: n_torsions_moved " << n_torsions_moved
                         << std::endl;
            }

            if (debug) { // debugging
               // where is restraints_residue_p now? What is its orientation?
               mmdb::Manager *mmm = create_mmdbmanager_from_residue(restraints_residue_p);
               mmm->WriteCIFASCII("restraints-residue-post-torsion-match.pdb");
            }

            // after torsion matching, let's superpose (again)

            mmdb::Manager *mol_from_restraints_residue_2 =
               create_mmdbmanager_from_residue(restraints_residue_p);
            rtop_info = get_lsq_matrix(mol, mol_from_restraints_residue_2, lsq_matchers, 1, true);
            if (rtop_info.first)
               transform_atoms(restraints_residue_p, rtop_info.second);

            if (false) { // debugging
               // where is restraints_residue_p now? What is its orientation?
               mmdb::Manager *mmm = create_mmdbmanager_from_residue(restraints_residue_p);
               mmm->WriteCIFASCII("restraints-residue-post-second-lsq.pdb");
            }

            // now copy or replace to coordinates the atoms of restraints_residue_p into residue_p
            // and remove atoms in residue_p that are not in restraints_residue_p

            reposition_copy_or_delete_atoms(residue_p, restraints_residue_p, false, false);

            if (debug) {
               mmdb::Atom **residue_atoms = 0;
               int n_residue_atoms = 0;
               residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
               for (int iat=0; iat<n_residue_atoms; iat++) {
                  mmdb::Atom *at = residue_atoms[iat];
                  if (! at->isTer()) {
                     std::cout << "residue atom " << iat << " " << at
                               << " " << atom_spec_t(at) << std::endl;
                  }
               }
            }

            std::string new_residue_name = restraints_new_type.residue_info.comp_id;
            if (! util::is_standard_amino_acid_name(new_residue_name))
               convert_to_hetatoms(residue_p);

            residue_p->SetResName(new_residue_name.c_str());

            mol->FinishStructEdit(); // needed?
            pdbcleanup_serial_residue_numbers(mol);

            delete mol_from_restraints_residue;
            delete mol_from_restraints_residue_2;

            status = 1;
         }

         if (is_nucl) {

            std::vector<std::string> purine_set = {" N9 ", " N7 ", " C5 ", " N1 ", " N3 "};
            std::vector<std::string> pyrimidine_set = {" N1 ", " C5 ", " N3 "};
            std::vector<std::string> purine_to_pyrimidine_set = {" N1 ", " C2 ", " N3 "};
            std::vector<std::string> pyrimidine_to_purine_set = {" N9 ", " C4 ", " N5 "};

            std::pair<bool, clipper::RTop_orth> rtop_info =
               nucleotide_to_nucleotide(residue_p, rr, false);
            if (rtop_info.first)
               transform_atoms(restraints_residue_p, rtop_info.second);

            const auto &tr_ligand  = restraints_new_type.torsion_restraint;
            const auto &tr_res_ref = restraints_current_type.torsion_restraint;
            match_torsions mt(restraints_residue_p, residue_p, restraints_new_type);
            int n_torsions_moved = mt.match(tr_ligand, tr_res_ref);

            // 20241115-PE this is more simple than the above block - do I need
            // to do the same sort of thing here too?

            reposition_copy_or_delete_atoms(residue_p, restraints_residue_p, true, true);

            status = 1;
         }
      }
   }
   delete restraints_residue_p;
   return status;
}

coot::acedrg_types_for_residue_t
coot::get_acedrg_types_for_residue(mmdb::Residue *residue_p, int imol_enc,
                                   const coot::protein_geometry &geom) {

   coot::acedrg_types_for_residue_t types;
   std::string residue_type = residue_p->GetResName();
   auto r = geom.get_monomer_restraints(residue_type, imol_enc);
   if (r.first) {
      const auto &restraints = r.second;
      for (unsigned int ib=0; ib<restraints.bond_restraint.size(); ib++) {
         const auto &bond_restraint = restraints.bond_restraint[ib];
         const std::string &atom_name_1 = bond_restraint.atom_id_1();
         const std::string &atom_name_2 = bond_restraint.atom_id_2();
         mmdb::Atom **residue_atoms = 0;
         int n_residue_atoms = 0;
         int idx_1 = -1;
         int idx_2 = -1;
         residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
         for (int iat=0; iat<n_residue_atoms; iat++) {
            mmdb::Atom *at = residue_atoms[iat];
            if (! at->isTer()) {
               std::string atom_name = at->GetAtomName();
               if (atom_name == atom_name_1) idx_1 = iat;
               if (atom_name == atom_name_2) idx_2 = iat;
            }
         }
         if (idx_1 != -1) {
            if (idx_2 != -1) {
               clipper::Coord_orth at_pos_1 = co(residue_atoms[idx_1]);
               clipper::Coord_orth at_pos_2 = co(residue_atoms[idx_2]);
               double bb = (at_pos_2 - at_pos_1).lengthsq();
               double bond_length = std::sqrt(bb);

               std::string type_1;
               std::string type_2;

               for (unsigned int ii=0; ii<restraints.atom_info.size(); ii++) {
                  const auto &atom = restraints.atom_info[ii];
                  if (atom.atom_id_4c == atom_name_1) type_1 = atom.acedrg_atom_type;
                  if (atom.atom_id_4c == atom_name_2) type_2 = atom.acedrg_atom_type;
               }

               if (! type_1.empty()) {
                  if (! type_2.empty()) {
                     bool in_same_ring_flag = restraints.in_same_ring(atom_name_1, atom_name_2);
                     acedrg_types_for_bond_t bt(atom_name_1, atom_name_2, type_1, type_2, bond_length, in_same_ring_flag);
                     types.bond_types.push_back(bt);
                  }
               }
            }
         }
      }
   }
   return types;
}
