/*
 * coot-utils/find-water-baddies.cc
 *
 * Copyright 2023 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */

#include "find-water-baddies.hh"
#include <string.h> // for strncmp()
#include "coot-map-utils.hh"
#include "utils/coot-utils.hh"

std::vector <coot::atom_spec_t>
coot::find_water_baddies_OR(atom_selection_container_t atom_sel,
                            float b_factor_lim, const clipper::Xmap<float> &xmap_in,
                            float map_in_sigma,
                            float outlier_sigma_level,
                            float min_dist, float max_dist,
                            short int ignore_part_occ_contact_flag,
                            short int ignore_zero_occ_flag) {

   std::vector <coot::atom_spec_t> v;

   std::vector<std::pair<mmdb::Atom *, float> > marked_for_display;

   bool this_is_marked;
   float den = 0.0;
   short int use_b_factor_limit_test = 1;
   short int use_map_sigma_limit_test = 1;
   short int use_min_dist_test = 1;
   short int use_max_dist_test = 1;
   bool sigma_warned = 0; // we only want to see this message once (at most!), not 184 times.

   if (b_factor_lim < 0.0)
      use_b_factor_limit_test = 0;
   if (outlier_sigma_level < -50.0)
      use_map_sigma_limit_test = 0;
   if (min_dist < 0.0)
      use_min_dist_test = 0;
   if ( max_dist < 0.0 )
      use_max_dist_test = 0;

//    std::cout << "DEBUG:: passed: b_factor_lim "
//              << b_factor_lim << " outlier_sigma_level "
//               << outlier_sigma_level << " min_dist "
//              << min_dist << " max_dist " << max_dist
//               << std::endl;

//     std::cout << "DEBUG:: Usage flags b-factor: " << use_b_factor_limit_test
//              << " map sigma: " << use_map_sigma_limit_test
//               << " min dist: " << use_min_dist_test
//               << " max_dist: " << use_max_dist_test << std::endl;

   if (atom_sel.n_selected_atoms > 0) {

      int n_models = atom_sel.mol->GetNumberOfModels();
      for (int imod=1; imod<=n_models; imod++) {

         mmdb::Model *model_p = atom_sel.mol->GetModel(imod);
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
                  mmdb::Atom *at;
                  for (int ires=0; ires<nres; ires++) {
                     residue_p = chain_p->GetResidue(ires);
                     int n_atoms = residue_p->GetNumberOfAtoms();

                     std::string resname = residue_p->GetResName();
                     if (resname == "WAT" || resname == "HOH") {

                        for (int iat=0; iat<n_atoms; iat++) {
                           at = residue_p->GetAtom(iat);
			   bool water_atom_is_hydrogen_atom = false;
			   // PDBv3 FIXME
			   if (! strncmp(at->name, " H", 2)) water_atom_is_hydrogen_atom = true;
			   if (! strncmp(at->name, " D", 2)) water_atom_is_hydrogen_atom = true;

			   if (water_atom_is_hydrogen_atom) continue;

                           this_is_marked = false;

                           if (! at->isTer()) {

                              if (false)
                                 std::cout << "found water atom " << coot::atom_spec_t(at) << std::endl;

                              // density check:
                              if (map_in_sigma > 0.0) { // it *should* be!
                                 clipper::Coord_orth a(at->x, at->y, at->z);
                                 den = coot::util::density_at_point(xmap_in, a);

                                 den /= map_in_sigma;
                                 if (den < outlier_sigma_level && use_map_sigma_limit_test) {
                                    this_is_marked = true;
                                    marked_for_display.push_back(std::pair<mmdb::Atom *, float>(at, den));
                                 }
                              } else {
                                 if (! sigma_warned) {
                                    std::cout << "Ooops! Map sigma is " << map_in_sigma << std::endl;
                                    sigma_warned = true;
                                 }
                              }

                              // B factor check:
                              if (! this_is_marked) {
                                 if (at->tempFactor > b_factor_lim && use_b_factor_limit_test) {
                                    marked_for_display.push_back(std::pair<mmdb::Atom *, float>(at, den));
				 }
                              }


                              // distance check
                              if (! this_is_marked) {

                                 // (ignoring things means less marked atoms)
                                 if (ignore_part_occ_contact_flag == 0) {

                                    if (ignore_zero_occ_flag == false || at->occupancy > 0.01) {

                                       double dist_to_atoms_min = 99999;
                                       double dc_sqrd = dist_to_atoms_min * dist_to_atoms_min;
                                       double d_sqrd_min = 999999999;
                                       clipper::Coord_orth a(at->x, at->y, at->z);
                                       for (int j=0; j<atom_sel.n_selected_atoms; j++) {
                                          if (at != atom_sel.atom_selection[j]) {
                                             bool is_H = false;
                                             // PDB v3 FIXME?
                                             if (! strncmp(atom_sel.atom_selection[j]->element, " H", 2))
                                                is_H = true;

                                             if (! is_H) {
                                                clipper::Coord_orth p(atom_sel.atom_selection[j]->x,
                                                                      atom_sel.atom_selection[j]->y,
                                                                      atom_sel.atom_selection[j]->z);
                                                double d_sqrd = (p-a).lengthsq();
                                                if (d_sqrd < d_sqrd_min) {
                                                   d_sqrd_min = d_sqrd;
                                                }
                                             }
                                          }
                                       }
                                       bool dist_to_atoms_min_is_set = false;
                                       if (d_sqrd_min < 999999998) {
                                          // should be!
                                          dist_to_atoms_min_is_set = true;
                                          dist_to_atoms_min = sqrt(d_sqrd_min);
                                       }
                                       bool failed_min_dist_test = false;
                                       bool failed_max_dist_test = false;

                                       if (dist_to_atoms_min_is_set && (dist_to_atoms_min < min_dist) && use_min_dist_test)
                                          failed_min_dist_test = true;

                                       if (dist_to_atoms_min_is_set && (dist_to_atoms_min > max_dist) && use_max_dist_test)
                                          failed_max_dist_test = true;

                                       if (failed_min_dist_test || failed_max_dist_test) {

                                          // std::cout << "mark for display " << coot::atom_spec_t(at) << std::endl;
                                          marked_for_display.push_back(std::pair<mmdb::Atom *, float>(at, den));
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

   std::cout << "marked_for_display size() " << marked_for_display.size() << std::endl;

   for (unsigned int i=0; i<marked_for_display.size(); i++) {
      std::string s = "B fac: ";
      s += coot::util::float_to_string(marked_for_display[i].first->tempFactor);
      if (map_in_sigma > 0.0) {
         s += "   ED: ";
         s += coot::util::float_to_string(marked_for_display[i].second);
         s += " rmsd";
      }
      coot::atom_spec_t as(marked_for_display[i].first, s);
      as.float_user_data = marked_for_display[i].first->occupancy;
      v.push_back(as);
   }
   return v;
}
