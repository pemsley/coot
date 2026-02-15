/* coords/phenix-geo-bonds.cc
 * 
 * Copyright 2014 by Medical Research Council
 * Author: Paul Emsley
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

#include "coot-utils/coot-coord-utils.hh"

#include "Bond_lines.hh"

// Phenix Geo
//
Bond_lines_container::Bond_lines_container(mmdb::Manager *mol,
                                           const coot::phenix_geo::phenix_geometry &pg) {

   struct bonded_atom_pair_t {
      mmdb::Atom *atom_1;
      mmdb::Atom *atom_2;
      double residual;
      bonded_atom_pair_t(mmdb::Atom *a1, mmdb::Atom *a2, double r) : atom_1(a1), atom_2(a2), residual(r) {}
   };

   auto atom_colour = [] (double residual) {
      // residual 0 -> colour 60, residual >= 10 -> colour 90, linear in between
      if (residual >= 10.0) return 90;
      if (residual <= 0.0) return 60;
      return static_cast<int>(60.0 + (residual / 10.0) * 30.0);
   };

   // fill this, then make bonds from it.
   std::vector<bonded_atom_pair_t> bonded_atom_pairs;

   coot::residue_spec_t previous_res_spec;
   mmdb::Atom *previous_atom_1 = NULL;
   mmdb::Residue *previous_residue = NULL;
   coot::atom_spec_t previous_atom_spec;
   for (unsigned int i=0; i<pg.geo_bonds.size(); i++) {
      mmdb::Residue *atom_1_res = NULL;
      mmdb::Atom *atom_1 = NULL;
      const coot::phenix_geo::phenix_geo_bond &gb = pg.geo_bonds[i];
      const coot::atom_spec_t &atom_1_spec = gb.atom_1;
      const coot::atom_spec_t &atom_2_spec = gb.atom_2;
      coot::residue_spec_t res_1_spec(atom_1_spec);

      if (res_1_spec == previous_res_spec) {

         if (atom_1_spec == previous_atom_spec) {

            if (previous_atom_1) {
               atom_1 = previous_atom_1;
               atom_1_res = atom_1->residue;
            }

            coot::residue_spec_t res_2_spec(gb.atom_2);

            if (atom_1) {
               mmdb::Atom *atom_2 = coot::util::get_atom(atom_2_spec, atom_1->residue);

               if (res_2_spec == res_1_spec) {

                  // OK both atoms in the same residue

                  if (atom_2) {
                     bonded_atom_pairs.push_back(bonded_atom_pair_t(atom_1, atom_2, gb.residual));
                  } else {
                     std::cout << "Null atom 2 - [A path] this should not happen " << std::endl;
                  }

               } else {

                  // atom_2 was in a different residue
                  //
                  mmdb::Residue *nr = coot::util::next_residue(atom_1_res);
                  mmdb::Residue *pr = coot::util::previous_residue(atom_1_res);

                  coot::residue_spec_t nr_spec(nr);
                  coot::residue_spec_t pr_spec(pr);

                  bool done = false;
                  if (res_2_spec == nr_spec) {
                     mmdb::Atom *atom_2 = coot::util::get_atom(atom_2_spec, nr);
                     done = true;
                     if (atom_2) {
                        bonded_atom_pairs.push_back(bonded_atom_pair_t(atom_1, atom_2, gb.residual));
                     } else {
                        std::cout << "Null atom 2 - [B path] this should not happen " << std::endl;
                     }
                  }

                  if (res_2_spec == pr_spec) {
                     mmdb::Atom *atom_2 = coot::util::get_atom(atom_2_spec, pr);
                     done = true;
                     if (atom_2) {
                        bonded_atom_pairs.push_back(bonded_atom_pair_t(atom_1, atom_2, gb.residual));
                     } else {
                        std::cout << "Null atom 2 - [C path] this should not happen " << std::endl;
                     }
                  }

                  if (! done) {

                     // debug
                     std::cout << "Residue for Atom 2 was not the same Residue as for Atom 1 "
                               << "and not the next or previous residues either  " << std::endl;

                     mmdb::Residue *residue_p = coot::util::get_residue(coot::residue_spec_t(atom_2_spec), mol);
                     mmdb::Atom *atom_2_l = coot::util::get_atom(atom_2_spec, residue_p);
                     done = true;
                     if (atom_2_l) {
                        bonded_atom_pairs.push_back(bonded_atom_pair_t(atom_1, atom_2_l, gb.residual));
                     } else {
                        std::cout << "Null atom 2 - [D path] this should not happen " << std::endl;
                     }
                  }

                  if (! done) {
                     std::cout << "Fail in atom 2 in a different residue - this should not happen"
                               << std::endl;
                  }
               }

            } else {
               std::cout << "Null atom_1 - this should not happen" << std::endl;
            }

         } else { // this atom_1 was not the same as the previous
                  // atom_1, but the residue of this atom_1 is the
                  // same as the residue of the previous atom_1.

            atom_1 = coot::util::get_atom(atom_1_spec, previous_residue);
            if (atom_1) {
               coot::residue_spec_t res_2_spec(gb.atom_2);
               if (res_2_spec == res_1_spec) {
                  mmdb::Atom *atom_2 = coot::util::get_atom(atom_2_spec, atom_1->residue);
                  if (atom_2) {
                     bonded_atom_pairs.push_back(bonded_atom_pair_t(atom_1, atom_2, gb.residual));
                  } else {
                     std::cout << "Null atom 2 - [E path] this should not happen " << std::endl;
                  }
               } else {
                  mmdb::Atom *atom_2 = coot::util::get_atom(atom_2_spec, mol);
                  if (atom_2) {
                     bonded_atom_pairs.push_back(bonded_atom_pair_t(atom_1, atom_2, gb.residual));
                  } else {
                     std::cout << "Null atom 2 - [E path] this should not happen " << std::endl;
                  }
               }

            } else {
               std::cout << "Null atom 2 - [F path] this should not happen " << std::endl;
            }

         }

      } else { // res of this atom_1 was not the same as the the res of the previous atom_1

         atom_1 = coot::util::get_atom(atom_1_spec, mol);

         if (atom_1) {

            res_1_spec = coot::residue_spec_t(atom_1->GetResidue());
            atom_1_res=atom_1->residue;

            coot::residue_spec_t res_2_spec(gb.atom_2);
            if (res_2_spec == res_1_spec) {
               mmdb::Atom *atom_2 = coot::util::get_atom(atom_2_spec, atom_1->residue);
               if (atom_2) {
                  bonded_atom_pairs.push_back(bonded_atom_pair_t(atom_1, atom_2, gb.residual));
               } else {
                  std::cout << "Null atom 2 - [G path] this should not happen " << std::endl;
               }
            } else {
               // the residue of atom_2 was different to the residue of atom_1
               mmdb::Atom *atom_2 = coot::util::get_atom(atom_2_spec, mol);
               if (atom_2) {
                  bonded_atom_pairs.push_back(bonded_atom_pair_t(atom_1, atom_2, gb.residual));
               } else {
                  std::cout << "Null atom 2 - [H path] this should not happen " << std::endl;
               }

            }

         } else {
            std::cout << "Null atom_1 - Path Z " << atom_1_spec << " mol: " << mol 
                      << " - this should not happen" << std::endl;
         }
      }

      // for next round
      previous_atom_spec = atom_1_spec;
      previous_atom_1 = atom_1;
      // do we need to do current_res_spec here too?
      previous_residue = atom_1_res;
      previous_res_spec = coot::residue_spec_t(previous_residue);
   }

   if (false)
      std::cout << "Made " << bonded_atom_pairs.size() << " bonded_atom_pairs"
                << " from " << pg.geo_bonds.size() << " geo-bonds" << std::endl;

   Bond_lines a;
   bonds.push_back(a); // bonded
   bonds.push_back(a); // unbonded

   int uddHnd = mol->RegisterUDInteger(mmdb::UDR_ATOM, "phenix-geo-bond-for-atom");
   bool have_udd_atoms = true;
   set_udd_unbonded(mol, uddHnd);

   int udd_atom_index_handle = mol->GetUDDHandle(mmdb::UDR_ATOM, "atom index");

   for (unsigned int i=0; i<bonded_atom_pairs.size(); i++) {
      const bonded_atom_pair_t &bap = bonded_atom_pairs[i];
      mmdb::Atom *a1 = bap.atom_1;
      mmdb::Atom *a2 = bap.atom_2;
      a1->PutUDData(uddHnd, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
      a2->PutUDData(uddHnd, graphical_bonds_container::BONDED_WITH_STANDARD_ATOM_BOND);
      graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;
      int model_number = a1->GetModel()->GetSerNum();
      int atom_index_1 = -1;
      int atom_index_2 = -1;
      a1->GetUDData(udd_atom_index_handle, atom_index_1);
      a2->GetUDData(udd_atom_index_handle, atom_index_2);
      int col = atom_colour(bap.residual);
      addBond(col, coot::Cartesian(a1->x, a1->y, a1->z),
                   coot::Cartesian(a2->x, a2->y, a2->z),
              cc, model_number, atom_index_1, atom_index_2);
   }

   stars_for_unbonded_atoms(mol, uddHnd);

}


void
Bond_lines_container::set_udd_unbonded(mmdb::Manager *mol, int uddHnd) {
   
   // set all atoms to unbonded initially
   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (! model_p) {
         std::cout << "Null model" << std::endl;
      } else { 
         mmdb::Chain *chain_p;
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            chain_p = model_p->GetChain(ichain);
            if (! chain_p) {
               std::cout << "Null chain" << std::endl;
            } else {
               int nres = chain_p->GetNumberOfResidues();
               mmdb::Residue *residue_p;
               mmdb::Atom *at;
               for (int ires=0; ires<nres; ires++) {
                  residue_p = chain_p->GetResidue(ires);
                  if (! residue_p) {
                     std::cout << "Null residue" << std::endl;
                  } else {
                     int n_atoms = residue_p->GetNumberOfAtoms();
                     for (int iat=0; iat<n_atoms; iat++) {
                        at = residue_p->GetAtom(iat);
                        if (at)
                           at->PutUDData(uddHnd, graphical_bonds_container::NO_BOND);
                     }
                  }
               }
            }
         }
      }
   }
}
   
void
Bond_lines_container::stars_for_unbonded_atoms(mmdb::Manager *mol, int uddHnd) {

   float star_size = 0.2;
   graphics_line_t::cylinder_class_t cc = graphics_line_t::SINGLE;

   coot::Cartesian small_vec_x(star_size, 0.1, 0.1);
   coot::Cartesian small_vec_y(0.1, star_size, 0.1);
   coot::Cartesian small_vec_z(0.1, 0.1, star_size);
   int ic;
   int col = 2;

   // set all atoms to unbonded initially
   for(int imod = 1; imod<=mol->GetNumberOfModels(); imod++) {
      mmdb::Model *model_p = mol->GetModel(imod);
      if (! model_p) {
         std::cout << "Null model" << std::endl;
      } else { 
         mmdb::Chain *chain_p;
         int n_chains = model_p->GetNumberOfChains();
         for (int ichain=0; ichain<n_chains; ichain++) {
            chain_p = model_p->GetChain(ichain);
            if (! chain_p) {
               std::cout << "Null chain" << std::endl;
            } else { 
               int nres = chain_p->GetNumberOfResidues();
               mmdb::Residue *residue_p;
               mmdb::Atom *at;
               for (int ires=0; ires<nres; ires++) { 
                  residue_p = chain_p->GetResidue(ires);
                  if (! residue_p) {
                     std::cout << "Null residue" << std::endl;
                  } else { 
                     int n_atoms = residue_p->GetNumberOfAtoms();
                     for (int iat=0; iat<n_atoms; iat++) {
                        at = residue_p->GetAtom(iat);
                        if (at) { 
                           if (at->GetUDData(uddHnd, ic) == mmdb::UDDATA_Ok) {
                              if (ic == graphical_bonds_container::NO_BOND) {
                                 coot::Cartesian atom_pos(at->x, at->y, at->z);
                                 addBond(col, atom_pos+small_vec_x, atom_pos-small_vec_x, cc, -1, -1, -1, true, true);
                                 addBond(col, atom_pos+small_vec_y, atom_pos-small_vec_y, cc, -1, -1, -1, true, true);
                                 addBond(col, atom_pos+small_vec_z, atom_pos-small_vec_z, cc, -1, -1, -1, true, true);
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
   
