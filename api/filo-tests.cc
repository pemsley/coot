/*
* @author Filomeno Sanchez
* 
* Copyright 2023, University of York
* All rights reserved.
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published
* by the Free Software Foundation; either version 3 of the License, or (at
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

#include "filo-tests.hh"

int test_something_filo(molecules_container_t &mc) {

   starting_test(__FUNCTION__);
   int status = 0;

   std::cout << "BEFORE " << mc.map_sampling_rate << "\n";

   mc.set_map_sampling_rate(1.8);

   std::cout << "AFTER " << mc.map_sampling_rate << "\n";

   int imol     = mc.read_pdb(reference_data("moorhen-tutorial-structure-number-1.pdb"));
   int imol_map = mc.read_mtz(reference_data("moorhen-tutorial-map-number-1.mtz"), "FWT", "PHWT", "W", false, false);

   if (mc.is_valid_model_molecule(imol)) {
      coot::atom_spec_t atom_spec("A", 270, "", " O  ","");
      mmdb::Atom *at_1 = mc.get_atom(imol, atom_spec);
      if (at_1) {
         coot::Cartesian atom_pos = atom_to_cartesian(at_1);
         double dd = coot::Cartesian::lengthsq(atom_pos, atom_pos);
         double d = std::sqrt(dd);
         std::cout << "test_ d " << d << std::endl;

         coot::validation_information_t dca = mc.density_correlation_analysis(imol, imol_map);
         for (const auto &chain : dca.cviv) {
            for (const auto &res : chain.rviv) {
               if (res.residue_spec.res_no == 62) {
                  std::cout << "function value " << res.function_value << std::endl;
                  if (res.function_value > 0.6) {
                     status = 1;
                  }
               }
            }
         }
      }
   }
   mc.close_molecule(imol);
   return status;
}

