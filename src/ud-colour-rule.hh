/*
 * src/ud-colour-rule.hh
 *
 * Copyright 2022 by Medical Research Council
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
#ifndef UD_COLOUR_RULE_HH
#define UD_COLOUR_RULE_HH

#include <MoleculesToTriangles/CXXClasses/ColorRule.h>
#include "utils/colour-holder.hh"

class ud_colour_rule : public ColorRule {
public:
   // this is the the ideal container!
   std::vector<std::pair<unsigned int, coot::colour_holder> > user_defined_colours;
   int udd_handle;
   int max_colour_index;
   mmdb::Manager *mol;
   ud_colour_rule(int udd_handle_in, mmdb::Manager *mol_in, const std::vector<std::pair<unsigned int, coot::colour_holder> > &user_defined_colours_in) :
      ColorRule(), user_defined_colours(user_defined_colours_in) {

      user_defined_colours = user_defined_colours_in;
      udd_handle = udd_handle_in;
      max_colour_index = -1;
      // size = user_defined_colours.size();
      // size = -1;
      for (unsigned int i=0; i<user_defined_colours.size(); i++) {
         int colour_index = user_defined_colours[i].first; // change type
         if (colour_index > max_colour_index) max_colour_index = user_defined_colours[i].first;
      }
      mol = mol_in;
      type = 1;
   }
   FCXXCoord colorForAtom(const mmdb::Atom *atom) {
      mmdb::Atom *nc_atom = const_cast<mmdb::Atom *> (atom);
      int idx = -1;
      // std::cout << "debug:: here in colorForAtom(): A " << std::endl;
      if (nc_atom->GetUDData(udd_handle, idx) == mmdb::UDDATA_Ok) {
         // std::cout << "debug:: here in colorForAtom(): B with idx " << idx << std::endl;
         if (idx>=0 && idx<=max_colour_index) {
            // std::cout << "debug:: here in colorForAtom(): C colorForAtom(): extracted idx " << idx << std::endl;
            bool found = false;
            for (unsigned int i=0; i<user_defined_colours.size(); i++) {
               if (static_cast<int>(user_defined_colours[i].first) == idx) {
                  const coot::colour_holder &col = user_defined_colours[i].second;
                  // std::cout << "debug:: here in colorForAtom(): D colorForAtom(): found idx " << idx
                  // << " col " << col<< std::endl;
                  FCXXCoord fxcol(col.red, col.green, col.blue, 1.0);
                  return fxcol;
               }
            }
            if (! found) {
               // std::cout << "no match for atom " << std::endl;
               return FCXXCoord(0.3, 0.1, 0.1, 1.0);
            }
         } else {
            return FCXXCoord(0.3, 0.3, 0.5, 1.0);
         }
      } else {
         return FCXXCoord(0.3, 0.3, 0.3, 1.0);
      }
      return FCXXCoord(0.3, 0.2, 0.2, 1.0);
   }
};



#endif // UD_COLOUR_RULE_HH
