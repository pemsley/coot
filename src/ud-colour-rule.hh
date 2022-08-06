#ifndef UD_COLOUR_RULE_HH
#define UD_COLOUR_RULE_HH

#include <MoleculesToTriangles/CXXClasses/ColorRule.h>
#include "utils/colour-holder.hh"

class ud_colour_rule : public ColorRule {
public:
   std::vector<coot::colour_holder> user_defined_colours;
   int udd_handle;
   int size; // because int vs unsigned int.
   mmdb::Manager *mol;
   ud_colour_rule(int udd_handle_in, mmdb::Manager *mol_in, const std::vector<coot::colour_holder> &user_defined_colours_in) : ColorRule(), user_defined_colours(user_defined_colours_in) {
      // user_defined_colours = user_defined_colours_in;
      udd_handle = udd_handle_in;
      size = user_defined_colours.size();
      mol = mol_in;
      type = 1;
   }
   FCXXCoord colorForAtom(const mmdb::Atom *atom) {
      mmdb::Atom *nc_atom = const_cast<mmdb::Atom *> (atom);
      int idx = -1;
      if (nc_atom->GetUDData(udd_handle, idx) == mmdb::UDDATA_Ok) {
         if (idx>=0 && idx<size) {
            coot::colour_holder col = user_defined_colours[idx];
            FCXXCoord fxcol(col.red, col.green, col.blue, 1.0);
            return fxcol;
         } else {
            return FCXXCoord(0.3, 0.3, 0.5, 1.0);
         }
      } else {
         return FCXXCoord(0.3, 0.3, 0.3, 1.0);
      }
   }
};



#endif // UD_COLOUR_RULE_HH
