
#include "mol-utils.hh"


std::vector<std::string>
coot::util::get_residue_alt_confs(mmdb::Residue *res) {

   std::vector<std::string> v;
   mmdb::PPAtom residue_atoms;
   int nResidueAtoms;
   res->GetAtomTable(residue_atoms, nResidueAtoms);
   bool ifound = 0;
   for (int iat=0; iat<nResidueAtoms; iat++) {
      ifound = 0;
      for(unsigned int i=0; i<v.size(); i++) {
         if (std::string(residue_atoms[iat]->altLoc) == v[i]) {
            ifound = 1;
            break;
         }
      }
      if (! ifound)
         v.push_back(std::string(residue_atoms[iat]->altLoc));
   }
   return v;
}


