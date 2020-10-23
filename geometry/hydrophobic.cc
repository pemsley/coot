

#include "hydrophobic.hh"

bool
coot::is_hydrophobic_atom(const std::string &residue_name, const std::string &atom_name) {

   bool status = false;

   if (atom_name[1] == 'N') return false;
   if (atom_name[1] == 'O') return false;

   // test by resisdue type only for now

   // -: ASP GLU
   // +: LYS ARG HIS
   // H: ASN GLN SER THR TYR

   if (residue_name == "GLY")
      status = true;
   else
      if (residue_name == "ALA")
         status = true;
      else
         if (residue_name == "VAL")
            status = true;
         else
            if (residue_name == "LEU")
               status = true;
            else
               if (residue_name == "ILE")
                  status = true;
               else
                  if (residue_name == "PRO")
                     status = true;
                  else
                     if (residue_name == "PHE")
                        status = true;
                     else
                        if (residue_name == "MET")
                           status = true;
                        else
                           if (residue_name == "TRP")
                              status = true;
                           else
                              if (residue_name == "CYS")
                                 status = true;
                              else
                                 if (residue_name == "TYR")      // ?
                                    status = true;
   return status;
}

bool
coot::is_hydrophobic_atom(mmdb::Atom *at) {

   std::string atom_name(at->GetAtomName());
   std::string res_name(at->residue->GetResName());
   return is_hydrophobic_atom(res_name, atom_name);
}
