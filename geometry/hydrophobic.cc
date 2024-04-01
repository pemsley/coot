/*
 * geometry/hydrophobic.cc
 *
 * Copyright 2020 by Medical Research Council
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
