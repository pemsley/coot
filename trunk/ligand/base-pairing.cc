/* ligand/base-pairing.cc
 * 
 * Copyright 2006 by The University of York
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
 * 02110-1301, USA.
 */

#include "coot-coord-utils.hh"
#include "ideal-rna.hh"
#include "base-pairing.hh"

CResidue *
coot::watson_crick_partner(CResidue *res_ref, CMMDBManager *standard_residues) {

   CResidue *res = NULL;

   std::string resname = res_ref->GetResName();
   std::string seq;

   if (resname == "Gd")
      seq = "g";
   if (resname == "Ad")
      seq = "a";
   if (resname == "Td")
      seq = "t";
   if (resname == "Cd")
      seq = "c";
   
   coot::ideal_rna ir("DNA", "B", 0, seq, standard_residues);
   CMMDBManager *mol = ir.make_molecule();

   if (mol) { 
      int imod = 1;
      CModel *model_p = mol->GetModel(imod);
      CChain *chain_p;
      // run over chains of the existing mol
      int nchains = model_p->GetNumberOfChains();
      if (nchains == 2) { 
	 chain_p = model_p->GetChain(1);
	 int nres = chain_p->GetNumberOfResidues();
	 if (nres == 1) { 
	    res = chain_p->GetResidue(0);
	 }
      }
   }

   return res;
} 

