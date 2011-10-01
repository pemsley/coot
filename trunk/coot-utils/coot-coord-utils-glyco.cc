/* coot-utils/coot-coord-utils-glyco.cc
 * 
 * Copyright 2011 by The University of York
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

#include "coot-utils.hh"
#include "coot-coord-utils.hh"

// Note: this is a simple-minded hack.  The right way of doing this
// is to define a bonding tree that includes atoms from both
// residues.  Then we don't need reference structures - the
// "moving" residue atoms will get placed by internal coordinates.
//
// This can be viewed as starting (or test material) for the Proper
// Way code.
// 
// Scenario Simple Beam-in:
//    User has an ASN onto which they want to beam in a NAG.
// 
//    setup:
//    Make CResidue *s and/or molecule for the N-linked NAG reference residues.
// 
///   Get the CResidue * for the user residue ASN 
//    Get the CAtoms *s for the OD1, ND2, GC and CB in user residue [1]
//    Get the CAtoms *s for the OD1, ND2, GC and CB in N-linked ASN molecule [2]
//
//    LSQ fit the NAG residue from the reference ASN+NAG pair using
//    matrix that rotates [2] onto [1].  (We don't need the ASN from
//    that pair).  Now we can add that rotated NAG CResidue * to user
//    molecule.  we have N-linked-NAG template - consider renaming to
//    ASN-NAG-via-NAG-ASN i.e. the general case
//    {ResType1}-{ResType2}-via-{LinkName} where ResType1 and ResType2
//    are comp-ids, of course.  Actually, NAG-ASN is a pyranose-ASN
//    link (group to comp_id). Hmm...
coot::beam_in_linked_residue::beam_in_linked_residue(CResidue *residue_ref_in,
						     const std::string &link_type,
						     const std::string &new_residue_type) {
   
   // do we have a template structure for the given args?
   have_template = false;
   template_res_ref = NULL;
   template_res_mov = NULL;

   if (residue_ref_in) {
      residue_ref = residue_ref_in;
      comp_id_ref = residue_ref->GetResName();
      std::string file_name = comp_id_ref + "-" + new_residue_type ;
      file_name += "-via-";
      file_name += link_type;
      file_name += ".pdb";

      std::string pkgdatadir = coot::package_data_dir();
      std::string full_path_pdb_filename = pkgdatadir; // and then add to it...
      full_path_pdb_filename += "/";
      full_path_pdb_filename += file_name;

      std::cout << "DEBUG:: Looking for template link file-name: "
		<< full_path_pdb_filename << std::endl;

      if (coot::file_exists(full_path_pdb_filename)) {
	 // Cool.
	 CMMDBManager *t_mol = new CMMDBManager;
	 int status = t_mol->ReadPDBASCII(full_path_pdb_filename.c_str());
	 if (status != Error_NoError) {
	    std::cout << "ERROR:: on reading " << full_path_pdb_filename << std::endl;
	 } else { 
	    // More cool.
	    template_res_ref = get_residue(comp_id_ref, t_mol);
	    if (! template_res_ref) {
	       std::cout << "ERROR:: failed to find residue with comp_id " << comp_id_ref
			 << " in " << full_path_pdb_filename << std::endl;
	    } else { 

	       // should be this path
	       
	       template_res_mov = get_residue(new_residue_type, t_mol);

	       if (! template_res_mov) {
		  std::cout << "ERROR:: failed to find residue with comp_id " << new_residue_type
			    << " in " << full_path_pdb_filename << std::endl;
	       } else { 
		  // Happy path
		  have_template = 1; // template_res_mov and
				     // template_res_ref are correctly
				     // set.
	       } 
	    }
	 } 
      } 
   } 
}

// This can return NULL if we were unable to make the residue to be attached.
CResidue *
coot::beam_in_linked_residue::get_residue() const {

   CResidue *r = NULL;
   if (have_template) {
      // this mean that comp_id_ref, template_res_ref and
      // template_res_mov are correctly set.

      // depends on the comp_id of the residue_ref.
      // 
      std::vector<std::string> lsq_reference_atom_names = 
	 make_reference_atom_names(comp_id_ref);
      if (lsq_reference_atom_names.size() == 0) {
	 std::cout << "WARNING:: no reference atoms for LSQing" << std::endl;
      } else { 
	 std::vector<CAtom *> va_1 = get_atoms(template_res_ref, lsq_reference_atom_names);
	 std::vector<CAtom *> va_2 = get_atoms(residue_ref,      lsq_reference_atom_names);

	 if (va_1.size() != lsq_reference_atom_names.size()) {
	    std::cout << "Mismatch atoms length for " << comp_id_ref << " in "
		      << "get_residue() (c.f. reference atoms) "
		      << va_1.size() << " need " << lsq_reference_atom_names.size() << std::endl;
	 } else {
	    if (va_1.size() != va_2.size()) {
	       std::cout << "Mismatch atoms length for " << comp_id_ref << " in "
			 << "get_residue()" << std::endl;
	    } else {
	       // Happy path
	       int n = lsq_reference_atom_names.size();
	       std::vector<clipper::Coord_orth> co_1(n);
	       std::vector<clipper::Coord_orth> co_2(n);
	       for (unsigned int iat=0; iat<va_1.size(); iat++) { 
		  co_1[iat] = clipper::Coord_orth(va_1[iat]->x, va_1[iat]->y, va_1[iat]->z);
		  co_2[iat] = clipper::Coord_orth(va_2[iat]->x, va_2[iat]->y, va_2[iat]->z);
	       }
	       clipper::RTop_orth rtop(co_1, co_2);
	       coot::util::transform_atoms(template_res_mov, rtop);
	       r = template_res_mov;
	    } 
	 } 
      }
   } else {
      std::cout << "WARNING:: no template" << std::endl;
   }
   return r;
}

std::vector<CAtom *>
coot::beam_in_linked_residue::get_atoms(CResidue *residue_p,
					const std::vector<std::string> &names) const {

   std::vector<CAtom *> v;
   PPCAtom residue_atoms = 0;
   int n_residue_atoms;
   residue_p->GetAtomTable(residue_atoms, n_residue_atoms);
   for (unsigned int iname=0; iname<names.size(); iname++) {
      for (unsigned int iat=0; iat<n_residue_atoms; iat++) {
	 std::string atom_name(residue_atoms[iat]->GetAtomName());
	 if (atom_name == names[iname]) {
	    v.push_back(residue_atoms[iat]);
	    break;
	 } 
      }
   }
   return v;

} 




std::vector<std::string>
coot::beam_in_linked_residue::make_reference_atom_names(const std::string &comp_id) const {

   std::vector<std::string> lsq_reference_atom_names;
   if (comp_id == "ASN") {
      lsq_reference_atom_names.push_back(" ND2");
      lsq_reference_atom_names.push_back(" OD1");
      lsq_reference_atom_names.push_back(" CG ");
      lsq_reference_atom_names.push_back(" CB ");
   }
   return lsq_reference_atom_names;
} 

CResidue *
coot::beam_in_linked_residue::get_residue(const std::string &comp_id,
					  CMMDBManager*mol) const {

   CResidue *r = NULL;
   int imod = 1;
   CModel *model_p = mol->GetModel(imod);
   CChain *chain_p;
   int n_chains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<n_chains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      CResidue *residue_p;
      CAtom *at;
      for (int ires=0; ires<nres; ires++) { 
	 residue_p = chain_p->GetResidue(ires);
	 std::string res_name(residue_p->GetResName());
	 if (res_name == comp_id) {
	    r = residue_p;
	    break;
	 }
      }
      if (r)
	 break;
   }

   return r;
} 
