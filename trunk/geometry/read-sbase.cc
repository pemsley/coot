
#include <stdlib.h>

#include "coot-utils.hh"
#include "protein-geometry.hh"

#ifdef USE_SBASE
void
coot::protein_geometry::read_sbase_residues() {

   if (SBase) { 
      PCSBStructure SBS;
      std::vector<std::string> local_residue_codes;
      local_residue_codes.push_back("ASN");
      local_residue_codes.push_back("TRP");

      int nStructures = SBase->GetNofStructures();
      printf ( "  Total %i structures in SBASE\n",nStructures );
      std::cout << "INFO:: Total of " << nStructures << " in SBase\n";
      
      for (int i=0; i<local_residue_codes.size(); i++) {
	 SBS = SBase->GetStructure(local_residue_codes[i].c_str());
	 if (SBS) { 
	    std::cout << SBS->compoundID << " " << SBS->Name << std::endl;
	    for (int iat=0; iat<SBS->nAtoms; iat++) {
	       CSBAtom *at = SBS->Atom[iat];
	       std::cout << "    " << at->sca_name << " ("
			 << at->x << "," << at->y << ","
			 << at->z <<")\n";
	    } 
	 } else {
	    std::cout << "WARNING:: structure " << local_residue_codes[i]
		      << " not found in SBase " << std::endl;
	 } 
      }
   } else {
      std::cout << "WARNING:: SBase not initialised"  << std::endl;
   } 
}
#endif // USE_SBASE

#ifdef USE_SBASE
CResidue *
coot::protein_geometry::get_sbase_residue(const std::string &res_name) const {

   CResidue *residue_p = NULL;
   if (SBase) { 
      CSBStructure *SBS = SBase->GetStructure(res_name.c_str());
      if (SBS) {
	 residue_p = new CResidue;
	 for (int iat=0; iat<SBS->nAtoms; iat++) {
	    CSBAtom *at = SBS->Atom[iat];
	    CAtom *new_atom = new CAtom;
	    new_atom->SetCoordinates(at->x, at->y, at->z, 1.0, 30.0);
	    std::string new_atom_name = coot::atom_id_mmdb_expand(at->sca_name, at->element);
	    // std::cout << "Adding SBase atom with name :" << new_atom_name << ":\n";
	    new_atom->SetAtomName(new_atom_name.c_str());
	    new_atom->SetElementName(at->element);
	    residue_p->AddAtom(new_atom);
	 }
      }
   }
   return residue_p;
}
#endif // USE_SBASE


#ifdef USE_SBASE
std::vector<std::string>
coot::protein_geometry::matching_sbase_residues_names(const std::string &compound_name_frag) const {

   std::string compound_name = coot::util::upcase(compound_name_frag);
   std::vector<std::string> v;
   if (SBase) {
      int nStructures = SBase->GetNofStructures();
      CFile *sf = SBase->GetStructFile();
      for (int is=0; is<nStructures; is++) {
	 CSBStructure *SBS = SBase->GetStructure(is, sf);
	 if (SBS) {
	    if (SBS->Name) { 
	       std::string sbase_residue_name(SBS->Name);
	       std::string::size_type imatch = sbase_residue_name.find(compound_name);
	       if (imatch != std::string::npos) {
		  if (0) 
		     std::cout << "DEBUG:: " << imatch << " :" << compound_name
			       << ": found in " << sbase_residue_name << std::endl;
		  std::string sbase_comp_id = SBS->compoundID;
		  v.push_back(sbase_comp_id);
	       }
	    }
	 }
      }
   }
   return v;
}
#endif // USE_SBASE


#ifdef USE_SBASE
// return mmdb sbase return codes
//
// Try to use the MONOMER_DIR_STR, ie. COOT_SBASE_DIR first, if that
// fails then use the fallback directory sbase_monomer_dir_in
// 
int
coot::protein_geometry::init_sbase(const std::string &sbase_monomer_dir_in) {
      
   int RC = SBASE_FileNotFound; // initial status.
   const char *monomer_dir = getenv(MONOMER_DIR_STR);
   
   if (!monomer_dir) {
      if (coot::is_directory_p(sbase_monomer_dir_in)) {
	 monomer_dir = sbase_monomer_dir_in.c_str();
      } else { 
	 RC = SBASE_FileNotFound; // fail
      }
   }

   if (monomer_dir) { 

      std::cout << "sbase monomer dir: " << monomer_dir << std::endl;
      SBase = new CSBase;
      RC = SBase->LoadIndex(monomer_dir);

      if (RC != SBASE_Ok) {
         std::cout << "sbase files not found in " << monomer_dir << std::endl;
	 delete SBase;
	 SBase = NULL;
      } else { 
         // std::cout << "sbase files found" << std::endl; 
      }
   }
   return RC; 
}
#endif // USE_SBASE
