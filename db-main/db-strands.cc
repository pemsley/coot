/* db-main/db-strand.cc
 * 
 * Copyright 2007 The University of York
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
 * 02110-1301, USA
 */


// for directory reading
#include <sys/types.h>   
#include <dirent.h>

#include "utils/coot-utils.hh"
#include "coot-utils/coot-coord-utils.hh"

#include "db-main.hh" // for matches_pdb_name
#include "db-strands.hh"

coot::db_strands::db_strands() {
   
   const char *ds = getenv("COOT_REF_SEC_STRUCTS");
   if (!ds) {
      const char *d = getenv("COOT_REF_STRUCTS");
      if (d)
	 ref_str_dir_str = d;
      else {

	 // 20080403 recactivate this code
	 // fall back to the db-main reference structures
	 std::string d1(PKGDATADIR); // $prefix/share/coot
	 // std::string d2 = coot::util::append_dir_dir(d1, "coot"); // no need for this
	 std::string d3 = coot::util::append_dir_dir(d1, "reference-structures");
	 ref_str_dir_str = d3;

	 // Before 20080403
// 	 ref_str_dir_str = "Undefined";
// 	 std::cout << "WARNING:: reference structure directory undefined. COOT_REF_STRUCTS currently "
// 		   << "\n"
// 		   << "WARNING:: won't do\n";
      }
   } else {
      ref_str_dir_str = ds;
   } 
}


std::vector<coot::minimol::molecule>
coot::db_strands::get_reference_strands(int n_strands, int strand_length) {

   std::vector<coot::minimol::molecule> mv;
   std::vector<std::string> v = get_reference_pdb_list();

   for (unsigned int ipdb=0; ipdb<v.size(); ipdb++) {
      if (mv.size() >= n_strands) {
	 break;
      } else { 
	 std::string filename = v[ipdb];
	 CMMDBManager *mol = get_mol(filename);
	 if (mol) { 
	    CModel *model_p = mol->GetModel(1);
#ifdef HAVE_MMDB_WITH_CISPEP	    
	    int status = model_p->CalcSecStructure(1); // Hmm. Used to have an atomselhnd arg.
#else 	    
	    int status = model_p->CalcSecStructure(1);
	    // int status = SSERC_Ok;
#endif // HAVE_MMDB_WITH_CISPEP
	    if (status == SSERC_Ok) {
	       std::cout << "INFO:: SSE status was OK\n";
	       std::vector<coot::minimol::molecule> v_strand = 
		  strand_analysis(model_p, mol, filename, strand_length);
	       if (v_strand.size() > 0) {
		  for (unsigned int iv=0; iv<v_strand.size(); iv++) {
		     if (int(mv.size()) < n_strands) {
			mv.push_back(v_strand[iv]);
		  }
		  }
	       }
	    } else {
	       std::cout << "INFO:: SSE status was bad\n" << status << "\n";
	    }
	 }
      }
   } 
   return mv;
} 

CMMDBManager *
coot::db_strands::get_mol(const std::string &filename) const {

   CMMDBManager *MMDBManager = new CMMDBManager();
   int err = MMDBManager->ReadCoorFile((char *)filename.c_str());
   if (err) {
      std::cout << "Error reading " << filename << std::endl;
      delete MMDBManager;
      MMDBManager = 0;
   } else {
      std::cout << "Read OK: " << filename << std::endl;
   } 
   return MMDBManager;
} 

std::vector<std::string>
coot::db_strands::get_reference_pdb_list() const {

   std::vector<std::string> pdb_list;
   
   DIR *ref_struct_dir = opendir(ref_str_dir_str.c_str()); 
   
   if (ref_struct_dir == NULL) { 
      std::cout << "An error occured on opendir on "
		<< ref_str_dir_str << std::endl;
   } else { 
      
      std::cout << ref_str_dir_str << " successfully opened" << std::endl; 
      
      // loop until the end of the filelist (readdir returns NULL)
      // 
      struct dirent *dir_ent; 
      
      while(1) {
	 dir_ent = readdir(ref_struct_dir); 
	 
	 if (dir_ent != NULL) { 
	    
	    std::string file_str(dir_ent->d_name); 
	    if (coot::matches_pdb_name(file_str)) {
	       
	       // Construct a file in a unix directory - non-portable?
	       // 
	       std::string s(ref_str_dir_str);
	       s += "/"; 
	       s += file_str; 
	       
	       pdb_list.push_back(s); 
	    } else {
	       std::cout << "DEBUG:: " << file_str << " fails pdb extension test"
			 << std::endl;
	    } 
	 } else { 
	    break;
	 }
      }
      closedir(ref_struct_dir); 
   }

   std::cout << "INFO:: found " << pdb_list.size() << " PDB files in reference"
	     << " structure direcotry" << std::endl;
   return pdb_list; 
} 

// We need to pass the mol, because we do an atom selection and the
// mol gives us the handle and access to the Select() function
// 
std::vector<coot::minimol::molecule>
coot::db_strands::strand_analysis(CModel *model_p, CMMDBManager *mol,
				  const std::string &filename, int strand_length) const {

   std::vector<coot::minimol::molecule> rv;

   PCSheets sheets = model_p->GetSheets();
   std::cout << "has " << sheets->nSheets << " sheets" << std::endl;
   for (int isheet=0; isheet<sheets->nSheets; isheet++) {
      CSheet *sheet = sheets->Sheet[isheet];
      int nstrands = sheet->nStrands;
      std::cout << "   Sheet " << isheet << " has " << nstrands
		<< " strands " << std::endl;
      for (int istrand=0; istrand<nstrands; istrand++) {
	 CStrand *strand = sheet->Strand[istrand];
	 if (strand) { 
	    std::cout << "      strand " << strand->strandNo << " from "
		      << strand->initChainID << " " << strand->initSeqNum
		      << " " << strand->initICode << " to " 
		      << strand->endChainID << " " << strand->endSeqNum
		      << " " << strand->endICode << std::endl;

	    // Now do a residue selection for that strand.  Once only
	    // and pass around the handle, from which we do a
	    // GetSelIndex().
	    //
	    int SelHnd = mol->NewSelection();
	    mol->Select(SelHnd, STYPE_RESIDUE, 1,
			strand->initChainID,
			strand->initSeqNum, strand->initICode,
			strand->endSeqNum, strand->endICode,
			"*", "*", 
			"*", "*",
			SKEY_NEW);
	    PPCResidue SelResidues = 0;
	    int nSelResidues;
	    mol->GetSelIndex(SelHnd, SelResidues, nSelResidues);
	    if (nSelResidues == strand_length) { 
	       
	       std::pair<bool, clipper::RTop_orth> ori = orient_strand_on_z(SelHnd, mol);
	       if (ori.first) { 
		  apply_rtop_to_strand(SelHnd, mol, ori.second);
		  std::pair<CMMDBManager *, int> p =
		     coot::util::create_mmdbmanager_from_res_selection(mol, SelResidues, nSelResidues,
								       0, 0, "", strand->initChainID, 0);
		  if (p.second) {
		     trim_to_mainchain(p.first);
		     coot::minimol::molecule m(p.first);
		     rv.push_back(m);
		  }
	       }
	    }
	    mol->DeleteSelection(SelHnd);
	 }
      }
   }
   return rv;
}

// And remove hydrogens too
void
coot::db_strands::trim_to_mainchain(CMMDBManager *mol) const {

   int imod = 1;
      
   CModel *model_p = mol->GetModel(imod);
   CChain *chain_p;
   // run over chains of the existing mol
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      PCResidue residue_p;
      CAtom *at;
      for (int ires=0; ires<nres; ires++) { 
	 residue_p = chain_p->GetResidue(ires);
	 int n_atoms = residue_p->GetNumberOfAtoms();
	 
	 for (int iat=0; iat<n_atoms; iat++) {
	    at = residue_p->GetAtom(iat);
	    std::string ele = at->element;
	    if (!is_main_chain_or_cb_p(at)
		|| ele == " H" || ele == " D") {
	       residue_p->DeleteAtom(iat);
	    }
	 }
	 residue_p->TrimAtomTable();
      }
   }
   mol->PDBCleanup(PDBCLEAN_SERIAL|PDBCLEAN_INDEX);
   mol->FinishStructEdit();
}

// return the bool == 1 on good rtop.
std::pair<bool, clipper::RTop_orth>
coot::db_strands::orient_strand_on_z(int SelHnd, CMMDBManager *mol) const {

   clipper::Mat33<double> m_dum(1,0,0,0,1,0,0,0,1);
   clipper::Coord_orth pt_dum(0,0,0);
   clipper::RTop_orth rtop(m_dum, pt_dum);
   bool stat = 0;
   int nSelResidues = 0;
   PPCResidue SelResidues;
   mol->GetSelIndex(SelHnd, SelResidues, nSelResidues);
   //    std::cout << "      selected " << nSelResidues << " residues "
   // << std::endl;
   std::vector<clipper::Coord_orth> z_points = z_control_points(nSelResidues);
   std::vector<clipper::Coord_orth> atom_vec;

   for (int ires=0; ires<nSelResidues; ires++) {
      CResidue *residue_p = SelResidues[ires];
      int n_atoms;
      PPCAtom atoms;
      residue_p->GetAtomTable(atoms, n_atoms);
      for (int iat=0; iat<n_atoms; iat++) {
	 CAtom *at = atoms[iat];
	 std::string atom_name(at->name);
	 if (atom_name == " N  ") {
	    clipper::Coord_orth pt(at->x, at->y, at->z);
	    atom_vec.push_back(pt);
	 } 
	 if (atom_name == " CA ") {
	    clipper::Coord_orth pt(at->x, at->y, at->z);
	    atom_vec.push_back(pt);
	 }
	 if (atom_name == " C  ") {
	    clipper::Coord_orth pt(at->x, at->y, at->z);
	    atom_vec.push_back(pt);
	 } 
      } 
   }
   if (int(atom_vec.size()) != 3* nSelResidues) {
      std::cout << "skipping this strange strand with " << atom_vec.size()
		<< " atoms (should be " << 3*nSelResidues << ") and "
		<< nSelResidues << " residues" << std::endl;
   } else {
      // find the rtop from this clipper method
      rtop = clipper::RTop_orth(atom_vec, z_points);
      stat = 1;
   }
   return std::pair<bool, clipper::RTop_orth> (stat, rtop);
}


// fiddle with mol
void
coot::db_strands::apply_rtop_to_strand(int SelHnd, CMMDBManager *mol,
				       const clipper::RTop_orth &rtop) const {

   int nSelResidues = 0;
   PPCResidue SelResidues;
   mol->GetSelIndex(SelHnd, SelResidues, nSelResidues);
   for (int ires=0; ires<nSelResidues; ires++) {
      CResidue *residue_p = SelResidues[ires];
      int n_atoms;
      PPCAtom atoms;
      residue_p->GetAtomTable(atoms, n_atoms);
      for (int iat=0; iat<n_atoms; iat++) {
	 CAtom *at = atoms[iat];
	 clipper::Coord_orth pt(at->x, at->y, at->z);
	 clipper::Coord_orth n = pt.transform(rtop);
	 at->x = n.x();
	 at->y = n.y();
	 at->z = n.z();
      }
   }
}


std::vector<clipper::Coord_orth>
coot::db_strands::z_control_points(int nres) const {

   std::vector<clipper::Coord_orth> r;
   double stl = 3.8 * 0.866; // the angle subtended on the z-axis of
			     // the CA-CA distance.
//    double curve_len = 0.2; // how much the strand curves per residues:
//                            // used to set the control points.

   double base = -double(nres-1)/2.0 * stl;
   // double curve_base = double(nres-1)/2.0 * curve_len; 

   for (int ires=0; ires<nres; ires++) {
      // Add control points for N CA and C
      double d = double(ires) * stl + base;
      // double cb = curve_base = double(ires)*0.2 + curve_base;
      double cb = 0;
      clipper::Coord_orth ca( cb,       0, d);
      clipper::Coord_orth  n( cb + 0.4, 0, d-1.1);
      clipper::Coord_orth  c( cb + 0.4, 0, d+1.1);
      r.push_back(ca);
      r.push_back(n);
      r.push_back(c);
      if (0) { 
	 clipper::String s1("(to-generic-object-add-point 0 \"blue\" 5 ");
	 s1 += clipper::String(n.x());
	 s1 += " ";
	 s1 += clipper::String(n.y());
	 s1 += " ";
	 s1 += clipper::String(n.z());
	 s1 += ") ";
	 std::cout << s1 << std::endl;
	 clipper::String s2("(to-generic-object-add-point 0 \"green\" 5 ");
	 s2 += clipper::String(ca.x());
	 s2 += " ";
	 s2 += clipper::String(ca.y());
	 s2 += " ";
	 s2 += clipper::String(ca.z());
	 s2 += ") ";
	 std::cout << s2 << std::endl;
	 clipper::String s3("(to-generic-object-add-point 0 \"red\" 5 ");
	 s3 += clipper::String(c.x());
	 s3 += " ";
	 s3 += clipper::String(c.y());
	 s3 += " ";
	 s3 += clipper::String(c.z());
	 s3 += ") ";
	 std::cout << s3 << std::endl;
      }
   }
   return r;
}
