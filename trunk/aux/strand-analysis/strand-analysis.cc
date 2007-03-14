
#include <sys/types.h>
#include <dirent.h>

#include <unistd.h>
#include <iostream>
#include "strand-analysis.hh"

int main(int argc, char **argv) {

   if (argc > 1) { 
      std::string dirname = argv[1];
      read_dir(dirname);
   } else {
      std::cout << "Usage: " << argv[0] << " <reference-structure-dir-name>"
		<< std::endl;
   } 
   return 0;
}

void read_dir(const std::string &dir_name) {

   coot::strands_t strands;

   std::vector<std::string> v = get_reference_pdb_list(dir_name);
   std::cout << "" << std::endl;

   for (int i=0; i<v.size(); i++) {
      strands.analyse_pdb_file(v[i]);
   }
}

std::vector<std::string>
get_reference_pdb_list(const std::string &dir_name) {

   // stat "MAPVIEW_REF_STRUCTS", opendir
   // 
   // match filenames,  read them all.
   //
   // 
   std::vector<std::string> pdb_list; 
   
   std::string ref_str_dir_str(dir_name); 
      
   DIR *ref_struct_dir = opendir(ref_str_dir_str.c_str()); 
   
   if (ref_struct_dir == NULL) { 
      std::cout << "An error occured on opendir" << std::endl;
   } else { 
      
      std::cout << ref_str_dir_str << " successfully opened" << std::endl; 
      
      // loop until the end of the filelist (readdir returns NULL)
      // 
      struct dirent *dir_ent; 
      
      while(1) { 
	 dir_ent = readdir(ref_struct_dir); 
	 
	 if (dir_ent != NULL) { 
	    
	    std::string file_str(dir_ent->d_name); 
	    if (matches_pdb_name(file_str)) {
	       
	       // Construct a file in a unix directory - non-portable?
	       // 
	       std::string s(dir_name);
	       s += "/"; 
	       s += file_str; 
	       
	       pdb_list.push_back(s); 
	    }
	 } else { 
	    break;
	 }
      }
      closedir(ref_struct_dir); 
   } 
   return pdb_list; 
} 

bool
matches_pdb_name(const std::string &file_str) { 

   short int match_flag = 0; 

   // We have to cast to (char *) here because it seems to me a bug in
   // the library.  The function strchr accepts const char *, but the
   // function that it calls FirstOccurence (and similarly
   // LastOccurence) are typed with char *.  Crapness.
   // 
   char *f_pt = strrchr( (char *) file_str.c_str(), '.');

   if (f_pt == NULL) return 0;

   if (! (strncmp(f_pt, ".pdb",    4))) { 
      match_flag = 1; 
   } else { 
      
      char *f_pt1 = strchr((char *) file_str.c_str(), '.'); 

      if (! (strncmp(f_pt1, ".pdb.gz", 7))) { 
	 match_flag = 1;

      } else { 
	 
	 // was it a pdb CD entries thing?
	 // e.g. pdb1ge0.gz
	 //
	 
	 if ( (! strncmp(f_pt, ".gz", 3)) && 
	      (! strncmp((char *) file_str.c_str(), "pdb", 3))) { 
	    match_flag = 1; 
	 }
      }
   }
   return match_flag; 
}

void
coot::strands_t::analyse_pdb_file(const std::string &filename) {

   CMMDBManager *mol = get_mol(filename);
   if (mol) { 
      CModel *model_p = mol->GetModel(1);
      int status = model_p->CalcSecStructure(1);
      if (status == SSERC_Ok) {
	 std::cout << "INFO:: SSE status was OK\n";
	 strand_analysis(model_p, mol, filename);
      } else {
	 std::cout << "INFO:: SSE status was bad\n" << status << "\n";
      }
   }
}

// Running over all of Top500:
// oxygen distance
//   Strand stats: 4.5171  +/- 0.356533
//   Helix  stats: 3.56436 +/- 0.333756
// 
void
distance_checks(CModel *model_p) {

   int imod = 1;
      
   CChain *chain_p;
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      int nres = chain_p->GetNumberOfResidues();
      PCResidue residue_p_1;
      PCResidue residue_p_2;
      for (int ires=0; ires<(nres-1); ires++) {
	 residue_p_1 = chain_p->GetResidue(ires);
	 residue_p_2 = chain_p->GetResidue(ires+1);
	 int n_atoms_1 = residue_p_1->GetNumberOfAtoms();
	 int n_atoms_2 = residue_p_2->GetNumberOfAtoms();
	 if (n_atoms_1 && n_atoms_2) {
	    if (residue_p_1->GetSeqNum() == (residue_p_2->GetSeqNum()-1)) {
	       int sse_1 = SSE(residue_p_1);
	       int sse_2 = SSE(residue_p_2);
	       if ((sse_1 == SSE_Strand) && (sse_2 == SSE_Strand) ||
		   (sse_1 == SSE_Helix)  && (sse_2 == SSE_Helix)) {
		  float d = oxygen_check(residue_p_1, residue_p_2);
		  if (d > 0) {
		     std::string ss("Helix ");
		     if (sse_1 == SSE_Strand)
			ss = "Strand";
		     std::cout << ss << " " << d << "  " 
			       << residue_p_1->GetChainID()
			       << " "
			       << residue_p_1->GetSeqNum()
			       << std::endl;
		  }
	       }
	    }
	 }
      }
   }
}

// We need to pass the mol, because we do an atom selection and the
// mol gives us the handle and access to the Select() function
// 
void
coot::strands_t::strand_analysis(CModel *model_p, CMMDBManager *mol,
			       const std::string &filename) {

   PCSheets sheets = model_p->GetSheets();
   std::cout << "has " << sheets->nSheets << " sheets" << std::endl;
   for (int isheet=0; isheet<sheets->nSheets; isheet++) {
      CSheet *sheet = sheets->Sheet[isheet];
      int nstrands = sheet->nStrands;
      std::cout << "   Sheet " << isheet << " has " << nstrands
		<< " strands " << std::endl;
      for (int istrand=0; istrand<nstrands; istrand++) {
	 CStrand *strand = sheet->Strand[istrand];
	 std::cout << "      strand " << strand->strandNo << " from "
		   << strand->initChainID << " " << strand->initSeqNum
		   << " " << strand->initICode << " to " 
		   << strand->endChainID << " " << strand->endSeqNum
		   << " " << strand->endICode << std::endl;
	 // Now do a residue selection for that strand:
	 //
	 int SelHnd = mol->NewSelection();
	 mol->Select(SelHnd, STYPE_RESIDUE, 1,
		     strand->initChainID,
		     strand->initSeqNum, strand->initICode,
		     strand->endSeqNum, strand->endICode,
		     "*", "*", 
		     "*", "*",
		     SKEY_NEW);
	 std::pair<bool, clipper::RTop_orth> ori = orient_strand_on_z(SelHnd, mol);
	 if (ori.first)
	    apply_rtop_to_strand(SelHnd, mol, ori.second);

	 add_strand(filename, mol, SelHnd);
	 
	 mol->DeleteSelection(SelHnd);
      }
   }

   std::string::size_type islash = filename.find_last_of("/");
   std::string::size_type idot   = filename.find_last_of(".");
   std::string new_pdb_name = filename.substr(islash + 1, idot - islash - 1) + 
      "-a.pdb";
//    std::cout << "debug filename bits: "
// 	     << filename << " " 
// 	     << filename.find_last_of("/") << " "
// 	     << filename.find_last_of(".") << std::endl;
   mol->WritePDBASCII((char *)new_pdb_name.c_str());
}

void
coot::strands_t::add_strand(const std::string &filename,
			    CMMDBManager *mol, int SelectionHandle) {

} 


// return the bool == 1 on good rtop.
std::pair<bool, clipper::RTop_orth>
orient_strand_on_z(int SelHnd, CMMDBManager *mol) {

   clipper::RTop_orth rtop;
   bool stat = 0;
   int nSelResidues = 0;
   PPCResidue SelResidues;
   mol->GetSelIndex(SelHnd, SelResidues, nSelResidues);
   std::cout << "      selected " << nSelResidues << " residues "
	     << std::endl;
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
   if (atom_vec.size() != 3* nSelResidues) {
      std::cout << "skipping this strange strand " << atom_vec.size() << " "
		<< nSelResidues << std::endl;
   } else {
      // find the rtop from this clipper method
      rtop = clipper::RTop_orth(atom_vec, z_points);
      stat = 1;
   }
   return std::pair<bool, clipper::RTop_orth> (stat, rtop);
}

// fiddle with mol
void apply_rtop_to_strand(int SelHnd, CMMDBManager *mol,
			  const clipper::RTop_orth &rtop) {

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


std::vector<clipper::Coord_orth> z_control_points(int nres) {

   std::vector<clipper::Coord_orth> r;
   double stl = 3.8 * 0.866; // the angle subtended on the z-axis of
			     // the CA-CA distance.
   double curve_len = 0.2; // how much the strand curves per residues:
                           // used to set the control points.

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


int SSE(CResidue *res) {
   return res->SSE;
} 

CMMDBManager *get_mol(const std::string &filename) {

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



// return -1 on failure
// 
float oxygen_check(CResidue *residue_p_1, CResidue *residue_p_2) {

   int n_atoms_1 = residue_p_1->GetNumberOfAtoms();
   int n_atoms_2 = residue_p_2->GetNumberOfAtoms();
   CAtom *at, *at2;
   float d = -1;
   short int done_it = 0; 
   
   for (int iat=0; iat<n_atoms_1; iat++) {
      at = residue_p_1->GetAtom(iat);
      std::string atom_name (at->name);
      if (atom_name == " O  ") {
	 for (int iat2=0; iat2<n_atoms_2; iat2++) {
	    at2 = residue_p_2->GetAtom(iat2);
	    std::string atom_name_2 (at2->name);
	    if (atom_name_2 == " O  ") {
	       std::string altconf1 = at->altLoc;
	       std::string altconf2 = at2->altLoc;
	       if ((altconf1 == "") && (altconf2 == "")) { 
		  clipper::Coord_orth c1(at->x,  at->y,  at->z);
		  clipper::Coord_orth c2(at2->x, at2->y, at2->z);
		  d = clipper::Coord_orth::length(c1, c2);
		  done_it = 1;
		  break;
	       }
	    }
	 }
      }
      if (done_it)
	 break;
   }
   return d;
}
