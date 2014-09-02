
#include <sys/types.h>
#include <dirent.h>
#include <string.h>
#include <stdlib.h>

#include <unistd.h>
#include <iostream>
#include "strand-analysis.hh"

int main(int argc, char **argv) {

   if (argc > 1) { 
      std::string dirname = argv[1];
      int strand_length = 9; 
      if (argc > 2)
	 strand_length = atoi(argv[2]);
      read_dir(dirname, strand_length);
   } else {
      std::cout << "Usage: " << argv[0] << " <reference-structure-dir-name>"
		<< std::endl;
   } 
   return 0;
}

void read_dir(const std::string &dir_name, int strand_length) {

   coot::strands_t strands;

   std::vector<std::string> v = get_reference_pdb_list(dir_name);
   std::cout << "" << std::endl;

   for (int i=0; i<v.size(); i++) {
      strands.analyse_pdb_file(v[i]);
   }

   strands.post_read_analysis(strand_length);
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
	       if (((sse_1 == SSE_Strand) && (sse_2 == SSE_Strand)) ||
		   ((sse_1 == SSE_Helix)  && (sse_2 == SSE_Helix))) {
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
	    std::pair<bool, clipper::RTop_orth> ori = orient_strand_on_z(SelHnd, mol);
	    if (ori.first)
	       apply_rtop_to_strand(SelHnd, mol, ori.second);
	    
	    add_strand(filename, mol, strand, SelHnd);
	    
	    // mol->DeleteSelection(SelHnd);
	 }
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
			    CMMDBManager *mol,
			    CStrand *strand_p,
			    int SelectionHandle) {

   int nSelResidues = 0;
   PPCResidue SelResidues;
   mol->GetSelIndex(SelectionHandle, SelResidues, nSelResidues);

   coot::strand_info_t s;
   s.length = nSelResidues;
   s.strand = *strand_p;
   std::string name = filename;
   s.name = filename;
   std::string::size_type islash = filename.find_last_of("/");
   if (islash == std::string::npos) {
      name = filename;
   } else {
      name = filename.substr(islash+1, filename.length());
   }
   std::string::size_type ipdb = name.rfind(".pdb");
   if (ipdb == std::string::npos) { 
      // name == name; 
   } else { 
      name = name.substr(0,ipdb);
   } 
   s.name = name;
   s.mol = mol;
   s.SelectionHandle = SelectionHandle;
   strand_infos.push_back(s);
   
} 


// return the bool == 1 on good rtop.
std::pair<bool, clipper::RTop_orth>
orient_strand_on_z(int SelHnd, CMMDBManager *mol) {

   clipper::RTop_orth rtop = clipper::RTop_orth::identity();
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
   if (atom_vec.size() != 3* nSelResidues) {
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

void
coot::strands_t::post_read_analysis(int this_length) const {

   std::cout << "Captured " << strand_infos.size() << " strands" << std::endl;

   std::vector<int> strand_lengths(20, 0);
   
   for (unsigned int i=0; i<strand_infos.size(); i++) {
      int l = strand_infos[i].length;
      if (l >= 0)
	 if (l < 20)
	    strand_lengths[l]++;
   }

   for (unsigned int i=0; i<20; i++) {
      std::cout << "strand_length: " << i << "  " << strand_lengths[i] << std::endl;
   }

   int n_strands = strand_lengths[this_length];
   std::cout << "There are " << n_strands << " strands of length " << this_length << std::endl;
   std::vector<std::pair<std::string,std::vector<float> > > dist_arr(n_strands, std::pair<std::string, std::vector<float> > ("name", std::vector<float>(n_strands, 0.0)));
   int iarr_index = 0; 
   for (unsigned int istrand=0; istrand<strand_infos.size(); istrand++) {
      int l1 = strand_infos[istrand].length;
      if (l1 == this_length) {
	 // strand_names[iarr_index] = strand_infos[istrand].name;
	 int jarr_index = 0; 
	 for (unsigned int jstrand=0; jstrand<strand_infos.size(); jstrand++) {
	    int l2 = strand_infos[jstrand].length;
	    if (l2 == this_length) {
	       float v = residual_between_strands(istrand, jstrand, this_length);
// 	       std::cout << "setting dist_arr[" << iarr_index << "][" << jarr_index
// 			 << "]" << std::endl;
	       dist_arr[iarr_index].second[jarr_index] = v;
	       dist_arr[iarr_index].first = strand_infos[istrand].name;
	       jarr_index++;
	    }
	 }
	 iarr_index++; 
      }
   }

   std::vector<std::pair< std::string,std::vector<float> > > f_dist_arr = filter_bad_fits(dist_arr); 
   
   // istrand is indexing dist_arr now
   //
   if (0) { 
      for (int istrand=0; istrand<n_strands; istrand++) {
	 for (int jstrand=0; jstrand<n_strands; jstrand++) {
	    std::cout << istrand << " " << jstrand << " "
		      << dist_arr[istrand].second[jstrand] << std::endl;
	 }
      }
   }

   if (0) {
      std::cout << " " << n_strands << std::endl;
      for (int istrand=0; istrand<n_strands; istrand++) {
	 // the front of the line
	 std::cout << "strand" << istrand << "    ";
	 // all the row
	 std::cout.setf(std::ios::fixed);
	 std::cout.precision(7);
	 for (int jstrand=0; jstrand<n_strands; jstrand++) {
	    std::cout << dist_arr[istrand].second[jstrand] << " ";
	 }
	 std::cout << std::endl;
      }
   }

   if (1) {
      int f_n_strands = f_dist_arr.size();
      printf("   %d\n", f_n_strands);
      for (int istrand=0; istrand<f_n_strands; istrand++) {

	 // the front of the line
	 std::string fl = f_dist_arr[istrand].first;
	 if (fl.length() < 6)
	    fl += " ";
	 printf("%s     ", fl.c_str());

	 // all the row
	 for (int jstrand=0; jstrand<n_strands; jstrand++) {
	    float v = f_dist_arr[istrand].second[jstrand];
	    if (v < 0.0) 
	       printf("wierd! %f\n", v);
	    else
	       if (v > 2000) 
		  printf("wierd! %f\n", v);
	    
	    printf("%7.4f ", v);
	 }
	 printf("\n");
      }
   } 

} 


std::vector<std::pair<std::string, std::vector<float> > >
coot::strands_t::filter_bad_fits(const std::vector<std::pair<std::string, std::vector<float> > > &dist_arr) const {
   
   std::vector<std::pair<std::string, std::vector<float> > > r = dist_arr;

   // first get the sums of all the rows
   std::vector<std::pair<float,float> > stats = get_stats(dist_arr);
   float sum = 0.0;
   float sum_sq = 0.0;
   for (int i=0; i<stats.size(); i++) {
      float v = stats[i].first;
      sum += v;
      sum_sq += v*v;
   }
   float mean_of_means = sum/float(stats.size());
   float var_of_means = sum_sq/float(stats.size()) - mean_of_means*mean_of_means;
   float std_dev = -1; // undefined
   if (var_of_means >= 0)
      std_dev = sqrt(var_of_means);
   else
      std::cout << "Ooops can't calculate std dev from variance: " << var_of_means << std::endl;
   for (int i=0; i<stats.size(); i++) {
      std::cout << "column " << i << " " << stats[i].first << "\n";
   }

   float threash = mean_of_means + 4*std_dev;
   threash = 10.0; 
   std::cout << "   mean= " << mean_of_means << " std_dev: " << std_dev
	     << " threash= " << threash << std::endl;
   int bad_column = -1; // unset
   for (int i=0; i<stats.size(); i++) {
      if (stats[i].first > threash) { 
	 std::cout << "column " << i << " is an outlier: v= " << stats[i].first
		   << "   mean= " << mean_of_means << " std_dev: " << std_dev << std::endl;
	 threash = stats[i].first;
	 bad_column = i;
      }
   }

   // recur if we cleaned out a column
   if (bad_column > -1) {
      return filter_bad_fits(clear_out_column(bad_column, dist_arr));
   } else { 
      return r;
   }
}

std::vector<std::pair<std::string, std::vector<float> > >
coot::strands_t::clear_out_column(int bad_col, const std::vector<std::pair<std::string, std::vector<float> > > &d) const {

   std::vector<std::pair<std::string, std::vector<float> > > r1;
   std::vector<std::pair<std::string, std::vector<float> > > r2;
   for (unsigned int i=0; i<d.size(); i++) {
      if (i != bad_col)
	 r1.push_back(d[i]);
      else
	 std::cout << "DEBUG:: ommmited col " << bad_col << std::endl;
   }
   
   for (unsigned int i=0; i<r1.size(); i++) {
      std::pair<std::string, std::vector<float> > p(r1[i].first, std::vector<float>(0));
      for (unsigned int j=0; j<r1[i].second.size(); j++) {
	 if (j != bad_col) {
	    float v = r1[i].second[j];
	    if (v < 0.0) {
	       std::cout << "replacing neg  weird value " << v << " at " << i << " " << j << " with 0"
			 << std::endl;
	       p.second.push_back(0.0);
	    } else {
	       if (v > 2000.0) {
		  std::cout << "replacing high weird value " << v << " at " << i << " " << j << " with 8"
			    << std::endl;
		  p.second.push_back(8.0);
	       } else {
		  p.second.push_back(v); // normal case
	       } 
	    }
	 } else {
	    std::cout << "DEBUG:: omitted row " << bad_col << " with value "
		      << r1[i].second[j] << std::endl;
	 }
      }
      r2.push_back(p);
   }
   return r2;
}

// Return a vector of sums and standard deviations
std::vector<std::pair<float, float> >
coot::strands_t::get_stats(const std::vector<std::pair<std::string, std::vector<float> > > &dist_arr) const {

   std::vector<std::pair<float, float> > sums(dist_arr.size(), std::pair<float, float>(0.0, 0.0));
   for (unsigned int i=0; i<dist_arr.size(); i++) {
      for (unsigned int j=0; j<dist_arr[i].second.size(); j++) {
	 float v = dist_arr[i].second[j];
	 sums[i].first += v;
	 sums[i].second += v*v;
      }
   }
   for (unsigned int i=0; i<sums.size(); i++) {
      sums[i].first  /= float(sums.size());
      sums[i].second /= float(sums.size()) - sums[i].first*sums[i].first;
   }
   return sums; // now mean and var
}


float
coot::strands_t::residual_between_strands(int istrand, int jstrand, int strand_length) const {

   float r = -1; // bad value
   
   // get first set of atoms
   int nSelResidues1 = 0;
   PPCResidue SelResidues1;
   CMMDBManager *mol1 = strand_infos[istrand].mol;
   mol1->GetSelIndex(strand_infos[istrand].SelectionHandle, SelResidues1, nSelResidues1);
   int nSelResidues2 = 0;
   PPCResidue SelResidues2;
   CMMDBManager *mol2 = strand_infos[jstrand].mol;
   mol2->GetSelIndex(strand_infos[jstrand].SelectionHandle, SelResidues2, nSelResidues2);
   std::vector<clipper::Coord_orth> found_atoms_strand_1;
   std::vector<clipper::Coord_orth> found_atoms_strand_2;

   if (nSelResidues1 == nSelResidues2) {
      // which they should do!
      float sum_dist_sq = 0;
      int n_res_with_atoms = 0; 
      for (int ires=0; ires<nSelResidues1; ires++) {

	 int n_atoms_1;
	 PPCAtom residue_atoms_1;
	 int n_atoms_2;
	 PPCAtom residue_atoms_2;
	 
	 SelResidues1[ires]->GetAtomTable(residue_atoms_1, n_atoms_1);
	 SelResidues2[ires]->GetAtomTable(residue_atoms_2, n_atoms_2);
	 
	 std::pair<CAtom *, CAtom *> o(0,0);
	 std::pair<CAtom *, CAtom *> c(0,0);
	 std::pair<CAtom *, CAtom *> ca(0,0);
	 std::pair<CAtom *, CAtom *> n(0,0);

	 for (int iatom1=0; iatom1<n_atoms_1; iatom1++) {
	    std::string at_name1 = residue_atoms_1[iatom1]->name;
	    if (at_name1 == " CA ") {
	       ca.first = residue_atoms_1[iatom1];
	    }
	    if (at_name1 == " O  ") {
	       o.first = residue_atoms_1[iatom1];
	    }
	    if (at_name1 == " C  ") {
	       c.first = residue_atoms_1[iatom1];
	    }
	    if (at_name1 == " N  ") {
	       n.first = residue_atoms_1[iatom1];
	    }
	 }

	 for (int iatom2=0; iatom2<n_atoms_2; iatom2++) {
	    std::string at_name2 = residue_atoms_2[iatom2]->name;
	    if (at_name2 == " CA ") {
	       ca.second = residue_atoms_2[iatom2];
	    }
	    if (at_name2 == " O  ") {
	       o.second = residue_atoms_2[iatom2];
	    }
	    if (at_name2 == " C  ") {
	       c.second = residue_atoms_2[iatom2];
	    }
	    if (at_name2 == " N  ") {
	       n.second = residue_atoms_2[iatom2];
	    }
	 }

	 // OK are all the atoms filled?
	 if (ca.first && ca.second) { 
	    if (c.first && c.second) { 
	       if (n.first && n.second) { 
		  if (o.first && o.second) {

		     clipper::Coord_orth cat1;
		     clipper::Coord_orth cat2;
		     
		     cat1 = clipper::Coord_orth(ca.first->x,  ca.first->y,  ca.first->z);
		     cat2 = clipper::Coord_orth(ca.second->x, ca.second->y, ca.second->z);
		     found_atoms_strand_1.push_back(cat1);
		     found_atoms_strand_2.push_back(cat2);

		     cat1 = clipper::Coord_orth(c.first->x,  c.first->y,  c.first->z);
		     cat2 = clipper::Coord_orth(c.second->x, c.second->y, c.second->z);
		     found_atoms_strand_1.push_back(cat1);
		     found_atoms_strand_2.push_back(cat2);

		     cat1 = clipper::Coord_orth(o.first->x,  o.first->y,  o.first->z);
		     cat2 = clipper::Coord_orth(o.second->x, o.second->y, o.second->z);
		     found_atoms_strand_1.push_back(cat1);
		     found_atoms_strand_2.push_back(cat2);

		     cat1 = clipper::Coord_orth(n.first->x,  n.first->y,  n.first->z);
		     cat2 = clipper::Coord_orth(n.second->x, n.second->y, n.second->z);
		     found_atoms_strand_1.push_back(cat1);
		     found_atoms_strand_2.push_back(cat2);


		     

		     float x2 = ca.first->x - ca.second->x;
		     float y2 = ca.first->y - ca.second->y;
		     float z2 = ca.first->z - ca.second->z;
		     sum_dist_sq += x2*x2 + y2*y2 + z2*z2;
		     n_res_with_atoms++; 
		  }
	       }
	    }
	 }
      }
      float var = -1;
      if (n_res_with_atoms > 0) { 
	 float var = sum_dist_sq/float(4*n_res_with_atoms);
	 r = var;
      }
      coot::stats_t s =
	 get_rtop_and_apply(found_atoms_strand_1, found_atoms_strand_2, SelResidues2, nSelResidues2);
      r = s.mean;
      
   } else {
      std::cout << "ERROR:: Something bad in residual_between_strands " << istrand << " "
		<< jstrand << std::endl;
   } 
   return r;
}


coot::stats_t
coot::get_rtop_and_apply(const std::vector<clipper::Coord_orth> &found_atoms_strand_1,
			 const std::vector<clipper::Coord_orth> &found_atoms_strand_2,
			 PPCResidue SelResidues2, int nSelResidues2) {
   coot::stats_t s(-1,0);

   clipper::RTop_orth matching_rtop(found_atoms_strand_2, found_atoms_strand_1);
   float sum_d = 0.0; 
   float sum_d2 = 0.0;
   float sum_do = 0.0; 
   float sum_d2o = 0.0;
   
   if (found_atoms_strand_1.size() != found_atoms_strand_2.size()) {
      std::cout << " Not match on number of atoms in get_rtop_and_apply()!\n";
   } else { 
      if (found_atoms_strand_1.size() == 0) {
	 std::cout << "Empty coord vectors in get_rtop_and_apply()!\n";
      } else {
	 for (unsigned int i=0; i<found_atoms_strand_1.size(); i++) {
	    float d = clipper::Coord_orth::length(found_atoms_strand_1[i], found_atoms_strand_2[i].transform(matching_rtop));
	    float dor = clipper::Coord_orth::length(found_atoms_strand_1[i],found_atoms_strand_2[i]);
	    sum_d += d;
	    sum_do += dor;
	    sum_d2 += d*d;
	    sum_d2o += dor*dor;
	 }
	 float mean = sum_d/float(found_atoms_strand_1.size());
	 float var = sum_d2/float(found_atoms_strand_1.size()) - mean*mean;
	 // std::cout << "Distances pre-match:"
	 // << sum_do/float(found_atoms_strand_1.size())
	 // << "  post-match: " << mean << std::endl;
	 s = coot::stats_t(mean, var);
      }
   }
   return s;
}
