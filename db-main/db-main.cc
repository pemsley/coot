/* db-main/db-main.cc
 * 
 * Copyright 2002, 2003, 2006 The University of York
 * Copyright 2008 The University of York
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

#include <sys/types.h> // for stating
#include <sys/stat.h>
#include <string.h>

#include <algorithm>

#if !defined _MSC_VER
#include <unistd.h>
#else
#define DATADIR "C:/coot/share"
#endif

#include "utils/coot-utils.hh"
#include "db-main.hh"

coot::main_fragment_t::main_fragment_t(int i_start_res,
				       int mol_no,
				       std::vector<float> eigns_sqrt,
				       std::string seg_id,
				       std::pair<bool, clipper::Coord_orth> mid_o_pt,
				       int ilen) {
   
   i_start_res_ = i_start_res; 
   molecule_number = mol_no; 
   sqrt_eigen_values = eigns_sqrt; 
   segment_id = seg_id;
   ilength = ilen;
   middle_carbonyl_oxygen_position = mid_o_pt;
}

int 
coot::db_main::fill_with_fragments(int ilength) {

   std::vector<std::string> pdb_filename_list = get_reference_pdb_list();

   // we need molecule_number index to correspond to the
   // pdb_filename_list index (because i is stored in
   // mainchain_frag_db
   
   molecule_list.resize(pdb_filename_list.size());

   for  (unsigned int i=0; i< pdb_filename_list.size(); i++) {
      std::cout << "     Adding fragments from: " << pdb_filename_list[i] << std::endl;
      coot::minimol::molecule m;
      m.read_file(pdb_filename_list[i]);

      if (m.fragments.size() > 0) {

	 // we need molecule_number index to correspond to the
	 // pdb_filename_list index.
	 molecule_list[i] = m;
                        	     // molecule_list is used to
				     // retrieve the Calphas if we get
				     // a close match (by eigens) and
				     // other mainchain atoms if the
				     // overlap is close.
	 int istart = m[0].min_res_no();
	 int nfragment = 0; // count for this pdb (info purposes only)

	 while ((istart+ilength) < m[0].max_residue_number()) {

	    coot::minimol::fragment fragment = get_fragment_ca_atoms(istart, ilength, m);

	    // was that an acceptable fragment? (when we get to the end
	    // of the protein chain, it is not).

	    //  min: 2 max: 7:  has range-length: 6:
	    if ( (int(fragment.n_filled_residues()) == ilength)) {

	       // What is the position of the oxygen in the 2th (0 indexed) residue?
	       
	       // now get the eigen values/functions for that fragment:
	       clipper::Matrix<float> mat = make_cov_matrix(frag_to_coords(fragment));
	       std::vector<float> eigens = mat.eigen( true );
	       for(unsigned int j=0; j<eigens.size(); j++)
		  eigens[j] = sqrt(eigens[j]);
	       std::pair<bool, clipper::Coord_orth> mid_o_pt = get_middle_ox_pos(fragment);
	       mainchain_frag_db.push_back(main_fragment_t(istart,i,eigens,
							   fragment.fragment_id,
							   mid_o_pt, ilength));
	       istart++;
	       nfragment++;
	    } else {
	       break;
	    }
	 }
	 std::cout << nfragment << " fragments garnered from "
		   << pdb_filename_list[i] << std::endl;
      }
   }

   int n5 = 0;
   for (unsigned int i=0; i<mainchain_frag_db.size(); i++) {
      if (mainchain_frag_db[i].ilength == 5)
	 n5++;
   }
   
   std::cout << "INFO:: " << mainchain_frag_db.size() << " fragments found in total";
   std::cout << "       of which " << n5 << " were 5 peptides " << std::endl;
   return mainchain_frag_db.size();
}

std::pair<bool, clipper::Coord_orth>
coot::db_main::get_middle_ox_pos(const coot::minimol::fragment &fragment) const {

   bool have = 0;
   clipper::Coord_orth pt(0,0,0); // unset

   int residue_count = 0;
   for (int ires=fragment.min_res_no(); ires<=fragment.max_residue_number(); ires++) {
      if (residue_count == 2) {
	 if (fragment[ires].atoms.size() > 0) {
	    for (unsigned int iatom=0; iatom<fragment[ires].atoms.size(); iatom++) {
	       if (fragment[ires][iatom].name == " O  ") {
		  if (fragment[ires][iatom].altLoc == "") {
		     have = 1;
		     pt = fragment[ires][iatom].pos;
		  }
	       }
	    }
	 }
      }
      residue_count++;
   }
   return std::pair<bool, clipper::Coord_orth> (have, pt);
}


std::vector<clipper::Coord_orth>
coot::db_main::frag_to_coords(const coot::minimol::fragment &fragment) const {

   std::vector<clipper::Coord_orth> a;
   std::vector<coot::minimol::atom *> c = fragment.select_atoms_serial();

   for (unsigned int i=0; i< c.size(); i++) {
      a.push_back(c[i]->pos); 
   }
   return a; 
} 

// So, we have a molecule, we want to chop it into fragments,
// starting at residue istart.  We don't want to keep the actual
// fragments, we will just store just the eigen values of the
// fragment that we pull out (and the discard), along with
// enough other information to get the fragment back if/when we
// want to use it in recombination.
// 
coot::minimol::fragment
coot::db_main::get_fragment_ca_atoms(int istart, int ilength,
				     const coot::minimol::molecule &m) const {

   coot::minimol::fragment f;
   coot::minimol::residue residoo;
   f.residues.push_back(residoo); // the 0th (unused residue)

   if (m[0].max_residue_number() >= istart+ilength) { 
      for (int ires=istart; ires<istart+ilength; ires++) {
	 for (unsigned int iat=0; iat<m[0][ires].atoms.size(); iat++) {
	    if (m[0][ires][iat].name == " CA ") {
	       if (m[0][ires][iat].altLoc == "") { 
//  	       std::cout << "found a CA in " << ires << " "
//  			 << m[0][ires].name << " "
//  			 << " fragment_id  " << m[0].fragment_id
//  			 << " molecule: " << m.name << std::endl; 
		  coot::minimol::residue res(ires);
		  res.addatom(m[0][ires][iat]);
		  //	       f.residues.push_back(res);

		  try { 
		     f.addresidue(res, 0);
		  }
		  catch (std::runtime_error rte) {
		     std::cout << "ERROR:: get_fragment_ca_atoms() " << rte.what() << std::endl;
		  }
		  break;
	       }
	    }
	 }
      }
   }

   return f;
}

clipper::Matrix<float>
coot::db_main::make_cov_matrix(const std::vector<clipper::Coord_orth> &frag_coord) const {
   
   clipper::Matrix<float> mat(3,3); 
   float sum[3]; 
   for(int i=0; i<3; i++) sum[i] = 0; 

   for (unsigned int i=0; i<frag_coord.size(); i++) { 
      sum[0] += frag_coord[i].x(); 
      sum[1] += frag_coord[i].y(); 
      sum[2] += frag_coord[i].z(); 
   }

   clipper::Coord_orth mean_pos(sum[0]/float(frag_coord.size()), 
				sum[1]/float(frag_coord.size()), 
				sum[2]/float(frag_coord.size()));
   // this is symmetric, isn't it?
   //
//    std::cout << "Adding " << frag_coord.size() << " coords diffs to covmatrix "
// 	     << std::endl;
   for (unsigned int i=0; i<frag_coord.size(); i++) {
//       std::cout << "was " << mat(0,0) << " "; 
//       std::cout << mat(0,1) << " "; 
//       std::cout << mat(0,1) << std::endl; 
//       std::cout << "was " << mat(1,0) << " "; 
//       std::cout << mat(1,1) << " "; 
//       std::cout << mat(1,1) << std::endl; 
//       std::cout << "was " << mat(2,0) << " "; 
//       std::cout << mat(2,1) << " "; 
//       std::cout << mat(2,1) << std::endl; 
      mat(0,0) += (frag_coord[i].x() - mean_pos.x()) * 
                  (frag_coord[i].x() - mean_pos.x()); 
      mat(0,1) += (frag_coord[i].x() - mean_pos.x()) * 
                  (frag_coord[i].y() - mean_pos.y()); 
      mat(0,2) += (frag_coord[i].x() - mean_pos.x()) * 
                  (frag_coord[i].z() - mean_pos.z()); 
      mat(1,0) += (frag_coord[i].y() - mean_pos.y()) * 
                  (frag_coord[i].x() - mean_pos.x()); 
      mat(1,1) += (frag_coord[i].y() - mean_pos.y()) * 
                  (frag_coord[i].y() - mean_pos.y()); 
      mat(1,2) += (frag_coord[i].y() - mean_pos.y()) * 
                  (frag_coord[i].z() - mean_pos.z()); 
      mat(2,0) += (frag_coord[i].z() - mean_pos.z()) *  
                  (frag_coord[i].x() - mean_pos.x()); 
      mat(2,1) += (frag_coord[i].z() - mean_pos.z()) * 
                  (frag_coord[i].y() - mean_pos.y()); 
      mat(2,2) += (frag_coord[i].z() - mean_pos.z()) * 
                  (frag_coord[i].z() - mean_pos.z());
   }

//    for (int i=0; i<3; i++) {
//       std::cout << "| ";
//       for (int j=0; j<3; j++) {
// 	 std::cout << mat(i,j) << "  ";
//       }
//       std::cout << "|" << std::endl;
//    }
//    std::cout << std::endl;
   return mat;
}

std::vector<std::string>
coot::db_main::get_reference_pdb_list() const {

   // stat "MAPVIEW_REF_STRUCTS", opendir
   // 
   // match filenames,  read them all.
   //
   // 
   std::vector<std::string> pdb_list; 
   std::string mapview_ref_structs("COOT_REF_STRUCTS"); 
   
   char *dir = getenv(mapview_ref_structs.c_str()); 
   std::string d(DATADIR); // prefix/share
   std::string ref_struct_dir = coot::util::append_dir_dir(d, "coot");
   ref_struct_dir = coot::util::append_dir_dir(ref_struct_dir, "reference-structures");

   if (!dir) { 
      // Let's try the "standard" place DATADIR/coot/reference-structures

      struct stat buf;
      int istat = stat(ref_struct_dir.c_str(), &buf);
      if (istat == 0) { 
	 dir = (char *) ref_struct_dir.c_str();
      }
   }


   if (!dir) {
      std::cout << "ERROR! COOT_REF_STRUCTS is not defined.  \n"
		<< "       Can't find " << ref_struct_dir << ".\n"
		<< "       Cannot continue with mainchain building.\n";
   } else { 
      std::string ref_str_dir_str(dir); 
      
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
	       if (coot::matches_pdb_name(file_str)) {

		  // Construct a file in a unix directory - non-portable?
		  // 
		  std::string s(dir);
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
   } 
   return pdb_list; 
} 

bool
coot::matches_pdb_name(std::string file_str) { 

   bool match_flag = false; 

   // did it end in: .pdb .pdb.gz or was it pdb*.gz?

   std::string::size_type t1 = file_str.find(".pdb");
   std::string::size_type t2 = file_str.find(".pdb.gz");
   std::string::size_type t3 = file_str.find("pdb");
   std::string::size_type t4 = file_str.find(".gz");

   if (t1 != std::string::npos) match_flag = true;
   if (t2 != std::string::npos) match_flag = true;
   if (t3 != std::string::npos && t4 != std::string::npos) match_flag = true;
   
   return match_flag; 
}


std::vector<clipper::Coord_orth>
coot::db_main::get_target_ca_coords(int iresno, int ilength,
				    const coot::minimol::molecule &target) const {
      
   std::vector<clipper::Coord_orth> cs;

   //    std::cout << "get_target_ca_coords: iresno is " << iresno << std::endl; 

   if ( (iresno+ilength-1) <= target[0].max_residue_number()) {
      for (int ires=iresno; ires<(iresno+ilength); ires++) {
	 for (unsigned int iat=0; iat<target[0][ires].atoms.size(); iat++) {
// 	    std::cout << "DEBUG:: get_target_ca_coords " << ires << " :"
// 		      << target[0][ires][iat].name << ":" << std::endl;
	    if (target[0][ires][iat].name == " CA ") { 
	       cs.push_back(target[0][ires][iat].pos);
	    }
	 }
      }
   } 
   return cs; 
} 

// Fill big_results:
//
// First set the output fragment size.
//
// We will use residue numbers for the output fragment that match the
// input target fragment (target), i.e.: we start at
// iresno_target_start of target and go on to iresno_target_end.
// 
void
coot::db_main::match_target_fragment(const coot::minimol::molecule &target,
				     int iresno_target_start_in,  // typically 1
				     int iresno_target_end_in,
				     int ilength) {

   max_devi = 2.6; // the upper limit of the sum of the deviations
		   // between the overlayed test fragment and the
		   // target fragment Cas.

   int iresno = iresno_target_start_in;
   iresno_start = iresno_target_start_in;  // Fill the class members.
   iresno_end   = iresno_target_end_in;    // Needed when we merge.
   target_fragment_fragment_id = target[0].fragment_id; // needed for returned fragment

   // for a given residue number (in the the chain), what's the residue name?
   // 
   // std::map<int, std::string> sequence;
   //
   sequence = get_sequence(target,
			   target_fragment_fragment_id,
			   iresno_target_start_in,
			   iresno_target_end_in);

   if (target.fragments.size() > 0) { 
      
      while (1) {

	 std::vector<clipper::Coord_orth> target_ca =
	    get_target_ca_coords(iresno, ilength, target);

	 // Typically, the following test fails when we have reached the
	 // end of the target and there were no more Ca left to make up
	 // 6.  In future, this will need to be more sophisticated,
	 // because we actually do want to try to fit the last couple of
	 // residues (rather than (as we currently do) discard them).  We
	 // should step back istart so that it is 6 positions from the
	 // end, in that case.
	 //
	 if (int(target_ca.size()) != ilength) {
	    std::cout << "short target_ca " << target_ca.size()
		      << " with iresno " << iresno
		      << ", trying " << iresno -1 << std::endl;
	    iresno--;
	    target_ca = get_target_ca_coords(iresno, ilength, target);
	    if (int(target_ca.size()) != ilength) { 
	       std::cout << "short target_ca  " << target_ca.size()
			 << " with iresno " << iresno
			 << ", trying " << iresno -1 << std::endl;
	       if (iresno > iresno_target_start_in) { 
		  iresno--;
		  target_ca = get_target_ca_coords(iresno, ilength, target);
		  if (int(target_ca.size()) != ilength) {
		     std::cout << "short target_ca  " << target_ca.size()
			       << " with iresno " << iresno
			       << ", breaking " << std::endl;
		     break;
		  }
	       } else {
		  std::cout << "WARNING:: short target at start - bailing out now."
			    << std::endl;
		  break;
	       }
	    }
	 }

	 if (iresno + 5 > iresno_end) { // 41 + 5 = 46, OK.
	    std::cout << "breaking out of match_target_fragment 6 with "
		      << "iresno= " << iresno << std::endl; 
	    break; 
	 } 

	 std::vector<coot::db_fitting_result> fits =
	    fit_reference_structures(max_devi, target_ca, iresno, ilength);

	 std::cout << "target fragment start res: " << iresno
		   << " had " << fits.size() << " db fragment hits"
		   << ", ilength=" << ilength << std::endl;

	 if (fits.size() == 0) {
	    target_ca = get_target_ca_coords(iresno, ilength-1, target);
	    fits = fit_reference_structures(max_devi, target_ca, iresno, ilength-1);
	    big_results.push_back(fits);
	    std::cout << " retrying with ilength = " << ilength-1
		      << " produces " << fits.size() << " hits." << std::endl;
	    
	    target_ca = get_target_ca_coords(iresno+1, ilength-1, target);
	    fits = fit_reference_structures(max_devi, target_ca, iresno+1, ilength-1);
	    std::cout << " 2nd part with ilength = " << ilength-1
		      << ", iresno=" << iresno+1
		      << " produces " << fits.size() << " hits." << std::endl;
	    big_results.push_back(fits);
	 } else {
	    // we had hits for the normal size (6).
	    big_results.push_back(fits);
	 } 
	 iresno += ilength/2; // or some such
      }
   } else {
      std::cout << "There were no fragments in the target molecule"
		<< std::endl;
   }
}

std::map<int, std::string>
coot::db_main::get_sequence(const coot::minimol::molecule &target,
			    std::string target_fragment_fragment_id,
			    int iresno_target_start,
			    int iresno_target_end) const {

   std::map<int, std::string> m;

   if (iresno_target_start > iresno_target_end)
      std::swap(iresno_target_start, iresno_target_end);

   int chain_idx = 0; // is this the same as the calling function?
 
   minimol::fragment f = target[0];
   for (int ires=iresno_target_start; ires <= iresno_target_end; ires++) {
      try {
	 std::string rn = f[ires].name;
	 m[ires] = rn;
      }
      catch (const std::runtime_error &rte) {
	 std::cout << "ERROR:: " << rte.what() << std::endl;
      } 
   } 
   return m;
} 


std::vector<coot::db_fitting_result>
coot::db_main::fit_reference_structures(float max_devi,
					std::vector<clipper::Coord_orth> target_ca,
					int iresno_start,
					int ilength) const {

   std::vector<coot::db_fitting_result> fits;
   float devi; 
   
   clipper::Matrix<float> mat = make_cov_matrix(target_ca);
   std::vector<float> target_eigens = mat.eigen( true );
   for (unsigned int j=0; j<target_eigens.size(); j++)
      target_eigens[j] = sqrt(target_eigens[j]);

   for (unsigned int i=0; i<mainchain_frag_db.size(); i++) { // several 1000s.

      if (similar_eigens(0.4, target_eigens,
			 mainchain_frag_db[i].sqrt_eigen_values)) {
	 std::vector<clipper::Coord_orth> mcfca =
	    mainchain_ca_coords_of_db_frag(i, ilength);
	    
	 // and get devi
	 if (int(mcfca.size()) == ilength ) {
	    if (int(target_ca.size()) == ilength) {
	       clipper::RTop_orth rtop(mcfca, target_ca);
	       devi = deviance(mcfca, target_ca, rtop);
	       if (devi < max_devi) { 
		  fits.push_back(db_fitting_result(rtop, i, devi,
						   iresno_start, ilength));
	       } 
	    } else {
	       std::cout << "wrong target ca size" << std::endl;
	    }
	 } else {
	    // too many residues were included, due to ins codes, I guess.
// 	    std::cout << "wrong mcfca size: mcfca: " << mcfca.size()
// 	             << " ilength: " << ilength << " for fragment number "
// 		      << i << " of " << mainchain_frag_db.size() << std::endl;
	 }
      }
   }

   return fits; 
}
					

std::vector<clipper::Coord_orth>
coot::db_main::mainchain_ca_coords_of_db_frag(int i, int ilength) const {

   int ifraglen = mainchain_frag_db[i].ilength ; // typically 6
   ifraglen = ilength; 
   
   std::vector<clipper::Coord_orth> a;
   
   int imol_no = mainchain_frag_db[i].molecule_number;
   int ires_st = mainchain_frag_db[i].i_start_res();

   int ifrag = 0;

   if (molecule_list[imol_no][ifrag].max_residue_number() < ires_st+ifraglen) {
      std::cout << "ERROR: Trapped bad residue index : "
		<< "imol_no: " << imol_no << " fragment: " << ifrag << " wanted "
		<< ires_st << "+" << ifraglen << " but short residues: "
		<< molecule_list[imol_no][ifrag].max_residue_number() << " " 
		<< molecule_list[imol_no][ifrag].fragment_id << " "
		<< std::endl;
   } 
   for (int ires=ires_st; ires<ires_st+ifraglen; ires++) {
      if (molecule_list[imol_no][ifrag][ires].atoms.size() == 0) {
	 std::cout << "oops: zero atoms for residue " << ires
		   << " in molecule number " << imol_no
		   << std::endl; 
      }
      for (unsigned int iat=0; iat<molecule_list[imol_no][ifrag][ires].atoms.size(); iat++) {
	 if (molecule_list[imol_no][ifrag][ires][iat].name == " CA ") {
	    a.push_back(molecule_list[imol_no][ifrag][ires][iat].pos);
	 }
      }
   }
   return a;
}

short int
coot::db_main::similar_eigens(float tolerance_frac, 
			const std::vector<float> &target, 
			const std::vector<float> &frag) const { 

//     std::cout << "comparing " << target[0] << "  " << target[1] << "  "
//  	     << target[2] << std::endl; 
//     std::cout << "with      " << frag[0] << "  " << frag[1] << "  "
//   	     << frag[2] << std::endl; 
   
   short int sim = 1;
   for (unsigned int i=0; i<target.size(); i++) {
      if (target[i] > frag[i]*(1+tolerance_frac) ||
	  target[i] < frag[i]*(1-tolerance_frac)) { 
	 sim = 0; 
	 break; 
      }
   }
   return sim;
}

// we need to apply rtop to frag1
float
coot::db_main::deviance(const std::vector<clipper::Coord_orth> &frag1, 
			const std::vector<clipper::Coord_orth> &frag2, 
			const clipper::RTop_orth &rtop) const { 
   float devi = 0; 
   for (unsigned int i=0; i<frag1.size(); i++)
      devi += frag1[i].length(frag2[i], frag1[i].transform(rtop)); 
   return devi; 
}


void
coot::db_main::merge_fragments() {
   
   // There are 3 weights: 
   // 
   // weight_pos_in_frag (positions 2 and 3 get higher weights)
   // 
   // weight_target_devi (fragments that match the Ca target fragment
   //                    more closely get higher weights)
   // 
   // weight_prev_min (fragments that more closely match any element
   //                 of the previous fragment get higher weights)

   // ifrag is the fragment counter
   // i_frag_pos is the position within the fragment (0,1,2,3,4,5)
   // i_out_res is the "output" residue number
   // 
   // All the atoms of a single residue have the same weight applied
   // to them.
   //

   // weight_prev_min seems particularly difficult to implement.
   // Let's implement merge_fragment without it.  Lets call this
   // function merge_fragment_s and give _s suffix to all the specific
   // functions that that this simple (_s) implementation needs.  This
   // allows us easily to add a more complex merging later, when we
   // have thought how to do it (if we ever do (and we decide that it
   // is the right thing to do anyway)).

   // float weight_pos_in_frag;
   float weight_target_devi;
   float devi;
   minimol::residue fragment_res;
   output_fragment.resize(iresno_end-iresno_start+1); // recall: std::vector<weighted_residue>

   // we presume that output_fragment is sorted on istart_res
   // 
   std::cout << "merge fragments " << iresno_start << " to "
	     << iresno_end << " with " << big_results.size()
	     << " fit sets to merge" << std::endl;

//    std::cout << "---- fit set sizes: ----- " << std::endl;
//    for (int ibr=0; ibr<big_results.size(); ibr++) {
//       std::cout << ibr << "   " << big_results[ibr].size() << std::endl;
//    }


// show us where the fits were:   
//    for (int ibr=0; ibr<big_results.size(); ibr++) {
//       std::cout << "fit set: " << ibr << ": ";
//       for (int j=0; j<big_results[ibr].size(); j++) {
// 	 std::cout << big_results[ibr][j].istart_res_of_ca_target << " ";
//       }
//       std::cout << std::endl;
//    }
	    

   for (unsigned int ibr=0; ibr<big_results.size(); ibr++) {
      
      if (big_results[ibr].size() > 0) {
	 
	 for (unsigned int i=0; i<big_results[ibr].size(); i++) {
	    devi = big_results[ibr][i].deviance;

	    weight_target_devi = 1/(devi*devi+0.001); // stabilization

	    // std::cout << "weight_target_devi is " << weight_target_devi << std::endl;

	    int ilength = big_results[ibr][i].ilength;
	    // assign fragment_res and pos_in_frag
	    //
	    for (int ipos=0; ipos<ilength; ipos++) {
	       fragment_res = pull_db_residue(big_results[ibr][i], ipos);
	       float w = weight_pos_in_frag(ipos,ilength) * weight_target_devi;
	       int j = big_results[ibr][i].istart_res_of_ca_target + ipos;
	       if (j<=iresno_end) {
		  // std::cout << "merge_fragments adding to residue " << j << std::endl;
		  int of_idx = j - iresno_start;
		  output_fragment[of_idx].add_residue_pos(fragment_res,
							  big_results[ibr][i].rtop, w);
	       }
	    }
	 }
      }
   }
   std::cout << "The merging is complete" << std::endl;
}

coot::minimol::residue
coot::db_main::pull_db_residue(const coot::db_fitting_result &fit, int ipos) const {
   int imolno = mainchain_frag_db[fit.db_frag_index].molecule_number;
   int frag_start_res = mainchain_frag_db[fit.db_frag_index].i_start_res();
   
   coot::minimol::residue res = molecule_list[imolno][0][frag_start_res+ipos];
//    if (res.atoms.size() == 0) {
//       std::cout << "WARNING: no atoms in fragment_res, ipos ="
// 		<< ipos << ", frag_start_res=" << frag_start_res << std::endl;
//    }
   
   return res;
}

coot::minimol::fragment
coot::db_main::pull_db_fragment(const coot::main_fragment_t &dbfit, int ilength) {
   
   int imolno = dbfit.molecule_number;
   int frag_start_res = dbfit.i_start_res();
   std::string chain_id = dbfit.segment_id;

   coot::minimol::fragment f(chain_id);
      

   for (int ipos=0; ipos<ilength; ipos++) {
      coot::minimol::residue res = molecule_list[imolno][0][frag_start_res+ipos];
      try { 
	 f.addresidue(res, 0);
      }
      catch (std::runtime_error rte) {
	 std::cout << "ERROR:: pull_db_fragment() " << rte.what() << std::endl;
      } 
   }
   return f;
}


coot::minimol::residue
coot::weighted_residue::pull_residue() const {

   minimol::residue res;
   minimol::atom at(std::string("tmp"), std::string("tmp"), 0,0,0, "", 1.0, 30.0); // overridden

   for (unsigned int iat=0; iat<atoms.size(); iat++) {
      at.name = atoms[iat].name;
      if (at.name == " CB " ) { 
	 at.pos  = 1.0/weight_sum_cb * atoms[iat].pos;
      } else {
//  	 std::cout << "pulling: weight_sum: " << weight_sum << ", pos "
//  		   << atoms[iat].pos.format() << std::endl; 
	 at.pos  = 1.0/weight_sum * atoms[iat].pos;
      }
      at.element = atoms[iat].element; 
      res.atoms.push_back(at);
   }
   return res; 
}

coot::minimol::fragment
coot::db_main::mainchain_fragment() const {

   minimol::fragment frag;
   frag.fragment_id = target_fragment_fragment_id; 

   frag.residues.resize(iresno_end +1);

   std::cout << "mainchain_fragment from " << iresno_start << " to "
	     << iresno_end << std::endl;

   for (int ires=iresno_start; ires<=iresno_end; ires++) {
      int of_idx = ires - iresno_start;
      minimol::residue res = output_fragment[of_idx].pull_residue();
      for (unsigned int i=0; i<res.atoms.size(); i++) {
	 if (res.atoms[i].name == " CA " ||
	     res.atoms[i].name == " CB " ||
	     res.atoms[i].name == " C  " ||
	     res.atoms[i].name == " N  " ||
	     res.atoms[i].name == " O  " ) { 
	    frag[ires].addatom(res.atoms[i]);
	 }
      }

      std::string res_name = "ALA";
      if (false) {
	 // 20180706-PE - why would I want to do this?
	 std::map<int, std::string>::const_iterator it = sequence.find(ires);
	 if (it != sequence.end())
	    res_name = it->second;
      }
      frag[ires].name = res_name;
      frag[ires].seqnum = ires;
   }

   std::cout << "done mainchain_fragment" << std::endl; 
   return frag;
} 

// 3,6 -> 3
// 1,6 -> 2
// 0,6 -> 1
// 
// (l+1)/2 - |(l-1)/2 - i|
float
coot::db_main::weight_pos_in_frag(int ipos, int ilength) const {

   return float(ilength+1)*0.5 - fabs( float(ilength-1)*0.5 - ipos ); 
} 

// Given a residue with the weight_target_devi.  We apply
// weight_pos_in_frag... 
//
// Note that we implicitly presume that all in_res has CA, C, N and O,
// but only some have CB.  So we need to keep a separate track of the
// sum of the CB weights (which is different from the sum of all
// weights).
//
// Notice that when we start, the vector atoms is empty, so we do a
// check that only fails for the first time, where we do an addatom,
// rather than adding to the position of already esisting atoms.
//
// For the first time, we check that CA, C, N, O exist after we have
// added in_res.
// 
void
coot::weighted_residue::add_residue_pos(const minimol::residue &in_res,
					const clipper::RTop_orth &rtop,
					float weight) {

   clipper::Coord_orth contrib; 

   weight_sum += weight;

   // Now, do we need to add to the CB weight too (ie does the in_res
   // have a Cb?)
   // 
   // Normally, we just add to the atoms already here, but the first
   // time we need to addatom()s.


   if (atoms.size() > 0) {

      for (unsigned int jat=0; jat<in_res.atoms.size(); jat++) {
	 if (in_res[jat].name == " CB ") {
	    weight_sum_cb += weight;
	    if (have_cb_flag == 0) {
	       minimol::atom tmp_atom = in_res[jat];
	       tmp_atom.pos = weight * in_res[jat].pos.transform(rtop);
	       addatom(tmp_atom);
	       cb_index = atoms.size() -1; // index of previous atom
	       // std::cout << "info: cb_index is " << cb_index << std::endl;
	       have_cb_flag = 1;
	    } else {
	       contrib = weight * in_res.atoms[jat].pos.transform(rtop);
	       atoms[cb_index].pos += contrib;
	    }
	 }
      }
      
      for (unsigned int iat=0; iat<atoms.size(); iat++) {
	 if (atoms[iat].name != " CB ") { 
	    for (unsigned int jat=0; jat<in_res.atoms.size(); jat++) {
	       if (atoms[iat].name == in_res.atoms[jat].name) {
		  contrib = weight * in_res.atoms[jat].pos.transform(rtop);
// 		  std::cout << "residue " << in_res.seqnum << " adding "
// 			    << in_res.atoms[jat].pos.transform(rtop).format()
// 			    << " with weight " << weight << std::endl;
		     
		  atoms[iat].pos += contrib;
	       }
	    }
	 }
      }
   } else {
      minimol::atom tmp_atom("tmp", "tmp", 0,0,0, "", 1.0, 30.0); // overridden
      for (unsigned int jat=0; jat<in_res.atoms.size(); jat++) {
	 if (in_res[jat].name != " CB ") { 
	    tmp_atom     =  in_res[jat];
	    tmp_atom.pos = weight * in_res[jat].pos.transform(rtop);
	    addatom(tmp_atom);
	 }
      }
      
      // now check that at least N, C, CA, O exist.
      int n_atom_count = 0; 
      for (unsigned int iat=0; iat<atoms.size(); iat++) {
	 if (atoms[iat].name == " CA ") n_atom_count++; 
	 if (atoms[iat].name == " C  ") n_atom_count++; 
	 if (atoms[iat].name == " N  ") n_atom_count++; 
	 if (atoms[iat].name == " O  ") n_atom_count++; 
      }
      // Programmer Error:
      if (n_atom_count != 4) {
	 std::cout << "ERROR: DISASTER! wrong number of mainchain"
		   << " atoms initially added." << std::endl;
      }
   }
}


bool
coot::db_main::is_empty() const {

   if (mainchain_frag_db.size() == 0)
      return 1;
   else
      return 0;
}

bool
coot::db_main::is_empty_of_pepflips() const {

   return 0;
}

      

void
coot::db_main::clear_results() {
   big_results.clear();
   output_fragment.clear();
   pepflip_fragments.clear();
}



// --------------------------------------------------------------------------
//        ----------------------- Pepflip extras --------------------
// --------------------------------------------------------------------------
// 
void
coot::db_main::match_targets_for_pepflip(const minimol::fragment &target_ca_coords_5_res_frag) {

   std::vector<peptide_match_fragment_info_t> fits;
   float devi;
   int ilength = 5;
   std::vector<clipper::Coord_orth> target_ca;
//    std::cout << "Running over target fragment, getting CA vector from "
// 	     << target_ca_coords_5_res_frag.min_res_no() << " to "
// 	     << target_ca_coords_5_res_frag.max_residue_number()
// 	     << std::endl;
   
   for (int ires=target_ca_coords_5_res_frag.min_res_no();
	ires<=target_ca_coords_5_res_frag.max_residue_number();
	ires++) {
//       std::cout << "     in CA vector " << ires << " "
// 		<< target_ca_coords_5_res_frag[ires] << std::endl;
      for (unsigned int iat=0; iat<target_ca_coords_5_res_frag[ires].atoms.size();
	   iat++) {
	 if (target_ca_coords_5_res_frag[ires][iat].name == " CA ") {
	    target_ca.push_back(target_ca_coords_5_res_frag[ires][iat].pos);
	 }
      }
   }

   // std::cout << "atom_count for vect: " << target_ca.size() << std::endl;
   if (target_ca.size() == 5) { 
      clipper::Matrix<float> mat = make_cov_matrix(target_ca);
      std::vector<float> target_eigens = mat.eigen( true );
//       std::cout << "Eigens: "
// 		<< target_eigens[0] << " "
// 		<< target_eigens[1] << " "
// 		<< target_eigens[2] << " " << std::endl;
      for (unsigned int j=0; j<target_eigens.size(); j++)
	 target_eigens[j] = sqrt(target_eigens[j]);

      float best_devi = 9999999.9;

      
      assign_eigen_similarity_scores(target_eigens);
      sort_mainchain_fragments_by_eigens(target_eigens);
      const unsigned int max_frag_count = 100;
      int frag_count = 0;
      
      for (unsigned int i=0; (i<mainchain_frag_db.size() && i<max_frag_count); i++) { // several 1000s.
	 if (mainchain_frag_db[i].ilength == 5) { 
// 	 std::cout << "eigen similarilty score: " << mainchain_frag_db[i].eigen_similarity_score
// 		   << std::endl;
	    frag_count++;

	    if (similar_eigens(0.2, target_eigens,
			       mainchain_frag_db[i].sqrt_eigen_values)) {
	       // std::cout << "eigens OK" << std::endl;
	       std::vector<clipper::Coord_orth> mcfca =
		  mainchain_ca_coords_of_db_frag(i, ilength);

	       // and get devi
	       if (int(mcfca.size()) == ilength ) {
		  if (int(target_ca.size()) == ilength) {
		     clipper::RTop_orth rtop(mcfca, target_ca);
		     coot::minimol::fragment db_fragment =
			pull_db_fragment(mainchain_frag_db[i], ilength);
		     db_fragment.transform(rtop);
		     float this_devi = deviance(mcfca, target_ca, rtop);
		     // std::cout << "d: " << this_devi << "\n";
		     if (this_devi < best_devi)
			best_devi = this_devi;
		     fits.push_back(coot::peptide_match_fragment_info_t(db_fragment, this_devi));
		  } else {
		     std::cout << "wrong target ca size" << std::endl;
		  }
	       }
	    }
	 }
      }
      //       std::cout << "best devi: " << best_devi << std::endl;
   }
//    std::cout << "Generated " << fits.size() << " fits from "
// 	     << mainchain_frag_db.size() << " candidates" << std::endl;
   pepflip_fragments = fits;
}


// Check that internal vector against the oxygen position.
//
// Return the fraction fo the best fitting peptides that have their
// oxygen closer than d_crit
float 
coot::db_main::mid_oxt_outliers(const clipper::Coord_orth &my_peptide_oxt_pos, int resno_oxt, float d_crit) {

   float ret_frac = -1.;
   int n_better = 0;
   const int frag_count_max = 70;
   int frag_count = 0;

   std::sort(pepflip_fragments.begin(), pepflip_fragments.end(), pepflip_sorter);
   std::vector<clipper::Coord_orth> matched_oxt_pos;
//    for (unsigned int ifit=0; (ifit<pepflip_fragments.size()) && (ifit<5) ; ifit++) {
//       std::cout << ifit << " " << pepflip_fragments[ifit].devi << "\n";
//    }

   for (unsigned int ifit=0; (ifit<pepflip_fragments.size()) && (frag_count<frag_count_max) ; ifit++) {
      //       std::cout << "pepflip_fragments devi: " << pepflip_fragments[ifit].devi << std::endl;

      int rescount=0;
      for (int ires=pepflip_fragments[ifit].fragment.min_res_no();
	   ires<=pepflip_fragments[ifit].fragment.max_residue_number();
	   ires++) {
	 if (pepflip_fragments[ifit].fragment[ires].atoms.size() > 0) {
	       rescount++;
	 }
	 if (rescount==3) {
	    for (unsigned int iat=0; iat<pepflip_fragments[ifit].fragment[ires].atoms.size();
		 iat++) {
	       if (pepflip_fragments[ifit].fragment[ires][iat].name == " O  ") {

		  float l =
		     clipper::Coord_orth::length(my_peptide_oxt_pos,
						 pepflip_fragments[ifit].fragment[ires][iat].pos);
// 		  std::cout << "devis_o_and_fit: " << l << "   "
// 			    << pepflip_fragments[ifit].devi/5.0
// 			    << std::endl;
		  matched_oxt_pos.push_back(pepflip_fragments[ifit].fragment[ires][iat].pos);
		  frag_count++;
		  if (l < d_crit)
		     n_better++;
		  break;
	       }
	    }
	 }
      } 
   }

   if (matched_oxt_pos.size() > 0) { 
      // get mid point of matched oxygens and dist from our oxt post to
      // that mid point.
      clipper::Coord_orth sum_mid_point(0.0,0.0,0.0);
      for (unsigned int ipt=0; ipt<matched_oxt_pos.size(); ipt++) {
	 sum_mid_point += matched_oxt_pos[ipt];
      }
      float div = 1.0/float(matched_oxt_pos.size());
      clipper::Coord_orth mid_point = clipper::Coord_orth(sum_mid_point.x()*div,
							  sum_mid_point.y()*div,
							  sum_mid_point.z()*div);
      double sum_d_sq = 0.0;
      double d_sq = 0.0;
      for (unsigned int ipt=0; ipt<matched_oxt_pos.size(); ipt++) {
	 double d = clipper::Coord_orth::length(matched_oxt_pos[ipt], mid_point);
	 sum_d_sq += d*d;
      }
      double rmsd = sqrt(sum_d_sq/float(matched_oxt_pos.size()));
      double midpt_to_my_oxt = clipper::Coord_orth::length(mid_point, my_peptide_oxt_pos);
      double z = 1.0;
      if (rmsd > 0.0) 
	 z = midpt_to_my_oxt/rmsd;

      // Now, Jones et al. (1991) use the rmsd of the distances from
      // the matched point to the model point and suggest values more
      // than 2.5 are worth investigating.
      sum_d_sq = 0.0;
      d_sq = 0.0;
      for (unsigned int ipt=0; ipt<matched_oxt_pos.size(); ipt++) {
	 double d = clipper::Coord_orth::length(matched_oxt_pos[ipt], my_peptide_oxt_pos);
	 sum_d_sq += d*d;
      }
      double jones_rms = sqrt(sum_d_sq/float(matched_oxt_pos.size()));
      
      double z_jones = 1;
      if (rmsd > 0.0) 
	 z_jones = jones_rms/rmsd;
   
      ret_frac = float(n_better)/float(frag_count);
      std::cout << "z_jones_frac: " << resno_oxt << " " << z << " " << z_jones << " "
		<< jones_rms << " " << rmsd << " " << 100.0*ret_frac <<  "% ["
		<< pepflip_fragments.size() << " samples]" << std::endl;
   }
   return ret_frac;
}

// static
bool
coot::db_main::pepflip_sorter(const coot::peptide_match_fragment_info_t &fit_a,
			      const coot::peptide_match_fragment_info_t &fit_b) {
   return (fit_a.devi < fit_b.devi);
}

// static
bool
coot::db_main::mainchain_fragment_sorter(const coot::main_fragment_t &fit_a,
					 const coot::main_fragment_t &fit_b) {

   return (fit_a.eigen_similarity_score < fit_b.eigen_similarity_score);
}

void
coot::db_main::assign_eigen_similarity_scores(const std::vector<float> &target_eigens) {

   for (unsigned int i=0; i<mainchain_frag_db.size(); i++) { // several 1000s.
      float score = 0.0;
      for (unsigned int j=0; j<target_eigens.size(); j++) {
	 score += fabs(target_eigens[j]-mainchain_frag_db[i].sqrt_eigen_values[j]);
      }
//       std::cout << i << " ["
// 		<< mainchain_frag_db[i].sqrt_eigen_values[0] << " "
// 		<< mainchain_frag_db[i].sqrt_eigen_values[1] << " "
// 		<< mainchain_frag_db[i].sqrt_eigen_values[2] << "] vs ["
// 		<< target_eigens[0] << " "
// 		<< target_eigens[1] << " "
// 		<< target_eigens[2] << "]  ---> score: " << score << "\n";
	 
      mainchain_frag_db[i].eigen_similarity_score = score;
   }
}

void
coot::db_main::sort_mainchain_fragments_by_eigens(std::vector<float> target_eigens) {

   tmp_target_eigens = target_eigens;
   std::sort(mainchain_frag_db.begin(), mainchain_frag_db.end(), mainchain_fragment_sorter);

}
