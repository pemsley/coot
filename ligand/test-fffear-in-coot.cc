/* ligand/test-fffear-in-coot.cc
 * 
 * Copyright 2005 by The University of York
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

// Portability gubbins
#include <unistd.h> // for getopt(3)

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#ifdef __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#else
#define __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#undef __GNU_LIBRARY__
#endif

#include <iostream>

#include <sys/types.h> // for stating
#include <sys/stat.h>
#include <unistd.h>

#include "ligand.hh"
#include "utils/coot-utils.hh"
#include "coot-utils/coot-map-utils.hh"
#include "clipper/core/map_utils.h"
#include "clipper/ccp4/ccp4_map_io.h"
#include "coot-utils/coot-map-heavy.hh"
#include "coot-utils/atom-selection-container.hh"

int
main(int argc, char **argv) {

   if (argc < 6) { 
      std::cout << "Usage: " << argv[0] 
		<< " --pdbin pdb-in-filename" << " --hklin mtz-filename"
		<< " --f f_col_label"
		<< " --phi phi_col_label"
		<< " --pdbout fitted-molecule-filename"
		<< "\n"
		<< "        --mapin ccp4-map-name can be used"
		<< " instead of --hklin --f --phi\n";
      std::cout << "        where pdbin is the protein (typically)\n"
		<< "        and pdbout is file for the fitted molecule.\n";

   } else { 

      std::string pdb_file_name;
      std::string  mtz_filename;
      std::string         f_col;
      std::string       phi_col;
      std::string    output_pdb;
      std::string map_file_name;
      std::string map_out_file_name;
      std::string resolution_string;
      float resolution = 1; // set later if needed.
      float angle_step = 15;
      bool trans_only = false;

      const char *optstr = "i:h:f:p:r:o:m";
      struct option long_options[] = {
	 {"pdbin",  1, 0, 0},
	 {"hklin",  1, 0, 0},
	 {"f",      1, 0, 0},
	 {"phi",    1, 0, 0},
	 {"pdbout", 1, 0, 0},
	 {"resolution", 1, 0, 0},
	 {"mapin",  1, 0, 0},
	 {"mapout", 1, 0, 0},
	 {"angle-step", 1, 0, 0},
	 {"trans-only", 0, 0, 0},
	 {0, 0, 0, 0}
      };

      int ch;
      int option_index = 0;
      while ( -1 != 
	      (ch = coot_getopt_long(argc, argv, optstr, long_options, &option_index))) {

	 switch(ch) {

	 case 0:
	    if (coot_optarg) { 
	       std::string arg_str = long_options[option_index].name;

	       if (arg_str == "pdbin") { 
		  pdb_file_name = coot_optarg;
	       } 
	       if (arg_str == "pdbout") { 
		  output_pdb = coot_optarg;
	       } 
	       if (arg_str == "hklin") { 
		  mtz_filename = coot_optarg;
	       } 
	       if (arg_str == "f") { 
		  f_col = coot_optarg;
	       } 
	       if (arg_str == "phi") {
		  phi_col = coot_optarg;
	       } 
	       if (arg_str == "mapin") {
		  map_file_name = coot_optarg;
	       }
	       if (arg_str == "mapout") {
		  map_out_file_name = coot_optarg;
	       }
	       if (arg_str == "resolution") {
		  resolution_string = coot_optarg;
	       }
	       if (arg_str == "angle-step") {
		  angle_step = atof(coot_optarg);
	       }
	       
	    } else { 
	       std::string arg_str = long_options[option_index].name;
	       if (arg_str == "trans-only")
		  trans_only=true;
	    }
	    break;

	 case 'i':
	    pdb_file_name = coot_optarg;
	    break;
	    
	 case 'o':
	    output_pdb = coot_optarg;
	    break;
	    
	 case 'h':
	    mtz_filename = coot_optarg;
	    break;
	    
	 case 'f':
	    f_col = coot_optarg;
	    break;
	    
	 case 'p':
	    phi_col = coot_optarg;
	    break;

	 case 'r':
	    resolution = atof(coot_optarg);
	    
	 default:
	    break;
	 }
      }

      short int do_it = 0;
      short int do_it_with_map = 0;
      short int have_map = 0;
      short int have_reso_limit = 0;
      if (map_file_name.length() > 0)
	 have_map = 1;
      
      if (pdb_file_name.length() == 0) {
	 std::cout << "Missing input PDB file\n";
	 exit(1);
      } else { 
	 if (output_pdb.length() == 0) { 
	    std::cout << "Missing output PDB file\n";
	    exit(1);
	 } else { 
	    if (map_out_file_name.length() == 0) {
	       std::cout << "Missing output map name\n";
	       exit(1);
	    } else {
	       if (have_map) {
		  do_it_with_map = 1;
		  do_it = 1;
	       } else { 
		  if (mtz_filename.length() == 0) { 
		     std::cout << "Missing MTZ file\n";
		     exit(1);
		  } else { 
		     if (f_col.length() == 0) { 
			std::cout << "Missing F column name\n";
			exit(1);
		     } else {
			if (phi_col.length() == 0) { 
			   std::cout << "Missing PHI column name\n";
			   exit(1);
			} else {
			   if (resolution_string.length() > 0) {
			      resolution = atof(resolution_string.c_str());
			      have_reso_limit = 1;
			   }
			   do_it = 1;
			}
		     }
		  }
	       }
	    }
	 }
      }
      
      if (do_it) { 

	 short int use_weights = 0;
	 short int is_diff_map = 0; 
	 clipper::Xmap<float> xmap;
	 if (have_map) {

	    struct stat buf;
	    int status =  stat(map_file_name.c_str(), &buf);
	    if (status == 0) { 
	       clipper::CCP4MAPfile file;
	       file.open_read(map_file_name);
	       file.import_xmap(xmap);
	       file.close_read();
	    } else {
	       std::cout << " Error: Can't find map file: " << map_file_name << std::endl;
	       exit(1);
	    } 
	 } else {
	    coot::util::map_fill_from_mtz(&xmap, mtz_filename,
					  f_col, phi_col, "", 
					  use_weights, is_diff_map,
					  resolution, have_reso_limit);
	 }

	 atom_selection_container_t atom_sel =
	    get_atom_selection(pdb_file_name, true, false, false);

	 coot::util::fffear_search f(atom_sel.mol,
				     atom_sel.SelectionHandle,
				     xmap, angle_step, trans_only);
	 clipper::Xmap<float> results_map = f.get_results_map();
	 clipper::CCP4MAPfile mapout;
	 mapout.open_write(map_out_file_name);
	 mapout.export_xmap(results_map);
	 mapout.close_write();

	 clipper::RTop_orth mid_point_transformation =
	    f.mid_point_transformation();
	 std::vector<std::pair<float, clipper::RTop_orth> > p =
	    f.scored_orientations();
	 std::cout << "mid point transformation: \n" << mid_point_transformation.format()
		   << std::endl;

	 if (p.size() > 0) {

	    clipper::RTop_orth rtop = p[0].second;
	    std::vector<clipper::RTop_orth> possible_rtops;

	    std::cout << "Looking at permutations of \n" << rtop.format()
		      << " with peak height " << p[0].first << std::endl;

	    possible_rtops.push_back(rtop);
	    possible_rtops.push_back(clipper::RTop_orth(rtop.rot(),           -rtop.trn()));
	    possible_rtops.push_back(clipper::RTop_orth(rtop.rot().inverse(),  rtop.trn()));
	    possible_rtops.push_back(clipper::RTop_orth(rtop.rot().inverse(), -rtop.trn()));

	    std::vector<clipper::Coord_orth> peaks;
	    for (unsigned int iop=0; iop<p.size(); iop++) { 
	       std::cout << "Translation match at position : " << p[iop].second.trn().format()
			 << " peak value " << p[iop].first << std::endl;
	       peaks.push_back(clipper::Coord_orth(p[iop].second.trn()));
	    }

	    coot::minimol::molecule peaks_mol(peaks, "HOH", " OW1", "A");
	    peaks_mol.write_file("fffear-search-peaks.pdb", 20.0); 
	    
	    for (unsigned int iop=0; iop<possible_rtops.size(); iop++) { 
	       clipper::RTop_orth rtop_local = possible_rtops[iop];
	       mmdb::Manager *new_mol = new mmdb::Manager;
	       new_mol->Copy(atom_sel.mol, mmdb::MMDBFCM_All);
	       // now apply rtop to all atoms of new_mol:
	       int imod = 1;
	       
	       mmdb::Model *model_p = new_mol->GetModel(imod);
	       mmdb::Chain *chain_p;
	       // run over chains of the existing mol
	       int nchains = model_p->GetNumberOfChains();
	       for (int ichain=0; ichain<nchains; ichain++) {
		  chain_p = model_p->GetChain(ichain);
		  int nres = chain_p->GetNumberOfResidues();
		  mmdb::PResidue residue_p;
		  mmdb::Atom *at;
		  for (int ires=0; ires<nres; ires++) { 
		     residue_p = chain_p->GetResidue(ires);
		     int n_atoms = residue_p->GetNumberOfAtoms();
		     
		     for (int iat=0; iat<n_atoms; iat++) {
			at = residue_p->GetAtom(iat);
			clipper::Coord_orth co(at->x, at->y, at->z);
			co += mid_point_transformation.trn();
			clipper::Coord_orth pt = co.transform(rtop_local);
			at->x = pt.x();
			at->y = pt.y();
			at->z = pt.z();
		     }
		  }
	       }
	       // write the file:
	       std::string filename("fffear-ori-test-");
	       filename += coot::util::int_to_string(iop);
	       filename += ".pdb";
	       new_mol->WritePDBASCII(filename.c_str());
	    }

	    for (unsigned int iop=0; iop<p.size(); iop++) {
	       clipper::RTop_orth rtop_local = p[iop].second;
	       mmdb::Manager *new_mol = new mmdb::Manager;
	       new_mol->Copy(atom_sel.mol, mmdb::MMDBFCM_All);
	       // now apply rtop to all atoms of new_mol:
	       int imod = 1;
	       
	       mmdb::Model *model_p = new_mol->GetModel(imod);
	       mmdb::Chain *chain_p;
	       // run over chains of the existing mol
	       int nchains = model_p->GetNumberOfChains();
	       for (int ichain=0; ichain<nchains; ichain++) {
		  chain_p = model_p->GetChain(ichain);
		  int nres = chain_p->GetNumberOfResidues();
		  mmdb::PResidue residue_p;
		  mmdb::Atom *at;
		  for (int ires=0; ires<nres; ires++) { 
		     residue_p = chain_p->GetResidue(ires);
		     int n_atoms = residue_p->GetNumberOfAtoms();
		     
		     for (int iat=0; iat<n_atoms; iat++) {
			at = residue_p->GetAtom(iat);
			clipper::Coord_orth co(at->x, at->y, at->z);
			co += mid_point_transformation.trn();
			clipper::Coord_orth pt = co.transform(rtop_local);
			at->x = pt.x();
			at->y = pt.y();
			at->z = pt.z();
		     }
		  }
	       }

	       // write the file:
	       std::string filename("fffear-results-");
	       filename += coot::util::int_to_string(iop);
	       filename += ".pdb";
	       new_mol->WritePDBASCII(filename.c_str());
	    }
	    
	 } else { 
	    std::cout << "No peaks found, sadly" << std::endl;
	 }
      }
   }
   return 0;
}

