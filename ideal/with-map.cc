/* ideal/with-map.cc
 * 
 * Copyright 2002, 2003 The University of York
 * Author: Paul Emsley
 * Copyright 2007, 2008 The University of Oxford
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

 
// York
// #define MTZFILENAME "/h/paule/data/rnase/rnasa-1.8-all_refmac1.mtz"
// Cambridge
// #define MTZFILENAME "/h/paul/data/rnase/rnasa-1.8-all_refmac1.mtz"
// Glasgow
#define MTZFILENAME "/home/paule/data/rnase/rnasa-1.8-all_refmac1.mtz"

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#ifdef __GNU_LIBRARY__
#include "coot-getopt.h"
#else
#define __GNU_LIBRARY__
#include "coot-getopt.h"
#undef __GNU_LIBRARY__
#endif

#include <sys/types.h> // for stating
#include <sys/stat.h>
#ifndef _MSC_VER
#include <unistd.h>
#endif

#include <iostream>
#include <string.h>
#include <math.h>

#include "clipper/core/xmap.h"
#include "clipper/core/map_utils.h"
#include "clipper/ccp4/ccp4_map_io.h"
#include "clipper/ccp4/ccp4_mtz_io.h"
#include "clipper/core/hkl_compute.h"
#include "clipper/core/map_utils.h" // Map_stats

#define UNSET -9999

#ifndef  __MMDB_MMCIF__
#include "mmdb_mmcif.h"
#endif
  
#include <iostream>
#include <string>
#include <vector>

#include "mmdb_manager.h"
#include "mmdb-extras.h"
#include "mmdb.h"

#include "simple-restraint.hh"

// for debugging (the density_around_point function)
#include "coot-map-utils.hh"

std::pair<bool,clipper::Xmap<float> >
map_from_mtz(std::string mtz_file_name,
	     std::string f_col,
	     std::string phi_col,
	     std::string weight_col,
	     int use_weights,
	     int is_diff_map,
	     bool is_debug_mode);

class input_data_t {
public:
   bool is_good;
   bool is_debug_mode; 
   bool given_map_flag;
   bool use_rama_targets;
   bool use_torsion_targets;
   int resno_start;
   int resno_end;
   float map_weight; 
   std::string chain_id;
   std::string mtz_file_name;
   std::string f_col;
   std::string phi_col;
   std::string map_file_name;
   std::string  input_pdb_file_name;
   std::string output_pdb_file_name;

};

input_data_t get_input_details(int argc, char **argv);


std::vector<std::pair<bool,CResidue *> >
fill_residues(const std::string &chain_id, int resno_start, int resno_end, CMMDBManager *mol) {

   std::vector<std::pair<bool,CResidue *> > v;

   
   int imod = 1;
   CModel *model_p = mol->GetModel(imod);
   CChain *chain_p;
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      std::string this_chain_id = chain_p->GetChainID();
      if (this_chain_id == chain_id) { 
	 int nres = chain_p->GetNumberOfResidues();
	 CResidue *residue_p;
	 CAtom *at;
	 for (int ires=0; ires<nres; ires++) {
	    residue_p = chain_p->GetResidue(ires);
	    int this_res_no = residue_p->GetSeqNum();
	    if (this_res_no >= resno_start) {
	       if (this_res_no <= resno_end) { 
		  std::pair<bool, CResidue *> p(0, residue_p);
		  v.push_back(p);
	       }
	    }
	 }
      }
   }
   return v;
} 



int
main(int argc, char **argv) {

#ifndef HAVE_GSL
   std::cout << "We don't have GSL, this program does nothing" << std::endl;
#else
   
   std::string dict_filename;
   coot::protein_geometry geom;

   std::string mtz_filename(MTZFILENAME);
   std::string f_col("FWT");
   std::string phi_col("PHWT");

   if (argc < 2) {
      std::cout << "Usage: " << argv[0] << "\n"
		<< "       --pdbin pdb-in-filename\n"
		<< "       --hklin mtz-filename\n"
		<< "       --f f_col_label\n"
		<< "       --phi phi_col_label\n"
		<< "       --pdbout output-filename\n"
		<< "       --resno-start resno_low\n"
		<< "       --resno-end resno_high\n"
		<< "       --chain-id chain-id\n"
		<< "       --weight w (weight of map gradients, default 60)\n"
		<< "       --rama\n"
		<< "\n"
		<< "     --mapin ccp4-map-name can be used\n"
		<< "       instead of --hklin --f --phi\n"
		<< std::endl;

   } else {

      geom.set_verbose(0);
      geom.init_standard();

      //coot::restraints_container_t restraints(asc);

      // So, we provide easy(?) access to the atoms of next and
      // previous residues (to those in the atom selection
      // moving_residue_atoms).  It is also possible to select "fixed"
      // atoms in the graphics (so that they don't move).  Let's
      // provide a vector of indices in the the moving_residue_atoms
      // array to define those (lovely mixture of styles - heh).
      //
      std::vector<coot::atom_spec_t> fixed_atom_specs;

      // This interface has been withdrawn because we need the whole
      // molecule (acutally, a pointer to it) to do some atom selection.
      // 
//       coot::restraints_container_t restraints(asc.atom_selection, // moving_residue_atoms,
// 					      asc.n_selected_atoms,
// 					      previous_residue,
// 					      next_residue,
// 					      fixed_atoms);
      
      input_data_t inputs = get_input_details(argc, argv);

      if (inputs.is_good) { 

	 string pdb_file_name(inputs.input_pdb_file_name);
	 bool map_is_good = 0;

	 // if pdb_file_name does not exist -> crash?
	 atom_selection_container_t asc = get_atom_selection(pdb_file_name);
      
	 const char *chain_id = inputs.chain_id.c_str();
	 clipper::Xmap<float> xmap;

	 if (! inputs.given_map_flag) {
	    std::pair<bool, clipper::Xmap<float> > xmap_pair = 
	    map_from_mtz(inputs.mtz_file_name, inputs.f_col,
			 inputs.phi_col, "", 0, 0, inputs.is_debug_mode);
	    xmap = xmap_pair.second;
	    map_is_good = xmap_pair.first;
	    if (inputs.is_debug_mode) { 
	       clipper::Coord_orth pt(51.148,   8.121,  -1.418);
	       coot::util::density_stats_info_t ds =
	       coot::util::density_around_point(pt, xmap, 10.0);
	       
	       std::pair<float, float> mv = ds.mean_and_variance();
	       std::cout << "INFO:: density stats: N:" << ds.n << " sum: " << ds.sum
		      << " mean: " << mv.first << "  variance: "
			 << mv.second << std::endl;
	    }
	    
	 } else { 
	    clipper::CCP4MAPfile file;
	    file.open_read(inputs.map_file_name);
	    file.import_xmap(xmap);
	    file.close_read();
	    map_is_good = 1;
	 }

	 float map_weight = 60.0;
	 if (inputs.map_weight > 0) { 
	    map_weight = inputs.map_weight;
	    std::cout << "INFO:: using weight " << map_weight << std::endl;
	 }

	 if (map_is_good) { 
	    std::string altloc("");

	    std::vector<std::pair<bool,CResidue *> > local_residues =
	       fill_residues(chain_id, inputs.resno_start, inputs.resno_end, asc.mol);

	    // 	 int have_flanking_residue_at_start = 0;
	    // 	 int have_flanking_residue_at_end   = 0;
	    // 	 int have_disulfide_residues        = 0;
	 
	    // 	 coot::restraints_container_t restraints(inputs.resno_start,
	    // 						 inputs.resno_end,
	    // 						 have_flanking_residue_at_start,
	    // 						 have_flanking_residue_at_end,
	    // 						 have_disulfide_residues,
	    // 						 altloc,
	    // 						 chain_id,
	    // 						 asc.mol,
	    // 						 fixed_atom_specs);

	    coot::restraints_container_t restraints(local_residues, geom,
						    asc.mol,
						    fixed_atom_specs);
      
	    restraints.add_map(xmap, map_weight);

	    // coot::restraint_usage_Flags flags = coot::NO_GEOMETRY_RESTRAINTS;
	    // coot::restraint_usage_Flags flags = coot::BONDS;
	    // coot::restraint_usage_Flags flags = coot::BONDS_AND_ANGLES;
	    // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_AND_PLANES;
	    // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_TORSIONS_AND_PLANES; 
	    // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_AND_NON_BONDED;
	    // flags = coot::NON_BONDED;
	    coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
	    if (inputs.use_torsion_targets) {
	       flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS;
	       if (inputs.use_rama_targets)
		  flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_RAMA;
	    } else { 
	       if (inputs.use_rama_targets)
		  flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_CHIRALS_AND_RAMA;
	    }

	    coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
	    bool do_rama_plot_restraints = 0;
	    if (inputs.use_rama_targets)
	       do_rama_plot_restraints = 1;
	    restraints.make_restraints(geom, flags, 1, 1.0, do_rama_plot_restraints, pseudos);

	    restraints.minimize(flags);
	    restraints.write_new_atoms(inputs.output_pdb_file_name);
	 }
      }
   }

#endif // HAVE_GSL
   return 0; 
} 

std::pair<bool, clipper::Xmap<float> > 
map_from_mtz(std::string mtz_file_name,
	     std::string f_col,
	     std::string phi_col,
	     std::string weight_col,
	     int use_weights,
	     int is_diff_map,
	     bool is_debug_mode) {


   bool status = 0; // not filled
   
   clipper::HKL_info myhkl; 
   clipper::MTZdataset myset; 
   clipper::MTZcrystal myxtl; 
   clipper::Xmap<float> xmap;

   try { 
      cout << "reading mtz file..." << endl; 
      clipper::CCP4MTZfile mtzin; 
      mtzin.open_read( mtz_file_name );       // open new file 
      mtzin.import_hkl_info( myhkl );         // read sg, cell, reso, hkls
      clipper::HKL_data< clipper::datatypes::F_sigF<float> >   f_sigf_data(myhkl, myxtl);
      clipper::HKL_data< clipper::datatypes::Phi_fom<float> > phi_fom_data(myhkl, myxtl);
      clipper::HKL_data< clipper::datatypes::F_phi<float> >       fphidata(myhkl, myxtl); 
   
      std::string mol_name = mtz_file_name + " "; 
      mol_name += f_col; 
      mol_name += " "; 
      mol_name += phi_col; 
   
      if (use_weights) { 
	 mol_name += " ";
	 mol_name += weight_col; 
      } 

      if ( use_weights ) {
	 clipper::String dataname = "/*/*/[" + f_col + " " + f_col + "]";
	 std::cout << dataname << "\n";
	 mtzin.import_hkl_data(  f_sigf_data, myset, myxtl, dataname ); 
	 dataname = "/*/*/[" + phi_col + " " + weight_col + "]";
	 std::cout << dataname << "\n";
	 mtzin.import_hkl_data( phi_fom_data, myset, myxtl, dataname );
	 mtzin.close_read(); 
	 cout << "We should use the weights: " << weight_col << endl;
	 // it seems to me that we should make 2 data types, an F_sigF and a phi fom
	 // and then combine them using a Convert_fsigf_phifom_to_fphi();

	 fphidata.compute(f_sigf_data, phi_fom_data,
			  clipper::datatypes::Compute_fphi_from_fsigf_phifom<float>());

      } else {
	 clipper::String dataname = "/*/*/[" + f_col + " " + phi_col + "]";
	 mtzin.import_hkl_data(     fphidata, myset, myxtl, dataname );
	 mtzin.close_read(); 
      }
      std::cout << "Number of reflections: " << myhkl.num_reflections() << "\n"; 
      xmap.init( myhkl.spacegroup(), myhkl.cell(),
		 clipper::Grid_sampling( myhkl.spacegroup(),
					 myhkl.cell(),
					 myhkl.resolution()) );
      cout << "Grid..." << xmap.grid_sampling().format() << "\n";
      cout << "doing fft..." << endl;

      if (is_debug_mode) { 
	 int count = 0; 
	 clipper::HKL_info::HKL_reference_index hri;
	 for (hri=fphidata.first(); !hri.last(); hri.next()) {
	    if (count == 10)
	       break;
	    std::cout << "sample data " << " "
		      << hri.hkl().h() << " " 
		      << hri.hkl().k() << " " 
		      << hri.hkl().l() << " : " 
		      << fphidata[hri].f() << " " << fphidata[hri].phi()*180/M_PI << std::endl;
	    count++;
	 }
      }
  
  
      xmap.fft_from( fphidata );                  // generate map
      cout << "done fft..." << endl;
      status = 1;
   }
   catch (clipper::Message_base exception) {
      std::cout << "Failed to read mtz file " << mtz_file_name << std::endl;
   }

   return std::pair<bool, clipper::Xmap<float> > (status, xmap);
}


input_data_t
get_input_details(int argc, char **argv) {

   input_data_t d;
   d.is_debug_mode = 0;
   d.is_good = 0;
   d.given_map_flag = 0;
   d.resno_start = UNSET;
   d.resno_end = UNSET;
   d.map_weight = UNSET;
   d.use_rama_targets = 0;
   d.use_torsion_targets = 0;
   
   int ch;
   int option_index = 0;

   const char *optstr = "i:h:f:p:o:m:1:2:c:w";

   struct option long_options[] = {
      {"pdbin",  1, 0, 0}, 
      {"hklin",  1, 0, 0}, 
      {"f",      1, 0, 0}, 
      {"phi",    1, 0, 0}, 
      {"pdbout", 1, 0, 0}, 
      {"mapin",  1, 0, 0}, 
      {"resno-start", 1, 0, 0}, 
      {"resno-end",   1, 0, 0}, 
      {"chain-id",    1, 0, 0},
      {"weight",    1, 0, 0},
      {"rama",      0, 0, 0},
      {"torsions",      0, 0, 0},
      {"torsion",      0, 0, 0},
      {"debug",      0, 0, 0},  // developer option
      {0, 0, 0, 0}
   };

   while (-1 !=
	  (ch = getopt_long(argc, argv, optstr, long_options, &option_index))) {

      switch (ch) { 
      case 0:
	 if (optarg) {
	    std::string arg_str = long_options[option_index].name;

	    if (arg_str == "pdbin") { 
	       d.input_pdb_file_name = optarg;
	    } 
	    if (arg_str == "pdbout") { 
	       d.output_pdb_file_name = optarg;
	    } 
	    if (arg_str == "hklin") { 
	       d.mtz_file_name = optarg;
	    } 
	    if (arg_str == "f") { 
	       d.f_col = optarg;
	    } 
	    if (arg_str == "phi") { 
	       d.phi_col = optarg;
	    } 
	    if (arg_str == "mapin") { 
	       d.map_file_name = optarg;
	    }
	    if (arg_str == "chain-id") { 
	       d.chain_id = optarg;
	    }
	    if (arg_str == "chain") { 
	       d.chain_id = optarg;
	    }
	    if (arg_str == "resno-start") { 
	       d.resno_start = atoi(optarg);
	    }
	    if (arg_str == "resno-end") { 
	       d.resno_end = atoi(optarg);
	    }
	    if (arg_str == "weight") { 
	       d.map_weight = atof(optarg);
	    }
	 } else {
	    // long argument without parameter:
	    std::string arg_str(long_options[option_index].name);
	 
	    if (arg_str == "rama") {
	       d.use_rama_targets = 1;
	    }
	    if (arg_str == "torsions") {
	       d.use_torsion_targets = 1;
	    }
	    if (arg_str == "torsion") {
	       d.use_torsion_targets = 1;
	    }
	    if (arg_str == "debug") {
	       d.is_debug_mode = 1;
	    }
	 } 
	 break;

      case 'i':
	 d.input_pdb_file_name = optarg;
	 break;

      case 'o':
	 d.output_pdb_file_name = optarg;
	 break;

      case 'h':
	 d.mtz_file_name = optarg;
	 break;

      case 'f':
	 d.f_col = optarg;
	 break;

      case 'p':
	 d.phi_col = optarg;
	 break;

      case 'm':
	 d.map_file_name = optarg;
	 d.given_map_flag = 1;
	 break;

      case '1':
	 d.resno_start = atoi(optarg);
	 break;

      case '2':
	 d.resno_end = atoi(optarg);
	 break;

      case 'c':
	 d.chain_id = optarg;
	 break;

      case 'w':
	 d.map_weight = atof(optarg);
	 break;
      }
   }

   if (d.input_pdb_file_name != "" && d.output_pdb_file_name != "") { 
      if (d.resno_start != UNSET && d.resno_end != UNSET) { 
	 if (d.chain_id != "") { 
	    if ( (d.map_file_name != "") || (d.mtz_file_name != "" && d.f_col != "" && d.phi_col != "") ) { 
	       d.is_good = 1;
	    } else {
	       std::cout << "error in input map or [HKLIN, f, or phi]" << std::endl;
	    }
	 } else {
	    std::cout << "error in input chain-id" << std::endl;
	 }
      } else {
	 std::cout << "error in input resno start and end" << std::endl;
      }
   } else {
      std::cout << "error in input pdb or out put pdb file names" << std::endl;
   } 
	    
   return d;
} 


