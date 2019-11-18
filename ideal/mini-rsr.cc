/* ideal/with-map.cc
 * 
 * Copyright 2002, 2003 The University of York
 * Author: Paul Emsley
 * Copyright 2007, 2008 The University of Oxford
 * Copyright 2013, 2015, 2016 by Medical Research Council
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
#include "compat/coot-getopt.h"
#undef __GNU_LIBRARY__
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

#define UNSET mmdb::MinInt4

#ifndef  __MMDB_MMCIF__
#include <mmdb2/mmdb_utils.h>
#endif
  
#include <iostream>
#include <string>
#include <vector>
#include <mmdb2/mmdb_manager.h>

#include "utils/coot-utils.hh"
#include "coords/mmdb-extras.h"
#include "coords/mmdb.h"

#include "simple-restraint.hh"
#include "crankshaft.hh"

#include "model-bond-deltas.hh"

// for debugging (the density_around_point function)
#include "coot-utils/coot-map-utils.hh"

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
   input_data_t() {
      is_good = false;
      is_debug_mode = false;
      given_map_flag = false;
      bond_length_deltas = false;
      radius = 4.2;
      residues_around = mmdb::MinInt4;
      resno_start = mmdb::MinInt4;
      resno_end = mmdb::MinInt4;
      f_col = "FWT";
      phi_col = "PHWT";
      use_rama_targets = false;
      use_planar_peptide_restraints = true;
      use_trans_peptide_restraints = true;
      use_torsion_targets = false;
      tabulate_distortions_flag = false;
      no_refine = false; // if we want distortions only, say
      correlations = false;
      do_crankshaft = false;
      crankshaft_n_peptides = 7;
   }
   bool is_good;
   bool is_debug_mode;
   bool bond_length_deltas;
   bool given_map_flag;
   bool use_rama_targets;
   bool use_torsion_targets;
   bool use_planar_peptide_restraints;
   bool use_trans_peptide_restraints;
   bool tabulate_distortions_flag;
   bool no_refine;
   bool correlations;
   bool do_crankshaft;
   int resno_start;
   int resno_end;
   int residues_around;
   int crankshaft_n_peptides;
   float map_weight; 
   float radius; 
   std::string chain_id;
   std::string mtz_file_name;
   std::string f_col;
   std::string phi_col;
   std::string map_file_name;
   std::string  input_pdb_file_name;
   std::string output_pdb_file_name;
   std::vector<std::string> dictionary_file_names;
   std::vector<std::string> extra_restraints_file_names;
   std::vector<coot::residue_spec_t> single_residue_specs; // not used as yet.
};

input_data_t get_input_details(int argc, char **argv);


std::vector<std::pair<bool,mmdb::Residue *> >
fill_residues(const std::string &chain_id, int resno_start, int resno_end, mmdb::Manager *mol) {

   std::vector<std::pair<bool,mmdb::Residue *> > v;
   
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   mmdb::Chain *chain_p;
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      bool done = false;
      chain_p = model_p->GetChain(ichain);
      std::string this_chain_id = chain_p->GetChainID();
      if (this_chain_id == chain_id) { 
	 int nres = chain_p->GetNumberOfResidues();
	 mmdb::Residue *residue_p;
	 mmdb::Atom *at;
	 for (int ires=0; ires<nres; ires++) {
	    residue_p = chain_p->GetResidue(ires);
	    int this_res_no = residue_p->GetSeqNum();
	    if (this_res_no >= resno_start) {
	       if (this_res_no <= resno_end) { 
		  std::pair<bool, mmdb::Residue *> p(false, residue_p);
		  v.push_back(p);
		  done = true;
	       }
	    }

	    if (!done) {
	       if ((resno_start == mmdb::MinInt4) && (resno_end == mmdb::MinInt4)) {
		  std::pair<bool, mmdb::Residue *> p(false, residue_p);
		  v.push_back(p);
	       }
	    }
	 }
      }
   }
   return v;
}

std::vector<std::pair<bool,mmdb::Residue *> >
fill_residues_all_atoms(mmdb::Manager *mol) {

   std::vector<std::pair<bool, mmdb::Residue *> > v;
   int imod = 1;
   mmdb::Model *model_p = mol->GetModel(imod);
   mmdb::Chain *chain_p;
   int nchains = model_p->GetNumberOfChains();
   for (int ichain=0; ichain<nchains; ichain++) {
      chain_p = model_p->GetChain(ichain);
      if (chain_p) {
	 int nres = chain_p->GetNumberOfResidues();
	 mmdb::Residue *residue_p;
	 for (int ires=0; ires<nres; ires++) {
	    residue_p = chain_p->GetResidue(ires);
	    std::pair<bool, mmdb::Residue *> p(false, residue_p);
	    v.push_back(p);
	 }
      }
   }
   return v;
}

void
print_version() {

   std::cout << "coot-mini-rsr version " << VERSION << std::endl;

}

void show_usage() {

   std::string prog_name = "coot-mini-rsr";
   
   std::cout << "Usage: " << prog_name << "\n"
	     << "       --pdbin pdb-in-filename\n"
	     << "       --hklin mtz-filename\n"
	     << "       --dictin cif-dictionary-filename\n"
	     << "       --f f_col_label\n"
	     << "       --phi phi_col_label\n"
	     << "       --pdbout output-filename\n"
	     << "       [ --resno-start resno_low\n"
	     << "         --resno-end   resno_high] \n"
	     << "       or [ --residues-around res_no ]\n"
	     << "       or [ --residue-number resno ]\n"
	     << "       or [ --res-no resno ]\n"
	     << "       --chain-id chain-id\n"
	     << "       --weight w (weight of map gradients, default 60)\n"
	     << "       --radius (default 4.2)\n"
	     << "       --rama\n"
	     << "       --torsions\n"
	     << "       --no-planar-peptide-restraints\n"
	     << "       --no-trans-peptide-restraints\n"
	     << "       --tabulate-distortions\n"
	     << "       --extra-restraints extra-restraints-file-name\n"
	     << "       --bond-length-deltas\n"
	     << "       --correlations\n"
	     << "       --crankshaft\n"
	     << "       --no-refine\n"
	     << "       --version\n"
	     << "       --debug\n"
	     << "\n"
	     << "     --mapin ccp4-map-name can be used\n"
	     << "       instead of --hklin --f --phi\n"
	     << std::endl;
}

void
execute_crankshaft(const coot::residue_spec_t &rs, int n_peptides, const clipper::Xmap<float> &xmap,
	   mmdb::Manager *mol, float map_weight, int n_samples, const std::string &pdb_out_file_name) {

   int n_solutions = 1; // just the best

   ctpl::thread_pool thread_pool;
   int n_threads = 4; // fixme
   std::vector<mmdb::Manager *> v =
      coot::crankshaft::crank_refine_and_score(rs, n_peptides, xmap, mol, map_weight, n_samples,
					       n_solutions, &thread_pool, n_threads);

   if (v.size() == 1) {
      int err = v[0]->WritePDBASCII(pdb_out_file_name.c_str());
      if (! err)
	 std::cout << "INFO:: wrote " << pdb_out_file_name << std::endl;
   } else {
      std::cout << "WARNING:: No crankshaft solutions" << std::endl;
   }

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
      show_usage();
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
      // molecule (actually, a pointer to it) to do some atom selection.
      // 
//       coot::restraints_container_t restraints(asc.atom_selection, // moving_residue_atoms,
// 					      asc.n_selected_atoms,
// 					      previous_residue,
// 					      next_residue,
// 					      fixed_atoms);
      
      input_data_t inputs = get_input_details(argc, argv);

      // Testing/debugging command line parsing
      // std::cout << "Here with good inputs " << std::endl;
      // inputs.input_pdb_file_name = "../src/test.pdb";
      // inputs.bond_length_deltas = true;
      // inputs.is_good = true;

      if (inputs.is_good) {

	 std::string pdb_file_name(inputs.input_pdb_file_name);
	 bool map_is_good = false; // currently

	 // if pdb_file_name does not exist -> crash?
	 atom_selection_container_t asc = get_atom_selection(pdb_file_name, true, true);

	 if (! asc.read_success) {
	    std::cout << "coordinates read failure " << pdb_file_name << std::endl;
	 }

	 if (inputs.bond_length_deltas) {
	    int imol = 0; // dummy
	    std::cout << "----------- bond length delta --------" << std::endl;
	    coot::model_bond_deltas mbd(asc.mol, imol, &geom);
	    mbd.resolve(); // not a good function name
	    return 0;
	 }

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

	    if (coot::file_exists(inputs.map_file_name)) {
	       // should we try/catch around this read?
	       clipper::CCP4MAPfile file;
	       file.open_read(inputs.map_file_name);
	       file.import_xmap(xmap);
	       file.close_read();
	       map_is_good = true;
	    } else {
	       std::cout << "WARNING:: map file " << inputs.map_file_name << " not found"
			 << std::endl;
	    }
	 }

	 float map_weight = 60.0;
	 if (inputs.map_weight > 0) { 
	    map_weight = inputs.map_weight;
	    std::cout << "INFO:: using weight " << map_weight << std::endl;
	 }

	 if (map_is_good) {
	    std::string altloc("");

	    // load up the protein_geometry with user-supplied dictionaries
	    int read_number = 55;
	    for (std::size_t i=0; i<inputs.dictionary_file_names.size(); i++) {
	       coot::read_refmac_mon_lib_info_t mli = geom.init_refmac_mon_lib(inputs.dictionary_file_names[i],
									       read_number++,
									       coot::protein_geometry::IMOL_ENC_ANY);
	    }


	    coot::extra_restraints_t extra_restraints;
	    for (unsigned int jj=0; jj<inputs.extra_restraints_file_names.size(); jj++) {
	       const std::string &file_name = inputs.extra_restraints_file_names[jj];
	       coot::extra_restraints_t r;
	       r.read_refmac_extra_restraints(file_name);
	       extra_restraints.add_restraints(r);
	    }

	    // the bool is to annotate for "fixed" residue - here nothing is fixed.
	    // 
	    std::vector<std::pair<bool,mmdb::Residue *> > local_residues;

	    if (inputs.chain_id.empty())
	       local_residues = fill_residues_all_atoms(asc.mol);
	    else
	       local_residues = fill_residues(chain_id, inputs.resno_start, inputs.resno_end, asc.mol);

	    if (inputs.residues_around != mmdb::MinInt4) {
	       coot::residue_spec_t res_spec(inputs.chain_id, inputs.residues_around);
	       mmdb::Residue *res_ref = coot::util::get_residue(res_spec, asc.mol);
	       if (res_ref) {
		  std::vector<mmdb::Residue *> residues =
		     coot::residues_near_residue(res_ref, asc.mol, inputs.radius);
		  local_residues.push_back(std::pair<bool,mmdb::Residue *>(false, res_ref));
		  for (unsigned int i=0; i<residues.size(); i++) {
		     std::pair<bool, mmdb::Residue *> p(false, residues[i]);
		     local_residues.push_back(p);
		  }
	       }
	    }

	    // now, for the correlations, make the residue specs from
	    // the local_residues (somewhat convoluted, I agree)
	    // 
	    std::vector<coot::residue_spec_t> residue_specs;
	    for (unsigned int ii=0; ii<local_residues.size(); ii++) { 
	       if (local_residues[ii].second) {
		  coot::residue_spec_t spec(local_residues[ii].second);
		  residue_specs.push_back(spec);
	       }
	    }

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

	    // 20150919: Request from Robbie Joosten: add in the planar peptide restraints.
	    //           (this should be user-configurable).
	    // 
	    if (inputs.use_planar_peptide_restraints)
	       geom.add_planar_peptide_restraint();

	    if (inputs.do_crankshaft) {
	       coot::residue_spec_t res_spec(inputs.chain_id, inputs.resno_start);
	       execute_crankshaft(res_spec, inputs.crankshaft_n_peptides, xmap, asc.mol,
			  inputs.map_weight, -1, inputs.output_pdb_file_name);
	    } else {

	       std::vector<mmdb::Link> links;
	       coot::restraints_container_t restraints(local_residues,
						       links,
						       geom,
						       asc.mol,
						       fixed_atom_specs, &xmap);
      
	       restraints.add_map(map_weight);

	       // coot::restraint_usage_Flags flags = coot::NO_GEOMETRY_RESTRAINTS;
	       // coot::restraint_usage_Flags flags = coot::BONDS;
	       // coot::restraint_usage_Flags flags = coot::BONDS_AND_ANGLES;
	       // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_AND_PLANES;
	       // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_TORSIONS_AND_PLANES;
	       // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_AND_NON_BONDED;
	       // flags = coot::NON_BONDED;
	       coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
	       flags = coot::TYPICAL_RESTRAINTS;

	       if (inputs.use_torsion_targets) {
		  flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_AND_CHIRALS;
		  if (inputs.use_rama_targets)
		     flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_RAMA;
	       } else {
		  if (inputs.use_rama_targets)
		     flags = coot::BONDS_ANGLES_TORSIONS_PLANES_NON_BONDED_CHIRALS_AND_RAMA;
	       }

	       coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
	       bool do_rama_plot_restraints = 0;
	       if (inputs.use_rama_targets)
		  do_rama_plot_restraints = 1;
	       // this should be a user-settable parameter.
	       bool make_trans_peptide_restraints = false;

	       if (inputs.use_trans_peptide_restraints)
		  make_trans_peptide_restraints = true;

	       int n_threads = coot::get_max_number_of_threads();
	       ctpl::thread_pool thread_pool(n_threads);
	       restraints.thread_pool(&thread_pool, n_threads);

	       int imol = 0; // dummy
	       restraints.make_restraints(imol, geom, flags, 1, make_trans_peptide_restraints,
	       				  1.0, do_rama_plot_restraints, true, true, pseudos);

               // std::vector<std::pair<bool,mmdb::Residue *> > residues;
               // coot::protein_geometry geom;
               // mmdb::Manager *mol = 0;
               // clipper::Xmap<float> xmap; // don't go out of scope until refinement is over
               // coot::restraints_container_t restraints_2(residues, geom, mol, &xmap);

	       if (inputs.extra_restraints_file_names.size() > 0) {
		  restraints.add_extra_restraints(imol, "user-defined restraints", extra_restraints, geom);
	       }

	       if (! inputs.no_refine) {
		  int nsteps_max = 4000;
		  short int print_chi_sq_flag = 1;
		  restraints.minimize(flags, nsteps_max, print_chi_sq_flag);
		  restraints.write_new_atoms(inputs.output_pdb_file_name);
	       }

	       if (inputs.tabulate_distortions_flag) {
		  coot::geometry_distortion_info_container_t gd = restraints.geometric_distortions();
		  gd.print();
	       }

	       if (inputs.correlations) {

		  std::vector<coot::residue_spec_t> neighbours;
		  int ATOM_MASK_MAINCHAIN           =  1; // hideous hack
		  int ATOM_MASK_NOT_MAINCHAIN       =  2;
		  int ATOM_MASK_NOT_MAINCHAIN_OR_CB =  3;
		  int ATOM_MASK_ALL_ATOM_B_FACTOR   = 10;
		  unsigned short int atom_mask_mode = ATOM_MASK_ALL_ATOM_B_FACTOR;
		  coot::map_stats_t map_stats_flag = coot::SIMPLE;
	       
		  coot::util::density_correlation_stats_info_t stats =
		     coot::util::map_to_model_correlation_stats(asc.mol,
								residue_specs,
								neighbours,
								atom_mask_mode,
								2.5, // dummy for this mode
								xmap,
								map_stats_flag);

		  std::vector<std::pair<coot::residue_spec_t, float> > correls =
		     coot::util::map_to_model_correlation_per_residue(asc.mol,
								      residue_specs,
								      atom_mask_mode,
								      2.5, // dummy
								      xmap);


		  std::cout << " Residue Correlation Table: " << std::endl;
		  for (unsigned int j=0; j<correls.size(); j++) {
		     std::string res_name;
		     mmdb::Residue *r = coot::util::get_residue(correls[j].first, asc.mol);
		     if (r) res_name = r->GetResName();
		     std::cout << "     "
			       << correls[j].first.chain_id << " "
			       << correls[j].first.res_no   << " "
			       << correls[j].first.ins_code << " "
			       << res_name
			       << "  " << correls[j].second << std::endl;
		  }
	       }
	    } 
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
      std::cout << "reading mtz file..." << std::endl;
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
	 std::cout << "We should use the weights: " << weight_col << std::endl;
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
      std::cout << "Grid..." << xmap.grid_sampling().format() << "\n";
      std::cout << "doing fft..." << std::endl;

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
      std::cout << "done fft..." << std::endl;
      status = 1;
   }
   catch (const clipper::Message_base &exc) {  // "exception" is a protected word, it seems.
      std::cout << "Failed to read mtz file " << mtz_file_name << std::endl;
   }

   return std::pair<bool, clipper::Xmap<float> > (status, xmap);
}


input_data_t
get_input_details(int argc, char **argv) {

   input_data_t d;
   d.is_debug_mode = 0;
   d.is_good = 0;
   d.bond_length_deltas = false;
   d.given_map_flag = 0;
   d.resno_start = UNSET;
   d.resno_end   = UNSET;
   d.map_weight  = UNSET;
   d.use_rama_targets = 0;
   d.use_torsion_targets = 0;
   d.use_planar_peptide_restraints = true;
   d.use_trans_peptide_restraints = true;
   d.correlations = false;
   
   int ch;
   int option_index = 0;

   const char *optstr = "i:h:f:p:o:m:1:2:c:w";

   struct option long_options[] = {
      {"version",0, 0, 0},
      {"pdbin",  1, 0, 0},
      {"hklin",  1, 0, 0},
      {"f",      1, 0, 0},
      {"phi",    1, 0, 0},
      {"pdbout", 1, 0, 0},
      {"dictin", 1, 0, 0},
      {"dictionary", 1, 0, 0},
      {"mapin",  1, 0, 0},
      {"extra-restraints",  1, 0, 0},
      {"resno-start", 1, 0, 0},
      {"resno-end",   1, 0, 0},
      {"residue-number",    1, 0, 0}, // for a single residue refinement of crankshafting
      {"resno",             1, 0, 0}, // same
      {"res-no",            1, 0, 0}, // same
      // {"residue",     2, 0, 0}, // can't give 2 args to one keyword, it seems?
      {"residues-around",   1, 0, 0},
      {"chain-id",    1, 0, 0},
      {"weight",    1, 0, 0},
      {"radius",    1, 0, 0},
      {"rama",      0, 0, 0},
      {"torsions",  0, 0, 0},
      {"torsion",   0, 0, 0},
      {"no-planar-peptide-restraints", 0, 0, 0},
      {"no-trans-peptide-restraints",  0, 0, 0},
      {"tabulate-distortions", 0, 0, 0},
      {"bond-length-deltas", 0, 0, 0},
      {"no-refine", 0, 0, 0},
      {"correlations", 0, 0, 0},
      {"crankshaft", 0, 0, 0},
      {"n-peptides", 1, 0, 0},
      {"crank", 0, 0, 0},
      {"help", 0, 0, 0},
      {"debug",     0, 0, 0},  // developer option
      {0, 0, 0, 0}
   };

   while (-1 !=
	  (ch = coot_getopt_long(argc, argv, optstr, long_options, &option_index))) {

      // std::cout << "here 1 ch " << ch << std::endl;

      switch (ch) {
      case 0:
	 if (! optarg) std::cout << "No optarg" << std::endl;
	 if (optarg) {

	    // std::cout << "here 2 optarg " << optarg << std::endl;

	    std::string arg_str = long_options[option_index].name;

	    if (arg_str == "pdbin") {
	       // std::cout << "found pdbin " << optarg << std::endl;
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
               d.given_map_flag = 1;
	    }
	    if (arg_str == "dictin") {
	       d.dictionary_file_names.push_back(optarg);
	    }
	    if (arg_str == "dictionary") {
	       d.dictionary_file_names.push_back(optarg);
	    }
	    if (arg_str == "extra-restraints") {
	       d.extra_restraints_file_names.push_back(optarg);
	    }
	    if (arg_str == "chain-id") {
	       d.chain_id = optarg;
	    }
	    if (arg_str == "chain") {
	       d.chain_id = optarg;
	    }
	    if (arg_str == "resno-start") {
	       d.resno_start = coot::util::string_to_int(optarg);
	    }
	    if (arg_str == "resno-end") {
	       d.resno_end = coot::util::string_to_int(optarg);
	    }
	    if (arg_str == "resno" || arg_str == "residue-number" || arg_str == "res-no") {
	       try {
		  d.resno_start = coot::util::string_to_int(optarg);
		  d.resno_end = d.resno_start;
	       }
	       catch (const std::runtime_error &rte) {
		  std::cout << "WARNING::" << rte.what() << std::endl;
	       }
	    }
	    if (arg_str == "n-peptides") {
	       try {
		  d.crankshaft_n_peptides = coot::util::string_to_int(optarg);
	       }
	       catch (const std::runtime_error &rte) {
		  std::cout << "WARNING::" << rte.what() << std::endl;
	       }
	    }
	    if (arg_str == "residues-around") {
	       try {
		  d.residues_around = coot::util::string_to_int(optarg);
	       }
	       catch (const std::runtime_error &rte) {
		  std::cout << "WARNING::" << rte.what() << std::endl;
	       }
	    }
	    if (arg_str == "weight") {
	       d.map_weight = coot::util::string_to_float(optarg);
	    }
	    if (arg_str == "radius") {
	       d.radius = coot::util::string_to_float(optarg);
	    }
	 } else {
	    // long argument without parameter:
	    std::string arg_str(long_options[option_index].name);
	 
	    if (arg_str == "help") {
	       show_usage();
	       exit(0);
	    }
	    if (arg_str == "version") {
	       print_version();
	       exit(0);
	    }
	    if (arg_str == "bond-length-deltas") {
	       d.bond_length_deltas = true;
	    }
	    if (arg_str == "rama") {
	       d.use_rama_targets = 1;
	    }
	    if (arg_str == "crank") {
	       d.do_crankshaft = 1;
	    }
	    if (arg_str == "crankshaft") {
	       d.do_crankshaft = 1;
	    }
	    if (arg_str == "torsions") {
	       d.use_torsion_targets = 1;
	    }
	    if (arg_str == "torsion") {
	       d.use_torsion_targets = 1;
	    }
	    if (arg_str == "no-planar-peptide-restraints") {
	       d.use_planar_peptide_restraints = false;
	    }
	    if (arg_str == "no-trans-peptide-restraints") {
	       d.use_trans_peptide_restraints = false;
	    }
	    if (arg_str == "tabulate-distortions") {
	       d.tabulate_distortions_flag = true;
	    }
	    if (arg_str == "no-refine") {
	       d.no_refine = true;
	    }
	    if (arg_str == "correlations") {
	       d.correlations = true;
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

      case 'v':
	 print_version();
	 exit(0);
	 break;

      case 'w':
	 d.map_weight = atof(optarg);
	 break;
      }
   }

   if (!d.input_pdb_file_name.empty()){
      if (!d.output_pdb_file_name.empty()) {
	 if ((d.resno_start != UNSET && d.resno_end != UNSET) || (d.residues_around != mmdb::MinInt4)) {
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
	 std::cout << "error in output pdb file name" << std::endl;
      }
   } else {
      std::cout << "error in input pdb file name" << std::endl;
   }
	    
   return d;
}
