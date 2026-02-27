/*
 * ideal/crankshaft-test.cc
 *
 * Copyright 2020 by Medical Research Council
 * Author: Paul Emsley
 *
 * This file is part of Coot
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 */



#ifdef __GNU_LIBRARY__
#include "coot-getopt.h"
#else
#define __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#undef __GNU_LIBRARY__
#endif

#include <iomanip>

#include <mmdb2/mmdb_manager.h>

#include <clipper/core/xmap.h>
#include <clipper/ccp4/ccp4_map_io.h>
#include <clipper/ccp4/ccp4_mtz_io.h>
#include <clipper/core/hkl_compute.h>

#include "coot-utils/coot-map-utils.hh"

#include "simple-restraint.hh"
#include "crankshaft.hh"
#include "utils/logging.hh"
extern logging logger;

class molecule_score_t {
public:
   molecule_score_t() { density_score = -1; model_score = -1;}
   float density_score; // ~50, big is good - needs thought with cryo-EM maps
   float model_score;   // ~50, big is bad

};

// abstract out for use with min-rsr also.

class input_data_t {
public:
   input_data_t() {
      is_good = false;
      given_map_flag = false;
      radius = 4.2;
      residues_around = mmdb::MinInt4;
      res_no = mmdb::MinInt4;
      f_col = "FWT";
      phi_col = "PHWT";
      use_rama_targets = false;
      use_planar_peptide_restraints = true;
      use_trans_peptide_restraints = true;
      use_torsion_targets = false;
      tabulate_distortions_flag = false;
      n_samples = 30;
      map_weight = -1;
   }
   bool is_good;
   bool given_map_flag;
   bool use_rama_targets;
   bool use_torsion_targets;
   bool use_planar_peptide_restraints;
   bool use_trans_peptide_restraints;
   bool tabulate_distortions_flag;
   int res_no;
   int residues_around;
   int n_samples;
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
   std::vector<coot::residue_spec_t> single_residue_specs; // not used as yet.
};

input_data_t
get_input_details(int argc, char **argv) {

   input_data_t inputs;

   int ch;
   int option_index = 0;

   const char *optstr = "i:h:f:p:o:m:1:2:c:w";

   struct option long_options[] = {
      {"pdbin",  1, 0, 0},
      {"hklin",  1, 0, 0},
      {"f",      1, 0, 0},
      {"phi",    1, 0, 0},
      {"pdbout", 1, 0, 0},
      {"dictin", 1, 0, 0},
      {"dictionary", 1, 0, 0},
      {"mapin",  1, 0, 0},
      {"res-no", 1, 0, 0},
      {"residues-around",   1, 0, 0},
      {"chain-id",    1, 0, 0},
      {"weight",    1, 0, 0},
      {"n_samples", 1, 0, 0},
      {"radius",    1, 0, 0},
      {"rama",      0, 0, 0},
      {"help", 0, 0, 0},
      {"debug",     0, 0, 0},  // developer option
      {0, 0, 0, 0}
   };

   while (-1 !=
	  (ch = coot_getopt_long(argc, argv, optstr, long_options, &option_index))) {

      switch (ch) {
      case 0:
	 if (coot_optarg) {
	    std::string arg_str = long_options[option_index].name;

	    if (arg_str == "pdbin") {
	       inputs.input_pdb_file_name = coot_optarg;
	    }
	    if (arg_str == "pdbout") {
	       inputs.output_pdb_file_name = coot_optarg;
	    }
	    if (arg_str == "mapin") {
	       inputs.map_file_name = coot_optarg;
	    }
	    if (arg_str == "hklin") {
	       inputs.mtz_file_name = coot_optarg;
	    }
	    if (arg_str == "chain-id") {
	       inputs.chain_id = coot_optarg;
	    }
	    if (arg_str == "resno") {
	       inputs.res_no = coot::util::string_to_int(coot_optarg);
	    }
	    if (arg_str == "res-no") {
	       inputs.res_no = coot::util::string_to_int(coot_optarg);
	    }
	    if (arg_str == "weight") {
	       inputs.map_weight = coot::util::string_to_float(coot_optarg);
	    }
	    if (arg_str == "n-samples") {
	       try {
		  inputs.map_weight = coot::util::string_to_int(coot_optarg);
	       }
	       catch (const std::runtime_error &rte) {
		  std::cout << "WARNING:: parsing n-samples " << rte.what() << std::endl;
	       }
	    }
	 }
      }
   }
   return inputs;
}

molecule_score_t
refine_and_score_mol(mmdb::Manager *mol,
		     const coot::residue_spec_t &res_spec_mid,
		     const std::vector<coot::residue_spec_t> &residue_specs_for_scoring,
		     const coot::protein_geometry &geom,
		     const clipper::Xmap<float> &xmap,
		     float map_weight,
		     const std::string &output_pdb_file_name) {

   molecule_score_t score;

   if (mol) {
      float radius = 2.2; // for the moment
      std::vector<std::pair<bool, mmdb::Residue *> > residues;
      std::vector<mmdb::Link> links;
      std::vector<coot::atom_spec_t> fixed_atom_specs;
      coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_NON_BONDED_AND_CHIRALS;
      flags = coot::TYPICAL_RESTRAINTS;
      int imol = 0; // dummy
      int nsteps_max = 4000;
      bool make_trans_peptide_restraints = true;
      short int print_chi_sq_flag = 1;
      bool do_rama_plot_restraints = true;
      coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
      int restraints_rama_type = coot::restraints_container_t::RAMA_TYPE_ZO;

      mmdb::Residue *res_ref = coot::util::get_residue(res_spec_mid, mol);
      if (res_ref) {
	 std::vector<mmdb::Residue *> residues_near = coot::residues_near_residue(res_ref, mol, radius);
	 residues.push_back(std::pair<bool,mmdb::Residue *>(true, res_ref));
	 for (unsigned int i=0; i<residues_near.size(); i++) {
	    std::pair<bool, mmdb::Residue *> p(true, residues_near[i]);
	    residues.push_back(p);
	 }
      }

      std::vector<coot::residue_spec_t> residue_specs_vec(residues.size());
      for (std::size_t ir=0; ir<residues.size(); ir++)
	 residue_specs_vec[ir] = coot::residue_spec_t(residues[ir].second);

      auto tp_0 = std::chrono::high_resolution_clock::now();
      coot::restraints_container_t restraints(residues, links, geom, mol, fixed_atom_specs, &xmap);
      restraints.add_map(map_weight);
      restraints.set_rama_type(restraints_rama_type);
      restraints.set_rama_plot_weight(1);
      restraints.make_restraints(imol, geom, flags, 1, make_trans_peptide_restraints,
				 1.0, do_rama_plot_restraints, true, true, false, pseudos);
      restraints.minimize(flags, nsteps_max, print_chi_sq_flag);
      if (! output_pdb_file_name.empty())
	 restraints.write_new_atoms(output_pdb_file_name);

      coot::geometry_distortion_info_container_t gdic = restraints.geometric_distortions();
      for (std::size_t id=0; id<gdic.geometry_distortion.size(); id++) {
	 // std::cout << "   " << gdic.geometry_distortion[id] << std::endl;
      }
      // std::cout << "   total-distortion: " << gdic.print() << std::endl;

      auto tp_1 = std::chrono::high_resolution_clock::now();
      bool mainchain_only_flag = true;
      float score_map = coot::util::map_score_by_residue_specs(mol, residue_specs_for_scoring,
							       xmap, mainchain_only_flag);
      auto tp_2 = std::chrono::high_resolution_clock::now();
      auto d10 = std::chrono::duration_cast<std::chrono::microseconds>(tp_1 - tp_0).count();
      auto d21 = std::chrono::duration_cast<std::chrono::microseconds>(tp_2 - tp_1).count();
      if (true) {
	 std::cout << "refine mol : d10  " << d10 << " microseconds\n";
	 std::cout << "map_score_by_residue_specs: d21  " << d21 << " microseconds\n";
	 std::cout << "scores map: " << score_map << " distortion " << gdic.distortion()
		   << " " << output_pdb_file_name << std::endl;
      }
      score.density_score = score_map;
      score.model_score = gdic.distortion();
   }
   return score;
}



// for use the threads:
void
refine_and_score_mols(std::vector<mmdb::Manager *> mols,
		      const std::vector<unsigned int> &mols_thread_vec,
		      const coot::residue_spec_t &res_spec_mid,
		      const std::vector<coot::residue_spec_t> &residue_specs_for_scoring,
		      const coot::protein_geometry &geom,
		      const clipper::Xmap<float> &xmap,
		      float map_weight,
		      std::vector<molecule_score_t> *mol_scores) {

   for (std::size_t i=0; i<mols_thread_vec.size(); i++) {
      molecule_score_t ms = refine_and_score_mol(mols[mols_thread_vec[i]],
						 res_spec_mid, residue_specs_for_scoring,
						 geom, xmap, map_weight, "");
      (*mol_scores)[mols_thread_vec[i]] = ms;
   }
}

// for testing threads
void
dummy_func(std::vector<mmdb::Manager *> mols) {

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
      xmap.init(myhkl.spacegroup(), myhkl.cell(),
		clipper::Grid_sampling(myhkl.spacegroup(),
				       myhkl.cell(),
				       myhkl.resolution()));
      std::cout << "Grid..." << xmap.grid_sampling().format() << "\n";
      xmap.fft_from( fphidata );                  // generate map
      status = 1;
   }
   catch (const clipper::Message_base &exc) {
      std::cout << "Failed to read mtz file " << mtz_file_name << std::endl;
   }

   return std::pair<bool, clipper::Xmap<float> > (status, xmap);
}

void
crank_refine_and_score(const coot::residue_spec_t &rs, // mid-residue
		       const clipper::Xmap<float> &xmap,
		       atom_selection_container_t asc,
		       float map_weight,
		       int n_samples) {
   
   // cranshaft takes the first of 3 residues, so we want to call cranshaft with
   // the residue before rs:
   mmdb::Residue *prev_res = coot::util::get_previous_residue(rs, asc.mol);
   if (! prev_res) {
      std::cout << "WARNING:: No residue previous to " << rs << std::endl;
   } else {
      coot::crankshaft cs(asc.mol);
      zo::rama_table_set zorts;

      coot::residue_spec_t prev_residue_spec(prev_res);
      std::cout << "INFO:: using residue specifier: " << rs << std::endl;

      // consider changing the API to find_maxima to be the middle of the 3 residues
      //
      //
      std::vector<coot::crankshaft::scored_triple_angle_set_t> sas =
	 cs.find_maxima_from_triples(prev_residue_spec, zorts, n_samples);
      // where did they move to?
      std::cout << "sas size: " << sas.size() << std::endl;

      coot::protein_geometry geom;
      geom.set_verbose(0);
      geom.init_standard();

      if (false) { // debug models
	 for (std::size_t i=0; i<sas.size(); i++) {
	    std::string output_pdb_file_name = "crankshaft-pre-ref-" +
	       coot::util::int_to_string(i) + ".pdb";
	    mmdb::Manager *mol = cs.new_mol_with_moved_atoms(sas[i]); // delete mol
	    mol->WritePDBASCII(output_pdb_file_name.c_str());
	    delete mol;
	 }
      }

      std::vector<coot::residue_spec_t> local_residue_specs;
      mmdb::Residue *res_ref = coot::util::get_residue(rs, asc.mol);
      if (res_ref) {
	 float radius = 5;
	 std::vector<mmdb::Residue *> residues_near =
	    coot::residues_near_residue(res_ref, asc.mol, radius);
	 local_residue_specs.push_back(rs);
	 for (unsigned int i=0; i<residues_near.size(); i++)
	    local_residue_specs.push_back(coot::residue_spec_t(residues_near[i]));
      }

      std::vector<mmdb::Manager *> mols(sas.size(), 0);
      std::vector<molecule_score_t> mol_scores(sas.size());
      for (std::size_t i=0; i<sas.size(); i++) {
	 // delete mol when finished with it
	 mmdb::Manager *mol = cs.new_mol_with_moved_atoms(sas[i]); // d
	 mols[i] = mol;
      }

#ifdef HAVE_CXX_THREAD
      unsigned int n_threads = coot::get_max_number_of_threads();
      auto tp_start = std::chrono::high_resolution_clock::now();

      if (n_threads > 0) {

	 std::vector<std::thread> threads;
	 unsigned int thread_idx = 0;
	 std::cout << "debug:: making mols_thread_vec " << n_threads << std::endl;
	 std::vector<std::vector<unsigned int> > mols_thread_vec(n_threads);
	 std::cout << "debug:: done making mols_thread_vec " << std::endl;
	 // put the mols in mols_thread_vec
	 for (unsigned int ii=0; ii<mols.size(); ii++) {
	    std::cout << "adding mol " << ii << " to thread " << thread_idx << std::endl;
	    mols_thread_vec[thread_idx].push_back(ii);
	    thread_idx++;
	    if (thread_idx == n_threads) thread_idx = 0;
	 }
      
	 std::cout << "debug:: the mols_thread_vecs are: " << std::endl;
	 for (std::size_t j=0; j<mols_thread_vec.size(); j++) {
	    std::cout << "thread " << j << " of " << mols_thread_vec.size() << " threads " << std::endl;
	    for (std::size_t k=0; k<mols_thread_vec[j].size(); k++) {
	       std::cout << " " << mols_thread_vec[j][k];
	    }
	    std::cout << std::endl;
	 }

         unsigned int n_pushed_threads = 0;
	 for (unsigned int i_thread=0; i_thread<n_threads; i_thread++) {

	    auto tp_0 = std::chrono::high_resolution_clock::now();
	    const std::vector<unsigned int> &i_mols = mols_thread_vec[i_thread];

	    // Don't bother to push a thread if there are no molecules to refine
	    // (this happens when there are more threads than there
	    // are crankshaft probability maxima/solutions).
	    //
	    if (mols_thread_vec[i_thread].size() > 0) {
	       threads.push_back(std::thread(refine_and_score_mols, mols, mols_thread_vec[i_thread],
					     std::cref(rs), std::cref(local_residue_specs),
					     std::cref(geom), std::cref(xmap), map_weight,
					     &mol_scores));
               ++n_pushed_threads;
	       // threads.push_back(std::thread(dummy_func, mols)); // testing
	       auto tp_1 = std::chrono::high_resolution_clock::now();
	       auto d10 = std::chrono::duration_cast<std::chrono::microseconds>(tp_1 - tp_0).count();
	       if (true)
		  std::cout << "time to push thread: " << d10 << " microseconds\n";
	    }
	 }


	 auto tp_2 = std::chrono::high_resolution_clock::now();
	 for (unsigned int i_thread=0; i_thread<n_pushed_threads; i_thread++)
	    threads.at(i_thread).join();
	 auto tp_3 = std::chrono::high_resolution_clock::now();
	 auto d32 = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_2).count();
	 auto d3start = std::chrono::duration_cast<std::chrono::milliseconds>(tp_3 - tp_start).count();
	 if (true)
	    std::cout << "time to join threads: " << d32 << " milliseconds\n";
	 if (true)
	    std::cout << "time to do all refinements: " << d3start << " milliseconds\n";

      } else {
	 std::cout << "ERROR:: No threads" << std::endl;
      }
#endif // HAVE_CXX_THREAD      

      for (std::size_t i=0; i<sas.size(); i++) {

	 // the higher the combi_score the better
	 const molecule_score_t &ms = mol_scores[i];
	 float combi_score = 0.01 * map_weight * ms.density_score - ms.model_score - sas[i].minus_log_prob;
	 std::cout << "scores: " << i << " "
		   << std::setw(9) << sas[i].minus_log_prob << " "
		   << std::setw(9) << ms.density_score << " "
		   << std::setw(9) << ms.model_score << " combi-score "
		   << std::setw(9) << combi_score << std::endl;
      }

      for (std::size_t i=0; i<sas.size(); i++)
	 delete mols[i];
      
   }
}



int main(int argc, char **argv) {

   int status = 0;

   input_data_t inputs = get_input_details(argc, argv);
   atom_selection_container_t asc = get_atom_selection(inputs.input_pdb_file_name, false, true, false);
   if (! asc.read_success) {
      std::cout << "fail on read pdb" << inputs.input_pdb_file_name << std::endl;
   } else {

      // this is the middle residue
      coot::residue_spec_t rs(inputs.chain_id, inputs.res_no);

      // happy path

      float map_weight = 60;
      if (inputs.map_weight > 0)
	 map_weight = inputs.map_weight;

      if (inputs.mtz_file_name.empty()) {

	 // were we passed a map?
	 if (! inputs.map_file_name.empty()) {

	    if (coot::file_exists(inputs.map_file_name)) {

	       // should we try/catch around this read?
	       clipper::CCP4MAPfile file;
	       clipper::Xmap<float> xmap;
	       file.open_read(inputs.map_file_name);
	       file.import_xmap(xmap);
	       file.close_read();
	       crank_refine_and_score(rs, xmap, asc, map_weight, inputs.n_samples);
	    } else {
	       std::cout << "ERROR: map " << inputs.map_file_name << " does not exist" << std::endl;
	    }
	 }

      } else {
	 std::pair<bool, clipper::Xmap<float> > xmap_pair =
	    map_from_mtz(inputs.mtz_file_name, inputs.f_col,
			 inputs.phi_col, "", 0, 0, false);
	 const clipper::Xmap<float> &xmap = xmap_pair.second;
	 bool map_is_good = xmap_pair.first;

	 if (map_is_good) {
	    crank_refine_and_score(rs, xmap, asc, map_weight, inputs.n_samples);
	 } else {
	    std::cout << "ERROR:: bad map " << std::endl;
	 }
      }
   }
   return status;
};
