

#ifdef __GNU_LIBRARY__
#include "coot-getopt.h"
#else
#define __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#undef __GNU_LIBRARY__
#endif

#include <iomanip>

#include <mmdb2/mmdb_manager.h>
#include "coords/mmdb-extras.h"
#include "coords/mmdb.h"

#include "crankshaft.hh"

#include "clipper/core/xmap.h"
#include "clipper/core/map_utils.h"
#include "clipper/ccp4/ccp4_map_io.h"
#include "clipper/ccp4/ccp4_mtz_io.h"
#include "clipper/core/hkl_compute.h"
#include "clipper/core/map_utils.h" // Map_stats

#include "simple-restraint.hh"
#include "coot-utils/coot-map-utils.hh"

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
      n_samples = -1;
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
      {"n-samples", 1, 0, 0},
      {"radius",    1, 0, 0},
      {"rama",      0, 0, 0},
      {"help", 0, 0, 0},
      {"debug",     0, 0, 0},  // developer option
      {0, 0, 0, 0}
   };

   while (-1 !=
	  (ch = getopt_long(argc, argv, optstr, long_options, &option_index))) {

      switch (ch) {
      case 0:
	 if (optarg) {
	    std::string arg_str = long_options[option_index].name;

	    if (arg_str == "pdbin") {
	       inputs.input_pdb_file_name = optarg;
	    }
	    if (arg_str == "pdbout") {
	       inputs.output_pdb_file_name = optarg;
	    }
	    if (arg_str == "mapin") {
	       inputs.map_file_name = optarg;
	    }
	    if (arg_str == "hklin") {
	       inputs.mtz_file_name = optarg;
	    }
	    if (arg_str == "chain-id") {
	       inputs.chain_id = optarg;
	    }
	    if (arg_str == "resno") {
	       inputs.res_no = coot::util::string_to_int(optarg);
	    }
	    if (arg_str == "res-no") {
	       inputs.res_no = coot::util::string_to_int(optarg);
	    }
	    if (arg_str == "weight") {
	       inputs.map_weight = coot::util::string_to_float(optarg);
	    }
	    if (arg_str == "n-samples") {
	       try {
		  inputs.n_samples = coot::util::string_to_int(optarg);
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

      if (use_weights) {
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



int main(int argc, char **argv) {

   int status = 0;

   input_data_t inputs = get_input_details(argc, argv);
   atom_selection_container_t asc = get_atom_selection(inputs.input_pdb_file_name, 1, 0);
   if (! asc.read_success) {
      std::cout << "fail on read pdb" << inputs.input_pdb_file_name << std::endl;
   } else {

      // this is the middle residue
      coot::residue_spec_t rs(inputs.chain_id, inputs.res_no);

      // happy path

      float map_weight = 60;
      unsigned int n_peptides = 7;
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
	       int n_solutions = 6;
	       std::vector<mmdb::Manager *> mols =
		  coot::crankshaft::crank_refine_and_score(rs, n_peptides, xmap, asc.mol, map_weight, inputs.n_samples, n_solutions);
	       for (std::size_t i=0; i<mols.size(); i++) {
		  std::string file_name = "crankshafted-";
		  file_name += coot::util::int_to_string(i);
		  file_name += ".pdb";
		  mols[i]->WritePDBASCII(file_name.c_str());
	       }
	       for (std::size_t i=0; i<mols.size(); i++)
		  delete mols[i];
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
	    int n_solutions = 6;
	    std::vector<mmdb::Manager *> mols =
	       coot::crankshaft::crank_refine_and_score(rs, n_peptides, xmap, asc.mol,
							map_weight, inputs.n_samples, n_solutions);
	    for (std::size_t i=0; i<mols.size(); i++) {
	       std::string file_name = "crankshafted-";
		  file_name += coot::util::int_to_string(i);
		  file_name += ".pdb";
		  mols[i]->WritePDBASCII(file_name.c_str());
	    }
	    for (std::size_t i=0; i<mols.size(); i++)
	       delete mols[i];
	 } else {
	    std::cout << "ERROR:: bad map " << std::endl;
	 }
      }
   }
   return status;
};
