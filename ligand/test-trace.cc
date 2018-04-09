/* ligand/test-trace.cc
 * 
 * Copyright 2016 by Medical Research Council
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

#include <clipper/ccp4/ccp4_map_io.h>

#include "../utils/coot-utils.hh"
#include "../geometry/residue-and-atom-specs.hh"
#include "../coot-utils/coot-coord-utils.hh"
#include "../coot-utils/coot-map-utils.hh"
#include "multi-peptide.hh"
#include "trace.hh"

int main(int argc, char **argv) {

   bool debug = false;

   std::string map_file_name = "trace-test.map";

   if (argc > 1) map_file_name = argv[1];

   if (coot::file_exists(map_file_name)) {

      try {
	 clipper::CCP4MAPfile file;
	 clipper::Xmap<float> xmap;
	 file.open_read(map_file_name);
	 file.import_xmap(xmap);
	 file.close_read();

	 coot::trace t(xmap);

	 if (argc > 3) {

	    // test-trace map-name pdb-name x x x

	    std::string pdb_name = argv[2];
	    mmdb::Manager *mol = new mmdb::Manager;
	    std::cout << "Reading coordinate file: " << pdb_name.c_str() << "\n";
	    mmdb::ERROR_CODE err = mol->ReadCoorFile(pdb_name.c_str());
	    if (err) {
	       std::cout << "There was an error reading " << pdb_name.c_str() << ". \n";
	       std::cout << "ERROR " << err << " READ: "
			 << mmdb::GetErrorDescription(err) << std::endl;
	    } else {

	       // coot::residue_spec_t spec("A", 425, "");

	       // going forward from here in the tutorial model fails
	       // at a glycine - needs investigation (it's not the number of trials)
	       coot::residue_spec_t spec("A", 28, "");

	       if (argc > 4) {
		  std::string chain_id(argv[3]);
		  std::string res_no_str(argv[4]);
		  try {
		     int res_no = coot::util::string_to_int(res_no_str);
		     spec = coot::residue_spec_t(chain_id, res_no, "");
		  }
		  catch (const std::runtime_error &rte) {
		     std::cout << "ERROR::" << rte.what() << std::endl;
		  }
	       }
	       mmdb::Residue *r = coot::util::get_residue(spec, mol);

	       if (r) {

		  bool debugging = false;
		  if (argc == 6) {
		     std::string a5 = argv[5];
		     if (a5 == "debug")
			debugging = true;
		  }
		  coot::protein_geometry geom;
		  geom.set_verbose(0);
		  geom.init_standard();
		  geom.remove_planar_peptide_restraint();
		  std::pair<float, float> mv = coot::util::mean_and_variance(xmap);
		  mmdb::Residue *r_next = coot::util::next_residue(r);
		  mmdb::Residue *r_prev = coot::util::previous_residue(r);

		  bool try_2018 = true;
		  if (try_2018) {

		     float b_fact = 30.0;
		     int n_trials = 10000;

		     coot::minimol::fragment fC =
			coot::multi_build_terminal_residue_addition_forwards_2018(r, r_prev, "A", b_fact, n_trials,
										  geom, xmap, mv, debugging);

		     std::string file_name = "trace-frag-C-build.pdb";
		     coot::minimol::molecule mmm(fC);
		     mmm.write_file(file_name, 10);

		  } else {
		     if (false)
			coot::minimol::fragment fN =
			   coot::multi_build_N_terminal_ALA(r, r_next, "A", 20, 500, geom, xmap, mv, debugging);

		     coot::minimol::fragment fC =
			coot::multi_build_C_terminal_ALA(r, r_prev, "A", 20, 500, geom, xmap, mv, debugging);
		  }
	       }
	    }
	 } else {

	    // not testing rama multi-build


	    // test from a pdb file
	    std::string test_pdb_file_name = "test-trace-template.pdb";

	    if (argc > 2) {
	       std::string fn = argv[2];
	       if (coot::file_exists(fn)) 
		  test_pdb_file_name = fn;
	    }

	 
	    // test with a null moll or flood mol
	    if (coot::file_exists(test_pdb_file_name)) {

	       mmdb::Manager *mol = new mmdb::Manager;
	       mmdb::ERROR_CODE err = mol->ReadCoorFile(test_pdb_file_name.c_str());
	       if (! err) {

		  if (coot::file_exists("test-scales")) {
		     t.optimize_weights(mol);

		     std::cout << "------------- Done optimize_weights() " << std::endl;

		  } else { 
		     std::cout << "----------------------------------------------------\n";
		     std::cout << "----------------------------------------------------\n";
		     std::cout << "running with test mol: " << std::endl;
		     std::cout << "----------------------------------------------------\n";
		     std::cout << "----------------------------------------------------\n";
		     t.test_model(mol);
		  }
	       }
	    } else {

	       t.action();

	    }
	 }
      }
      
      // problem reading the map, perhaps?
      // 
      catch (const clipper::Message_fatal &mess) {
	 std::cout << "ERROR:: " << mess.text() << std::endl;
      }
   }
   return 0;
}
