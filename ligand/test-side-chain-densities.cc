/* coot-utils/residue-and-atom-specs.hh
 * 
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

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <fstream>

#include <clipper/ccp4/ccp4_map_io.h>

#include "utils/coot-utils.hh"
#include "coot-utils/atom-selection-container.hh"
#include "coot-utils/xmap-stats.hh"
#include "side-chain-densities.hh"

void
make_useable_grid_points(int n_steps, float grid_box_radius,
			 int res_no, const std::string &file_name) {
   
   std::string pdb_file_name("test.pdb");
   atom_selection_container_t asc = get_atom_selection(pdb_file_name, true, false);
   coot::side_chain_densities scd(n_steps, grid_box_radius, "");
   coot::residue_spec_t spec_this("A", res_no, "");
   coot::residue_spec_t spec_next("A", res_no+1, "");
   mmdb::Residue *residue_this_p = coot::util::get_residue(spec_this, asc.mol);
   mmdb::Residue *residue_next_p = coot::util::get_residue(spec_next, asc.mol);
   if (residue_this_p && residue_next_p) {
      // residue_next_p is only needed when we are making exclusion points
      // (or vice versa) the acceptable grid points. In that case the N
      // is needed to exclude some grid points.
      scd.gen_useable_grid_points(residue_this_p, residue_next_p, n_steps, grid_box_radius, file_name);
   }
}

void
check_useable_grid_points(int n_steps, float grid_box_radius,
			  int res_no,
			  const std::string &useable_grid_points_file_name,
			  const std::string &useable_grid_points_mapped_to_residue) {

  std::string pdb_file_name("test.pdb");
   atom_selection_container_t asc = get_atom_selection(pdb_file_name, true, false);
   coot::residue_spec_t spec_this("A", res_no, "");
   mmdb::Residue *residue_p = coot::util::get_residue(spec_this, asc.mol);
   if (residue_p) {
      coot::side_chain_densities scd(n_steps, grid_box_radius, useable_grid_points_file_name);
      scd.check_useable_grid_points(residue_p, useable_grid_points_mapped_to_residue);
   }
}

void
check_stats(int n_steps, float grid_box_radius, const std::string &res_name,
	    const std::string &rot_name,
	    const std::string &pdb_file_name, const std::string &chain_id, int res_no,
	    const std::string &file_name) {

   // do the stats look like the side chain that they are supposed to be?

   // Show me a dotted grids where the size is proportional to the mean

   atom_selection_container_t asc = get_atom_selection(pdb_file_name, true, false);
   coot::residue_spec_t spec_this(chain_id, res_no, "");
   mmdb::Residue *residue_p = coot::util::get_residue(spec_this, asc.mol);
   if (residue_p) {
      std::string grid_points_file_name = file_name;
      coot::side_chain_densities scd(n_steps, grid_box_radius, grid_points_file_name);
      scd.check_stats(residue_p, res_name, rot_name);
   }
}

void
test_residue_vs_likelihoods(int n_steps, float grid_box_radius, int res_no, const std::string &file_name) {

   // Does this look like a TYR m-85 (say)?

   std::string pdb_file_name("test.pdb");
   std::string map_file_name("blurred-test.map");
   std::string grid_points_file_name = file_name;
   clipper::CCP4MAPfile file;
   try {
      file.open_read(map_file_name);
      clipper::Xmap<float> xmap;
      file.import_xmap(xmap);
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, true, false);
      if (asc.read_success) {
	 std::string id_1 = coot::util::name_sans_extension(pdb_file_name);
	 std::string id = coot::util::file_name_non_directory(id_1);

	 coot::side_chain_densities test_ll(n_steps, grid_box_radius, grid_points_file_name);

	 test_ll.set_data_dir("side-chain-data");
	 coot::residue_spec_t spec_this("A", res_no, "");
	 mmdb::Residue *residue_this_p = coot::util::get_residue(spec_this, asc.mol);
	 if (residue_this_p) {
	    // residue_next_p is only needed when we are making exclusion points
	    // (or vice versa) the acceptable grid points. In that case the N
	    // is needed to exclude some grid points.
	    test_ll.get_rotamer_likelihoods(residue_this_p, xmap);
	 } else {
	    std::cout << "residue " << spec_this << " not found in molecule from  "
		      << pdb_file_name << std::endl;
	 }
      }
   }
   catch (const clipper::Message_base &exc) {
      std::cout << "WARNING:: failed to open " << pdb_file_name << std::endl;
      std::cout << "WARNING:: failed to open " << map_file_name << std::endl;
   }
}

void find_probabilities_of_rotamers(int n_steps, float grid_box_radius,
				    const std::string &useable_grid_points_file_name,
				    const std::string &pdb_file_name,
				    const std::string &chain_id,
				    int resno_start,
				    int resno_end,
				    const std::string &map_file_name) {

   try {
      clipper::CCP4MAPfile file;
      file.open_read(map_file_name);
      clipper::Xmap<float> xmap;
      file.import_xmap(xmap);
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, true, false);
      if (asc.read_success) {
	 // "analysis" constructor
	 coot::side_chain_densities scd(n_steps, grid_box_radius, useable_grid_points_file_name);
	 scd.set_data_dir("side-chain-data");
	 scd.probability_of_each_rotamer_at_each_residue(asc.mol, chain_id, resno_start, resno_end, xmap);
      }
   }
   catch (const clipper::Message_base &exc) {
      std::cout << "WARNING:: failed to open " << pdb_file_name << std::endl;
      std::cout << "WARNING:: failed to open " << map_file_name << std::endl;
   }
}

void test_sequence(int n_steps, float grid_box_radius,
		   const std::string &useable_grid_points_file_name,
		   const std::string &pdb_file_name,
		   const std::string &chain_id,
		   int resno_start,
		   int resno_end,
		   const std::string &map_file_name,
		   const std::string &sequence) {

   std::cout << "test sequence " << sequence << std::endl;
   try {
      clipper::CCP4MAPfile file;
      file.open_read(map_file_name);
      clipper::Xmap<float> xmap;
      file.import_xmap(xmap);
      atom_selection_container_t asc = get_atom_selection(pdb_file_name, true, false);
      if (asc.read_success) {
	 // "analysis" constructor
	 coot::side_chain_densities scd(n_steps, grid_box_radius, useable_grid_points_file_name);
	 scd.set_data_dir("side-chain-data");
	 scd.test_sequence(asc.mol, chain_id, resno_start, resno_end, xmap, sequence);
      }
   }
   catch (const clipper::Message_base &exc) {
      std::cout << "WARNING:: failed to open " << pdb_file_name << std::endl;
      std::cout << "WARNING:: failed to open " << map_file_name << std::endl;
   }
}

// generate "stats.table" for every rotamer
void
combine(int n_steps) {

   // Set these from parsing the command line.
   // When optimizing, we need only re-run this combine stage - the grid data
   // generation stage need not be run multiple times.
   //
   unsigned int mn_unreliable_minimum_counts;
   unsigned int mn_unreliable_minimum_counts_for_low_variance;
   double mn_unreliable_minimum_variance;
   double mn_use_this_variance_for_unreliable;
   mn_unreliable_minimum_counts = 10;
   mn_unreliable_minimum_counts_for_low_variance = 20;
   mn_unreliable_minimum_variance = 0.1;
   mn_use_this_variance_for_unreliable = 4.0;

   std::cout << "combine() " << std::endl;
   std::string dir = "side-chain-data";

   // put the globbing in the side_chain_densities class
   //
   std::string glob_pattern = "*";
   std::vector<std::string> dirs = coot::util::glob_files(dir, glob_pattern);
   std::cout << "found " << dirs.size() << " files in " << dir << std::endl;

   for (std::size_t idir=0; idir<dirs.size(); idir++) {
      const std::string &dir = dirs[idir];
      // std::cout << "dir " << dir << std::endl;

      std::vector<std::string> rot_dirs = coot::util::glob_files(dir, glob_pattern);
      std::cout << "found " << rot_dirs.size() << " files in " << dir << std::endl;
      for (std::size_t jdir=0; jdir<rot_dirs.size(); jdir++) {
	 const std::string &dir = rot_dirs[jdir];
	 std::cout << "combine(): rot_dir: " << dir << std::endl;

	 coot::side_chain_densities::combine_directory(dir, n_steps,
						       mn_unreliable_minimum_counts,
						       mn_unreliable_minimum_counts_for_low_variance,
						       mn_unreliable_minimum_variance,
						       mn_use_this_variance_for_unreliable);
      }
   }
}

void
check_map_mean(const std::string &map_file_name) {

   clipper::CCP4MAPfile file;
   try {
      clipper::Xmap<float> xmap;
      file.open_read(map_file_name);
      file.import_xmap(xmap);
      mean_and_variance<float> mv_0 = map_density_distribution(xmap, 400, false, false);
      mean_and_variance<float> mv_1 = map_density_distribution(xmap, 400, false, true);
      float n_sds_0 = mv_0.mean/sqrt(mv_0.variance);
      float n_sds_1 = mv_1.mean/sqrt(mv_1.variance);
      std::cout << map_file_name << " "
                << std::fixed << std::right << std::setprecision(12)
		<< mv_0.mean << " " << sqrt(mv_0.variance) << " " << n_sds_0 << "    "
		<< mv_1.mean << " " << sqrt(mv_1.variance) << " " << n_sds_1
		<< std::endl;
   }
   catch (const clipper::Message_base &exc) {
      std::cout << "WARNING:: failed to open " << map_file_name << std::endl;
   }

}

int main(int argc, char **argv) {

   bool done = false;

   int n_steps = 5;
   float grid_box_radius = 5.0; // half the width/length of the box (not diagonal)

   if (argc == 2) {
      // generate the stats from the sampled maps.
      std::string a1(argv[1]);
      if (a1 == "combine") {
	 combine(n_steps);
	 done = true;
      }
   }

   if (argc == 4) {
      std::string a1(argv[1]);
      // make a file that includes grid_idx x y z
      if (a1 == "generate-useable-grid-points") {
	 int res_no = coot::util::string_to_int(argv[2]);
	 std::string grid_points_file_name(argv[3]);
	 make_useable_grid_points(n_steps, grid_box_radius, res_no, grid_points_file_name);
	 done = true;
      }
   }

   if (argc == 5) {
      std::string a1(argv[1]);
      if (a1 == "check-generated-useable-grid-points") {
	 int res_no = coot::util::string_to_int(argv[2]);
	 std::string grid_points_file_name(argv[3]);
	 std::string useable_grid_points_mapped_to_residue(argv[4]);
	 check_useable_grid_points(n_steps, grid_box_radius, res_no,
				   grid_points_file_name, useable_grid_points_mapped_to_residue);
	 done = true;
      }
   }

   if (argc == 4) {
      std::string a1(argv[1]);
      if (a1 == "test-residue") {
	 // what are the probabilities that this residue is any of the rotamers
	 // using test.pdb and blurred-test.map
	 int res_no = coot::util::string_to_int(argv[2]);
	 std::string grid_points_file_name(argv[3]);
	 test_residue_vs_likelihoods(n_steps, grid_box_radius, res_no, grid_points_file_name);
	 done = true;
      }
   }

   if (argc == 8) {
      // Make a 3D dot plot - representing the stats for the given rotamer
      // point (size) represent values in stats.table
      // different rotamers should look different
      std::string a1(argv[1]);
      if (a1 == "check-stats") {
	 std::string res_name(argv[2]);
	 std::string rot_name(argv[3]);
	 std::string pdb_file_name(argv[4]);
	 std::string chain_id(argv[5]);
	 std::string res_no_str(argv[6]);
	 std::string grid_points_file_name(argv[7]);
	 int res_no = coot::util::string_to_int(res_no_str); // protect with try/catch
	 // the residue to which the grid is matched is set in the function
	 check_stats(n_steps, grid_box_radius, res_name, rot_name,
		     pdb_file_name, chain_id, res_no,
		     grid_points_file_name);
	 done = true;
      }
   }

   if (argc == 8) {
      std::string a1(argv[1]);
      if (a1 == "find-probabilities-of-rotamers") {
	 try {
	    std::string pdb_file_name(argv[2]); // poly-ALA model
	    std::string chain_id(argv[3]);
	    int resno_start = coot::util::string_to_int(argv[4]);
	    int resno_end   = coot::util::string_to_int(argv[5]);
	    std::string map_file_name(argv[6]);
	    std::string useable_grid_points_file_name(argv[7]);
	    find_probabilities_of_rotamers(n_steps, grid_box_radius, useable_grid_points_file_name,
					   pdb_file_name, chain_id, resno_start, resno_end,
					   map_file_name);
	 }
	 catch (const std::runtime_error &rte) {
	    std::cout << "ERROR:: " << rte.what() << std::endl;
	 }
	 done = true;
      }
   }

   if (argc == 9) {
      std::string a1(argv[1]);
      if (a1 == "test-sequence") {
	 try {
	    std::string pdb_file_name(argv[2]); // poly-ALA model
	    std::string chain_id(argv[3]);
	    int resno_start = coot::util::string_to_int(argv[4]);
	    int resno_end   = coot::util::string_to_int(argv[5]);
	    std::string map_file_name(argv[6]);
	    std::string sequence(argv[7]);
	    std::string useable_grid_points_file_name(argv[8]);
	    test_sequence(n_steps, grid_box_radius, useable_grid_points_file_name,
			  pdb_file_name, chain_id, resno_start, resno_end,
			  map_file_name, sequence);
	 }
	 catch (const std::runtime_error &rte) {
	    std::cout << "ERROR:: " << rte.what() << std::endl;
	 }
	 done = true;
      }
   }

   if (argc == 3) {
      std::string a1(argv[1]);
      if (a1 == "test-grid-interconversion") {
	 std::string useable_grid_points_file_name(argv[2]);
	 coot::side_chain_densities scd(n_steps, grid_box_radius, useable_grid_points_file_name);
	 scd.set_data_dir("side-chain-data");
	 bool success = scd.test_grid_point_to_coords_interconversion();
	 if (success)
	    std::cout << "Correct" << std::endl;
	 done = true;
      }
   }

   if (argc == 3) {
      std::string a1(argv[1]);
      if (a1 == "calc-mean") {
	 std::string map_file_name(argv[2]);
	 std::cout << "### check means for map " << map_file_name << std::endl;
	 check_map_mean(map_file_name);
	 done = true;
      }
   }

   if (! done) {
      // The engine of this program.
      //
      // Make the grid samples for each residue of each structure and map
      //
      if (argc == 5) {
	 std::string a1(argv[1]);
	 if (a1 == "sample") {
	    std::string pdb_file_name(argv[2]);
	    std::string map_file_name(argv[3]);
	    std::string useable_grid_points_file_name(argv[4]);
	    clipper::CCP4MAPfile file;
	    try {
	       file.open_read(map_file_name);
	       clipper::Xmap<float> xmap;
	       file.import_xmap(xmap);
	       atom_selection_container_t asc = get_atom_selection(pdb_file_name, true, false);
	       if (asc.read_success) {
		  std::string id_1 = coot::util::name_sans_extension(pdb_file_name);
		  std::string id = coot::util::file_name_non_directory(id_1);
		  coot::side_chain_densities(id, asc.mol, n_steps, grid_box_radius, xmap,
					     useable_grid_points_file_name);
	       }
	    }
	    catch (const clipper::Message_base &exc) {
	       std::cout << "WARNING:: failed to open " << pdb_file_name << std::endl;
	       std::cout << "WARNING:: failed to open " << map_file_name << std::endl;
	    }
	    done = true;
	 }
      }
   }

   if (! done) std::cout << "Nothing done\n";
   return 0;
} 
