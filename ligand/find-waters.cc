/* ligand/find-waters.cc
 * 
 * Copyright 2004, 2005 The University of York
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
#ifndef _MSC_VER
#include <unistd.h> // for getopt(3)
#endif

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

#include <stdlib.h>

// no dependency on coords files

// #include "coords/mmdb-extras.h"
// #include "coords/mmdb.h"
#include "ligand.hh"
#include "utils/coot-utils.hh"
#include "coot-utils/coot-map-utils.hh"
#include "coot-utils/coot-trim.hh"
#include "clipper/core/map_utils.h"
#include "clipper/ccp4/ccp4_map_io.h"
#include "clipper/ccp4/ccp4_mtz_io.h"
#include <iostream>

void show_usage(std::string pname) {
   std::cout << "Usage: " << pname
	     << " --pdbin pdb-in-filename" << " --hklin mtz-filename"
	     << " --f f_col_label"
	     << " --phi phi_col_label"
	     << " --pdbout waters-filename"
	     << " --sigma sigma-level"
	     << " --min-dist min-dist-to-protein"
	     << " --max-dist max-dist-to-protein"
	     << " --flood"
	     << " --flood-atom-radius"
	     << " --chop"
	     << "\n"
	     << "        --mapin ccp4-map-name can be used"
	     << " instead of --hklin --f --phi\n";
   std::cout << "        where pdbin is the protein (typically)\n"
	     << "        and pdbout is file for the waters.\n"
	     << "        The default sigma level is 2.0\n"
	     << "        Use --chop to remove waters below given sigma-level\n"
	     << "            In this case, pdbout is the modified input coordinates\n"
	     << "        Use --flood to fill everything with waters "
 	     << "(not just water peaks)\n"
	     << "        and --flood-atom-radius to adjust contact distance\n"
	     << "           (default 1.4A).\n";
} 

int
main(int argc, char **argv) {
   
   if (argc < 3) {
      show_usage(argv[0]); 
   } else { 

      std::string pdb_file_name;
      std::string  mtz_filename;
      std::string         f_col;
      std::string       phi_col;
      std::string    output_pdb;
      std::string     sigma_str;
      std::string map_file_name;
      std::string min_dist_str;
      std::string max_dist_str;
      short int do_flood_flag = 0;
      short int do_chop_flag = 0;
      float flood_atom_mask_radius = 1.4;
      float water_to_protein_max_dist = 99.9;
      float water_to_protein_min_dist = 1.5;

      const char *optstr = "i:h:f:p:o:s:e:c:m";
      struct option long_options[] = {
	 {"pdbin",  1, 0, 0},
	 {"hklin",  1, 0, 0},
	 {"f",      1, 0, 0},
	 {"phi",    1, 0, 0},
	 {"pdbout", 1, 0, 0},
	 {"sigma",  1, 0, 0},
	 {"mapin",  1, 0, 0},
	 {"min-dist",  1, 0, 0},
	 {"max-dist",  1, 0, 0},
	 {"flood-atom-radius",  1, 0, 0},
	 {"water-to-protein-max-dist",  1, 0, 0},
	 {"water-to-protein-min-dist",  1, 0, 0},
	 {"flood",  0, 0, 0},
	 {"chop",   0, 0, 0},
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
	       if (arg_str == "sigma") {
		  sigma_str = coot_optarg;
	       }
	       if (arg_str == "mapin") {
		  map_file_name = coot_optarg;
	       }
	       if (arg_str == "min-dist") {
		  min_dist_str = coot_optarg;
	       }
	       if (arg_str == "max-dist") {
		  max_dist_str = coot_optarg;
	       }
	       if (arg_str == "flood-atom-radius") {
		  try {
		     flood_atom_mask_radius =
			coot::util::string_to_float(coot_optarg);
		  }
		  catch (const std::exception &e) {
		     std::cout << "argument for --flood_atom_mask_radius"
			       << " is not a number: " << coot_optarg
			       << std::endl;
		     exit(1);
		  } 
	       }
               if (arg_str == "water-to-protein-max-dist") {
		  try {
		     water_to_protein_max_dist =
			coot::util::string_to_float(coot_optarg);
		  }
		  catch (const std::exception &e) {
		     std::cout << "argument for --water-to-protein-max-dist"
			       << " is not a number: " << coot_optarg
			       << std::endl;
		     exit(1);
		  } 
               }
               if (arg_str == "water-to-protein-min-dist") {
		  try {
		     water_to_protein_min_dist =
			coot::util::string_to_float(coot_optarg);
		  }
		  catch (const std::exception &e) {
		     std::cout << "argument for --water-to-protein-min-dist"
			       << " is not a number: " << coot_optarg
			       << std::endl;
		     exit(1);
		  } 
               }
	       
	    } else { 
	       std::string arg_str = long_options[option_index].name;
	       if (arg_str == "flood") {
		  do_flood_flag = 1;
	       } else {
		  if (arg_str == "chop") {
		     std::cout << "Removing waters";
		     do_chop_flag = 1;
		  } else { 
		     std::cout << "Malformed option: "
			       << long_options[option_index].name << std::endl;
		  }
	       }
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
	    
	 case 's':
	    sigma_str = coot_optarg;
	    break;

	 case 'e':
	    do_flood_flag = 1;
	    break;

	 case 'c':
	    do_chop_flag = 1;
	    break;
	    
	 default:
	    break;
	 }
      }

      short int do_it = 0;  // vs do flood
      short int do_it_with_map = 0;
      bool have_map = false; // meaning map was specified on the command line (was-given map)
      if (map_file_name.length() > 0)
	 have_map = true;

      if ( (pdb_file_name.length() == 0) && !do_flood_flag) {
	 std::cout << "Missing input PDB file\n";
	 exit(1);
      } else { 
	 if (output_pdb.length() == 0) { 
	    std::cout << "Missing output PDB file\n";
	    exit(1);
	 } else {
	    if (have_map) {
	       do_it_with_map = 1;
	       do_it = 1;
	       if (sigma_str.length() == 0) {
		  sigma_str = "2.0";
	       }
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
			if (sigma_str.length() == 0) {
			   sigma_str = "2.0";
			}
			do_it = 1;
		     }
		  }
	       }
	    }
	 }
      }
      
      if (! do_it) {

	 std::cout << "Bad command line" << std::endl;

      } else { 

	 float input_sigma_level = atof(sigma_str.c_str());
	 short int use_weights = 0;
	 short int is_diff_map = 0;
	 bool set_wpdl = false;
	 float wpdl_max = 3.2; // defaults in ligand() constructor
	 float wpdl_min = 2.4; //            ""

	 if (! min_dist_str.empty()) {
	    try {
	       wpdl_min = coot::util::string_to_float(min_dist_str);
	       set_wpdl = true;
	    }
	    catch (const std::exception &e) {
	       std::cout << "WARNING:: unable to convert: " << e.what() << std::endl;
	    }
	 }
	 if (! max_dist_str.empty()) {
	    try {
	       wpdl_max = coot::util::string_to_float(max_dist_str);
	       set_wpdl = true;
	    }
	    catch (const std::exception &e) {
	       std::cout << "WARNING:: unable to convert: " << e.what() << std::endl;
	    }
	 }

	 coot::ligand lig;
	 if (have_map) {
	    clipper::CCP4MAPfile file;
	    clipper::Xmap<float> xmap;
	    if (coot::file_exists(map_file_name)) {
	       file.open_read(map_file_name);
	       file.import_xmap(xmap);
	       file.close_read();
	       // this can take a long time for a cryo-em map
	       // is it needed for cryoem flooding?
	       clipper::Xmap<float> fine_xmap =
		  coot::util::reinterp_map_fine_gridding(xmap);
	       lig.import_map_from(fine_xmap);
	       // lig.output_map("fine-map.map");
	       xmap = fine_xmap;
	    } else {
	       std::cout << "ERROR:: Missing map " << map_file_name << "\n";
	       exit(1);
	    }
	 } else {
	    clipper::Xmap<float> xmap;
	    bool stat = coot::util::map_fill_from_mtz(&xmap,
						      mtz_filename,
						      f_col, phi_col, "", 
						      use_weights, is_diff_map);

	    if (! stat) {
	       std::cout << "ERROR: in filling map from mtz file: " << mtz_filename
			 << std::endl;
	       exit(1);
	    } else {
	       clipper::Xmap<float> fine_xmap =
		  coot::util::reinterp_map_fine_gridding(xmap);
	       lig.import_map_from(fine_xmap);
	    } 
	 }

	 if (! do_flood_flag) { 
	    if (pdb_file_name.length() == 0) {
	       std::cout << "confused input: no pdb input name and no --flood specified...\n";
	       exit(1);
	    } else { 
	       if (do_chop_flag == 0) { 

		  // mask_by_atoms() [which calls mask_map(flag)] here in
		  // find-waters currently behaves differently to
		  // c-interface-waters.cc's mask_map(mol, water_flag)
		  // resulting in move_atom_to_peak() moving wats to the
		  // wrong place
		  //
		  lig.set_map_atom_mask_radius(1.9);
		  lig.mask_by_atoms(pdb_file_name);
		  if (lig.masking_molecule_has_atoms()) { 
		     // lig.output_map("find-waters-masked.map");
		     // std::cout << "DEBUG:: in findwaters: using input_sigma_level: "
		     // 			    << input_sigma_level << std::endl;

		     if (set_wpdl)
			lig.set_water_to_protein_distance_limits(wpdl_max, wpdl_min); // max min
		     
		     lig.water_fit(input_sigma_level, 3); // e.g. 2.0 sigma for 3 cycles 
		     coot::minimol::molecule water_mol = lig.water_mol();
		     water_mol.write_file(output_pdb, 20.0);
		  } else {
		     std::cout << "No atoms found in masking molecule: " 
			       << pdb_file_name << std::endl;
		  }

	       } else {

		  // chop[py] waters!
		  clipper::Xmap<float> xmap;
		  coot::util::map_fill_from_mtz(&xmap, mtz_filename, f_col, phi_col, "", 
						use_weights, is_diff_map);

		  // we can't use atom_selection_container_t here
		  // libcoot-coords depends on libcoot-ligand - not the other way round.
		  //
		  // atom_selection_container_t asc = get_atom_selection(pdb_file_name, true, true);

		  mmdb::Manager *mol = new mmdb::Manager;
		  std::cout << "Reading coordinate file: " << pdb_file_name.c_str() << "\n";
		  mmdb::ERROR_CODE err = mol->ReadCoorFile(pdb_file_name.c_str());

		  if (err) {
		     std::cout << "There was an error reading " << pdb_file_name.c_str() << ". \n";
		     std::cout << "ERROR " << err << " READ: "
			       << mmdb::GetErrorDescription(err) << std::endl;
		  } else {
		     short int waters_only_flag = 1;
		     short int remove_or_zero_occ_flag = coot::util::TRIM_BY_MAP_DELETE;
		     clipper::Map_stats stats(xmap);
		     float map_level = stats.mean() + input_sigma_level * stats.std_dev();
		     int n_atoms = coot::util::trim_molecule_by_map(mol, xmap, map_level,
								    remove_or_zero_occ_flag,
								    waters_only_flag);
		     std::cout << "INFO:: " << n_atoms << " waters removed" << std::endl;
		     mol->WritePDBASCII(output_pdb.c_str());
		  }
	       }
	    }
	 } else {

	    std::cout << "===================== Flood mode ======================= "
		      << std::endl;
	    std::cout << "input_sigma_level:: " << input_sigma_level << std::endl;
	    // if a pdb file was defined, let's mask it
	    if (pdb_file_name.length() > 0) {
	       std::cout << "INFO:: masking map by coords in " << pdb_file_name
			 << std::endl;
	       // lig.set_map_atom_mask_radius();
	       lig.mask_by_atoms(pdb_file_name);
	    } 
	    lig.set_cluster_size_check_off();
	    lig.set_chemically_sensible_check_off();
	    lig.set_sphericity_test_off();

	    lig.set_map_atom_mask_radius(flood_atom_mask_radius);
	    // lig.set_water_to_protein_distance_limits(10.0, 1.5); // should not be 
                                                                    // used in lig.
            lig.set_water_to_protein_distance_limits(water_to_protein_max_dist,
                                                     water_to_protein_min_dist);
	    lig.flood2(input_sigma_level); // with atoms
	    coot::minimol::molecule water_mol = lig.water_mol();
	    water_mol.write_file(output_pdb, 30.0);
	    // lig.output_map("find-waters-masked-flooded.map");
	 }
      }
   }
   return 0;
}

