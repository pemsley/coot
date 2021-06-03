/* ligand/find-ligands.cc
 * 
 * Copyright 2004, 2005, 2006, 2007 The University of York
 * Author: Paul Emsley
 * Copyright 2007 The University of Oxford
 * Copyright 2012, 2014 by Medical Research Council
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

// Portability (getopt) gubbins
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

#include <iostream>
#include <string>
#include <stdexcept>

#include "utils/coot-utils.hh"
#include "coot-utils/coot-coord-utils.hh"
#include "wligand.hh"

int
main(int argc, char **argv) {

   if (argc < 6) { 
      std::cout << "Usage: " << argv[0] << "\n     "
		<< " --pdbin pdb-in-filename" << " --hklin mtz-filename" << "\n     "
		<< " --f f_col_label" << " --phi phi_col_label" << "\n     "
		<< " --clusters nclust"   << "\n     "
		<< " --sigma sigma-level" << "\n     "
		<< " --absolute level"    << "\n     "
		<< " --fit-fraction frac" << "\n     "
		<< " --flexible"          << "\n     "
		<< " --samples nsamples"  << "\n     "
		<< " --sampling-rate map-sampling-rate"  << "\n     "
		<< " --dictionary cif-dictionary-name" << "\n     "
		<< " --script script-file-name\n"
		<< "   ligand-pdb-file-name(s)\n";
      std::cout << "     where pdbin is the protein (typically)\n"
		<< "           nclust is the number of clusters to fit [default 10]\n"
		<< "           sigma is the search level of the map (default 2.0)\n"
		<< "              or specify --absolute and pass level in e/A^3\n"
		<< "           --flexible means use torsional conformation ligand search\n"
		<< "           --dictionary file containing the CIF ligand dictionary description\n"
		<< "           nsamples is the number of flexible conformation samples [default 30]\n"
		<< "           frac is the minimum fraction of atoms in density allowed after fit [default 0.75]\n"
		<< "           script-file-name is a file name of helper script suitable for use in Coot\n"
		<< "               (default: coot-ligands.scm)."
		<< std::endl;

   } else { 

      int n_used_args = 1;
      std::string pdb_file_name;
      std::string  mtz_filename;
      std::string         f_col;
      std::string       phi_col;
      std::string     sigma_str;
      std::string     n_cluster_string;
      std::vector<std::string> lig_files;
      short int use_wiggly_ligand = 0;
      int wiggly_ligand_n_samples = 30;
      std::string cif_file_name = "";
      std::string fit_frac_str = "";
      std::string coot_ligands_script_file_name = "coot-ligands.scm";
      float fit_frac = -1.0; // tested for positivity
      bool set_absolute = 0;
      float absolute_level = 0.0;
      std::string absolute_string = "";
      bool blobs_mode = false;
      float map_sampling_factor = 1.5;
      std::string map_sampling_factor_str;

      // These hold the coordinates of a particular position.  The
      // bool is whether they were set in the input or not.  Only of
      // all 3 are set will findligand search in a particular position
      // (pretty obviously).
      // 
      std::pair<bool,float> pos_x(0, 0.0), pos_y(0, 0.0), pos_z(0, 0.0);

      const char *optstr = "i:h:f:p:s:c:w:n:d";
      struct option long_options[] = {
	 {"pdbin",         1, 0, 0},
	 {"hklin",         1, 0, 0},
	 {"f",             1, 0, 0},
	 {"phi",           1, 0, 0},
	 {"sigma",         1, 0, 0},
	 {"absolute",      1, 0, 0},
	 {"clusters",      1, 0, 0},
	 {"samples",       1, 0, 0},
	 {"dictionary",    1, 0, 0},
	 {"fit-fraction",  1, 0, 0},
	 {"flexible",      0, 0, 0},
	 {"script",        0, 0, 0},
	 {"blobs",         0, 0, 0},
	 {"sampling-factor", 1, 0, 0},
	 {"pos-x",         1, 0, 0},
	 {"pos-y",         1, 0, 0},
	 {"pos-z",         1, 0, 0},
	 {0, 0, 0, 0}
      };

      int ch;
      int option_index = 0;
      while ( -1 != 
	      (ch = coot_getopt_long(argc, argv, optstr, long_options, &option_index))) {

	 switch(ch) {
	    
	 case 0:
	    if (coot_optarg) { 

// 	       std::cout << "DEBUG:: " << option_index << " " << strlen(coot_optarg) << std::endl;
// 	       std::cout << " " << coot_optarg << std::endl;
// 	       std::cout << "   ch:: " << ch << std::endl;
	 
	       std::string arg_str = long_options[option_index].name;

	       // std::cout << " considering arg_str :" << arg_str << ":\n";

	       if (arg_str == "pdbin") { 
		  pdb_file_name = coot_optarg;
		  n_used_args += 2;
	       } 
	       if (arg_str == "pdb") { 
		  pdb_file_name = coot_optarg;
		  n_used_args += 2;
	       } 
	       if (arg_str == "hklin") { 
		  mtz_filename = coot_optarg;
		  n_used_args += 2;
	       } 
	       if (arg_str == "f") { 
		  f_col = coot_optarg;
		  n_used_args += 2;
	       } 
	       if (arg_str == "phi") {
		  phi_col = coot_optarg;
		  n_used_args += 2;
	       } 
	       if (arg_str == "sigma") {
		  sigma_str = coot_optarg;
		  n_used_args += 2;
	       } 
	       if (arg_str == "absolute") {
		  absolute_string = coot_optarg;
		  set_absolute = 1;
		  n_used_args += 2;
	       } 
	       if (arg_str == "clusters") {
		  n_cluster_string = coot_optarg;
		  n_used_args += 2;
	       }

	       if (arg_str == "samples") { 
		  wiggly_ligand_n_samples = atoi(coot_optarg);
		  n_used_args += 2;
	       }

	       if (arg_str == "dictionary") { 
		  cif_file_name = coot_optarg;
		  n_used_args += 2;
	       }

	       if (arg_str == "sampling-factor") { 
		  map_sampling_factor_str = coot_optarg;
		  n_used_args += 2;
	       }

	       if (arg_str == "fit-fraction") { 
		  fit_frac_str = coot_optarg;
		  n_used_args += 2;
	       }
	       
	       if (arg_str == "script") { 
		  coot_ligands_script_file_name = coot_optarg;
		  n_used_args += 2;
	       }

	       if (arg_str == "pos-x" || arg_str == "pos-y" || arg_str == "pos-z") {
		  try {
		     float v = coot::util::string_to_float(coot_optarg);
		     if (arg_str == "pos-x")
			pos_x = std::pair<bool, float> (1, v);
		     if (arg_str == "pos-y")
			pos_y = std::pair<bool, float> (1, v);
		     if (arg_str == "pos-z")
			pos_z = std::pair<bool, float> (1, v);
		  }
		  catch (const std::runtime_error &rte) {
		     std::cout << "WARNING:: Failed to convert " << arg_str << " to a number"
			       << std::endl;
		  } 
		  n_used_args += 2;
	       } 
	       
	    } else {

	       // options without arguments:
	       
	       // long argument without parameter:
	       std::string arg_str(long_options[option_index].name);
	       
	       if (arg_str == "flexible") {
		  use_wiggly_ligand = 1;
		  n_used_args += 1;
	       }

	       if (arg_str == "blobs") { 
		 blobs_mode = true;
		 std::cout << "set blobs mode " << std::endl;
		 n_used_args++;
	       } 

	    }
	    break;

	 case 'i':
	    pdb_file_name = coot_optarg;
	    n_used_args += 2;
	    break;
	    
	 case 'h':
	    mtz_filename = coot_optarg;
	    n_used_args += 2;
	    break;
	    
	 case 'f':
	    f_col = coot_optarg;
	    n_used_args += 2;
	    break;
	    
	 case 'p':
	    phi_col = coot_optarg;
	    n_used_args += 2;
	    break;
	    
	 case 's':
	    sigma_str = coot_optarg;
	    n_used_args += 2;
	    break;
	    
	 case 'c':
	    n_cluster_string = coot_optarg;
	    n_used_args += 2;
	    break;
	    
	 default:
	    std::cout << "default coot_optarg: " << coot_optarg << std::endl;
	    break;
	 }
      }

      bool do_it = false;
      if (pdb_file_name.length() == 0) { 
	 std::cout << "Missing input PDB file\n";
	 exit(1);
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
		  if (n_cluster_string.length() == 0) { 
		     n_cluster_string = "10";
		  }

 		  // std::cout << "------------ argc: " << argc << " n_used_args: " 
		  // << n_used_args << std::endl;

		  // 20140520 What was I thinking?  The last argument
		  // is the one where the ligand pdb file is.  Now we
		  // check all the arguments.
		  // 
		  // for (int i=n_used_args; i<(argc-1); i++)
		  
		  for (int i=n_used_args; i<argc; i++) {
		     std::cout << "----- pushing back ligand file name \"" << argv[i] << "\"" 
			       << std::endl;
		     lig_files.push_back(argv[i]);
		  }
		  do_it = true;
	       }
	    }
	 }
      }
      
      if (do_it) { 

	if ((lig_files.size() == 0) && ! blobs_mode) { 
	     std::cout << "No ligand pdb files specified\n";
	     exit(1);

	 } else {

	    if (fit_frac_str != "") { // it was set
	       fit_frac = atof(fit_frac_str.c_str());
	    } else {
	       fit_frac = 0.75;
	    } 
	    std::cout << "INFO:: Using acceptable fit fraction of: " << fit_frac << std::endl;
	    float input_sigma_level = 2.0;
	    try {
	       input_sigma_level = coot::util::string_to_float(sigma_str);
	    }
	    catch (const std::runtime_error &rte) {
	       std::cout << rte.what() << " - using default 2.0 sigma";
	    }

	    short int use_weights = 0;
	    short int is_diff_map = 0; 
	    coot::wligand wlig;
	    // wlig.set_verbose_reporting();

	    // sampling rate
	    if (! map_sampling_factor_str.empty()) {
	       try {
		  map_sampling_factor = coot::util::string_to_float(map_sampling_factor_str);
	       }
	       catch(const std::exception &e) {
		  std::cout << "failed to convert " << map_sampling_factor_str << std::endl;
	       }
	    }
	       

	    if (! use_wiggly_ligand) {

	       // rigid ligands path

	       short int map_stat = wlig.map_fill_from_mtz(mtz_filename, f_col, phi_col, "", 
							   use_weights, is_diff_map,
							   map_sampling_factor);

	       // Now we have the map, we can the convert absolute level (if it was set by
	       // the user to absolute) which is what the interface to the ligand class
	       // want to be passed).
	       // 
	       if (set_absolute) {
		  try { 
		     float ab = coot::util::string_to_float(absolute_string);
		     clipper::Map_stats stats = wlig.map_statistics();
		     input_sigma_level = ab/stats.std_dev();
		     std::cout << "Masking with level " << absolute_string << "("
			       << input_sigma_level << " sigma)" << std::endl;
		  }
		  catch (const std::runtime_error &rte) {
		     std::cout << rte.what() << std::endl;
		  }
	       } 

	       if (map_stat == 0) { 
		     std::cout << "Map making failure." << std::endl;
	       } else {
		  wlig.mask_by_atoms(pdb_file_name);
		  if (wlig.masking_molecule_has_atoms()) {
		     wlig.set_acceptable_fit_fraction(fit_frac); 
		     // lig.output_map("find-ligands-masked.map"); // debugging

		     if (pos_x.first && pos_y.first && pos_z.first) {

			clipper::Coord_orth pt(pos_x.second, pos_y.second, pos_z.second);
			wlig.cluster_from_point(pt, input_sigma_level);
			wlig.fit_ligands_to_clusters(1); // just this cluster.
			
		     } else { 
			
			if (blobs_mode) { 
			   int n_cycles = 1;
			   wlig.water_fit(input_sigma_level, n_cycles);
			   std::vector<std::pair<clipper::Coord_orth, double> > big_blobs = wlig.big_blobs();
			   unsigned int n_big_blobs = big_blobs.size();
			   if (n_big_blobs) {
			      std::cout << "=============== start blob-table ==========\n";
			      for (unsigned int i=0; i<n_big_blobs; i++) {
				 std::cout << "  blob " << i << " " << big_blobs[i].first.format()
					   << " sum: " << big_blobs[i].second << std::endl;
			      } 
			      std::cout << "=============== end blob-table ==========\n";
			   } 
			} else { 
			   wlig.find_clusters(input_sigma_level);
			}
		     }
		     // install ligands:
		     for (unsigned int ilig=0; ilig<lig_files.size(); ilig++)
			wlig.install_ligand(lig_files[ilig]);
		     wlig.fit_ligands_to_clusters(10); 
		  } else {
		     std::cout << "No atoms found in masking molecule: " 
			       << pdb_file_name << std::endl;
		  }
	       }

	    } else {

	       // wiggly ligands path

	       coot::protein_geometry geom;
	       // wlig.set_verbose_reporting();
	       
	       // this might be a pain if the flexible ligand is in the standard
	       // refmac dictionary...
	       if (cif_file_name.length() == 0) {
		  std::cout << "No cif dictionary file given\n";
	       } else { 
		  // geom.init_standard();
		  coot::read_refmac_mon_lib_info_t rmit = geom.init_refmac_mon_lib(cif_file_name, 10);
		  if (rmit.n_bonds == 0) {
		     std::cout << "Critical cif dictionary reading failure." << std::endl;
		  } else { 
		     short int map_stat = wlig.map_fill_from_mtz(mtz_filename, f_col, phi_col, "", 
								 use_weights, is_diff_map,
								 map_sampling_factor);

		     // Now we have the map, we can the convert absolute level (if it was set by
		     // the user to absolute) which is what the interface to the ligand class
		     // want to be passed).
		     // 
		     if (set_absolute) {
			try { 
			   float ab = coot::util::string_to_float(absolute_string);
			   clipper::Map_stats stats = wlig.map_statistics();
			   input_sigma_level = ab/stats.std_dev();
			}
			catch (const std::runtime_error &rte) {
			   std::cout << rte.what() << std::endl;
			}
		     } 

		     
		     if (map_stat == 0) {
			std::cout << "Map making failure." << std::endl;
		     } else { 
			
			wlig.mask_by_atoms(pdb_file_name);
			if (wlig.masking_molecule_has_atoms()) { 
			   wlig.output_map("find-ligands-masked.map");
			   wlig.set_acceptable_fit_fraction(fit_frac);

			   
			   if (pos_x.first && pos_y.first && pos_z.first) {
			      
			      clipper::Coord_orth pt(pos_x.second, pos_y.second, pos_z.second);
			      wlig.cluster_from_point(pt, input_sigma_level);
			      wlig.fit_ligands_to_clusters(1); // just this cluster.
			      
			   } else { 
			      wlig.find_clusters(input_sigma_level);

			   }

			   // install wiggly ligands...
			   // 
			   // wiggly ligands currently have to be minimols
			   for (unsigned int ilig=0; ilig<lig_files.size(); ilig++) { 
			      coot::minimol::molecule mmol;
			      mmol.read_file(lig_files[ilig]);

			      bool optim_geom = 1;
			      bool fill_vec = 0;
			      try {
				 int imol_ligand = 0;
				 std::vector<coot::installed_wiggly_ligand_info_t>
				    wiggled_ligands =
				    wlig.install_simple_wiggly_ligands(&geom, mmol,
								       imol_ligand,
								       wiggly_ligand_n_samples,
								       optim_geom,
								       fill_vec);
			      }
			      catch (const std::runtime_error &mess) {
				 std::cout << "Failed to install flexible ligands\n" << mess.what()
					   << std::endl;
			      } 
			   }
			   wlig.fit_ligands_to_clusters(10); 
			}
		     }
		  }
	       }
	    }

	    // what are the results?
	    
	    // now add in the solution ligands:
	    int n_clusters = wlig.n_clusters_final();

	    int n_new_ligand = 0;
	    coot::minimol::molecule m;
	    for (int iclust=0; iclust<n_clusters; iclust++) {

	       // frac_lim is the fraction of the score of the best solutions
	       // that we should consider as solutions. 0.5 is generous, I
	       // think.
	       float frac_lim = 0.7;
	       float correl_frac_lim = 0.9;
	       
	       // nino-mode
	       unsigned int nlc = 1;
	       wlig.score_and_resort_using_correlation(iclust, nlc);

	       correl_frac_lim = 0.975;

	       if (nlc > 12) nlc = 12; // arbitrary limit of max 12 solutions per cluster
	       float tolerance = 20.0;
	       // limit_solutions should be run only after a post-correlation sort.
	       //
	       wlig.limit_solutions(iclust, correl_frac_lim, nlc, tolerance, true);
			   
	       for (unsigned int isol=0; isol<nlc; isol++) { 

		  m = wlig.get_solution(isol, iclust);
		  if (! m.is_empty()) {
		     float bf = 30;
		     mmdb::Manager *ligand_mol = m.pcmmdbmanager();
		     coot::hetify_residues_as_needed(ligand_mol);

		     std::string file_name = "fitted-ligand-" +
			coot::util::int_to_string(iclust) + std::string("-") + 
			coot::util::int_to_string(isol) + std::string(".pdb");
		     ligand_mol->WritePDBASCII(file_name.c_str());
		  }
	       }
	    }
	 }
      }
   }
   return 0; 
}
