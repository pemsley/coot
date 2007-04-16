/* ideal/with-map.cc
 * 
 * Copyright 2002, 2003 The University of York
 * Author: Paul Emsley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
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

#include <sys/types.h> // for stating
#include <sys/stat.h>
#include <unistd.h>

#include <iostream>
#include <string.h>
#include <math.h>

#include "clipper/core/xmap.h"
#include "clipper/ccp4/ccp4_mtz_io.h"
#include "clipper/core/hkl_compute.h"
#include "clipper/core/map_utils.h" // Map_stats

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

clipper::Xmap<float>
map_from_mtz(std::string mtz_file_name,
	     std::string f_col,
	     std::string phi_col,
	     std::string weight_col,
	     int use_weights,
	     int is_diff_map);


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
      std::cout << "usage: " << argv[0] << " pdb_filename " << std::endl;

   } else {

      geom.init_standard();

      string pdb_file_name(argv[1]);

      // if pdb_file_name does not exist -> crash?
      atom_selection_container_t asc = get_atom_selection(pdb_file_name); 
      //coot::restraints_container_t restraints(asc);

      // So, we provide easy(?) access to the atoms of next and
      // previous residues (to those in the atom selection
      // moving_residue_atoms).  It is also possible to select "fixed"
      // atoms in the graphics (so that they don't move).  Let's
      // provide a vector of indices in the the moving_residue_atoms
      // array to define those (lovely mixture of styles - heh).
      //
      std::vector<CAtom *> fixed_atoms;

      // This interface has been withdrawn because we need the whole
      // molecule (acutally, a pointer to it) to do some atom selection.
      // 
//       coot::restraints_container_t restraints(asc.atom_selection, // moving_residue_atoms,
// 					      asc.n_selected_atoms,
// 					      previous_residue,
// 					      next_residue,
// 					      fixed_atoms);
      
      // int istart_res = 72;
      // int iend_res   = 72;  // ropey.pdb 

      int istart_res = 89;
      int iend_res   = 89;  // ropey.pdb
      char *chain_id   = asc.atom_selection[0]->residue->GetChainID();

      clipper::Xmap<float> xmap = map_from_mtz(mtz_filename,
					       f_col,
					       phi_col,
					       "", 0, 0);
					       
      std::string altloc("");

      coot::restraints_container_t restraints(istart_res,
					      iend_res,
					      1, 1, 0,
					      altloc,
					      chain_id,
					      asc.mol,
					      fixed_atoms,
					      xmap, 1000.0);
      
      // coot::restraint_usage_Flags flags = coot::NO_GEOMETRY_RESTRAINTS;
      //coot::restraint_usage_Flags flags = coot::BONDS;
      // coot::restraint_usage_Flags flags = coot::BONDS_AND_ANGLES;
      coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_AND_PLANES;
      // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_TORSIONS_AND_PLANES; 
      // coot::restraint_usage_Flags flags = coot::BONDS_ANGLES_PLANES_AND_NON_BONDED;
      flags = coot::NON_BONDED;
      
      coot::pseudo_restraint_bond_type pseudos = coot::NO_PSEUDO_BONDS;
      restraints.make_restraints(geom, flags, 1, 0, pseudos);


      restraints.minimize(flags);
      restraints.write_new_atoms("new.pdb");
   }

#endif // HAVE_GSL
   return 0; 
} 

clipper::Xmap<float>
map_from_mtz(std::string mtz_file_name,
		  std::string f_col,
		  std::string phi_col,
		  std::string weight_col,
		  int use_weights,
		  int is_diff_map) {

   
   clipper::HKL_info myhkl; 
   clipper::MTZdataset myset; 
   clipper::MTZcrystal myxtl; 

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
   clipper::Xmap<float> xmap;
   

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

  cout << "finding ASU unique map points..." << endl;
  xmap.init( myhkl.spacegroup(), myhkl.cell(),
	     clipper::Grid_sampling( myhkl.spacegroup(),
				     myhkl.cell(),
				     myhkl.resolution()) );
  cout << "Grid..." << xmap.grid_sampling().format() << "\n";


  cout << "doing fft..." << endl;
  xmap.fft_from( fphidata );                  // generate map
  cout << "done fft..." << endl;

  clipper::Map_stats stats(xmap);

  std::cout << "map mean " << stats.mean() << ", std dev: " << stats.std_dev() << std::endl;

  return xmap;
    
}
