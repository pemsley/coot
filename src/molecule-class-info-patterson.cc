/* src/molecule-class-info-maps.cc
 * 
 * Copyright 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
 * Copyright 2007 by Paul Emsley
 * Copyright 2007, 2008, 2009, 2010 by The University of Oxford
 * Copyright 2016 by Medical Research Council
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

#ifdef USE_PYTHON
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif


#include "compat/coot-sysdep.h"

// Having to set up the include files like this so that
// molecule-class-info.h can be parsed, is silly.

// For stat, mkdir:
#include <sys/types.h>
#include <sys/stat.h>

#include <string>
#include <mmdb2/mmdb_manager.h>
#include "coords/mmdb-extras.h"
#include "coords/Cartesian.h"
#include "coords/mmdb-crystal.h"
#include "molecule-class-info.h"
#include "coot-utils/coot-coord-utils.hh"
#include "density-contour/CIsoSurface.h"

#include "clipper/ccp4/ccp4_mtz_io.h"
// #include "clipper/cns/cns_hkl_io.h"
// #include "clipper/cns/cns_map_io.h"
// #include "clipper/core/hkl_compute.h"
#include "clipper/core/map_utils.h" // Map_stats
// #include "clipper/core/resol_basisfn.h"
// #include "clipper/core/resol_targetfn.h"
// #include "clipper/mmdb/clipper_mmdb.h"
// #include "clipper/clipper-phs.h"
// #include "clipper/contrib/sfcalc_obs.h"
// #include "clipper/contrib/sfscale.h"
// #include "clipper/contrib/sfweight.h"

// #include "clipper/clipper-cif.h"
// #include "clipper/contrib/sfcalc.h"

// #include "graphics-info.h"

#include "xmap-utils.h"
#include "coot-utils/xmap-stats.hh"



bool
molecule_class_info_t::make_patterson(std::string mtz_file_name,
				      std::string f_col,
				      std::string sigf_col,
				      float map_sampling_rate) {

   bool ret_val = 0;
   if (coot::file_exists(mtz_file_name)) {

      try {
	 // KDC Code (examples/cpatterson.cpp)
	 // 
	 clipper::CCP4MTZfile mtzin;
	 mtzin.open_read(mtz_file_name);
	 clipper::Spacegroup spgr = mtzin.spacegroup();
	 clipper::Cell       cell = mtzin.cell();
	 clipper::Resolution reso = mtzin.resolution();
	 clipper::HKL_info hkls, hklp;
	 typedef clipper::HKL_data_base::HKL_reference_index HRI;

	 clipper::HKL_data<clipper::data32::F_sigF> fsig(hkls);

	 // If use weights, use both strings, else just use the first
	 int use_weights = 0;        // dummy value
	 std::string weight_col = ""; // dummy value
	 std::pair<std::string, std::string> p =
	    make_import_datanames(f_col, sigf_col, weight_col, use_weights);
	 // std::cout << "======= import_hkl_data using " << p.first << std::endl;
	 mtzin.import_hkl_data(fsig, p.first);
	 mtzin.close_read();
      
	 // get Patterson spacegroup
	 clipper::Spacegroup pspgr(clipper::Spgr_descr(spgr.generator_ops().patterson_ops()));
	 hklp.init( pspgr, cell, reso, true );
      
	 // make patterson coeffs
	 int count = 0;
	 clipper::HKL_data<clipper::data32::F_phi> fphi( hklp );
	 for ( HRI ih = fphi.first(); !ih.last(); ih.next() ) {
	    clipper::data32::F_sigF f = fsig[ih.hkl()];
	    fphi[ih].phi() = 0.0 ;
	    if ( !f.missing() ) {
	       fphi[ih].f() = f.f()*f.f();
	       count++;
	    }
	 }
	 std::cout << "DEBUG:: read " << count << " reflections" << std::endl;
	 if (count) {

	    // This is a problem for Kevin
	    // original_fphis_filled = true;
	    // original_fphis.init(fphi.spacegroup(),fphi.cell(),fphi.hkl_sampling());
	    // original_fphis = fphi;
	 }

	 clipper::Grid_sampling grid;
	 grid.init( pspgr, cell, reso );
      
	 // make xmap
	 clipper::Xmap<float> xmap(pspgr, cell, grid);
	 xmap.fft_from(fphi);
	 is_patterson = 1; // needed for contour level protection
			   // (Pattersons are off-scale, but we should
			   // be able to change contour).

	 std::string map_name = "Patterson ";
	 map_name += mtz_file_name;
	 map_name += " ";
	 map_name += f_col;
	 new_map(xmap, map_name);
	 update_map_in_display_control_widget();
	 mean_and_variance<float> mv = map_density_distribution(xmap, 40, false);
	 map_mean_  = mv.mean; 
	 map_sigma_ = sqrt(mv.variance);
	 map_max_   = mv.max_density;
	 map_min_   = mv.min_density;
	 
	 set_initial_contour_level();
	 update_map();
	 ret_val = 1;
      }

      catch (const clipper::Message_base &rte) {
	 std::cout << "WARNING:: bad read of HKL file " << mtz_file_name << std::endl;
      }
   } else {
      std::cout << "No such file: " << mtz_file_name << std::endl;
   }
   return ret_val;
}


bool
molecule_class_info_t::make_patterson_using_intensities(std::string mtz_file_name,
							std::string i_col,
							std::string sigi_col,
							float map_sampling_rate) {

   bool ret_val = false;
   if (coot::file_exists(mtz_file_name)) {

      try {
	 clipper::HKL_info hkls, hklp;
	 clipper::CCP4MTZfile mtzin;
	 mtzin.open_read(mtz_file_name);
	 clipper::Spacegroup spgr = mtzin.spacegroup();
	 clipper::Resolution reso = mtzin.resolution();
	 clipper::Cell       cell = mtzin.cell();
	 clipper::HKL_data<clipper::data32::I_sigI> isig(hkls);
	 std::string data_name = "/*/*/[" + i_col + "," + sigi_col + "]";
	 mtzin.import_hkl_data(isig, data_name);
	 mtzin.close_read();
	 std::cout << "INFO:: make_patterson_using_intensities() Got "
		   << isig.num_obs() << " Is" << std::endl;

	 // get Patterson spacegroup
	 clipper::Spacegroup pspgr(clipper::Spgr_descr(spgr.generator_ops().patterson_ops()));
	 hklp.init(pspgr, cell, reso, true);
      
	 // make patterson coeffs
	 int count = 0;
	 clipper::HKL_data<clipper::data32::F_phi> fphi(hklp);
	 typedef clipper::HKL_data_base::HKL_reference_index HRI;
	 for (HRI ih = fphi.first(); !ih.last(); ih.next()) {
	    if (!isig[ih.hkl()].missing()) {
	       fphi[ih].f() =  isig[ih.hkl()].I();
	    } else {
	       fphi[ih].f() = 0.0;
	    }
	    fphi[ih].phi() = 0.0;
	 }
	 clipper::Grid_sampling grid;
	 grid.init(pspgr, cell, reso);
      
	 // make xmap
	 clipper::Xmap<float> xmap(pspgr, cell, grid);
	 xmap.fft_from(fphi);
	 is_patterson = 1; // needed for contour level protection
			   // (Pattersons are off-scale, but we should
			   // be able to change contour).

	 std::string map_name = "Patterson ";
	 map_name += mtz_file_name;
	 map_name += " ";
	 map_name += i_col;
	 new_map(xmap, map_name);
	 update_map_in_display_control_widget();
	 mean_and_variance<float> mv = map_density_distribution(xmap, 40, false);
	 map_mean_  = mv.mean; 
	 map_sigma_ = sqrt(mv.variance);
	 map_max_   = mv.max_density;
	 map_min_   = mv.min_density;
	 
	 set_initial_contour_level();
	 update_map();
	 ret_val = 1;
	 
      }
      catch (const std::exception &rte) {
	 std::cout << "WARNING:: something got caught: " << rte.what() << std::endl;
      }
      catch (const clipper::Message_fatal &rte) {
	 std::cout << "WARNING:: something got caught: " << rte.text() << std::endl;
      }
   } else {
      std::cout << "WARNING:: No such file: " << mtz_file_name << std::endl;
   }

   return ret_val;
}

