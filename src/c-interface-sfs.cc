/* src/c-interface-build.cc
 * 
 * Copyright 2012 The University of Oxford
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

#include <sys/types.h>
#include <sys/stat.h>

#include <gtk/gtk.h>
#include "c-interface.h"

#include "clipper/cif/cif_data_io.h"
#include "clipper/ccp4/ccp4_mtz_io.h"

#include "coot-utils.hh"

/* ------------------------------------------------------------------------- */
/*                      mmCIF sfs -> MTZ conversion                          */
/* ------------------------------------------------------------------------- */

/*! \brief convert the structure factors in cif_file_name to an mtz
  file.  Return 1 on success. */
int mmcif_sfs_to_mtz(const char *cif_file_name, const char *mtz_out_file_name) {

   int r = -1;

   clipper::HKL_data<clipper::datatypes::F_sigF    <float> > f_sigf;
   clipper::HKL_data<clipper::datatypes::Phi_fom   <float> > phi_fom;
   clipper::HKL_data<clipper::datatypes::F_sigF_ano<float> > f_sigf_anom;
   clipper::HKL_data<clipper::datatypes::I_sigI    <float> > i_sigi;
   clipper::HKL_data<clipper::datatypes::I_sigI_ano<float> > i_sigi_anom;
   clipper::HKL_data<clipper::datatypes::D_sigD    <float> > d_sigd;
   clipper::HKL_data<clipper::datatypes::ABCD      <float> > ABCD;
   clipper::HKL_data<clipper::datatypes::Flag> free_r_flags;
   clipper::HKL_info hkl_info;
   clipper::CIFfile cif;

   try { 
      cif.open_read (cif_file_name);
      cif.import_hkl_info(hkl_info);
      cif.import_hkl_data(f_sigf);
      cif.import_hkl_data(f_sigf_anom);
      cif.import_hkl_data(i_sigi);
      cif.import_hkl_data(i_sigi_anom);
      cif.import_hkl_data(d_sigd);
      cif.import_hkl_data(ABCD);
      cif.import_hkl_data(free_r_flags);
      cif.import_hkl_data(phi_fom);
      cif.close_read();

      if (hkl_info.is_null()) {
	 std::cout << "ERROR:: null HKL info for " << cif_file_name << std::endl;
      } else { 
	 clipper::Spacegroup sg = f_sigf.spacegroup();
	 if (sg.is_null()) {
	    std::cout << "ERROR:: null space group for "
		      << cif_file_name << std::endl;
	 } else {
	    clipper::CCP4MTZfile mtz;
	    mtz.open_write(mtz_out_file_name);
	    mtz.export_hkl_info(hkl_info);
	    if (i_sigi.num_obs())
	       mtz.export_hkl_data(i_sigi,       "/*/*/I" );
	    if (i_sigi_anom.num_obs())
	       mtz.export_hkl_data(i_sigi_anom,  "/*/*/I_anom" );
	    if (f_sigf.num_obs())
	       mtz.export_hkl_data(f_sigf,       "/*/*/F" );
	    if (f_sigf_anom.num_obs())
	       mtz.export_hkl_data(f_sigf_anom,  "/*/*/F_anom" );
	    if (phi_fom.num_obs())
	       mtz.export_hkl_data(phi_fom,      "/*/*/phi_fom" );
	    if (d_sigd.num_obs())
	       mtz.export_hkl_data(d_sigd,       "/*/*/D" );
	    if (ABCD.num_obs())
	       mtz.export_hkl_data(ABCD,         "/*/*/HL_ABCD" );
	    if (free_r_flags.num_obs())
	       mtz.export_hkl_data(free_r_flags, "/*/*/Rfree" );
	    mtz.close_write();
	    
	    // OK, except if we don't have R-free (that happens on
	    // broken files (rare) and old versions of clipper).
	    // 
	    if (free_r_flags.num_obs())
	       r = 1;
	    else
	       r = 0;
	 }
      }
   }
   
   catch (clipper::Message_fatal m) {
      // the error message has already been printed
   } 
   return r;
} 
