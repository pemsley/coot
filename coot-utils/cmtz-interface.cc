/*
 * coot-utils/cmtz-interface.cc
 *
 * Copyright 2023 by Medical Research Council
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
 * Lesser General Public License for more details.
 *
 * You should have received a copies of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,
 * Fifth Floor, Boston, MA, 02110-1301, USA.
 * See http://www.gnu.org/licenses/
 *
 */


#include <clipper/ccp4/ccp4_mtz_io.h>
#include "cmtz-interface.hh"
#include "utils/logging.hh"

extern logging logger;

coot::mtz_column_types_info_t
coot::get_mtz_columns(const std::string &filename) {

   coot::mtz_column_types_info_t a;
   a.read_success = 0;
   a.selected_f_col      = 0;
   a.selected_phi_col    = -1; /* unset */
   a.selected_weight_col = 0;
   a.selected_refmac_fobs_col = 0;
   a.selected_refmac_sigfobs_col = 0;
   a.selected_refmac_r_free_col = 0;
   a.selected_refmac_phi_col = -1; /* unset */
   a.selected_refmac_fom_col = 0;
   a.selected_refmac_hla_col = -1; /* unset */
   a.selected_refmac_hlb_col = 0;
   a.selected_refmac_hlc_col = 0;
   a.selected_refmac_hld_col = 0;
   // new ones for SAD and twin refinement
   // I, I+/-, F+/- and appropriate sigs
   a.selected_refmac_fp_col = 0;
   a.selected_refmac_sigfp_col = 0;
   a.selected_refmac_fm_col = 0;
   a.selected_refmac_sigfm_col = 0;
   a.selected_refmac_iobs_col = -1; /* unset */
   a.selected_refmac_sigiobs_col = 0;
   a.selected_refmac_ip_col = 0;
   a.selected_refmac_sigip_col = 0;
   a.selected_refmac_im_col = 0;
   a.selected_refmac_sigim_col = 0;

   clipper::CCP4MTZfile f;
   short int is_mtz_file = 1;
   // new try catch here

   try {
      f.open_read(filename);
   }

   catch (...) {
      logger.log(log_t::INFO, logging::ltw("not an mtz file: "), logging::ltw(filename));
      // std::cout << "INFO:: not an mtz file: " << filename << std::endl;
      is_mtz_file = 0;
   }

   if (is_mtz_file) {
      std::vector<clipper::String> v = f.column_labels();
      // std::cout << "INFO:: found " << v.size() << " column labels in " << filename << "\n";
      if (v.size() > 1) {
         a.read_success = 1;
         a.mtz_filename = filename;
         for (unsigned int i=0; i<v.size(); i++) {
            // std::cout << i << " " << v[i] << "\n";
            std::string label;
            std::string type;
            std::string::size_type ispace = v[i].find_last_of(" ");
            if (ispace == std::string::npos) {
               std::cout <<  "WARNING:: uninterprettable label \"" << v[i] << "\" of "
                         << filename << "\n";
            } else {
               label = v[i].substr(0, ispace);
               type  = v[i].substr(ispace+1);
               //std::cout << "Got label :" << label << ": and type :" << type << ":\n";
               if (type == "F")
                 a.f_cols.push_back(coot::mtz_type_label(label, 'F', i));
               if (type == "G")
                 a.fpm_cols.push_back(coot::mtz_type_label(label, 'G', i));
               if (type == "L")
                 a.sigfpm_cols.push_back(coot::mtz_type_label(label, 'L', i));
               if (type == "Q")
                 a.sigf_cols.push_back(coot::mtz_type_label(label, 'Q', i));
               if (type == "P")
                 a.phi_cols.push_back(coot::mtz_type_label(label, 'P', i));
               if (type == "D")
                 a.d_cols.push_back(coot::mtz_type_label(label, 'D', i));
               if (type == "W")
                 a.weight_cols.push_back(coot::mtz_type_label(label, 'W', i));
               if (type == "I")
                 a.r_free_cols.push_back(coot::mtz_type_label(label, 'I', i));
               if (type == "A")
                 a.hl_cols.push_back(coot::mtz_type_label(label, 'A', i));
               if (type == "J")
                 a.i_cols.push_back(coot::mtz_type_label(label, 'J', i));
               // for completeness; not used yet
               if (type == "K")
                 a.ipm_cols.push_back(coot::mtz_type_label(label, 'K', i));
               if (type == "M")
                 a.sigipm_cols.push_back(coot::mtz_type_label(label, 'M', i));
            }
         }
      }
   }
   return a;
}
