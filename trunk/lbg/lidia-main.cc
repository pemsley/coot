/* lbg/lidia-main.cc
 * 
 * Author: Paul Emsley
 * Copyright 2010 by The University of Oxford
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

#include <sys/types.h>  // for stating
#include <sys/stat.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdexcept>
#include <fstream>
#include <iomanip>
#include <algorithm>


#include <cairo.h>
#if CAIRO_HAS_PDF_SURFACE
#include <cairo-pdf.h>
#endif
#include "lbg.hh"

#ifdef MAKE_ENTERPRISE_TOOLS
void set_show_additional_representation(int imol, int representation_number, int on_off_flag) {}
void set_show_all_additional_representations(int imol, int on_off_flag) {}
void orient_view(int imol,
		 const coot::residue_spec_t &central_residue_spec,
		 const coot::residue_spec_t &neighbour_residue_spec) {}
void all_additional_representations_off_except(int imo, int addrep) {}
#endif


int
main(int argc, char *argv[]) {

   InitMatType(); // mmdb program. 
   
   gtk_init (&argc, &argv);
   std::string molecule_file_name = "";
   lig_build::molfile_molecule_t mm;
   CMMDBManager *mol = NULL; // no atom names to transfer

   if (argc > 1) {
      std::string file_name(argv[1]);
      mm.read(file_name);
   }

   bool stand_alone_flag = 1;
   std::string view_name = "";
   std::pair<bool, coot::residue_spec_t> p(0, coot::residue_spec_t());
   int imol = -1; // dummy/unset
   if (lbg(mm, p, mol, view_name, molecule_file_name, imol, stand_alone_flag)) {
       gtk_main ();
   } 
   return 1;
}
