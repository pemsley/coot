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
   
   if (lbg(mm, mol, molecule_file_name)) {
       gtk_main ();
   } 
   return 1;
}
