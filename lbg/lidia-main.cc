/* lbg/lidia-main.cc
 * 
 * Author: Paul Emsley
 * Copyright 2010, 2011, 2012 by The University of Oxford
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
#include <Python.h>
#endif

#include <iostream>

#ifdef HAVE_GOOCANVAS

#include <sys/types.h>  // for stating
#include <sys/stat.h>
#include <unistd.h>
#include <sys/time.h>
#include <stdexcept>
#include <fstream>
#include <iomanip>
#include <algorithm>

#include "compat/coot-sysdep.h"

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#include <RDGeneral/versions.h>
#include <RDGeneral/FileParseException.h>
#include <RDGeneral/BadFileException.h>
#include <GraphMol/FileParsers/FileParsers.h>
#ifdef USE_PYTHON
#include <boost/python.hpp>
#endif
#endif

#include <cairo.h>
#if CAIRO_HAS_PDF_SURFACE
#include <cairo-pdf.h>
#endif
#include "lbg.hh"

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
void set_show_additional_representation(int imol, int representation_number, int on_off_flag) {}
void set_show_all_additional_representations(int imol, int on_off_flag) {}
void orient_view(int imol,
		 const coot::residue_spec_t &central_residue_spec,
		 const coot::residue_spec_t &neighbour_residue_spec) {}
void all_additional_representations_off_except(int imol, int addrep,
					       short int ball_and_sticks_off_too_flag) {}
void set_rotation_centre(const clipper::Coord_orth &pos) {}
#endif


int
main(int argc, char *argv[]) {

   mmdb::InitMatType(); // mmdb program. 

#ifdef USE_PYTHON
   std::cout << "FIXME! main()" << std::endl;
#if 0
   Py_Initialize();
   PySys_SetArgv(argc, argv);
   PyRun_SimpleString("global user_defined_alert_smarts ; user_defined_alert_smarts = []");
#endif   
#endif
   
   gtk_init (&argc, &argv);
   std::string molecule_file_name = "";
   lig_build::molfile_molecule_t mm;
   mmdb::Manager *mol = NULL; // no atom names to transfer

   // don't use my mol reader here, use the RDKit one after we have started lbg.
   // Currently command line file reading is disabled then.

#if MAKE_ENHANCED_LIGAND_TOOLS
   if (argc > 1) {
      std::string file_name(argv[1]);
      try {
	 RDKit::RWMol *rdkm = RDKit::MolFileToMol(file_name);
	 if (rdkm) {

	    RDKit::MolOps::Kekulize(*rdkm);
	    double weight_for_3d_distances = 0.1;
	    int iconf = coot::add_2d_conformer(rdkm, weight_for_3d_distances);
	    mm = coot::make_molfile_molecule(*rdkm, iconf);
	 } else {
	    std::cout << "WARNING:: file " << file_name << " null molecule" << std::endl;
	 }
      }
      catch (const RDKit::MolSanitizeException &e) {
	 std::cout << "WARNING:: file " << file_name << " "  << e.what() << std::endl;
      }
      catch (const RDKit::BadFileException &e) {
	 std::cout << "WARNING:: Bad file " << file_name << " "  << e.what() << std::endl;
      }
      catch (const RDKit::FileParseException &e) {
	 std::cout << "WARNING:: File Parse error " << file_name << " "  << e.what() << std::endl;
      }
   }
#else
    if (argc > 1) {
       std::string file_name(argv[1]);
       mm.read(file_name);
    }
#endif

   bool stand_alone_flag = 1;
   bool use_graphics = 1;
   std::string view_name = "";
   std::pair<bool, coot::residue_spec_t> p(0, coot::residue_spec_t());
   int imol = -1; // dummy/unset
   int (*func) (const char *s1, const char *s2) = NULL;
   void (*prodrg_import_func_ptr) (std::string file_name, std::string comp_id) = NULL;
   void (*sbase_import_func_ptr)  (std::string file_name) = NULL;
   std::string (*get_drug_mdl_via_wikipedia_and_drugbank) (std::string drug_name) = NULL;
   coot::protein_geometry *geom_p = NULL;
   
#if ( ( (GTK_MAJOR_VERSION == 2) && (GTK_MINOR_VERSION > 11) ) || GTK_MAJOR_VERSION > 2)
   if (lbg(mm, p, mol, view_name, molecule_file_name, imol,
	   geom_p,
	   use_graphics, stand_alone_flag,
	   func,
	   prodrg_import_func_ptr,
	   sbase_import_func_ptr,
	   get_drug_mdl_via_wikipedia_and_drugbank)) {
       gtk_main ();
   } 
#endif // GTK_VERSION
   return 1;
}


#else

int
main(int argc, char *argv[]) {

   std::cout << "No goo canvas at compile-time, no lidia " << std::endl;
   return 0;
}



#endif 
