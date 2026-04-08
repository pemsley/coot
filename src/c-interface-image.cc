/*
 * src/c-interface-image.cc
 *
 * Copyright 2015 by Medical Research Council
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

#include <fstream>
#ifdef USE_PYTHON
#include <Python.h>  // before system includes to stop "POSIX_C_SOURCE" redefined problems
#endif

// #define RDKIT_HAS_CAIRO_SUPPORT

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#ifdef RDKIT_HAS_CAIRO_SUPPORT

#include <cairo.h>
#include <GraphMol/MolDraw2D/MolDraw2DCairo.h>
#include "lidia-core/rdkit-interface.hh"
#else
#include "lidia-core/rdkit-interface.hh"
#endif // RDKIT_HAS_CAIRO_SUPPORT


#include "utils/coot-utils.hh"
#endif // MAKE_ENHANCED_LIGAND_TOOLS

#include "graphics-info.h"
#include "globjects.h" //includes gtk/gtk.h

#include "c-interface-image-widget.hh"

#include "utils/logging.hh"
extern logging logger;

GtkWidget *test_get_image_widget_for_comp_id(const std::string &comp_id) {

   GtkWidget *r = 0;


#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#ifdef RDKIT_HAS_CAIRO_SUPPORT_xyz // doesn't compile now

   std::string smiles="CO[C@@H](O)C1=C(O[C@H](F)Cl)C(C#N)=C1ONNC[NH3+]";
   RDKit::ROMol *m_local = RDKit::SmilesToMol(smiles);
   TEST_ASSERT(m_local);
   RDDepict::compute2DCoords(*m_local);
   WedgeMolBonds(*m_local,&(m_local->getConformer()));
   std::string png_file_name = "image-" + comp_id + ".png";

   {
      RDKit::MolDraw2DCairo drawer(200,200);
      drawer.drawMolecule(*m_local);
      drawer.finishDrawing();
      std::string dt = drawer.getDrawingText();
      // std::cout << "PE-debug drawing-text :" << dt << ":" << std::endl;
      std::cout << "drawingtext is of length " << dt.length() << std::endl;
      drawer.writeDrawingText(png_file_name.c_str());
   }
#endif // RDKIT_HAS_CAIRO_SUPPORT
#endif // RDKIT

   return r;
}


GtkWidget *get_image_widget_for_comp_id(const std::string &comp_id, int imol, unsigned int image_size) {

   auto write_image_data = [] (const std::string &dt, const std::string &comp_id) {

      std::string fn = comp_id + ".png";
      std::ofstream f(fn);
      if (f) {
         f << dt;
         f.close();
      }
   };

   GtkWidget *r = 0;

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
#ifdef RDKIT_HAS_CAIRO_SUPPORT

   if (false)
      std::cout << "Here in get_image_widget_for_comp_id() with image size " << image_size << std::endl;

   graphics_info_t g;
   g.Geom_p()->try_dynamic_add(comp_id, g.cif_dictionary_read_number++);
   std::pair<bool, coot::dictionary_residue_restraints_t> dict =
      g.Geom_p()->get_monomer_restraints(comp_id, imol);

   if (dict.first) {

      try {
         RDKit::RWMol rdk_m = rdkit_mol(dict.second);
         coot::assign_formal_charges(&rdk_m);
         coot::rdkit_mol_sanitize(rdk_m);
         RDKit::RWMol rdk_mol_with_no_Hs = coot::remove_Hs_and_clean(rdk_m);

         int iconf_2d = RDDepict::compute2DCoords(rdk_mol_with_no_Hs);
         WedgeMolBonds(rdk_mol_with_no_Hs, &(rdk_mol_with_no_Hs.getConformer(iconf_2d)));
         bool debug_make_mol_file = false;
         if (debug_make_mol_file) {
            std::string smb = RDKit::MolToMolBlock(rdk_mol_with_no_Hs, true, iconf_2d);
            std::string fn = "test-" + comp_id + ".mol";
            std::ofstream f(fn.c_str());
            if (f)
               f << smb << std::endl;
            f.close();
         }

	 int n_conf = rdk_mol_with_no_Hs.getNumConformers();
	 // std::cout << "debug:: n_conf for " << comp_id << " is " << n_conf << std::endl;
	 if (n_conf > 0) {

            RDKit::MolDraw2DCairo drawer(image_size, image_size);
            RDKit::MolDrawOptions& draw_opts = drawer.drawOptions();
            // 20251211-PE we need to turn this off so that png can be parsed:
            // see https://greglandrum.github.io/rdkit-blog/posts/2025-02-15-storing-molecules-in-images.html
            // and https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html
            // eog and Image Viewer (on Arch Linux currently) don't display these images
            // if I don't remove metadata.
            // Offset	Keyword	  Content
            // 0x00ac9	rdkitPKL Binary RDKit pickle data
            // 0x00b6c	MOL	 RDKit MOL data
            // 0x00c5f	SMILES	 RDKit SMILES data
            draw_opts.includeMetadata = false;
            drawer.drawMolecule(rdk_mol_with_no_Hs);
            drawer.finishDrawing();
            std::string dt = drawer.getDrawingText();
            GError *error = NULL;

            // --- START: GTK4 Recommended Method ---
            // 1. Create a GBytes buffer from the std::string data
            //    g_bytes_new_static is efficient and uses the exact size of the string.
            GBytes *bytes = g_bytes_new_static(dt.c_str(), dt.length());

            // 2. Create the GdkTexture from the buffer
            GdkTexture *texture = gdk_texture_new_from_bytes(bytes, &error);
            g_bytes_unref(bytes); // We are done with the GBytes object

            if (texture) {
               // 3. Create the GtkImage widget from the texture
               r = gtk_image_new_from_paintable(GDK_PAINTABLE(texture));
               g_object_unref(texture); // We are done with the texture object
            } else {
               // Check the error if the texture failed to load
               // std::cout << "ERROR:: loading PNG texture: " << error->message << " for " << comp_id << std::endl;
               logger.log(log_t::ERROR, "loading PNG texture: " + std::string(error->message) + " for " + comp_id);
               g_error_free(error);
            }
         }
      }
      catch (...) {
         std::cout << "WARNING:: hack caught a ... exception " << std::endl;
      }

   } else {
      // std::cout << "No dictionary for rdkit_mol from " << comp_id << std::endl;
      // too noisy now the we don't distribute the monomer library with coot.
      //   logger.log(log_t::WARNING, "No dictionary for rdkit_mol from " + comp_id);
   }

#endif   // RDKIT_HAS_CAIRO_SUPPORT
#endif   // MAKE_ENHANCED_LIGAND_TOOLS
   return r;
}

