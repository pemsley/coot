/*
 * src/graphics-info-init.cc
 *
 * Copyright 2018 by Medical Research Council
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

#include "utils/coot-utils.hh"
#ifdef USE_PYTHON
#include <Python.h>
#endif

#include "coords/cos-sin.h"
#include "graphics-info.h"

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

void
graphics_info_t::load_freetype_font_textures() {

   // std::cout << "------------------------------- load_freetype_font_textures() -------" << std::endl;

   // ----------------------------- font test -----------------
   FT_Library ft;
   if (FT_Init_FreeType(&ft))
   std::cout << "ERROR::FREETYPE: Could not init FreeType Library" << std::endl;

   FT_Face face;
   std::string pkgdatadir = coot::package_data_dir();
   std::string font_dir  = coot::util::append_dir_dir(pkgdatadir, "fonts");
   std::string font_path = coot::util::append_dir_file(font_dir, "Vera.ttf");
   if (! coot::file_exists(font_path))
       font_path = "fonts/Vera.ttf";
   if (FT_New_Face(ft, font_path.c_str(), 0, &face)) {
      std::cout << "ERROR::FREETYPE: Failed to load font" << std::endl;
   } else {
      vera_font_loaded = true;
      // FT_Set_Pixel_Sizes(face, 0, 24); too big for labels
      // FT_Set_Pixel_Sizes(face, 0, 16);
      FT_Set_Pixel_Sizes(face, 0, 40);

      glPixelStorei(GL_UNPACK_ALIGNMENT, 1); // Disable byte-alignment restriction
      // only using one byte.

      for (GLubyte ic = 0; ic < 255; ic++) {
         // Load character glyph
         if (FT_Load_Char(face, ic, FT_LOAD_RENDER)) {
            std::cout << "ERROR::FREETYTPE: Failed to load Glyph " << ic << std::endl;
            continue;
         }
         // Generate texture
         GLuint texture;
         glGenTextures(1, &texture);
         GLenum err = glGetError();
         if (err) std::cout << "Loading characture textures glGenTextures err " << err << std::endl;
         // if (!err) std::cout << "OK loading texture " << texture << " for ic " << ic << std::endl;
         glBindTexture(GL_TEXTURE_2D, texture);
         err = glGetError(); if (err) std::cout << "Loading characture textures glBindTexture err " << err << std::endl;
         glTexImage2D( GL_TEXTURE_2D,
                       0,
                       GL_RED,
                       face->glyph->bitmap.width,
                       face->glyph->bitmap.rows,
                       0,
                       GL_RED,
                       GL_UNSIGNED_BYTE,
                       face->glyph->bitmap.buffer);
         // std::cout << "load_freetype_font_textures(): character " << ic << " " << face->glyph->bitmap.width
         // << " " << face->glyph->bitmap.rows << std::endl;
         // Set texture options
         err = glGetError(); if (err) std::cout << "Loading characture textures glTexImage2D err " << err << std::endl;
         glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
         glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
         glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
         glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
         err = glGetError(); if (err) std::cout << "Loading characture textures glTexParameteri err " << err << std::endl;

         // std::cout << "Storing characture with texture id " << texture << std::endl;
         // Now store the character
         GLuint face_glyph_advance_x = face->glyph->advance.x;
         FT_character character = {
                                   texture,
                                   glm::ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
                                   glm::ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
                                   face_glyph_advance_x
         };
         ft_characters[ic] = character; // ic type change
      }
      glBindTexture(GL_TEXTURE_2D, 0);
      // Destroy FreeType once we're finished
      FT_Done_Face(face);
      FT_Done_FreeType(ft);
   }
}

void
graphics_info_t::init() {

#ifdef WINDOWS_MINGW
   prefer_python = 1;
#endif
   // The cosine->sine lookup table, used in picking.
   //
   // The data in it are static, so we can get to them anywhere
   // now that we have run this
   cos_sin cos_sin_table(1000);


      // transform = Transform(glm::vec3(0.0, 0.0, 0.0),
      //                       glm::vec3(0.0, 0.0, 0.0),
      //                       glm::vec3(1.0, 1.0, 1.0));

      // camera = benny::Camera(glm::vec3(0.0f,0.0f,-3.5f), 70, aspect, 0.1f, 500.0f);

      on_going_updating_map_lock = false;

      find_ligand_ligand_mols_ = new std::vector<std::pair<int, bool> >;
      geom_p = new coot::protein_geometry;
      geom_p->set_verbose(true); // was false

      use_gemmi = false;
      cif_dictionary_read_number = geom_p->init_standard();

      geom_p->add_planar_peptide_restraint();
      convert_dictionary_planes_to_improper_dihedrals_flag = false;

      use_harmonic_approximation_for_NBCs = false; // Hard-mode by default

      geom_p->init_ccp4srs("srsdata"); // overridden by COOT_CCP4SRS_DIR and CCP4_LIB

      // rotamer probabilitiles
      // guess we shall rather use COOT_DATA_DIR and only as fallback PKGDATADIR?!
      // maybe only for windows!?
      //
      // 20090920-PE, no, not only windows.  If they set
      // COOT_DATA_DIR, let's use that instead of PKGDATADIR (useful
      // for Justin Lecher and Gentoo who test before installing (and
      // they need a way to specify the data dir (before installing
      // it's not in PKGDATADIR)).
      //
      std::string tables_dir = coot::package_data_dir();

      char *data_dir = getenv("COOT_DATA_DIR");
      if (data_dir) {
	      tables_dir = data_dir;
      }

      tables_dir += "/rama-data";
      rot_prob_tables.set_tables_dir(tables_dir);

      moving_atoms_asc = new atom_selection_container_t;
      moving_atoms_asc->mol = NULL;
      moving_atoms_asc->atom_selection = NULL;
      moving_atoms_asc->n_selected_atoms = 0;

      standard_residues_asc.read_success = 0;
      standard_residues_asc.n_selected_atoms = 0;
      read_standard_residues(); // updates read_success

      symmetry_colour_merge_weight = 0.5; // 0.0 -> 1.0

      // use_graphics_interface_flag = 1;  don't (re)set this here,
      // it is set as a static and possibly modified by immediate
      // handling of command line data in main.cc

      // moving_atoms_asc gets filled in copy_mol_and_regularize, not
      // here.

      // db_main = NULL;

      // command line scripts:
      // command_line_scripts = new std::vector<std::string>; no longer 20210923-PE

      console_display_commands.display_commands_flag = false;

      // LSQ matching info
      lsq_matchers = new std::vector<coot::lsq_range_match_info_t>;

      // LSQ Plane
      lsq_plane_atom_positions = new std::vector<clipper::Coord_orth>;

      directory_for_fileselection = "";
      directory_for_filechooser = "";

      // rotamer distortion graph scale
      rotamer_distortion_scale = 0.3;

      // cif dictionary
      cif_dictionary_filename_vec = new std::vector<std::string>;

      std::string pdd = coot::package_data_dir();
      std::filesystem::path pdd_path(pdd);
      std::filesystem::path ptm_database_json = pdd_path / "data" / "ptm_database.json";
      ptm_database.read(ptm_database_json.string());

      // 2026-02-06-PE
      // put the inchikey reader here.

      // 2026-02-06-PE is this still a thing?
      // ligand blobs:
      ligand_big_blobs = new std::vector<clipper::Coord_orth>;

      // rot_trans adjustments:
      for (int i=0; i<6; i++)
	 previous_rot_trans_adjustment[i] = -10000;

      // merging molecules
      merge_molecules_merging_molecules = new std::vector<int>;

      // generic display objects
      // generic_objects_p = new std::vector<coot::old_generic_display_object_t>;
      generic_objects_dialog = NULL;

      // views
      // views = new std::vector<coot::view_info_t>; // not a pointer any more - hooray.

      // glob extensions:
      coordinates_glob_extensions = new std::vector<std::string>;
      data_glob_extensions = new std::vector<std::string>;
      map_glob_extensions = new std::vector<std::string>;
      dictionary_glob_extensions  = new std::vector<std::string>;

      coordinates_glob_extensions->push_back(".pdb");
      coordinates_glob_extensions->push_back(".pdb.gz");
      coordinates_glob_extensions->push_back(".brk");
      coordinates_glob_extensions->push_back(".brk.gz");
      coordinates_glob_extensions->push_back(".ent");
      coordinates_glob_extensions->push_back(".ent.gz");
      coordinates_glob_extensions->push_back(".ent.Z");
      coordinates_glob_extensions->push_back(".cif");
      coordinates_glob_extensions->push_back(".mmcif");
      coordinates_glob_extensions->push_back(".mmCIF");
      coordinates_glob_extensions->push_back(".cif.gz");
      coordinates_glob_extensions->push_back(".mmcif.gz");
      coordinates_glob_extensions->push_back(".mmCIF.gz");
      coordinates_glob_extensions->push_back(".res");  // SHELX
      coordinates_glob_extensions->push_back(".ins");  // SHELX
      coordinates_glob_extensions->push_back(".pda");  // SHELX

      data_glob_extensions->push_back(".mtz");
      data_glob_extensions->push_back(".hkl");
      data_glob_extensions->push_back(".data");
      data_glob_extensions->push_back(".phs");
      data_glob_extensions->push_back(".pha");
      data_glob_extensions->push_back(".cif");
      data_glob_extensions->push_back(".fcf"); // SHELXL
      data_glob_extensions->push_back(".mmcif");
      data_glob_extensions->push_back(".mmCIF");
      data_glob_extensions->push_back(".cif.gz");
      data_glob_extensions->push_back(".mmcif.gz");
      data_glob_extensions->push_back(".mmCIF.gz");

      map_glob_extensions->push_back(".map");
      map_glob_extensions->push_back(".mrc");
      map_glob_extensions->push_back(".map.gz");
      map_glob_extensions->push_back(".mrc.gz");
      map_glob_extensions->push_back(".ext");
      map_glob_extensions->push_back(".msk");
      map_glob_extensions->push_back(".ccp4");
      map_glob_extensions->push_back(".cns");

      dictionary_glob_extensions->push_back(".cif");
      dictionary_glob_extensions->push_back(".mmcif");
      dictionary_glob_extensions->push_back(".mmCIF");
      dictionary_glob_extensions->push_back(".cif.gz");
      dictionary_glob_extensions->push_back(".mmcif.gz");
      dictionary_glob_extensions->push_back(".mmCIF.gz");
      dictionary_glob_extensions->push_back(".lib");

      /* things for preferences */

      // -------------these are all frame ------------------------------

      //preferences_internal = new std::vector<coot::preference_info_t>;

      // preferences_general_tabs.push_back("preferences_file_selection");
      // preferences_general_tabs.push_back("preferences_dock_accept_dialog");
      preferences_general_tabs.push_back("preferences_hid");
      preferences_general_tabs.push_back("preferences_noughties_physics");
      preferences_general_tabs.push_back("preferences_recentre_pdb");
      preferences_general_tabs.push_back("preferences_smooth_scroll");
      // preferences_general_tabs.push_back("preferences_model_toolbar_style");
      // preferences_general_tabs.push_back("preferences_main_toolbar_style");

      preferences_bond_tabs.push_back("preferences_bond_parameters");
      preferences_bond_tabs.push_back("preferences_bond_colours");

      preferences_map_tabs.push_back("preferences_map_parameters");
      preferences_map_tabs.push_back("preferences_map_colours");
      preferences_map_tabs.push_back("preferences_map_drag");

      preferences_geometry_tabs.push_back("preferences_cis_peptides");
      preferences_geometry_tabs.push_back("preferences_default_b_factor_frame");

      preferences_colour_tabs.push_back("preferences_background_colour");
      preferences_colour_tabs.push_back("preferences_bond_colours");
      preferences_colour_tabs.push_back("preferences_map_colours");

      preferences_other_tabs.push_back("preferences_console");
      preferences_other_tabs.push_back("preferences_tips");
      preferences_other_tabs.push_back("preferences_speed");
      preferences_other_tabs.push_back("preferences_antialias");
      preferences_other_tabs.push_back("preferences_font");
      preferences_other_tabs.push_back("preferences_pink_pointer");

      // for toolbar icons in preferences
      model_toolbar_icons = new std::vector<coot::preferences_icon_info_t>;
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(0, "refine-1.svg",
								   "Real Space Refine",
								   "model_toolbar_refine_togglebutton",
								   1, 1));
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(1, "regularize-1.svg",
								   "Regularize",
								   "model_toolbar_regularize_togglebutton",
								   1, 1));
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(2, "anchor.svg",
								   "Fixed Atoms...",
								   "model_toolbar_fixed_atoms_button",
								   1, 1));
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(3, "rigid-body.svg",
								   "Rigid Body Fit Zone",
								   "model_toolbar_rigid_body_fit_togglebutton",
								   1, 1));
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(4, "rtz.svg",
								   "Rotate/Translate Zone",
								   "model_toolbar_rot_trans_toolbutton",
								   1, 1));
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(5, "auto-fit-rotamer.svg",
								   "Auto Fit Rotamer",
								   "model_toolbar_auto_fit_rotamer_togglebutton",
								   1, 1));
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(6, "rotamers.svg",
								   "Rotamers",
								   "model_toolbar_rotamers_togglebutton",
								   1,1 ));
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(7, "edit-chi.svg",
								   "Edit Chi Angles",
								   "model_toolbar_edit_chi_angles_togglebutton",
								   1, 1));
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(8, "torsion-general.svg",
								   "Torsion General",
								   "model_toolbar_torsion_general_toggletoolbutton",
								   1, 0));
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(9, "flip-peptide.svg",
								   "Flip Peptide",
								   "model_toolbar_flip_peptide_togglebutton",
								   1, 1));
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(10, "side-chain-180.svg",
								   "Side Chain 180 Degree Flip",
								   "model_toolbar_sidechain_180_togglebutton",
								   1, 1));
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(11, "edit-backbone.svg",
								   "Edit Backbone Torsions",
								   "model_toolbar_edit_backbone_torsions_toggletoolbutton",
								   1, 1));
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(12, "",
								   "---------------------",
								   "model_toolbar_hsep_toolitem",
								   1, 1));
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(13, "mutate-auto-fit.svg",
								   "Mutate and Auto-Fit...",
								   "model_toolbar_mutate_and_autofit_togglebutton",
								   1, 1));
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(14, "mutate.svg",
								   "Simple Mutate...",
								   "model_toolbar_simple_mutate_togglebutton",
								   1, 1));
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(15, "add-peptide-1.svg",
								   "Add Terminal Residue...",
								   "model_toolbar_add_terminal_residue_togglebutton",
								   1, 1));
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(16, "add-alt-conf.svg",
								   "Add Alt Conf...",
								   "model_toolbar_add_alt_conf_toolbutton",
								   1, 1));
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(17, "atom-at-pointer.svg",
								   "Place Atom at Pointer",
								   "model_toolbar_add_atom_button",
								   1, 1));
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(18, "gtk-clear",
								   "Clear Pending Picks",
								   "model_toolbar_clear_pending_picks_button",
								   1, 1));
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(19, "gtk-delete",
								   "Delete...",
								   "model_toolbar_delete_button",
								   1, 1));
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(20, "gtk-undo",
								   "Undo",
								   "model_toolbar_undo_button",
								   1, 1));
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(21, "gtk-redo",
								   "Redo",
								   "model_toolbar_redo_button",
								   1, 1));
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(22, "",
								   "---------------------",
								   "model_toolbar_hsep_toolitem2",
								   1, 0));
      model_toolbar_icons->push_back(coot::preferences_icon_info_t(23, "azerbaijan.svg",
								   "Run Refmac...",
								   "model_toolbar_refmac_button",
								   1, 0));
      // for main icons in preferences
      main_toolbar_icons = new std::vector<coot::preferences_icon_info_t>;
      main_toolbar_icons->push_back(coot::preferences_icon_info_t(0, "gtk-open",
								  "Open Coords...",
								  "coords_toolbutton",
								  1, 1));
      main_toolbar_icons->push_back(coot::preferences_icon_info_t(1, "gtk-zoom-fit",
								  "Reset View",
								  "reset_view_toolbutton",
								  1, 1));
      main_toolbar_icons->push_back(coot::preferences_icon_info_t(2, "display-manager.png",
								  "Display Manager",
								  "display_manager_toolbutton",
								  1, 1));
      main_toolbar_icons->push_back(coot::preferences_icon_info_t(3, "go-to-atom.svg",
								  "Go To Atom...",
								  "go_to_atom_toolbutton",
								  1, 1));
      main_toolbar_icons->push_back(coot::preferences_icon_info_t(4, "go-to-ligand.svg",
								  "Go To Ligand",
								  "go_to_ligand_toolbutton",
								  1, 1));

      do_expose_swap_buffers_flag = 1;
      vera_font_loaded = false;

      regenerate_bonds_needs_make_bonds_type_checked_flag = true;

      refmac_dialog_mtz_file_label = NULL;
      /* set no of refmac cycles */
      preset_number_refmac_cycles = new std::vector<int>;
      preset_number_refmac_cycles->push_back(0);
      preset_number_refmac_cycles->push_back(1);
      preset_number_refmac_cycles->push_back(2);
      preset_number_refmac_cycles->push_back(3);
      preset_number_refmac_cycles->push_back(5);
      preset_number_refmac_cycles->push_back(7);
      preset_number_refmac_cycles->push_back(10);
      preset_number_refmac_cycles->push_back(15);
      preset_number_refmac_cycles->push_back(20);
      preset_number_refmac_cycles->push_back(50);

      validation_graph_model_list = gtk_list_store_new(2,G_TYPE_STRING,G_TYPE_INT);
      active_validation_graph_model_idx = -1;

      ramachandran_plot_model_list = gtk_list_store_new(2,G_TYPE_STRING,G_TYPE_INT);

   }

