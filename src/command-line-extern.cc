/* src/command-line.cc
 * 
 * Copyright 2001, 2002, 2003, 2004, 2005, 2006 by The University of York
 * Copyright 2015 by Medical Research Council
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */


#if defined USE_PYTHON 
#include "Python.h"  // before system includes to stop "POSIX_C_SOURCE" redefined problems
// BL says:: need to exclude in Windows as conflicting getopt definitions.
//           do we need python here at all?!
#endif


#ifdef _MSC_VER
#include <windows.h>
#else
#include <unistd.h> // for getopt(3)
#endif

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif


#include "graphics-info.h"
// Including python needs to come after graphics-info.h, because
// something in Python.h (2.4 - chihiro) is redefining FF1 (in
// ssm_superpose.h) to be 0x00004000 (Grrr).
// BL says:: and (2.3 - dewinter), i.e. is a Mac - Python issue
// since the follwing two include python graphics-info.h is moved up

#include "c-interface.h"
#include "cc-interface.hh"
#include "coot-version.hh"
#include "command-line.hh"
#include "get-monomer.hh"



extern "C"
void
handle_command_line_data_argc_argv(int argc, char **argv) {

   command_line_data cld = parse_command_line(argc, argv);
   handle_command_line_data(cld);
} 


extern "C"
void
handle_command_line_data(command_line_data cld) {


   // We *should* run scripts first and they can make setting that
   // affect the other command line options (e.g. column labels for
   // auto-read).
   //
   // But scripts are run in c_inner_main().

   // script

   // We need to setup guile before we run the script.
   //
   // i.e. handle_command_line args needs to go "in" scm_boot_guile
   // 
   // OK we do that using run_command_line_scripts() which runs
   // commands stored in graphics_info_t::command_line_scripts.  So
   // store stuff there.
   // 
   for (unsigned int i=0; i< cld.script.size(); i++) {
      graphics_info_t::command_line_scripts->push_back(cld.script[i]);
   }

   // command line scripting (direct using -c)
#ifdef USE_GUILE
   if (cld.script_is_python_flag)
      graphics_info_t::command_line_commands.is_python = 1;
#else
#ifdef USE_PYTHON
   // Use python automatically if there is only python available
   graphics_info_t::command_line_commands.is_python = 1;
#endif //USE_PYTHON
#endif // USE_GUILE
   for (unsigned int i=0; i<cld.command.size(); i++) {
      graphics_info_t::command_line_commands.commands.push_back(cld.command[i]);
   }

   // getting map model by passing the accession code:
   for (unsigned int i=0; i<cld.accession_codes.size(); i++) {
      graphics_info_t::command_line_accession_codes.push_back(cld.accession_codes[i]);
   }


   if (cld.em_mode) {
      graphics_info_t::box_radius_xray = graphics_info_t::box_radius_em;
      // what else?
   }

   // coordinates

   for (unsigned int i=0; i< cld.coords.size(); i++) {
      // don't slide around for 100 ligands
      short int smooth_scroll_on_state_pre = graphics_info_t::smooth_scroll_on;
      short int smooth_scroll_state_pre = graphics_info_t::smooth_scroll;
      graphics_info_t::smooth_scroll_on = 0;
      handle_read_draw_molecule(cld.coords[i].c_str());
      graphics_info_t::smooth_scroll_on = smooth_scroll_on_state_pre;
      graphics_info_t::smooth_scroll    = smooth_scroll_state_pre;
   }

   // datasets

   for (unsigned int i=0; i< cld.datasets.size(); i++) {
      std::cout << "debug: manage_column_selector for file: "
	   << cld.datasets[i].c_str() << std::endl;
      manage_column_selector(cld.datasets[i].c_str());
   }

   // auto-datasets

   for (unsigned int i=0; i< cld.auto_datasets.size(); i++) { 
      auto_read_make_and_draw_maps(cld.auto_datasets[i].c_str()); 
   }

   // maps

   for (unsigned int i=0; i< cld.maps.size(); i++) { 
      handle_read_ccp4_map(cld.maps[i], 0); // not difference map
   }

   // emdb codes
   for (unsigned int i=0; i< cld.emdb_codes.size(); i++) { 
      handle_read_emdb_data(cld.emdb_codes[i]); // not difference map
   }

   // cif dictionaries
   
   for (unsigned int i=0; i< cld.dictionaries.size(); i++) {
      read_cif_dictionary(cld.dictionaries[i].c_str());
   }

   for (unsigned int i=0; i<cld.comp_ids.size(); i++) { 
      get_monomer(cld.comp_ids[i].c_str());
   }

   // ccp4 project directory given?
   if (cld.ccp4_project != "") {
      std::string dir = ccp4_project_directory(cld.ccp4_project);
      if (dir != "") {
	 graphics_info_t g;
	 g.set_directory_for_fileselection_string(dir);
      }
   }

   // title
   if (cld.title.length() > 0)
      set_main_window_title(cld.title.c_str());

   // --no-guano used?
   if (cld.disable_state_script_writing)
      graphics_info_t::disable_state_script_writing = 1;

   //
   if (cld.try_listener) { 
      std::cout << "INFO:: setting port and host "
		<< cld.port << " " << cld.hostname << std::endl;
      graphics_info_t::try_port_listener = 1;
      graphics_info_t::remote_control_port_number = cld.port;
      graphics_info_t::remote_control_hostname = cld.hostname;
   }

   // more settings for small screen
   if (cld.small_screen_display && graphics_info_t::use_graphics_interface_flag) {
     std::cout << "INFO:: setting only main icons for small screen" << std::endl;
     show_model_toolbar_main_icons();
     set_graphics_window_size(400, 400);
   }

}
