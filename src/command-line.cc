/* src/command-line.cc
 * 
 * Copyright 2001, 2002, 2003, 2004, 2005, 2006 by The University of York
 * Copyright 2015, 2016 by Medical Research Council
 * Author: Paul Emsley
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
 * You should have received a copy of the GNU General Public License and
 * the GNU Lesser General Public License along with this program; if not,
 * write to the Free Software Foundation, Inc., 51 Franklin Street,  02110-1301, USA
 */


#if defined USE_PYTHON && !defined WINDOWS_MINGW
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

#ifdef __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#else
#define __GNU_LIBRARY__
#include "compat/coot-getopt.h"
#undef __GNU_LIBRARY__
#endif

#ifdef WINDOWS_MINGW
// for whatever reason (getopt) include python here for windows
#include "Python.h"
#endif //WINDOWS_MINGW

#include <iostream>
#include <string>
#include <vector>
#include <gtk/gtk.h>

// #include "mtz-bits.h" use cmtz-interface

#include "graphics-info.h"
// Including python needs to come after graphics-info.h, because
// something in Python.h (2.4 - chihiro) is redefining FF1 (in
// ssm_superpose.h) to be 0x00004000 (Grrr).
// BL says:: and (2.3 - dewinter), i.e. is a Mac - Python issue
// since the follwing two include python graphics-info.h is moved up
#include "c-interface.h"
#include "coot-version.hh"
#include "command-line.hh"


command_line_data
parse_command_line(int argc, char ** argv ) {

   command_line_data cld;

   coot_optind = 0; // reset the (global) extern because
	            // parse_command_line() is/can be called more once.

   int ch;

 /*
   * the option string will be passed to getopt(3), the format
   * of our string "sf:v:" will allow us to accept -s as a flag,
   * and -f or -v with an argument, the colon suffix tells getopt(3)
   * that we're expecting an argument.  Eg:  optest -s -f this -v8
   *
   */

   const char *optstr = "p:m:d:s:c:";

     /*
   * getopt(3) takes our argc, and argv, it also takes
   * the option string we set up earlier.  It will assign
   * the switch character to ch, and -1 when there are no more
   * command line options to parse.
   *
   */

   static struct option long_options[] = {
      {"pdb" ,   1, 0, 0},
      {"coords", 1, 0, 0},
      {"xyzin",  1, 0, 0},
      {"map" ,   1, 0, 0},
      {"data",   1, 0, 0},
      {"hklin",  1, 0, 0},
      {"auto",   1, 0, 0},
      {"script", 1, 0, 0},
      {"buster", 1, 0, 0},
      {"command", 1, 0, 0},
      {"ccp4-project", 1, 0, 0},
      {"dictionary", 1, 0, 0},
      {"code",       1, 0, 0},
      {"comp_id",    1, 0, 0},
      {"comp-id",    1, 0, 0},
      {"em",         0, 0, 0},
      {"emdb",       1, 0, 0},
      {"title",      1, 0, 0},
      {"port",       1, 0, 0},
      {"host",       1, 0, 0},
      {"hostname",   1, 0, 0}, // alternate for host
      {"help",       0, 0, 0},
      {"python",     0, 0, 0},
      {"run-state-script",   0, 0, 0},
      {"splash-screen",      1, 0, 0}, // alternate splash screen
      {"self-test",          0, 0, 0},
      {"show-ccp4i2-save-button", 0, 0, 0},
      {"opengl-es",          0, 0, 0},
      {"no-state-script",    0, 0, 0},
      {"no-startup-scripts", 0, 0, 0},
      {"no-graphics",      0, 0, 0},
      {"no-splash-screen", 0, 0, 0},
      {"stereo",           0, 0, 0},        // no arguments
      {"side-by-side",     0, 0, 0},  //
      {"zalman-stereo",    0, 0, 0},  //
      {"version",          0, 0, 0},        //
      {"version-full",     0, 0, 0},     //
      {"no-guano",         0, 0, 0},        //
      {"small-screen", 0, 0, 0},      // no arguments (setting for small screens)
      {"update-self", 0, 0, 0},
      {0, 0, 0, 0}	       // must have blanks at end
   };

   int option_index = 0;

   bool found_no_graphics_in_the_command_line = false;

   while( -1 !=
	  (ch = coot_getopt_long(argc, argv, optstr, long_options, &option_index) )) {

      switch(ch) {

      case 0:

	 if (coot_optarg) {

	    // options that need an argument:

	    std::string arg_str = long_options[option_index].name;

	    if (arg_str == "pdb") {
	       cld.coords.push_back(coot_optarg);
	    }
	    if (arg_str == "coords") {
	       cld.coords.push_back(coot_optarg);
	    }
	    if (arg_str == "xyzin") {
	       cld.coords.push_back(coot_optarg);
	    }
	    if (arg_str == "map") {
	       cld.maps.push_back(coot_optarg);
	    }
	    if (arg_str == "data") {
	       cld.datasets.push_back(coot_optarg);
	    }
	    if (arg_str == "hklin") {
	       cld.datasets.push_back(coot_optarg);
	    }
	    if (arg_str == "script") {
	       cld.script.push_back(coot_optarg);
	    }
	    if (arg_str == "command") {
	       cld.command.push_back(coot_optarg);
	    }
	    if (arg_str == "port") {
	       cld.port = atoi(coot_optarg);
	    }
	    if (arg_str == "host") {
	       cld.hostname = coot_optarg;
	    }
	    if (arg_str == "hostname") {
	       cld.hostname = coot_optarg;
	    }
	    if (arg_str == "auto") {
	       cld.auto_datasets.push_back(coot_optarg);
	    }
	    if (arg_str == "dictionary") {
	       cld.dictionaries.push_back(coot_optarg);
	    }
	    if (arg_str == "dictionary-with-mol") {
	       cld.dictionaries_with_mol.push_back(coot_optarg);
	    }
	    if (arg_str == "ccp4-project") {
	       cld.ccp4_project = coot_optarg;
	    }
	    if (arg_str == "code") {
	       cld.accession_codes.push_back(coot_optarg);
	    }
	    if (arg_str == "emdb") {
	       cld.emdb_codes.push_back(coot_optarg);
	    }
	    if (arg_str == "comp_id") {
	       cld.comp_ids.push_back(coot_optarg);
	    }
	    if (arg_str == "comp-id") {
	       cld.comp_ids.push_back(coot_optarg);
	    }
	    if (arg_str == "show-ccp4i2-save-button") {
	       cld.show_ccp4i2_save_button = true;
	    }
	    if (arg_str == "buster") {
	       cld.open_buster_output_files = true;
	    }
	    if (arg_str == "title") {
	       cld.title = coot_optarg;
	    }
	    if (arg_str == "splash-screen") {
	       cld.alternate_splash_screen_file_name = coot_optarg;
	    }

	 } else {

	    // long argument without parameter:
	    std::string arg_str(long_options[option_index].name);

	    if (arg_str == "stereo") {
	       cld.hardware_stereo_flag = 1;
	    } else {
	       if (arg_str == "zalman-stereo") {
		  cld.hardware_stereo_flag = 5;
	       } else {

		  // Thanks for suggesting this Ezra.
		  //
		  if (arg_str == "help") {
		     std::cout << std::endl
			       << "Usage: coot [--pdb pdb-file-name]\n"
			       << "            [--coords pdb/cif/shelx-filename]\n"
			       << "            [--map ccp4-map-file-name]\n"
			       << "            [--data mtz-file-name]\n"
			       << "            [--hklin mtz-file-name]\n"
			       << "            [--auto mtz-file-name]\n"
			       << "            [--dictionary cif-dictionary-file-name]\n"
			       << "            [--dictionary-with-mol cif-dictionary-file-name]\n"
			       << "            [--script script-file-name]\n"
			       << "            [--em]\n"
			       << "            [--title some-title]\n"
			       << "            [--command command-script]\n"
			       << "            [--small-screen]\n"
			       << "            [--splash-screen]\n"
			       << "            [--stereo]\n"
                        //			       << "            [--zalman-stereo]\n"
                        //             << "            [--side-by-side]\n"
			       << "            [--version]\n"
			       << "            [--show-ccp4i2-save-button]\n"
			       << "            [--self-test]\n"
			       << "            [--no-state-script]\n"
			       << "            [--no-startup-scripts]\n"
			       << "            [--no-splash-screen]\n"
			       << "            [--opengl-es]\n"
			       << "            [--no-graphics]\n"
			       << "            [--no-guano]\n"
			       << std::endl;
		     coot_no_state_real_exit(0); // merge conflict resolved 4c1ace414
		  } else {

			if (arg_str == "version") {
			   std::cout  << VERSION << " " << coot_version_extra_info();
			   // this is in coot_version_extra_info() now
			   // std::cout << "Binary type: " << COOT_SYS_BUILD_TYPE << std::endl;
			   std::vector<std::string> enableds;
#ifdef MAKE_ENHANCED_LIGAND_TOOLS
			   enableds.push_back("Enhanced-ligand-tools");
#endif
#ifdef HAVE_CXX11
			   enableds.push_back("C++-11");
#endif
#ifdef HAVE_CXX_THREAD
			   enableds.push_back("Threads");
#endif
#ifdef HAVE_BOOST_BASED_THREAD_POOL_LIBRARY
			   enableds.push_back("Boost-based-thread-pool");
#endif
#ifdef USE_MOLECULES_TO_TRIANGLES
			   enableds.push_back("MoleculesToTriangles");
#endif
#ifdef HAVE_GOOCANVAS
			   enableds.push_back("Goocanvas");
#endif
#ifdef USE_GEMMI
            enableds.push_back("GEMMI");
#endif
#ifdef USE_SQLITE3
			   enableds.push_back("SQLite3");
#endif
#ifdef HAVE_CCP4SRS
			   enableds.push_back("CCP4SRS");
#endif
#ifdef USE_LIBCURL
			   enableds.push_back("LibCurl");
#endif
			   if (enableds.size()) {
			      std::cout << "Enabled: ";
			      for (unsigned int i=0; i<enableds.size(); i++)
				 std::cout << enableds[i] << " ";
			      std::cout << std::endl;
			   }
			   std::string s = COOT_BUILD_INFO_STRING;
			   if (s.length() > 0)
			      std::cout << "Builder_info: " << s << std::endl;
			   exit(0);
			} else {
			   if (arg_str == "python") {
			      cld.script_is_python_flag = 1;
			   } else {
			      if (arg_str == "no-state-script") {
				 graphics_info_t::run_state_file_status = 0;
			      } else {
                                 if (arg_str == "run-state-script") {
                                    graphics_info_t::run_state_file_status = 2;
                                 }  else {
                                    if (arg_str == "no-startup-scripts") {
                                       graphics_info_t::run_startup_scripts_flag = false;
                                    } else {
                                       if (arg_str == "no-graphics") {
                                          cld.do_graphics = false; // 20220409-PE it's already false now
                                          found_no_graphics_in_the_command_line = true;
                                       } else {
                                          if (arg_str == "em") {
                                             cld.em_mode = true;
                                          } else {
                                             if (arg_str == "side-by-side") {
                                                cld.hardware_stereo_flag = 2;
                                             } else {
                                                if (arg_str == "no-guano") {
                                                   cld.disable_state_script_writing = 1;
                                                } else {
                                                   if (arg_str == "small-screen") {
                                                      cld.small_screen_display = 1;
                                                   } else {
                                                      if (arg_str == "no-splash-screen") {
                                                         cld.use_splash_screen = 0;
                                                      } else {
                                                         if (arg_str == "self-test") {
                                                            cld.run_internal_tests_and_exit = 1;
                                                         } else {
                                                            if (arg_str == "opengl-es") {
                                                               cld.use_opengl_es = true;
                                                            } else {
                                                               if (arg_str == "show-ccp4i2-save-button") {
                                                                  cld.show_ccp4i2_save_button = true;
                                                               } else {
                                                                  if (arg_str == "update-self") {
                                                                     cld.update_self = 1;
                                                                     cld.do_graphics = 0;
                                                                  } else {
                                                                     std::cout << "WARNING! Malformed option - needs an argument: "
                                                                               << long_options[option_index].name
                                                                               << std::endl << std::endl;
                                                                  }
                                                               }
                                                            }
                                                         }
                                                      }
                                                   }
                                                }
                                             }
                                          }
                                       }
                                    }
                                 }
                              }
			   }
			}

		  }
	       }
	    }
	 }
	 break;

	 // try short options then...

      case 'p':
	 cld.coords.push_back(coot_optarg);
	 break;

      case 's':
	 cld.script.push_back(coot_optarg);
	 break;

      case 'd':
	 cld.datasets.push_back(coot_optarg);
	 break;

      case 'a':
	 cld.auto_datasets.push_back(coot_optarg);
	 break;

      case 'm':
	 cld.maps.push_back(coot_optarg);
	 break;

      case 'c':
         if (coot_optarg) {
            // std::cout << "command coot_optarg: " << coot_optarg << std::endl;
	    cld.command.push_back(coot_optarg);
         } else {
            std::cout << "command coot_optarg is NULL " << std::endl;
         }
	 break;

      case '?':
	 std::cout << "Unrecognised option: " << optopt << std::endl;
	 break;

      default:
	 std::cout << "Unaccounted for coot_optarg condition " << std::endl;
	 break;
      }
   }

   if (cld.hostname != "" && cld.port != 0)
      cld.try_listener = true;

   if (! found_no_graphics_in_the_command_line) // we do this for "import coot" in Python, which doesn't use main()
      cld.do_graphics = true;                   // and doesn't run this function.

   cld.roberto_pdbs(argc, argv);

   return cld;

}

void
command_line_data::handle_immediate_settings() {

   // start the graphics?
   if (do_graphics) {
      graphics_info_t::use_graphics_interface_flag = true;
   } else {
      graphics_info_t::use_graphics_interface_flag = false;
   }

   if (script_is_python_flag)
      graphics_info_t::python_at_prompt_flag = 1;

   if (update_self)
      graphics_info_t::update_self = 1;

   // small screen
   if (small_screen_display && graphics_info_t::use_graphics_interface_flag) {
      std::cout <<"INFO:: set labels and icons for small screens" <<std::endl;

      std::cout << "Fix small screen parsing in handle_immediate_settings() " << std::endl;
      // gtk_rc_parse_string("gtk-icon-sizes=\"gtk-large-toolbar=10,10:gtk-button=10,10\"");
      // gtk_rc_parse_string("class \"GtkLabel\" style \"small-font\"");
      graphics_info_t::graphics_x_size = 400;
      graphics_info_t::graphics_y_size = 400;
   }
}


// add any pdb files not alread added with --pdb/coords/xyzin
void
command_line_data::roberto_pdbs(int argc, char **argv) {

   for (int i=1; i<argc; i++) {
      std::string file(argv[i]);
      if (coot::util::extension_is_for_coords(coot::util::file_name_extension(file)) ||
          coot::util::extension_is_for_shelx_coords(coot::util::file_name_extension(file)))
         if (std::find(coords.begin(), coords.end(), file) == coords.end())
            coords.push_back(file);
   }
}


void
command_line_data::add(const std::string &file_name) {

   // add to command_line_scripts or command_line_commands
   // so that it behaves like it did before GtkApplication usage
   //

   auto is_dictionary = [] (const std::string &file_name) {
      bool state = false;
      if (coot::file_exists(file_name)) {
         std::ifstream f(file_name);
         if (f) {
            unsigned int line_count = 0;
            std::string line;
            while (std::getline(f, line)) {
               if (line.find("_chem_comp_atom.x") != std::string::npos) {
                  state = true;
                  break;
               }
               line_count++;
               if (line_count> 1000) break;
            };
         }
      }
      return state;
   };

   std::string extension = coot::util::file_name_extension(file_name);
   if (extension == ".pdb")
      coords.push_back(file_name);
   if (extension == ".ent")
      coords.push_back(file_name);
   if (extension == ".mmcif") // mmcifs are always coordinates, right?
      coords.push_back(file_name);
   if (extension == ".map")
      maps.push_back(file_name);
   if (extension == ".mrc")
      maps.push_back(file_name);
   if (extension == ".mtz")
      auto_datasets.push_back(file_name);
   if (extension == ".py")
      script.push_back(file_name);
   if (extension == ".scm")
      script.push_back(file_name);

   if (false) // debugging
      if (extension == ".cif")
         std::cout << ":::::::::::::::::::::::::: is_dictionary " << file_name
                   << " " << is_dictionary(file_name) << std::endl;

   if (extension == ".cif") {
      if (is_dictionary(file_name))
         dictionaries.push_back(file_name);
      else
         coords.push_back(file_name);
   }

}

