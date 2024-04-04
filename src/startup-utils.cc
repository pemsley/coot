/*
 * src/startup-utils.cc
 *
 * Copyright 2022 by Medical Research Council
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
#include "startup-utils.hh"
extern "C" int git_revision_count(); // 20220612-PE doesn this need to be extern now?

std::string
make_main_window_title() {

   std::string version_string = VERSION;
   std::string main_title = "Coot " + version_string;

#ifdef MAKE_ENHANCED_LIGAND_TOOLS
   // main_title += " EL";
#endif

#ifdef COOT_MAIN_TITLE_EXTRA
   main_title += COOT_MAIN_TITLE_EXTRA;
#else

   // if this is a pre-release, stick in the revision number too
   if (version_string.find("-pre") != std::string::npos) {
      main_title += " (revision count ";
      main_title += coot::util::int_to_string(git_revision_count());
      main_title += ")";
   }
#endif

#ifdef WINDOWS_MINGW
   main_title = "Win" + main_title;
#endif

   return main_title;
}
