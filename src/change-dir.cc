/* src/change-dir.cc
 * 
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#include <sys/types.h>
#include <sys/stat.h>

#ifndef _MSC_VER 
#include <unistd.h>
#else
#define S_IRUSR S_IREAD
#define S_IWUSR S_IWRITE
#define S_IXUSR S_IEXEC
#include <windows.h>
#include <direct.h>
#endif // _MSC_VER

#include <stdlib.h>

#include <iostream>
#include "change-dir.hh"
#include "utils/coot-utils.hh"

// Mac users often start somewhere where thy can't write files
// 
void change_directory_maybe() {

   struct stat buf;

   int i = stat(".", &buf);
   if (i == 0) {
      if (! S_ISDIR(buf.st_mode) ) {
	 std::cout << "INFO:: in change_directory_maybe() strange " << std::endl;
      } else {
	 if (buf.st_mode & S_IWUSR) {
	    // OK
	 } else {
	    // No write permission.
            std::string s = coot::get_home_dir();
            if (!s.empty()) {
	       std::cout << "INFO:: changing working directory to " << s << std::endl;
	       int state = chdir(s.c_str());
               if (state != 0)
                  std::cout << "Faked to change dir to " << s << std::endl;
	    } else {
	       const char *cs = getenv("COOT_HOME");
	       if (cs) {
		  std::cout << "INFO:: changing working directory to " << s << std::endl;
		  int state = chdir(cs);
                  if (state != 0)
                     std::cout << "Faked to change dir to " << s << std::endl;
	       }
	    }
         }
      }
   } else {
      // can't stat .
   }

}
