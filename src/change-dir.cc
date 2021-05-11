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
	    const char *s = getenv("HOME");
	    if (s) {
	       std::cout << "INFO:: changing working directory to " << s << std::endl;
	       int status = chdir(s);
               if (status != 0)
                  std::cout << "ERROR:: failed to change working directory to " << s << std::endl;
	    } else {
	       s = getenv("COOT_HOME");
	       if (s) {
		  std::cout << "INFO:: changing working directory to " << s << std::endl;
		  int status = chdir(s);
                  if (status != 0)
                     std::cout << "ERROR:: failed to change working directory to " << s << std::endl;
	       }
	    }
	 } 
      } 
   } else {
      // can't stat .
   } 

} 
