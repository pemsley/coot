/* lbg/some-coot-utils.cc
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

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>

#ifdef MAKE_ENTERPRISE_TOOLS
// no stand-alone functions need be added
#else 

#include "some-coot-utils.hh"

bool coot::is_directory_p(const std::string &filename) {

   bool st = 0;
   struct stat s; 
   int fstat = stat(filename.c_str(), &s);
   if ( fstat == -1 ) { // file not exist
      return 0;
   } else {
      if (S_ISDIR(s.st_mode)) {
	 return 1;
      } else {
	 return 0;
      }
   }
   return st;
}



std::string
coot::util::int_to_string(int i) {
   char s[100];
   snprintf(s,99,"%d",i);
   return std::string(s);
}


std::string
coot::util::float_to_string(float f) {
   char s[100];
   snprintf(s,99,"%5.2f",f);
   return std::string(s);
}

std::vector<std::string>
coot::util::split_string_no_blanks(const std::string &string_in,
				   const std::string &splitter) {

  
   std::vector<std::string> v;
   std::string s=string_in;

   while (1) {
      std::string::size_type isplit=s.find_first_of(splitter);
      if (isplit != std::string::npos) {
         std::string f = s.substr(0, isplit);
	 if (f.length() > 0)
	    v.push_back(f);
         if (s.length() >= (isplit+splitter.length())) { 
            s = s.substr(isplit+splitter.length());
         } else {
            break;
         }
      } else {
         v.push_back(s);
         break;
      } 
   }
   return v;
}



std::string
coot::util::downcase(const std::string &s) {

   std::string r = s;
   std::string::iterator it=r.begin();

   while ( (it!=r.end()) ) { 
      *it = ::tolower(*it);
      it++;
   }
   return r;
}

#endif // MAKE_ENTERPRISE_TOOLS
