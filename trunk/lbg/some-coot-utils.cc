
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>

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
