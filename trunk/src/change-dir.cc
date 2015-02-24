

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
	       chdir(s);
	    } else {
	       s = getenv("COOT_HOME");
	       if (s) { 
		  std::cout << "INFO:: changing working directory to " << s << std::endl;
		  chdir(s);
	       }
	    }
	 } 
      } 
   } else {
      // can't stat .
   } 

} 
