
#include <sys/types.h>  // stating
#include <sys/stat.h>

#if !defined _MSC_VER
#include <unistd.h>
#endif
#include <stdlib.h>  // for malloc
#include <stdio.h>  // for printf

#include "utils.h"

int 
does_file_exist(const char *file_name) { 

   int iout;
   int i;

   struct stat buf;

   i = stat(file_name, &buf);

   if (i == 0) { 
      printf("%s is statable\n", file_name); 
      iout = 1;
   } else { 
      printf("%s is not statable\n", file_name); 
      iout = 0;
   } 
   return iout;
} 

