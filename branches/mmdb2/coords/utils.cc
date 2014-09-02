/* coords/utils.h
 * 
 * Copyright 2004, 2005 by The University of York
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

