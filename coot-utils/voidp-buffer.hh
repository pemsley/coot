/* coot-utils/voidp-buffer.hh
 * 
 * Copyright 2023 by Global Phasing Ltd.
 * Author: Jakub Smulski
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

#ifndef VOIDP_BUFFER_HH
#define VOIDP_BUFFER_HH

#include <algorithm>

#include <zlib.h>
#include <stdlib.h>
#include <memory.h>

class voidp_buffer_t {
   voidp buffer;
   size_t buffer_size;

public:

   explicit voidp_buffer_t(size_t inital_size) : buffer(calloc(inital_size, sizeof(char))) {
      buffer_size = inital_size;
   }
   ~voidp_buffer_t() {
      free(buffer);
   }
   inline size_t size() const noexcept {
      return buffer_size;
   }
   void resize(size_t new_size) {
      voidp new_buffer = calloc(new_size, sizeof(char));
      memcpy(new_buffer, buffer, std::min(buffer_size, new_size) * sizeof(char));
      free(buffer);
      buffer = new_buffer;
      buffer_size = new_size;
   }
   inline voidp get() {
      return buffer;
   }
};


#endif // VOIDP_BUFFER_HH
