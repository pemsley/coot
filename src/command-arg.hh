/*
 * src/command-arg.hh
 *
 * Copyright 2011 by University of York
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


#ifndef COMMAND_ARG_HH
#define COMMAND_ARG_HH

#include <clipper/core/clipper_types.h>

#include "utils/coot-utils.hh"

namespace coot {
   class command_arg_t {
   public:
      enum coot_script_arg_type{UNSET, INT, FLOAT, STRING, BOOL};
      command_arg_t(int iin) {
         i = iin;
         type = INT;
         b = false;
         f = -1;
      }
      command_arg_t(float fin) {
         f = fin;
         type = FLOAT;
         b = false;
         i = -1;
      }
      command_arg_t(const clipper::String &sin) : s(sin) {
         type = STRING;
         b = false;
         f = -1;
         i = -1;
      }
      command_arg_t(const std::string &sin) : s(sin) {
         type = STRING;
         b = false;
         f = -1;
         i = -1;
      }
      command_arg_t(const char *sin) : s(sin) {
         type = STRING;
         b = false;
         f = -1;
         i = -1;
      }
      command_arg_t(bool bin) {
         b = bin;
         type = BOOL;
         f = -1;
         i = -1;
      }
      command_arg_t() {
         type = UNSET;
         b = false;
         f = -1;
         i = -1;
      }
      coot_script_arg_type type;
      bool b;
      float f;
      int i;
      std::string s;
      std::string as_string() const {
         std::string os("unknown-arg-type");
         if (type == INT)
            os = coot::util::int_to_string(i);
         if (type == FLOAT)
            os = coot::util::float_to_string(f);
         if (type == STRING)
            os = std::string("'") + s + std::string("'");
         if (type == BOOL) {
            if (b)
               os = "True";
            else
               os = "False";
         }
         return os;
      }
   };
}

#endif // COMMAND_ARG_HH
