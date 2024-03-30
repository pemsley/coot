/*
 * utils/logging.cc
 *
 * Copyright 2009 by University of Oxford
 * Author: Paul Emsley
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 3 of the License, or (at
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

#ifndef _MSC_VER
#include <sys/time.h>
#else
#include <windows.h>
#endif

#include <iostream>
#include "logging.hh"

#if defined(_MSC_VER)
#define EPOCHFILETIME (116444736000000000i64)

int
gettimeofday (struct timeval *tv, void *tz)
{
  union
  {
    __int64 ns100;              /*time since 1 Jan 1601 in 100ns units */
    FILETIME ft;
  } now;

  GetSystemTimeAsFileTime (&now.ft);
  tv->tv_usec = (long) ((now.ns100 / 10LL) % 1000000LL);
  tv->tv_sec = (long) ((now.ns100 - EPOCHFILETIME) / 10000000LL);
  return (0);
}
#endif /* _MSC_VER */

void
logging::log(logging::type_t type, const std::string &s) {

   // use:  add_log_item_to_history(log_time(type, s));
   log_item l(type, s);

   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   
   history.push_back(l);

}

void
logging::log(const std::string &s) {

   // extract the type from s:
   log_item l(WARNING, s);
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   
   history.push_back(l);

}


void
logging::show() const {

   for (std::size_t i=0; i<history.size(); i++) {
      const log_item &h = history[i];
      // maybe use localtime()
      std::cout << ctime(&h.t) << " " << h.type << " " << h.message << std::endl;
   }

}

