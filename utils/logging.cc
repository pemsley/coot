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
gettimeofday (struct timeval *tv, void *tz) {

   union {
      __int64 ns100;              /*time since 1 Jan 1601 in 100ns units */
      FILETIME ft;
  } now;

  GetSystemTimeAsFileTime(&now.ft);
  tv->tv_usec = (long) ((now.ns100 / 10LL) % 1000000LL);
  tv->tv_sec = (long) ((now.ns100 - EPOCHFILETIME) / 10000000LL);
  return (0);
}
#endif /* _MSC_VER */

logging::ltw::ltw(const std::string &s_in) {
   type = type_t::STRING_TYPE;
   s = s_in;
}

logging::ltw::ltw(bool b_in) {
   type = type_t::BOOL_TYPE;
   b = b_in;
}

logging::ltw::ltw(const int &i_in) {
   type = type_t::INT_TYPE;
   i = i_in;
}

logging::ltw::ltw(float f_in) {
   type = type_t::FLOAT_TYPE;
   f = f_in;
}

logging::ltw::ltw(const double &d_in) {
   type = type_t::DOUBLE_TYPE;
   d = d_in;
}

std::string
logging::ltw::to_string() const {
   if (type == type_t::STRING_TYPE)
      return s;
   if (type == type_t::INT_TYPE)
      return std::to_string(i);
   if (type == type_t::FLOAT_TYPE)
      return std::to_string(f);
   if (type == type_t::DOUBLE_TYPE)
      return std::to_string(d);
   return "---unknown-type-in-to-string---";
}

void
logging::log(log_t type, const std::string &s) {

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
   log_item l(log_t::WARNING, s);
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;

   history.push_back(l);

}

void
logging::log(log_t type_in, const function_name_t &fn, const std::string &s1) {
   log_item l(type_in, fn, s1);
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   history.push_back(l);
}


void
logging::log(log_t type_in, const std::string &s1, bool v1, const std::string &s2) {

   log_item l(type_in);
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   std::string b = "False";
   if (v1) b = "True";
   l.message = s1;
   l.message += " ";
   l.message += b;
   l.message += " ";
   l.message += s2;
   history.push_back(l);

}

void
logging::log(log_t type_in, const std::string &s1, const int &i) {

   log_item l(type_in);
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   l.message = s1;
   l.message += " ";
   l.message += std::to_string(i);
   history.push_back(l);
}


void
logging::log(log_t type_in, const std::string &s1, bool v1, const std::string &s2, const std::string &s3) {

   log_item l(type_in);
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   std::string b = "False";
   if (v1) b = "True";
   l.message = s1;
   l.message += " ";
   l.message += b;
   l.message += " ";
   l.message += s2;
   l.message += " ";
   l.message += s3;
   history.push_back(l);
}

void
logging::log(log_t type_in, ltw l1, ltw l2) {

   log_item l;
   l.type = type_in;
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   l.message += l1.to_string();
   l.message += " ";
   l.message += l2.to_string();
   history.push_back(l);

}

void
logging::log(log_t type_in, ltw l1, ltw l2, ltw l3) {

   log_item l;
   l.type = type_in;
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   l.message += l1.to_string();
   l.message += " ";
   l.message += l2.to_string();
   l.message += " ";
   l.message += l3.to_string();
   history.push_back(l);

}

void
logging::log(log_t type_in, const ltw &l1, const ltw &l2, const ltw &l3, const ltw &l4) {

   log_item l;
   l.type = type_in;
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   l.message += l1.to_string();
   l.message += " ";
   l.message += l2.to_string();
   l.message += " ";
   l.message += l3.to_string();
   l.message += " ";
   l.message += l4.to_string();
   history.push_back(l);

}

void
logging::log(log_t type_in, const std::vector<ltw> &ls) {

   log_item l;
   l.type = type_in;
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   for (const auto &item : ls) {
      l.message += item.to_string();
      l.message += " ";
   }
   history.push_back(l);
}

void
logging::show() const {

   for (std::size_t i=0; i<history.size(); i++) {
      const log_item &h = history[i];
      // maybe use localtime()
      std::string ctime_str = ctime(&h.t);
      if (! ctime_str.empty()) {
         ctime_str.pop_back();
         ctime_str += ":";
      }
      std::cout << ctime_str << " " << h.message << std::endl;
   }
}

void
logging::show_last() const {

   if (! history.empty()) {
      const log_item &h = history.back();

      bool do_it = false;
      if (output_type == output_t::TERMINAL)
         if (h.type != log_t::DEBUG)
            do_it = true;
      if (output_type == output_t::TERMINAL_WITH_DEBUGGING)
         do_it = true;

      if (do_it) {
         // make a show function for a log_item!
         std::string type_as_string;
         if (h.type == log_t::INFO)        type_as_string = "INFO::";
         if (h.type == log_t::DEBUG)       type_as_string = "DEBUG::";
         if (h.type == log_t::ERROR)       type_as_string = "ERROR::";
         if (h.type == log_t::WARNING)     type_as_string = "WARNING::";
         if (h.type == log_t::UNSPECIFIED) type_as_string = "UNSPECIFIED::";
         std::string ctime_str = ctime(&h.t);
         if (! ctime_str.empty()) {
            ctime_str.pop_back();
            ctime_str += ":";
         }
         std::cout << type_as_string << " " << ctime_str << " " << h.message << std::endl;
      }
   }
}
