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
#include <filesystem>
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

void
logging::set_log_file(const std::string &file_name) {

   std::filesystem::path fp(file_name);
   std::filesystem::file_status status = std::filesystem::status(fp);
   output_file.open(file_name);

}


logging::ltw::ltw(const std::string &s_in) {
   type = type_t::STRING_TYPE;
   s = s_in;
}

logging::ltw::ltw(const char *s_in) {
   type = type_t::STRING_TYPE;
   if (s_in)
      s = std::string(s_in);
}

logging::ltw::ltw(bool b_in) {
   type = type_t::BOOL_TYPE;
   b = b_in;
}

logging::ltw::ltw(const int &i_in) {
   type = type_t::INT_TYPE;
   i = i_in;
}

logging::ltw::ltw(const unsigned int &i_in) {
   // slightly naughty
   type = type_t::INT_TYPE;
   i = i_in;
}

logging::ltw::ltw(const unsigned long &i_in) {
   // slightly naughty
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
   if (type == type_t::BOOL_TYPE)
      return std::to_string(b);
   return "--- unknown-type in ltw::to_string() ---";
}

void
logging::notify() {

   if (update_notifier_function)
      update_notifier_function();

   if (output_file) {
      if (! history.empty()) {
	 const auto &last_item = history.back();
	 std::string s = last_item.to_string();
	 output_file << s << "\n";
	 output_file.flush();
      }
   }
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
   notify();

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
   output_to_terminal_maybe();
   notify();
}

void
logging::log(log_t type, const std::string &s1, const std::string &s2) {

   // use:  add_log_item_to_history(log_time(type, s));
   log_item l(type);

   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   l.add_to_message(s1);
   l.add_to_message(" ");
   l.add_to_message(s2);
   history.push_back(l);
   output_to_terminal_maybe();
   notify();
}

void
logging::log(log_t type_in, const function_name_t &fn, const std::string &s1) {

   log_item l(type_in, fn, s1);
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   history.push_back(l);
   output_to_terminal_maybe();
   notify();
}

void
logging::log(log_t type_in, const function_name_t &fn, const std::string &s1, const std::string &s2) {

   log_item l(type_in, fn);
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   l.add_to_message(s1);
   l.add_to_message(" ");
   l.add_to_message(s2);
   history.push_back(l);
   output_to_terminal_maybe();
   notify();
}

void
logging::log(log_t type_in, const function_name_t &fn,
             const std::string &s1, const std::string &s2, const std::string &s3) {

   log_item l(type_in, fn);
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   l.add_to_message(s1);
   l.add_to_message(" ");
   l.add_to_message(s2);
   l.add_to_message(" ");
   l.add_to_message(s3);
   history.push_back(l);
   output_to_terminal_maybe();
   notify();
}

void
logging::log(log_t type_in, const function_name_t &fn, const std::string &s1, int i) {

   log_item l(type_in, fn);
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   l.add_to_message(s1);
   l.add_to_message(" ");
   l.add_to_message(std::to_string(i));
   history.push_back(l);
   output_to_terminal_maybe();
   notify();
}

void
logging::log(log_t type_in, const function_name_t &fn, const std::string &s1, int i, int j) {

   log_item l(type_in, fn);
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   l.add_to_message(s1);
   l.add_to_message(" ");
   l.add_to_message(std::to_string(i));
   l.add_to_message(std::to_string(j));
   history.push_back(l);
   output_to_terminal_maybe();
   notify();
}


void
logging::log(log_t type_in, const function_name_t &fn, const std::string &s1, int i, const std::string &s2) {

   log_item l(type_in, fn);
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   l.add_to_message(s1);
   l.add_to_message(" ");
   l.add_to_message(std::to_string(i));
   l.add_to_message(" ");
   l.add_to_message(s2);
   history.push_back(l);
   output_to_terminal_maybe();
   notify();
}


void
logging::log(log_t type_in, const function_name_t &fn, const std::vector<ltw> &v) {

   log_item l(type_in, fn);
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   for (const auto &item : v) {
      l.add_to_message(item.to_string());
      l.add_to_message(" ");
   }
   history.push_back(l);
   output_to_terminal_maybe();
   notify();
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
   output_to_terminal_maybe();
   notify();

}

void
logging::log(log_t type_in, const std::string &s1, std::size_t s, const std::string &s2) {

   log_item l(type_in);
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   l.message = s1;
   l.message += " ";
   l.message += std::to_string(s);
   l.message += " ";
   l.message += s2;
   history.push_back(l);
   output_to_terminal_maybe();
   notify();

}

void
logging::log(log_t type_in, const std::string &s1, int i, const std::string &s2) {

   log_item l(type_in);
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   l.message = s1;
   l.message += " ";
   l.message += std::to_string(i);
   l.message += " ";
   l.message += s2;
   history.push_back(l);
   output_to_terminal_maybe();
   notify();

}

void
logging::log(log_t type_in, const std::string &s1, unsigned int i, const std::string &s2) {

   log_item l(type_in);
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   l.message = s1;
   l.message += " ";
   l.message += std::to_string(i);
   l.message += " ";
   l.message += s2;
   history.push_back(l);
   output_to_terminal_maybe();
   notify();

}
void
logging::log(log_t type_in, const std::string &s1, double d, const std::string &s2) {

   log_item l(type_in);
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   l.message = s1;
   l.message += " ";
   l.message += std::to_string(d);
   l.message += " ";
   l.message += s2;
   history.push_back(l);
   output_to_terminal_maybe();
   notify();
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
   output_to_terminal_maybe();
   notify();
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
   output_to_terminal_maybe();
   notify();
}


void
logging::log(log_t type_in, const std::string &s1, const double &d1) {

   log_item l(type_in);
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   l.add_to_message(s1);
   l.add_to_message(" ");
   l.add_to_message(std::to_string(d1));
   history.push_back(l);
   output_to_terminal_maybe();
   notify();
}

void
logging::log(log_t type_in, const std::string &s1, const float &f1) {

   log_item l(type_in);
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   l.add_to_message(s1);
   l.add_to_message(" ");
   l.add_to_message(std::to_string(f1));
   history.push_back(l);
   output_to_terminal_maybe();
   notify();
}

void logging::log(log_t type_in, const std::string &s1, const float &f1, const float &f2) {

   log_item l(type_in);
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   l.add_to_message(s1);
   l.add_to_message(" ");
   l.add_to_message(std::to_string(f1));
   l.add_to_message(" ");
   l.add_to_message(std::to_string(f2));
   history.push_back(l);
   output_to_terminal_maybe();
   notify();
}

void logging::log(log_t type_in, ltw l1, ltw l2) {

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
   output_to_terminal_maybe();
   notify();
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
   output_to_terminal_maybe();
   notify();
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
   output_to_terminal_maybe();
   notify();

}

void
logging::log(log_t type_in, ltw l1, ltw l2, ltw l3, ltw l4, ltw l5) {

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
   l.message += " ";
   l.message += l5.to_string();
   history.push_back(l);
   output_to_terminal_maybe();
   notify();
}

void
logging::log(log_t type_in, ltw l1, ltw l2, ltw l3, ltw l4, ltw l5, ltw l6) {

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
   l.message += " ";
   l.message += l5.to_string();
   l.message += " ";
   l.message += l6.to_string();
   history.push_back(l);
   output_to_terminal_maybe();
   notify();
}

void
logging::log(log_t type_in,
             logging::ltw const &l1, logging::ltw const &l2,
             logging::ltw const &l3, logging::ltw const &l4,
             logging::ltw const &l5, logging::ltw const &l6,
             logging::ltw const &l7, logging::ltw const &l8,
             logging::ltw const &l9, logging::ltw const &l10) {

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
   l.message += " ";
   l.message += l5.to_string();
   l.message += " ";
   l.message += l6.to_string();
   l.message += " ";
   l.message += l7.to_string();
   l.message += " ";
   l.message += l8.to_string();
   l.message += " ";
   l.message += l9.to_string();
   l.message += " ";
   l.message += l10.to_string();
   history.push_back(l);
   output_to_terminal_maybe();
   notify();
}

void
logging::log(log_t type_in, const std::string &s1, const double &d1, const std::string &s2, const double &d2) {

   log_item l;
   l.type = type_in;
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   l.add_to_message(s1);
   l.add_to_message(" ");
   l.add_to_message(std::to_string(d1));
   l.add_to_message(" ");
   l.add_to_message(s2);
   l.add_to_message(" ");
   l.add_to_message(std::to_string(d2));
   l.add_to_message(" ");
   history.push_back(l);
   output_to_terminal_maybe();
   notify();

}

void
logging::log(log_t type_in, const function_name_t &fn, const std::string &s1, bool b, const std::string &s2, float f) {

   log_item l(type_in, fn);
   timeval current_time;
   int success = gettimeofday(&current_time, NULL);
   if (success == 0) // was successful
      l.t = current_time.tv_sec;
   l.add_to_message(s1);
   l.add_to_message(" ");
   l.add_to_message(std::to_string(b));
   l.add_to_message(" ");
   l.add_to_message(s2);
   l.add_to_message(" ");
   l.add_to_message(std::to_string(f));
   l.add_to_message(" ");
   history.push_back(l);
   output_to_terminal_maybe();
   notify();

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
   output_to_terminal_maybe();
   notify();
}

void
logging::output_to_terminal_maybe() {
   show_last();
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
      std::cout << h.to_string() << std::endl;
   }
}

std::string
logging::log_item::type_as_string() const {

   std::string tt;
   if (type == log_t::INFO)        tt = "INFO::      ";
   if (type == log_t::DEBUG)       tt = "DEBUG::  ";
   if (type == log_t::ERROR)       tt = "ERROR::  ";
   if (type == log_t::WARNING)     tt = "WARNING::";
   if (type == log_t::GL_ERROR)    tt = "GL_ERROR::";
   if (type == log_t::UNSPECIFIED) tt = "UNSPECIFIED::";
   return tt;
}

std::string
logging::log_item::to_string(bool include_datetime, bool use_markup) const {

   std::string ctime_str;
   if (include_datetime) {
      std::string ctime_str = ctime(&t);
      if (! ctime_str.empty()) {
	 ctime_str.pop_back();
	 ctime_str += ":";
      }
   }
   std::string tas = type_as_string();
   std::string o;
   if (ctime_str.empty())
      o = tas + " ";
   else
      o = tas + " " + ctime_str + " ";
   if (! function_name.empty())
      o += function_name.fn + "(): ";
   o += message;
   return o;
}

std::string
logging::log_item::get_date_string() const {

   std::string ctime_str = ctime(&t);
   if (! ctime_str.empty()) {
      ctime_str.pop_back();
      ctime_str += ":";
   }
   return ctime_str;
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
         std::cout << h.to_string() << std::endl;
      }
   }
}

std::vector<logging::log_item>
logging::get_log_history_from(unsigned int idx_start) const {

   std::vector<log_item> v;
   if (idx_start < history.size()) {
      for (unsigned int i=idx_start; i<history.size(); i++) {
	 v.push_back(history[i]);
      }
   }
   return v;
}
