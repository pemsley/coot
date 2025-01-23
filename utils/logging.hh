/*
 * utils/logging.hh
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


#ifndef LOGGING_HH
#define LOGGING_HH

#include <ctime>
#include <string>
#include <vector>

#ifdef _MSC_VER
#undef WARNING
#undef ERROR
#endif

class logging {

public:
   enum type_t { INFO, WARNING, DEBUG, ERROR, UNSPECIFIED};

   class function_name_t {
   public:
      explicit function_name_t() { fn = "unspecified"; }
      explicit function_name_t(const std::string &s) { fn = s; }
      std::string fn;
   };

   class log_item {
   public:
      time_t t;
      type_t type;
      function_name_t function_name;
      std::string message;
      log_item(type_t type_in, const std::string &message_in) : t(0), type(type_in), message(message_in) {}
      log_item(type_t type_in, const function_name_t &fn, const std::string &message_in) :
         t(0), type(type_in), function_name(fn), message(message_in) {}
      log_item(const std::string &message_in) : t(0), type(UNSPECIFIED), message(message_in) {}
      friend std::ostream& operator<<(std::ostream &o, const log_item &li);
   };

   class ltw { // logging type wrapper
   public:
      ltw(const std::string &s);
      ltw(bool &b);
      ltw(const int &i);
      ltw(float f);
      ltw(const double &d);
   };

private:
   std::vector<log_item> history;
   void operator<<(const std::string &s);
 public:
   logging() {}
   void log(const std::string &s);
   void log(type_t type_in, const std::string &s1, bool v1, const std::string &s2);
   void log(type_t type_in, const std::string &s1, bool v1, const std::string &s2, const std::string &s3);
   void log(type_t type_in, const std::string &s1);
   void log(type_t type_in, const std::string &s1, const std::string &s2);
   void log(type_t type_in, const std::string &s1, const int &i);
   void log(type_t type_in, const std::string &s1, const float &f1);
   void log(type_t type_in, const std::string &s1, const double &d1);
   void log(type_t type_in, const std::string &s1, const float &f1, const float &f2);
   void log(type_t type_in, const std::string &s1, const double &d1, const std::string &s2, const double &d2);
   void log(type_t type_in, const std::string &s1, const std::string &s2, const std::string &s3, const float &f1, const std::string &s4, const float &f2);
   void log(type_t type_in, const std::string &s1, const float &f1, const std::string &s2, const std::string &s3, const std::string &s4, const float &f2);
   void log(type_t type_in, const ltw &l1, const ltw &l2, const ltw &l3, const ltw &l4);
   void log(type_t type_in, ltw l1, ltw l2);
   void log(type_t type_in, ltw l1, ltw l2, ltw l3, ltw l4, ltw l5);
   void log(type_t type_in,
            const ltw &l1, const ltw &l2, const ltw &l3, const ltw &l4, const ltw &l5,
            const ltw &l6, const ltw &l7, const ltw &l8, const ltw &l9, const ltw &l10);
   void log(type_t type_in,
            const ltw &l1, const ltw &l2, const ltw &l3, const ltw &l4, const ltw &l5,
            const ltw &l6, const ltw &l7, const ltw &l8, const ltw &l9, const ltw &l10, const ltw &l11);

   // 20241211-PE now I want to pass the function name also - let's do those...
   void log(type_t type_in, const function_name_t &fn, const std::string &s1);

   void show() const;
   friend std::ostream& operator<<(std::ostream &o, const log_item &li);
};

// logging logging;

#endif // LOGGING_HH
