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
#include <fstream>

#ifdef _MSC_VER
#undef WARNING
#undef ERROR
#endif

namespace coot {
   class residue_spec_t;
   class atom_spec_t;
}

enum class log_t { INFO, WARNING, DEBUG, ERROR, GL_ERROR, UNSPECIFIED};

class logging {

public:

   class function_name_t {
   public:
      explicit function_name_t() { fn = "unspecified"; }
      explicit function_name_t(const std::string &s) { fn = s; }
      std::string fn;
      std::string get_name() const { return fn + "()"; }
      bool empty() const { return (fn == std::string("unspecified")); }
   };

   class ltw { // logging type wrapper
      enum class type_t { STRING_TYPE, BOOL_TYPE, INT_TYPE, FLOAT_TYPE, DOUBLE_TYPE};
      type_t type;
      std::string s;
      bool b;
      int i;
      double d;
      float f;
   public:
      ltw(const std::string &s);
      ltw(const char *s);
      ltw(bool b);
      ltw(const int &i);
      ltw(const unsigned int &i);
      ltw(const unsigned long &i);
      ltw(float f);
      ltw(const double &d);
      std::string to_string() const;
   };

   class log_item {
   public:
      time_t t;
      log_t type;
      function_name_t function_name;
      std::string message;
      log_item() : t(0), type(log_t::UNSPECIFIED) {}
      log_item(log_t type_in) : type(type_in) {}
      log_item(log_t type_in, const std::string &message_in) : t(0), type(type_in), message(message_in) {}
      log_item(log_t type_in, const function_name_t &fn) : t(0), type(type_in), function_name(fn) {}
      log_item(log_t type_in, const function_name_t &fn, const std::string &message_in) :
         t(0), type(type_in), function_name(fn), message(message_in) {}
      log_item(const std::string &message_in) : t(0), type(log_t::UNSPECIFIED), message(message_in) {}
      log_item(const std::vector<ltw> &ls);
      void add_to_message(const std::string &s) { message += s; }
      std::string to_string(bool include_datetime=false, bool use_markup=true) const;
      std::string type_as_string() const;
      std::string get_date_string() const;
      friend std::ostream& operator<<(std::ostream &o, const log_item &li);
   };

private:
   std::vector<log_item> history;
   void operator<<(const std::string &s);
   void output_to_terminal_maybe();
   void (*update_notifier_function)();
   void notify();
   std::ofstream output_file;

public:
   enum class output_t {TERMINAL, TERMINAL_WITH_DEBUGGING, INTERNAL, BOTH};
   logging() : update_notifier_function(nullptr), output_type(output_t::TERMINAL) {}
   output_t output_type;
   //! if output_type is TERMINAL, show_last() writes to the terminal. Else
   //! show_last() does nothing.
   void set_output_type(output_t ot) {output_type = ot; }
   void set_log_file(const std::string &file_name);
   void log(const std::string &s);
   void log(log_t type_in, const std::string &s1, bool v1, const std::string &s2);
   void log(log_t type_in, const std::string &s1, bool v1, const std::string &s2, const std::string &s3);
   void log(log_t type_in, const std::string &s1, int, const std::string &s2);
   void log(log_t type_in, const std::string &s1, unsigned int, const std::string &s2);
   void log(log_t type_in, const std::string &s1, std::size_t s, const std::string &s2);
   void log(log_t type_in, const std::string &s1, double d, const std::string &s2);
   void log(log_t type_in, const std::string &s1);
   void log(log_t type_in, const std::string &s1, const std::string &s2);
   void log(log_t type_in, const std::string &s1, const int &i);
   void log(log_t type_in, const std::string &s1, const float &f1);
   void log(log_t type_in, const std::string &s1, const double &d1);
   void log(log_t type_in, const std::string &s1, const float &f1, const float &f2);
   void log(log_t type_in, const std::string &s1, const double &d1, const std::string &s2, const double &d2);
   void log(log_t type_in, const std::string &s1, const std::string &s2, const std::string &s3, const float &f1, const std::string &s4, const float &f2);
   void log(log_t type_in, const std::string &s1, const float &f1, const std::string &s2, const std::string &s3, const std::string &s4, const float &f2);
   void log(log_t type_in, ltw l1, ltw l2);
   void log(log_t type_in, ltw l1, ltw l2, ltw l3);
   void log(log_t type_in, const ltw &l1, const ltw &l2, const ltw &l3, const ltw &l4);
   void log(log_t type_in, ltw l1, ltw l2, ltw l3, ltw l4, ltw l5);
   void log(log_t type_in, ltw l1, ltw l2, ltw l3, ltw l4, ltw l5, ltw l6);
   void log(log_t type_in,
            const ltw &l1, const ltw &l2, const ltw &l3, const ltw &l4, const ltw &l5,
            const ltw &l6, const ltw &l7, const ltw &l8, const ltw &l9, const ltw &l10);
   void log(log_t type_in,
            const ltw &l1, const ltw &l2, const ltw &l3, const ltw &l4, const ltw &l5,
            const ltw &l6, const ltw &l7, const ltw &l8, const ltw &l9, const ltw &l10, const ltw &l11);

   void log(log_t type_in, const std::vector<ltw> &ls);

   // 20241211-PE now I want to pass the function name also - let's do those...
   void log(log_t type_in, const function_name_t &fn, const std::string &s1);

   void log(log_t type_in, const function_name_t &fn, const std::string &s1, const std::string &s2);

   void log(log_t type_in, const function_name_t &fn, const std::string &s1, int);

   void log(log_t type_in, const function_name_t &fn, const std::vector<ltw> &v);

   void log(log_t type_in, const function_name_t &fn, const std::string &s1, int i, const std::string &s2);

   void log(log_t type_in, const function_name_t &fn, const std::string &s1, bool b, const std::string &s2, float f);

   void show() const;
   void show_last() const;
   void set_update_notifier_function(void (*func)()) { update_notifier_function = func; }
   std::vector<log_item> get_log_history_from(unsigned int idx_start) const;

   friend std::ostream& operator<<(std::ostream &o, const log_item &li);
};

// logging logging;

#endif // LOGGING_HH
