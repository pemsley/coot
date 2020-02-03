
#ifndef COMMAND_ARG_HH
#define COMMAND_ARG_HH

#include "utils/coot-utils.hh"

namespace coot {
   class command_arg_t {
   public:
      enum coot_script_arg_type{UNSET, INT, FLOAT, STRING, BOOL};
      command_arg_t(int iin) {
         i = iin;
         type = INT;
         b = false;
      }
      command_arg_t(float fin) {
         f = fin;
         type = FLOAT;
         b = false;
      }
      command_arg_t(const clipper::String &sin) {
         s = sin;
         type = STRING;
         b = false;
      }
      command_arg_t(const std::string &sin) {
         s = sin;
         type = STRING;
         b = false;
      }
      command_arg_t(const char *sin) {
         s = sin;
         type = STRING;
         b = false;
      }
      command_arg_t(bool bin) {
         b = bin;
         type = BOOL;
      }
      command_arg_t() {
         type = UNSET;
         b = false;
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
         os = s;
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
