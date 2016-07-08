

#include "win-compat.hh"

std::string
coot::get_fixed_font() {

   std::string fixed_font_str;
#if defined(WINDOWS_MINGW) || defined(_MSC_VER)
   fixed_font_str = "monospace";
#else
   fixed_font_str = "fixed";
   fixed_font_str = "Sans 9";
#endif
   return fixed_font_str;
}
