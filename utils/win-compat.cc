
#include <sys/stat.h>

#include "win-compat.hh"

#if defined _MSC_VER
#define S_ISDIR(m)  (((m) & S_IFMT) == S_IFDIR)
#define S_ISREG(m)  (((m) & S_IFMT) == S_IFREG)
#endif

std::string
coot::get_fixed_font() {

   std::string fixed_font_str;
#if defined(WINDOWS_MINGW) || defined(_MSC_VER)
   fixed_font_str = "monospace";
#else
   // fixed_font_str = "fixed";
   fixed_font_str = "Sans 9";
#endif
   return fixed_font_str;
}

bool
coot::is_dir_or_link(const std::string & file_name) {

   bool r = false;
   struct stat buf;
   stat(file_name.c_str(), &buf);
   if (S_ISDIR(buf.st_mode))
       r = true;
#if defined(WINDOWS_MINGW) || defined(_MSC_VER)
#else
   if (S_ISLNK(buf.st_mode))
      r = true;
#endif   
   return r;
}

bool
coot::is_regular_file(const std::string & file_name) {

   struct stat buf;
   stat(file_name.c_str(), &buf);
   bool r = false;
   if (S_ISREG(buf.st_mode))
       r = true;
   return r;
}

std::string
coot::uri_to_file_name(const std::string &uri) {

   std::string r = uri;

   // Why this? https://en.wikipedia.org/wiki/File_URI_scheme
   
#ifdef WINDOWS_MINGW
	 r = uri.substr(8);
#else
	 r = uri.substr(7);
#endif

   return r;

}
