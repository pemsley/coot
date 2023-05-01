
#include <sys/stat.h>
#include <iostream>

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

#ifdef WINDOWS_MINGW
int
coot::rename_win(const char *old_filename, const char *new_filename) {
    // BL says:: on windows (non POSIX) rename wont overwrite, so
    // need to remove first.
    // Could just use coot::rename for all... maybe needs some rewrite then
    // and Windows check here.
    // return 0 on success

    int ret = -1;
    if (access(new_filename, F_OK) != 0) {
        // normal rename if new_file does not exist.
        ret = rename(old_filename, new_filename);
    } else {
        // file exists, make backup then remove file to be replaced
        std::string backup_fn = new_filename;
        backup_fn += ".bak";
        if (access(backup_fn.c_str(), F_OK) == 0) {
            // backup file exsists - shouldnt, remove first
            int status_rm = remove(backup_fn.c_str());
            // assume ok - for now FIXME
        }
        int status_back = rename(new_filename,
                                backup_fn.c_str());
        if (status_back != 0) {
           std::cout << "WARNING:: rename failed. Cannot make backup of "
                     << new_filename << std::endl;
           ret = status_back;
        } else {
            // now can do a rename (with old file out of the way
            int status = rename(old_filename, new_filename);
            if (status != 0) {
                std::cout << "WARNING:: rename status " << status
                          << " failed to rename to " << new_filename << std::endl;
                std::cout << "restore from backup" << std::endl;
                // restore
                int status_rest = rename(backup_fn.c_str(),
                                         new_filename);
                if (status_rest != 0) {
                    std::cout << "WARNING:: oh dear, failed to restore from backup..." << std::endl;
                }
            } else {
                std::cout << "debug:: renaming successful" << std::endl;
                // simple remove backup file (assume works since just created)
                int status_backrm = remove(backup_fn.c_str());
                std::cout << "INFO:: remove backup file status " << status_backrm << std::endl;
            }
            ret = status;
        }
    }

    return ret;

}
#endif // WINDOWS_MINGW
