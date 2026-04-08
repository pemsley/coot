#ifndef FILE_SYSTEM_UTILS_HH
#define FILE_SYSTEM_UTILS_HH

#include <string>

/* section File System Functions (string) */
/*!  \name File System Functions (string) */
/*! \{ */

/*! \brief make a directory dir (if it doesn't exist) and return error code

   If it can be created, create the directory dir, return the success status
   like mkdir: mkdir

   @return zero on success, or -1 if an  error  occurred.
   If dir already exists as a directory, return 0 of course.
 */
int make_directory_maybe(const std::string &dir);

/*! \} */


#endif // FILE_SYSTEM_UTILS_HH

