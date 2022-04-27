// vim: sw=4 ts=4 sts=4 expandtab

/*
 * This provides C wrappers for platform-specific functions
 * from coot-sysdep.h.
 */

#ifndef _COMPAT_COOT_SYSDEP_C_H
#define _COMPAT_COOT_SYSDEP_C_H

#include <stddef.h>

typedef struct {
    char ** found;
    size_t n_found;
} Coot_C_Found;

extern "C" {
    void coot_c_free_found(Coot_C_Found *found);
    Coot_C_Found coot_c_gather_files_by_patterns(const char *dir_path, const char **patterns, size_t n_patterns);
}

#endif // _COMPAT_COOT_SYSDEP_C_H
