/* compat/coot-sysdep.h
 * 
 * Copyright 2008, The University of York
 * Author: Bernhard Lohkamp
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
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

#ifndef _COMPAT_COOT_SYSDEP_H
#define _COMPAT_COOT_SYSDEP_H

#include <string>
#include <vector>

namespace coot {
    std::string current_working_dir();
    std::vector<std::string> gather_files_by_patterns(const std::string &dir_path, const std::vector<std::string> &pattern);
    std::string get_fixed_font();
    bool is_dir(const std::string &file_path);
    bool is_link(const std::string &file_path);
    bool is_regular_file(const std::string &file_path);
    bool rename(const char *old_file_path, const char *new_file_path, std::string &error_message);
    void set_os_error_mode();
} // namespace coot

#if defined(COOT_BUILD_WINDOWS)
# include "coot-win32.h"
#elif defined(COOT_BUILD_POSIX)
# include "coot-posix.h"
#else
# error "Misdetected or unsupported platform"
#endif // COOT_

#endif // _COMPAT_COOT_SYSDEP_H
