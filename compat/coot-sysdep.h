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

// window magic jiggery pokery
#if defined(WINDOWS_MINGW) || defined(_MSC_VER)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif

// for timeval; needs to be moved before windows.h nowadays...
#include <winsock2.h>

#define HKL HKL_RENAMED
#include <windows.h>
#undef HKL

#define AddAtomA AddAtom
#define GetAtomNameA GetAtomName
#undef V_UNKNOWN
#define V_UNKNOWNA V_UNKNOWN
#undef small
#undef near
#undef far
#ifdef PLANE
#undef PLANE
#endif
// for nomenclature errors
#ifdef IGNORE
#undef IGNORE
#endif
#endif //Windows

// some redefinitions for Visual C++
#if defined _MSC_VER
#define PKGDATADIR "C:/coot/share"
#define snprintf _snprintf
#define S_ISDIR(m)  (((m) & S_IFMT) == S_IFDIR)
#define S_ISREG(m)  (((m) & S_IFMT) == S_IFREG)
#define S_IWUSR S_IWRITE
#define S_IXUSR S_IEXEC
#define S_IRUSR S_IREAD
#endif

