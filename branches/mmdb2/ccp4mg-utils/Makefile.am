# ccp4mg-utils/Makefile.am
# 
# Copyright 2005 The University of York
# Author: Paul Emsley
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA

lib_LTLIBRARIES=libccp4mg-utils.la

pkginclude_HEADERS = geomutil.h mgtree.h cartesian.h matrix.h quat.h ppmutil.h
pkgincludedir = $(includedir)/coot/ccp4mg-utils

EXTRA_DIST = convert.sh get-new-utils.sh

INCLUDES =  \
	$(MMDB_CXXFLAGS) \
	$(LIBPNG_CFLAGS)

libccp4mg_utils_la_SOURCES = geomutil.cc mgtree.cc cartesian.cc \
	matrix.cc quat.cc ppmutil.cc


libccp4mg_utils_la_LIBADD = \
   $(LIBPNG_LIBS)

# -lgomp



