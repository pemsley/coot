# build/Makefile.am
# 
# Copyright 2001, 2002, 2003, 2006 The University of York
# Author: Paul Emsley
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
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

lib_LTLIBRARIES=libcoot-cabuild.la

pkginclude_HEADERS = CalphaBuild.hh
pkgincludedir = $(includedir)/coot/build

INCLUDES =  \
	-I$(top_srcdir) \
	$(CLIPPER_CXXFLAGS) \
	$(MMDB_CXXFLAGS)

noinst_PROGRAMS = testcabuild

testcabuild_SOURCES = testcabuild.cc

libcoot_cabuild_la_SOURCES = CalphaBuild.cc

libcoot_cabuild_la_LIBADD =  \
	$(top_builddir)/angles/libcoot-AngleInfo.la \
        $(CLIPPER_LIBS)

testcabuild_LDADD = \
	./libcoot-cabuild.la \
	$(top_builddir)/angles/libcoot-AngleInfo.la \
	$(CLIPPER_LIBS) $(MMDB_LIBS)

EXTRA_DIST = build.txt
