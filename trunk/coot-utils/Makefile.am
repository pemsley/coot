# coot-utils/Makefile.am
# 
# Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 The University of York
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

if COOT_USE_GUILE
WITH_GUILE_CFLAGS=-DUSE_GUILE
else
if COOT_USE_PYTHON
WITH_PYTHON_CFLAGS=-DUSE_PYTHON
endif
endif

lib_LTLIBRARIES=libcoot-utils.la libcoot-map-utils.la libcoot-coord-utils.la  libcoot-coord-extras.la libcoot-map-heavy.la libcoot-shelx.la

include_HEADERS = coot-utils.hh peak-search.hh coot-trim.hh \
	coot-coord-utils.hh coot-map-utils.hh coot-coord-extras.hh \
	coot-fasta.hh coot-map-heavy.hh coot-lsq-types.h coot-shelx.hh

INCLUDES = \
	-I$(top_srcdir)/coords  	\
	-I$(top_srcdir)/geometry	\
	$(CLIPPER_CXXFLAGS) 		\
	$(MMDB_CXXFLAGS) 		\
	$(GSL_CFLAGS)			\
	$(WITH_GUILE_CFLAGS)		\
	$(WITH_PYTHON_CFLAGS)

AM_CXXFLAGS = $(COOT_GSL_FLAGS) $(COOT_WINDOWS_CFLAGS) \
	-DPKGDATADIR='"$(pkgdatadir)"'

# for simple things:
libcoot_utils_la_SOURCES = coot-utils.cc coot-fasta.cc atomic-numbers.cc \
	gdk-keysyms.cc

# for pure coordinates (mmdb) things
libcoot_coord_utils_la_SOURCES = coot-coord-utils.cc \
	coot-coord-utils-nucleotides.cc

# shelx
libcoot_shelx_la_SOURCES = coot-shelx-ins.cc

# for non-simple coordinates things
libcoot_coord_extras_la_SOURCES = coot-coord-extras.cc coot-coord-lsq.cc

# for (map) things that need to link mmdb and clipper:
libcoot_map_utils_la_SOURCES = coot-map-utils.cc peak-search.cc trim.cc 

# for map things that need mmdb, clipper, gsl and other things
libcoot_map_heavy_la_SOURCES = coot-map-heavy.cc coot-fffear.cc

noinst_PROGRAMS = testshelxreader testcootutils

testshelxreader_SOURCES = test-shelx-reader.cc
testshelxreader_LDADD = \
	./libcoot-shelx.la \
	./libcoot-coord-utils.la \
	./libcoot-utils.la \
	$(CLIPPER_LIBS) $(MMDB_LIBS)

testcootutils_SOURCES = test-utils.cc
testcootutils_LDADD = \
	./libcoot-shelx.la \
	./libcoot-coord-utils.la \
	./libcoot-utils.la \
	$(CLIPPER_LIBS) $(MMDB_LIBS)


