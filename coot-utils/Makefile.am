# coot-utils/Makefile.am
# 
# Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 The University of York
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


WITH_GUILE_CFLAGS=$(COOT_USE_GUILE)

if COOT_USE_PYTHON
WITH_PYTHON_CFLAGS=-DUSE_PYTHON
endif

lib_LTLIBRARIES= libcoot-coord-utils.la  libcoot-shelx.la libcoot-coord-extras.la libcoot-map-utils.la libcoot-map-heavy.la

include_HEADERS = peak-search.hh coot-trim.hh \
	coot-coord-utils.hh coot-map-utils.hh coot-coord-extras.hh \
	coot-map-heavy.hh coot-lsq-types.h coot-shelx.hh \
	coot-rama.hh xmap-stats.hh coot-h-bonds.hh lsq-improve.hh coot-hole.hh read-sm-cif.hh
INCLUDES = \
	-I$(top_srcdir)/utils	        \
	-I$(top_srcdir)/geometry	\
	-I$(top_srcdir)/compat		\
	-I$(top_srcdir)/mini-mol	\
	$(CLIPPER_CXXFLAGS) 		\
	$(MMDB_CXXFLAGS) 		\
	$(GSL_CFLAGS)			\
	$(WITH_GUILE_CFLAGS)		\
	$(WITH_PYTHON_CFLAGS)

AM_CXXFLAGS = $(COOT_GSL_FLAGS) $(COOT_WINDOWS_CFLAGS) \
	-DPKGDATADIR='"$(pkgdatadir)"'


# for pure coordinates (mmdb) things
libcoot_coord_utils_la_SOURCES = coot-coord-utils.cc \
	coot-coord-utils-nucleotides.cc coot-rama.cc coot-h-bonds.cc hole.cc read-sm-cif.cc \
	coot-coord-utils-glyco.cc

# clean
libcoot_coord_utils_la_LIBADD = \
    $(top_builddir)/geometry/libcoot-geometry.la \
    $(top_builddir)/utils/libcoot-utils.la \
    $(CLIPPER_LIBS) $(MMDB_LIBS)

# shelx
libcoot_shelx_la_SOURCES = coot-shelx-ins.cc

# clean
libcoot_shelx_la_LIBADD =  \
    libcoot-coord-utils.la \
    $(top_builddir)/utils/libcoot-utils.la

# for non-simple coordinates things
libcoot_coord_extras_la_SOURCES = coot-coord-extras.cc coot-coord-lsq.cc \
	coot-tree-extras.cc water-coordination.cc lsq-improve.cc

# needs more, mini-mol, geometry
libcoot_coord_extras_la_LIBADD = \
   libcoot-coord-utils.la \
    $(top_builddir)/mini-mol/libcoot-mini-mol.la \
    $(top_builddir)/utils/libcoot-utils.la

# for (map) things that need to link mmdb and clipper:
libcoot_map_utils_la_SOURCES = coot-map-utils.cc peak-search.cc peak-search-from-list.cc trim.cc xmap-stats.cc

# clean
libcoot_map_utils_la_LIBADD = \
   libcoot-coord-utils.la

# for map things that need mmdb, clipper, gsl and other things
libcoot_map_heavy_la_SOURCES = coot-map-heavy.cc coot-fffear.cc

# and mini-mol (just one).
libcoot_map_heavy_la_LIBADD = \
   libcoot-map-utils.la \
   libcoot-coord-utils.la \
   $(top_builddir)/mini-mol/libcoot-mini-mol.la \
   $(top_builddir)/utils/libcoot-utils.la \
   $(GSL_LIBS)


noinst_PROGRAMS = test-h-bonds testshelxreader testcootutils test-hole

# testshelxreader fails to link for some reason! coot::protein_geometry::get_h_bond_type 
# testcootutils also.  20100812 Fixed now - added libcoot-geometry.la


testshelxreader_SOURCES = test-shelx-reader.cc
testshelxreader_LDADD =          \
	./libcoot-shelx.la       \
	./libcoot-coord-utils.la \
	$(top_builddir)/geometry/libcoot-geometry.la \
        $(top_builddir)/utils/libcoot-utils.la

testcootutils_SOURCES = test-utils.cc
testcootutils_LDADD =            \
	./libcoot-shelx.la       \
	./libcoot-coord-extras.la\
	./libcoot-coord-utils.la \
	$(top_builddir)/geometry/libcoot-geometry.la \
        $(top_builddir)/utils/libcoot-utils.la

test_h_bonds_SOURCES = test-h-bonds.cc
test_h_bonds_LDADD = \
	 ./libcoot-coord-utils.la \
	$(top_builddir)/geometry/libcoot-geometry.la \
        $(top_builddir)/utils/libcoot-utils.la

test_hole_SOURCES = test-hole.cc
test_hole_LDADD = ./libcoot-coord-utils.la \
	$(top_builddir)/geometry/libcoot-geometry.la \
	$(top_builddir)/utils/libcoot-utils.la


#    $(MMDB_LIBS)
