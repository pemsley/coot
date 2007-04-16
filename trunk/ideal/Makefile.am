# ideal/Makefile.am
# 
# Copyright 2004, 2005, 2006 The University of York
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



# should be set by autoconf
# GSL_LIBS = -lgsl -lgslcblas -lm

lib_LTLIBRARIES=libcoot-ideal.la

# is this needed? no.
# top_srcdir = ..

EXTRA_DIST = notes Makefile.simple molmantest.cc

include_HEADERS =  test-indexing.hh simple-restraint.hh pepflip.hh regularize-minimol.hh

INCLUDES = \
	-I$(top_srcdir)/coords     \
	-I$(top_srcdir)/coot-utils \
	-I$(top_srcdir)/mini-mol   \
	-I$(top_srcdir)/geometry   \
	$(CLIPPER_CXXFLAGS)        \
	$(MMDB_CXXFLAGS)

AM_CXXFLAGS = $(COOT_GSL_FLAGS) $(GSL_CFLAGS)

libcoot_ideal_la_SOURCES = simple-restraint.cc  pepflip.cc regularize-minimol.cc

noinst_PROGRAMS = test-indexing with-geometry with-map
# bin_PROGRAMS = 

with_geometry_SOURCES = with-geometry.cc

with_map_SOURCES = with-map.cc

test_indexing_SOURCES = test-indexing.cc

# note that we need libcoot-coords only because we use 
# atom_selection_container_t asc get_atom_selection(filename).  If we get rid of that 
# in the future the dependency on libcoot-coords goes away.
#
with_geometry_LDADD = \
        ./libcoot-ideal.la                                \
	$(top_builddir)/geometry/libcoot-geometry.la      \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la      \
	$(top_builddir)/coords/libcoot-coords.la          \
	$(top_builddir)/coot-utils/libcoot-shelx.la       \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/coot-utils/libcoot-utils.la       \
	$(CLIPPER_LIBS) $(MMDB_LIBS) $(GSL_LIBS)

with_map_LDADD =  \
        ./libcoot-ideal.la                                 \
	$(top_builddir)/geometry/libcoot-geometry.la       \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
	$(top_builddir)/coords/libcoot-coords.la           \
	$(top_builddir)/coot-utils/libcoot-shelx.la        \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
	$(top_builddir)/coot-utils/libcoot-utils.la        \
	$(CLIPPER_LIBS) $(MMDB_LIBS) $(GSL_LIBS)

test_indexing_LDADD = \
        ./libcoot-ideal.la                                \
	$(top_builddir)/geometry/libcoot-geometry.la      \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la      \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/coot-utils/libcoot-utils.la       \
	$(CLIPPER_LIBS) $(MMDB_LIBS) $(GSL_LIBS)

