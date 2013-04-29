# geometry/Makefile.am
# 
# Copyright 2002, 2003 The University of York
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

lib_LTLIBRARIES=libcoot-geometry.la

EXTRA_DIST = notes Makefile.simple

include_HEADERS =  protein-geometry.hh read-sbase.hh lbg-graph.hh srs-geometry.hh

INCLUDES = -I$(top_srcdir)/coords     \
	   -I$(top_srcdir)/utils      \
	   -I$(top_srcdir)/mini-mol   \
	   -I$(top_srcdir)/coot-utils \
	   -I$(top_srcdir)/compat     \
	 $(CLIPPER_CXXFLAGS) $(CCP4SRS_CXXFLAGS) $(MMDB_CXXFLAGS)

AM_CXXFLAGS = \
	-DDATADIR='"$(datadir)"'       \
	-DPKGDATADIR='"$(pkgdatadir)"' \
	$(COOT_WINDOWS_CFLAGS)


libcoot_geometry_la_SOURCES = protein-geometry.cc protein-geometry-mmdb.cc \
    read-sbase.cc energy-lib.cc 	 \
    chem_mods.cc \
    lbg-graph.cc  dreiding.cc


#libcoot_geometry_la_LIBADD = \
#        $(top_builddir)/utils/libcoot-utils.la \
#        $(MMDB_LIBS) $(CLIPPER_LIBS)

libcoot_geometry_la_LIBADD = \
        $(top_builddir)/mini-mol/libcoot-mini-mol.la \
        $(top_builddir)/utils/libcoot-utils.la       \
        $(CLIPPER_LIBS) $(CCP4SRS_LIBS) $(MMDB_LIBS)

# libcoot_geometry_la_CXXFLAGS = -g -O1

noinst_PROGRAMS = test-geometry test-read-components

test_geometry_SOURCES = test-geometry.cc

test_read_components_SOURCES = test-read-components.cc test-read-components.hh

test_geometry_LDADD = ./libcoot-geometry.la    \
        $(top_builddir)/mini-mol/libcoot-mini-mol.la \
        $(top_builddir)/utils/libcoot-utils.la       \
	$(CLIPPER_LIBS) $(CCP4SRS_LIBS) $(MMDB_LIBS)

test_read_components_LDADD =  \
	./libcoot-geometry.la \
        $(top_builddir)/mini-mol/libcoot-mini-mol.la \
        $(top_builddir)/utils/libcoot-utils.la       \
	$(CLIPPER_LIBS) $(CCP4SRS_LIBS) $(MMDB_LIBS)

