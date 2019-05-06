# coords/Makefile.am
# 
# Copyright 2001, 2002, 2003, 2006, The University of York
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

lib_LTLIBRARIES=libcoot-coords.la

pkginclude_HEADERS =  Bond_lines.h Cartesian.h cos-sin.h \
	mmdb.h mmdb-extras.h utils.h mmdb-crystal.h Bond_lines_ext.h \
	coot-close.hh phenix-geo.hh ramachandran-container.hh \
	loop-path.hh

pkgincludedir = $(includedir)/coot/coords

AM_CPPFLAGS = \
	-I$(top_srcdir)                 \
	$(ENHANCED_LIGAND_TOOLS_CFLAGS) \
	$(RDKIT_CXXFLAGS) $(BOOST_CPPFLAGS)	\
	$(CLIPPER_CXXFLAGS) $(MMDB_CXXFLAGS)

AM_CXXFLAGS = \
	-DPKGDATADIR='"$(pkgdatadir)"' \
	-DDATADIR='"$(datadir)"'   	\
	$(COOT_WINDOWS_CFLAGS)

libcoot_coords_la_SOURCES =  Bond_lines.cc Cartesian.cc cos-sin.cc \
		 mmdb-extras.cc utils.cc mmdb.cc mmdb-crystal.cc \
		 Bond_lines_ext.cc coot-close.cc phenix-geo.cc phenix-geo-bonds.cc \
		 rotamer-markup.cc \
		 loop-path.cc

libcoot_coords_la_LIBADD = \
   $(top_builddir)/utils/libcoot-utils.la            \
   $(top_builddir)/mini-mol/libcoot-mini-mol.la      \
   $(top_builddir)/geometry/libcoot-geometry.la      \
   $(top_builddir)/coot-utils/libcoot-coord-utils.la \
   $(top_builddir)/lidia-core/libcoot-lidia-core.la  \
   $(top_builddir)/ligand/libcoot-ligand.la


#   $(MMDB_LIBS) \
#   $(RDKIT_LIBS) $(BOOST_LDFLAGS) $(L_BOOST_PYTHON) $(PYTHON_LIBS) -lssl -lcrypto -lldap -lz

libcoot_coords_la_LDFLAGS = $(SHARED_LDFLAGS)

check_PROGRAMS = test-phenix-geo test-transform

test_transform_SOURCES = test-transform.cc

test_transform_LDADD = libcoot-coords.la

test_phenix_geo_SOURCES = test-phenix-geo.cc

test_phenix_geo_LDADD = libcoot-coords.la \
   $(top_builddir)/coot-utils/libcoot-coord-utils.la \
   $(top_builddir)/utils/libcoot-utils.la           
