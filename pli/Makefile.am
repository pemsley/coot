# pli/Makefile.am
# 
# Copyright 2017 by Medical Research Council
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

# Python is used by specs.cc/hh
AM_CPPFLAGS = \
	-DUSE_PYTHON \
	-I$(top_srcdir) \
	$(ENHANCED_LIGAND_TOOLS_CPPFLAGS) \
	$(CLIPPER_CXXFLAGS)             \
	$(MMDB_CXXFLAGS) \
	$(PYTHON_CPPFLAGS)  \
	$(RDKIT_CXXFLAGS) $(BOOST_CPPFLAGS) \
	$(COOT_WINDOWS_CFLAGS)

lib_LTLIBRARIES=libcoot-pli.la

# put specs into a python-loving environment when we need access to coot-pli from python

# remove (python-using) specs from coot-pli
libcoot_pli_la_SOURCES = protein-ligand-interactions.cc flev-annotations.cc pi-stacking.cc \
	sdf-interface.hh         \
	sdf-internal.hh          \
	sdf-interface.cc

libcoot_pli_la_LIBADD = $(top_builddir)/coot-utils/libcoot-coord-utils.la \
                        $(top_builddir)/geometry/libcoot-geometry.la      \
                        $(top_builddir)/lidia-core/libcoot-lidia-core.la  \
                        $(CLIPPER_LIBS) $(MMDB_LIBS)                      \
                        $(RDKIT_LIBS) $(BOOST_LDFLAGS)

libcoot_pli_la_LDFLAGS = $(SHARED_LDFLAGS)

pkginclude_HEADERS = flev-annotations.hh protein-ligand-interactions.hh pi-stacking.hh specs.hh \
	sdf-interface.hh sdf-internal.hh sdf-interface-for-export.hh

pkgincludedir = $(includedir)/coot/pli
