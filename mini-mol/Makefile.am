# mini-mol/Makefile.am
# 
# Copyright 2003 The University of York
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

lib_LTLIBRARIES=libcoot-mini-mol.la

pkginclude_HEADERS =  mini-mol.hh mini-mol-utils.hh atom-quads.hh
pkgincludedir = $(includedir)/coot/mini-mol

AM_CPPFLAGS =  \
	-I$(top_srcdir) \
	$(CLIPPER_CXXFLAGS)    \
	$(MMDB_CXXFLAGS)

AM_CXXFLAGS = $(COOT_WINDOWS_CFLAGS)

check_PROGRAMS = testminimol

testminimol_SOURCES = testminimol.cc

libcoot_mini_mol_la_SOURCES =  mini-mol.cc mini-mol-utils.cc atom-quads.cc

# Oh, oops, there is a circular dependency, minimol uses coot-util's
# min_resno_in_chain().  OK, so minimol should go into coot utils and
# coot-utils becomes one big utils lib.
#
libcoot_mini_mol_la_LIBADD = \
	$(top_builddir)/utils/libcoot-utils.la \
	$(CLIPPER_LIBS) $(MMDB_LIBS)

libcoot_mini_mol_la_LDFLAGS = $(SHARED_LDFLAGS)

testminimol_LDADD = libcoot-mini-mol.la \
	$(CLIPPER_LIBS) $(MMDB_LIBS)
