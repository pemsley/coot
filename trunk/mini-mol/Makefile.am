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

include_HEADERS =  mini-mol.hh mini-mol-utils.hh atom-quads.hh

INCLUDES =  \
	-I$(top_srcdir)/coot-utils \
	-I$(top_srcdir)/compat \
	$(CLIPPER_CXXFLAGS)    \
	$(MMDB_CXXFLAGS)

AM_CXXFLAGS = $(COOT_WINDOWS_CFLAGS)

noinst_PROGRAMS = testminimol

testminimol_SOURCES = testminimol.cc

libcoot_mini_mol_la_SOURCES =  mini-mol.cc mini-mol-utils.cc

testminimol_LDADD = libcoot-mini-mol.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/coot-utils/libcoot-utils.la \
	$(CLIPPER_LIBS) $(MMDB_LIBS)
