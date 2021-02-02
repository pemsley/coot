# high-res/Makefile.am
# 
# Copyright 2004 The University of York
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

lib_LTLIBRARIES=libcoot-high-res.la

pkginclude_HEADERS = high-res.hh coot-atom-graph.hh coot-node-info.hh \
	sequence-assignment.hh
pkgincludedir = $(includedir)/coot/high-res


AM_CXXFLAGS = $(COOT_GSL_FLAGS) $(GSL_CFLAGS)

AM_CPPFLAGS = \
	-I$(top_srcdir)   \
	$(MMDB_CXXFLAGS)  \
	$(CLIPPER_CXXFLAGS)

libcoot_high_res_la_SOURCES = high-res.cc coot-atom-graph.cc residue-distortions.cc \
	sequence-assignment.cc

libcoot_high_res_la_LIBADD = \
   $(top_builddir)/coot-utils/libcoot-coord-utils.la \
   $(top_builddir)/coords/libcoot-coords.la \
   $(top_builddir)/ligand/libcoot-ligand.la \
   $(top_builddir)/coot-utils/libcoot-map-utils.la   \
   $(top_builddir)/mini-mol/libcoot-mini-mol.la \
   $(top_builddir)/utils/libcoot-utils.la \
   $(CLIPPER_LIBS) $(MMDB_LIBS)

libcoot_high_res_la_LDFLAGS = $(SHARED_LDFLAGS)


#check_PROGRAMS = testatomgraph tracehighres testseqalign testseqassignment globularize_protein

check_PROGRAMS = testatomgraph tracehighres

testatomgraph_SOURCES = test-atom-graph.cc

testatomgraph_LDADD =  libcoot-high-res.la \
   $(top_builddir)/mini-mol/libcoot-mini-mol.la\
   $(MMDB_LIBS)

#testatomgraph_LDADD =  libcoot-high-res.la                \
#	$(top_builddir)/ligand/libcoot-ligand.la          \
#	$(top_builddir)/db-main/libcoot-db-main.la        \
#	$(top_builddir)/ideal/libcoot-ideal.la            \
#	$(top_builddir)/geometry/libcoot-geometry.la      \
#	$(top_builddir)/coords/libcoot-coords.la          \
#	$(top_builddir)/mini-mol/libcoot-mini-mol.la      \
#	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
#	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
#	$(top_builddir)/coot-utils/libcoot-map-utils.la   \
#	$(top_builddir)/coot-utils/libcoot-shelx.la       \
#	$(top_builddir)/utils/libcoot-utils.la       \
#	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la   \
#	$(top_builddir)/compat/libcoot-compat.la          \
#	$(GSL_LIBS)			                  \
#	$(CLIPPER_LIBS)

#testseqassignment_SOURCES = testseqassignment.cc
#testseqassignment_LDADD = libcoot-high-res.la             \
#	$(top_builddir)/ligand/libcoot-ligand.la          \
#	$(top_builddir)/db-main/libcoot-db-main.la        \
#	$(top_builddir)/ideal/libcoot-ideal.la            \
#	$(top_builddir)/geometry/libcoot-geometry.la      \
#	$(top_builddir)/coords/libcoot-coords.la          \
#	$(top_builddir)/mini-mol/libcoot-mini-mol.la      \
#	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
#	$(top_builddir)/coot-utils/libcoot-map-utils.la   \
#	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
#	$(top_builddir)/coot-utils/libcoot-shelx.la       \
#	$(top_builddir)/utils/libcoot-utils.la       \
#	$(top_builddir)/compat/libcoot-compat.la          \
#	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la   \
#	$(GSL_LIBS)				          \
#	$(CLIPPER_LIBS)

## tracehighres_LDFLAGS = -static 

tracehighres_SOURCES = trace-high-res.cc

tracehighres_LDADD = libcoot-high-res.la \
   $(top_builddir)/coot-utils/libcoot-map-utils.la   \
   $(top_builddir)/mini-mol/libcoot-mini-mol.la      \
   $(top_builddir)/utils/libcoot-utils.la \
   $(top_builddir)/compat/libcoot-compat.la \
   $(CLIPPER_LIBS)


#tracehighres_LDADD = libcoot-high-res.la                  \
#	$(top_builddir)/ligand/libcoot-ligand.la          \
#	$(top_builddir)/ideal/libcoot-ideal.la            \
#	$(top_builddir)/db-main/libcoot-db-main.la        \
#	$(top_builddir)/geometry/libcoot-geometry.la      \
#	$(top_builddir)/coords/libcoot-coords.la          \
#	$(top_builddir)/mini-mol/libcoot-mini-mol.la      \
#	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
#	$(top_builddir)/coot-utils/libcoot-map-utils.la   \
#	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
#	$(top_builddir)/coot-utils/libcoot-shelx.la       \
#	$(top_builddir)/utils/libcoot-utils.la       \
#	$(top_builddir)/compat/libcoot-compat.la          \
#	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la   \
#	$(GSL_LIBS)				          \
#	$(CLIPPER_LIBS)


#testseqalign_SOURCES = testseqalign.cc

## why do we need to link with all this (junk)?
#testseqalign_LDADD = libcoot-high-res.la                  \
#	$(top_builddir)/ligand/libcoot-ligand.la          \
#	$(top_builddir)/ideal/libcoot-ideal.la            \
#	$(top_builddir)/db-main/libcoot-db-main.la        \
#	$(top_builddir)/geometry/libcoot-geometry.la      \
#	$(top_builddir)/coords/libcoot-coords.la          \
#	$(top_builddir)/mini-mol/libcoot-mini-mol.la      \
#	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
#	$(top_builddir)/coot-utils/libcoot-map-utils.la   \
#	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
#	$(top_builddir)/coot-utils/libcoot-shelx.la       \
#	$(top_builddir)/utils/libcoot-utils.la       \
#	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la   \
#	$(top_builddir)/compat/libcoot-compat.la          \
#	$(GSL_LIBS)				          \
#	$(CLIPPER_LIBS)


#globularize_protein_SOURCES = globularize.cc

#globularize_protein_LDADD = libcoot-high-res.la           \
#	$(top_builddir)/ligand/libcoot-ligand.la          \
#	$(top_builddir)/db-main/libcoot-db-main.la        \
#	$(top_builddir)/ideal/libcoot-ideal.la            \
#	$(top_builddir)/geometry/libcoot-geometry.la      \
#	$(top_builddir)/coords/libcoot-coords.la          \
#	$(top_builddir)/mini-mol/libcoot-mini-mol.la      \
#	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
#	$(top_builddir)/coot-utils/libcoot-map-utils.la   \
#	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
#	$(top_builddir)/coot-utils/libcoot-shelx.la       \
#	$(top_builddir)/utils/libcoot-utils.la       \
#	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la   \
#	$(top_builddir)/compat/libcoot-compat.la          \
#	$(GSL_LIBS)				          \
#	$(CLIPPER_LIBS)


