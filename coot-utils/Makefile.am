# coot-utils/Makefile.am
# 
# Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 The University of York
# Copyright 2014 by Medical Research Council
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

lib_LTLIBRARIES= libcoot-coord-utils.la  libcoot-shelx.la libcoot-coord-extras.la \
	libcoot-map-utils.la libcoot-map-heavy.la

pkginclude_HEADERS = peak-search.hh coot-trim.hh \
	coot-coord-utils.hh coot-map-utils.hh coot-coord-extras.hh \
	coot-map-heavy.hh coot-lsq-types.h coot-shelx.hh \
	coot-rama.hh xmap-stats.hh coot-h-bonds.hh lsq-improve.hh coot-hole.hh \
	read-sm-cif.hh bonded-pairs.hh tree.hh helix-analysis.hh elastic.hh emma.hh \
	glyco-torsions.hh coot-least-squares.hh coot-density-stats.hh \
        atom-overlaps.hh reduce.hh secondary-structure-headers.hh plane-utils.hh

pkgincludedir = $(includedir)/coot/coot-utils

AM_CPPFLAGS = \
	-I$(top_srcdir)	\
	$(CLIPPER_CXXFLAGS) 		\
	$(MMDB_CXXFLAGS) 		\
	$(GSL_CFLAGS)			\
	$(WITH_GUILE_CFLAGS)		\
	$(WITH_PYTHON_CFLAGS)

AM_CXXFLAGS = $(COOT_GSL_FLAGS) $(COOT_WINDOWS_CFLAGS) \
	-DPKGDATADIR='"$(pkgdatadir)"'


# for pure coordinates (mmdb) things
# - note to self: don't some of these include usage of clipper::Coord_orth - which
# means that we should add $CLIPPER_LIBS to libcoot_coord_utils_la_LDADD
libcoot_coord_utils_la_SOURCES = coot-coord-utils.cc \
	coot-coord-utils-nucleotides.cc coot-rama.cc coot-h-bonds.cc hole.cc read-sm-cif.cc \
	bonded-pairs.cc helix-analysis.cc elastic.cc glyco-torsions.cc coot-least-squares.cc \
	atom-overlaps.cc reduce.cc secondary-structure-headers.cc plane-utils.cc


#    $(top_builddir)/mini-mol/libcoot-mini-mol.la 
#    $(top_builddir)/utils/libcoot-utils.la 
#    $(CLIPPER_LIBS) $(MMDB_LIBS)

# clean
libcoot_coord_utils_la_LIBADD = \
    $(top_builddir)/geometry/libcoot-geometry.la \
    $(GSL_LIBS)

libcoot_coord_utils_la_LDFLAGS = $(SHARED_LDFLAGS)

# shelx
libcoot_shelx_la_SOURCES = coot-shelx-ins.cc

# clean
libcoot_shelx_la_LIBADD =  \
    libcoot-coord-utils.la \
    $(top_builddir)/utils/libcoot-utils.la

libcoot_shelx_la_LDFLAGS = $(SHARED_LDFLAGS)

# for non-simple coordinates things
libcoot_coord_extras_la_SOURCES = coot-coord-extras.cc coot-coord-lsq.cc \
	coot-tree-extras.cc water-coordination.cc lsq-improve.cc coot-coord-utils-glyco.cc

# needs more, mini-mol, geometry
#     $(top_builddir)/mini-mol/libcoot-mini-mol.la
#    $(top_builddir)/utils/libcoot-utils.la

libcoot_coord_extras_la_LIBADD = \
   libcoot-coord-utils.la 

libcoot_coord_extras_la_LDFLAGS = $(SHARED_LDFLAGS)


# for (map) things that need to link mmdb and clipper:
libcoot_map_utils_la_SOURCES = coot-map-utils.cc peak-search.cc peak-search-from-list.cc trim.cc \
	xmap-stats.cc maps-spherical.cc

# clean
libcoot_map_utils_la_LIBADD = \
   libcoot-coord-utils.la \
   $(CLIPPER_LIBS)

libcoot_map_utils_la_LDFLAGS = $(SHARED_LDFLAGS)

# for map things that need mmdb, clipper, gsl and other things
libcoot_map_heavy_la_SOURCES = coot-map-heavy.cc coot-fffear.cc

#   $(top_builddir)/utils/libcoot-utils.la 

# and mini-mol (just one).
libcoot_map_heavy_la_LIBADD = \
   libcoot-map-utils.la \
   libcoot-coord-utils.la \
   $(GSL_LIBS)

libcoot_map_heavy_la_LDFLAGS = $(SHARED_LDFLAGS)


# test-agg-clust
check_PROGRAMS = test-h-bonds testshelxreader testcootutils test-hole test-elastic test-spherical glyco-torsions \
	test-tt-morph test-coot-probe

# test_zo_chi_table_SOURCES = test-zo-chi-table.cc

# test_agg_clust_SOURCES = test-agg-clust.cc
# test_agg_clust_LDADD =     $(top_builddir)/utils/libcoot-utils.la $(CLIPPER_LIBS)


#	$(top_builddir)/geometry/libcoot-geometry.la 
#        $(top_builddir)/mini-mol/libcoot-mini-mol.la 
#        $(top_builddir)/utils/libcoot-utils.la $(GLOB_LIBS)


test_coot_probe_SOURCES = test-coot-probe.cc
test_coot_probe_LDADD =          \
	./libcoot-coord-utils.la

testshelxreader_SOURCES = test-shelx-reader.cc
testshelxreader_LDADD =          \
	./libcoot-shelx.la       \
	./libcoot-coord-utils.la 

#	$(top_builddir)/geometry/libcoot-geometry.la 
#       $(top_builddir)/mini-mol/libcoot-mini-mol.la 
#       $(top_builddir)/utils/libcoot-utils.la

testcootutils_SOURCES = test-utils.cc
testcootutils_LDADD =            \
	./libcoot-shelx.la       \
	./libcoot-coord-extras.la\
	./libcoot-coord-utils.la \
	$(GSL_LIBS) \
        $(GLOB_LIBS)

#	$(top_builddir)/geometry/libcoot-geometry.la 
#       $(top_builddir)/mini-mol/libcoot-mini-mol.la 
#       $(top_builddir)/utils/libcoot-utils.la 

test_h_bonds_SOURCES = test-h-bonds.cc
test_h_bonds_LDADD = \
	 ./libcoot-coord-utils.la \
        $(GLOB_LIBS)

test_hole_SOURCES = test-hole.cc
test_hole_LDADD = ./libcoot-coord-utils.la

# \
#	$(top_builddir)/utils/libcoot-utils.la $(GLOB_LIBS)
#	$(top_builddir)/geometry/libcoot-geometry.la \
#        $(top_builddir)/mini-mol/libcoot-mini-mol.la 

test_tt_morph_SOURCES = test-morph.cc

test_tt_morph_LDADD = 	./libcoot-map-heavy.la

test_elastic_SOURCES = test-elastic.cc
test_elastic_LDADD   = ./libcoot-coord-utils.la $(GLOB_LIBS)

test_spherical_SOURCES = test-spherical.cc
test_spherical_LDADD   = ./libcoot-coord-utils.la $(GSL_LIBS)

glyco_torsions_SOURCES = test-glyco-torsions.cc
glyco_torsions_LDADD   = ./libcoot-coord-utils.la

#    $(MMDB_LIBS)
