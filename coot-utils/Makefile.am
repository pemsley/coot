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

lib_LTLIBRARIES= libcoot-coord-utils.la  libcoot-coord-extras.la \
	libcoot-map-utils.la libcoot-map-heavy.la

pkginclude_HEADERS = peak-search.hh coot-trim.hh \
	coot-coord-utils.hh coot-map-utils.hh coot-coord-extras.hh \
	coot-map-heavy.hh coot-lsq-types.h coot-shelx.hh \
	coot-rama.hh xmap-stats.hh coot-h-bonds.hh lsq-improve.hh coot-hole.hh \
	read-sm-cif.hh bonded-pairs.hh tree.hh helix-analysis.hh elastic.hh emma.hh \
	glyco-torsions.hh coot-least-squares.hh coot-density-stats.hh \
        atom-overlaps.hh reduce.hh atom-selection-container.hh contact-info.hh plane-utils.hh \
        secondary-structure-headers.hh c-beta-deviations.hh jed-flip.hh lidia-core-functions.hh \
	bonded-atoms.hh spherical-harmonics.hh helix-like.hh merge-atom-selections.hh \
	contacts-by-bricks.hh amp-reso.hh stack-and-pair.hh

pkgincludedir = $(includedir)/coot/coot-utils

AM_CPPFLAGS = \
	-I$(top_srcdir)	\
	$(ENHANCED_LIGAND_TOOLS_CFLAGS) \
	$(CLIPPER_CXXFLAGS) 		\
	$(MMDB_CXXFLAGS) 		\
	$(GSL_CFLAGS)			\
	$(WITH_GUILE_CFLAGS)		\
	$(RDKIT_CXXFLAGS) 		\
	$(WITH_PYTHON_CFLAGS)

AM_CXXFLAGS = $(COOT_GSL_FLAGS) $(COOT_WINDOWS_CFLAGS) \
	-DPKGDATADIR='"$(pkgdatadir)"'


# for pure coordinates (mmdb) things
# - note to self: don't some of these include usage of clipper::Coord_orth - which
# means that we should add $CLIPPER_LIBS to libcoot_coord_utils_la_LDADD
libcoot_coord_utils_la_SOURCES = coot-coord-utils.cc \
	coot-coord-utils-nucleotides.cc coot-rama.cc coot-h-bonds.cc hole.cc read-sm-cif.cc \
	bonded-pairs.cc helix-analysis.cc elastic.cc glyco-torsions.cc coot-least-squares.cc \
	atom-overlaps.cc reduce.cc atom-selection-container.cc contact-info.cc coot-shelx-ins.cc \
	secondary-structure-headers.cc plane-utils.cc c-beta-deviations.cc lidia-core-functions.cc \
	atom-tools.cc bonded-atoms.cc bricks.cc helix-like.cc merge-atom-selections.cc \
	contacts-by-bricks.cc stack-and-pair.cc


#    $(top_builddir)/mini-mol/libcoot-mini-mol.la 
#    $(top_builddir)/utils/libcoot-utils.la 
#    $(CLIPPER_LIBS) $(MMDB_LIBS)

libcoot_coord_utils_la_LIBADD = \
    $(top_builddir)/lidia-core/libcoot-lidia-core.la \
    $(top_builddir)/geometry/libcoot-geometry.la \
    $(GSL_LIBS) $(RDKIT_LIBS) $(BOOST_LDFLAGS) $(L_BOOST_PYTHON) $(PYTHON_LIBS)

libcoot_coord_utils_la_LDFLAGS = $(SHARED_LDFLAGS) $(CLIPPER_LIBS) $(MMDB_LIBS)

# shelx
# libcoot_shelx_la_SOURCES = coot-shelx-ins.cc

# clean
# libcoot_shelx_la_LIBADD =  \
# libcoot-coord-utils.la \
# $(top_builddir)/utils/libcoot-utils.la

# ibcoot_shelx_la_LDFLAGS = $(SHARED_LDFLAGS)

# for non-simple coordinates things
libcoot_coord_extras_la_SOURCES = coot-coord-extras.cc coot-coord-lsq.cc \
	coot-tree-extras.cc water-coordination.cc lsq-improve.cc coot-coord-utils-glyco.cc jed-flip.cc

# needs more, mini-mol, geometry
#     $(top_builddir)/mini-mol/libcoot-mini-mol.la
#    $(top_builddir)/utils/libcoot-utils.la

# should libcoot-coord-utils link $CLIPPER_LIBS?
libcoot_coord_extras_la_LIBADD = \
   libcoot-coord-utils.la 

libcoot_coord_extras_la_LDFLAGS = $(SHARED_LDFLAGS)


# for (map) things that need to link mmdb and clipper:
libcoot_map_utils_la_SOURCES = coot-map-utils.cc peak-search.cc peak-search-from-list.cc trim.cc \
	xmap-stats.cc maps-spherical.cc

# should libcoot-coord-utils link $CLIPPER_LIBS directly? (Hence CLIPPER_LIBS is not needed here)
libcoot_map_utils_la_LIBADD = \
   libcoot-coord-utils.la \
   $(CLIPPER_LIBS)

libcoot_map_utils_la_LDFLAGS = $(SHARED_LDFLAGS)

# for map things that need mmdb, clipper, gsl and other things
libcoot_map_heavy_la_SOURCES = coot-map-heavy.cc coot-fffear.cc spherical-harmonics.cc segmap.cc

#   $(top_builddir)/utils/libcoot-utils.la 

# and mini-mol (just one).
libcoot_map_heavy_la_LIBADD = \
   libcoot-map-utils.la \
   libcoot-coord-utils.la \
   $(GSL_LIBS)

libcoot_map_heavy_la_LDFLAGS = $(SHARED_LDFLAGS)

# lidia-core functions moved here because they use coot-utils functions
#
check_PROGRAMS = test-h-bonds testshelxreader test-coot-utils test-hole test-elastic test-spherical glyco-torsions \
	test-tt-morph test-coot-probe test-chem-feat test-cod-atom-types dictionary-atom-types cod-db-access \
	test-emma test-convolute test-segmap


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
testshelxreader_LDADD =  ./libcoot-coord-utils.la 

test_emma_SOURCES = test-emma.cc
test_emma_LDADD =  libcoot-map-utils.la \
	$(top_builddir)/coords/libcoot-coords.la

test_convolute_SOURCES = test-convolute.cc
test_convolute_LDADD = \
        ./libcoot-map-utils.la \
	./libcoot-map-heavy.la \
	$(top_builddir)/coords/libcoot-coords.la

test_segmap_SOURCES = test-segmap.cc
test_segmap_LDADD = \
        ./libcoot-map-utils.la \
	./libcoot-map-heavy.la \
	$(top_builddir)/coords/libcoot-coords.la

# other test binaries may need to be linked like this too
test_coot_utils_SOURCES = test-utils.cc
test_coot_utils_LDADD = \
	./libcoot-map-heavy.la    \
	./libcoot-map-utils.la    \
	./libcoot-coord-extras.la \
	./libcoot-coord-utils.la  \
	$(top_builddir)/coords/libcoot-coords.la \
	$(CLIPPER_LIBS) $(CCP4SRS_LIBS)

test_h_bonds_SOURCES = test-h-bonds.cc
test_h_bonds_LDADD = \
	 ./libcoot-coord-utils.la \
        $(GLOB_LIBS)

test_hole_SOURCES = test-hole.cc
test_hole_LDADD = ./libcoot-coord-utils.la

test_tt_morph_SOURCES = test-morph.cc
test_tt_morph_LDADD = 	./libcoot-map-heavy.la

test_elastic_SOURCES = test-elastic.cc
test_elastic_LDADD   = ./libcoot-coord-utils.la $(GLOB_LIBS)

test_spherical_SOURCES = test-spherical.cc
test_spherical_LDADD   = ./libcoot-coord-utils.la $(CLIPPER_LIBS) $(GSL_LIBS)

glyco_torsions_SOURCES = test-glyco-torsions.cc
glyco_torsions_LDADD   = ./libcoot-coord-utils.la

#    $(MMDB_LIBS)

test_chem_feat_SOURCES = test-chem-feat.cc

test_chem_feat_LDADD = \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(RDKIT_LIBS)

test_cod_atom_types_SOURCES = test-cod-atom-types.cc

test_cod_atom_types_LDADD = \
	libcoot-coord-utils.la \
	$(RDKIT_LIBS)

cod_db_access_SOURCES = cod-db-access.cc

cod_db_access_LDADD = $(top_builddir)/lidia-core/libcoot-lidia-core.la $(RDKIT_LIBS)

dictionary_atom_types_SOURCES = dictionary-atom-types.cc

dictionary_atom_types_LDADD = \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/lidia-core/libcoot-lidia-core.la \
	$(top_builddir)/geometry/libcoot-geometry.la \
        $(RDKIT_LIBS)
