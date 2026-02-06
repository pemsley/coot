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

WITH_PYTHON_CFLAGS=-DUSE_PYTHON

lib_LTLIBRARIES= libcoot-coord-utils.la  libcoot-coord-extras.la \
	libcoot-map-utils.la libcoot-map-heavy.la

pkginclude_HEADERS = peak-search.hh coot-trim.hh \
	coot-coord-utils.hh coot-map-utils.hh coot-coord-extras.hh \
	coot-map-heavy.hh coot-lsq-types.hh coot-shelx.hh \
	coot-rama.hh xmap-stats.hh coot-h-bonds.hh lsq-improve.hh coot-hole.hh \
	read-sm-cif.hh bonded-pairs.hh tree.hh helix-analysis.hh elastic.hh emma.hh \
	glyco-torsions.hh coot-least-squares.hh coot-density-stats.hh \
        atom-overlaps.hh reduce.hh atom-selection-container.hh contact-info.hh plane-utils.hh \
        secondary-structure-headers.hh c-beta-deviations.hh jed-flip.hh lidia-core-functions.hh \
	bonded-atoms.hh spherical-harmonics.hh helix-like.hh merge-atom-selections.hh \
	contacts-by-bricks.hh amp-reso.hh stack-and-pair.hh atom-tools.hh segmap.hh fast-eigens.hh strand-fragments.hh \
	sfcalc-genmap.hh cablam-markup.hh slurp-map.hh merge-molecules.hh arc-info.hh exp-fit.hh \
	coot_shiftfield.h shiftfield.h fib-sphere.hh cis-peptide-info.hh \
	pepflip-using-difference-map.hh blob-line.hh shape-types.hh shapes.hh \
        polar-atoms.hh 3d-texture.hh fragment-container.hh map-index.hh \
        atom-tree.hh jed-flip.hh atom-vertex.hh merge-C-and-N-terminii.hh simple-mesh.hh vertex.hh g_triangle.hh \
	cylinder.hh cylinder-with-rotation-translation.hh oct.hh prideout-octasphere.hh cmtz-interface.hh \
        mtz-column-auto-read.hh pugixml.hpp pugiconfig.hpp find-water-baddies.hh voidp-buffer.hh parse-pisa.hh \
        diff-diff-map-peaks.hh tiny_gltf.h json.hpp \
	stb_image.h stb_image_write.h glyco-tree.hh positron.hh mini-texture.hh texture-as-floats.hh pae.hh q-score.hh \
	acedrg-types-for-residue.hh gl-matrix.h ortep.hh coot-align.hh grid-balls.hh dict-link-info.hh \
        dirichlet-process.hh cfc.hh cremer-pople.hh born.hh inchikey-store.hh cylinder-utils.hh gmm.hh nef-utils.hh ptm-database.hh


pkgincludedir = $(includedir)/coot/coot-utils

AM_CPPFLAGS = \
	-I$(top_srcdir)	\
	$(ENHANCED_LIGAND_TOOLS_CPPFLAGS) \
	$(CLIPPER_CXXFLAGS)  \
	$(GEMMI_CXXFLAGS)    \
	$(MMDB_CXXFLAGS)     \
        $(CAIRO_CFLAGS)      \
	$(gsl_CFLAGS)	     \
	$(GLM_CXXFLAGS)	     \
	$(WITH_GUILE_CFLAGS) \
	$(RDKIT_CXXFLAGS)    \
	$(BOOST_CPPFLAGS)    \
	$(WITH_PYTHON_CFLAGS)

AM_CXXFLAGS = $(COOT_WINDOWS_CFLAGS) \
	-DPKGDATADIR='"$(pkgdatadir)"'

if USE_BACKWARD
   BACKWARD_FLAGS=-DUSE_BACKWARD
endif

# for pure coordinates (mmdb) things
# - note to self: don't some of these include usage of clipper::Coord_orth - which
# means that we should add $CLIPPER_LIBS to libcoot_coord_utils_la_LDADD
libcoot_coord_utils_la_SOURCES = coot-coord-utils.cc \
	coot-coord-utils-nucleotides.cc coot-rama.cc coot-h-bonds.cc hole.cc read-sm-cif.cc \
	bonded-pairs.cc helix-analysis.cc elastic.cc glyco-torsions.cc coot-least-squares.cc \
	atom-overlaps.cc reduce.cc atom-selection-container.cc contact-info.cc coot-shelx-ins.cc \
	secondary-structure-headers.cc plane-utils.cc c-beta-deviations.cc lidia-core-functions.cc \
	atom-tools.cc bonded-atoms.cc bricks.cc helix-like.cc merge-atom-selections.cc \
	contacts-by-bricks.cc stack-and-pair.cc fast-eigens.cc strand-fragments.cc sfcalc-genmap.cc \
	cablam-markup.cc coot_shiftfield.cpp shiftfield.cpp fib-sphere.cc polar-atoms.cc \
	3d-texture.cc merge-molecules.cc fragment-container.cc mutate.cc ptm-database.cc \
	cylinder.cc cylinder-with-rotation-translation.cc oct.cc prideout-octasphere.cc simple-mesh.cc \
	shapes.cc cmtz-interface.cc pugixml.cpp parse-pisa.cc gltf-export.cc pae.cc \
	gl-matrix.cc gl-matrix-eigensystem.cc coot-align.cc grid-balls.cc dict-link-info.cc \
	cremer-pople.cc inchikey-store.cc cylinder-utils.cc nef-utils.cc gaussian-atom-map-for-mask.cc

#    $(top_builddir)/mini-mol/libcoot-mini-mol.la
#    $(top_builddir)/utils/libcoot-utils.la
#    $(CLIPPER_LIBS) $(MMDB_LIBS)

# this needs to be linked with clipper the (ubuntu 18.04) linker tells me - what uses clipper
# in the source files for libcoot-coord-utils?
libcoot_coord_utils_la_LIBADD = \
    $(top_builddir)/lidia-core/libcoot-lidia-core.la \
    $(top_builddir)/geometry/libcoot-geometry.la     \
    $(top_builddir)/mini-mol/libcoot-mini-mol.la     \
    $(top_builddir)/utils/libcoot-utils.la           \
    $(top_builddir)/mini-mol/libcoot-mini-mol.la     \
    $(CLIPPER_LIBS) $(GEMMI_LIBS) $(MMDB_LIBS) \
    $(gsl_LIBS) $(RDKIT_LIBS) $(CAIRO_LIBS) $(BOOST_LDFLAGS) -lpthread

libcoot_coord_utils_la_LDFLAGS = $(SHARED_LDFLAGS)

# shelx
# libcoot_shelx_la_SOURCES = coot-shelx-ins.cc

# clean
# libcoot_shelx_la_LIBADD =  \
# libcoot-coord-utils.la \
# $(top_builddir)/utils/libcoot-utils.la

# ibcoot_shelx_la_LDFLAGS = $(SHARED_LDFLAGS)

# for non-simple coordinates things
libcoot_coord_extras_la_SOURCES = coot-coord-extras.cc coot-coord-lsq.cc \
	coot-tree-extras.cc water-coordination.cc lsq-improve.cc coot-coord-utils-glyco.cc jed-flip.cc cfc.cc

# needs more, mini-mol, geometry
#     $(top_builddir)/mini-mol/libcoot-mini-mol.la
#    $(top_builddir)/utils/libcoot-utils.la

# should libcoot-coord-utils link $CLIPPER_LIBS? Yes, I think so.
libcoot_coord_extras_la_LIBADD = \
   libcoot-coord-utils.la \
   $(top_builddir)/lidia-core/libcoot-lidia-core.la \
   $(top_builddir)/geometry/libcoot-geometry.la \
   $(top_builddir)/mini-mol/libcoot-mini-mol.la \
   $(top_builddir)/utils/libcoot-utils.la \
   $(CLIPPER_LIBS)

libcoot_coord_extras_la_LDFLAGS = $(SHARED_LDFLAGS)


# for (map) things that need to link mmdb and clipper:
libcoot_map_utils_la_SOURCES = coot-map-utils.cc peak-search.cc peak-search-from-list.cc trim.cc \
	xmap-stats.cc maps-spherical.cc pepflip-using-difference-map.cc exp-fit.cc blob-line.cc \
        slurp-map.cc find-water-baddies.cc diff-diff-map-peaks.cc mini-texture.cc positron.cc \
        texture-as-floats.cc merge-C-and-N-terminii.cc born.cc


# should libcoot-coord-utils link $CLIPPER_LIBS directly? (Hence CLIPPER_LIBS is not needed here)
libcoot_map_utils_la_LIBADD = \
   libcoot-coord-utils.la \
   $(top_builddir)/lidia-core/libcoot-lidia-core.la \
   $(top_builddir)/protein_db/libcoot-protein-db.la \
   $(top_builddir)/geometry/libcoot-geometry.la \
   $(top_builddir)/utils/libcoot-utils.la \
   $(CLIPPER_LIBS) $(gsl_LIBS) -lpthread -lz

libcoot_map_utils_la_LDFLAGS = \
   $(SHARED_LDFLAGS) -lpthread

# for map things that need mmdb, clipper, gsl and other things
libcoot_map_heavy_la_SOURCES = coot-map-heavy.cc coot-fffear.cc spherical-harmonics.cc segmap.cc

#   $(top_builddir)/utils/libcoot-utils.la

# and mini-mol (just one).
libcoot_map_heavy_la_LIBADD = \
   libcoot-map-utils.la \
   libcoot-coord-utils.la \
   $(top_builddir)/lidia-core/libcoot-lidia-core.la \
   $(top_builddir)/mini-mol/libcoot-mini-mol.la \
   $(top_builddir)/utils/libcoot-utils.la \
   $(CLIPPER_LIBS) $(gsl_LIBS) -lpthread

libcoot_map_heavy_la_LDFLAGS = $(SHARED_LDFLAGS)

# lidia-core functions moved here because they use coot-utils functions
#
check_PROGRAMS = test-h-bonds testshelxreader test-hole test-elastic \
		 test-spherical glyco-torsions test-tt-morph test-coot-probe test-chem-feat \
		 test-cod-atom-types dictionary-atom-types cod-db-access \
		 test-5-residue-fragments test-diff-diff-map-peaks \
		 make-cablam-like-stats test-cod-atom-types-from-sm-cif test-parse-pisa \
		 test-coot-utils test-pae servalcat-tracker test-zd test-grid-balls \
		 test-dirichlet-process test-born test-gmm test-nef-utils

# wrong build order - delete from make check for now.
# test-emma test-convolute test-coot-utils test-segmap

bin_PROGRAMS = coot-mmrrcc

# test_zo_chi_table_SOURCES = test-zo-chi-table.cc

# test_agg_clust_SOURCES = test-agg-clust.cc
# test_agg_clust_LDADD =     $(top_builddir)/utils/libcoot-utils.la $(CLIPPER_LIBS)


#	$(top_builddir)/geometry/libcoot-geometry.la
#        $(top_builddir)/mini-mol/libcoot-mini-mol.la
#        $(top_builddir)/utils/libcoot-utils.la $(GLOB_LIBS)


# Why do I need to add lidia-core here? It's already added to libcoot-coord-utils.
coot_mmrrcc_SOURCES = mmrrcc.cc
coot_mmrrcc_LDADD = ./libcoot-map-utils.la \
	./libcoot-coord-utils.la                         \
        $(top_builddir)/lidia-core/libcoot-lidia-core.la \
	$(top_builddir)/geometry/libcoot-geometry.la     \
	$(top_builddir)/lidia-core/libcoot-lidia-core.la \
	$(top_builddir)/utils/libcoot-utils.la           \
	$(CLIPPER_LIBS)

make_cablam_like_stats_SOURCES = make-cablam-like-stats.cc

make_cablam_like_stats_LDADD =   \
        ./libcoot-coord-utils.la \
	./libcoot-map-utils.la   \
	$(top_builddir)/geometry/libcoot-geometry.la \
	$(top_builddir)/utils/libcoot-utils.la \
	$(CLIPPER_LIBS)	$(MMDB_LIBS)

test_coot_probe_SOURCES = test-coot-probe.cc
test_coot_probe_LDADD =          \
	./libcoot-coord-utils.la \
	$(top_builddir)/geometry/libcoot-geometry.la \
	$(MMDB_LIBS)

servalcat_tracker_SOURCES = servalcat-tracker.cc

# servalcat_tracker_LDADD =

test_parse_pisa_SOURCES = test-parse-pisa.cc
test_parse_pisa_LDADD = \
	./libcoot-coord-utils.la

test_pae_SOURCES = test-pae.cc
test_pae_LDADD = \
	./libcoot-coord-utils.la

test_born_SOURCES = test-born.cc
test_born_LDADD = \
	./libcoot-map-utils.la \
	$(CLIPPER_LIBS)

test_nef_utils_SOURCES = test-nef-utils.cc
test_nef_utils_LDADD = ./libcoot-coord-utils.la

test_gmm_SOURCES = test-gmm.cc
test_gmm_LDADD =

test_dirichlet_process_SOURCES = test-dirichlet-process.cc
test_dirichlet_process_LDADD = $(top_builddir)/utils/libcoot-utils.la

testshelxreader_SOURCES = test-shelx-reader.cc
testshelxreader_LDADD =  ./libcoot-coord-utils.la    \
	$(top_builddir)/geometry/libcoot-geometry.la

# if OOPS_WRONG_BUILD_ORDER # - No - I can't use coords libraries - they have not been built yet - fix later.
# test_emma_SOURCES = test-emma.cc
# test_emma_LDADD =  libcoot-map-utils.la \
# 	./libcoot-coord-utils.la  \
# 	$(top_builddir)/coords/libcoot-coords.la \
# 	$(CLIPPER_LIBS) $(FFTW2_LIBS)

# test_convolute_SOURCES = test-convolute.cc
# test_convolute_LDADD = \
# 	./libcoot-coord-utils.la  \
# 	./libcoot-map-utils.la \
# 	./libcoot-map-heavy.la \
# 	$(top_builddir)/coords/libcoot-coords.la \
# 	$(top_builddir)/utils/libcoot-utils.la \
# 	$(CLIPPER_LIBS)

# # other test binaries may need to be linked like this too
test_coot_utils_SOURCES = test-coot-utils.cc
test_coot_utils_LDADD = \
	./libcoot-map-heavy.la    \
	./libcoot-map-utils.la    \
	./libcoot-coord-extras.la \
	./libcoot-coord-utils.la  \
	$(top_builddir)/geometry/libcoot-geometry.la \
	$(top_builddir)/utils/libcoot-utils.la \
	$(CLIPPER_LIBS) $(CCP4SRS_LIBS) -lpthread

# test_segmap_SOURCES = test-segmap.cc
# test_segmap_LDADD = \
# 	./libcoot-map-utils.la \
# 	./libcoot-map-heavy.la \
# 	$(top_builddir)/coords/libcoot-coords.la \
# 	$(CLIPPER_LIBS) -lpthread


test_h_bonds_SOURCES = test-h-bonds.cc
test_h_bonds_LDADD = \
	 ./libcoot-coord-utils.la \
	$(top_builddir)/geometry/libcoot-geometry.la \
	$(top_builddir)/utils/libcoot-utils.la \
	$(GLOB_LIBS) $(MMDB_LIBS)

test_hole_SOURCES = test-hole.cc
test_hole_LDADD = ./libcoot-coord-utils.la \
	$(top_builddir)/geometry/libcoot-geometry.la \
	$(MMDB_LIBS)

test_tt_morph_SOURCES = test-morph.cc
test_tt_morph_LDADD = ./libcoot-map-heavy.la \
	./libcoot-map-utils.la \
	./libcoot-coord-utils.la \
	$(top_builddir)/utils/libcoot-utils.la \
	$(top_builddir)/geometry/libcoot-geometry.la \
	$(CLIPPER_LIBS)

test_elastic_SOURCES = test-elastic.cc
test_elastic_LDADD   = ./libcoot-coord-utils.la      \
	$(top_builddir)/geometry/libcoot-geometry.la \
	$(GLOB_LIBS) $(MMDB_LIBS)

test_spherical_SOURCES = test-spherical.cc
test_spherical_LDADD   = ./libcoot-coord-utils.la $(CLIPPER_LIBS) $(gsl_LIBS)

glyco_torsions_SOURCES = test-glyco-torsions.cc
glyco_torsions_LDADD   = ./libcoot-coord-utils.la    \
	$(top_builddir)/geometry/libcoot-geometry.la \
	$(top_builddir)/utils/libcoot-utils.la       \
	$(MMDB_LIBS)

test_grid_balls_SOURCES = test-grid-balls.cc
test_grid_balls_LDADD   = ./libcoot-coord-utils.la    \
	$(CLIPPER_LIBS) $(MMDB_LIBS)


test_chem_feat_SOURCES = test-chem-feat.cc

test_chem_feat_LDADD = \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/geometry/libcoot-geometry.la \
	$(top_builddir)/lidia-core/libcoot-lidia-core.la \
	$(MMDB_LIBS) $(RDKIT_LIBS)

test_zd_SOURCES = test-zd.cc

test_zd_LDADD = libcoot-map-utils.la $(CLIPPER_LIBS)

test_cod_atom_types_SOURCES = test-cod-atom-types.cc

test_cod_atom_types_LDADD = \
	libcoot-coord-utils.la \
	$(top_builddir)/geometry/libcoot-geometry.la \
	$(top_builddir)/lidia-core/libcoot-lidia-core.la \
	$(top_builddir)/utils/libcoot-utils.la \
	$(RDKIT_LIBS)

test_cod_atom_types_from_sm_cif_SOURCES = test-cod-atom-types-from-sm-cif.cc

test_cod_atom_types_from_sm_cif_LDADD = \
	libcoot-coord-utils.la \
	$(top_builddir)/geometry/libcoot-geometry.la \
	$(top_builddir)/lidia-core/libcoot-lidia-core.la \
	$(top_builddir)/utils/libcoot-utils.la \
	$(CLIPPER_LIBS) $(RDKIT_LIBS)

cod_db_access_SOURCES = cod-db-access.cc

cod_db_access_LDADD = $(top_builddir)/lidia-core/libcoot-lidia-core.la $(RDKIT_LIBS)

test_5_residue_fragments_SOURCES = test-5-residue-fragments.cc

test_5_residue_fragments_LDADD = \
	libcoot-coord-utils.la \
	$(top_builddir)/geometry/libcoot-geometry.la \
	$(top_builddir)/utils/libcoot-utils.la \
	$(CLIPPER_LIBS) \
	$(RDKIT_LIBS)

# This cannot be here - it depends on libcoot-analysis which hasn't been built yet.
# test_q_score_SOURCES = test-q-score.cc
# 
# test_q_score_LDADD = libcoot-map-utils.la \
#                      libcoot-coord-utils.la \
#                      $(top_builddir)/geometry/libcoot-geometry.la \
#                      $(top_builddir)/analysis/libcoot-analysis.la \
#                      $(top_builddir)/utils/libcoot-utils.la \
#                      $(CLIPPER_LIBS)

test_diff_diff_map_peaks_SOURCES = test-diff-diff-map-peaks.cc

test_diff_diff_map_peaks_LDADD = libcoot-map-utils.la libcoot-coord-utils.la $(CLIPPER_LIBS) $(EMBEDDED_PYTHON_LIBS)

# This cannot be here - it depends on libcoot-analysis which hasn't been built yet.
# test_map_tubes_SOURCES = test-map-tubes.cc
# 
# test_map_tubes_LDADD = libcoot-map-utils.la $(top_builddir)/analysis/libcoot-analysis.la $(CLIPPER_LIBS)

dictionary_atom_types_SOURCES = dictionary-atom-types.cc

dictionary_atom_types_LDADD = \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/lidia-core/libcoot-lidia-core.la \
	$(top_builddir)/geometry/libcoot-geometry.la \
	$(RDKIT_LIBS)
