# ideal/Makefile.am
# 
# Copyright 2004, 2005, 2006 The University of York
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



# should be set by autoconf
# GSL_LIBS = -lgsl -lgslcblas -lm

lib_LTLIBRARIES=libcoot-ideal.la

EXTRA_DIST = notes Makefile.simple molmantest.cc test-indexing.hh test-indexing.cc

pkginclude_HEADERS =  test-indexing.hh simple-restraint.hh pepflip.hh \
	regularize-minimol.hh lograma.h zo-rama.hh phi-psi.hh process-dfs-in-range.hh \
	parallel-planes.hh torsion-bonds.hh extra-restraints.hh model-bond-deltas.hh \
        refinement-lights.hh refinement-results-t.hh rama-types.hh new-linked-residue-t.hh

pkgincludedir = $(includedir)/coot/ideal

AM_CPPFLAGS = \
	-I$(top_srcdir) \
	$(CLIPPER_CXXFLAGS)        \
	$(MMDB_CXXFLAGS)

AM_CXXFLAGS = $(COOT_GSL_FLAGS) $(GSL_CFLAGS) $(COOT_WINDOWS_CFLAGS) \
     $(CLIPPER_CXXFLAGS) $(MMDB_CXXFLAGS) $(COOT_GSL_FLAGS) 
AM_LDFLAGS = $(MMDB_LIBS) $(CLIPPER_LIBS) $(GLS_LIBS)


#libcoot_ideal_la_LIBADD = ../coot-utils/libcoot-utils.la                       \
#	../coot-utils/libcoot-coord-utils.la                                   \
#	../coot-utils/libcoot-coord-extras.la ../mini-mol/libcoot-mini-mol.la  \
#	../geometry/libcoot-geometry.la

# Tim Fenn suggestion to clear up the undefined-non-weak-symbols.
#
# From the book, libxxx_LDFLAGS is for flags passed to the linker
#                libxxx_LDADD extra flags passed to the linker
#                libxxx_LIBADD as above, but for the static lib target
# 
# However, for me, it doesn't change the unresolved symbols (it does
# change dependency_libs in the .la file though).
# 
# libcoot_ideal_la_LDFLAGS = $(GSL_LIBS) $(MMDB_LIBS) $(CLIPPER_LIBS)
# reworking library order  $(top_builddir)/coords/libcoot-coords.la           

libcoot_ideal_la_LIBADD =                 \
   $(top_builddir)/coot-utils/libcoot-map-heavy.la     \
   $(top_builddir)/coot-utils/libcoot-map-utils.la     \
   $(top_builddir)/geometry/libcoot-geometry.la        \
   $(top_builddir)/coot-utils/libcoot-coord-extras.la  \
   $(top_builddir)/coot-utils/libcoot-coord-utils.la   \
   $(top_builddir)/utils/libcoot-utils.la              \
   $(top_builddir)/mini-mol/libcoot-mini-mol.la        \
   $(top_builddir)/lidia-core/libcoot-lidia-core.la    \
   $(MMDB_LIBS) $(CLIPPER_LIBS) $(GSL_LIBS) $(LIBPNG_LIBS)


libcoot_ideal_la_SOURCES = simple-restraint.cc make-restraints.cc ng.cc chi-squareds.cc \
	mods.cc distortion.cc chirals.cc flanking.cc gradients.cc extra-restraints.cc \
	pepflip.cc regularize-minimol.cc lograma.cpp process-dfs-in-range.cc \
	crankshaft.cc extra-restraints-kk.cc simple-restraint-kk.cc torsion-bonds.cc \
	coot-compare-residues.hh coot-compare-residues.cc link-restraints.cc \
	parallel-planes.cc pull-restraint.cc trans-peptide.cc model-bond-deltas.cc \
	write-png.cc zo-rama.cc crankshaft.hh phi-psi.cc neutron.cc

libcoot_ideal_la_LDFLAGS = $(SHARED_LDFLAGS) -lpthread

check_PROGRAMS =  with-geometry logramatest crankshaft-test test-rama-plot
if OS_WIN32
bin_PROGRAMS = mini-rsr-bin
else
libexec_PROGRAMS = mini-rsr-bin coot-ligand-validation-bin
endif

logramatest_SOURCES = logramatest.cpp

crankshaft_test_SOURCES = crankshaft-test.cc

crankshaft_test_LDADD = ./libcoot-ideal.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/coot-utils/libcoot-map-utils.la   \
	$(top_builddir)/geometry/libcoot-geometry.la      \
	$(top_builddir)/utils/libcoot-utils.la            \
        $(top_builddir)/compat/libcoot-compat.la -lpthread

with_geometry_SOURCES = with-geometry.cc

mini_rsr_bin_SOURCES = mini-rsr.cc

test_rama_plot_SOURCES = test-rama-plot.cc

test_rama_plot_LDADD = ./libcoot-ideal.la

coot_ligand_validation_bin_SOURCES = lv.cc

coot_ligand_validation_bin_LDADD = ./libcoot-ideal.la      \
	$(top_builddir)/coot-utils/libcoot-map-heavy.la    \
	$(top_builddir)/coot-utils/libcoot-map-utils.la    \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
        $(top_builddir)/lidia-core/libcoot-lidia-core.la   \
	$(top_builddir)/geometry/libcoot-geometry.la       \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
	$(top_builddir)/utils/libcoot-utils.la             \
        $(top_builddir)/compat/libcoot-compat.la -lpthread

# test_indexing_SOURCES = test-indexing.cc

# note that we need libcoot-coords only because we use 
# atom_selection_container_t asc get_atom_selection(filename).  If we get rid of that 
# in the future the dependency on libcoot-coords goes away.
#
with_geometry_LDADD = \
        ./libcoot-ideal.la                                 \
	$(top_builddir)/coot-utils/libcoot-map-heavy.la    \
	$(top_builddir)/geometry/libcoot-geometry.la       \
	$(top_builddir)/coords/libcoot-coords.la           \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
	$(top_builddir)/utils/libcoot-utils.la             \
        $(top_builddir)/lidia-core/libcoot-lidia-core.la   \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
	$(GLOB_LIBS) \
	$(RDKIT_LIBS) \
        $(CLIPPER_LIBS) $(MMDB_LIBS) $(GSL_LIBS)

# add libcoot-map-utils for  the debugging density_around_point function
#
# 20191226-PE I have added -lpthread to mini_rsr linking and coot-ideal
#    linking - maybe that's a mistake - there *should* only be one (but
#    I tried adding it to just the library, but mini-rsr didn't link).
#    Not sure either that this is the canonical way of adding threads.
#
mini_rsr_bin_LDADD =  \
        ./libcoot-ideal.la                                 \
	$(top_builddir)/coot-utils/libcoot-map-heavy.la    \
	$(top_builddir)/geometry/libcoot-geometry.la       \
	$(top_builddir)/coot-utils/libcoot-map-utils.la   \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
        $(top_builddir)/lidia-core/libcoot-lidia-core.la  \
	$(top_builddir)/utils/libcoot-utils.la        \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
        $(top_builddir)/compat/libcoot-compat.la           \
	$(GLOB_LIBS) \
	$(RDKIT_LIBS) \
        $(CLIPPER_LIBS) $(MMDB_LIBS) $(GSL_LIBS) -lpthread

#test_indexing_LDADD = \
#        ./libcoot-ideal.la                                \
#	$(top_builddir)/geometry/libcoot-geometry.la      \
#	$(top_builddir)/mini-mol/libcoot-mini-mol.la      \
#	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
#	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
#	$(top_builddir)/utils/libcoot-utils.la       \
#	$(CLIPPER_LIBS) $(MMDB_LIBS) $(GSL_LIBS)

# Why do I need all these libraries to simply link a trivial test program?
# Badness....
logramatest_LDADD = \
	./libcoot-ideal.la \
	$(top_builddir)/coot-utils/libcoot-map-heavy.la    \
	$(top_builddir)/geometry/libcoot-geometry.la      \
	$(top_builddir)/coords/libcoot-coords.la           \
        $(top_builddir)/lidia-core/libcoot-lidia-core.la  \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la      \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/utils/libcoot-utils.la       \
	$(RDKIT_LIBS) \
	$(CLIPPER_LIBS) $(GSL_LIBS)

dist_bin_SCRIPTS = coot-mini-rsr coot-ligand-validation

