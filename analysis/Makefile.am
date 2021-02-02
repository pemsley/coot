
lib_LTLIBRARIES=libcoot-analysis.la

pkginclude_HEADERS=bfkurt.hh mogul-interface.hh stats.hh
pkgincludedir = $(includedir)/coot/analysis

# CXXFLAGS = 	$(COOT_GSL_FLAGS)

AM_CPPFLAGS = \
	-I$(top_srcdir) \
	$(CLIPPER_CXXFLAGS)		\
	$(COOT_GSL_FLAGS) $(GSL_CFLAGS)	\
	$(MMDB_CXXFLAGS)

libcoot_analysis_la_SOURCES = bfkurt.cc mogul.cc kolmogorov.hh kolmogorov.cc stats.cc \
	cablam.hh cablam.cc typed-distances.cc typed-distances.hh \
	b-factor-histogram.cc b-factor-histogram.hh daca.cc daca.hh

# note to self: does RDKIT_LIBS belong here or should it be added to the LDADD of
# one of these libraries?
libcoot_analysis_la_LIBADD = \
   $(top_builddir)/coot-utils/libcoot-coord-utils.la    \
   $(top_builddir)/geometry/libcoot-geometry.la    	\
   $(top_builddir)/mini-mol/libcoot-mini-mol.la    	\
   $(top_builddir)/utils/libcoot-utils.la 		\
   $(CLIPPER_LIBS) $(MMDB_LIBS)  \
   $(GSL_LIBS) $(RDKIT_LIBS)

libcoot_analysis_la_LDFLAGS = $(SHARED_LDFLAGS)

bin_PROGRAMS = coot-bfactan

check_PROGRAMS = test-mogul test-cablam peptide-projection test-align improper-dihedrals \
	test-stats coot-daca

coot_daca_SOURCES = coot-daca.cc

# what is goin on here!? 	$(CLIPPER_LIBS) is needed to resolve libRDKitPartialCharges.so.1
# at run time! I so don't understand!
coot_daca_LDADD = libcoot-analysis.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la 	\
	$(top_builddir)/geometry/libcoot-geometry.la    	\
	$(top_builddir)/utils/libcoot-utils.la 		\
	$(CLIPPER_LIBS)

coot_bfactan_SOURCES = bfactan.cc

test_mogul_SOURCES = test-mogul.cc

test_mogul_LDADD = libcoot-analysis.la

test_cablam_SOURCES = test-cablam.cc

test_stats_SOURCES = test-stats.cc

test_align_SOURCES = test-align.cc

test_cablam_LDADD = \
	libcoot-analysis.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la	\
	$(top_builddir)/coords/libcoot-coords.la \
	$(MMDB_LIBS)

peptide_projection_SOURCES  = peptide-projection.cc

peptide_projection_LDADD = libcoot-analysis.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la	\
	$(CLIPPER_LIBS) \
	$(MMDB_LIBS)

improper_dihedrals_SOURCES = improper-dihedrals.cc

improper_dihedrals_LDADD = \
	libcoot-analysis.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la	\
	$(top_builddir)/coords/libcoot-coords.la \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la	\
	$(top_builddir)/utils/libcoot-utils.la		\
	$(MMDB_LIBS)

test_align_LDADD = libcoot-analysis.la \
	$(top_builddir)/utils/libcoot-utils.la

coot_bfactan_LDADD = \
	libcoot-analysis.la					\
	$(top_builddir)/coot-utils/libcoot-coord-utils.la 	\
	$(top_builddir)/lidia-core/libcoot-lidia-core.la    	\
	$(top_builddir)/geometry/libcoot-geometry.la    	\
	$(top_builddir)/mini-mol/libcoot-mini-mol.la    	\
	$(top_builddir)/utils/libcoot-utils.la 		\
	$(GLOB_LIBS) 		\
	$(CLIPPER_LIBS) 	\
	$(MMDB_LIBS)

EXTRA_DIST = Makefile.simple
