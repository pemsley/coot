
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
	b-factor-histogram.cc b-factor-histogram.hh

libcoot_analysis_la_LIBADD = \
   $(top_builddir)/coot-utils/libcoot-coord-utils.la    \
   $(top_builddir)/geometry/libcoot-geometry.la    	\
   $(top_builddir)/mini-mol/libcoot-mini-mol.la    	\
   $(top_builddir)/utils/libcoot-utils.la 		\
   $(MMDB_LIBS)  \
   $(GSL_LIBS)

libcoot_analysis_la_LDFLAGS = $(SHARED_LDFLAGS)

bin_PROGRAMS = coot-bfactan

check_PROGRAMS = test-mogul test-cablam peptide-projection test-align improper-dihedrals test-stats

coot_bfactan_SOURCES = bfactan.cc

test_mogul_SOURCES = test-mogul.cc

test_mogul_LDADD = libcoot-analysis.la

test_cablam_SOURCES = test-cablam.cc

test_stats_SOURCES = test-stats.cc

test_align_SOURCES = test-align.cc

test_cablam_LDADD =    $(top_builddir)/coords/libcoot-coords.la    \
                       libcoot-analysis.la

peptide_projection_SOURCES  = peptide-projection.cc

peptide_projection_LDADD = libcoot-analysis.la

improper_dihedrals_SOURCES = improper-dihedrals.cc

improper_dihedrals_LDADD = $(top_builddir)/coords/libcoot-coords.la    \
			libcoot-analysis.la

test_align_LDADD = libcoot-analysis.la

coot_bfactan_LDADD = \
	libcoot-analysis.la 					\
	$(top_builddir)/coot-utils/libcoot-coord-utils.la 	\
	$(top_builddir)/geometry/libcoot-geometry.la    	\
        $(top_builddir)/mini-mol/libcoot-mini-mol.la    	\
	$(top_builddir)/utils/libcoot-utils.la 		\
	$(GLOB_LIBS) 		\
	$(CLIPPER_LIBS) 	\
	$(MMDB_LIBS)

EXTRA_DIST = Makefile.simple
