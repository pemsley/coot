
lib_LTLIBRARIES=libcoot-analysis.la

pkginclude_HEADERS=bfkurt.hh mogul-interface.hh
pkgincludedir = $(includedir)/coot/analysis

# CXXFLAGS = 	$(COOT_GSL_FLAGS)

INCLUDES = \
	-I$(top_srcdir) \
	$(CLIPPER_CXXFLAGS)		\
	$(MMDB_CXXFLAGS)

libcoot_analysis_la_SOURCES = bfkurt.cc mogul.cc 

libcoot_analysis_la_LIBADD = \
   $(top_builddir)/coot-utils/libcoot-coord-utils.la    \
   $(top_builddir)/geometry/libcoot-geometry.la    	\
   $(top_builddir)/mini-mol/libcoot-mini-mol.la    	\
   $(top_builddir)/utils/libcoot-utils.la 		\
   $(MMDB_LIBS)  \
   $(GSL_LIBS)

bin_PROGRAMS = bfactan

check_PROGRAMS = test-mogul

bfactan_SOURCES = bfactan.cc

test_mogul_SOURCES = test-mogul.cc

test_mogul_LDADD = libcoot-analysis.la

bfactan_LDADD = \
	libcoot-analysis.la 					\
	$(top_builddir)/coot-utils/libcoot-coord-utils.la 	\
	$(top_builddir)/geometry/libcoot-geometry.la    	\
        $(top_builddir)/mini-mol/libcoot-mini-mol.la    	\
	$(top_builddir)/utils/libcoot-utils.la 		\
	$(GLOB_LIBS) 		\
	$(CLIPPER_LIBS) 	\
	$(MMDB_LIBS)

EXTRA_DIST = Makefile.simple
