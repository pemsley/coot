
lib_LTLIBRARIES=libcoot-analysis.la

include_HEADERS=bfkurt.hh mogul-interface.hh

INCLUDES = \
	-I$(top_srcdir)/coot-utils   	\
	-I$(top_srcdir)/mini-mol   	\
	$(CLIPPER_CXXFLAGS)		\
	$(MMDB_CXXFLAGS)

libcoot_analysis_la_SOURCES = bfkurt.cc mogul.cc 

libcoot_analysis_la_LIBADD = \
   $(top_builddir)/coot-utils/libcoot-coord-utils.la    \
   $(top_builddir)/geometry/libcoot-geometry.la    	\
   $(top_builddir)/mini-mol/libcoot-mini-mol.la    	\
   $(top_builddir)/utils/libcoot-utils.la 		\
   $(MMDB_LIBS) 

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
	$(CLIPPER_LIBS) 	\
	$(MMDB_LIBS)

EXTRA_DIST = Makefile.simple
