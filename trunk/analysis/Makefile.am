
lib_LTLIBRARIES=libcoot-analysis.la

include_HEADERS=bfkurt.hh

INCLUDES = \
	-I$(top_srcdir)/coot-utils   	\
	-I$(top_srcdir)/mini-mol   	\
	$(CLIPPER_CXXFLAGS)		\
	$(MMDB_CXXFLAGS)

libcoot_analysis_la_SOURCES = bfkurt.cc

libcoot_analysis_la_LIBADD = \
   ../coot-utils/libcoot-coord-utils.la \
   $(MMDB_LIBS) 

bin_PROGRAMS = bfactan

bfactan_SOURCES = bfactan.cc

bfactan_LDADD = \
	libcoot-analysis.la 					\
	$(top_builddir)/coot-utils/libcoot-coord-utils.la 	\
	$(top_builddir)/geometry/libcoot-geometry.la    	\
	$(top_builddir)/utils/libcoot-utils.la 		\
	$(CLIPPER_LIBS) 	\
	$(MMDB_LIBS)

EXTRA_DIST = Makefile.simple
