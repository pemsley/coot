
lib_LTLIBRARIES=libcoot-utils.la libcoot-map-utils.la libcoot-coord-utils.la  libcoot-coord-extras.la libcoot-map-heavy.la libcoot-shelx.la

include_HEADERS = coot-utils.hh peak-search.hh coot-trim.hh \
	coot-coord-utils.hh coot-map-utils.hh coot-coord-extras.hh \
	coot-fasta.hh coot-map-heavy.hh coot-lsq-types.h coot-shelx.hh

INCLUDES = \
	-I$(top_srcdir)/coords  	\
	-I$(top_srcdir)/geometry	\
	$(CLIPPER_CXXFLAGS) 		\
	$(MMDB_CXXFLAGS) 		\
	$(GSL_CFLAGS)

AM_CXXFLAGS = $(COOT_GSL_FLAGS) $(COOT_WINDOWS_CFLAGS) \
	-DPKGDATADIR='"$(pkgdatadir)"'

# for simple things:
libcoot_utils_la_SOURCES = coot-utils.cc coot-fasta.cc atomic-numbers.cc

# for pure coordinates (mmdb) things
libcoot_coord_utils_la_SOURCES = coot-coord-utils.cc 

# shelx
libcoot_shelx_la_SOURCES = coot-shelx-ins.cc

# for non-simple coordinates things
libcoot_coord_extras_la_SOURCES = coot-coord-extras.cc coot-coord-lsq.cc

# for (map) things that need to link mmdb and clipper:
libcoot_map_utils_la_SOURCES = coot-map-utils.cc peak-search.cc trim.cc 

# for map things that need mmdb, clipper, gsl and other things
libcoot_map_heavy_la_SOURCES = coot-map-heavy.cc coot-fffear.cc

noinst_PROGRAMS = testshelxreader testcootutils

testshelxreader_SOURCES = test-shelx-reader.cc
testshelxreader_LDADD = \
	./libcoot-shelx.la \
	./libcoot-coord-utils.la \
	./libcoot-utils.la \
	$(CLIPPER_LIBS) $(MMDB_LIBS)

testcootutils_SOURCES = test-utils.cc
testcootutils_LDADD = \
	./libcoot-shelx.la \
	./libcoot-coord-utils.la \
	./libcoot-utils.la \
	$(CLIPPER_LIBS) $(MMDB_LIBS)


