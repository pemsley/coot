
lib_LTLIBRARIES=libcoot-geometry.la

EXTRA_DIST = read-sbase.cc read-sbase.hh notes Makefile.simple

include_HEADERS =  protein-geometry.hh

INCLUDES = -I$(top_srcdir)/coords -I$(top_srcdir)/coot-utils $(CLIPPER_CXXFLAGS) $(MMDB_CXXFLAGS)

AM_CXXFLAGS = \
	-DDATADIR='"$(datadir)"'       \
	-DPKGDATADIR='"$(pkgdatadir)"' \
	$(COOT_WINDOWS_CFLAGS)

libcoot_geometry_a_sources = protein-geometry.cc protein-geometry-mmdb.cc

libcoot_geometry_la_SOURCES = $(libcoot_geometry_a_sources)

# libcoot_geometry_la_CXXFLAGS = -g -O1

noinst_PROGRAMS = test-geometry

test_geometry_SOURCES = test-geometry.cc

test_geometry_LDADD = ./libcoot-geometry.la         \
        $(top_builddir)/coot-utils/libcoot-utils.la \
	$(CLIPPER_LIBS) $(MMDB_LIBS)
