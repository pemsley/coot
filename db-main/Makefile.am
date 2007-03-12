
lib_LTLIBRARIES=libcoot-db-main.la

# why doesn't this work?
# INCLUDES = -I${top_scrdir}/mini-mol $(MMDB_CXXFLAGS) $(CLIPPER_CXXFLAGS)

AM_CXXFLAGS = \
	-DPKGDATADIR='"$(pkgdatadir)"' \
	-DDATADIR='"$(datadir)"'

INCLUDES = -I$(top_srcdir)/coot-utils -I$(top_srcdir)/mini-mol $(MMDB_CXXFLAGS) $(CLIPPER_CXXFLAGS)

include_HEADERS = db-main.hh

libcoot_db_main_la_SOURCES  = db-main.cc

noinst_PROGRAMS = testdbmain

testdbmain_SOURCES = testdbmain.cc

testdbmain_LDADD = libcoot-db-main.la \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/coot-utils/libcoot-utils.la \
	$(CLIPPER_LIBS)      \
	$(MMDB_LIBS)


EXTRA_DIST = Makefile.simple
