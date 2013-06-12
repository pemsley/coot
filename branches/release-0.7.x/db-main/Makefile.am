
lib_LTLIBRARIES=libcoot-db-main.la

AM_CXXFLAGS = \
	-DPKGDATADIR='"$(pkgdatadir)"' \
	-DDATADIR='"$(datadir)"'

INCLUDES = -I$(top_srcdir)/utils -I$(top_srcdir)/coot-utils -I$(top_srcdir)/mini-mol \
   $(MMDB_CXXFLAGS) $(CLIPPER_CXXFLAGS)

include_HEADERS = db-main.hh db-strands.hh

libcoot_db_main_la_SOURCES  = db-main.cc db-strands.cc

libcoot_db_main_la_LIBADD  = \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/utils/libcoot-utils.la \
        $(CLIPPER_LIBS) $(MMDB_LIBS)

noinst_PROGRAMS = testdbmain train-pepflip

testdbmain_SOURCES = testdbmain.cc

testdbmain_LDADD = libcoot-db-main.la \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/geometry/libcoot-geometry.la \
	$(top_builddir)/utils/libcoot-utils.la \
	$(GLOB_LIBS)	     \
	$(CLIPPER_LIBS)      \
	$(MMDB_LIBS)

train_pepflip_SOURCES = train-pepflip.cc

train_pepflip_LDADD = libcoot-db-main.la \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/geometry/libcoot-geometry.la \
	$(top_builddir)/utils/libcoot-utils.la \
	$(GLOB_LIBS)	     \
	$(CLIPPER_LIBS)      \
	$(MMDB_LIBS)


EXTRA_DIST = Makefile.simple
