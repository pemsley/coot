
lib_LTLIBRARIES=libcoot-db-main.la

AM_CXXFLAGS = \
	-DPKGDATADIR='"$(pkgdatadir)"' \
	-DDATADIR='"$(datadir)"'

AM_CPPFLAGS = -I$(top_srcdir) \
   $(MMDB_CXXFLAGS) $(CLIPPER_CXXFLAGS) $(COOT_WINDOWS_CFLAGS)

pkginclude_HEADERS = db-main.hh db-strands.hh
pkgincludedir = $(includedir)/coot/db-main

libcoot_db_main_la_SOURCES  = db-main.cc db-strands.cc

libcoot_db_main_la_LIBADD  = \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la \
        $(top_builddir)/lidia-core/libcoot-lidia-core.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/geometry/libcoot-geometry.la \
	$(top_builddir)/utils/libcoot-utils.la \
        $(CLIPPER_LIBS) $(MMDB_LIBS)

libcoot_db_main_la_LDFLAGS = $(SHARED_LDFLAGS)

check_PROGRAMS = testdbmain train-pepflip

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
