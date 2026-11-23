# locres/Makefile.am
#
# libcoot-locres uses clipper and nothing else from Coot, so that it can be
# built and used without the rest of Coot. fsc-core.hh is also used by
# coot-utils, which includes it as "locres/fsc-core.hh".

lib_LTLIBRARIES=libcoot-locres.la

pkginclude_HEADERS = coot-locres.hh fsc-core.hh
pkgincludedir = $(includedir)/coot/locres

AM_CPPFLAGS = $(CLIPPER_CXXFLAGS)

AM_CXXFLAGS = $(COOT_WINDOWS_CFLAGS)

libcoot_locres_la_SOURCES = coot-locres.cc

libcoot_locres_la_LIBADD = $(CLIPPER_LIBS)

libcoot_locres_la_LDFLAGS = $(SHARED_LDFLAGS)

bin_PROGRAMS = coot-locres

coot_locres_SOURCES = coot-locres-main.cc

coot_locres_LDADD = libcoot-locres.la $(CLIPPER_LIBS)

check_PROGRAMS = test-locres

test_locres_SOURCES = test-locres.cc

test_locres_LDADD = libcoot-locres.la $(CLIPPER_LIBS)
