
lib_LTLIBRARIES=libcoot-cabuild.la

include_HEADERS = CalphaBuild.hh

INCLUDES =  \
	-I$(top_srcdir)/angles \
	$(CLIPPER_CXXFLAGS) \
	$(MMDB_CXXFLAGS)

noinst_PROGRAMS = testcabuild

testcabuild_SOURCES = testcabuild.cc

libcoot_cabuild_la_SOURCES = CalphaBuild.cc

testcabuild_LDADD = \
	./libcoot-cabuild.la \
	$(top_builddir)/angles/libcoot-AngleInfo.la \
	$(CLIPPER_LIBS) $(MMDB_LIBS)

EXTRA_DIST = build.txt
