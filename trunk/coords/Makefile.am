
lib_LTLIBRARIES=libcoot-coords.la

include_HEADERS =  Bond_lines.h Cartesian.h cos-sin.h \
	mmdb.h mmdb-extras.h utils.h mmdb-crystal.h Bond_lines_ext.h \
	coot-close.hh

INCLUDES = \
	-I$(top_builddir)/skeleton \
	-I$(top_srcdir)/coot-utils \
	-I$(top_srcdir)/mini-mol \
	$(CLIPPER_CXXFLAGS) $(MMDB_CXXFLAGS)

AM_CXXFLAGS = \
	-DPKGDATADIR='"$(pkgdatadir)"' \
	-DDATADIR='"$(datadir)"'    

libcoot_coords_la_SOURCES =  Bond_lines.cc Cartesian.cc cos-sin.cc \
		 mmdb-extras.cc utils.cc mmdb.cc mmdb-crystal.cc \
		 Bond_lines_ext.cc coot-close.cc

