
lib_LTLIBRARIES=libcoot-lidia.la

INCLUDES = \
	-I$(top_srcdir)/lidia-core \
	-I$(top_srcdir)/coot-utils \
	-I$(top_srcdir)/mini-mol \
	-I$(top_srcdir)/geometry \
	-I$(top_srcdir)/coot-utils \
	-I$(top_srcdir)/utils \
        $(GSL_CFLAGS) \
        $(GOOCANVAS_CFLAGS) \
        $(GTK_CFLAGS) \
	$(CLIPPER_CXXFLAGS) $(MMDB_CXXFLAGS) \
	$(COOT_WINDOWS_CFLAGS) \
	$(ENTERPRISE_TOOLS_CFLAGS)      \
	$(RDKIT_CXXFLAGS)

AM_CXXFLAGS = \
	-DPKGDATADIR='"$(pkgdatadir)"'


bin_PROGRAMS = lidia

libcoot_lidia_la_SOURCES = lbg.cc wmolecule.cc  \
    lbg-search.cc  residue-circles.cc \
    flev-annotations.cc

libcoot_lidia_la_LIBADD = \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/geometry/libcoot-geometry.la \
	$(top_builddir)/lidia-core/libcoot-lidia-core.la \
	$(top_builddir)/utils/libcoot-utils.la \
	$(CLIPPER_LIBS) $(MMDB_LIBS) $(GSL_LIBS) $(RDKIT_LIBS) $(GOOCANVAS_LIBS) $(GTK_LIBS)

lidia_SOURCES = lidia-main.cc lbg-callbacks-outer.cc

include_HEADERS = lbg.hh wmolecule.hh  flev-annotations.hh graphics-c-interface-functions.hh


# res spec output operator means we include libcoot-coord-utils
# (if that is annoying, then remove it and expand it by hand)
#
lidia_LDADD = \
	./libcoot-lidia.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/geometry/libcoot-geometry.la \
	$(top_builddir)/utils/libcoot-utils.la \
	$(top_builddir)/lidia-core/libcoot-lidia-core.la \
	$(CLIPPER_LIBS) $(MMDB_LIBS) $(GSL_LIBS) $(RDKIT_LIBS) $(GOOCANVAS_LIBS) $(GTK_LIBS)

EXTRA_DIST=lbg.glade $(appicon_DATA) lbg-callbacks.cc

subpkgdatadir = $(pkgdatadir)
subpkgdata_DATA = lbg.glade lbg-callbacks.cc

check_PROGRAMS = test_lbg_functions
test_lbg_functions_SOURCES = test_lbg_functions.cc  graphics-c-interface-functions-blanks.cc
test_lbg_functions_LDADD = \
       ./libcoot-lidia.la \
	$(top_builddir)/lidia-core/libcoot-lidia-core.la \
	$(top_builddir)/geometry/libcoot-geometry.la \
	$(top_builddir)/utils/libcoot-utils.la \
	$(RDKIT_LIBS)




appicondir=$(pkgdatadir)/pixmaps
appicon_DATA = \
	3-c-new.svg    6-c-shadow-experiment.svg	   single.svg                           \
	3-c-other.svg  6-c.svg				   sp3.svg                              \
	3-c.svg        7-c.svg				   stereo-out.svg                       \
	4-c.svg        8-c.svg				   TheStructorr_magnifying_glass-32.svg \
	5-c.svg        Anonymous_Architetto_--_Chiave.svg  TheStructorr_magnifying_glass.svg    \
	6-c-arom.svg   charge.svg			   triple.svg                           \
	6-c-exp-2.svg  double.svg Trash-hydrogens-Andy_Trash_Can.svg

## #else 

##libcoot_lidia_la_SOURCES = flev-annotations.cc
## include_HEADERS = flev-annotations.hh

## endif
