
lib_LTLIBRARIES=libcoot-lidia.la

if MAKE_ENTERPRISE_TOOLS

AM_CXXFLAGS = \
	-DPKGDATADIR='"$(pkgdatadir)"'


bin_PROGRAMS = lidia

libcoot_lidia_la_SOURCES = lbg.cc wmolecule.cc lbg-callbacks.cc \
    lbg-search.cc some-coot-utils.cc residue-circles.cc enterprise.cc

lidia_SOURCES = lidia-main.cc

# some-coot-utils.hh needed?
#
include_HEADERS = lbg.hh wmolecule.hh enterprise.hh

INCLUDES = $(GSL_CFLAGS) \
        $(GOOCANVAS_CFLAGS) \
        $(GTK_CFLAGS) \
	$(CLIPPER_CXXFLAGS) $(MMDB_CXXFLAGS) \
	$(ENTERPRISE_TOOLS_CFLAGS)      \
	$(RDKIT_CXXFLAGS)


lidia_LDADD = \
	./libcoot-lidia.la \
	$(top_builddir)/lidia-core/libcoot-lidia-core.la \
	$(CLIPPER_LIBS) $(MMDB_LIBS) $(GSL_LIBS) $(RDKIT_LIBS) $(GOOCANVAS_LIBS) $(GTK_LIBS)

EXTRA_DIST=lbg.glade

subpkgdatadir = $(pkgdatadir)
subpkgdata_DATA = lbg.glade

appicondir=$(pkgdatadir)/pixmaps
appicon_DATA = \
	3-c-new.svg    6-c-shadow-experiment.svg	   single.svg                           \
	3-c-other.svg  6-c.svg				   sp3.svg                              \
	3-c.svg        7-c.svg				   stereo-out.svg                       \
	4-c.svg        8-c.svg				   TheStructorr_magnifying_glass-32.svg \
	5-c.svg        Anonymous_Architetto_--_Chiave.svg  TheStructorr_magnifying_glass.svg    \
	6-c-arom.svg   charge.svg			   triple.svg                           \
	6-c-exp-2.svg  double.svg

else 

libcoot_lidia_la_SOURCES = 

endif
