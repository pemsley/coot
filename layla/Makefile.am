
AUTOMAKE_OPTIONS = subdir-objects

lib_LTLIBRARIES=libcoot-layla.la

AM_CPPFLAGS = \
	-I$(top_srcdir)         \
	$(COOT_GSL_FLAGS)	\
	$(GSL_CFLAGS)		\
	$(CLIPPER_CXXFLAGS)     \
	$(RDKIT_CXXFLAGS)       \
        $(BOOST_CPPFLAGS)       \
	-DUSE_PYTHON=1 $(PYTHON_CPPFLAGS)      \
	$(GTK_CFLAGS)           \
	$(GLIB_CFLAGS)          \
	-DUSE_LIBCURL=1 $(CURL_CFLAGS)          \
	$(MMDB_CXXFLAGS)

# needed? test make dist
#  EXTRA_DIST = layla.ui

subpkgdatadir = $(pkgdatadir)/ui
subpkgdata_DATA = layla.ui

pkginclude_HEADERS = \
	layla_embedded.hpp \
	generators.hpp \
	notifier.hpp \
	ui.hpp \
	utils.hpp \
	state.hpp \
	ligand_editor_canvas/core.hpp  \
	ligand_editor_canvas/model.hpp \
	ligand_editor_canvas/tools.hpp \
	ligand_editor_canvas/render.hpp \
	ligand_editor_canvas.hpp \
	python_utils.hpp

pkgincludedir = $(includedir)/coot/layla

libcoot_layla_la_SOURCES = \
	layla_embedded.cpp \
	generators.cpp \
	notifier.cpp \
	signals.cpp \
	state.cpp \
	ui.cpp \
	utils.cpp \
	ligand_editor_canvas/core.cpp \
	ligand_editor_canvas.cpp \
	ligand_editor_canvas/model.cpp \
	ligand_editor_canvas/tools.cpp \
	ligand_editor_canvas/render.cpp \
	python_utils.cpp

libcoot_layla_la_LIBADD = $(top_builddir)/geometry/libcoot-geometry.la \
                          $(top_builddir)/lidia-core/libcoot-lidia-core.la \
                          $(top_builddir)/utils/libcoot-utils.la \
                          $(GTK_LIBS) $(GLIB_LIBS) \
                          $(RDKIT_LIBS) \
                          -lpthread

libcoot_layla_la_LDFLAGS = $(SHARED_LDFLAGS)

EXTRA_DIST = layla.ui

AM_CXXFLAGS = \
	-DPKGDATADIR='"$(pkgdatadir)"' \
	$(SQLITE3_CFLAGS)              \
	$(COOT_WINDOWS_CFLAGS)         \
	$(CURL_LIBS)                   \
        $(LIBPNG_CFLAGS)

# maybe libexec_PROGRAMS if/when we use a shell-script wrapper
bin_PROGRAMS = layla

layla_SOURCES = main.cpp

layla_LDADD = libcoot-layla.la $(top_builddir)/utils/libcoot-utils.la $(GTK_LIBS) $(GLIB_LIBS)

