
lib_LTLIBRARIES=libcoot-layla.la

AM_CPPFLAGS = \
	-I$(top_srcdir)         \
	$(COOT_GSL_FLAGS)	\
	$(GSL_CFLAGS)		\
	$(CLIPPER_CXXFLAGS)     \
	$(RDKIT_CXXFLAGS)       \
        $(BOOST_CPPFLAGS)       \
	$(PYTHON_CPPFLAGS)      \
        $(GTK_CFLAGS)           \
        $(GLIB_CFLAGS)          \
        $(MMDB_CXXFLAGS)

# needed? test make dist
#  EXTRA_DIST = layla.ui

subpkgdatadir = $(pkgdatadir)/ui
subpkgdata_DATA = layla.ui

pkginclude_HEADERS = ligand_builder_generators.hpp  ligand_builder_notifier.hpp  ligand_builder_state.hpp \
                     ligand_builder_ui.hpp  ligand_builder_utils.hpp  ligand_editor_canvas.hpp

pkgincludedir = $(includedir)/coot/layla

libcoot_layla_la_SOURCES = ligand_builder_generators.cpp  ligand_builder_signals.cpp  ligand_builder_ui.cpp    \
                           ligand_editor_canvas.cpp       ligand_builder_notifier.cpp ligand_builder_state.cpp \
                           ligand_builder_utils.cpp       python-utils.cpp \
		           ligand_editor_canvas/core.cpp  \
		           ligand_editor_canvas/model.cpp \
		           ligand_editor_canvas/tools.cpp

libcoot_layla_la_LIBADD = $(top_builddir)/geometry/libcoot-geometry.la \
                          $(top_builddir)/utils/libcoot-utils.la \
                          $(top_builddir)/lidia-core/libcoot-lidia-core.la \
                          $(RDKIT_LIBS) -lpthread

libcoot_layla_la_LDFLAGS = $(SHARED_LDFLAGS)

AM_CXXFLAGS = \
	-DPKGDATADIR='"$(pkgdatadir)"' \
	$(SQLITE3_CFLAGS)              \
	$(COOT_WINDOWS_CFLAGS)         \
        $(LIBPNG_CFLAGS)

# maybe libexec_PROGRAMS if/when we use a shell-script wrapper
bin_PROGRAMS = layla

layla_SOURCES = main.cpp

layla_LDADD = libcoot-layla.la $(GTK_LIBS) $(GLIB_LIBS)

