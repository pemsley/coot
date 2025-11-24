

lib_LTLIBRARIES=libcoot-validation-graphs.la

AM_CPPFLAGS = -I$(top_srcdir) $(MMDB_CXXFLAGS) $(CLIPPER_CXXFLAGS) $(BOOST_CPPFLAGS) $(COOT_WINDOWS_CFLAGS) @GTK_CFLAGS@


pkginclude_HEADERS = \
	residue-validation-information.hh \
	validation-information.hh \
	validation-graphs.hh \
	validation-graph-widget.hh \
	sequence-view-widget.hh

pkgincludedir = $(includedir)/coot/test-validation

#    validation-information.cc \
#

libcoot_validation_graphs_la_SOURCES = \
	residue-validation-information.cc \
	validation-graph-widget.cc \
	validation-graphs.cc \
	validation-information.cc \
	sequence-view-widget.cc

libcoot_validation_graphs_la_LIBADD = \
   $(top_builddir)/ideal/libcoot-ideal.la             \
   $(top_builddir)/coot-utils/libcoot-map-heavy.la    \
   $(top_builddir)/lidia-core/libcoot-lidia-core.la   \
   $(top_builddir)/cootaneer/libcootaneer.la          \
   $(top_builddir)/db-main/libcoot-db-main.la         \
   $(top_builddir)/ccp4mg-utils/libccp4mg-utils.la    \
   $(top_builddir)/geometry/libcoot-geometry.la       \
   $(top_builddir)/mini-mol/libcoot-mini-mol.la       \
   $(top_builddir)/coot-utils/libcoot-coord-extras.la \
   $(top_builddir)/coot-utils/libcoot-coord-utils.la  \
   $(top_builddir)/coot-utils/libcoot-map-utils.la    \
   $(top_builddir)/analysis/libcoot-analysis.la	      \
   $(top_builddir)/utils/libcoot-utils.la             \
   $(CLIPPER_LIBS) $(MMDB_LIBS) $(GTK_LIBS) -lpthread

libcoot_validation_graphs_la_LDFLAGS = $(SHARED_LDFLAGS)

check_PROGRAMS = gtk4-test-validation-graph gtk4-test-sequence-view

gtk4_test_validation_graph_SOURCES = gtk4-test-validation-graph.cc

gtk4_test_validation_graph_LDADD = ./libcoot-validation-graphs.la \
        $(top_builddir)/ideal/libcoot-ideal.la            \
	$(top_builddir)/ligand/libcoot-ligand.la          \
	$(top_builddir)/coords/libcoot-coords.la          \
        $(top_builddir)/coot-utils/libcoot-map-utils.la   \
        $(top_builddir)/coot-utils/libcoot-coord-utils.la \
        $(top_builddir)/geometry/libcoot-geometry.la      \
	$(CLIPPER_LIBS) $(MMDB_LIBS) $(GTK_LIBS)

gtk4_test_sequence_view_SOURCES = gtk4-test-sequence-view.cc

gtk4_test_sequence_view_LDADD = ./libcoot-validation-graphs.la \
	$(top_builddir)/coords/libcoot-coords.la               \
        $(top_builddir)/coot-utils/libcoot-coord-utils.la      \
        $(top_builddir)/geometry/libcoot-geometry.la           \
	$(MMDB_LIBS) $(GTK_LIBS)
