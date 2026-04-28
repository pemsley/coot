
lib_LTLIBRARIES = libcoot-graphs.la

AM_CPPFLAGS = -I$(top_srcdir) $(COOT_WINDOWS_CFLAGS) @GTK_CFLAGS@

pkginclude_HEADERS = \
	graph2d.hh

pkgincludedir = $(includedir)/coot/graphs

libcoot_graphs_la_SOURCES = \
	graph2d.cc

libcoot_graphs_la_LIBADD = $(GTK_LIBS)

libcoot_graphs_la_LDFLAGS = $(SHARED_LDFLAGS)
