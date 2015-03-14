
AM_CPPFLAGS = -I$(top_srcdir) $(GLIB_CFLAGS) $(GTK_CFLAGS)

lib_LTLIBRARIES=libcoot-gtkgraph.la

libcoot_gtkgraph_la_SOURCES = \
    annotation.c  Doxyfile	  gtkgraph.h	       \
    gtkgraph_original.h  polar_util.c  trace.c         \
    axis.c	      gtkgraph.c  gtkgraph_internal.h  \
    polar.c		    smith.c

libcoot_gtkgraph_la_LIBADD = $(GTK_LIBS) -lm

libcoot_gtkgraph_la_LDFLAGS = $(SHARED_LDFLAGS)

