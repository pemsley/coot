

INCLUDES = -I$(top_srcdir)  $(GOOCANVAS_CFLAGS) $(GTK_CFLAGS)        

lib_LTLIBRARIES=libcoot-goograph.la 

libcoot_goograph_la_SOURCES = \
    goograph.cc goograph.h 

libcoot_goograph_la_LIBADD = $(top_builddir)/lidia-core/libcoot-lidia-core.la \
			     $(GOOCANVAS_LIBS) $(GTK_LIBS)

check_PROGRAMS = test-goograph

test_goograph_SOURCES = test-goograph.cc

test_goograph_LDADD = libcoot-goograph.la
