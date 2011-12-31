

INCLUDES = -I$(top_srcdir)/lidia-core \
           -I$(top_srcdir)/utils \
     $(GOOCANVAS_CFLAGS) $(GTK_CFLAGS)        

lib_LTLIBRARIES=libcoot-goograph.la 

libcoot_goograph_la_SOURCES = \
    goograph.cc goograph.h 

libcoot_goograph_la_LIBADD = $(top_builddir)/lidia-core/libcoot-lidia-core.la \
			     $(top_builddir)/utils/libcoot-utils.la \
			     $(GOOCANVAS_LIBS) $(GTK_LIBS)

check_PROGRAMS = test-goograph

test_goograph_SOURCES = test-goograph.cc

test_goograph_LDADD = libcoot-goograph.la \
			     $(top_builddir)/utils/libcoot-utils.la
