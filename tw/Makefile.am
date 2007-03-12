
# Do it properly with configure (not like this)
# OGG_VORBIS_CFLAGS=-I/sw/include/vorbis -I/usr/include/vorbis
# AO_CFLAGS=-I/sw/include -I/usr/include

AM_CXXFLAGS = \
	-DPKGDATADIR='"$(pkgdatadir)"'

INCLUDES = \
	$(GTKCANVAS_CFLAGS)  \
	$(GTK_CFLAGS)	     \
	$(OGG_VORBIS_CFLAGS) \
	$(AO_CFLAGS)

# don't make this yet, just carry it around.
EXTRA_DIST = cootsound.c

lib_LTLIBRARIES=libcoot-tw.la

include_HEADERS = coot-tw.hh

libcoot_tw_la_SOURCES = tw.cc

noinst_PROGRAMS = coot-tiddly

coot_tiddly_SOURCES = tiddly.cc

coot_tiddly_LDADD = libcoot-tw.la \
	$(GTKCANVAS_LIBS) \
	@GTK_LIBS@


# cootsound_LDADD = -L/sw/lib -lvorbisfile -logg -lao

