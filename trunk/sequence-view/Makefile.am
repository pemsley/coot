
lib_LTLIBRARIES=libcoot-sequence-view.la

include_HEADERS = sequence-view.hh seq-view-interface.h

CANVAS_CFLAGS = $(GTKCANVAS_CFLAGS) $(GNOME_CANVAS_CFLAGS)
CANVAS_LIBS   = $(GTKCANVAS_LIBS) $(GNOME_CANVAS_LIBS)

# we need these include because we use graphics-info for the
# recentering callback - we don't use c-interface because we dont 
# want the function that returns a sequence_view object to be included
# into c-interface.h (mmdb poisoning).
INCLUDES = \
	-I $(top_srcdir)/coords		\
	-I $(top_srcdir)/coot-utils	\
	-I $(top_srcdir)/src		\
	-I $(top_srcdir)/geometry	\
	-I $(top_srcdir)/mini-mol	\
	-I $(top_srcdir)/build		\
	-I $(top_srcdir)/db-main	\
	-I $(top_srcdir)/angles 	\
	-I $(top_srcdir)/ideal 		\
	-I $(top_srcdir)/coot-surface	\
	-I $(top_srcdir)/surface	\
	$(CLIPPER_CXXFLAGS) 		\
	$(MMDB_CXXFLAGS) 		\
	$(CANVAS_CFLAGS) 		\
	@GTK_CFLAGS@ 

# CPPFLAGS = -DSEQ_VIEW_STANDALONE
AM_CXXFLAGS = $(COOT_GSL_FLAGS) $(GSL_CFLAGS) $(COOT_WINDOWS_CFLAGS)

libcoot_sequence_view_la_SOURCES = sequence-view.cc  seq-view-interface.c

# bin_PROGRAMS = testsequenceview

# testsequenceview_SOURCES = testsequenceview.cc

# testsequenceview_LDADD = \
#	-L../coords -lcoot-coords	\
#	-L. -lcoot-sequence-view 	\
#	$(GTKCANVAS_LIBS) 	 	\
#	$(MMDB_LIBS)		 	\
#	@GTK_LIBS@	
