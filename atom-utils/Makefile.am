
lib_LTLIBRARIES=libcoot-atomutils.la

include_HEADERS = atom-utils.h

GTKGL_CONTEXT_CFLAGS = $(GTKGL_CFLAGS) $(GTKGLEXT_CFLAGS)

INCLUDES =  \
	-I$(top_srcdir)/coords 		\
	-I$(top_srcdir)/geometry 	\
	-I$(top_srcdir)/coot-utils	\
	-I$(top_srcdir)/angles		\
	-I$(top_srcdir)/build		\
	-I$(top_srcdir)/mini-mol	\
	-I$(top_srcdir)/db-main		\
	-I$(top_srcdir)/sequence-view 	\
	-I$(top_srcdir)/ideal	 	\
	-I$(top_srcdir)/coot-surface	\
	-I$(top_srcdir)/surface		\
	-I$(top_srcdir)/src		\
	$(MMDB_CXXFLAGS) 		\
	$(MMDBSSM_CXXFLAGS) 		\
	$(CLIPPER_CXXFLAGS) 		\
	$(GTKCANVAS_CFLAGS) 		\
	$(COOT_GSL_FLAGS) 		\
	$(GSL_CFLAGS) 		 	\
	$(GLUT_CFLAGS)			\
	$(GTKGL_CONTEXT_CFLAGS)		\
	@GTK_CFLAGS@ 

# we need GTK_CFLAGS because we include graphics_info.h which includes
# gtk/gtk.h (sigh)

libcoot_atomutils_la_SOURCES = atom-utils.cc

