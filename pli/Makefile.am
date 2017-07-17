

AM_CPPFLAGS = \
	-I$(top_srcdir) \
	$(COOT_GSL_FLAGS)		\
	$(GSL_CFLAGS)			\
	$(CLIPPER_CXXFLAGS) $(MMDB_CXXFLAGS) 

lib_LTLIBRARIES=libcoot-pli.la

libcoot_pli_la_SOURCES = protein-ligand-interactions.cc flev-annotations.cc

pkginclude_HEADERS = flev-annotations.hh protein-ligand-interactions.hh 

pkgincludedir = $(includedir)/coot/pli


