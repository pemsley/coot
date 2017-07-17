

AM_CPPFLAGS = \
	-I$(top_srcdir) \
	$(ENHANCED_LIGAND_TOOLS_CFLAGS) \
	$(COOT_GSL_FLAGS)		\
	$(GSL_CFLAGS)			\
	$(CLIPPER_CXXFLAGS) $(MMDB_CXXFLAGS) \
	$(RDKIT_CXXFLAGS)

lib_LTLIBRARIES=libcoot-pli.la

libcoot_pli_la_SOURCES = protein-ligand-interactions.cc flev-annotations.cc pi-stacking.cc

pkginclude_HEADERS = flev-annotations.hh protein-ligand-interactions.hh pi-stacking.hh

pkgincludedir = $(includedir)/coot/pli


