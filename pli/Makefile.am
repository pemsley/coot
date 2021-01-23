

AM_CPPFLAGS = \
	-I$(top_srcdir) \
	$(ENHANCED_LIGAND_TOOLS_CFLAGS) \
	$(COOT_GSL_FLAGS)		\
	$(GSL_CFLAGS)			\
	$(CLIPPER_CXXFLAGS)             \
        $(MMDB_CXXFLAGS) \
        $(PYTHON_CFLAGS)  \
	$(RDKIT_CXXFLAGS) $(BOOST_CPPFLAGS) \
	$(COOT_WINDOWS_CFLAGS)

lib_LTLIBRARIES=libcoot-pli.la

libcoot_pli_la_SOURCES = protein-ligand-interactions.cc flev-annotations.cc pi-stacking.cc specs.cc

libcoot_pli_la_LIBADD = $(top_builddir)/coot-utils/libcoot-coord-utils.la \
                        $(top_builddir)/geometry/libcoot-geometry.la      \
                        $(top_builddir)/lidia-core/libcoot-lidia-core.la  \
                        $(CLIPPER_LIBS) $(MMDB_LIBS)                      \
                        $(RDKIT_LIBS) $(BOOST_LDFLAGS) $(L_BOOST_PYTHON)  \
                        $(PYTHON_LIBS)


libcoot_pli_la_LDFLAGS = $(SHARED_LDFLAGS)

pkginclude_HEADERS = flev-annotations.hh protein-ligand-interactions.hh pi-stacking.hh specs.hh

pkgincludedir = $(includedir)/coot/pli


