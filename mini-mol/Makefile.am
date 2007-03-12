
lib_LTLIBRARIES=libcoot-mini-mol.la

include_HEADERS =  mini-mol.hh mini-mol-utils.hh

INCLUDES =  \
	-I$(top_srcdir)/coot-utils \
	$(CLIPPER_CXXFLAGS)    \
	$(MMDB_CXXFLAGS)

noinst_PROGRAMS = testminimol

testminimol_SOURCES = testminimol.cc

libcoot_mini_mol_la_SOURCES =  mini-mol.cc mini-mol-utils.cc

testminimol_LDADD = libcoot-mini-mol.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/coot-utils/libcoot-utils.la \
	$(CLIPPER_LIBS) $(MMDB_LIBS)
