
lib_LTLIBRARIES=libcoot-high-res.la

include_HEADERS = high-res.hh coot-atom-graph.hh coot-node-info.hh \
	sequence-assignment.hh


AM_CXXFLAGS = $(COOT_GSL_FLAGS) $(GSL_CFLAGS)

INCLUDES = \
	-I$(top_srcdir)/geometry     \
	-I$(top_srcdir)/coords       \
	-I$(top_srcdir)/ccp4mg-utils \
	-I$(top_srcdir)/mini-mol     \
	-I$(top_srcdir)/ligand       \
	-I$(top_srcdir)/coot-utils   \
	-I$(top_srcdir)/compat       \
	$(MMDB_CXXFLAGS)  \
	$(CLIPPER_CXXFLAGS)

libcoot_high_res_la_SOURCES = high-res.cc coot-atom-graph.cc residue-distortions.cc \
	sequence-assignment.cc

noinst_PROGRAMS = testatomgraph tracehighres testseqalign testseqassignment globularize_protein

testatomgraph_SOURCES = test-atom-graph.cc

testatomgraph_LDADD =  libcoot-high-res.la                \
	$(top_builddir)/ligand/libcoot-ligand.la          \
	$(top_builddir)/ideal/libcoot-ideal.la            \
	$(top_builddir)/geometry/libcoot-geometry.la      \
	$(top_builddir)/coords/libcoot-coords.la          \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la      \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/coot-utils/libcoot-map-utils.la   \
	$(top_builddir)/coot-utils/libcoot-shelx.la       \
	$(top_builddir)/coot-utils/libcoot-utils.la       \
	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la   \
	$(top_builddir)/compat/libcoot-compat.la          \
	$(GSL_LIBS)			                  \
	$(CLIPPER_LIBS)

testseqassignment_SOURCES = testseqassignment.cc
testseqassignment_LDADD = libcoot-high-res.la             \
	$(top_builddir)/ligand/libcoot-ligand.la          \
	$(top_builddir)/ideal/libcoot-ideal.la            \
	$(top_builddir)/geometry/libcoot-geometry.la      \
	$(top_builddir)/coords/libcoot-coords.la          \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la      \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
	$(top_builddir)/coot-utils/libcoot-map-utils.la   \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/coot-utils/libcoot-shelx.la       \
	$(top_builddir)/coot-utils/libcoot-utils.la       \
	$(top_builddir)/compat/libcoot-compat.la          \
	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la   \
	$(GSL_LIBS)				          \
	$(CLIPPER_LIBS)

# tracehighres_LDFLAGS = -static 
tracehighres_SOURCES = trace-high-res.cc
tracehighres_LDADD = libcoot-high-res.la                  \
	$(top_builddir)/ligand/libcoot-ligand.la          \
	$(top_builddir)/ideal/libcoot-ideal.la            \
	$(top_builddir)/geometry/libcoot-geometry.la      \
	$(top_builddir)/coords/libcoot-coords.la          \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la      \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
	$(top_builddir)/coot-utils/libcoot-map-utils.la   \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/coot-utils/libcoot-shelx.la       \
	$(top_builddir)/coot-utils/libcoot-utils.la       \
	$(top_builddir)/compat/libcoot-compat.la          \
	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la   \
	$(GSL_LIBS)				          \
	$(CLIPPER_LIBS)


testseqalign_SOURCES = testseqalign.cc

# why do we need to link with all this (junk)?
testseqalign_LDADD = libcoot-high-res.la                  \
	$(top_builddir)/ligand/libcoot-ligand.la          \
	$(top_builddir)/ideal/libcoot-ideal.la            \
	$(top_builddir)/geometry/libcoot-geometry.la      \
	$(top_builddir)/coords/libcoot-coords.la          \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la      \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
	$(top_builddir)/coot-utils/libcoot-map-utils.la   \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/coot-utils/libcoot-shelx.la       \
	$(top_builddir)/coot-utils/libcoot-utils.la       \
	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la   \
	$(top_builddir)/compat/libcoot-compat.la          \
	$(GSL_LIBS)				          \
	$(CLIPPER_LIBS)


globularize_protein_SOURCES = globularize.cc

globularize_protein_LDADD = libcoot-high-res.la           \
	$(top_builddir)/ligand/libcoot-ligand.la          \
	$(top_builddir)/ideal/libcoot-ideal.la            \
	$(top_builddir)/geometry/libcoot-geometry.la      \
	$(top_builddir)/coords/libcoot-coords.la          \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la      \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
	$(top_builddir)/coot-utils/libcoot-map-utils.la   \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/coot-utils/libcoot-shelx.la       \
	$(top_builddir)/coot-utils/libcoot-utils.la       \
	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la   \
	$(top_builddir)/compat/libcoot-compat.la          \
	$(GSL_LIBS)				          \
	$(CLIPPER_LIBS)


