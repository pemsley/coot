#
AM_CPPFLAGS = \
	-I$(top_srcdir) \
	$(ENHANCED_LIGAND_TOOLS_CPPFLAGS) \
	$(COOT_GSL_FLAGS)		\
	$(GSL_CFLAGS)			\
	$(CLIPPER_CXXFLAGS)             \
	$(MMDB_CXXFLAGS) \
	$(PYTHON_CPPFLAGS)  \
	$(RDKIT_CXXFLAGS) $(BOOST_CPPFLAGS) \
	$(COOT_WINDOWS_CFLAGS)

pkgincludedir = $(includedir)/coot/api

pkginclude_HEADERS = \
coot_molecule.hh	 make-a-dodec.hh	                                           \
colour-functions.hh	 gensurf.hh		    validation-information.hh              \
coot-colour.hh		 ghost-molecule-display.hh  residue-validation-information.hh      \
molecules_container.hh	 ncs.hh			    geometry-graphs.hh                     \
json.hpp 		 rama-plot-phi-psi.hh	    superpose-results.hh                   \
add-terminal-residue.hh  interfaces.hh              moorhen-h-bonds.hh                     \
atom-pull.hh             rigid-body-fit.hh          svg-store-key.hh                       \
phi-psi-prob.hh          saved-strand-info.hh       simple-molecule.hh                     \
generic-3d-lines.hh      cell.hh                    coot-simple-molecule.hh

lib_LTLIBRARIES=libcoot-api.la

libcoot_api_la_SOURCES = \
coot_molecule.cc	     coot-colour.cc	     add-terminal-residue.cc  rama-plot-phi-psi.cc   \
colour-functions.cc          coot_molecule_maps.cc                            coot_molecule_merge_molecules.cc \
coot_molecule_bonds.cc	     coot_molecule_validation.cc                      molecules_container_modelling.cc \
molecules_container.cc       make-a-dodec.cc	                              molecules_container_ligand_fitting.cc \
phi-psi-prob.cc              coot_molecule_moltris.cc                         coot_molecule_replace_fragment.cc \
coot_molecule_backup.cc      coot_molecule_bonds_instanced.cc                 molecules_container_mmrrcc.cc \
coot_molecule_refine.cc      molecules_container_superpose.cc                 rigid-body-fit.cc \
moorhen-h-bonds.cc

libcoot_api_la_LDFLAGS = $(SHARED_LDFLAGS)            \
   $(top_builddir)/MoleculesToTriangles/CXXClasses/libMoleculesToTrianglesCXXClasses.la \
   $(top_builddir)/MoleculesToTriangles/CXXSurface/libMoleculesToTrianglesCXXSurface.la \
   $(top_builddir)/ligand/libcoot-ligand.la           \
   $(top_builddir)/pli/libcoot-pli.la                 \
   $(top_builddir)/ideal/libcoot-ideal.la             \
   $(top_builddir)/lidia-core/libcoot-lidia-core.la   \
   $(top_builddir)/cootaneer/libcootaneer.la          \
   $(top_builddir)/db-main/libcoot-db-main.la         \
   $(top_builddir)/ccp4mg-utils/libccp4mg-utils.la    \
   $(top_builddir)/analysis/libcoot-analysis.la	      \
   $(top_builddir)/density-contour/libcoot-isosurface.la \
   $(top_builddir)/coot-utils/libcoot-map-heavy.la    \
   $(top_builddir)/coot-utils/libcoot-coord-extras.la \
   $(top_builddir)/coot-utils/libcoot-coord-utils.la  \
   $(top_builddir)/coot-utils/libcoot-map-utils.la    \
   $(top_builddir)/geometry/libcoot-geometry.la       \
   $(top_builddir)/mini-mol/libcoot-mini-mol.la       \
   $(top_builddir)/coords/libcoot-coords.la	      \
   $(top_builddir)/utils/libcoot-utils.la             \
   $(CLIPPER_LIBS) $(MMDB_LIBS) $(SQLITE3_LIBS) $(GSL_LIBS) $(RDKIT_LIBS) -lpthread

check_PROGRAMS = test-molecules-container

test_molecules_container_SOURCES = test_molecules_container.cc

test_molecules_container_LDFLAGS = ./libcoot-api.la \
   $(top_builddir)/MoleculesToTriangles/CXXClasses/libMoleculesToTrianglesCXXClasses.la \
   $(top_builddir)/MoleculesToTriangles/CXXSurface/libMoleculesToTrianglesCXXSurface.la \
   $(top_builddir)/ligand/libcoot-ligand.la           \
   $(top_builddir)/pli/libcoot-pli.la                 \
   $(top_builddir)/density-contour/libcoot-isosurface.la \
   $(top_builddir)/coot-utils/libcoot-map-utils.la    \
   $(top_builddir)/coot-utils/libcoot-coord-extras.la \
   $(top_builddir)/coot-utils/libcoot-coord-utils.la \
   $(top_builddir)/geometry/libcoot-geometry.la      \
   $(top_builddir)/coords/libcoot-coords.la	     \
   $(top_builddir)/utils/libcoot-utils.la            \
   $(CLIPPER_LIBS) $(MMDB_LIBS) 

# ---------- pythonic ------------------

dist_python_PYTHON = coot_headless_api.py

pyexec_LTLIBRARIES = _coot_headless_api.la

wrap-coot-api.cc: coot-headless-api.i molecules_container.hh
	swig -o $@ -python -c++ $<

_coot_headless_api_la_SOURCES = wrap-coot-api.cc molecules_container_swig.cc

_coot_headless_api_la_CXXFLAGS = -DUSE_PYTHON -DSWIG

_coot_headless_api_la_LIBADD = libcoot-api.la

_coot_headless_api_la_LDFLAGS = -module $(SHARED_LDFLAGS) --no-undefined
