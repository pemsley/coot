#
# Makefile.am
# 
# Copyright 2020 Medical Research Council
# Author: Paul Emsley
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA

## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = \
	-I$(top_srcdir)                     \
	$(ENHANCED_LIGAND_TOOLS_CPPFLAGS)   \
	$(gsl_CFLAGS)			    \
	$(LIBSSM_CXXFLAGS)                  \
	$(CLIPPER_CXXFLAGS)                 \
	$(MMDB_CXXFLAGS)                    \
	$(GEMMI_CXXFLAGS)                   \
	$(PYTHON_CPPFLAGS)                  \
	$(CAIRO_CFLAGS)                     \
	$(RDKIT_CXXFLAGS) $(BOOST_CPPFLAGS) \
	$(COOT_WINDOWS_CFLAGS)

pkgincludedir = $(includedir)/coot/api

pkginclude_HEADERS = \
coot-molecule.hh	 make-a-dodec.hh	    lsq-results.hh                         \
colour-functions.hh	 gensurf.hh		    validation-information.hh              \
coot-colour.hh		 ghost-molecule-display.hh  residue-validation-information.hh      \
molecules-container.hh	 ncs.hh			    geometry-graphs.hh                     \
header-info.hh 		 rama-plot-phi-psi.hh	    superpose-results.hh                   \
add-terminal-residue.hh  interfaces.hh              moorhen-h-bonds.hh                     \
atom-pull.hh             rigid-body-fit.hh          svg-store-key.hh                       \
phi-psi-prob.hh          saved-strand-info.hh       simple-molecule.hh                     \
generic-3d-lines.hh      cell.hh                    coot-simple-molecule.hh                \
api-cell.hh              bond-colour.hh             make-instanced-graphical-bonds.hh      \
symmetry-info.hh         instancing.hh              merge-molecule-results-info-t.hh       \
coot-molecule-bonds.hh   clipper-ccp4-map-file-wrapper.hh blender-mesh.hh                  \
ncs-ghosts.hh            moved-atom.hh              moved-residue.hh \
filo-tests.hh            plain-atom-overlap.hh      test-molecules-container.hh \
lucrezia-tests.hh        user-defined-colour-table.hh positioned-atom-spec.hh

lib_LTLIBRARIES=libcoot-api.la

libcoot_api_la_SOURCES = \
coot-molecule.cc	     coot-colour.cc	     add-terminal-residue.cc  rama-plot-phi-psi.cc   \
colour-functions.cc          coot-molecule-maps.cc                            coot-molecule-merge-molecules.cc \
coot-molecule-bonds.cc	     coot-molecule-validation.cc                      molecules-container-modelling.cc \
molecules-container.cc       make-a-dodec.cc	                              molecules-container-ligand-fitting.cc \
phi-psi-prob.cc              coot-molecule-moltris.cc                         coot-molecule-replace-fragment.cc \
coot-molecule-backup.cc      coot-molecule-bonds-instanced.cc                 molecules-container-mmrrcc.cc \
molecules-container-gltf.cc  mmcif-tests.cc                                   ncs-ghosts.cc \
coot-molecule-refine.cc      molecules-container-superpose.cc                 rigid-body-fit.cc \
moorhen-h-bonds.cc           molecules-container-bucca-ml-growing.cc          instancing.cc  \
molecules-container-maps.cc  coot-molecule-change-chain-id.cc                 coot-molecule-modelling.cc \
blender-mesh.cc              coot-molecule-blender.cc                         molecules-container-blender.cc \
ghost-molecule-display.cc    model-analysis.cc                                coot-molecule-json.cc \
molecules-container-pyrogen.cc

libcoot_api_la_LDFLAGS = $(SHARED_LDFLAGS)            \
   $(top_builddir)/MoleculesToTriangles/CXXClasses/libMoleculesToTrianglesCXXClasses.la \
   $(top_builddir)/MoleculesToTriangles/CXXSurface/libMoleculesToTrianglesCXXSurface.la \
   $(top_builddir)/ligand/libcoot-ligand.la           \
   $(top_builddir)/pli/libcoot-pli.la                 \
   $(top_builddir)/ideal/libcoot-ideal.la             \
   $(top_builddir)/lidia-core/libcoot-lidia-core.la   \
   $(top_builddir)/pyrogen/libcoot_py_restraints.la   \
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
   $(CLIPPER_LIBS) $(LIBSSM_LIBS) $(MMDB_LIBS) $(SQLITE3_LIBS) $(gsl_LIBS) $(RDKIT_LIBS) -lpthread

check_PROGRAMS = test-molecules-container

test_molecules_container_SOURCES = test-molecules-container.cc filo-tests.cc lucrezia-tests.cc

test_molecules_container_LDFLAGS = ./libcoot-api.la \
   $(top_builddir)/MoleculesToTriangles/CXXClasses/libMoleculesToTrianglesCXXClasses.la \
   $(top_builddir)/MoleculesToTriangles/CXXSurface/libMoleculesToTrianglesCXXSurface.la \
   $(top_builddir)/pyrogen/libcoot_py_restraints.la    \
   $(top_builddir)/ligand/libcoot-ligand.la           \
   $(top_builddir)/ideal/libcoot-ideal.la             \
   $(top_builddir)/pli/libcoot-pli.la                 \
   $(top_builddir)/density-contour/libcoot-isosurface.la \
   $(top_builddir)/coot-utils/libcoot-map-utils.la    \
   $(top_builddir)/coot-utils/libcoot-coord-extras.la \
   $(top_builddir)/coot-utils/libcoot-coord-utils.la \
   $(top_builddir)/geometry/libcoot-geometry.la      \
   $(top_builddir)/coords/libcoot-coords.la	     \
   $(top_builddir)/utils/libcoot-utils.la            \
   $(CLIPPER_LIBS) $(GEMMI_LIBS) $(MMDB_LIBS) $(RDKIT_LIBS)

# 2025-07-03-PE swig stuff was here.
