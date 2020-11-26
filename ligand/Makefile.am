# ligand/Makefile.am
# 
# Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 The University of York
# Copyright 2009 by the University of Oxford
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

lib_LTLIBRARIES=libcoot-ligand.la

EXTRA_DIST = Makefile.simple read-dunbrack.awk

AM_CPPFLAGS = \
	-I$(top_srcdir) \
	$(COOT_GSL_FLAGS)		\
	$(GSL_CFLAGS)			\
	$(CLIPPER_CXXFLAGS) $(MMDB_CXXFLAGS) 

pkginclude_HEADERS = ligand.hh residue_by_phi_psi.hh wligand.hh monomer-utils.hh \
	chi-angles.hh dunbrack.hh helix-placement.hh fast-ss-search.hh        \
	rigid-body.hh ideal-rna.hh rotamer.hh richardson-rotamer.hh           \
	torsion-general.hh primitive-chi-angles.hh dipole.hh    \
	base-pairing.hh backrub-rotamer.hh trace.hh multi-peptide.hh ligands-db.hh \
   scored-coord.hh ligand-score-card.hh scored-molecule.hh map-point-cluster.hh \
	side-chain-densities.hh

pkgincludedir = $(includedir)/coot/ligand

libcoot_ligand_la_SOURCES = ligand.cc ligand-extras.cc  residue_by_phi_psi.cc \
	wiggly-ligand.cc monomer-utils.cc dunbrack.cc chi-angles-autogen.cc   \
	chi-angles.cc helix-placement.cc fast-ss-search.cc rigid-body.cc      \
	ideal-rna.cc richardson-rotamers.cc richardson-rotamer.cc rotamer.cc  \
	torsion-general.cc rotamer-tables.cc primitive-chi-angles.cc dipole.cc\
        base-pairing.cc backrub-rotamer.cc trace.cc trace-2.cc multi-peptide.cc \
        ligands-db.cc joining-residues.cc ligand-score-card.cc side-chain-densities.cc

# chi-angles, monomer-utils, rotamer torsion-general use coot::Cartesian
# from libcoot-coords
#   $(top_builddir)/coords/libcoot-coords.la           # library re-order

# clean
libcoot_ligand_la_LIBADD = \
   $(top_builddir)/ideal/libcoot-ideal.la             \
   $(top_builddir)/coot-utils/libcoot-map-heavy.la    \
   $(top_builddir)/cootaneer/libcootaneer.la \
   $(top_builddir)/db-main/libcoot-db-main.la         \
   $(top_builddir)/ccp4mg-utils/libccp4mg-utils.la    \
   $(top_builddir)/geometry/libcoot-geometry.la       \
   $(top_builddir)/mini-mol/libcoot-mini-mol.la       \
   $(top_builddir)/coot-utils/libcoot-coord-extras.la \
   $(top_builddir)/coot-utils/libcoot-coord-utils.la  \
   $(top_builddir)/coot-utils/libcoot-map-utils.la    \
   $(top_builddir)/analysis/libcoot-analysis.la	      \
   $(top_builddir)/utils/libcoot-utils.la             \
   $(CLIPPER_LIBS) $(MMDB_LIBS) $(SQLITE3_LIBS) $(RDKIT_LIBS)

libcoot_ligand_la_LDFLAGS = $(SHARED_LDFLAGS)

AM_CXXFLAGS = \
         $(OPENMP_CFLAGS) \
	-DPKGDATADIR='"$(pkgdatadir)"' \
	$(SQLITE3_CFLAGS) \
	$(COOT_WINDOWS_CFLAGS)

# ligands-2016.db is made and updated by make-ligands-db
#
cootdatadir = $(pkgdatadir)/data

dist_cootdata_DATA = \
	ligands-2016.db


if OS_WIN32
bin_PROGRAMS = findligand-bin findwaters-bin coot-make-ligands-db
else
libexec_PROGRAMS = findligand-bin findwaters-bin coot-make-ligands-db
endif

# check_PROGRAMS = test-helix-placement test-fffear-in-coot testdunbrack 
#                  test-random test-ligand test-rotamer-tables-sql

check_PROGRAMS = test-fffear-in-coot test-trace ligand-percentiles test-side-chain-densities \
		test-ligand

findligand_bin_SOURCES = find-ligand.cc
findwaters_bin_SOURCES = find-waters.cc

test_trace_SOURCES = test-trace.cc
test_trace_LDADD = ./libcoot-ligand.la $(top_builddir)/compat/libcoot-compat.la

coot_make_ligands_db_SOURCES = make-ligands-db.cc
coot_make_ligands_db_LDADD = ./libcoot-ligand.la $(top_builddir)/compat/libcoot-compat.la


ligand_percentiles_SOURCES = ligand-percentiles.cc
ligand_percentiles_LDADD = ./libcoot-ligand.la $(top_builddir)/compat/libcoot-compat.la

test_side_chain_densities_SOURCES = test-side-chain-densities.cc
test_side_chain_densities_LDADD = \
  ./libcoot-ligand.la \
   $(top_builddir)/coot-utils/libcoot-map-utils.la \
   $(top_builddir)/coot-utils/libcoot-coord-utils.la \
   $(top_builddir)/coords/libcoot-coords.la        \
   $(top_builddir)/geometry/libcoot-geometry.la \
   $(top_builddir)/utils/libcoot-utils.la \
   $(CLIPPER_LIBS) $(MMDB_LIBS)

# test_rotamer_tables_sql = test-rotamer-tables-sql.cc
#
# test_rotamer_tables_sql_LDADD = ./libcoot-ligand.la  $(MMDB_LIBS) $(install_dir)/lib/libsqlite3.la

test_ligand_SOURCES = test-ligand.cc

# test_helix_placement_SOURCES = test-helix-placement.cc

# test_random_SOURCES = test-random.cc

#test_random_LDADD = \
#	./libcoot-ligand.la                                \
#	$(top_builddir)/ideal/libcoot-ideal.la             \
#	$(top_builddir)/coords/libcoot-coords.la           \
#	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
#	$(top_builddir)/coot-utils/libcoot-map-utils.la    \
#	$(top_builddir)/db-main/libcoot-db-main.la         \
#	$(top_builddir)/geometry/libcoot-geometry.la       \
#	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
#	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
#	$(top_builddir)/utils/libcoot-utils.la             \
#	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la    \
#	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
#	$(top_builddir)/compat/libcoot-compat.la           \
#	$(CLIPPER_LIBS) $(MMDB_LIBS) $(GSL_LIBS)

test_ligand_LDADD = \
	./libcoot-ligand.la

#test_helix_placement_LDADD = \
#	./libcoot-ligand.la                                \
#	$(top_builddir)/ideal/libcoot-ideal.la             \
#	$(top_builddir)/coords/libcoot-coords.la           \
#	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
#	$(top_builddir)/coot-utils/libcoot-map-utils.la    \
#	$(top_builddir)/geometry/libcoot-geometry.la       \
#	$(top_builddir)/db-main/libcoot-db-main.la         \
#	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
#	$(top_builddir)/utils/libcoot-utils.la             \
#	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
#	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la    \
#	$(top_builddir)/compat/libcoot-compat.la           \
#	$(CLIPPER_LIBS) $(MMDB_LIBS) $(GSL_LIBS)

test_fffear_in_coot_SOURCES = test-fffear-in-coot.cc

# 	$(top_builddir)/coords/libcoot-coords.la

test_fffear_in_coot_LDADD = \
	./libcoot-ligand.la \
	$(top_builddir)/ideal/libcoot-ideal.la             \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
	$(top_builddir)/geometry/libcoot-geometry.la       \
	$(top_builddir)/coot-utils/libcoot-map-heavy.la    \
	$(top_builddir)/db-main/libcoot-db-main.la         \
	$(top_builddir)/coot-utils/libcoot-map-utils.la    \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
	$(top_builddir)/utils/libcoot-utils.la             \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la    \
	$(top_builddir)/compat/libcoot-compat.la           \
	$(CLIPPER_LIBS) $(MMDB_LIBS) $(GSL_LIBS)

# testwlig_SOURCES = testwlig.cc
#	-L$(top_builddir)/high-res     -lcoot-high-res    
#	-I$(top_srcdir)/high-res

# -L$(top_builddir)/ideal        -lcoot-ideal       

#	$(top_builddir)/coords/libcoot-coords.la

findligand_bin_LDADD = \
	./libcoot-ligand.la                                \
	$(top_builddir)/ideal/libcoot-ideal.la             \
        $(top_builddir)/coot-utils/libcoot-map-heavy.la    \
	$(top_builddir)/geometry/libcoot-geometry.la       \
	$(top_builddir)/db-main/libcoot-db-main.la         \
	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la    \
	$(top_builddir)/coot-utils/libcoot-map-utils.la    \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
	$(top_builddir)/utils/libcoot-utils.la             \
        $(top_builddir)/lidia-core/libcoot-lidia-core.la   \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
	$(top_builddir)/compat/libcoot-compat.la           \
        $(GLOB_LIBS) \
        $(RDKIT_LIBS) \
	$(CLIPPER_LIBS) $(MMDB_LIBS) $(GSL_LIBS)

# 	$(top_builddir)/coords/libcoot-coords.la

findwaters_bin_LDADD = \
	./libcoot-ligand.la                                \
	$(top_builddir)/ideal/libcoot-ideal.la             \
        $(top_builddir)/coot-utils/libcoot-map-heavy.la    \
	$(top_builddir)/geometry/libcoot-geometry.la       \
	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la    \
	$(top_builddir)/db-main/libcoot-db-main.la         \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
	$(top_builddir)/coot-utils/libcoot-map-utils.la    \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
	$(top_builddir)/utils/libcoot-utils.la             \
        $(top_builddir)/lidia-core/libcoot-lidia-core.la   \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
	$(top_builddir)/compat/libcoot-compat.la           \
        $(GLOB_LIBS) \
        $(RDKIT_LIBS) \
	$(CLIPPER_LIBS) $(MMDB_LIBS) $(GSL_LIBS)

#testwlig_LDADD = -L. -lcoot-ligand      \
#	-L../mini-mol -lcoot-mini-mol    \
#	-L../geometry -lcoot-geometry    \
#	-L../coords -lcoot-coords        \
#	-L../ccp4mg-utils -lccp4mg-utils \
#	$(CLIPPER_LIBS) $(MMDB_LIBS)


#  -L../coords -lcoot-coords -L../atom-utils -lcoot-atomutils

# testori_SOURCES = testori.cc

dist_bin_SCRIPTS = findwaters findligand

#testdunbrack_SOURCES = testdunbrack.cc

#testdunbrack_LDADD = ./libcoot-ligand.la                   \
#	$(top_builddir)/ideal/libcoot-ideal.la             \
#	$(top_builddir)/geometry/libcoot-geometry.la       \
#	$(top_builddir)/coords/libcoot-coords.la           \
#	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la    \
#	$(top_builddir)/coot-utils/libcoot-map-utils.la    \
#	$(top_builddir)/db-main/libcoot-db-main.la         \
#	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
#	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
#	$(top_builddir)/utils/libcoot-utils.la             \
#	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
#        $(CLIPPER_LIBS) $(MMDB_LIBS) $(GSL_LIBS)

