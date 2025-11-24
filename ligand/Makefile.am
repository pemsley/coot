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

lib_LTLIBRARIES=libcoot-ligand.la libcoot-res-tracer.la

EXTRA_DIST = Makefile.simple read-dunbrack.awk side-chain-data.tar.gz

AM_CPPFLAGS = \
	-I$(top_srcdir)         \
	$(gsl_CFLAGS)		\
	$(CLIPPER_CXXFLAGS)     \
	$(RDKIT_CXXFLAGS)       \
        $(BOOST_CPPFLAGS)       \
        $(MMDB_CXXFLAGS)

pkginclude_HEADERS = ligand.hh residue_by_phi_psi.hh wligand.hh monomer-utils.hh \
	chi-angles.hh dunbrack.hh helix-placement.hh fast-ss-search.hh        \
	rigid-body.hh ideal-rna.hh rotamer.hh richardson-rotamer.hh           \
	torsion-general.hh primitive-chi-angles.hh dipole.hh    \
	base-pairing.hh backrub-rotamer.hh trace.hh multi-peptide.hh ligands-db.hh \
	scored-coord.hh ligand-score-card.hh scored-molecule.hh map-point-cluster.hh \
	side-chain-densities.hh tube-finder.hh smith-waterman.hh side-chain.hh \
	new-residue-by-3-phi-psi.hh refine-isolated-chain.hh libres-tracer.hh \
	dSFMT.h dSFMT-params.h dSFMT-params19937.h dSFMT-common.h scored-node.hh \
        rama-rsr-extend-fragments.hh

pkgincludedir = $(includedir)/coot/ligand

libcoot_ligand_la_SOURCES = ligand.cc ligand-extras.cc  residue_by_phi_psi.cc \
	wiggly-ligand.cc monomer-utils.cc dunbrack.cc chi-angles-autogen.cc   \
	chi-angles.cc helix-placement.cc fast-ss-search.cc rigid-body.cc      \
	ideal-rna.cc richardson-rotamers.cc richardson-rotamer.cc rotamer.cc  \
	torsion-general.cc rotamer-tables.cc primitive-chi-angles.cc dipole.cc\
        base-pairing.cc backrub-rotamer.cc trace.cc trace-2.cc multi-peptide.cc \
        ligands-db.cc joining-residues.cc ligand-score-card.cc side-chain-densities.cc \
	tube-finder.cc smith-waterman.cc side-chain.cc refine-isolated-chain.cc

# chi-angles, monomer-utils, rotamer torsion-general use coot::Cartesian
# from libcoot-coords
#   $(top_builddir)/coords/libcoot-coords.la           # library re-order

# clean
libcoot_ligand_la_LIBADD = \
   $(top_builddir)/ideal/libcoot-ideal.la             \
   $(top_builddir)/coot-utils/libcoot-map-heavy.la    \
   $(top_builddir)/lidia-core/libcoot-lidia-core.la   \
   $(top_builddir)/cootaneer/libcootaneer.la          \
   $(top_builddir)/db-main/libcoot-db-main.la         \
   $(top_builddir)/ccp4mg-utils/libccp4mg-utils.la    \
   $(top_builddir)/geometry/libcoot-geometry.la       \
   $(top_builddir)/mini-mol/libcoot-mini-mol.la       \
   $(top_builddir)/coot-utils/libcoot-coord-extras.la \
   $(top_builddir)/coot-utils/libcoot-coord-utils.la  \
   $(top_builddir)/coot-utils/libcoot-map-utils.la    \
   $(top_builddir)/analysis/libcoot-analysis.la	      \
   $(top_builddir)/utils/libcoot-utils.la             \
   $(CLIPPER_LIBS) $(MMDB_LIBS) $(SQLITE3_LIBS) $(gsl_LIBS) -lpthread

libcoot_ligand_la_LDFLAGS = $(SHARED_LDFLAGS)

# 20210818-PE consolidating libcoot-res-tracer into libcoot-ligand created various head-scratching memory errors.
libcoot_res_tracer_la_SOURCES = libres-tracer.cc rama-rsr-extend-fragments.cc new-residue-by-3-phi-psi.cc dSFMT.c

libcoot_res_tracer_la_LDFLAGS = $(SHARED_LDFLAGS)

libcoot_res_tracer_la_LIBADD = libcoot-ligand.la           \
	$(top_builddir)/ideal/libcoot-ideal.la             \
	$(top_builddir)/coot-utils/libcoot-map-utils.la    \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
	$(top_builddir)/geometry/libcoot-geometry.la       \
	$(top_builddir)/utils/libcoot-utils.la             \
	$(CLIPPER_LIBS) $(MMDB_LIBS)

AM_CXXFLAGS = \
         $(OPENMP_CFLAGS) \
	-DPKGDATADIR='"$(pkgdatadir)"' \
	$(SQLITE3_CFLAGS) \
	$(COOT_WINDOWS_CFLAGS) \
        $(LIBPNG_CFLAGS)

# ligands-2016.db is made and updated by make-ligands-db
#
cootdatadir = $(pkgdatadir)/data

dist_cootdata_DATA = ligands-2016.db

cootdata_assign_side_chain_datadir = $(cootdatadir)/assign-side-chains

dist_cootdata_assign_side_chain_data_DATA = useable-grid-points-nstep=5,box_radius=5.0-charybdis.data

# now think about side-chain-data dir.  Hmm!


# this block is a bit dodgy, I feel. It uses tar to install a directory of data files
.PHONY: install_side_chain_data

# thanks John Bollinger
install_side_chain_data: install-dist_cootdataDATA
	$(MKDIR_P) $(DESTDIR)$(cootdata_assign_side_chain_datadir) && $(AMTAR) -x -C $(DESTDIR)$(cootdata_assign_side_chain_datadir) -f $(top_srcdir)/ligand/side-chain-data.tar.gz

install-data-local: install_side_chain_data


# check_PROGRAMS = test-helix-placement test-fffear-in-coot testdunbrack 
#                  test-random test-ligand test-rotamer-tables-sql

check_PROGRAMS = test-fffear-in-coot test-trace ligand-percentiles  \
		test-ligand test-tube-finder res-tracer

# removed for now wrong build order: test-side-chain-densities

coot_findligand_bin_SOURCES = find-ligand.cc
coot_findwaters_bin_SOURCES = find-waters.cc

test_tube_finder_SOURCES = test-tube-finder.cc
test_tube_finder_LDADD = ./libcoot-ligand.la          \
   $(top_builddir)/coot-utils/libcoot-coord-utils.la  \
   $(CLIPPER_LIBS)

res_tracer_SOURCES = res-tracer.cc

res_tracer_LDADD = ./libcoot-res-tracer.la            \
   libcoot-ligand.la \
   $(top_builddir)/ideal/libcoot-ideal.la             \
   $(top_builddir)/cootaneer/libcootaneer.la          \
   $(top_builddir)/coot-utils/libcoot-coord-utils.la  \
   $(top_builddir)/coot-utils/libcoot-map-utils.la    \
   $(top_builddir)/geometry/libcoot-geometry.la       \
   $(top_builddir)/mini-mol/libcoot-mini-mol.la       \
   $(top_builddir)/utils/libcoot-utils.la             \
   $(CLIPPER_LIBS) $(MMDB_LIBS)

if OS_WIN32
bin_PROGRAMS = coot-findligand-bin coot-findwaters-bin coot-make-ligands-db
else
libexec_PROGRAMS = coot-findligand-bin coot-findwaters-bin coot-make-ligands-db coot-identify-protein-bin
endif

test_trace_SOURCES = test-trace.cc
test_trace_LDADD = ./libcoot-ligand.la                \
   $(top_builddir)/ideal/libcoot-ideal.la             \
   $(top_builddir)/coot-utils/libcoot-coord-utils.la  \
   $(top_builddir)/coot-utils/libcoot-map-utils.la    \
   $(top_builddir)/geometry/libcoot-geometry.la       \
   $(top_builddir)/mini-mol/libcoot-mini-mol.la       \
   $(top_builddir)/utils/libcoot-utils.la             \
   $(top_builddir)/compat/libcoot-compat.la           \
   $(CLIPPER_LIBS) $(MMDB_LIBS)

coot_make_ligands_db_SOURCES = make-ligands-db.cc
coot_make_ligands_db_LDADD = ./libcoot-ligand.la \
	$(top_builddir)/ideal/libcoot-ideal.la             \
        $(top_builddir)/coot-utils/libcoot-map-heavy.la    \
	$(top_builddir)/geometry/libcoot-geometry.la       \
	$(top_builddir)/db-main/libcoot-db-main.la         \
	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la    \
	$(top_builddir)/coot-utils/libcoot-map-utils.la    \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
        $(top_builddir)/cootaneer/libcootaneer.la          \
	$(top_builddir)/utils/libcoot-utils.la             \
        $(top_builddir)/lidia-core/libcoot-lidia-core.la   \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
        $(top_builddir)/analysis/libcoot-analysis.la	   \
	$(top_builddir)/compat/libcoot-compat.la

coot_identify_protein_bin_SOURCES = identify-protein.cc


# Why do I need to add lidia-core here? It's already added to libcoot-coord-utils.
coot_identify_protein_bin_LDADD = ./libcoot-ligand.la      \
	$(top_builddir)/ideal/libcoot-ideal.la             \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
        $(top_builddir)/lidia-core/libcoot-lidia-core.la   \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
	$(top_builddir)/geometry/libcoot-geometry.la       \
	$(top_builddir)/utils/libcoot-utils.la             \
        $(CLIPPER_LIBS) $(MMDB_LIBS)


ligand_percentiles_SOURCES = ligand-percentiles.cc
ligand_percentiles_LDADD = ./libcoot-ligand.la       \
	$(top_builddir)/ideal/libcoot-ideal.la       \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
	$(top_builddir)/geometry/libcoot-geometry.la       \
	$(top_builddir)/utils/libcoot-utils.la             \
	$(top_builddir)/compat/libcoot-compat.la

# wrong build order - coords depends on this directory (because it uses rotamer info)
# solve by changing the directory for test-side-chain-densities to be in coords
# test_side_chain_densities_SOURCES = test-side-chain-densities.cc
# test_side_chain_densities_LDADD = \
#   ./libcoot-ligand.la \
#    $(top_builddir)/coot-utils/libcoot-map-utils.la \
#    $(top_builddir)/coot-utils/libcoot-coord-utils.la \
#    $(top_builddir)/coords/libcoot-coords.la        \
#    $(top_builddir)/geometry/libcoot-geometry.la \
#    $(top_builddir)/utils/libcoot-utils.la \
#    $(CLIPPER_LIBS) $(MMDB_LIBS)

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
	./libcoot-ligand.la \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
	$(top_builddir)/geometry/libcoot-geometry.la       \
	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la    \
	$(CLIPPER_LIBS)

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
	$(CLIPPER_LIBS) $(MMDB_LIBS) $(gsl_LIBS)

# testwlig_SOURCES = testwlig.cc
#	-L$(top_builddir)/high-res     -lcoot-high-res    
#	-I$(top_srcdir)/high-res

# -L$(top_builddir)/ideal        -lcoot-ideal       

#	$(top_builddir)/coords/libcoot-coords.la

coot_findligand_bin_LDADD = \
	./libcoot-ligand.la                                \
	$(top_builddir)/ideal/libcoot-ideal.la             \
        $(top_builddir)/coot-utils/libcoot-map-heavy.la    \
	$(top_builddir)/geometry/libcoot-geometry.la       \
	$(top_builddir)/db-main/libcoot-db-main.la         \
        $(top_builddir)/cootaneer/libcootaneer.la          \
	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la    \
	$(top_builddir)/coot-utils/libcoot-map-utils.la    \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
	$(top_builddir)/utils/libcoot-utils.la             \
        $(top_builddir)/lidia-core/libcoot-lidia-core.la   \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
        $(top_builddir)/analysis/libcoot-analysis.la	   \
	$(top_builddir)/compat/libcoot-compat.la           \
        $(GLOB_LIBS) \
	$(CLIPPER_LIBS) $(MMDB_LIBS) $(gsl_LIBS)

# 	$(top_builddir)/coords/libcoot-coords.la

coot_findwaters_bin_LDADD = \
	./libcoot-ligand.la                                \
	$(top_builddir)/ideal/libcoot-ideal.la             \
        $(top_builddir)/coot-utils/libcoot-map-heavy.la    \
	$(top_builddir)/coot-utils/libcoot-map-utils.la    \
	$(top_builddir)/geometry/libcoot-geometry.la       \
	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la    \
        $(top_builddir)/cootaneer/libcootaneer.la          \
	$(top_builddir)/db-main/libcoot-db-main.la         \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
	$(top_builddir)/utils/libcoot-utils.la             \
        $(top_builddir)/lidia-core/libcoot-lidia-core.la   \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
        $(top_builddir)/analysis/libcoot-analysis.la	   \
	$(top_builddir)/compat/libcoot-compat.la           \
        $(GLOB_LIBS) \
	$(CLIPPER_LIBS) $(MMDB_LIBS) $(gsl_LIBS)

#testwlig_LDADD = -L. -lcoot-ligand      \
#	-L../mini-mol -lcoot-mini-mol    \
#	-L../geometry -lcoot-geometry    \
#	-L../coords -lcoot-coords        \
#	-L../ccp4mg-utils -lccp4mg-utils \
#	$(CLIPPER_LIBS) $(MMDB_LIBS)


#  -L../coords -lcoot-coords -L../atom-utils -lcoot-atomutils

# testori_SOURCES = testori.cc

dist_bin_SCRIPTS = coot-findwaters coot-findligand coot-identify-protein

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

