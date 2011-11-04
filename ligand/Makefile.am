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

INCLUDES = \
	-I$(top_srcdir)/ccp4mg-utils \
	-I$(top_srcdir)/mini-mol \
	-I$(top_srcdir)/geometry \
	-I$(top_srcdir)/ideal \
	-I$(top_srcdir)/coords \
	-I$(top_srcdir)/lidia-core \
	-I$(top_srcdir)/utils \
	-I$(top_srcdir)/coot-utils \
	-I$(top_srcdir)/compat \
	-I$(top_srcdir)/db-main \
	$(COOT_GSL_FLAGS)		\
	$(GSL_CFLAGS)			\
	$(CLIPPER_CXXFLAGS) $(MMDB_CXXFLAGS) 

include_HEADERS = ligand.hh residue_by_phi_psi.hh wligand.hh monomer-utils.hh \
	chi-angles.hh dunbrack.hh helix-placement.hh fast-ss-search.hh        \
	rigid-body.hh ideal-rna.hh rotamer.hh richardson-rotamer.hh           \
	torsion-general.hh primitive-chi-angles.hh dipole.hh    \
	base-pairing.hh backrub-rotamer.hh

libcoot_ligand_la_SOURCES = ligand.cc ligand-extras.cc  residue_by_phi_psi.cc \
	wiggly-ligand.cc monomer-utils.cc dunbrack.cc chi-angles-autogen.cc   \
	chi-angles.cc helix-placement.cc fast-ss-search.cc rigid-body.cc      \
	ideal-rna.cc richardson-rotamers.cc richardson-rotamer.cc rotamer.cc  \
	torsion-general.cc rotamer-tables.cc primitive-chi-angles.cc dipole.cc\
        base-pairing.cc backrub-rotamer.cc

# clean
libcoot_ligand_la_LIBADD = \
   $(top_builddir)/ideal/libcoot-ideal.la             \
   $(top_builddir)/db-main/libcoot-db-main.la         \
   $(top_builddir)/ccp4mg-utils/libccp4mg-utils.la    \
   $(top_builddir)/geometry/libcoot-geometry.la       \
   $(top_builddir)/mini-mol/libcoot-mini-mol.la       \
   $(top_builddir)/coords/libcoot-coords.la           \
   $(top_builddir)/coot-utils/libcoot-coord-extras.la \
   $(top_builddir)/coot-utils/libcoot-coord-utils.la  \
   $(top_builddir)/coot-utils/libcoot-map-utils.la    \
   $(top_builddir)/utils/libcoot-utils.la             \
   $(CLIPPER_LIBS) $(MMDB_LIBS)

AM_CXXFLAGS = \
	-DPKGDATADIR='"$(pkgdatadir)"' \
	$(COOT_WINDOWS_CFLAGS)

bin_PROGRAMS = findligand-real findwaters-real

# check_PROGRAMS = test-helix-placement test-fffear-in-coot testdunbrack test-random test-ligand

findligand_real_SOURCES = find-ligand.cc
findwaters_real_SOURCES = find-waters.cc

# test_ligand_SOURCES = test-ligand.cc

# test_helix_placement_SOURCES = test-helix-placement.cc

# test_random_SOURCES = test-random.cc

#test_random_LDADD = \
#	./libcoot-ligand.la                                \
#	$(top_builddir)/ideal/libcoot-ideal.la             \
#	$(top_builddir)/coords/libcoot-coords.la           \
#	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
#	$(top_builddir)/coot-utils/libcoot-map-utils.la    \
#	$(top_builddir)/db-main/libcoot-db-main.la         \
#	$(top_builddir)/coot-utils/libcoot-shelx.la        \
#	$(top_builddir)/geometry/libcoot-geometry.la       \
#	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
#	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
#	$(top_builddir)/utils/libcoot-utils.la             \
#	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la    \
#	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
#	$(top_builddir)/compat/libcoot-compat.la           \
#	$(CLIPPER_LIBS) $(MMDB_LIBS) $(GSL_LIBS)

#test_ligand_LDADD = \
#	./libcoot-ligand.la                                \
#	$(top_builddir)/ideal/libcoot-ideal.la             \
#	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
#	$(top_builddir)/coords/libcoot-coords.la           \
#	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
#	$(top_builddir)/coot-utils/libcoot-map-utils.la    \
#	$(top_builddir)/db-main/libcoot-db-main.la         \
#	$(top_builddir)/coot-utils/libcoot-shelx.la        \
#	$(top_builddir)/geometry/libcoot-geometry.la       \
#	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
#	$(top_builddir)/utils/libcoot-utils.la             \
#	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
#	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la    \
#	$(top_builddir)/compat/libcoot-compat.la           \
#	$(CLIPPER_LIBS) $(MMDB_LIBS) $(GSL_LIBS)

#test_helix_placement_LDADD = \
#	./libcoot-ligand.la                                \
#	$(top_builddir)/ideal/libcoot-ideal.la             \
#	$(top_builddir)/coords/libcoot-coords.la           \
#	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
#	$(top_builddir)/coot-utils/libcoot-map-utils.la    \
#	$(top_builddir)/coot-utils/libcoot-shelx.la        \
#	$(top_builddir)/geometry/libcoot-geometry.la       \
#	$(top_builddir)/db-main/libcoot-db-main.la         \
#	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
#	$(top_builddir)/utils/libcoot-utils.la             \
#	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
#	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la    \
#	$(top_builddir)/compat/libcoot-compat.la           \
#	$(CLIPPER_LIBS) $(MMDB_LIBS) $(GSL_LIBS)

#test_fffear_in_coot_SOURCES = test-fffear-in-coot.cc

#test_fffear_in_coot_LDADD = \
#	./libcoot-ligand.la \
#	$(top_builddir)/ideal/libcoot-ideal.la             \
#	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
#	$(top_builddir)/coords/libcoot-coords.la           \
#	$(top_builddir)/coot-utils/libcoot-shelx.la        \
#	$(top_builddir)/geometry/libcoot-geometry.la       \
#	$(top_builddir)/coot-utils/libcoot-map-heavy.la    \
#	$(top_builddir)/db-main/libcoot-db-main.la         \
#	$(top_builddir)/coot-utils/libcoot-map-utils.la    \
#	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
#	$(top_builddir)/utils/libcoot-utils.la             \
#	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
#	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la    \
#	$(top_builddir)/compat/libcoot-compat.la           \
#	$(CLIPPER_LIBS) $(MMDB_LIBS) $(GSL_LIBS)

# testwlig_SOURCES = testwlig.cc
#	-L$(top_builddir)/high-res     -lcoot-high-res    
#	-I$(top_srcdir)/high-res

# -L$(top_builddir)/ideal        -lcoot-ideal       

findligand_real_LDADD = \
	./libcoot-ligand.la                                \
	$(top_builddir)/ideal/libcoot-ideal.la             \
	$(top_builddir)/geometry/libcoot-geometry.la       \
	$(top_builddir)/coords/libcoot-coords.la           \
	$(top_builddir)/db-main/libcoot-db-main.la         \
	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la    \
	$(top_builddir)/coot-utils/libcoot-map-utils.la    \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
	$(top_builddir)/utils/libcoot-utils.la             \
        $(top_builddir)/lidia-core/libcoot-lidia-core.la   \
	$(top_builddir)/coot-utils/libcoot-shelx.la        \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
	$(top_builddir)/compat/libcoot-compat.la           \
	$(CLIPPER_LIBS) $(MMDB_LIBS) $(GSL_LIBS)

findwaters_real_LDADD = \
	./libcoot-ligand.la                                \
	$(top_builddir)/ideal/libcoot-ideal.la             \
	$(top_builddir)/geometry/libcoot-geometry.la       \
	$(top_builddir)/coords/libcoot-coords.la           \
	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la    \
	$(top_builddir)/db-main/libcoot-db-main.la         \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
	$(top_builddir)/coot-utils/libcoot-map-utils.la    \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
	$(top_builddir)/utils/libcoot-utils.la             \
        $(top_builddir)/lidia-core/libcoot-lidia-core.la   \
	$(top_builddir)/coot-utils/libcoot-shelx.la        \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
	$(top_builddir)/compat/libcoot-compat.la           \
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
#	$(top_builddir)/coot-utils/libcoot-shelx.la        \
#	$(top_builddir)/mini-mol/libcoot-mini-mol.la       \
#        $(CLIPPER_LIBS) $(MMDB_LIBS) $(GSL_LIBS)

