# Makefile.am
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

lib_LTLIBRARIES=libcoot-lidia-core.la

pkginclude_HEADERS = lbg-shared.hh lbg-molfile.hh lig-build.hh \
	rdkit-interface.hh use-rdkit.hh cod-types.hh neighbour-sorter.hh \
	bond-record-container-t.hh bond-table-record-t.hh \
	cod-atom-types.hh cod-atom-type-t.hh primes.hh \
	chemical-feature-clusters.hh third-neighbour-info-t.hh get-residue.hh \
	svg-molecule.hh svg-container.hh atom-types.hh

pkgincludedir = $(includedir)/coot/lidia-core


AM_CPPFLAGS = \
	 -I$(top_srcdir) \
	$(ENHANCED_LIGAND_TOOLS_CPPFLAGS) \
	$(RDKIT_CXXFLAGS) 		  \
	$(BOOST_CPPFLAGS) 	          \
	$(SQLITE3_CFLAGS)                 \
        $(CLIPPER_CXXFLAGS) $(MMDB_CXXFLAGS) $(COOT_WINDOWS_CFLAGS) $(BOOST_CPPFLAGS)

libcoot_lidia_core_la_SOURCES = lbg-shared.cc lbg-molfile.cc lig-build.cc \
	rdkit-interface.cc cod-types.cc atom-types.cc bond-record-container-t.cc \
	cod-atom-type-t.cc cod-atom-types.cc bond-table-record-t.cc \
	chemical-feature-clusters.cc get-residue.cc svg-molecule.cc

# Are these needed?
#   $(top_builddir)/coot-utils/libcoot-coord-utils.la  
#   $(top_builddir)/utils/libcoot-utils.la             
#
# bond-record-container-t has sqlite3 code in it (to create the atoms database currently)
libcoot_lidia_core_la_LIBADD =                    \
   $(top_builddir)/geometry/libcoot-geometry.la   \
   $(top_builddir)/utils/libcoot-utils.la         \
   $(top_builddir)/mini-mol/libcoot-mini-mol.la   \
   $(CLIPPER_LIBS)                                \
   $(RDKIT_LIBS)                                  \
   $(BOOST_LDFLAGS)                               \
   $(SQLITE3_LIBS)

#   -lRDKitmaeparser

# interesting link failure when using -no-undefined.
#
# libcoot_lidia_core_la_LDFLAGS = $(SHARED_LDFLAGS) -no-undefined
libcoot_lidia_core_la_LDFLAGS = $(SHARED_LDFLAGS)

# noinst_PROGRAMS = test-molfile

# test_molfile_SOURCES: test-molfile.cc

# test_molfile_LDADD: libcoot-lidia-core.la

check_PROGRAMS = test-lig-build

test_lig_build_SOURCES = test-lig-build.cc

test_lig_build_CXXFLAGS = $(CLIPPER_CXXFLAGS)

test_lig_build_LDFLAGS = -lcoot-lidia-core      \
   $(top_builddir)/geometry/libcoot-geometry.la \
   $(top_builddir)/utils/libcoot-utils.la       \
   $(CLIPPER_LIBS)                              \
   $(RDKIT_LIBS)

