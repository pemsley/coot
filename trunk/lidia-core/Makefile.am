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

include_HEADERS = lbg-shared.hh lbg-molfile.hh lig-build.hh rdkit-interface.hh use-rdkit.hh

#  20110408 I add a CResidue,dictionary_residue_restraints_t constructor to molfile_molecule_t
# simple
#INCLUDES = $(CLIPPER_CXXFLAGS)
#
INCLUDES = \
	 -I$(top_srcdir)/geometry       \
	 -I$(top_srcdir)/mini-mol       \
	 -I$(top_srcdir)/coot-utils     \
	 -I$(top_srcdir)/utils          \
	$(ENTERPRISE_TOOLS_CFLAGS)      \
	$(RDKIT_CXXFLAGS) 		\
        $(CLIPPER_CXXFLAGS) $(MMDB_CXXFLAGS)

libcoot_lidia_core_la_SOURCES = lbg-shared.cc lbg-molfile.cc lig-build.cc rdkit-interface.cc \
	cod-types.cc

libcoot_lidia_core_la_LIBADD =                        \
   $(top_builddir)/geometry/libcoot-geometry.la       \
   $(top_builddir)/utils/libcoot-utils.la             \
   $(RDKIT_LIBS)

# noinst_PROGRAMS = test-molfile

# test_molfile_SOURCES: test-molfile.cc

# test_molfile_LDADD: libcoot-lidia-core.la

noinst_PROGRAMS = test-cod-atom-types
test_cod_atom_types_SOURCES = test-cod-atom-types.cc

test_cod_atom_types_LDADD = libcoot-lidia-core.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la 

