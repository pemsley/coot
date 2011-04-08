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

include_HEADERS = lbg-shared.hh lbg-molfile.hh lig-build.hh

#  20110408 I add a CResidue,dictionary_residue_restraints_t constructor to molfile_molecule_t
# simple
#INCLUDES = $(CLIPPER_CXXFLAGS)
#
INCLUDES = \
	 -I$(top_srcdir)/geometry \
	 $(CLIPPER_CXXFLAGS) $(MMDB_CXXFLAGS)

libcoot_lidia_core_la_SOURCES = lbg-molfile.cc lig-build.cc

