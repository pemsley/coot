# protein_db/Makefile.am
# 
# Copyright 2011 The University of Oxford
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

lib_LTLIBRARIES=libcoot-protein-db.la

INCLUDES= $(CLIPPER_CXXFLAGS)  $(MMDB_CXXFLAGS)

include_HEADERS=protein_db.h  protein_db_utils.h rna_db.h

libcoot_protein_db_la_SOURCES=protein_db.cpp  protein_db_utils.cpp rna_db.cpp

libcoot_protein_db_la_LIBADD = $(CLIPPER_LIBS) $(MMDB_LIBS)

EXTRA_DIST=protein.db

protein_db_datadir=$(pkgdatadir)/protein_db

protein_db_data_DATA=$(EXTRA_DIST)

