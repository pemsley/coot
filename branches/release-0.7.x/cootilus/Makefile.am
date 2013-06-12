# Makefile.am
# 
# Copyright 2007 by Paul Emsley
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

EXTRA_DIST = 		\
	nautilus_lib.pdb

INCLUDES = 	$(CLIPPER_CXXFLAGS) $(MMDB_CXXFLAGS) 

cootilus_demo_SOURCES = cootilus-demo.cpp

cootilus_demo_LDADD = ./libcootilus.la $(CLIPPER_LIBS) $(MMDB_LIBS) 

noinst_PROGRAMS = cootilus-demo

lib_LTLIBRARIES=libcootilus.la

libcootilus_la_SOURCES = \
	nucleicacid_db.cpp nautilus-tools.cpp nautilus-ss-find.cpp \
	nautilus-join.cpp nautilus-target.cpp cootilus-build.cpp \
	nucleicacid_db.h nautilus-tools.h nautilus-ss-find.h \
	nautilus-join.h nautilus-target.h cootilus-build.h

libcootilus_la_LIBADD = $(CLIPPER_LIBS)

subpkgdatadir = $(pkgdatadir)
subpkgdata_DATA = nautilus_lib.pdb
