# tw/Makefile.am
# 
# Copyright 2001, 2002, 2003 The University of York
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

# Do it properly with configure (not like this)
# OGG_VORBIS_CFLAGS=-I/sw/include/vorbis -I/usr/include/vorbis
# AO_CFLAGS=-I/sw/include -I/usr/include

AM_CXXFLAGS = \
	-DPKGDATADIR='"$(pkgdatadir)"'

AM_CPPFLAGS = \
	$(GTKCANVAS_CFLAGS)  \
	$(GTK_CFLAGS)	     \
	$(OGG_VORBIS_CFLAGS) \
	$(AO_CFLAGS)

# don't make this yet, just carry it around.
EXTRA_DIST = cootsound.c

lib_LTLIBRARIES=libcoot-tw.la

pkginclude_HEADERS = coot-tw.hh
pkgincludedir = $(includedir)/coot/tw

libcoot_tw_la_SOURCES = tw.cc

libcoot_tw_la_LDFLAGS = $(SHARED_LDFLAGS)

check_PROGRAMS = coot-tiddly

coot_tiddly_SOURCES = tiddly.cc

coot_tiddly_LDADD = libcoot-tw.la \
	$(GTKCANVAS_LIBS) \
	@GTK_LIBS@


# cootsound_LDADD = -L/sw/lib -lvorbisfile -logg -lao

