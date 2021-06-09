# src/Makefile.am
# 
# Copyright 2011, 2012 by The University of Oxford
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


AM_CPPFLAGS = -I$(top_srcdir) \
     $(GOOCANVAS_CFLAGS) $(GTK_CFLAGS)        

lib_LTLIBRARIES=libcoot-goograph.la 

libcoot_goograph_la_SOURCES = \
    goograph.cc goograph.hh

libcoot_goograph_la_LIBADD = $(top_builddir)/lidia-core/libcoot-lidia-core.la \
			     $(top_builddir)/utils/libcoot-utils.la \
			     $(GOOCANVAS_LIBS) $(GTK_LIBS)

libcoot_goograph_la_LDFLAGS = $(SHARED_LDFLAGS)

check_PROGRAMS = test-goograph

test_goograph_SOURCES = test-goograph.cc

test_goograph_LDADD = libcoot-goograph.la \
	$(top_builddir)/lidia-core/libcoot-lidia-core.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/utils/libcoot-utils.la \
        $(GTK_LIBS)
