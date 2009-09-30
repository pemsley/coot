# Makefile.am
# 
# Copyright 2002 The University of York
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

INCLUDES = -I$(top_srcdir)/surface \
	$(MMDB_CXXFLAGS) $(CLIPPER_CXXFLAGS) $(GTKGL_CFLAGS) $(GLUT_CFLAGS)

lib_LTLIBRARIES = libcoot-surface.la

libcoot_surface_la_SOURCES = \
	 coot-surface.cc coot-surface.hh rgbreps.cc rgbreps.h

libcoot_surface_la_LIBADD = \
	$(top_builddir)/surface/libccp4mg-surface.la \
	$(GLUT_LIBS) $(CLIPPER_LIBS) $(MMDB_LIBS) 
