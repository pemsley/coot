# Makefile.am
# 
# Copyright 2002  The University of York
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

AM_CPPFLAGS = -I$(top_srcdir) $(CLIPPER_CXXFLAGS) 

lib_LTLIBRARIES = libcoot-isosurface.la

pkginclude_HEADERS = Vectors.h CIsoSurface.h density-contour-triangles.hh occlusion.hh transfer-occlusions.hh
pkgincludedir = $(includedir)/coot/density-contour

check_PROGRAMS = test-occlusion

libcoot_isosurface_la_SOURCES = Vectors.cpp CIsoSurface.cpp density-contour-triangles.cc occlusion.cc

#libcoot_isosurface_la_LIBADD = \
#	$(top_builddir)/coords/libcoot-coords.la \
#        $(CLIPPER_LIBS)

libcoot_isosurface_la_LIBADD = $(top_builddir)/coords/libcoot-coords.la    \
                               $(CLIPPER_LIBS)

libcoot_isosurface_la_LDFLAGS = $(SHARED_LDFLAGS)

test_occlusion_SOURCES = test-occlusion.cc
# I need to add $(CLIPPER_LIBS) here, even though libcoot-coord-utils is correctly linked because
# I use clipper functions in the test-occlusion program
test_occlusion_LDADD =  libcoot-isosurface.la $(top_builddir)/coot-utils/libcoot-coord-utils.la $(CLIPPER_LIBS)


