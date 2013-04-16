# src/Makefile.am
# 
# Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 The University of York
# Author: Paul Emsley
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or (at
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

## Process this file with automake to produce Makefile.in

# and similarly for the GL context:
GTKGL_CONTEXT_CFLAGS = $(GTKGL_CFLAGS) $(GTKGLEXT_CFLAGS)
GTKGL_CONTEXT_LIBS   = $(GTKGL_LIBS) $(GTKGLEXT_LIBS)

# remove this from coot dependencies:
# 	$(top_builddir)/mmtz-extras/libcoot-mmtz-extras.la 

lib_LTLIBRARIES=libccp4mg-opengl.la	\
	libccp4mg-font_cache.la		\
	libccp4mg-cprim.la		\
	libccp4mg-symmetry.la		\
	libccp4mg-ppmutil.la            \
	libccp4mg-psutil.la

include_HEADERS = cbuild.h 	\
	sphere_0.h sphere_1.h sphere_2.h sphere_3.h sphere_4.h 		\
	sphere_arrays.h glew.h wglew.h atom_texture_coords.h

libccp4mg_opengl_la_SOURCES = \
	help.cc splineset.cc texture.cc help.h 	\
	texture.h 
if BUILD_GLEW
libccp4mg_opengl_la_SOURCES+= glew.c glewUtils.cc 
GLEW_FLAG=-DGLEW_BUILD
endif

libccp4mg_opengl_la_LDFLAGS = -L../util -lccp4mg-pygl_coord

libccp4mg_font_cache_la_SOURCES = \
	font_info.cc x11_font.cc win_font.cc freetype_font.cc	\
	font_info.h x11_font.h win_font.h freetype_font.h	\
	font_util.cc freetype_dl.c font_util.h freetype_dl.h

libccp4mg_font_cache_la_LDFLAGS = -L../util -lccp4mg-mgutil

libccp4mg_cprim_la_SOURCES = \
	sphere.cc subdivide.cc cdisplayobject.cc	\
	sphere.h subdivide.h cdisplayobject.h	\
	build_tree_primitives.cc cprimitive.cc text.cc	\
	billboard.cc billboard_primitive.cc cprimitive.h \
	cprimitive_vbo.cc atom_texture_coords.cc

libccp4mg_cprim_la_LDFLAGS = -L../mmut -lccp4mg-mmut -L../util -lccp4mg-mginterrupt -lccp4mg-pygl_coord -lccp4mg-mgutil -lccp4mg-paramsmanager -lccp4mg-rgbreps -lccp4mg-atom_util -lccp4mg-bezier -L. -lccp4mg-opengl -lccp4mg-font_cache

libccp4mg_symmetry_la_SOURCES = \
	symmetry.cc symmetry.h

libccp4mg_symmetry_la_LDFLAGS = -L../util -ccp4mg-pygl_coord -L../mmut -lccp4mg-mmut -L../util -lccp4mg-mginterrupt -lccp4mg-mgutil

libccp4mg_psutil_la_SOURCES = \
	psutil.cc zsortps.cc psutil.h 

libccp4mg_ppmutil_la_SOURCES = \
	ppmutil.cc


# must have GUILE and GUILE_GTK flags before $(INCLUDES), and this 
# does not happen because LTCXXCOMPILE has AM_CXXFLAGS after INCLUDES.
# But DEFAULT_INCLUDES is *before* INCLUDES
# DEFAULT_INCLUDES = ... snip.  Rethink...

# temporary hack FIXME
AM_CXXFLAGS = \
	-DPKGDATADIR='"$(pkgdatadir)"' \
	-DPACKAGE_LOCALE_DIR='"$(pkgdatadir)/locale"' \
	-DDATADIR='"$(datadir)"'    \
	$(GLEW_FLAG)


# SWIG_FLAGS is something of a misnomer.  It actually defines
# COOT_USE_GTK2_INTERFACE, which is needed to correctly declare the
# correct callbacks in both GTK1 and GTK2 regimes.
#
INCLUDES = \
	-I$(top_srcdir)/ccp4mg-utils/util	\
	-I$(top_srcdir)/ccp4mg-utils/mmut	\
	-I$(top_srcdir)/ccp4mg-utils/pygl	\
	-I$(top_srcdir)/ccp4mg-utils	\
	$(GTKGL_CONTEXT_CFLAGS)		\
	$(LIBGLADE_CFLAGS)		\
	@GTK_CFLAGS@ 			\
	$(GLUT_CFLAGS)			\
	$(MMDBSSM_CXXFLAGS) 		\
	$(MMDB_CXXFLAGS) 		\
	$(CLIPPER_CXXFLAGS) 		

#	$(WITH_SCRIPT_CFLAGS) 

