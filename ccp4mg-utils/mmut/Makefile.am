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

lib_LTLIBRARIES=libccp4mg-mmut.la

include_HEADERS=splineinfo.h

libccp4mg_mmut_la_SOURCES = \
	mmut_manager.cc mmut_util.cc mman_manager.cc	\
	mmut_manager.h mmut_util.h mman_manager.h	\
	mman_base.cc mmut_secstr.cc mmut_bonds.cc	\
	mman_base.h mmut_secstr.h mmut_bonds.h		\
	mmut_sbase.cc mmut_sasarea.cc mmut_connectivity.cc	\
	mmut_sbase.h mmut_sasarea.h mmut_connectivity.h	\
	mmut_hbond.cc mmut_contact.cc mmut_basepairs.cc	\
	mmut_hbond.h mmut_contact.h mmut_basepairs.h	\
	mmut_spline.cc mmut_lipids.cc mmut_morph.cc	\
	mmut_lipids.h mmut_morph.h	\
	mmut_nma.cc mmut_colour.cc mmut_nma.h mmut_colour.h


# must have GUILE and GUILE_GTK flags before $(INCLUDES), and this 
# does not happen because LTCXXCOMPILE has AM_CXXFLAGS after INCLUDES.
# But DEFAULT_INCLUDES is *before* INCLUDES
# DEFAULT_INCLUDES = ... snip.  Rethink...

AM_CXXFLAGS = \
	-DPKGDATADIR='"$(pkgdatadir)"' \
	-DPACKAGE_LOCALE_DIR='"$(pkgdatadir)/locale"' \
	-DDATADIR='"$(datadir)"'    \
	$(COOT_WINDOWS_CFLAGS) 


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

