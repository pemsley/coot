# sequence-view/Makefile.am
# 
# Copyright 2003, 2004, 2005 The University of York
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



lib_LTLIBRARIES=libcoot-sequence-view.la

include_HEADERS = sequence-view.hh seq-view-interface.h

CANVAS_CFLAGS = $(GTKCANVAS_CFLAGS) $(GNOME_CANVAS_CFLAGS)
CANVAS_LIBS   = $(GTKCANVAS_LIBS) $(GNOME_CANVAS_LIBS)
GTKGL_CONTEXT_CFLAGS = $(GTKGL_CFLAGS) $(GTKGLEXT_CFLAGS)

# we need these include because we use graphics-info for the
# recentering callback - we don't use c-interface because we dont 
# want the function that returns a sequence_view object to be included
# into c-interface.h (mmdb poisoning).
INCLUDES = \
	-I$(top_srcdir)/compat 		\
	-I$(top_srcdir)/coords 		\
	-I$(top_srcdir)/ccp4mg-utils	\
	-I$(top_srcdir)/coot-utils	\
	-I$(top_srcdir)/intl 		\
	-I$(top_srcdir)/skeleton	\
	-I$(top_srcdir)/angles		\
	-I$(top_srcdir)/atom-utils	\
	-I$(top_srcdir)/db-main		\
	-I$(top_srcdir)/mini-mol	\
	-I$(top_srcdir)/ligand		\
	-I$(top_srcdir)/sequence-view   \
	-I$(top_srcdir)/geometry	\
	-I$(top_srcdir)/ideal   	\
	-I$(top_srcdir)/coot-surface	\
	-I$(top_srcdir)/surface		\
	-I$(top_srcdir)/density-contour	\
	-I$(top_srcdir)/build		\
	-I$(top_srcdir)/high-res	\
	-I$(top_srcdir)/analysis	\
	-I$(top_srcdir)/cootaneer	\
	-I$(top_srcdir)/src     	\
	$(CLIPPER_CXXFLAGS) 		\
	$(MMDB_CXXFLAGS) 		\
	$(MMDBSSM_CXXFLAGS) 		\
	$(GLUT_CFLAGS)			\
	$(CANVAS_CFLAGS) 		\
	$(GTKGL_CONTEXT_CFLAGS)		\
	@GTK_CFLAGS@ 

# CPPFLAGS = -DSEQ_VIEW_STANDALONE
AM_CXXFLAGS = 	-DPKGDATADIR='"$(pkgdatadir)"' \
	$(COOT_GSL_FLAGS) $(GSL_CFLAGS) $(COOT_WINDOWS_CFLAGS)

libcoot_sequence_view_la_SOURCES = sequence-view.cc  seq-view-interface.c

# hmm... still a few undefined symbols from src/ code.
libcoot_sequence_view_la_LIBADD = \
   $(top_builddir)/coot-utils/libcoot-coord-utils.la \
   $(MMDB_LIBS) $(CANVAS_LIBS) $(GTK_LIBS)

# bin_PROGRAMS = testsequenceview

# testsequenceview_SOURCES = testsequenceview.cc

# testsequenceview_LDADD = \
#	-L../coords -lcoot-coords	\
#	-L. -lcoot-sequence-view 	\
#	$(GTKCANVAS_LIBS) 	 	\
#	$(MMDB_LIBS)		 	\
#	@GTK_LIBS@	
