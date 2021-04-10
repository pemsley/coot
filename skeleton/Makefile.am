# skeleton/Makefile.am
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


lib_LTLIBRARIES=libcoot-skeleton.la

AM_CPPFLAGS = \
	  -I$(top_srcdir) \
	 $(MMDB_CXXFLAGS) $(CLIPPER_CXXFLAGS) 

AM_CXXFLAGS = $(COOT_WINDOWS_CFLAGS)

pkginclude_HEADERS = graphical_skel.h dm-utils.h BuildCas.h
pkgincludedir = $(includedir)/coot/skeleton

libcoot_skeleton_la_SOURCES =  \
	dm-utils.cc \
	graphical_skel.cc \
	BuildCas.cc

libcoot_skeleton_la_LIBADD = \
   $(top_builddir)/angles/libcoot-AngleInfo.la  \
   $(top_builddir)/ideal/libcoot-ideal.la       \
   $(top_builddir)/geometry/libcoot-geometry.la \
   $(top_builddir)/coords/libcoot-coords.la     \
   $(CLIPPER_LIBS) $(MMDB_LIBS)

libcoot_skeleton_la_LDFLAGS = $(SHARED_LDFLAGS)



# bin_PROGRAMS = ncssearch test_various test_torsion test_theta_2

# ncssearch_SOURCES = ncs-search.cc ncs-search.h

# test_various_SOURCES = test_various.cc 

# test_torsion_SOURCES = test-torsion.cc

# test_theta_2_SOURCES = test_theta_2.cc

# ncssearch_LDADD = -L. -lcoot-skeleton \
# 	-L../coords -lcoot-coords	\
# 	-L../angles -lcoot-AngleInfo \
# 	$(CLIPPER_LIBS) 	\
# 	$(MMDB_LIBS)

# test_torsion_LDADD = -L. -lcoot-skeleton 	\
# 		-L../coords -lcoot-coords	\
# 		-L../angles -lcoot-AngleInfo \
# 		$(CLIPPER_LIBS) 	\
# 		$(MMDB_LIBS)

# test_theta_2_LDADD = -L../angles -lcoot-AngleInfo

# sktest_SOURCES = sktest.cpp
# sktest_LDADD = -L. -lskeleton $(CLIPPER_LIBS) 

# test_various_LDADD = -L. -lcoot-skeleton \
# 	-L../coords -lcoot-coords	\
# 	-L../angles -lcoot-AngleInfo \
# 	$(CLIPPER_LIBS) 	\
# 	$(MMDB_LIBS)

