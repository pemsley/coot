# skeleton/Makefile.am
# 
# Copyright 2001, 2002, 2003 The University of York
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


lib_LTLIBRARIES=libcoot-skeleton.la

INCLUDES =-I$(top_srcdir)/coords -I$(top_srcdir)/angles 	\
	 -I$(top_srcdir)/atom-utils  -I$(top_srcdir)/compat	\
	 $(MMDB_CXXFLAGS) $(CLIPPER_CXXFLAGS) 

AM_CXXFLAGS = $(COOT_WINDOWS_CFLAGS)

include_HEADERS = graphical_skel.h dm-utils.h BuildCas.h

libskeleton_a_sources =  \
	dm-utils.cc \
	graphical_skel.cc \
	BuildCas.cc

# No more fc_skeleton
#	fc_skeleton.cc fc_skeleton.h


libcoot_skeleton_la_SOURCES = $(libskeleton_a_sources)


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

