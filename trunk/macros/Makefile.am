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

## Please update this variable if any new macros are created

# 20091111 Actually, probably not necessary, neither MACROS or
# EXTRA_DIST because there are rules when Makefile (or Makefile.in is
# built that finds .m4 files and puts them in the
# dist). c.f. clipper's config/Makefile.am

MACROS=						\
  aclocal-include.m4				\
  compiler-flags.m4				\
  gtkgl.m4					\
  gtk-canvas.m4					\
  clipper.m4					\
  clipper.m4.old				\
  mmdb.m4					\
  mmdb-ssm.m4					\
  with-guile.m4					\
  with-python.m4				\
  guile.m4					\
  coot-windows.m4				\
  lib-prefix.m4					\
  glut.m4				

EXTRA_DIST=$(MACROS)
