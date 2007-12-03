# Makefile.am
# 
# Copyright 2005, 2006, 2007 The University of York
# Copyright 2005, 2006, 2007 Bernhard Lohkamp
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

pythondir = $(pkgdatadir)/python

python_DATA =  $(EXTRA_DIST) coot_load_modules.py

if COOT_USE_PYTHON
TEXI_FILE =  coot-python-functions.texi
coot-python-functions.texi:
	python python_coot_docs.py 
endif

EXTRA_DIST =    \
   Makefile.am  \
   americanisms.py       fitting.py    refmac.py      \
   background_demo.py    gap.py        shelx.py       \
   get_ebi.py            mutate.py     coot_lsq.py    \
   hello.py              povray.py     coot_utils.py  \
   raster3d.py           extensions.py generic_objects.py \
   ncs.py                tips_gui.py   tips.py          \
   coot_load_modules.py  filter.py     local_code.py	\
   coot_gui.py           libcheck.py   cns2coot.py      \
   group_settings.py     ../src/coot.py $(TEXI_FILE)

