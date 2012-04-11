# Makefile.am
# 
# Copyright 2005, 2006, 2007 The University of York
# Copyright 2005, 2006, 2007 Bernhard Lohkamp
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

pythondir = $(pkgdatadir)/python


# how are these generated?  AFAICS, they are not.  comment out for now
# 
#if COOT_USE_PYTHON
#TEXI_FILE =  coot-python-functions.texi


#coot-python-functions.texi:
#	python $(top_srcdir)/python/python_coot_docs.py 

#endif

dist_python_DATA = \
                        coot_utils.py        get_ebi.py      mutate.py        \
   background_demo.py   hello.py             raster3d.py                      \
   brute_lsqman.py                                           refmac.py        \
   coot_gui.py          filter.py            libcheck.py                      \
   coot_load_modules.py                      local_code.py                    \
   gap.py               remote_control.py    povray.py       coot_lsq.py      \
   shelx.py             group_settings.py    fitting.py      americanisms.py  \
   tips.py              tips_gui.py          generic_objects.py  cns2coot.py  \
   extensions.py        ncs.py               create_server.py                 \
   python_coot_docs.py  clear_backup.py      bones_reader.py                  \
   shelx_extensions.py  parse_pisa_xml.py    		        		          \
   check_for_updates.py prodrg_import.py     								  \
   get_recent_pdbe.py   jligand.py	     	 jligand_gui.py	              \
   coot_toolbuttons.py              \
   $(TEXI_FILE)

