# Makefile.am
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

# Tell autoreconf where to find macros
ACLOCAL_AMFLAGS=-I macros

# see configure.in for po intl 
SUBDIRS = auxil scheme python macros compat coot-utils surface \
	skeleton ccp4mg-utils ccp4mg-utils/util ccp4mg-utils/mmut \
	ccp4mg-utils/pygl ccp4mg-utils/mgapp coot-surface density-contour \
	coords angles mini-mol db-main  geometry  \
	ideal ligand high-res build greg-tests python-tests \
	analysis sequence-view tw atom-utils cootaneer src rama-data doc 

# po, no po today.

cootrc: $(top_srcdir)/cootrc.template
	$(MAKE_COOTRC_CMD) $(top_srcdir)/cootrc.template > $@


EXTRA_DIST = 		\
	autogen.sh 	\
	make-scheme-tar \
	mapview.glade 	\
	cootrc.template	\
	standard-residues.pdb \
	BUGS		\
	depcomp 	\
	libtool 	\
	config.sub 	\
	config.guess 	\
	install-sh 	\
	missing 	\
	TODO 		\
	count-loc 	\
	syminfo.lib	\
	colours.def	\
	function-notes	\
	my.supp		\
	guile.supp	\
	make-patches	\
	ifdef-for-create-aboutdialog.awk \
	mtrix-to-ncs-matrix.awk  \
	theor-helix-70.pdb       \
	RELEASE-NOTES	         \
	copyright-disclaimer.txt \
	intltool-extract.in      \
	intltool-merge.in        \
	intltool-update.in 

install-data-local:
	@$(NORMAL_INSTALL)
	if test -d $(srcdir)/pixmaps; then \
	  $(mkinstalldirs) $(DESTDIR)$(pkgdatadir)/pixmaps; \
	  for pixmap in $(srcdir)/pixmaps/*; do \
	    if test -f $$pixmap; then \
	      $(INSTALL_DATA) $$pixmap $(DESTDIR)$(pkgdatadir)/pixmaps; \
	    fi \
	  done \
	fi

dist-hook:
	if test -d $(srcdir)/pixmaps; then \
	  mkdir $(distdir)/pixmaps; \
	  for pixmap in $(srcdir)/pixmaps/*; do \
	    if test -f $$pixmap; then \
	      cp -p $$pixmap $(distdir)/pixmaps; \
	    fi \
	  done \
	fi

subpkgdatadir = $(pkgdatadir)
subpkgdata_DATA = cootrc standard-residues.pdb syminfo.lib colours.def \
	mtrix-to-ncs-matrix.awk theor-helix-70.pdb


