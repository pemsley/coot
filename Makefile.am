# Makefile.am
# 
# Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 The University of York
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

## Process this file with automake to produce Makefile.in

# Tell autoreconf where to find macros
ACLOCAL_AMFLAGS=-I macros

# see configure.in for po intl 
#
# coot-utils depends on mini-mol.
# geometry depends on utils
#

# also: coords depends on ligand because now we do rotamer markup.
#       so ligand functions must not use coot::Cartesian
if COND_LIDIA_ONLY
   SUBDIRS = compat utils protein_db mini-mol geometry lidia-core coot-utils lbg
else
if COND_LIBCOOTAPI_ONLY
   SUBDIRS = auxil scheme python pdb-templates macros compat utils protein_db mini-mol geometry  \
        lidia-core coot-utils angles ccp4mg-utils db-main analysis cootaneer ideal ligand coords \
	density-contour skeleton  high-res build MoleculesToTriangles pli \
	python-tests cootilus fonts pyrogen lbg data \
	man rama-data glTF pixmaps textures api
else
   # rcrane is not included yet.
   SUBDIRS = auxil scheme python pdb-templates macros compat utils protein_db mini-mol geometry  \
        lidia-core coot-utils angles ccp4mg-utils db-main analysis cootaneer ideal ligand coords \
	density-contour skeleton high-res build MoleculesToTriangles pli api greg-tests \
	python-tests cootilus fonts pyrogen lbg data \
	man rama-data doc glTF pixmaps textures validation-graphs src
endif
endif


# po, no po today.

cootrc: $(top_srcdir)/cootrc.template
	$(MAKE_COOTRC_CMD) $(top_srcdir)/cootrc.template > $@


EXTRA_DIST = 		\
	autogen.sh 	\
	make-scheme-tar \
	cootrc.template	\
	py-compile      \
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
	a6.glade        \
	coot-gtk4.ui    \
	single-map-properties-dialog.glade \
	syminfo.lib	\
	colours.def	\
	my.supp		\
	guile.supp	\
	windows/syminfo.lib	\
	fix-cif-loops.awk       \
	ifdef-for-create-aboutdialog.awk \
	mtrix-to-ncs-matrix.awk  \
	RELEASE-NOTES	         \
	copyright-disclaimer.txt \
	intltool-extract.in      \
	intltool-merge.in        \
	intltool-update.in 

# install-data-local:
# 	@$(NORMAL_INSTALL)
# 	if test -d $(srcdir)/pixmaps; then \
# 	  $(mkinstalldirs) $(DESTDIR)$(pkgdatadir)/pixmaps; \
# 	  for pixmap in $(srcdir)/pixmaps/*; do \
# 	    if test -f $$pixmap; then \
# 	      $(INSTALL_DATA) $$pixmap $(DESTDIR)$(pkgdatadir)/pixmaps; \
# 	    fi \
# 	  done \
# 	fi
# 	if test -z "$(COOT_WINDOWS_CFLAGS)" ; then \
# 	  $(INSTALL_DATA) $(srcdir)/syminfo.lib $(DESTDIR)$(pkgdatadir); \
# 	  else \
# 	  $(INSTALL_DATA) $(srcdir)/windows/syminfo.lib $(DESTDIR)$(pkgdatadir); \
# 	fi

# 20140523 Surely this is not needed now?  I have a Makefile in pixmaps.
#
#dist-hook:
#	if test -d $(srcdir)/pixmaps; then \
#	  mkdir $(distdir)/pixmaps; \
#	  for pixmap in $(srcdir)/pixmaps/*; do \
#	    if test -f $$pixmap; then \
#	      cp -p $$pixmap $(distdir)/pixmaps; \
#	    fi \
#	  done \
#	fi

# packdatadir is xxx/share/coot
subpkgdatadir = $(pkgdatadir)
subpkgdata_DATA = cootrc standard-residues.pdb colours.def \
	mtrix-to-ncs-matrix.awk syminfo.lib

cootgladedir = $(pkgdatadir)/glade
dist_cootglade_DATA = coot-gtk4.ui a6.glade preferences-gtk3.glade single-map-properties-dialog-gtk4.ui

appdatadir = ${datarootdir}/appdata
dist_appdata_DATA = coot.appdata.xml

applicationdir = ${datarootdir}/applications
dist_application_DATA = coot.desktop


