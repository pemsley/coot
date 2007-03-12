## Process this file with automake to produce Makefile.in

# Tell autoreconf where to find macros
ACLOCAL_AMFLAGS=-I macros

# see configure.in for po intl 
SUBDIRS = scheme python macros compat coot-utils surface coot-surface \
	skeleton ccp4mg-utils density-contour coords angles \
	mini-mol db-main  geometry  \
	ideal ligand high-res build \
	analysis sequence-view tw atom-utils src doc

EXTRA_DIST = 		\
	autogen.sh 	\
	make-scheme-tar \
	mapview.glade 	\
	cootrc		\
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
	mtrix-to-ncs-matrix.awk \
	theor-helix-70.pdb \
	RELEASE-NOTES	\
	copyright-disclaimer.txt

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


