# src/Makefile.am
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

if COOT_USE_PYTHON
pythondir = $(pkgdatadir)/python
python_SCRIPTS = coot.py
endif


######## 20091001 ########################################################


# COOT_STDTYPES is something like: DHAVE_SYS_STDTYPES_H=1

# USE_GUILE_DIRECTIVE=`if [ -z "$(COOT_USE_GUILE)"]; then echo ""; else echo -DUSE_GUILE; fi`

WITH_GUILE_CFLAGS=$(COOT_STDTYPES) $(COOT_USE_GUILE)


# GUILE_VERSION=$(GUILE_CONFIG)

MAKE_COOT_WRAP_GUILE=$(SHELL) $(top_srcdir)/src/post-process-coot-wrap-guile.sh "$(GUILE_CONFIG)"

MAKE_COOT_WRAP_PYTHON_XXX=$(SHELL) $(top_srcdir)/src/post-process-coot-wrap-python.sh 


##########################################################################

CANVAS_CFLAGS = $(GTKCANVAS_CFLAGS) $(GNOME_CANVAS_CFLAGS)
CANVAS_LIBS   = $(GTKCANVAS_LIBS) $(GNOME_CANVAS_LIBS)

# and similarly for the GL context:
GTKGL_CONTEXT_CFLAGS = $(GTKGL_CFLAGS) $(GTKGLEXT_CFLAGS)
GTKGL_CONTEXT_LIBS   = $(GTKGL_LIBS) $(GTKGLEXT_LIBS)

# remove this from coot dependencies:
# 	$(top_builddir)/mmtz-extras/libcoot-mmtz-extras.la 

COOT_EXTRA_LIBS = \
	$(top_builddir)/ligand/libcoot-ligand.la      		\
	$(top_builddir)/lidia-core/libcoot-lidia-core.la      	\
	$(top_builddir)/density-contour/libcoot-isosurface.la  	\
	$(top_builddir)/build/libcoot-cabuild.la		\
	$(top_builddir)/skeleton/libcoot-skeleton.la    	\
	$(top_builddir)/angles/libcoot-AngleInfo.la   		\
	$(top_builddir)/db-main/libcoot-db-main.la     		\
	$(top_builddir)/high-res/libcoot-high-res.la  		\
	$(top_builddir)/ideal/libcoot-ideal.la    		\
	$(top_builddir)/geometry/libcoot-geometry.la    	\
	$(top_builddir)/coords/libcoot-coords.la      		\
	$(top_builddir)/coot-surface/libcoot-surface.la  	\
	$(top_builddir)/surface/libccp4mg-surface.la  		\
	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la 	\
	$(top_builddir)/coot-utils/libcoot-map-heavy.la  	\
	$(top_builddir)/coot-utils/libcoot-map-utils.la  	\
	$(top_builddir)/coot-utils/libcoot-coord-extras.la  	\
	$(top_builddir)/coot-utils/libcoot-coord-utils.la  	\
	$(top_builddir)/coot-utils/libcoot-shelx.la	  	\
	$(top_builddir)/utils/libcoot-utils.la  		\
	$(top_builddir)/mini-mol/libcoot-mini-mol.la    	\
	$(top_builddir)/analysis/libcoot-analysis.la 		\
	$(top_builddir)/cootaneer/libcootaneer.la 		\
	$(top_builddir)/compat/libcoot-compat.la  		



# Remember to regenerate coot_wrap_python_pre.cc before making dist, eh?
#
EXTRA_DIST = \
	coot.i post-glade add-c-extern.sh blank.cc fix-gettext.sh \
	mol-hierachy.cc \
	coot.py coot_gtk1.py coot_gtk2.py \
	junk.cc coot.defs coot-dox.cfg post-gtk2-glade.sh \
        post-process-coot-wrap-guile.sh \
        post-process-coot-wrap-python.sh



# use -DHAVE_MMDBSSM or use (expanded from) $(MMDBSSM) works.
# SWIG_FLAGS is either blank or defines COOT_USE_GTK2_INTERFACE usage
#
coot_wrap_guile_pre_gtk1.cc: $(top_srcdir)/src/c-interface.h $(top_srcdir)/src/cc-interface.hh $(top_srcdir)/src/coot.i
	swig -o $@ $(WITH_GUILE_CFLAGS) $(USE_LIBCURL) -I$(top_srcdir)/src -guile -c++ $(top_srcdir)/src/coot.i

coot_wrap_guile_pre_gtk2.cc: $(top_srcdir)/src/c-interface.h $(top_srcdir)/src/cc-interface.hh $(top_srcdir)/src/coot.i
	swig -o $@ -DCOOT_USE_GTK2_INTERFACE $(WITH_GUILE_CFLAGS)  $(USE_LIBCURL) -I$(top_srcdir)/src -guile -c++ $(top_srcdir)/src/coot.i

coot_wrap_guile_gtk2.cc: coot_wrap_guile_pre_gtk2.cc
	$(MAKE_COOT_WRAP_GUILE) gtk2 $^ $@

coot_wrap_guile_gtk1.cc: coot_wrap_guile_pre_gtk1.cc
	$(MAKE_COOT_WRAP_GUILE) gtk1 $^ $@

# python

coot_wrap_python_pre_gtk1.cc: $(top_srcdir)/src/c-interface.h $(top_srcdir)/src/cc-interface.hh $(top_srcdir)/src/coot.i
	swig -o $@ $(PYTHON_CFLAGS)  $(USE_LIBCURL) -I$(top_srcdir)/src -shadow -python -c++ $(top_srcdir)/src/coot.i && mv coot.py coot_gtk1.py

coot_wrap_python_pre_gtk2.cc: $(top_srcdir)/src/c-interface.h $(top_srcdir)/src/cc-interface.hh $(top_srcdir)/src/coot.i
	swig -o $@ -DCOOT_USE_GTK2_INTERFACE $(PYTHON_CFLAGS)  $(USE_LIBCURL) -I$(top_srcdir)/src -shadow -python -c++ $(top_srcdir)/src/coot.i && mv coot.py coot_gtk2.py

coot_wrap_python_gtk1.cc: coot_wrap_python_pre_gtk1.cc
	$(MAKE_COOT_WRAP_PYTHON_XXX) gtk1 $^ $@

coot_wrap_python_gtk2.cc: coot_wrap_python_pre_gtk2.cc
	$(MAKE_COOT_WRAP_PYTHON_XXX) gtk2 $^ $@

coot_gtk1.py: coot_wrap_python_pre_gtk1.cc

coot_gtk2.py: coot_wrap_python_pre_gtk2.cc

# coot.py is different for the GTK1 and GTK2 versions.
#
if USING_GTK2_VERSION
coot.py: coot_gtk2.py
	cp coot_gtk2.py coot.py
else
coot.py: coot_gtk1.py
	cp coot_gtk1.py coot.py
endif



#libcoot_wrap.so: coot_wrap.o
#	-DVERSION=$(VERSION) already in config.h
#	g++ --shared -o libcoot_wrap.so coot_wrap.o

AM_CFLAGS = \
	-DPKGDATADIR='"$(pkgdatadir)"' \
	$(WITH_GUILE_CFLAGS)	\
	$(USE_GUILE_GTK) 	\
	$(GUILE_CFLAGS)		\
	$(PYTHON_CFLAGS)

# must have GUILE and GUILE_GTK flags before $(INCLUDES), and this 
# does not happen because LTCXXCOMPILE has AM_CXXFLAGS after INCLUDES.
# But DEFAULT_INCLUDES is *before* INCLUDES
# DEFAULT_INCLUDES = ... snip.  Rethink...

AM_CXXFLAGS = \
	-DPKGDATADIR='"$(pkgdatadir)"' \
	-DPACKAGE_LOCALE_DIR='"$(pkgdatadir)/locale"' \
	-DDATADIR='"$(datadir)"'    \
        -DCOOT_SYS_BUILD_TYPE='"$(COOT_SYS_BUILD_TYPE)"' \
	$(USE_LIBCURL) 	        \
	$(WITH_GUILE_CFLAGS)	\
	$(USE_GUILE_GTK) 	\
	$(GUILE_CFLAGS)		\
	$(PYTHON_CFLAGS)	\
	$(COOT_WINDOWS_CFLAGS) 	\
	$(COOT_GSL_FLAGS) $(OPENMP_CFLAGS)




# SWIG_FLAGS is something of a misnomer.  It actually defines
# COOT_USE_GTK2_INTERFACE, which is needed to correctly declare the
# correct callbacks in both GTK1 and GTK2 regimes.
#
INCLUDES = \
	-I$(top_srcdir)/compat 		\
	-I$(top_srcdir)/coords 		\
	-I$(top_srcdir)/ccp4mg-utils	\
	-I$(top_srcdir)/utils	        \
	-I$(top_srcdir)/coot-utils	\
	-I$(top_srcdir)/intl 		\
	-I$(top_srcdir)/skeleton	\
	-I$(top_srcdir)/angles		\
	-I$(top_srcdir)/db-main		\
	-I$(top_srcdir)/mini-mol	\
	-I$(top_srcdir)/ligand		\
	-I$(top_srcdir)/geometry	\
	-I$(top_srcdir)/ideal   	\
	-I$(top_srcdir)/coot-surface	\
	-I$(top_srcdir)/surface		\
	-I$(top_srcdir)/density-contour	\
	-I$(top_srcdir)/build		\
	-I$(top_srcdir)/high-res	\
	-I$(top_srcdir)/analysis	\
	-I$(top_srcdir)/cootaneer	\
	-I$(top_srcdir)/lidia-core      \
	$(WITH_GUILE_CFLAGS)		\
	$(USE_GUILE_GTK) 		\
	$(GUILE_GTK_CFLAGS)		\
	$(GUILE_CFLAGS)			\
	$(CANVAS_CFLAGS) 		\
	$(SWIG_FLAGS)			\
	$(GTKGL_CONTEXT_CFLAGS)		\
	$(LIBGLADE_CFLAGS)		\
	@GTK_CFLAGS@ 			\
	$(GLUT_CFLAGS)			\
	$(MMDBSSM_CXXFLAGS) 		\
	$(MMDB_CXXFLAGS) 		\
	$(GSL_CFLAGS)			\
	$(CLIPPER_CXXFLAGS) 		\
	$(PYGTK_CFLAGS)                 \
	$(CURL_CFLAGS)                 

#	$(WITH_SCRIPT_CFLAGS) 


bin_PROGRAMS = coot-real density-score-by-residue-real fix-nomenclature-errors
dist_bin_SCRIPTS = coot density-score-by-residue
noinst_PROGRAMS = test-coot-align test-cholesky test-lbg-graph
#check_PROGRAMS = test-coot-align test-choleksy
# check_PROGRAMS = test-lbg-graph
# testmap write_map testcartesian graphics_crash_test test_cholesky test-command-line

#coot: coot.in
#	cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@

fix_nomenclature_errors_SOURCES = \
	fix-nomenclature.cc    \
	coot-nomenclature.hh   \
	coot-nomenclature.cc

coot_real_SOURCES = \
	main.cc 		 \
	svn-revision.cc 	 \
	c-inner-main.c		 \
	c-inner-main.h		 \
	c-interface.h            \
	c-interface-build.cc	 \
	c-interface-build-gui.cc \
	c-interface-database.hh  \
	c-interface-database.cc  \
	c-interface-pisa.cc	 \
	pisa-interface.hh	 \
	c-interface-pisa-gui.cc	 \
	c-interface-python.hh	 \
	c-interface.cc		 \
	c-interface-gui.cc	 \
	c-interface-gtk2.cc	 \
	c-interface-info.cc	 \
	c-interface-kk.cc	 \
	c-interface-ligands.hh   \
	c-interface-ligands.cc   \
	c-interface-ligands-gui.cc \
	c-interface-mmdb.cc	 \
	c-interface-mmdb.hh      \
	c-interface-mutate.cc    \
	c-interface-ncs.cc	 \
	c-interface-network.cc	 \
	c-interface-nucleotides.cc \
	c-interface-python.cc	 \
	c-interface-preferences.cc \
	c-interface-refine.cc	 \
	c-interface-scm.hh	 \
	c-interface-scm.cc	 \
	c-interface-superpose.cc \
	c-interface-validate.cc	 \
	c-interface-widgets.cc   \
	c-interface-waters.cc	 \
	c-interface-waters-gui.cc\
	cc-interface.hh		 \
	callbacks.c callbacks.h	 \
        canvas-fixes.hh          \
	cmtz-interface.hh	 \
	cmtz-interface.cc	 \
	command-line.cc		 \
	command-line.hh		 \
	coot-database.hh	 \
	coot-align.cc 		 \
	coot-align.hh		 \
	coot-nomenclature.cc	 \
	coot-fileselections.h    \
	coot-glue.hh             \
	coot-glue.c              \
	coot-glue-gtk2.c	 \
	coot-glue-gtk2-python.cc \
	coot-render.hh           \
	coot_pythonmodule.cc     \
	coot-preferences.h       \
	coot-references.h        \
	coot_wrap_guile_gtk1.cc  \
	coot_wrap_guile_gtk2.cc  \
	generic-display-object.hh\
	geometry-graphs.cc	 \
	geometry-graphs.hh	 \
	get-residue.hh           \
	get-residue.cc           \
	flev.cc                  \
	gl-bits.hh		 \
	globjects.cc 		 \
	globjects.h gl-matrix.h  \
	graphics-info.cc	 \
	graphics-info-defines.cc \
	graphics-info-refmac.cc  \
	graphics-info-graphs.cc	 \
	graphics-info-gui.cc	 \
	graphics-info-modelling.cc \
	graphics-info-navigation.cc	\
	graphics-info-navigation-gtk1.cc\
	graphics-info-navigation-gtk2.cc\
	graphics-info-preferences.cc	\
	graphics-info-network.cc \
	graphics-info-mutate.cc	 \
	graphics-info-lsq.cc	 \
	graphics-info-pick.cc	 \
	graphics-info-render.cc	 \
	graphics-info-state.cc	 \
	graphics-info-superpose.cc	\
	gl-matrix.cc 		 \
	graphics-info.h 	 \
	gtk-manual.cc		 \
	gtk-manual.h		 \
	gtk2-interface.c	 \
	gtk-manual.hh		 \
	guile-fixups.h		 \
	history_list.hh		 \
	interface.c	 	 \
	interface.h              \
	lbg-graph.hh   		 \
	lbg-graph.cc   		 \
        manipulation-modes.hh    \
	molecule-class-info.cc 	 \
	molecule-class-info-maps.cc 	\
	molecule-class-info-other.cc 	\
	molecule-class-info.h 	        \
	molecule-class-info-mutate.cc 	\
	molecule-class-info-residues.cc \
	molecule-class-info-refine.cc   \
	molecule-class-info-surface.cc 	\
	molecule-class-info-ncs.cc 	\
	molecule-class-info-widget-work.cc \
	nsv.cc nsv.hh    	 \
        pick.h pick.cc		 \
	positioned-widgets.h     \
	read-phs.h		 \
	read-phs.c		 \
	read-cif.h		 \
	read-cif.c		 \
	rama_mousey.cc		 \
	rama_mousey.hh		 \
	rama_plot.cc		 \
	rama_plot.hh		 \
	restraints-editor.cc     \
	restraints-editor.hh     \
	restraints-editor-c.h    \
	rotate-translate-modes.hh \
	rotamer-search-modes.hh  \
	support.c support.h 	 \
	support-gtk2-extras.c    \
	testing.cc testing.hh    \
	trackball.c trackball.h  \
	triangle-extras.h 	 \
	select-atom-info.hh      \
	sequence-view.hh         \
	sequence-view.cc         \
	seq-view-interface.c     \
	seq-view-interface.h     \
	user-mods.hh user-mods.cc\
	validation-graphs.hh	 \
	xmap-utils.cc 		 \
	xmap-utils.h		 

if COOT_USE_PYTHON
if USING_GTK2_VERSION
   coot_real_SOURCES += coot_wrap_python_gtk2.cc
else
   coot_real_SOURCES += coot_wrap_python_gtk1.cc
endif
endif


# testcartesian_SOURCES = testcartesian.cc

# graphics_crash_test_SOURCES = graphics-crash-test.cc graphics-crash-test.h


# testcartesian_LDADD = $(COOT_EXTRA_LIBS)

# graphics_crash_test_LDADD = @GTK_LIBS@ $(GTKGL_CONTEXT_LIBS) $(INTLLIBS) $(GLUT_LIBS) -lGL -lGLU

fix_nomenclature_errors_LDADD = \
	$(top_builddir)/ligand/libcoot-ligand.la      		\
	$(top_builddir)/db-main/libcoot-db-main.la      	\
	$(top_builddir)/ideal/libcoot-ideal.la    		\
	$(top_builddir)/geometry/libcoot-geometry.la    	\
	$(top_builddir)/coords/libcoot-coords.la      		\
	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la 	\
	$(top_builddir)/coot-utils/libcoot-coord-extras.la	\
	$(top_builddir)/coot-utils/libcoot-coord-utils.la	\
	$(top_builddir)/coot-utils/libcoot-map-utils.la		\
	$(top_builddir)/coot-utils/libcoot-shelx.la	  	\
	$(top_builddir)/utils/libcoot-utils.la	  	        \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la 		\
	$(GSL_LIBS)		\
	$(CLIPPER_LIBS) 	\
	$(MMDB_LIBS)		


# coot_real_LDFLAGS = -static -lstdc++
coot_real_LDFLAGS =

coot_real_LDADD =  $(OPENMP_CFLAGS)     \
                $(COOT_EXTRA_LIBS)      \
		$(GUILE_GTK_LIBS)       \
		$(GUILE_LDFLAGS)	\
		$(CLIPPER_LIBS) 	\
		$(MMDBSSM_LIBS)		\
		$(MMDB_LIBS)		\
		$(PYTHON_LIBS)		\
		$(GTKGL_CONTEXT_LIBS) 	\
		$(GLUT_LIBS)		\
		$(CANVAS_LIBS) 		\
		$(GSL_LIBS)		\
		$(MYSQL_LIBS)		\
		$(CURL_LIBS)		\
		$(LIBGLADE_LIBS)	\
		@GTK_LIBS@	

#                -lglade-2.0


#		@COOT_WRAP_LIBS@	do it a different way

#		$(INTLLIBS) 		Fix when GETTEXT workds

test_cholesky_SOURCES = test-chol.cc gl-matrix.cc
test_cholesky_LDADD =	\
	$(top_builddir)/coords/libcoot-coords.la          \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la      \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/coot-utils/libcoot-shelx.la       \
	$(top_builddir)/geometry/libcoot-geometry.la      \
	$(top_builddir)/utils/libcoot-utils.la            \
	$(CLIPPER_LIBS) $(MMDB_LIBS) $(GSL_LIBS)

# testmap_SOURCES = testmap.cc xmap-interface.cc
# testmap_LDADD = $(COOT_EXTRA_LIBS) $(SURFACE_LIB) $(CLIPPER_LIBS) 


#write_map_SOURCES = write-map.cc xmap-interface.cc
#write_map_LDADD = $(COOT_EXTRA_LIBS) $(SURFACE_LIB)  $(CLIPPER_LIBS)
# example_LDADD =$(GNOME_LIBDIR)$(GNOMEUI_LIBS) $(INTLLIBS)


# test_command_line_SOURCES = test-command-line.cc
# test_command_line_LDADD = -L$(top_srcdir)/compat      -lcoot-compat    

# if you like, you can investigate why coot-geometry and clipper is needed here:
test_coot_align_SOURCES = align.cc coot-align.cc
test_coot_align_LDADD =	\
	$(top_builddir)/geometry/libcoot-geometry.la    	\
	$(top_builddir)/coords/libcoot-coords.la 		\
	$(top_builddir)/mini-mol/libcoot-mini-mol.la 		\
	$(top_builddir)/coot-utils/libcoot-coord-utils.la 	\
	$(top_builddir)/coot-utils/libcoot-shelx.la  		\
	$(top_builddir)/utils/libcoot-utils.la  		\
	$(CLIPPER_LIBS) 	\
	$(MMDB_LIBS)

test_lbg_graph_SOURCES = test-lbg-graph.cc lbg-graph.hh lbg-graph.cc
test_lbg_graph_LDADD = $(MMDB_LIBS)


density_score_by_residue_real_SOURCES = density-score-by-residue-src.cc
density_score_by_residue_real_LDADD  =	\
	$(top_builddir)/geometry/libcoot-geometry.la    	\
	$(top_builddir)/coords/libcoot-coords.la 		\
	$(top_builddir)/mini-mol/libcoot-mini-mol.la 		\
	$(top_builddir)/coot-utils/libcoot-map-utils.la  	\
	$(top_builddir)/coot-utils/libcoot-coord-utils.la  	\
	$(top_builddir)/coot-utils/libcoot-shelx.la  		\
	$(top_builddir)/utils/libcoot-utils.la  		\
	$(top_builddir)/compat/libcoot-compat.la  		\
	$(CLIPPER_LIBS) 	\
	$(MMDB_LIBS)		\
	$(GSL_LIBS)
