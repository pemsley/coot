# src/Makefile.am
#
# Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
# Copyright 2007, 2008, 2009, 2010, 2011, 2012 by The University of Oxford
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

# if coot.py is called __init.py and goes in site-packages/coot/
#dist_pkgpython_SCRIPTS = __init__.py
#dist_python_SCRIPTS = coot.pth

#__init__.py: coot.py
#	cp coot.py __init__.py


dist_python_SCRIPTS = coot.py coot.pth



######## 20091001 ########################################################


# COOT_STDTYPES is something like: DHAVE_SYS_STDTYPES_H=1

# USE_GUILE_DIRECTIVE=`if [ -z "$(COOT_USE_GUILE)"]; then echo ""; else echo -DUSE_GUILE; fi`

WITH_GUILE_CFLAGS=$(COOT_STDTYPES) $(COOT_USE_GUILE)


# GUILE_VERSION=$(GUILE_CONFIG)

MAKE_COOT_WRAP_GUILE=$(SHELL) $(top_srcdir)/src/post-process-coot-wrap-guile.sh "$(GUILE_CONFIG)"

MAKE_COOT_WRAP_PYTHON_XXX=$(SHELL) $(top_srcdir)/src/post-process-coot-wrap-python.sh


##########################################################################


# put these into configure
GLM_CXXFLAGS = -I$(HOME)/glm/include
EPOXY_LIBS = -lepoxy

CANVAS_CFLAGS = $(GTKCANVAS_CFLAGS) $(GNOME_CANVAS_CFLAGS)
CANVAS_LIBS   = $(GTKCANVAS_LIBS) $(GNOME_CANVAS_LIBS)

# and similarly for the GL context:
GTKGL_CONTEXT_CFLAGS = $(GTKGL_CFLAGS) $(GTKGLEXT_CFLAGS)
GTKGL_CONTEXT_LIBS   = $(GTKGL_LIBS) $(GTKGLEXT_LIBS)

# MOLECULESTOTRIANGLES=$(install_top_dir)
# MOLECULESTOTRIANGLES_LIBS = -L$(MOLECULESTOTRIANGLES)/lib -lMoleculesToTrianglesCXXClasses
# MOLECULESTOTRIANGLES_CXXFLAGS = -I$(MOLECULESTOTRIANGLES)/include

# remove this from coot dependencies:
# 	$(top_builddir)/mmtz-extras/libcoot-mmtz-extras.la
#

COOT_EXTRA_LIBS = \
        $(top_builddir)/lbg/libcoot-lidia.la      		\
	$(top_builddir)/pli/libcoot-pli.la      		\
	$(top_builddir)/ligand/libcoot-res-tracer.la  		\
	$(top_builddir)/ligand/libcoot-ligand.la      		\
	$(top_builddir)/lidia-core/libcoot-lidia-core.la      	\
	$(top_builddir)/goograph/libcoot-goograph.la      	\
	$(top_builddir)/density-contour/libcoot-isosurface.la  	\
	$(top_builddir)/build/libcoot-cabuild.la		\
	$(top_builddir)/skeleton/libcoot-skeleton.la    	\
	$(top_builddir)/angles/libcoot-AngleInfo.la   		\
	$(top_builddir)/db-main/libcoot-db-main.la     		\
	$(top_builddir)/high-res/libcoot-high-res.la  		\
	$(top_builddir)/ideal/libcoot-ideal.la    		\
	$(top_builddir)/geometry/libcoot-geometry.la    	\
	$(top_builddir)/coords/libcoot-coords.la      		\
	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la 	\
	$(top_builddir)/coot-utils/libcoot-map-heavy.la  	\
	$(top_builddir)/coot-utils/libcoot-map-utils.la  	\
	$(top_builddir)/coot-utils/libcoot-coord-extras.la  	\
	$(top_builddir)/coot-utils/libcoot-coord-utils.la  	\
	$(top_builddir)/utils/libcoot-utils.la  		\
	$(top_builddir)/mini-mol/libcoot-mini-mol.la    	\
	$(top_builddir)/analysis/libcoot-analysis.la 		\
	$(top_builddir)/protein_db/libcoot-protein-db.la      	\
	$(top_builddir)/cootilus/libcootilus.la 		\
	$(top_builddir)/cootaneer/libcootaneer.la 		\
	$(top_builddir)/MoleculesToTriangles/CXXClasses/libMoleculesToTrianglesCXXClasses.la     		\
	$(top_builddir)/MoleculesToTriangles/CXXSurface/libMoleculesToTrianglesCXXSurface.la     		\
	$(top_builddir)/compat/libcoot-compat.la



# Remember to regenerate coot_wrap_python_pre.cc before making dist, eh?
#
EXTRA_DIST = \
	coot.i post-glade add-c-extern.sh blank.cc fix-gettext.sh \
	mol-hierarchy.cc \
	junk.cc coot-gtk2.defs coot-dox.cfg post-gtk2-glade.sh \
        post-process-coot-wrap-guile.sh \
        post-process-coot-wrap-python.sh \
        fixup-interface.h.sh \
	coot_wrap_guile_pre.cc dynarama.glade


# subpkgdatadir = $(pkgdatadir)
# subpkgdata_DATA = dynarama.glade geometry_graphs_builder.glade

cootgladedir = $(pkgdatadir)/glade
dist_cootglade_DATA = dynarama.glade geometry_graphs_builder.glade

cootdatadir = $(pkgdatadir)/shaders
dist_cootdata_DATA = atom-label.shader blur.shader central-cube.shader \
                     hud-text.shader  hud-bars.shader hud-labels.shader \
                     map.shader model.shader screen.shader \
		     moleculestotriangles.shader \
                     instanced-objects.shader particles.shader lines.shader \
                     draw-map-cap.shader lines-pulse.shader rama-balls.shader \
		     hud-geometry-tooltip-text.shader \
		     residue-markers.shader extra-distance-restraints.shader \
		     residue-markers-for-ssao.shader \
		     model-as-mesh.shader ligand-view.shader symmetry-atoms-lines.shader \
		     rama-plot-phi-psi-markers.shader rama-plot-axes-and-ticks.shader \
		     hud-image-texture.shader hud-lines.shader \
		     depth-of-field.shader blur-x.shader blur-y.shader texture-meshes.shader \
		     camera-facing-quad-shader-for-testing.shader outline-of-active-residue.shader \
                     meshes.shader meshes-shadow-map.shader meshes-with-shadows.shader \
                     meshes-for-ssao.shader effects.shader 9.ssao_blur.shader 9.ssao_geometry.shader 9.ssao.shader \
		     texture-meshes-for-ssao.shader texture-meshes-with-shadows.shader \
		     texture-meshes-shadow-map.shader shadow-map-image-texture.shader \
		     rotation-centre-cross-hairs-for-ssao.shader background-image.shader



# use -DHAVE_MMDBSSM or use (expanded from) $(MMDBSSM) works.
# SWIG_FLAGS is either blank or defines COOT_USE_GTK2_INTERFACE usage
#
coot_wrap_guile_pre.cc:      $(top_srcdir)/src/c-interface.h    \
                             $(top_srcdir)/src/cc-interface.hh  \
                             $(top_srcdir)/src/cc-interface-ncs.hh  \
                             $(top_srcdir)/src/cc-interface-molecular-representation.hh  \
                             $(top_srcdir)/src/cc-interface-scripting.hh \
                             $(top_srcdir)/src/coot.i           \
                             $(top_srcdir)/src/lbg-interface.hh \
                             $(top_srcdir)/src/c-interface-refine.hh \
                             $(top_srcdir)/src/c-interface-ligands-swig.hh \
                             $(top_srcdir)/src/c-interface-generic-objects.h \
                             $(top_srcdir)/src/sdf-interface.hh \
                             $(top_srcdir)/src/cfc.hh \
                             $(top_srcdir)/src/curlew.h \
                             $(top_srcdir)/src/curlew.hh \
                             $(top_srcdir)/src/validation.hh \
                             $(top_srcdir)/src/cmtz-interface.hh \
                             $(top_srcdir)/src/c-interface-mogul.hh
	swig -o $@ $(CPPFLAGS) $(MOLTRIANGLES_CPPFLAGS) -DCOOT_USE_GTK2_INTERFACE $(ENHANCED_LIGAND_TOOLS_CPPFLAGS) $(WITH_GUILE_CFLAGS)  $(USE_LIBCURL) -I$(top_srcdir)/src -I$(top_srcdir) -guile -c++  $(top_srcdir)/src/coot.i

# For GOOPS support: -proxy -Linkage module (or passive), however, doing that makes SWIG_init() go away.
# -proxy  exports the GOOPS wrapper definitions into module.scm

coot_wrap_guile.cc: coot_wrap_guile_pre.cc
	$(MAKE_COOT_WRAP_GUILE) gtk2 $^ $@

# python

# you can't use PYTHON_CFLAGS for swig because that introduces
# -iwithsysroot dir - and swig barfs on that.
# so make a variable that gives us python flags without -iwithsysroot
coot_wrap_python_pre.cc: $(top_srcdir)/src/c-interface.h      \
                         $(top_srcdir)/src/cc-interface.hh    \
                         $(top_srcdir)/src/coot.i             \
                         $(top_srcdir)/src/cc-interface-ncs.hh  \
                         $(top_srcdir)/src/cc-interface-molecular-representation.hh  \
	                      $(top_srcdir)/src/cc-interface-scripting.hh \
                         $(top_srcdir)/src/lbg-interface.hh   \
                         $(top_srcdir)/src/c-interface-refine.hh \
                         $(top_srcdir)/src/c-interface-ligands-swig.hh \
                         $(top_srcdir)/src/sdf-interface.hh   \
                         $(top_srcdir)/src/cfc.hh \
                         $(top_srcdir)/src/curlew.h \
                         $(top_srcdir)/src/curlew.hh \
                         $(top_srcdir)/src/validation.hh \
                         $(top_srcdir)/src/c-interface-mogul.hh
	swig -o $@ $(CPPFLAGS) $(MOLTRIANGLES_CPPFLAGS) -DCOOT_USE_GTK2_INTERFACE $(ENHANCED_LIGAND_TOOLS_CPPFLAGS) $(PYTHON_CPPFLAGS)  -DUSE_PYTHON $(USE_LIBCURL) -I$(top_srcdir)/src  -I$(top_srcdir) -shadow -python -py3 -c++ $(top_srcdir)/src/coot.i

coot_wrap_python.cc: coot_wrap_python_pre.cc
	$(MAKE_COOT_WRAP_PYTHON_XXX) gtk2 $^ $@

# we need to run swig to make coot.py
coot.py: coot_wrap_python.cc

#libcoot_wrap.so: coot_wrap.o
#	-DVERSION=$(VERSION) already in config.h
#	g++ --shared -o libcoot_wrap.so coot_wrap.o

AM_CFLAGS = \
	-DPKGDATADIR='"$(pkgdatadir)"' \
	$(WITH_GUILE_CFLAGS)	\
	$(USE_GUILE_GTK) 	\
	$(GUILE_CFLAGS)		\
	$(PYTHON_CPPFLAGS)	\
	$(PYTHON_CFLAGS)

# must have GUILE and GUILE_GTK flags before $(INCLUDES), and this
# does not happen because LTCXXCOMPILE has AM_CXXFLAGS after INCLUDES.
# But DEFAULT_INCLUDES is *before* INCLUDES
# DEFAULT_INCLUDES = ... snip.  Rethink...

# get this from configure
# /Users/pemsley/python-3.8/lib/python3.8/site-packages/PyGObject-3.36.1-py3.8-macosx-11.5-x86_64.egg/include/pygobject-3.0/pygobject.h
# /Users/pemsley/test/ccp4-20210929/include/pygobject-3.0/pygobject.h
PYGOBJECT_CFLAGS=-I$(HOME)/python-3.8/lib/python3.8/site-packages/PyGObject-3.36.1-py3.8-macosx-11.5-x86_64.egg/include



AM_CXXFLAGS = \
        -DUSE_MOLECULES_TO_TRIANGLES \
        $(MOLECULESTOTRIANGLES_CXXFLAGS) \
	-DPKGDATADIR='"$(pkgdatadir)"' \
	-DPACKAGE_LOCALE_DIR='"$(pkgdatadir)/locale"' \
	-DPKGPYTHONDIR='"$(pkgpythondir)"' \
	-DPYTHONDIR='"$(pythondir)"' \
	-DDATADIR='"$(datadir)"'    \
        -DCOOT_SYS_BUILD_TYPE='"$(COOT_SYS_BUILD_TYPE)"' \
        -DCOOT_BUILD_INFO_STRING='"$(COOT_BUILD_INFO_STRING)"' \
	$(USE_LIBCURL) 	        \
	$(WITH_GUILE_CFLAGS)	\
	$(USE_GUILE_GTK) 	\
	$(GUILE_CFLAGS)		\
        $(PYGOBJECT_CFLAGS)     \
	$(PYTHON_CPPFLAGS)	\
	$(PYTHON_CFLAGS)	\
	$(COOT_WINDOWS_CFLAGS) 	\
	$(COOT_GSL_FLAGS) $(OPENMP_CFLAGS)


# SWIG_FLAGS is something of a misnomer.  It actually defines
# COOT_USE_GTK2_INTERFACE, which is needed to correctly declare the
# correct callbacks in both GTK1 and GTK2 regimes.
#
#
AM_CPPFLAGS = \
        -DUSE_PYTHON \
	@GTK_CFLAGS@ 			\
	-I$(top_srcdir)                 \
	$(RDKIT_CXXFLAGS) 		\
	$(BOOST_CPPFLAGS) 		\
	$(ENHANCED_LIGAND_TOOLS_CPPFLAGS) \
	$(GLM_CXXFLAGS)                 \
	$(SWIG_FLAGS)			\
	$(GTKGL_CONTEXT_CFLAGS)		\
	$(GOOCANVAS_CFLAGS)             \
	$(WITH_GUILE_CFLAGS)		\
	$(USE_GUILE_GTK) 		\
	$(GUILE_GTK_CFLAGS)		\
	$(GUILE_CFLAGS)			\
	$(CANVAS_CFLAGS) 		\
	$(MOLTRIANGLES_CPPFLAGS)	\
	$(LIBSSM_CXXFLAGS) 		\
	$(MMDB_CXXFLAGS) 		\
	$(CCP4SRS_CXXFLAGS) 		\
	$(GSL_CFLAGS)			\
	$(CLIPPER_CXXFLAGS) 		\
	$(PYGTK_CFLAGS)                 \
	$(GLOB_CFLAGS)                  \
	$(WITH_SOUND)                   \
	$(VORBISFILE_CFLAGS)            \
	$(SQLITE3_CFLAGS)               \
	$(GEMMI_CXXFLAGS)               \
	$(CURL_CFLAGS)

#	$(WITH_SCRIPT_CFLAGS)

pkginclude_HEADERS = crunch-model.hh

_coot_embedded_la_SOURCES = coot-embedded.cc
_coot_extended_la_SOURCES = coot-extended.cc

# use -Wl,--no-undefined or  -Wl,-undefined,error (clang) ?
_coot_embedded_la_LDFLAGS = -module $(SHARED_LDFLAGS) --no-undefined
_coot_embedded_la_LIBADD = \
                       $(top_builddir)/lidia-core/libcoot-lidia-core.la   \
                       $(RDKIT_LIBS) $(BOOST_LDFLAGS) $(L_BOOST_PYTHON)   \
                       $(PYTHON_LIBS) \
                       $(MMDB_LIBS)

_coot_extended_la_LDFLAGS = -module $(SHARED_LDFLAGS) --no-undefined
_coot_extended_la_LIBADD = \
                        $(top_builddir)/pli/libcoot-pli.la                 \
                        $(top_builddir)/coot-utils/libcoot-coord-utils.la  \
                        $(top_builddir)/lidia-core/libcoot-lidia-core.la   \
                        $(top_builddir)/geometry/libcoot-geometry.la       \
                        $(top_builddir)/utils/libcoot-utils.la             \
                        $(RDKIT_LIBS) $(BOOST_LDFLAGS) $(L_BOOST_PYTHON)   \
                        $(PYTHON_LIBS) \
                        $(MMDB_LIBS)

pyexec_LTLIBRARIES = _coot_embedded.la _coot_extended.la _coot.la

_coot_la_SOURCES = coot.cc coot_wrap_python.cc

_coot_la_LDFLAGS = -module $(SHARED_LDFLAGS) --no-undefined  \
           libcootsumo.la                \
           $(top_builddir)/coot-utils/libcoot-coord-utils.la \
           $(top_builddir)/coords/libcoot-coords.la          \
	   $(top_builddir)/geometry/libcoot-geometry.la      \
	   $(top_builddir)/pli/libcoot-pli.la                \
	   $(top_builddir)/lidia-core/libcoot-lidia-core.la  \
           $(top_builddir)/utils/libcoot-utils.la            \
           $(PYTHON_LIBS) $(CLIPPER_LIBS)


if OS_WIN32
# no more libexec for windows
bin_PROGRAMS     = coot-bin coot-density-score-by-residue-bin dynarama-bin
else
libexec_PROGRAMS = coot-bin coot-density-score-by-residue-bin dynarama-bin
bin_PROGRAMS = 
endif

dist_bin_SCRIPTS = coot coot-density-score-by-residue dynarama
check_PROGRAMS = test-coot-align test-cholesky

#coot_fix_nomenclature_errors_SOURCES = \
#	fix-nomenclature.cc    \
#	coot-nomenclature.cc

coot_bin_SOURCES = main.cc glade-callbacks.cc glade-callbacks-main-window.cc glade-callbacks-preferences.cc glade-callbacks-shader-settings.cc

coot_bin_LDADD =  \
        libcootsumo.la \
	$(top_builddir)/utils/libcoot-utils.la \
        $(top_builddir)/lidia-core/libcoot-lidia-core.la  \
	$(CLIPPER_LIBS) $(CURL_LIBS) $(GTK_LIBS)


lib_LTLIBRARIES=libcootsumo.la

libcootsumo_la_SOURCES = \
	atom-pull.hh		 \
        atom-label-info.hh       \
	array-2d.hh		 \
	init_coot_as_python_module.cc \
	init_coot_as_python_module.hh \
	boids.hh boids.cc	 \
        bond-colour-mode.hh      \
	boot-python.hh           \
	boot-python.cc           \
	atom-overlaps.cc         \
	git-revision-count.cc 	 \
	scm-boot-guile.hh 	 \
	scm-boot-guile.cc 	 \
	clipper-ccp4-map-file-wrapper.hh \
	c-interface.h            \
	c-interface-analysis.cc	 \
	c-interface-bonds.hh	 \
	c-interface-bonds.cc	 \
	c-interface-build.cc	 \
	c-interface-build-morph.cc	 \
	c-interface-build-main-chain.cc	 \
	c-interface-build-side-chain.cc	 \
	c-interface-build-symmetry.cc	 \
	c-interface-build-gui.cc \
	c-interface-curlew.cc    \
	c-interface-curlew.hh    \
	c-interface-database.hh  \
	c-interface-database.cc  \
	c-interface-generic-objects.h	 \
	c-interface-geometry-distortion.hh \
	c-interface-gui.cc	 \
	c-interface-gui.hh	 \
	c-interface-gtk2.cc	 \
	c-interface-gtk-widgets.h \
	c-interface-info.cc	 \
	c-interface-kk.cc	 \
	c-interface-image.cc     \
	c-interface-image-widget.hh    \
	c-interface-ligand-search.hh   \
	c-interface-ligands.hh         \
	c-interface-ligands-widgets.hh \
	c-interface-ligands-swig.hh    \
	c-interface-ligands-gui.cc \
	c-interface-ligands.cc   \
	c-interface-maps.cc	 \
	c-interface-mmdb.cc	 \
	c-interface-mmdb.hh      \
	c-interface-mogul.cc     \
	c-interface-mogul.hh     \
	c-interface-mogul-scripting.cc     \
	c-interface-mutate.cc    \
	c-interface-ncs.cc	 \
	c-interface-network.cc	 \
	c-interface-nucleotides.cc \
	c-interface-pisa.cc	 \
	c-interface-pisa-gui.cc	 \
	c-interface-preferences.cc \
	c-interface-preferences.h  \
	c-interface-python.cc	 \
	c-interface-python.hh	 \
	c-interface-refine.hh	 \
	c-interface-refine.h	 \
	c-interface-refine.cc	 \
	c-interface-residues.cc	 \
	c-interface-refmac.cc	 \
	simple-refmac.cc	 \
	c-interface-refmac.h	 \
	c-interface-scm.hh	 \
	c-interface-scm.cc	 \
	c-interface-sequence.cc  \
	c-interface-sequence.hh  \
	c-interface-sfs.cc	 \
	c-interface-sharpen.cc   \
	c-interface-superpose.cc \
	c-interface-validate.cc	 \
	c-interface-view.cc	 \
	c-interface-widgets.cc   \
	c-interface-widgets.hh   \
	c-interface-waters.cc	 \
	c-interface-waters-gui.cc\
	c-interface.cc		 \
	c-interface-test.cc      \
	cc-interface.hh		 \
	cc-interface-alignment.cc \
	cc-interface-graphics.cc \
	cc-interface-merge.cc    \
	cc-interface-alignment.hh \
	cc-interface-ncs.hh      \
	cc-interface-molecular-representation.cc \
	cc-interface-molecular-symmetry.cc \
	cc-interface-network.hh  \
	cc-interface-mogul.hh    \
	cc-interface-scripting.cc\
	cc-interface-scripting.hh\
	cc-interface-map-utils.cc\
	cc-interface-image.hh    \
	cc-interface-png.cc      \
	cc-interface-user-defined-atom-colours.hh \
	cc-interface-user-defined-atom-colours.cc \
	cc-interface-molecular-representation.hh  \
        canvas-fixes.hh          \
        cfc.cc cfc.hh            \
        cfc-widgets.hh           \
        cfc-widgets-c-interface.h \
	change-dir.cc		 \
	change-dir.hh		 \
	cmtz-interface-gui.hh	 \
	cmtz-interface.hh	 \
	cmtz-interface.cc	 \
	clipper-to-glm.hh        \
	command-arg.hh           \
	command-line.cc		 \
	command-line-extern.cc   \
	command-line-extern.hh   \
	command-line.hh		 \
	coot-database.hh	 \
	coot-colour.hh	 	 \
	coot-colour.cc	 	 \
	coot-align.cc 		 \
	coot-align.hh		 \
	coot-nomenclature.cc	 \
	coot-nomenclature.hh     \
	coot-fileselections.h    \
	coot-glue.hh             \
	coot-glue-gtk2.c	 \
	coot-glue-gtk2-python.cc \
	coot-init-glue.hh        \
	coot-hydrogens.hh        \
	coot-hydrogens.cc        \
	coot-render.hh           \
	coot_pythonmodule.cc     \
	coot-preferences.h       \
	coot-references.h        \
	coot-setup-python.cc     \
	coot-setup-python.hh     \
	coot-version.hh          \
	coot_wrap_guile.cc       \
        curlew.h                 \
        curlew.hh                \
	cylinder.hh		 \
	cylinder.cc		 \
	cylinder-with-rotation-translation.hh	 \
	cylinder-with-rotation-translation.cc	 \
	data-pair-remover.hh     \
        density-results-container-t.hh \
        dots-representation.hh   \
	draw.cc                  \
	draw-2.cc                \
	draw-2.hh                \
	draw.hh                  \
	drag-and-drop.hh         \
	drag-and-drop.cc         \
	draw-generic-display-objects.cc \
	dynamic-menus.hh         \
	dynamic-menus.cc         \
	eyes.cc eyes.hh		 \
	fit-loop-gui.cc          \
	fit-loop-gui.hh          \
	flev.hh			 \
	flev.cc                  \
        framebuffer.cc           \
        framebuffer.hh           \
        fresnel-settings.hh      \
	ft-character.hh \
	old-generic-display-object.hh \
	old-generic-display-object.cc \
	gaussian-surface.hh      \
	gaussian-surface.cc      \
	generic-objects.cc   	 \
        generic-display-objects-c.h \
	generic-vertex.hh	 \
	gensurf.hh               \
	geometry-graphs.cc	 \
	geometry-graphs.hh	 \
	get-residue.hh           \
	get-residue.cc           \
	get-monomer.hh           \
	get-monomer.cc           \
	gl-rama-plot.cc          \
	gl-rama-plot.hh          \
	gl-bits.hh		 \
	globjects.cc 		 \
	globjects.h gl-matrix.h  \
	globularize.cc           \
	glarea_tick_function.hh  \
	globularize.hh           \
	graphics-info.cc	 \
	graphics-info-defines.cc \
	graphics-info-draw.cc    \
	graphics-info-draw-models.cc \
	graphics-info-graphs.cc	 \
	graphics-info-gui.cc	 \
	graphics-info-init.cc	 \
	graphics-info-intermediate.cc    \
	graphics-info-modelling.cc       \
	graphics-info-navigation.cc	 \
	graphics-info-navigation-gtk2.cc \
	graphics-info-preferences.cc	 \
	graphics-info-network.cc \
	graphics-info-mutate.cc	 \
	graphics-info-lsq.cc	 \
	graphics-info-opengl.cc	 \
	graphics-info-pick.cc	 \
	graphics-info-refmac.cc  \
	graphics-info-render.cc	 \
	graphics-info-render-scene.cc \
	graphics-info-residues.cc     \
	graphics-info-refine.cc	 \
	graphics-info-scripting.cc	\
	graphics-info-state.cc	 \
	graphics-info-statics.cc \
	graphics-info-superpose.cc	\
	graphics-info-trackball.cc	\
	graphics-ligand-view.cc	 \
	graphics-ligand-view.hh	 \
	graphics-ligand-mesh-molecule.cc \
	graphics-ligand-mesh-molecule.hh \
	gl-matrix.cc 		 \
	graphics-info.h 	 \
	gtk-manual.cc		 \
	gtk-manual.h		 \
	gtk-widget-conversion-utils.h  \
	gtk-widget-conversion-utils.cc \
	gtk-manual.hh		 \
	guile-fixups.h		 \
	history_list.hh		 \
	HUDMesh.cc               \
	HUDMesh.hh               \
	HUDTextureMesh.cc        \
	HUDTextureMesh.hh        \
	Instanced-Markup-Mesh.hh \
	Instanced-Markup-Mesh.cc \
        idles.hh                 \
	interface.h              \
	key-bindings.cc		 \
	key-bindings.hh		 \
        manipulation-modes.hh    \
	LinesMesh.cc LinesMesh.hh \
	LigandViewMesh.cc        \
	LigandViewMesh.hh        \
        lights-info.hh           \
	lsq-dialog-values.hh     \
	lbg-interface.hh         \
	lbg-callbacks-outer.cc   \
	lbg-interface.cc         \
	ligand-check.hh          \
	Model.cc Model.hh        \
	make-a-dodec.cc make-a-dodec.hh \
	matrix-utils.hh          \
	matrix-utils.cc          \
        map-statistics.hh        \
        Material.hh              \
        map_triangle.hh          \
	Mesh.cc Mesh.hh          \
	Mesh-from-glTF.cc        \
	Mesh-from-graphical-bonds-instanced.cc     \
	Mesh-from-graphical-bonds.cc     \
	meshed-generic-display-object.hh \
	meshed-generic-display-object.cc \
	map-statistics.hh        \
	merge-molecule-results-info-t.hh \
	model-composition-statistics.hh  \
        molecular-triangles-mesh.hh      \
	molecule-class-info.cc 	         \
	molecule-class-info-alignment.cc \
	molecule-class-info-analysis.cc \
	molecule-class-info-bonds.cc    \
	molecule-class-info-build.cc 	\
	molecule-class-info-maps.cc 	\
	molecule-class-info-maps-2.cc 	\
	molecule-class-info-molecular-symmetry.cc \
	molecule-class-info-kdc.cc      \
	molecule-class-info-dots.cc 	\
	molecule-class-info-dots-internal.cc 	\
	molecule-class-info-other.cc 	\
	molecule-class-info.h 	        \
	molecule-class-info-mutate.cc 	\
	molecule-class-info-residues.cc \
	molecule-class-info-refine.cc   \
	molecule-class-info-kk.cc       \
	molecule-class-info-surface.cc 	\
	molecule-class-info-string.cc 	\
	molecule-class-info-patterson.cc   \
	molecule-class-info-ncs.cc         \
	molecule-class-info-mol-tris.cc	   \
	molecule-class-info-widget-work.cc \
	molecule-class-info-geometry.cc    \
	molecular-mesh-generator.cc \
	molecular-mesh-generator-mol-tris.cc \
	molecular-mesh-generator.hh \
	mtz-column-auto-read.hh  \
	named-rotamer-score.hh   \
	nsv.cc nsv.hh    	 \
	oct.cc oct.hh            \
	obj_loader.h obj_loader.cc \
	pisa-interface.hh	 \
        pick.hh pick.cc		 \
	positioned-widgets.h     \
	pre-load.hh              \
	prideout-octasphere.hh   \
	prideout-octasphere.cc   \
	probe-clash-score.cc     \
	probe-clash-score.hh     \
	protein_db-interface.hh  \
	protein_db-interface.cc  \
	python-classes.cc	 \
	python-classes.hh	 \
	python-3-interface.hh 	 \
	Particle.cc Particle.hh  \
	pulse-data.hh            \
	read-phs.h		 \
	read-phs.c		 \
	read-cif.h		 \
	read-cif.c		 \
	rama_mousey.cc		 \
	rama_mousey.hh		 \
	rama_plot.cc		 \
	rama_plot.hh		 \
        text-rendering-utils.cc  \
        text-rendering-utils.hh  \
	restraints-editor.cc     \
	restraints-editor.hh     \
	restraints-editor-c.h    \
	rotate-translate-modes.hh \
	rotamer-search-modes.hh  \
	rsr-functions.cc         \
	rsr-functions.hh         \
	screen-utils.cc		 \
	screendump-tga.hh        \
	screendump-tga.cc        \
	sound.cc \
	sdf-interface.hh         \
	sdf-internal.hh          \
	sdf-interface.cc         \
	select-atom-info.hh      \
	sequence-view.hh         \
	sequence-view.cc         \
	seq-view-interface.c     \
	seq-view-interface.h     \
	Shader.cc Shader.hh      \
        simple-distance-object.hh \
        sound.cc                 \
	sound.hh                 \
	startup-scripts.hh       \
	support.c support.h      \
	stb_image.h              \
	stb_image_write.h        \
	startup-scripts.hh       \
	tiny_gltf.h              \
	support-gtk2-extras.c    \
	testing.cc testing-2.cc  \
	testing-data.hh          \
	testing.hh               \
	Texture.cc Texture.hh    \
	TextureMesh.cc TextureMesh.hh \
	trackball.c trackball.h  \
	triangle-extras.h        \
	user-mods.hh user-mods.cc\
	validation-graphs.hh	 \
	xmap-utils.cc 		 \
	xmap-utils.h             \
	remarks-browser-gtk-widgets.hh \
	rama-plot-phi-psi.cc     \
	rama-plot-phi-psi.hh     \
	rama-rota-score.hh 	 \
	spherical-density.cc     \
	single-map-properties-dialog.cc \
	single-map-properties-dialog.hh \
	fragment-info.hh 	 \
	atom-name-bits.hh 	 \
	additional-representation.hh 	   \
	extra-restraints-representation.hh \
	atom-attribute.hh 	 \
	animated-ligand.hh 	 \
	atom-selection.hh 	 \
	ncs.hh 		 \
	model-view.hh 		 \
	validation.hh            \
	validation.cc            \
	updating-map-params.hh   \
	updating-coordinates-molecule-parameters.hh \
	ud-colour-rule.hh        \
	view.hh 		 \
	view.cc 		 \
	new-centre.hh            \
        widget-headers.hh        \
	widget-from-builder.hh   \
	widget-from-builder.cc   \
        json.hpp


# 20220228-PE this appeared in a merge - not sure if it's needed (or should be included directly)
# coot_bin_SOURCES += coot_wrap_python.cc

EXTRA_DIST += coot_wrap_python_pre.cc
dist_python_PYTHON = coot.py coot.pth

if OS_WIN32
# Python .pyd modules are simply DLLs, but they have to be called .pyd for
# python to find them, and libtool only creates .dll.
install-exec-hook:
	mv $(pyexecdir)/_coot-0.dll $(pyexecdir)/_coot.pyd
endif

# testcartesian_SOURCES = testcartesian.cc

# graphics_crash_test_SOURCES = graphics-crash-test.cc graphics-crash-test.h


# testcartesian_LDADD = $(COOT_EXTRA_LIBS)

# graphics_crash_test_LDADD = @GTK_LIBS@ $(GTKGL_CONTEXT_LIBS) $(INTLLIBS) $(GLUT_LIBS) -lGL -lGLU

# merge problems I think. fix another time.
#
# coot_fix_nomenclature_errors_LDADD = \
# 	$(top_builddir)/ligand/libcoot-ligand.la      		\
# 	$(top_builddir)/cootaneer/libcootaneer.la 		\
# 	$(top_builddir)/db-main/libcoot-db-main.la      	\
# 	$(top_builddir)/ideal/libcoot-ideal.la    		\
#         $(top_builddir)/coot-utils/libcoot-map-heavy.la         \
# 	$(top_builddir)/geometry/libcoot-geometry.la    	\
# 	$(top_builddir)/coords/libcoot-coords.la      		\
# 	$(top_builddir)/coot-utils/libcoot-coord-extras.la	\
# 	$(top_builddir)/coot-utils/libcoot-coord-utils.la	\
# 	$(top_builddir)/coot-utils/libcoot-map-utils.la		\
# 	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la 	\
#         $(top_builddir)/lidia-core/libcoot-lidia-core.la        \
# 	$(top_builddir)/analysis/libcoot-analysis.la 		\
# 	$(top_builddir)/utils/libcoot-utils.la	  	        \
# 	$(top_builddir)/mini-mol/libcoot-mini-mol.la 		\
# 	$(CLIPPER_LIBS) 	\
# 	$(MMDB_LIBS)		\
# 	$(PYTHON_LIBS)		\
# 	$(RDKIT_LIBS)	        \
# 	$(BOOST_LDFLAGS)        \
# 	$(L_BOOST_PYTHON)

# removed from fix_nomenclature_errors_LDADD:
#        $(RDKIT_LIBS)
#        $(L_BOOST_PYTHON)	# also needs PYTHON_LIBS, I guess?

# coot_bin_LDFLAGS = -static -lstdc++
#
coot_bin_LDFLAGS = $(EXPORT_DYNAMIC_FLAG_SPEC)

# coot_bin_CXXFLAGS=-Wl,--export-dynamic -dynamiclib
#  -Wl,-undefined -Wl,dynamic_lookup

# 20101220 Note: when libtool/ltmain.sh is fixed, OPENMP_CFLAGS will
# not be needed here. So says Justin Lecher (see email).
#
# the -lblas here is a bandaid - the GSL usage of BLAS needs to be
# fixed in configure.
#
libcootsumo_la_LIBADD = $(OPENMP_CFLAGS) \
      $(COOT_EXTRA_LIBS) \
		$(VORBISFILE_LIBS)           \
		$(CLIPPER_LIBS) 	\
		$(LIBSSM_LIBS)		\
		$(CCP4SRS_LIBS) 	\
		$(MMDB_LIBS)		\
		$(MOLTRIANGLES_LIBS)	\
		$(PYTHON_LIBS)		\
		$(GUILE_GTK_LIBS) \
		$(GUILE_LDFLAGS)	\
		$(GLOB_LIBS)		\
		$(GTKGL_CONTEXT_LIBS) \
		$(EPOXY_LIBS)        \
		$(CANVAS_LIBS) 		\
		$(GSL_LIBS) -lblas   \
		$(MYSQL_LIBS)		   \
		$(CURL_LIBS)		   \
	  	$(RDKIT_LIBS)	      \
		$(BOOST_LDFLAGS)     \
		$(L_BOOST_PYTHON)	   \
		$(GOOCANVAS_LIBS)    \
		-lfreetype           \
		$(SQLITE3_LIBS)      \
		@GTK_LIBS@           \
		$(SHARED_LDFLAGS)    \
		-lpthread

# make freetype a configure variable FREETYPE_LIBS

#		$(INTLLIBS) 		Fix when GETTEXT workds

# the filename for libcootsumo.la in LIBADD must be direct, or the Makefile doesn't know how to make it.
test_cholesky_SOURCES = test-chol.cc
test_cholesky_LDADD =	\
        libcootsumo.la                                    \
	$(top_builddir)/coords/libcoot-coords.la          \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la      \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/geometry/libcoot-geometry.la      \
        $(top_builddir)/lidia-core/libcoot-lidia-core.la  \
	$(top_builddir)/utils/libcoot-utils.la            \
        $(RDKIT_LIBS) \
	$(CLIPPER_LIBS) $(MMDB_LIBS) $(GSL_LIBS)

# testmap_SOURCES = testmap.cc xmap-interface.cc
# testmap_LDADD = $(COOT_EXTRA_LIBS) $(SURFACE_LIB) $(CLIPPER_LIBS)


#write_map_SOURCES = write-map.cc xmap-interface.cc
#write_map_LDADD = $(COOT_EXTRA_LIBS) $(SURFACE_LIB)  $(CLIPPER_LIBS)
# example_LDADD =$(GNOME_LIBDIR)$(GNOMEUI_LIBS) $(INTLLIBS)


# test_command_line_SOURCES = test-command-line.cc
# test_command_line_LDADD = -L$(top_srcdir)/compat      -lcoot-compat

# if you like, you can investigate why coot-geometry and clipper is needed here:
test_coot_align_SOURCES = align.cc
test_coot_align_LDADD =	\
        libcootsumo.la                                          \
	$(top_builddir)/geometry/libcoot-geometry.la    	\
	$(top_builddir)/coords/libcoot-coords.la 		\
	$(top_builddir)/mini-mol/libcoot-mini-mol.la 		\
	$(top_builddir)/coot-utils/libcoot-coord-utils.la 	\
	$(top_builddir)/utils/libcoot-utils.la  		\
        $(top_builddir)/lidia-core/libcoot-lidia-core.la        \
        $(GLOB_LIBS)		\
        $(RDKIT_LIBS) 		\
	$(CLIPPER_LIBS) 	\
	$(MMDB_LIBS)

# test_lbg_graph_SOURCES = test-lbg-graph.cc lbg-graph.hh lbg-graph.cc
# test_lbg_graph_LDADD = $(MMDB_LIBS)


coot_density_score_by_residue_bin_SOURCES = density-score-by-residue-src.cc
coot_density_score_by_residue_bin_LDADD  =	\
        libcootsumo.la                                          \
	$(top_builddir)/geometry/libcoot-geometry.la    	\
	$(top_builddir)/coot-utils/libcoot-coord-extras.la  	\
	$(top_builddir)/coords/libcoot-coords.la 		\
	$(top_builddir)/mini-mol/libcoot-mini-mol.la 		\
	$(top_builddir)/ligand/libcoot-ligand.la      		\
	$(top_builddir)/analysis/libcoot-analysis.la 		\
	$(top_builddir)/db-main/libcoot-db-main.la     		\
	$(top_builddir)/ideal/libcoot-ideal.la    		\
	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la 	\
	$(top_builddir)/coot-utils/libcoot-map-utils.la  	\
	$(top_builddir)/coot-utils/libcoot-coord-utils.la  	\
        $(top_builddir)/lidia-core/libcoot-lidia-core.la        \
	$(top_builddir)/utils/libcoot-utils.la  		\
	$(top_builddir)/cootaneer/libcootaneer.la 		\
	$(top_builddir)/compat/libcoot-compat.la  		\
        $(GLOB_LIBS)  		 \
        $(RDKIT_LIBS) 		\
	$(CLIPPER_LIBS) 	\
	$(MMDB_LIBS)		\
	$(GSL_LIBS)

dynarama_bin_SOURCES = dynarama-main.cc
dynarama_bin_LDADD = \
        libcootsumo.la                                          \
        $(top_builddir)/coords/libcoot-coords.la                \
	$(top_builddir)/ligand/libcoot-ligand.la      		\
	$(top_builddir)/coot-utils/libcoot-map-heavy.la  	\
	$(top_builddir)/coot-utils/libcoot-map-utils.la  	\
	$(top_builddir)/coot-utils/libcoot-coord-extras.la  	\
        $(top_builddir)/coot-utils/libcoot-coord-utils.la       \
	$(top_builddir)/ideal/libcoot-ideal.la    		\
	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la 	\
        $(top_builddir)/lidia-core/libcoot-lidia-core.la        \
	$(top_builddir)/mini-mol/libcoot-mini-mol.la    	\
        $(top_builddir)/geometry/libcoot-geometry.la            \
	$(top_builddir)/pli/libcoot-pli.la                      \
	$(top_builddir)/db-main/libcoot-db-main.la     		\
	$(top_builddir)/analysis/libcoot-analysis.la 		\
	$(top_builddir)/cootaneer/libcootaneer.la 		\
        $(top_builddir)/compat/libcoot-compat.la                \
        $(top_builddir)/utils/libcoot-utils.la                  \
        $(CLIPPER_LIBS)         \
        $(MMDB_LIBS)            \
        $(GLOB_LIBS)            \
        $(GOOCANVAS_LIBS)       \
        @GTK_LIBS@
