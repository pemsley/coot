# src/Makefile.am
#
# Copyright 2001, 2002, 2003, 2004, 2005, 2006, 2007 by The University of York
# Copyright 2007, 2008, 2009, 2010, 2011, 2012 by The University of Oxford
# Copyright 2013, 2014, 2015, 2018, 2020, 2021 by Medical Research Council
# Author: Paul Emsley
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation; either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License and
# the GNU Lesser General Public License along with this program; if not,
# write to the Free Software Foundation, Inc., 51 Franklin Street,
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

MAKE_COOT_WRAP_GUILE=$(SHELL) $(top_srcdir)/src/post-process-coot-wrap-guile.sh

MAKE_COOT_WRAP_PYTHON_XXX=$(SHELL) $(top_srcdir)/src/post-process-coot-wrap-python.sh


##########################################################################


# remove this from coot dependencies:
# 	$(top_builddir)/mmtz-extras/libcoot-mmtz-extras.la
#
# we don't have lidia yet.
#        $(top_builddir)/lbg/libcoot-lidia.la

COOT_EXTRA_LIBS = \
	$(top_builddir)/api/libcoot-api.la      		\
	$(top_builddir)/pli/libcoot-pli.la      		\
	$(top_builddir)/layla/libcoot-layla.la                  \
	$(top_builddir)/ligand/libcoot-res-tracer.la  		\
	$(top_builddir)/ligand/libcoot-ligand.la      		\
	$(top_builddir)/validation-graphs/libcoot-validation-graphs.la \
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
	coot_wrap_guile_pre.cc


# subpkgdatadir = $(pkgdatadir)
# subpkgdata_DATA = dynarama.glade geometry_graphs_builder.glade

cootgladedir = $(pkgdatadir)/glade
dist_cootglade_DATA = geometry_graphs_builder.glade

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
		     meshes.shader meshes-for-shadow-map.shader instanced-meshes-for-shadow-map.shader \
		     meshes-with-shadows.shader instanced-meshes-with-shadows.shader \
		     meshes-for-ssao.shader effects.shader 9.ssao_blur.shader 9.ssao_geometry.shader 9.ssao.shader \
		     instanced-meshes-for-ssao.shader texture-meshes-for-ssao.shader texture-meshes-with-shadows.shader \
		     texture-meshes-shadow-map.shader shadow-map-image-texture.shader \
		     rotation-centre-cross-hairs-for-ssao.shader background-image.shader

MOLTRIANGLES_CPPFLAGS=-DUSE_MOLECULES_TO_TRIANGLES

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
                             $(top_srcdir)/src/pythonic-cfc.hh \
                             $(top_srcdir)/src/curlew.h \
                             $(top_srcdir)/src/curlew.hh \
                             $(top_srcdir)/src/validation.hh \
                             $(top_srcdir)/src/cmtz-interface.hh \
                             $(top_srcdir)/src/c-interface-mogul.hh
	swig -o $@ $(MOLTRIANGLES_CPPFLAGS) -DCOOT_USE_GTK2_INTERFACE $(ENHANCED_LIGAND_TOOLS_CPPFLAGS) $(WITH_GUILE_CFLAGS) $(USE_LIBCURL) -I$(top_srcdir)/src -I$(top_srcdir) -guile -c++  $(top_srcdir)/src/coot.i

# For GOOPS support: -proxy -Linkage module (or passive), however, doing that makes SWIG_init() go away.
# -proxy  exports the GOOPS wrapper definitions into module.scm

coot_wrap_guile.cc: coot_wrap_guile_pre.cc
	$(MAKE_COOT_WRAP_GUILE) $^ $@

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
                         $(top_srcdir)/src/c-interface-generic-objects.h \
                         $(top_srcdir)/src/pythonic-cfc.hh \
                         $(top_srcdir)/src/curlew.h \
                         $(top_srcdir)/src/curlew.hh \
                         $(top_srcdir)/src/validation.hh \
                         $(top_srcdir)/src/c-interface-mogul.hh
	swig -o $@ $(MOLTRIANGLES_CPPFLAGS) -DCOOT_USE_GTK2_INTERFACE $(ENHANCED_LIGAND_TOOLS_CPPFLAGS) $(PYTHON_CPPFLAGS)  -DUSE_PYTHON $(USE_LIBCURL) -I$(top_srcdir)/src  -I$(top_srcdir) -shadow -python -c++ $(top_srcdir)/src/coot.i

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
	$(GUILE_CFLAGS)		\
	$(PYTHON_CPPFLAGS)	\
	$(PYTHON_CFLAGS)

# must have GUILE flags before $(INCLUDES), and this
# does not happen because LTCXXCOMPILE has AM_CXXFLAGS after INCLUDES.
# But DEFAULT_INCLUDES is *before* INCLUDES
# DEFAULT_INCLUDES = ... snip.  Rethink...

# get this from configure
# /Users/pemsley/python-3.8/lib/python3.8/site-packages/PyGObject-3.36.1-py3.8-macosx-11.5-x86_64.egg/include/pygobject-3.0/pygobject.h
# /Users/pemsley/test/ccp4-20210929/include/pygobject-3.0/pygobject.h
# PYGOBJECT_CFLAGS=-I$(HOME)/python-3.8/lib/python3.8/site-packages/PyGObject-3.36.1-py3.8-macosx-11.5-x86_64.egg/include

# the value from pkgconfig is wrong: -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -I/home/paule/autobuild/Linux-penelope-gtk4/include/pygobject-3.0'
# the include files are not there. Use python3 install_headers to instal the headers
#                  
# MyPyGObject_CFLAGS=-I/home/paule/autobuild/Linux-penelope-gtk4/lib/python3.9/site-packages/PyGObject-3.42.2-py3.9-linux-x86_64.egg/include
MyPyGObject_CFLAGS=-I$(pyexecdir)/PyGObject-3.42.2-py3.9-linux-x86_64.egg/include


# On merging gtk3 branch I had a merge conflict for MyPyGObject_CFLAGS and PYGOBJECT_CFLAGS - I kept both
AM_CXXFLAGS = \
        -DUSE_MOLECULES_TO_TRIANGLES \
        $(MOLECULESTOTRIANGLES_CXXFLAGS) \
	-DPKGDATADIR='"$(pkgdatadir)"' \
	-DPACKAGE_LOCALE_DIR='"$(pkgdatadir)/locale"' \
	-DPKGPYTHONDIR='"$(pkgpythondir)"' \
	-DPYTHONDIR='"$(pythondir)"' \
	-DXDATADIR='"$(datadir)"'    \
	-DCOOT_SYS_BUILD_TYPE='"$(COOT_SYS_BUILD_TYPE)"' \
	-DCOOT_BUILD_INFO_STRING='"$(COOT_BUILD_INFO_STRING)"' \
    $(openal_CFLAGS)        \
	$(RSVG_CFLAGS)		\
	$(USE_LIBCURL)          \
	$(WITH_GUILE_CFLAGS)	\
	$(GUILE_CFLAGS)		\
        $(MyPyGObject_CFLAGS)   \
	$(PYGOBJECT_CFLAGS)	\
	$(PYTHON_CPPFLAGS)	\
	$(PYTHON_CFLAGS)	\
	$(COOT_WINDOWS_CFLAGS)	\
	$(OPENMP_CFLAGS)


if USE_BACKWARD
   BACKWARD_FLAGS=-DUSE_BACKWARD
endif


# SWIG_FLAGS is something of a misnomer.  It actually defines
# COOT_USE_GTK2_INTERFACE, which is needed to correctly declare the
# correct callbacks in both GTK1 and GTK2 regimes.
# 20230218-PE the webassembly merge: SWIG_FLAGS is removed from AM_CPPFLAGS
#             because we want operator== to work for residue specs.
#
AM_CPPFLAGS = \
	-DUSE_PYTHON \
	@GTK_CFLAGS@			\
	-I$(top_srcdir)                 \
	$(RDKIT_CXXFLAGS)		\
	$(BOOST_CPPFLAGS)		\
	$(ENHANCED_LIGAND_TOOLS_CPPFLAGS) \
	$(GLM_CXXFLAGS)                 \
	$(WITH_GUILE_CFLAGS)		\
	$(GUILE_CFLAGS)		        \
	$(MOLTRIANGLES_CPPFLAGS)	\
	$(LIBSSM_CXXFLAGS)		\
	$(MMDB_CXXFLAGS)		\
	$(CCP4SRS_CXXFLAGS)		\
	$(gsl_CFLAGS)			\
	$(CLIPPER_CXXFLAGS)		\
	$(PYGTK_CFLAGS)                 \
	$(GLOB_CFLAGS)                  \
	$(WITH_SOUND)                   \
	$(VORBISFILE_CFLAGS)            \
	$(SOUND_FLAGS)                  \
	$(SQLITE3_CFLAGS)               \
	$(GEMMI_CXXFLAGS)               \
	$(BACKWARD_FLAGS)               \
        $(CURL_CFLAGS)

#	$(WITH_SCRIPT_CFLAGS)

pkginclude_HEADERS = crunch-model.hh

pyexec_LTLIBRARIES = _coot.la

# the python libraries depend on libcootsumo - don't race
BUILT_SOURCES = libcootsumo.la

# Does this stop the race condition?
# Well, it might do, but if I add it then _coot.so is not updated.
# So, how about I added libcootsumo.la to the SOURCES, will that work?
# _coot.la: libcootsumo.la

_coot_la_SOURCES = coot.cc coot_wrap_python.cc libcootsumo.la

_coot_la_LDFLAGS = -module $(SHARED_LDFLAGS) --no-undefined  \
	   libcootsumo.la                                    \
	   $(top_builddir)/coot-utils/libcoot-coord-utils.la \
	   $(top_builddir)/coords/libcoot-coords.la          \
	   $(top_builddir)/geometry/libcoot-geometry.la      \
	   $(top_builddir)/pli/libcoot-pli.la                \
	   $(top_builddir)/lidia-core/libcoot-lidia-core.la  \
	   $(top_builddir)/utils/libcoot-utils.la            \
	   $(PYTHON_LIBS) $(CLIPPER_LIBS)


if OS_WIN32
# no more libexec for windows
bin_PROGRAMS     = WinCoot coot-density-score-by-residue-bin
else
if OS_DARWIN
libexec_PROGRAMS = MacCoot coot-density-score-by-residue-bin
bin_PROGRAMS =
else
libexec_PROGRAMS = coot-1 coot-density-score-by-residue-bin
bin_PROGRAMS =
endif
endif

dist_bin_SCRIPTS = coot coot-density-score-by-residue

check_PROGRAMS = test-coot-align test-cholesky test-tubular-spline

#coot_fix_nomenclature_errors_SOURCES = \
#	fix-nomenclature.cc    \
#	coot-nomenclature.cc

if OS_WIN32
WinCoot_SOURCES = main.cc glade-callbacks.cc glade-callbacks-main-window.cc \
        glade-callbacks-preferences.cc glade-callbacks-shader-settings.cc \
        new-startup.cc \
        setup-gui-components.cc \
        setup-gui-components.hh \
	create-menu-item-actions.cc \
	create-menu-item-actions.hh
else
if OS_DARWIN
MacCoot_SOURCES = main.cc glade-callbacks.cc glade-callbacks-main-window.cc \
        glade-callbacks-preferences.cc glade-callbacks-shader-settings.cc \
        new-startup.cc \
        setup-gui-components.cc \
        setup-gui-components.hh \
	create-menu-item-actions.cc \
	create-menu-item-actions.hh

else
coot_1_SOURCES = main.cc glade-callbacks.cc glade-callbacks-main-window.cc \
        glade-callbacks-preferences.cc glade-callbacks-shader-settings.cc \
        new-startup.cc \
        setup-gui-components.cc \
        setup-gui-components.hh \
	create-menu-item-actions.cc \
	create-menu-item-actions.hh
endif
endif

# is it a good idea to move functions from coot-bin dependencies into libcootsumo?
# e.g. atoms_with_zero_occupancy()?
if OS_WIN32
WinCoot_LDADD =  \
	libcootsumo.la                                     \
	$(top_builddir)/pli/libcoot-pli.la                 \
	$(top_builddir)/density-contour/libcoot-isosurface.la  	\
	$(top_builddir)/layla/libcoot-layla.la             \
	$(top_builddir)/coot-utils/libcoot-map-utils.la	   \
	$(top_builddir)/coords/libcoot-coords.la	   \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
	$(top_builddir)/lidia-core/libcoot-lidia-core.la   \
	$(top_builddir)/geometry/libcoot-geometry.la       \
	$(top_builddir)/utils/libcoot-utils.la             \
	$(CLIPPER_LIBS) $(CURL_LIBS) $(GUILE_LDFLAGS) $(GTK_LIBS) $(EPOXY_LIBS) $(EMBEDDED_PYTHON_LIBS)
else
if OS_DARWIN
MacCoot_LDADD =  \
	libcootsumo.la                                     \
	$(top_builddir)/pli/libcoot-pli.la                 \
	$(top_builddir)/density-contour/libcoot-isosurface.la  	\
	$(top_builddir)/layla/libcoot-layla.la             \
	$(top_builddir)/coot-utils/libcoot-map-utils.la	   \
	$(top_builddir)/coords/libcoot-coords.la	   \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
	$(top_builddir)/lidia-core/libcoot-lidia-core.la   \
	$(top_builddir)/geometry/libcoot-geometry.la       \
	$(top_builddir)/utils/libcoot-utils.la             \
	$(CLIPPER_LIBS) $(CURL_LIBS) $(GUILE_LDFLAGS) $(GTK_LIBS) $(EPOXY_LIBS) $(EMBEDDED_PYTHON_LIBS)
else
coot_1_LDADD =  \
	libcootsumo.la                                     \
	$(top_builddir)/pli/libcoot-pli.la                 \
	$(top_builddir)/density-contour/libcoot-isosurface.la  	\
	$(top_builddir)/layla/libcoot-layla.la             \
	$(top_builddir)/coot-utils/libcoot-map-utils.la	   \
	$(top_builddir)/coords/libcoot-coords.la	   \
	$(top_builddir)/coot-utils/libcoot-coord-extras.la \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la  \
	$(top_builddir)/lidia-core/libcoot-lidia-core.la   \
	$(top_builddir)/geometry/libcoot-geometry.la       \
	$(top_builddir)/utils/libcoot-utils.la             \
	$(CLIPPER_LIBS) $(CURL_LIBS) $(GUILE_LDFLAGS) $(GTK_LIBS) $(EPOXY_LIBS) $(EMBEDDED_PYTHON_LIBS)
endif
endif


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
	git-revision-count.cc	 \
	scm-boot-guile.hh	 \
	scm-boot-guile.cc	 \
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
	cc-interface-utils.cc    \
	cc-interface-alignment.cc \
	cc-interface-graphics.cc \
	cc-interface-graphics.hh \
	cc-interface-merge.cc    \
	cc-interface-alignment.hh \
	cc-interface-ncs.hh      \
	cc-interface-molecular-representation.cc \
	cc-interface-molecular-symmetry.cc \
	cc-interface-network.hh  \
	cc-interface-mogul.hh    \
	cc-interface-scripting.cc  \
	cc-interface-scripting.hh  \
	cc-interface-subprocess.cc \
	cc-interface-map-utils.cc\
	cc-interface-image.hh    \
	cc-interface-png.cc      \
	cc-interface-validation.cc      \
	cc-interface-user-defined-atom-colours.hh \
	cc-interface-user-defined-atom-colours.cc \
	cc-interface-molecular-representation.hh  \
	canvas-fixes.hh          \
	pythonic-cfc.cc pythonic-cfc.hh            \
	cfc-widgets.hh           \
	cfc-widgets-c-interface.h \
	cfc-2025.cc              \
	cfc-2025.hh              \
	cfc-gui.cc               \
	cfc-gui.hh               \
	change-dir.cc		 \
	change-dir.hh		 \
	cmtz-interface-gui.hh	 \
	cmtz-interface.hh	 \
	cmtz-interface.cc	 \
	clip-hybrid-atom.hh      \
	clipper-to-glm.hh        \
	colour-functions.hh      \
	colour-holder-to-glm.hh  \
	command-arg.hh           \
	command-line.cc	         \
	command-line-extern.cc   \
	command-line-extern.hh   \
	command-line.hh		 \
	coot-application.hh      \
	coot-application.cc      \
	coot-database.hh	 \
	coot-colour.hh		 \
	coot-colour.cc		 \
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
	coot-phi-psi.cc          \
	coot-phi-psi.hh          \
	coot_pythonmodule.cc     \
	coot-preferences.h       \
	coot-references.h        \
	coot-setup-python.cc     \
	coot-setup-python.hh     \
	coot-version.hh          \
	coot_wrap_guile.cc       \
        curl-utils.hh            \
	curlew.h                 \
	curlew.hh                \
	curlew-gtk4.cc           \
	curlew-gtk4.hh           \
	data-pair-remover.hh     \
	density-results-container-t.hh \
	dots-representation.hh   \
	draw.cc                  \
	draw.hh                  \
	draw-2.hh                \
	drag-and-drop.hh         \
	drag-and-drop.cc         \
	draw-generic-display-objects.cc \
	drawn-ghost-molecule-display.hh \
	dynamic-validation.cc    \
	dynamic-validation.hh    \
	eyes.cc eyes.hh		 \
	event-controller-callbacks.hh \
	event-controller-callbacks.cc \
	file-system-utils.cc     \
	file-system-utils.hh     \
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
	gaussian-surface.cc      \
	generic-objects.cc	 \
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
	globjects.h              \
	globularize.cc           \
	glarea_tick_function.hh  \
	globularize.hh           \
	gphl-screen.cc           \
	graphics-info.cc	 \
	graphics-info-defines.cc \
	graphics-info-defines-gtk4.cc \
	graphics-info-draw.cc    \
	graphics-info-draw-utils.cc \
	graphics-info-draw-models.cc \
	graphics-info-enums.hh \
	graphics-info-geometry-graphs.cc \
	graphics-info-graphs.cc	 \
	graphics-info-gui.cc	 \
	graphics-info-init.cc	 \
	graphics-info-intermediate.cc    \
	graphics-info-keybindings.cc     \
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
	graphics-info-restraints-editor.cc \
	graphics-info-scripting.cc	\
	graphics-info-state.cc	 \
	graphics-info-statics.cc \
	graphics-info-superpose.cc	\
	graphics-info-tick-function.cc \
	graphics-info-trackball.cc	\
	graphics-info-utils.cc	 \
	graphics-ligand-view.cc	 \
	graphics-ligand-view.hh	 \
	graphics-ligand-mesh-molecule.cc \
	graphics-ligand-mesh-molecule.hh \
	graphics-info.h 	 \
	gtkglarea-rama-plot.cc   \
	gtkglarea-rama-plot.hh   \
	gtk-manual.cc		 \
	gtk-manual.h		 \
	gtk-widget-conversion-utils.h  \
	gtk-widget-conversion-utils.cc \
	gtk-manual.hh		 \
	gtk-utils.cc \
	gtk-utils.hh \
	guile-fixups.h		 \
	history_list.hh		 \
	intermediate-atom.hh     \
	HUDMesh.cc               \
	HUDMesh.hh               \
	HUDTextureMesh.cc        \
	HUDTextureMesh.hh        \
	Instanced-Markup-Mesh.hh \
	Instanced-Markup-Mesh.cc \
	init-from-gtkbuilder.hh  \
	init-from-gtkbuilder.cc  \
        idles.hh                 \
	interface.h              \
	key-bindings.cc	         \
	key-bindings.hh	         \
	manipulation-modes.hh    \
        labelled-button-info.hh  \
	LinesMesh.cc LinesMesh.hh \
	LigandViewMesh.cc        \
	LigandViewMesh.hh        \
	lights-info.hh           \
	lsq-dialog-values.hh     \
	lbg-interface.hh         \
	lbg-callbacks-outer.cc   \
	lbg-interface.cc         \
	ligand-check.hh          \
        M2T-interface.hh         \
        M2T-interface.cc         \
	Model.cc Model.hh        \
	make-a-dodec.cc make-a-dodec.hh \
	map-statistics.hh               \
	matrix-utils.hh          \
	matrix-utils.cc          \
	map-statistics.hh        \
	Material.hh              \
	map_triangle.hh          \
	Mesh.cc Mesh.hh          \
	Mesh-from-glTF.cc        \
	Mesh-from-graphical-bonds-instanced.cc \
	Mesh-from-graphical-bonds.cc     \
	meshed-generic-display-object.hh \
	meshed-generic-display-object.cc \
	meshed-particle-container.hh     \
	model-composition-statistics.hh  \
	model-molecule-meshes.hh         \
	model-molecule-meshes.cc         \
	molecular-triangles-mesh.hh      \
	molecule-class-info.cc		 \
	molecule-class-info-alignment.cc \
	molecule-class-info-analysis.cc \
	molecule-class-info-bonds.cc    \
	molecule-class-info-build.cc	\
	molecule-class-info-maps.cc	\
	molecule-class-info-maps-2.cc	\
	molecule-class-info-molecular-symmetry.cc \
	molecule-class-info-kdc.cc      \
	molecule-class-info-dots.cc	\
	molecule-class-info-dots-internal.cc	\
	molecule-class-info-other.cc	\
	molecule-class-info.h		\
	molecule-class-info-mutate.cc	\
	molecule-class-info-residues.cc \
	molecule-class-info-refine.cc   \
	molecule-class-info-kk.cc       \
	molecule-class-info-surface.cc	\
	molecule-class-info-string.cc	\
	molecule-class-info-patterson.cc   \
	molecule-class-info-ncs.cc         \
	molecule-class-info-mol-tris.cc \
	molecule-class-info-widget-work.cc \
	molecule-class-info-geometry.cc    \
	molecular-mesh-generator.cc \
	molecular-mesh-generator-mol-tris.cc \
	molecular-mesh-generator.hh \
	mtz-column-auto-read.hh  \
	named-rotamer-score.hh   \
	new-molecule-by-symmetry-matrix.hh \
	nsv.hh                   \
	obj_loader.h obj_loader.cc \
	pandda-inspect.cc        \
	pae-interface.cc         \
	pisa-interface.hh	 \
	pisa-internal.hh         \
	pick.hh pick.cc	         \
	positioned-widgets.h     \
	positron-plot.cc         \
	pre-load.hh              \
	probe-clash-score.cc     \
	probe-clash-score.hh     \
	protein_db-interface.hh  \
	protein_db-interface.cc  \
	pumpkin.cc pumpkin.hh    \
	python-classes.cc	 \
	python-classes.hh	 \
	python-3-interface.hh	 \
	python-scripting-gui.hh  \
	Particle.cc Particle.hh  \
	pulse-data.hh            \
	rail-points.hh           \
        read-molecule.hh         \
	read-phs.h		 \
	read-phs.c		 \
	read-cif.h		 \
	read-cif.cc		 \
	rama_mousey.cc		 \
	rama_mousey.hh		 \
	rama_plot_with_canvas.cc \
	rama_plot_with_canvas.hh \
	text-rendering-utils.cc  \
	text-rendering-utils.hh  \
	restraints-editor.cc     \
	restraints-editor.hh     \
	restraints-editor-c.h    \
	rotate-translate-modes.hh \
	rotamer-search-modes.hh  \
	rsr-functions.cc         \
	rsr-functions.hh         \
	screendump-tga.hh        \
	screendump-tga.cc        \
	sdf-interface.cc         \
	sdf-interface.hh         \
	sdf-internal.hh          \
	select-atom-info.hh      \
	Shader.cc Shader.hh      \
	simple-distance-object.hh \
	sound.cc                 \
	sound.hh                 \
	startup-utils.hh         \
	startup-utils.cc         \
	startup-scripts.hh       \
	startup-scripts.hh       \
	stringify-error-code.hh  \
	testing.cc testing-2.cc  \
	testing-data.hh          \
	testing.hh               \
	Texture.cc Texture.hh    \
	TextureMesh.cc TextureMesh.hh \
	tomo-tools.cc            \
	translation-gizmo.cc     \
	translation-gizmo.hh     \
	triangle-extras.h        \
	user-mods.hh user-mods.cc \
	xmap-utils.cc 		 \
	xmap-utils.h             \
	remarks-browser-gtk-widgets.hh \
	rama-plot-phi-psi.cc     \
	rama-plot-phi-psi.hh     \
	rama-rota-score.hh	 \
	spherical-density.cc     \
	single-map-properties-dialog.cc \
	single-map-properties-dialog.hh \
	fragment-info.hh	 \
	atom-name-bits.hh	 \
	additional-representation.hh	   \
	extra-restraints-representation.hh \
	extra-distance-restraint-markup.hh \
	atom-attribute.hh	 \
	animated-ligand.hh	 \
	atom-selection.hh	 \
	ncs.hh			 \
	model-view.hh		 \
	validation.hh            \
	validation.cc            \
	updating-map-params.hh   \
	updating-coordinates-molecule-parameters.hh \
	ud-colour-rule.hh        \
	vnc-vertex-to-generic-vertex.hh \
	view.hh		 \
	view.cc		 \
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

# graphics_crash_test_LDADD = @GTK_LIBS@ $(INTLLIBS) $(GLUT_LIBS) -lGL -lGLU

# merge problems I think. fix another time.
#
# coot_fix_nomenclature_errors_LDADD = \
# 	$(top_builddir)/ligand/libcoot-ligand.la      		\
# 	$(top_builddir)/cootaneer/libcootaneer.la 		\
# 	$(top_builddir)/db-main/libcoot-db-main.la      	\
# 	$(top_builddir)/ideal/libcoot-ideal.la    		\
#       $(top_builddir)/coot-utils/libcoot-map-heavy.la         \
# 	$(top_builddir)/geometry/libcoot-geometry.la    	\
# 	$(top_builddir)/coords/libcoot-coords.la      		\
# 	$(top_builddir)/coot-utils/libcoot-coord-extras.la	\
# 	$(top_builddir)/coot-utils/libcoot-coord-utils.la	\
# 	$(top_builddir)/coot-utils/libcoot-map-utils.la		\
# 	$(top_builddir)/ccp4mg-utils/libccp4mg-utils.la 	\
#       $(top_builddir)/lidia-core/libcoot-lidia-core.la        \
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
if OS_DARWIN
MacCoot_LDFLAGS = $(EXPORT_DYNAMIC_FLAG_SPEC)
else
coot_1_LDFLAGS = $(EXPORT_DYNAMIC_FLAG_SPEC)
endif

# coot_bin_CXXFLAGS=-Wl,--export-dynamic -dynamiclib
#  -Wl,-undefined -Wl,dynamic_lookup

# 20101220 Note: when libtool/ltmain.sh is fixed, OPENMP_CFLAGS will
# not be needed here. So says Justin Lecher (see email).
#
# the -lblas here is a bandaid - the GSL usage of BLAS needs to be
# fixed in configure.
#
libcootsumo_la_LIBADD = $(OPENMP_CFLAGS) \
                $(COOT_EXTRA_LIBS)      \
		$(VORBISFILE_LIBS)      \
		$(GTK_LIBS)             \
		$(CLIPPER_LIBS) 	\
		$(LIBSSM_LIBS)		\
		$(CCP4SRS_LIBS) 	\
		$(MMDB_LIBS)		\
		$(MOLTRIANGLES_LIBS)	\
		$(PYTHON_LIBS)		\
		$(GUILE_LDFLAGS)      \
		$(GLOB_LIBS)	      \
		$(EPOXY_LIBS)         \
		$(gsl_LIBS) -lblas    \
		$(MYSQL_LIBS)	      \
		$(CURL_LIBS)	      \
		$(RSVG_LIBS)	      \
	  	$(RDKIT_LIBS)	      \
		$(BOOST_LDFLAGS)      \
		$(L_BOOST_PYTHON)     \
		-lfreetype            \
		$(SQLITE3_LIBS)       \
		$(openal_LIBS)        \
		$(SHARED_LDFLAGS)     \
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
	$(CLIPPER_LIBS) $(MMDB_LIBS) $(gsl_LIBS)

test_tubular_spline_SOURCES = test-tubular-spline.cc
test_tubular_spline_LDADD = \
	$(top_builddir)/coot-utils/libcoot-coord-utils.la \
	$(top_builddir)/MoleculesToTriangles/CXXClasses/libMoleculesToTrianglesCXXClasses.la

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
	$(gsl_LIBS)

