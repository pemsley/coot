# pyrogen/Makefile.am
# 
# Copyright 2012 The University of Oxford
# Copyright 2014 by Medical Research Council
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


# we need files for make dist configured without MAKE_ENHANCED_LIGAND_TOOLS and USE_PYTHON

# pyexec_LTLIBRARIES allows us to make libtool libraries that get the
# file names we like and installed where we like, not forced to be
# libxxx.la in $(libdir) as lib_LTLIBRARIES.
#
pyexec_LTLIBRARIES = 

# rdkit includes are in the rdkit directory, boost includes are the (normal) prefix above that
#
# ENHANCED_LIGAND_TOOLS_CFLAGS is needed to compile with pli pi-stacking header
AM_CPPFLAGS = $(BOOST_CPPFLAGS) $(RDKIT_CXXFLAGS) $(PYTHON_CPPFLAGS) \
              $(ENHANCED_LIGAND_TOOLS_CPPFLAGS) -I$(top_srcdir)

# we don't want "-pthread" which comes in BOOST_CPPFLAGS now
SWIG_CPPFLAGS = $(RDKIT_CXXFLAGS) $(PYTHON_CXXFLAGS) \
                $(ENHANCED_LIGAND_TOOLS_CPPFLAGS) -I$(top_srcdir)

AM_CXXFLAGS = -DPKGDATADIR='"$(pkgdatadir)"' $(CAIRO_CFLAGS) $(COOT_WINDOWS_CFLAGS)

pyrogen_swig.py: pyrogen-swig.cc

mmdb2.py: mmdb2.cc

#_mmdb2_la_SOURCES = mmdb2.cc
#
#_mmdb2_la_LIBADD = ../lidia-core/libcoot-lidia-core.la $(MMDB_LIBS)
#
#_mmdb2_la_LDFLAGS = -module $(SHARED_LDFLAGS)

# EXTRA_DIST = pyrogen-swig.i mmdb2.i mmdb-typemaps.i pointers.i mmdb2-interface.cc
EXTRA_DIST = pyrogen-swig.i

# we can't use AM_CPPFLAGS because that uses PYTHON_CFLAGS which in turn adds
# -iwithsysroot /System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7
# and swig barfs on that. So instead of AM_CPPFLAGS let's make a variable for just swig
pyrogen-swig.cc: $(top_srcdir)/pyrogen/pyrogen-swig.i $(top_srcdir)/pyrogen/restraints.hh
	swig $(SWIG_CPPFLAGS) -o $@ -python -shadow -c++ $<

# weird:
# mmdb2-interface.cc doesn't have its own object code - it is compiled at the end
# of mmdb2.cc - because it relies on mmdb2.cc "header" code that is not in a header.
#    SWIG_ConvertPtr(swigged_lig_res, (void **) &residue_p,
#		   SWIGTYPE_p_mmdb__Residue, SWIG_POINTER_EXCEPTION);
#
# We need to consider where to put this when mmdb2-interface.cc depends on the RDKit.
#
#mmdb2.cc: mmdb2.i mmdb2-interface.cc
#	swig $(AM_CPPFLAGS) -o $@ -python -shadow -c++ $<
#	echo '#include "mmdb2-interface.cc"' >> $@

pyexec_LTLIBRARIES += _pyrogen_swig.la pyrogen_boost.la
# add these when it works geometry-store.cc geometry-store.hh geometry-store-private.hh
_pyrogen_swig_la_SOURCES = pyrogen-swig.cc

## _pyrogen_swig_la_CFLAGS = $(PYTHON_CFLAGS) -I$(top_srcdir)/include $(PYTHON_CPPFLAGS) 

_pyrogen_swig_la_LIBADD = libcoot_py_restraints.la \
                          $(top_builddir)/ideal/libcoot-ideal.la             \
                          $(top_builddir)/coot-utils/libcoot-coord-utils.la  \
                          $(top_builddir)/lidia-core/libcoot-lidia-core.la   \
			  $(top_builddir)/analysis/libcoot-analysis.la       \
			  $(top_builddir)/geometry/libcoot-geometry.la       \
			  $(top_builddir)/utils/libcoot-utils.la             \
                          $(RDKIT_LIBS) $(BOOST_LDFLAGS) $(L_BOOST_PYTHON)   \
                          $(CLIPPER_LIBS) $(MMDB_LIBS) \
			  $(PYTHON_LIBS) $(CAIRO_LIBS)

# What does this do? If we add it, we don't get this message:
#  `_pyrogen_swig.la' is not a standard libtool library name
# pyrogen/Makefile.am:30: did you mean `lib_pyrogen_swig.la'?
# I doubt the need for the other args, so just use -module for now.
#
# _pyrogen_swig_la_LDFLAGS = -avoid-version -module -shared -export-dynamic
#
_pyrogen_swig_la_LDFLAGS = -module $(SHARED_LDFLAGS)

pyrogen_boost_la_LDFLAGS = -module $(SHARED_LDFLAGS)


# These get installed in libdir.
#
lib_LTLIBRARIES=libcoot_py_restraints.la

# this is a bona fide library
#
libcoot_py_restraints_la_SOURCES = mmff-restraints.cc mmff-restraints.hh \
	py-restraints.cc py-restraints.hh restraints-private.hh cairo-molecule.hh \
	cairo-molecule.cc restraints.cc restraints.hh

# 20141014 Add $RDKIT_LIBS for Bill
libcoot_py_restraints_la_LIBADD =  \
                                $(top_builddir)/geometry/libcoot-geometry.la \
                                $(top_builddir)/lidia-core/libcoot-lidia-core.la \
                                $(top_builddir)/utils/libcoot-utils.la       \
                                $(CLIPPER_LIBS) $(MMDB_LIBS) \
                                $(RDKIT_LIBS) $(BOOST_LDFLAGS) $(L_BOOST_PYTHON) \
                                $(PYTHON_LIBS) $(CAIRO_LIBS)

# -no-undefined causes strange/interesting link errors with clang/Mac OS X 10.10.3
# libcoot_py_restraints_la_LDFLAGS = -no-undefined 
libcoot_py_restraints_la_LDFLAGS = $(SHARED_LDFLAGS)

# The name of this starts with lib so it fits with boost::python wrapper
pyrogen_boost_la_SOURCES = restraints-boost.cc restraints.hh

pyrogen_boost_la_LIBADD = libcoot_py_restraints.la \
                                $(top_builddir)/pli/libcoot-pli.la                 \
                                $(top_builddir)/ideal/libcoot-ideal.la             \
                                $(top_builddir)/coot-utils/libcoot-coord-utils.la  \
                                $(top_builddir)/geometry/libcoot-geometry.la       \
                                $(top_builddir)/coords/libcoot-coords.la           \
			        $(top_builddir)/lidia-core/libcoot-lidia-core.la   \
                                $(MMDB_LIBS) $(RDKIT_LIBS) $(L_BOOST_PYTHON)

libpyrogen_boost_la_LIBDFLAGS = -no-undefined

if OS_WIN32
# need to dist pyrogen.bat too
dist_bin_SCRIPTS = pyrogen.bat pyrogen
else
dist_bin_SCRIPTS = pyrogen
endif

dist_pkgpython_PYTHON = pyrogen.py atom_types.py jay_util.py tautomer.py pyrogen_swig.py \
	coot_git.py

# Something here is broken, I think.
# on $ bash autogen.sh, I get:
# pyrogen/Makefile.am:79: OS_WIN32 does not appear in AM_CONDITIONAL
#
if OS_WIN32
# Python .pyd modules are simply DLLs, but they have to be called .pyd for
# python to find them, and libtool only creates .dll.
install-exec-hook:
	mv $(pyexecdir)/_pyrogen_swig-0.dll $(pyexecdir)/_pyrogen_swig.pyd
	mv $(pyexecdir)/pyrogen_boost-0.dll $(pyexecdir)/pyrogen_boost.pyd
endif

# clean needs to remove _pyrogen_swig.so
#
clean-local: 
	rm -f _pyrogen_swig.so

