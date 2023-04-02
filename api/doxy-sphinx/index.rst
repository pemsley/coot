.. libcootapi documentation master file, created by
   sphinx-quickstart on Tue Nov 15 01:24:41 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Libcootapi's Documentation
==========================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

Libcootapi is a new (as of 2022) interface to the functionality of
*Coot*. It is designed to be (is intended to be) a clean, consistent and easy-to-use
interface primarily targeting WebAssembly/JavaScript and secondarily, Python.

The previous API for Coot still exists and is available for use with Python and
Guile/Scheme (it provides the ``coot`` module). That interface is considerably more
extensive than this one (consisting of several thousand API functions), but is
significantly less easy to use as it embeds OpenGL and GTK libraries (and all of their
dependencies, of course).

The Molecules Container
=======================

*Coot*'s molecules are referred to by index and should only be accessed *via* member
functions of the ``molecules_container_t`` class.

The functions in this class for the most part return simple types. There are some functions
that return a more complex type, such as ``coot::simple_mesh_t`` or ``coot::validation_information_t``.

Note: The private types, functions, attributes and members are listed here, but they
are not, for the most part, useful for exporting. Which is not to rule out that there
may be *something* there that could usefully be declared as public.

.. doxygenfile:: molecules_container.hh
   :project: libcootapi


Mesh Objects
============

.. doxygenfile:: coot-utils/g_triangle.hh
     :project: libcootapi

.. doxygenfile:: coot-utils/vertex.hh
     :project: libcootapi

.. doxygenfile:: coot-utils/simple-mesh.hh
     :project: libcootapi

.. doxygenfile:: instancing.hh
     :project: libcootapi

Atom and Residue Specifiers
===========================

.. doxygenfile:: geometry/residue-and-atom-specs.hh
      :project: libcootapi

Validation Information
======================

.. doxygenfile:: validation-information.hh
      :project: libcootapi

.. doxygenfile:: residue-validation-information.hh
      :project: libcootapi

.. doxygenfile:: coot-utils/coot-density-stats.hh
      :project: libcootapi

Superposition
=============

.. doxygenfile:: superpose-results.hh
      :project: libcootapi


3D lines
========

.. doxygenfile:: generic-3d-lines.hh
      :project: libcootapi

Coordinates
===========

.. doxygenfile:: coords/Cartesian.h
      :project: libcootapi

Symmetry
========

.. doxygenfile:: coords/mmdb-crystal.h
      :project: libcootapi

Structure Factors
=================

.. doxygenfile:: coot-utils/sfcalc-genmap.hh
      :project: libcootapi


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
