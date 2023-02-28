---
layout: post
title:  "libcootapi"
date: Sun 26 Feb 2023 02:27:51 GMT
---

It's time to talk about libcootapi.

For about a year now, I've been working with a team from CCP4 to
bring Coot to the web browser - i.e. WebCoot (this project is known
as Moorhen and will have its own post soon).

As you may know, Coot already has an API, and that was used to create
the python module "coot" as well as having a scheme/guile interface.
Howver, to make that interface useful in the web browser (i.e. compilable
to WebAssembly using emscripten) then it needed to be stripped of GTK and OpenGL.

After some testing, I decided to abandon that idea (it was too ugly -
too many #ifdefs) and to create a new API from scratch - one where GTK
and OpenGL are not part of the dependencies. A new directory was
created ("api") and the resulting library is called libcootapi. It is
available in the `gtk3` branch of Coot.

The libcootapi API has been compiled to WebAssembly and has JavaScript
bindings so that it is useful in the web browser and thus is the basis
of Moorhen.

Additionally, this API now has Python bindings (currently generated
using SWIG) so that it is useful from command-line python3 or plugging
into other programs, such as... oh, I don't know... Blender perhaps.

libcootapi is the API that I am now recommending to Python scripters
that want to use Coot functions. Although it does not contain many of
the functions of the `coot` module, I am keen to transfer them into
libcootapi if/when they are needed.

However, with real space refinement, additions, deletions, validation
information, superpositions, rotamers, ligand-fitting and mutations
already in place, you can get quite a lot done as it stands.

```
$ python3
>>> import coot_headless_api
>>> coot = coot_headless_api.molecules_container_t
```

Many useful functions are now available from `coot`:

```
>>> imol     = coot.read_pdb("test.pdb")
>>> imol_map = coot.read_mtz("test.mtz", "FWT", "PHWT", "W", 0, 0)
>>> dca      = coot.density_correlation_analysis(imol, imol_map)
```

The libcootapi documentation is here:

[https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/docs/api/html](https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/docs/api/html)

libcootapi uses mmdb coordinate identifiers (cids) extensively. The mmdb documentation is here:

[https://www.ebi.ac.uk/pdbe/docs/cldoc/object/cl_obj_selfnc.html](https://www.ebi.ac.uk/pdbe/docs/cldoc/object/cl_obj_selfnc.html)

Notes:

1. Compared to the Python 2 version of the API, the arguments to `auto_fit_best_rotamer()` have been changed
2. To run real-space refinement one would now use `refine_residues_using_atom_cid()` or `refine_residues()` or perhaps `refine_residue_range()`


