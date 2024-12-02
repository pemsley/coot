---
layout: post
title:  "Doxygen documentation to Python Docstrings"
date: Thu 31 Oct 09:38:05 GMT 2024
---

[Libcootapi](https://www2.mrc-lmb.cam.ac.uk/personal/pemsley/coot/docs/api/html/) is written
in C++ and documented with [Doxygen](https://www.doxygen.nl/).

I wanted this functionality to be available via Python and I chose
[nanobind](https://nanobind.readthedocs.io/en/latest/) as the means to do that.

Nanobind provides access to practically all of the many hundreds of
functions in the C++ class and it would be useful is the Python API
was as well documented as was the C++ api.  Transfering the
documentation from C++ to Python "by hand" would be tedious. So
Lucrezia (@lulu_catapano) and I wrote a [Python
script](https://github.com/pemsley/coot/blob/main/api/doxy-sphinx/xml-to-python.py)
to transfer the Doxygen documentation into Python docstrings in reStructured Text format,
which then get processed by Sphinx so that the Python API looks like read-the-docs
documentation from native Python code.

As well as attractive HTML output, the intermediate documented python
function stubs work with LSP and intellisense to provide api documentation
and completions.

The script works pretty well for our [chapi documentation](https://www.mrc-lmb.cam.ac.uk/lucrezia/libcootapi-documentation/api.html).
It might work for other nanobinders in a similar situation.

