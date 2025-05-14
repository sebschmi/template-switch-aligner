# Python Bindings for [`lib_tsalign`](https://crates.io/crates/lib_tsalign)

[![PyPI](https://img.shields.io/pypi/v/tsalign)](https://pypi.org/project/tsalign/)

These bindings are still very minimal and are subject to improvement and/or breaking changes with future versions.

## Usage
Install with `pip install tsalign`.

The most important function is `tsalign.align(reference, query, **settings)`. On the object that is returned, you can e.g. call `.stats()` or `.cigar()`.
