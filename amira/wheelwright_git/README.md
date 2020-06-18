# Wheelwright {#mainpage}

A C++17 header-only library offering algorithms and data structures related to
Rademacher averages.

Examples of algorithms are the procedures for computing upper bounds to the
maximum deviation of sample averages from their expectations, and upper bounds
to the empirical Rademacher averages.

Example of data structures are data structures for the projection (set of
vectors) of a family of functions onto a sample, and for the vectors used in the
projection.

For any question or to report bugs, contact Matteo Riondato <riondato@acm.org>.

## Name

_Wheelwright_ is the English translation of the German word _rademacher_. A
wheelwright is a craftsman who builds or repairs wooden wheels.

## Requirements

MatteoUtils requires [Boost][boost] and [MatteoUtils][matteoutils]. The latter
is included as a submodule, so if you have a git clone of the Wheelwright
repository, remember to run

```sh
git submodule update --init
```
to init and fetch the submodule.

Some components optionally leverage parallelism using [OpenMP][openmp].

Some optional but important components require [NLOpt][nlopt], possibly the [git
version][nloptgit].

Wheelwright uses [CMake][cmake] for installation and to build and run the test.
A C++1z compiler is needed to build the tests.

[Doxygen][doxygen] is needed to build the documentation.

## Installation

TODO README.md : write installation instructions

### Building the documentation

Run the following command to build and install the HTML documentation under
`./docs/wheelwright/html/`:

```sh
mkdir build-docs && cd build-docs && cmake -DWHEELWRIGHT_BUILD_DOCS=On ../ && make docs
```

Point a browser to `./docs/html/index.html` to start browsing the docs.

### Building and running the tests

Run the following command to build the tests:

```sh
mkdir build-tests && cd build-tests && cmake -DWHEELWRIGHT_BUILD_TESTS=On ../ && make
```

Run `make test` to run the tests.

# Using Wheelwright

TODO (README.md): write information about using the library (including options)

### Disabling parallelism

Parallelism is enabled when specifying appropriate OpenMP flags when
compiling sources that use Wheelwright.

If you include Wheelwright as a subproject of a CMake project, Wheelwright
will use OpenMP for parallelism by default when your compiler supports it. You
can disable the OpenMP parallelism using the option `WHEELWRIGHT_USE_OPENMP`
(set it to `Off`).

## License and Copyright

Copyright 2017 Matteo Riondato <riondato@acm.org>

See the [LICENSE](./LICENSE) and the [NOTICE](./NOTICE) files for licensing and
copyright information.

[boost]: http://www.boost.org/
[cmake]: https://www.cmake.org/
[doxygen]: http://www.stack.nl/~dimitri/doxygen/
[matteoutils]: https://github.com/rionda/matteoutils/
[nlopt]: http://ab-initio.mit.edu/wiki/index.php/NLopt
[nloptgit]: https://github.com/stevengj/nlopt/
[openmp]: http://www.openmp.org/
