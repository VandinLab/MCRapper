# MatteoUtils {#mainpage}

A C++17 header-only library offering some utility functions.

For any question or to report bugs, contact Matteo Riondato <riondato@acm.org>.

## Requirements

MatteoUtils requires [Boost][boost]. Some components optionally leverage
parallelism using [OpenMP][openmp].

MatteoUtils uses [CMake][cmake] for installation and to build and run the test.
A C++1z compiler is needed to build the tests.

[Doxygen][doxygen] is needed to build the documentation.

## Installation

MatteoUtils is a header-only library, so there is no compilation step. Run the
following command to install the header files under `./include/matteoutils`:

```sh
mkdir build && cd build && cmake ../ && make install
```

See the CMake documentation for how to specify a different installation prefix,
or just move the include files wherever you want.

### Building the documentation

Run the following command to build and install the HTML documentation under
`./docs/matteoutils/html/`:

```sh
mkdir build-docs && cd build-docs && cmake -DMATTEOUTILS_BUILD_DOCS=On ../ && make docs
```

Point a browser to `./docs/html/index.html` to start browsing the docs.

### Building and running the tests

Run the following command to build the tests:

```sh
mkdir build-tests && cd build-tests && cmake -DBUILD_MATTEOUTILS_TESTS=On ../ && make
```

Run `make test` to run the tests.

## Using MatteoUtils

Just include the MatteoUtils headers in your source file and go!

### Enabling invariants checking

By default MatteoUtils does _not_ check for the validity of function arguments.
Checking for arguments invariants can be done by defining the preprocessor
identifier `MATTEO_INVARIANTS`.

The CMake option `MATTEOUTILS_ENABLE_INVARIANTS` can be set to `On` to enable
invariant checking when MatteoUtils is a subproject of another CMake project.

### Disabling parallelism

Parallelism is enabled when specifying appropriate OpenMP flags when
compiling sources that use MatteoUtils.

If you include MatteoUtils as a subproject of a CMake project, MatteoUtils
will use OpenMP for parallelism by default when your compiler supports it. You
can disable the OpenMP parallelism using the option `MATTEOUTILS_USE_OPENMP`
(set it to `Off`).

## License and Copyright

Copyright 2017 Matteo Riondato <riondato@acm.org>

See the [LICENSE](./LICENSE) and the [NOTICE](./NOTICE) files for licensing and
copyright information.

[boost]: http://www.boost.org/
[cmake]: http://www.cmake.org/
[doxygen]: http://www.stack.nl/~dimitri/doxygen/
[openmp]: http://www.openmp.org/
