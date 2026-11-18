# Classical Coulomb Solver

A C++20 solver that computes, for N charged particles, the Coulomb force on
each particle and, optionally, the potential at each particle and the total
Coulomb energy of the system. All pairwise interactions are summed directly.

## Build

Requires CMake 3.20 or newer and a C++20 compiler (tested with GCC 14.2).

```sh
cmake -S . -B build -G Ninja
cmake --build build
```

The first configure downloads Catch2 for the tests. The default build type is
`Release` (`-O3`).

## Run

```sh
build/coulomb_solver input.txt output.txt
```

The computation time, and the Coulomb energy if enabled, are printed to the
console. Only the computation is timed, not reading or writing the files.

**Input:** one particle per line, `q x y z` (charge, then position).
Blank lines are ignored. A malformed line stops the program with its line number.

**Output:** one particle per line, `Fx Fy Fz` followed by `phi` if the
potential is enabled. Values are written with enough digits to be read back
exactly.

## Configuration

The settings are chosen at compile time with CMake options:

| Option                   | Values             | Default  |
|--------------------------|--------------------|----------|
| `COULOMB_PRECISION`      | `float` / `double` | `double` |
| `COULOMB_WITH_POTENTIAL` | `ON` / `OFF`       | `ON`     |
| `COULOMB_WITH_ENERGY`    | `ON` / `OFF`       | `ON`     |
| `COULOMB_BUILD_TESTS`    | `ON` / `OFF`       | `ON`     |

For example, a single-precision build without potentials:

```sh
cmake -S . -B build-float -G Ninja -DCOULOMB_PRECISION=float -DCOULOMB_WITH_POTENTIAL=OFF
cmake --build build-float
```

CMake writes the chosen settings into a generated header, `coulomb/config.hpp`,
as a C++ type alias and `constexpr` constants, so no preprocessor macros are
needed. The solver is a template over the floating-point type, and disabled
features are removed at compile time with `if constexpr`.

## Tests

```sh
ctest --test-dir build
```

The Catch2 tests run every physics check in both `float` and `double`:
hand-calculated systems, the forces adding up to zero, the energy against the
direct double sum, independence of particle order, the compile-time switches,
and the input/output format.

## Project layout

```
include/coulomb/solver.hpp      Particle, Results and compute()
include/coulomb/io.hpp          reading particles and writing results
include/coulomb/config.hpp.in   template for the generated settings header
src/main.cpp                    command-line program
tests/                          Catch2 tests
input.txt                       small example input (3 particles)
```
