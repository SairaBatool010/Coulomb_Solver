# Build file for the Coulomb solver.
#
# Configure and build (from the project folder):
#   cmake -S . -B build -G Ninja
#   cmake --build build
#
# Run the tests:
#   ctest --test-dir build
#
# Choose the configuration with -D options when configuring, for example:
#   cmake -S . -B build-float -G Ninja -DCOULOMB_PRECISION=float -DCOULOMB_WITH_POTENTIAL=OFF

cmake_minimum_required(VERSION 3.20)
project(CoulombSolver LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)  # plain C++20, no GNU extensions

# Build with optimization unless another build type is requested
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# ------------------------------ Options ---------------------------------
set(COULOMB_PRECISION "double" CACHE STRING "Floating-point type: float or double")
set_property(CACHE COULOMB_PRECISION PROPERTY STRINGS float double)
option(COULOMB_WITH_POTENTIAL "Compute and write the potentials" ON)
option(COULOMB_WITH_ENERGY "Compute and print the Coulomb energy" ON)
option(COULOMB_BUILD_TESTS "Build the Catch2 tests" ON)

if(NOT COULOMB_PRECISION MATCHES "^(float|double)$")
    message(FATAL_ERROR "COULOMB_PRECISION must be float or double, not '${COULOMB_PRECISION}'")
endif()

# CMake options are ON/OFF, but the generated C++ header needs true/false
foreach(feature POTENTIAL ENERGY)
    if(COULOMB_WITH_${feature})
        set(COULOMB_WITH_${feature}_CXX true)
    else()
        set(COULOMB_WITH_${feature}_CXX false)
    endif()
endforeach()

# Writes the chosen settings into a real C++ header (no preprocessor macros)
configure_file(
    include/coulomb/config.hpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/generated/coulomb/config.hpp
    @ONLY)

# ------------------------------ Compiler --------------------------------
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    set(COULOMB_WARNINGS -Wall -Wextra -Wpedantic -Wconversion -Wdouble-promotion)
elseif(MSVC)
    set(COULOMB_WARNINGS /W4)
endif()

# MinGW: put the C++ runtime inside each .exe, so a different libstdc++ DLL
# found on PATH can never be loaded by mistake
if(MINGW)
    add_link_options(-static)
endif()

# ------------------------------ Targets ---------------------------------
# The solver is header-only (templates), so this library has no .cpp files.
# Linking to it gives a target the include path and C++20.
add_library(coulomb INTERFACE)
target_include_directories(coulomb INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_features(coulomb INTERFACE cxx_std_20)

add_executable(coulomb_solver src/main.cpp)
target_include_directories(coulomb_solver PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/generated)
target_link_libraries(coulomb_solver PRIVATE coulomb)
target_compile_options(coulomb_solver PRIVATE ${COULOMB_WARNINGS})

# ------------------------------- Tests ----------------------------------
if(COULOMB_BUILD_TESTS)
    include(FetchContent)
    FetchContent_Declare(
        Catch2
        URL https://github.com/catchorg/Catch2/archive/refs/tags/v3.8.1.tar.gz
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
    FetchContent_MakeAvailable(Catch2)

    enable_testing()

    # The tests check both float and double, whatever the options above say
    add_executable(coulomb_tests tests/test_solver.cpp tests/test_io.cpp)
    target_link_libraries(coulomb_tests PRIVATE coulomb Catch2::Catch2WithMain)
    target_compile_options(coulomb_tests PRIVATE ${COULOMB_WARNINGS})

    list(APPEND CMAKE_MODULE_PATH ${catch2_SOURCE_DIR}/extras)
    include(Catch)
    catch_discover_tests(coulomb_tests)
endif()
