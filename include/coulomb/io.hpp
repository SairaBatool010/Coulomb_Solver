// Reading particles and writing results as text.
//
// Input:  one particle per line:  q x y z
// Output: one particle per line:  Fx Fy Fz [phi]
//
// The functions that do the work take any stream, not only files, so the
// tests can feed them text directly. The *_file functions wrap them for files.

#pragma once

#include <coulomb/solver.hpp>

#include <concepts>   // std::floating_point
#include <cstddef>    // std::size_t
#include <fstream>    // std::ifstream, std::ofstream
#include <iomanip>    // std::setprecision, std::scientific
#include <istream>    // std::istream
#include <limits>     // std::numeric_limits
#include <ostream>    // std::ostream
#include <sstream>    // std::istringstream
#include <stdexcept>  // std::runtime_error
#include <string>     // std::string, std::getline, std::to_string
#include <vector>     // std::vector

namespace coulomb {

// Reads all particles from a stream. source_name is only used in error
// messages. Throws if a line is not exactly four numbers.
template <std::floating_point Real>
std::vector<Particle<Real>> read_particles(std::istream& in, const std::string& source_name) {
    std::vector<Particle<Real>> particles;
    std::string line;
    std::size_t line_number = 0;

    while (std::getline(in, line)) {
        ++line_number;

        // Skip blank lines (\r handles Windows line endings)
        if (line.find_first_not_of(" \t\r") == std::string::npos) {
            continue;
        }

        // Treat this one line as a small stream and read its four numbers
        std::istringstream fields(line);
        Particle<Real> p{};
        std::string extra;
        if (!(fields >> p.q >> p.x >> p.y >> p.z) || (fields >> extra)) {
            throw std::runtime_error("line " + std::to_string(line_number) +
                                     " of " + source_name + " is not 'q x y z'");
        }
        particles.push_back(p);
    }

    if (in.bad()) {
        throw std::runtime_error("failed while reading " + source_name);
    }
    return particles;
}

template <std::floating_point Real>
std::vector<Particle<Real>> read_particles_from_file(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("cannot open input file: " + path);
    }
    return read_particles<Real>(in, path);
}

// Writes one line per particle: Fx Fy Fz, plus phi if it is enabled.
template <std::floating_point Real, bool WithPotential>
void write_results(std::ostream& out, const Results<Real>& results) {
    // Scientific notation with enough digits that reading a value back gives
    // exactly the same number. One digit sits before the decimal point, which
    // is why this is max_digits10 - 1. The digit count follows the type.
    out << std::scientific
        << std::setprecision(std::numeric_limits<Real>::max_digits10 - 1);

    for (std::size_t i = 0; i < results.Fx.size(); ++i) {
        out << results.Fx[i] << ' ' << results.Fy[i] << ' ' << results.Fz[i];
        if constexpr (WithPotential) {
            out << ' ' << results.phi[i];
        }
        out << '\n';
    }
}

template <std::floating_point Real, bool WithPotential>
void write_results_to_file(const std::string& path, const Results<Real>& results) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("cannot open output file: " + path);
    }

    write_results<Real, WithPotential>(out, results);

    out.close();
    if (!out) {
        throw std::runtime_error("failed while writing output file: " + path);
    }
}

}  // namespace coulomb
