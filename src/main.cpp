// Command-line program for the Coulomb solver.
//
// Usage: coulomb_solver <input file> <output file>
//
// Workflow (section 2.1 of the assignment): read the input, start the clock,
// compute, stop the clock, write the output, print the time.
//
// The precision and features are chosen when configuring with CMake (see
// CMakeLists.txt), which generates coulomb/config.hpp.

#include <coulomb/config.hpp>
#include <coulomb/io.hpp>
#include <coulomb/solver.hpp>

#include <chrono>     // std::chrono (timing)
#include <concepts>   // std::floating_point
#include <exception>  // std::exception
#include <iomanip>    // std::setprecision
#include <iostream>   // std::cout, std::cerr
#include <limits>     // std::numeric_limits
#include <string>     // std::string
#include <vector>     // std::vector

namespace {

// The whole workflow for one configuration.
template <std::floating_point Real, bool WithPotential, bool WithEnergy>
void run(const std::string& input_path, const std::string& output_path) {
    // 1. Read the input (not timed)
    const std::vector<coulomb::Particle<Real>> particles =
        coulomb::read_particles_from_file<Real>(input_path);

    // 2. Compute (timed)
    const auto start = std::chrono::steady_clock::now();
    const coulomb::Results<Real> results =
        coulomb::compute<Real, WithPotential, WithEnergy>(particles);
    const auto stop = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsed = stop - start;

    // 3. Write the output (not timed)
    coulomb::write_results_to_file<Real, WithPotential>(output_path, results);

    // 4. Report on the console
    std::cout << "Precision:        " << sizeof(Real) * 8 << " bit\n";
    std::cout << "Potential:        " << (WithPotential ? "on" : "off") << "\n";
    std::cout << "Particles:        " << particles.size() << "\n";
    if constexpr (WithEnergy) {
        std::cout << "Coulomb energy:   "
                  << std::setprecision(std::numeric_limits<Real>::max_digits10)
                  << results.Ec << "\n";
    }
    std::cout << "Computation time: " << std::setprecision(6)
              << elapsed.count() << " s\n";
}

}  // namespace

int main(int argc, char* argv[]) {
    // argv[0] is the program name, so two arguments means argc == 3
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input file> <output file>\n";
        return 1;
    }

    try {
        using namespace coulomb::config;
        run<Precision, kWithPotential, kWithEnergy>(argv[1], argv[2]);
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << "\n";
        return 1;
    }

    return 0;
}
