// Tests of the input and output file format.

#include <coulomb/io.hpp>
#include <coulomb/solver.hpp>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cstddef>  // std::size_t
#include <sstream>  // std::istringstream, std::ostringstream
#include <string>   // std::string, std::getline
#include <vector>   // std::vector

using Catch::Matchers::ContainsSubstring;
using coulomb::Particle;
using coulomb::Results;

namespace {

// Splits text into lines and counts the values on each line.
std::vector<std::size_t> values_per_line(const std::string& text) {
    std::vector<std::size_t> counts;
    std::istringstream lines(text);
    std::string line;
    while (std::getline(lines, line)) {
        std::istringstream fields(line);
        std::string value;
        std::size_t count = 0;
        while (fields >> value) {
            ++count;
        }
        counts.push_back(count);
    }
    return counts;
}

}  // namespace

TEMPLATE_TEST_CASE("Reading valid input", "[io]", float, double) {
    // Blank lines, extra spaces, tabs and Windows line endings are allowed
    std::istringstream input("1 0 0 0\r\n\r\n  -2.5\t1e-3  2 3  \r\n");
    const auto particles = coulomb::read_particles<TestType>(input, "test input");

    REQUIRE(particles.size() == 2);
    CHECK(particles[1].q == TestType(-2.5));
    CHECK(particles[1].x == TestType(1e-3));
    CHECK(particles[1].y == TestType(2));
    CHECK(particles[1].z == TestType(3));
}

TEST_CASE("Reading invalid input reports the line", "[io]") {
    SECTION("too few values") {
        std::istringstream input("1 0 0 0\n1 1 0\n");
        CHECK_THROWS_WITH(coulomb::read_particles<double>(input, "in.txt"),
                          ContainsSubstring("line 2") && ContainsSubstring("in.txt"));
    }
    SECTION("too many values") {
        std::istringstream input("1 0 0 0 9\n");
        CHECK_THROWS_WITH(coulomb::read_particles<double>(input, "in.txt"),
                          ContainsSubstring("line 1"));
    }
    SECTION("text instead of a number") {
        std::istringstream input("1 0 0 0\n\nabc 0 0 0\n");
        CHECK_THROWS_WITH(coulomb::read_particles<double>(input, "in.txt"),
                          ContainsSubstring("line 3"));
    }
    SECTION("a number too large for float") {
        std::istringstream input("1e40 0 0 0\n");
        CHECK_THROWS(coulomb::read_particles<float>(input, "in.txt"));
    }
}

TEST_CASE("A missing input file is an error", "[io]") {
    CHECK_THROWS_WITH(coulomb::read_particles_from_file<double>("this_file_does_not_exist.txt"),
                      ContainsSubstring("cannot open input file"));
}

TEST_CASE("Output has 4 columns with the potential and 3 without", "[io]") {
    Results<double> results;
    results.Fx = {1, 2};
    results.Fy = {3, 4};
    results.Fz = {5, 6};
    results.phi = {7, 8};

    std::ostringstream with_potential;
    coulomb::write_results<double, true>(with_potential, results);
    CHECK(values_per_line(with_potential.str()) == std::vector<std::size_t>{4, 4});

    std::ostringstream without_potential;
    coulomb::write_results<double, false>(without_potential, results);
    CHECK(values_per_line(without_potential.str()) == std::vector<std::size_t>{3, 3});
}

TEMPLATE_TEST_CASE("Written values read back exactly", "[io]", float, double) {
    // Values with many digits, which a default 6-digit output would round
    Results<TestType> results;
    results.Fx = {TestType(1) / TestType(3)};
    results.Fy = {TestType(-2) / TestType(7)};
    results.Fz = {TestType(1e-20)};
    results.phi = {TestType(123456.789)};

    std::ostringstream out;
    coulomb::write_results<TestType, true>(out, results);

    std::istringstream in(out.str());
    TestType fx{}, fy{}, fz{}, phi{};
    REQUIRE(static_cast<bool>(in >> fx >> fy >> fz >> phi));
    CHECK(fx == results.Fx[0]);  // exactly equal, not just close
    CHECK(fy == results.Fy[0]);
    CHECK(fz == results.Fz[0]);
    CHECK(phi == results.phi[0]);
}
