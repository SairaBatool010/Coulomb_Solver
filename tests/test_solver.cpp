// Tests of the physics in coulomb::compute.
//
// Every TEMPLATE_TEST_CASE runs twice: once with float and once with double.

#include <coulomb/solver.hpp>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>  // std::max, std::reverse
#include <cmath>      // std::abs, std::sqrt
#include <cstddef>    // std::size_t
#include <limits>     // std::numeric_limits
#include <random>     // std::mt19937, std::uniform_real_distribution
#include <vector>     // std::vector

using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;
using coulomb::Particle;

namespace {

// Relative tolerance for results that should agree up to rounding:
// about 1.2e-4 for float and 2.2e-13 for double.
template <typename Real>
constexpr double tolerance() {
    return 1000.0 * static_cast<double>(std::numeric_limits<Real>::epsilon());
}

// Checks a computed value against an expected value from a hand calculation.
template <typename Real>
void check_close(Real actual, double expected) {
    if (expected == 0.0) {
        CHECK_THAT(static_cast<double>(actual), WithinAbs(0.0, tolerance<Real>()));
    } else {
        CHECK_THAT(static_cast<double>(actual), WithinRel(expected, tolerance<Real>()));
    }
}

// A reproducible random system: the same seed always gives the same particles.
template <typename Real>
std::vector<Particle<Real>> random_system(std::size_t n, Real min_charge) {
    std::mt19937 rng(12345);
    std::uniform_real_distribution<Real> charge(min_charge, Real{1});
    std::uniform_real_distribution<Real> coordinate(Real{-10}, Real{10});

    std::vector<Particle<Real>> particles;
    for (std::size_t i = 0; i < n; ++i) {
        // Braced lists are evaluated left to right, so this order is fixed
        particles.push_back({charge(rng), coordinate(rng), coordinate(rng), coordinate(rng)});
    }
    return particles;
}

// Equation 5 exactly as written (double sum), computed in double precision.
// This is an independent reference for the solver's energy.
template <typename Real>
double energy_double_sum(const std::vector<Particle<Real>>& particles) {
    double Ec = 0.0;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        for (std::size_t j = 0; j < particles.size(); ++j) {
            if (i == j) {
                continue;
            }
            const double dx = static_cast<double>(particles[i].x) - static_cast<double>(particles[j].x);
            const double dy = static_cast<double>(particles[i].y) - static_cast<double>(particles[j].y);
            const double dz = static_cast<double>(particles[i].z) - static_cast<double>(particles[j].z);
            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            Ec += 0.5 * static_cast<double>(particles[i].q) * static_cast<double>(particles[j].q) / r;
        }
    }
    return Ec;
}

}  // namespace

// ------------------------ Hand-calculated systems ------------------------

TEMPLATE_TEST_CASE("Two equal charges repel", "[solver]", float, double) {
    const std::vector<Particle<TestType>> particles = {
        {1, 0, 0, 0},
        {1, 1, 0, 0},
    };
    const auto results = coulomb::compute<TestType, true, true>(particles);

    check_close(results.Fx[0], -1.0);
    check_close(results.Fy[0], 0.0);
    check_close(results.Fz[0], 0.0);
    check_close(results.Fx[1], 1.0);
    check_close(results.Fy[1], 0.0);
    check_close(results.Fz[1], 0.0);
    check_close(results.phi[0], 1.0);
    check_close(results.phi[1], 1.0);
    check_close(results.Ec, 1.0);
}

TEMPLATE_TEST_CASE("Opposite charges attract", "[solver]", float, double) {
    const std::vector<Particle<TestType>> particles = {
        {1, 0, 0, 0},
        {-1, 1, 0, 0},
    };
    const auto results = coulomb::compute<TestType, true, true>(particles);

    check_close(results.Fx[0], 1.0);   // pulled towards +x
    check_close(results.Fx[1], -1.0);  // pulled towards -x
    check_close(results.phi[0], -1.0);
    check_close(results.phi[1], 1.0);
    check_close(results.Ec, -1.0);
}

TEMPLATE_TEST_CASE("Worked example: charges 2 and 3 along z", "[solver]", float, double) {
    const std::vector<Particle<TestType>> particles = {
        {2, 0, 0, 0},
        {3, 0, 0, 2},
    };
    const auto results = coulomb::compute<TestType, true, true>(particles);

    check_close(results.Fx[0], 0.0);
    check_close(results.Fy[0], 0.0);
    check_close(results.Fz[0], -1.5);
    check_close(results.Fz[1], 1.5);
    check_close(results.phi[0], 1.5);
    check_close(results.phi[1], 1.0);
    check_close(results.Ec, 3.0);
}

TEMPLATE_TEST_CASE("Three particles on a line", "[solver]", float, double) {
    const std::vector<Particle<TestType>> particles = {
        {1, 0, 0, 0},
        {1, 1, 0, 0},
        {2, 2, 0, 0},
    };
    const auto results = coulomb::compute<TestType, true, true>(particles);

    check_close(results.Fx[0], -1.5);
    check_close(results.Fx[1], -1.0);
    check_close(results.Fx[2], 2.5);
    check_close(results.phi[0], 2.0);
    check_close(results.phi[1], 3.0);
    check_close(results.phi[2], 1.5);
    check_close(results.Ec, 4.0);
}

TEMPLATE_TEST_CASE("Empty and single-particle systems", "[solver]", float, double) {
    SECTION("no particles") {
        const auto results = coulomb::compute<TestType, true, true>({});
        CHECK(results.Fx.empty());
        CHECK(results.phi.empty());
        CHECK(results.Ec == TestType{0});
    }
    SECTION("one particle feels nothing") {
        const auto results = coulomb::compute<TestType, true, true>({{5, 1, 2, 3}});
        REQUIRE(results.Fx.size() == 1);
        CHECK(results.Fx[0] == TestType{0});
        CHECK(results.Fy[0] == TestType{0});
        CHECK(results.Fz[0] == TestType{0});
        CHECK(results.phi[0] == TestType{0});
        CHECK(results.Ec == TestType{0});
    }
}

// ------------------------- Physical properties --------------------------

TEMPLATE_TEST_CASE("Forces add up to zero (Newton's third law)", "[solver]", float, double) {
    const auto particles = random_system<TestType>(200, TestType{-1});
    const auto results = coulomb::compute<TestType, true, true>(particles);

    double sum_x = 0, sum_y = 0, sum_z = 0, sum_magnitude = 0;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        sum_x += static_cast<double>(results.Fx[i]);
        sum_y += static_cast<double>(results.Fy[i]);
        sum_z += static_cast<double>(results.Fz[i]);
        sum_magnitude += std::abs(static_cast<double>(results.Fx[i])) +
                         std::abs(static_cast<double>(results.Fy[i])) +
                         std::abs(static_cast<double>(results.Fz[i]));
    }

    // Zero up to rounding, measured against the size of the forces
    const double allowed = tolerance<TestType>() * sum_magnitude;
    CHECK_THAT(sum_x, WithinAbs(0.0, allowed));
    CHECK_THAT(sum_y, WithinAbs(0.0, allowed));
    CHECK_THAT(sum_z, WithinAbs(0.0, allowed));
}

TEMPLATE_TEST_CASE("Energy from potentials matches the double sum", "[solver]", float, double) {
    // Positive charges only, so the energy terms cannot cancel each other
    const auto particles = random_system<TestType>(200, TestType(0.1));
    const auto results = coulomb::compute<TestType, true, true>(particles);

    CHECK_THAT(static_cast<double>(results.Ec),
               WithinRel(energy_double_sum(particles), tolerance<TestType>()));
}

TEMPLATE_TEST_CASE("Results do not depend on particle order", "[solver]", float, double) {
    const auto particles = random_system<TestType>(100, TestType{-1});
    auto reversed = particles;
    std::reverse(reversed.begin(), reversed.end());

    const auto a = coulomb::compute<TestType, true, true>(particles);
    const auto b = coulomb::compute<TestType, true, true>(reversed);

    // Compare against the largest value, because single components can be
    // close to zero after cancellation
    double scale = 0;
    for (std::size_t i = 0; i < particles.size(); ++i) {
        scale = std::max({scale, std::abs(static_cast<double>(a.Fx[i])),
                          std::abs(static_cast<double>(a.Fy[i])),
                          std::abs(static_cast<double>(a.Fz[i])),
                          std::abs(static_cast<double>(a.phi[i]))});
    }
    const double allowed = tolerance<TestType>() * scale;

    const std::size_t n = particles.size();
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t k = n - 1 - i;  // where particle i ended up
        CHECK_THAT(static_cast<double>(b.Fx[k]), WithinAbs(static_cast<double>(a.Fx[i]), allowed));
        CHECK_THAT(static_cast<double>(b.Fy[k]), WithinAbs(static_cast<double>(a.Fy[i]), allowed));
        CHECK_THAT(static_cast<double>(b.Fz[k]), WithinAbs(static_cast<double>(a.Fz[i]), allowed));
        CHECK_THAT(static_cast<double>(b.phi[k]), WithinAbs(static_cast<double>(a.phi[i]), allowed));
    }
    CHECK_THAT(static_cast<double>(b.Ec),
               WithinAbs(static_cast<double>(a.Ec), tolerance<TestType>() * std::abs(static_cast<double>(a.Ec)) + allowed));
}

// --------------------------- Compile-time switches ----------------------

TEMPLATE_TEST_CASE("Switches change what is computed, not the results", "[solver][switches]",
                   float, double) {
    const auto particles = random_system<TestType>(50, TestType{-1});
    const auto all = coulomb::compute<TestType, true, true>(particles);

    SECTION("forces only") {
        const auto r = coulomb::compute<TestType, false, false>(particles);
        CHECK(r.Fx == all.Fx);  // same calculation, so exactly equal
        CHECK(r.Fy == all.Fy);
        CHECK(r.Fz == all.Fz);
        CHECK(r.phi.empty());
        CHECK(r.Ec == TestType{0});
    }
    SECTION("potential without energy") {
        const auto r = coulomb::compute<TestType, true, false>(particles);
        CHECK(r.Fx == all.Fx);
        CHECK(r.phi == all.phi);
        CHECK(r.Ec == TestType{0});
    }
    SECTION("energy without potential") {
        const auto r = coulomb::compute<TestType, false, true>(particles);
        CHECK(r.Fx == all.Fx);
        CHECK(r.phi.empty());
        CHECK(r.Ec == all.Ec);
    }
}
