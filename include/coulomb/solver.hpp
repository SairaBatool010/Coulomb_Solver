// The Coulomb solver: data types and the computation of forces, potentials
// and energy.
//
// Everything here is a template, and templates must be visible wherever they
// are used, so the whole solver lives in this header.

#pragma once

#include <cmath>     // std::sqrt
#include <concepts>  // std::floating_point
#include <cstddef>   // std::size_t
#include <vector>    // std::vector

namespace coulomb {

template <std::floating_point Real>
struct Particle {
    Real q;        // charge
    Real x, y, z;  // position
};

// Everything the solver computes.
template <std::floating_point Real>
struct Results {
    std::vector<Real> Fx, Fy, Fz;  // force components per particle
    std::vector<Real> phi;         // potential per particle (empty if disabled)
    Real Ec{0};                    // total energy (stays 0 if disabled)
};

// Computes the forces, plus the potentials and energy if they are enabled.
//
//   F_i   = q_i * sum_{j != i} q_j * (x_i - x_j) / |x_i - x_j|^3   (equation 3)
//   phi_i = sum_{j != i} q_j / |x_i - x_j|                         (equation 4)
//   Ec    = 1/2 * sum_i q_i * phi_i                                 (equation 5)
//
// The energy needs each particle's potential even when the potentials are not
// requested. For that reason the potential is summed in a local variable
// (phi_i) whenever either feature is on.
template <std::floating_point Real, bool WithPotential, bool WithEnergy>
Results<Real> compute(const std::vector<Particle<Real>>& particles) {
    constexpr bool need_potential = WithPotential || WithEnergy;
    const std::size_t n = particles.size();

    Results<Real> results;
    results.Fx.assign(n, Real{0});
    results.Fy.assign(n, Real{0});
    results.Fz.assign(n, Real{0});
    if constexpr (WithPotential) {
        results.phi.assign(n, Real{0});
    }

    for (std::size_t i = 0; i < n; ++i) {
        const Particle<Real>& pi = particles[i];
        [[maybe_unused]] Real phi_i{0};

        for (std::size_t j = 0; j < n; ++j) {
            if (i == j) {
                continue;  // no self-interaction
            }

            const Particle<Real>& pj = particles[j];

            const Real dx = pi.x - pj.x;
            const Real dy = pi.y - pj.y;
            const Real dz = pi.z - pj.z;
            const Real r = std::sqrt(dx * dx + dy * dy + dz * dz);
            const Real r3 = r * r * r;

            results.Fx[i] += pi.q * pj.q * dx / r3;  // equation 3
            results.Fy[i] += pi.q * pj.q * dy / r3;
            results.Fz[i] += pi.q * pj.q * dz / r3;

            if constexpr (need_potential) {
                phi_i += pj.q / r;  // equation 4
            }
        }

        if constexpr (WithPotential) {
            results.phi[i] = phi_i;
        }
        if constexpr (WithEnergy) {
            // A plain 0.5 is a double, which would turn float math into double math
            results.Ec += static_cast<Real>(0.5) * pi.q * phi_i;  // equation 5
        }
    }
    return results;
}

}  // namespace coulomb
