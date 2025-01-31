//
// Created by jkr on 12/23/24.
//
#ifndef HARMONICFORCE_H
#define HARMONICFORCE_H
#include <cmath>

#include "SingularForce.h"

/**
 * @brief Harmonic force for inter particle forces in a membrane
 */
class HarmonicForce final : public SingularForce {
 private:
  int k_{};
  double r_0_{};
  /**
   * sqrt(2) * bond_length == diag. bond length
   */
  double sr_0_{};

 public:
  /**
   * @brief Instantiates Harmonic Force
   * @param k Spring constant
   * @param r_0 average bond length
   */
  explicit HarmonicForce(const int k, const double r_0) : k_(k), r_0_(r_0) {
    sr_0_ = std::sqrt(2) * r_0;
  }

  /**
   * @brief Calculate force applied by the neighbours
   * @param p particle
   * @return force vector
   */
  [[nodiscard]] dvec3 applyForce(const Particle& p) const override;
};

#endif  // HARMONICFORCE_H
