//
// Created by jkr on 11/17/24.
//

#ifndef SPHEREGENERATOR_H
#define SPHEREGENERATOR_H
#include "ParticleGenerator.h"

/**
 * @brief generates particles in a 1 or 2 sphere
 */
class SpheroidGenerator final : public ParticleGenerator {
 private:
  dvec3 origin_{};
  /**
   * @brief this is an integer as it's the number of particles possible along
   * the radius of the spheroid, so we move along the radii in increments of h
   */
  const int radius_;
  double h_{};
  double m{};
  const dvec3 initial_velocity;
  double epsilon{};
  double sigma{};
  const int type{};
  double mv{};
  /**
   * only if this is passed with true, spheres will be two Dimensional
   */
  const bool two_d{};

 public:
  SpheroidGenerator(const dvec3 &origin, int radius, double h, double m,
                    const dvec3 &initial_velocity, double epsilon, double sigma,
                    int type, double mv, bool two_d);

  /**
   * @brief generates a 1 or 2 sphere of particles, vector size is approximated
   * by volume of 2-sphere as disk is contained.
   * @param particles particles into which the spheroid is placed
   */
  void generate(std::vector<Particle> &particles) override;
};

#endif  // SPHEREGENERATOR_H
