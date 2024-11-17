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
  dvec3 origin{};
  /**
   * @brief this is an integer as it's the number of particles possible along
   * the radius of the spheroid
   */
  const int radius;
  double h{};
  double m{};
  const dvec3 initialVelocity;
  double mv{};
  double epsilon{};
  double sigma{};
  const int type{};
  const bool threeD{};

 public:
  SpheroidGenerator(const dvec3 &origin, int radius, double h, double m,
                    const dvec3 &initialVelocity, double epsilon, double sigma,
                    int type, bool threeD);

  /**
   * @brief generates a 1 or 2 sphere of particles, vector size is approximated
   * by volume of 2-sphere as disk is contained.
   * @param particles particles into which the spheroid is place
   */
  void generate(std::vector<Particle> &particles) override;
};

#endif  // SPHEREGENERATOR_H
