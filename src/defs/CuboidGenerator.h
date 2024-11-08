//
// Created by jkr on 10/31/24.
//

#ifndef CUBOIDGENERATOR_H
#define CUBOIDGENERATOR_H
#include "ParticleGenerator.h"

/**
 * Generates Particles to a container in the shape of a cuboid [dim1, dim2,
 * dim3]
 */
class CuboidGenerator final : public ParticleGenerator {
 private:
  dvec3 corner;
  std::array<int, 3> dimensions;
  double h;
  double m;
  const dvec3 initialVelocity;
  double mv;
  double epsilon;
  double sigma;
  const int type{};

 public:
  /**
   * @brief Constructor for the particle generator
   * @param corner Lower left corner / origin of the shape
   * @param dimensions number of particles in each unit direction
   * @param h distance between pairwise particles
   * @param m mass of the particles in the cuboid
   * @param initialVelocity velocity (imagine this as a net-zero movement of all
   * particles) given by the predefined maxwell-boltzmann generator
   * @param epsilon lj - epsilon
   * @param sigma lj - sigma
   * @param mv temperature of our system
   * @param type type of the particle in the system
   */
  CuboidGenerator(const dvec3 &corner, const std::array<int, 3> &dimensions,
                  const double h, const double m,
                  const std::array<double, 3> &initialVelocity, const double mv,
                  const double epsilon, const double sigma, const int type);

  /**
   * @brief generates particles in the shape of a cuboid
   * @param particles vector from container which contains the vector in which
   * this cuboid is saved in
   */
  void generate(std::vector<Particle> &particles) override;
};

#endif  // CUBOIDGENERATOR_H
