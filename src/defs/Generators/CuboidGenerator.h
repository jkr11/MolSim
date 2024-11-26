//
// Created by jkr on 10/31/24.
//

#ifndef CUBOIDGENERATOR_H
#define CUBOIDGENERATOR_H
#include "ParticleGenerator.h"

/**
 * @brief Generates Particles to a container in the shape of a cuboid [dim1,
 * dim2, dim3]
 */
class CuboidGenerator final : public ParticleGenerator {
 private:
  dvec3 corner;
  ivec3 dimensions;
  double h;
  double m;
  const dvec3 initialVelocity;
  double mv;
  /**
   * technically these are only relevant for week 4 but here they dont really
   * matter right now
   */
  double epsilon;
  double sigma;
  const int type{};
  /**
   * here this just describes the behaviour of the brownian motion
   * initialization
   */
  const bool twoD{};

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
   * @param twoD dimension of velocity vector of brownian motion
   */
  CuboidGenerator(const dvec3 &corner, const std::array<int, 3> &dimensions,
                  double h, double m,
                  const std::array<double, 3> &initialVelocity, double mv,
                  double epsilon, double sigma, int type, bool twoD);

  /**
   * @brief generates particles in the shape of a cuboid
   * @param particles vector from container which contains the vector in which
   * this cuboid is saved in
   */
  void generate(std::vector<Particle> &particles) override;
};

#endif  // CUBOIDGENERATOR_H
