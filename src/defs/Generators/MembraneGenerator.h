//
// Created by jkr on 10/31/24.
//

#ifndef MEMBRANEGENERATOR_H
#define MEMBRANEGENERATOR_H
#pragma once
#include "ParticleGenerator.h"
#include "utils/SpdWrapper.h"

/**
 * @brief Generates Particles to a container in the shape of a cuboid [dim1,
 * dim2, dim3]
 */
class MembraneGenerator final : public ParticleGenerator {
 private:
  dvec3 corner_;
  ivec3 dimensions_;
  double h_;
  double m_;
  const dvec3 initial_velocity_;
  double mv_;
  double epsilon_;
  double sigma_;
  const int type_;
  /**
   * here this just describes the behaviour of the brownian motion
   * initialization
   */
  const bool two_d_;
  std::vector<int> ids_{};
  std::vector<ivec3> indices_{};

 public:
  /**
   * @brief Constructor for the particle generator
   * @param corner Lower left corner / origin of the shape
   * @param dimensions number of particles in each unit direction
   * @param h distance between pairwise particles
   * @param m mass of the particles in the cuboid
   * @param initial_velocity velocity (imagine this as a net-zero movement of
   * all particles) given by the predefined maxwell-boltzmann generator
   * @param epsilon lj - epsilon
   * @param sigma lj - sigma
   * @param mv temperature of our system
   * @param type type of the particle in the system
   * @param two_d dimension of velocity vector of brownian motion
   * @param indices
   */
  MembraneGenerator(const dvec3 &corner, const std::array<int, 3> &dimensions,
                    double h, double m,
                    const std::array<double, 3> &initial_velocity, double mv,
                    double epsilon, double sigma, int type, bool two_d,
                    const std::vector<ivec3> &indices);

  /**
   * @brief generates particles in the shape of a cuboid, and connects each
   * particle to its direct neighbours. Furthermore, writes if each neighbour
   * is diagonal or not. Neighbours are directly saved by uint_ptr
   * @param particles vector from container which contains the vector in which
   * this cuboid is saved in
   */
  void generate(std::vector<Particle> &particles) override;

  /**
   * @return particle ids (so matching by p.getId())
   */
  [[nodiscard]] std::vector<int> getIndices() const;

  /**
   * @brief calculates the particle ids from the cuboid coordinates
   * @param indices the cuboid coordinates of the membrane
   */
  void setTargetIndices(const std::vector<ivec3> &indices) {
    this->indices_ = indices;
  }
};

#endif  // MEMBRANEGENERATOR_H
