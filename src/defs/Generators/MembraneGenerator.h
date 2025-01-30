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
  std::vector<int> ids{};
  std::vector<ivec3> indeces{};

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
   * @param indeces
   */
  MembraneGenerator(const dvec3 &corner, const std::array<int, 3> &dimensions,
                    double h, double m,
                    const std::array<double, 3> &initialVelocity, double mv,
                    double epsilon, double sigma, int type, bool twoD,
                    const std::vector<ivec3> &indeces);

  /**
   * @brief generates particles in the shape of a cuboid
   * @param particles vector from container which contains the vector in which
   * this cuboid is saved in
   */
  void generate(std::vector<Particle> &particles) override;

  /**
   * @return particle ids
   */
  [[nodiscard]] std::vector<int> getIndeces() const;

  void setTargetIndeces(const std::vector<ivec3> &_indeces) {
    SpdWrapper::get()->info("Setting indeces");

    this->indeces = _indeces;
  }
};

#endif  // MEMBRANEGENERATOR_H
