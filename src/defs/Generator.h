//
// Created by jkr on 10/31/24.
//
#pragma once
#ifndef GENERATOR_H
#define GENERATOR_H
#include "Particle.h"
#include "ParticleContainer.h"
class ParticleGenerator {
private:
  dvec3 corner;
  std::array<int, 3> dimensions;
  double h;
  double m;
  const dvec3 initialVelocity;
  double temp;
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
   * @param temperature temperature of our system
   * @param type type of the particle in the system
   */
  ParticleGenerator(const dvec3 &corner, const std::array<int, 3>& dimensions,
                    const double h, const double m,
                    const std::array<double, 3> &initialVelocity,
                    const double temperature, const int type);

  /**
   * @brief spawns particles in a given shape
   * @param container container holding the vector of particles
   */
  void generate(ParticleContainer &container);

};
#endif // GENERATOR_H
