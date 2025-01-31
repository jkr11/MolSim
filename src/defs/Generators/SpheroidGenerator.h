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
  dvec3 origin_;
  /**
   * @brief this is an integer as it's the number of particles possible along
   * the radius of the spheroid, so we move along the radii in increments of h
   */
  const int radius_;
  /**
   * distance between particles
   */
  double h_;
  /**
   * mass of particles
   */
  double m_;
  /**
   * initial velocity of particles
   */
  const dvec3 initial_velocity_;
  /**
   * epsilon of particles
   */
  double epsilon_;
  /**
   * sigma of particles
   */
  double sigma_;
  /**
   * type of particles
   */
  const int type_;
  /**
   * TODO
   */
  double mv_;
  /**
   * only if this is passed with true, spheres will be two Dimensional
   */
  const bool two_d_;

 public:
  /**
   * @brief Instantiate SpheroidGenerator
   * @param origin Origin of the generator
   * @param radius Radius
   * @param h Distance between particles
   * @param m mass of particles
   * @param initial_velocity initial velocity of generated particles
   * @param epsilon epsilon of generated particles
   * @param sigma sigma of generated particles
   * @param type type of generated particles
   * @param mv TODO of generated particles
   * @param two_d whether generator is 2D or 3D
   */
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
