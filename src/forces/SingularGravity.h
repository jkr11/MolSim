//
// Created by jkr on 11/24/24.
//

#ifndef SINGULARGRAVITY_H
#define SINGULARGRAVITY_H
#include "SingularForce.h"

/**
 * @brief Globally acts on each particle using the gravitational constant g in
 * the y direction
 */
class SingularGravity final : public SingularForce {
 private:
  double g{};
  int axis{};

 public:
  /**
   *
   * @param g the gravitational coefficient
   * @param axis the axis 0 - 2 x - z on which gravity acts
   */
  explicit SingularGravity(double g, int axis);

  // TODO: do we need to do this modular over the axis? or is it only in y
  // direction?
  /**
   * @calculates f_i = p_i.mass * [0,g,0]
   * @param p the particle to act on
   * @return the force to be added to the particles force
   */
  [[nodiscard]] dvec3 applyForce(const Particle &p, ParticleContainer& container) const override;
};

#endif  // SINGULARGRAVITY_H
