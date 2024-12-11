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

 public:
  /**
   *
   * @param g the gravitational coefficient
   */
  explicit SingularGravity(double g);

  // TODO: do we need to do this modular over the axis? or is it only in y
  // direction?
  /**
   * @calculates f_i = p_i.mass * [0,g,0]
   * @param p the particle to act on
   * @return the force to be added to the particles force
   */
  [[nodiscard]] inline dvec3 applyForce(const Particle &p) const override {
    return {0.0, g * p.getM(), 0.0};
  }
};

#endif  // SINGULARGRAVITY_H
