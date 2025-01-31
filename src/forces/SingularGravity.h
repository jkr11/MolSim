//
// Created by jkr on 11/24/24.
//

#ifndef SINGULARGRAVITY_H
#define SINGULARGRAVITY_H
#include "SingularForce.h"

/**
 * @brief Globally acts on each particle using the gravitational constant g in
 * the specified axis direction [0 : x, 1 : y, 2 : z]
 */
class SingularGravity final : public SingularForce {
 private:
  /**
   * The gravitational scaling constant
   */
  double g_{};
  /**
   * The axis on which the force is applied {0 : x, 1 : y, 2 : z}
   */
  int axis_{};

 public:
  /**
   * @param g the gravitational coefficient
   * @param axis the axis 0 - 2 x - z on which gravity acts
   */
  explicit SingularGravity(double g, int axis);

  /**
   * @calculates f_i = p_i.mass * g * 1_{axis}
   * @param p the particle to act on
   * @return the force to be added to the particles force
   */
  [[nodiscard]] dvec3 applyForce(const Particle &p) const override;
};

#endif  // SINGULARGRAVITY_H
