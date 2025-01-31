//
// Created by jkr on 12/26/24.
//

#ifndef INDEXFORCE_H
#define INDEXFORCE_H

#include <vector>

#include "defs/containers/ParticleContainer.h"

/**
 * @brief Holds the Index Force. Index force acts on cuboid coordinates specified
 * before the simulation
 */
class IndexForce {
  /**
   * The particle ids that are ultimately calculated in the membrane generator
   */
  std::vector<int> indeces_{};
  /**
   * The time until which the force is active
   */
  double time_{};
  /**
   * The value of the actual force that is then added
   */
  dvec3 force_values_{};

 public:
  /**
   * @brief Instantiate Index Force
   * @param ids Vector of ids, which the force acts upon
   * @param time Time until force is applied
   * @param force_values force vector which is added to targeted particles
   */
  explicit IndexForce(const std::vector<int>& ids, const double time,
                      const dvec3& force_values)
      : indeces_(ids), time_(time), force_values_(force_values) {}
  IndexForce() = default;

  /**
   * @brief Calculates the index force on the specified particle
   * @param p the particle on where the force is applied
   * @param sim_time the current time of the simulation
   * @return the force vector thats applied (0 if sime_time > time_)
   */
  dvec3 applyForce(Particle& p, double sim_time) const;

  /**
   * @brief Returns the ids on which the force has to act
   * @return the particle ID indices of the force
   */
  [[nodiscard]] std::vector<int> getIndices() const { return indeces_; }
};

#endif  // INDEXFORCE_H
