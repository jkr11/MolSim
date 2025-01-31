//
// Created by jkr on 12/25/24.
//

#ifndef TRUNCATEDLENNARDJONES_H
#define TRUNCATEDLENNARDJONES_H
#include "InteractiveForce.h"

/**
 * @brief Repulsive Lennard Jones force that is truncated at 2^(1/6) * sigma
 */
class TruncatedLennardJones final : public InteractiveForce {
 public:
  /**
   * @brief Default construtor for truncated Lennard Jones force
   */
  explicit TruncatedLennardJones() = default;

  /**
   * @brief calculates the repulsive part of the lennard jones force between two
   * particles if they are close enough, 0 else.
   * @param p1 source particle
   * @param p2 target particle
   * @return force from particle 1 to particle 2
   */
  dvec3 directionalForce(Particle& p1, Particle& p2) const override;
};

#endif  // TRUNCATEDLENNARDJONES_H
