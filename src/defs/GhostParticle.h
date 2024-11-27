//
// Created by maximilian on 27.11.24.
//

#pragma once

#include <array>
#include <string>

#include "defs/types.h"
#include "defs/Particle.h"

class GhostParticle final {
 private:
  /**
   * Position of the particle
   */
  std::array<double, 3> x{};

  // no force or velocity needed as it is deleted and newly created every time to ensure consistency

  /**
   * Mass of this particle
   */
  double m{};

  /**
   * depth of the potential well
   * used in the caluclation of lennard-jones force
   */
  double epsilon{};

  /**
   * distance from the particle at which the potential is zero \n
   * used in the calculation of lennard-jones force
   */
  double sigma{};

  Particle *ghostParticleOf{};

 public:
  GhostParticle(
      // for visualization, we need always 3 coordinates
      // -> in case of 2d, we use only the first and the second
      const std::array<double, 3> &x_arg,
      double m_arg, double _epsilon, double _sigma, Particle *_ghostParticleOf);

  ~GhostParticle();

  [[nodiscard]] const std::array<double, 3> &getX() const;

  [[nodiscard]] double getM() const;

  [[nodiscard]] double getEpsilon() const;

  [[nodiscard]] double getSigma() const;

  [[nodiscard]] static std::string toString();
};

