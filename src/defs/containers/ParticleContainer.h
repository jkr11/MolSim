//
// Created by jkr on 10/18/24.
//
#pragma once

#include <cstddef>
#include <functional>
#include <vector>

#include "defs/Particle.h"

/**
 * @brief Interface for an object storing particles while providing single and
 * pairwise iterators.
 */
class ParticleContainer {
 public:
  /**
   * @brief Virtual destructor to allow proper cleanup of derived classes.
   */
  virtual ~ParticleContainer() = default;

  /**
   * @brief Add a particle to the particle system.
   * @param p Particle to be added.
   */
  virtual void addParticle(const Particle& p) = 0;

  virtual void addParticles(const std::vector<Particle>& particles) = 0;
  /**
   * @brief Remove a particle from the particle system.
   * @param p Particle to be removed.
   */
  virtual void removeParticle(const Particle& p) = 0;

  /**
   * @brief Get a vector containing references to all particles in the
   * container.
   * @return Vector of all particles.
   */
  [[nodiscard]] virtual std::vector<Particle> getParticles() const = 0;

  /**
   * @brief Get the number of particles in the container.
   * @return Number of particles in the container.
   */
  [[nodiscard]] virtual std::size_t size() const = 0;

  /**
   * @brief Apply a function to each single particle.
   * @param f Lambda to be applied to each particle (p).
   */
  virtual void singleIterator(const std::function<void(Particle&)>& f) = 0;

  /**
   * @brief Apply a function to each unique pair of particles (i.e., (p1, p2) ==
   * (p2, p1)).
   * @param f Lambda to be applied to each particle pair (p1, p2).
   */
  virtual void pairIterator(
      const std::function<void(Particle&, Particle&)>& f) = 0;
};
