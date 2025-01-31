//
// Created by jkr on 10/18/24.
//
#pragma once

#include <cstddef>
#include <functional>
#include <vector>

#include "defs/Particle.h"
#include "forces/InteractiveForce.h"
#include "forces/SingularForce.h"

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

  virtual void addParticles(const std::vector<Particle>& particles) = 0;
  /**
   * @brief Remove a particle from the particle system.
   * DO NOT USE, JUST FOR TESTING!
   * @param p Particle to be removed.
   */
  virtual void removeParticle(const Particle& p) = 0;

  /**
   * @brief Get a vector containing references to all particles in the
   * container.
   * @return Vector of all particles.
   */
  [[nodiscard]] virtual std::vector<Particle*> getParticles() = 0;

  /**
   * @brief Get a vector containing particles
   * @return Vector of all particles.
   */
  [[nodiscard]] virtual std::vector<Particle> getParticlesObjects() = 0;

  /**
   * @brief Get the number of particles in the container.
   * @return Number of particles in the container.
   */
  [[nodiscard]] virtual std::size_t size() const = 0;

  /**
   * @brief Impose the invariant, that the particles are spatially sorted into
   * the correct vectors and apply boundary condition
   */
  virtual void imposeInvariant() = 0;

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

  /**
   * @brief Compute interactive forces using the C18 coloring scheme
   * (Parallelization strategy 1)
   */
  virtual void computeInteractiveForcesC18(
      const std::vector<std::unique_ptr<InteractiveForce>>&
          interactive_forces) = 0;

  /**
   * @brief Compute interactive forces using Force Buffers (Parallalization
   * strategy 2)
   */
  virtual void computeInteractiveForcesForceBuffer(
      const std::vector<std::unique_ptr<InteractiveForce>>&
          interactive_forces) = 0;

  /**
   * @brief Compute singular forces
   */
  virtual void computeSingularForces(
      const std::vector<std::unique_ptr<SingularForce>>& singular_forces) = 0;

  /**
   * @brief Returns the kinetic energy of the system, E_kin = 1/2 \sum_i^n m_i *
   * <v_i | v_i>
   * @return kinetic energy of the system
   */
  virtual double getKineticEnergy() = 0;

  /**
   * @brief Increment Time for index forces
   */
  virtual void incrementTime() { this->current_time++; }

  /**
   * Current time for index forces
   */
  double current_time = 0;

  /**
   * @brief the exact number of current particles, updated accordingly
   * @return the current count of particles left in the simulation
   */
  [[nodiscard]] virtual size_t getParticleCount() = 0;

  /**
   * @brief the exact number of current special particles, updated accordingly
   * @return the current count of special particles left in the simulation
   */
  [[nodiscard]] virtual size_t getSpecialParticleCount() = 0;

  /**
   * @brief returns the domain of the container
   * @return the domain of the container
   */
  virtual ivec3 getDomain() = 0;
};
