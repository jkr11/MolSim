//
// Created by jkr on 10/18/24.
//
#pragma once
#include <cstddef>
#include <list>
#include <vector>

#include "Particle.h"

/**
 * @brief A wrapper for a std::vector<Particle> providing single and pairwise
 * iterators
 */
class ParticleContainer {
 private:
  std::vector<Particle> particles;

 public:
  /**
   * @brief Initialize empty particle container
   */
  ParticleContainer();

  /**
   * @brief Create particle container from vector of particles
   * @param particles is the vector of Particles to be used
   * @note A new vector is initialized to hold the references to the passed
   * particles
   */
  explicit ParticleContainer(const std::vector<Particle>& particles);

  /**
   * @brief Create particle container from list of particles
   * @param particles is the list of Particles to be used
   * @note A new vector is initialized to hold the references to the passed
   * particles
   */
  explicit ParticleContainer(const std::list<Particle>& particles);

  /**
   * @brief Destructor
   */
  ~ParticleContainer();

  /**
   * @brief Add particle into the particle-system
   * @param p particle to be added
   */
  void addParticle(const Particle& p);

  /**
   * @brief Get the beginning of the iterator
   * @return Iterator start
   */
  auto begin();

  /**
   * @brief Get the end of the iterator
   * @return Iterator end
   */
  auto end();

  /**
   * @brief Access particles using the array operator
   * @param n index of particle to be accessed
   * @return particle at index
   */
  Particle operator[](int n);

  /**
   * @brief Get vector of all particles in the container
   * @return vector of all particles
   */
  std::vector<Particle> getParticles();

  /**
   * @brief Get reference to vector of all particles in the container
   * @return reference to vector of all particles
   */
  std::vector<Particle>& getParticlesReference();

  /**
   * @brief Get the number of particles in the container
   * @return number of particles in the container
   */
  [[nodiscard]] std::size_t size() const;

  /**
   * @brief setter for the vector of particles
   * @param particles particles to be set in the container
   */
  void setParticles(const std::vector<Particle>& particles);

  /**
   * @brief resizes particles to the desires size, allows for faster additions
   * @param size desired size of the vector
   */
  void resize(std::size_t size);

  /**
   * @brief Iterator over single particles p
   * @tparam UnOp Single parameter lambda taking a particle as input
   * @param f Lambda that's applied to (p)
   */
  // this is defined in the .h in order for it to be instantiable everywhere
  template <typename UnOp>
  void single_iterator(const UnOp& f) {
    for (auto& p : particles) {
      f(p);
    }
  }

  /**
   * @brief Iterator over unique pairs (i.e. (p1, p2) == (p2, p1))
   * @tparam BinOp Double parameter lambda taking two particles as input
   * @param f Lambda that's applied to (p1, p2)
   */
  // this is defined in the .h in order for it to be instantiable everywhere
  template <typename BinOp>
  void pairIterator(const BinOp& f) {
    // note that the upper tri-diag matrix is iterated over
    for (size_t i = 0; i < particles.size(); ++i) {
      for (size_t j = i + 1; j < particles.size(); ++j) {
        f(particles[i], particles[j]);
      }
    }
  }
};
