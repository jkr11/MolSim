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
  ParticleContainer();
  explicit ParticleContainer(const std::vector<Particle>& particles);
  explicit ParticleContainer(const std::list<Particle>& particles);

  ~ParticleContainer();

  /**
   * Add particle into the particle-system
   * @param p particle to be added
   */
  void addParticle(const Particle& p);

  /**
   * @return Iterator start
   */
  auto begin();

  /**
   * @return Iterator end
   */
  auto end();

  /**
   * @param n index of particle to be accessed
   * @return particle at index
   */
  Particle operator[](int n);

  /**
   * @return returns all particles
   */
  std::vector<Particle> getParticles();

  /**
   * @return number of particles in the container
   */
  [[nodiscard]] std::size_t size() const;

  /**
   * Iterator over single particles p
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
   * Iterator over unique pairs (i.e. (p1, p2) == (p2, p1))
   * @tparam BinOp Double parameter lambda taking two particles as input
   * @param f Lambda that's applied to (p1, p2)
   */
  // this is defined in the .h in order for it to be instantiable everywhere
  template <typename BinOp>
  void pairIterator(const BinOp& f) {
    // note that the upper tri-diag matrix is iterated over
    for (size_t i = 0; i < particles.size();
         ++i) {  // TODO: I assume the last iteration of the outer loop is
                 // unnecessary
      for (size_t j = i + 1; j < particles.size(); ++j) {
        f(particles[i], particles[j]);
      }
    }
  }
};
