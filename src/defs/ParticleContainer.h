//
// Created by jkr on 10/18/24.
//
#pragma once
#include "Particle.h"
#include <cstddef>
#include <list>
#include <vector>


/**
 * @brief Particle Container managing a vector of Particles and providing two functions for single and pair iteration
 */
class ParticleContainer {
private:
  std::vector<Particle> particles;

public:
  ParticleContainer();

  explicit ParticleContainer(const std::vector<Particle>& particles);

  explicit ParticleContainer(const std::list<Particle>& particles);

  ~ParticleContainer();

  void addParticle(const Particle& p) {
    particles.push_back(p);
  }

  auto begin() {
    return particles.begin();
  }

  auto end() {
    return particles.end();
  }

  auto operator[](const int n) {
    return particles[n];
  }

  auto getParticles() {
    return particles;
  }

  [[nodiscard]] std::size_t size() const {
    return particles.size();
  }

  /**
   * Templating is faster than directly passing a function or a reference to one as it avoids copying and allowins inlining
   * @tparam UnOp Single parameter function taking a Particle as input
   * @param f instantiation of the UnOp Type
   */
  template <typename UnOp>
  void single_iterator(const UnOp &f) {
    for (auto& p : particles) {
      f(p);
    }
  }

  // From what I understand, passing a std::function<void(Particle&, Particle&)> &f here prevents inlining, so we send an anonymous callable thats able to be inlined
  // probably doesn't even matter
  template <typename BinOp>
  void pairIterator(const BinOp &f) {
    for (size_t i = 0; i < particles.size(); ++i) {
      for (size_t j = i + 1; j < particles.size(); ++j) {
        f(particles[i], particles[j]);
      }
    }
  }
};
